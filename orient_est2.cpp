// functions governing orientation estimation and interpolation
#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
// MARKED - check later
void est2orient::init_orient(int setblksz, int setprec, int setmaxshift)
{
  ofield.init_orient(setblksz,setprec,setmaxshift);
  edgeintensity = new double[h*w];
  edgeangle = new double[h*w];
  blk_intensity = new double[ofield.numblks];
  blk_angle = new double[ofield.numblks];
  return;
}
// MARKED - check later
void est2orient::transpose()
{
  dwtnode::transpose();
  return;
}
void est2orient::calc_raw_edges()
{
  int x,y;
  double *horzedges = new double[h*w];
  double *vertedges = new double[h*w];
  const double sqrt2 = 1.41421356;
  double *tmp = new double[h*w];
  // find horizontal edges
  for (y=1;y<h-1;y++)
    for (x=0;x<w;x++)
    { // vertical filter [1,0,-1] for edge detection
      tmp[y*w+x] = pixels[(y-1)*w+x] - pixels[(y+1)*w+x];
    }
  for (y=1;y<h-1;y++)
    for (x=1;x<w-1;x++)
    { // horizontal smoothing filter [1,sqrt2,1]
      horzedges[y*w+x] = tmp[y*w+x-1] + sqrt2*tmp[y*w+x] + tmp[y*w+x+1];
    }
  // find vertical edges
  for (y=0;y<h;y++)
    for (x=1;x<w-1;x++)
      tmp[y*w+x] = pixels[y*w+x-1] - pixels[y*w+x+1];
  for (y=1;y<h-1;y++)
    for (x=1;x<w-1;x++)
      vertedges[y*w+x] = tmp[(y-1)*w+x] + sqrt2*tmp[y*w+x] + tmp[(y+1)*w+x];
  delete[] tmp;
  for (x=0;x<w;x++)
  { // set top and bottom edges to 0
    edgeintensity[x] = 0;
    edgeintensity[(h-1)*w+x] = 0;
  }
  for (y=0;y<h;y++)
  { // set left and bottom edges to 0
    edgeintensity[y*w] = 0;
    edgeintensity[y*w+(w-1)] = 0;
  }
  // calc edge intensity
  for (y=1;y<h-1;y++)
    for (x=1;x<w-1;x++)
      edgeintensity[y*w+x] = sqrt(horzedges[y*w+x]*horzedges[y*w+x]
        + vertedges[y*w+x]*vertedges[y*w+x]);
  // calculate edge angle
  const double halfpi = 1.5707963268;
  for (int n=0;n<h*w;n++)
  {
    if (vertedges[n]==0)
      edgeangle[n] = 0;
    else if (horzedges[n]==0)
      edgeangle[n] = halfpi;
    else
      edgeangle[n] = atan(vertedges[n]/horzedges[n]);
  }
  delete[] horzedges;
  delete[] vertedges;
  return;
}
void est2orient::calc_edge_blocks()
{
  const int fieldw = w/ofield.blksz;
  for (int nblk=0;nblk<ofield.numblks;nblk++)
  {
    // location of the topleft pixel in the block
    int xstart = (nblk%fieldw)*ofield.blksz;
    int ystart = (nblk/fieldw)*ofield.blksz;
    // find average edge angle over the block
    blk_intensity[nblk] = 0;
    blk_angle[nblk] = 0;
    for (int x=xstart; x<xstart+ofield.blksz; x++)
      for (int y=ystart; y<ystart+ofield.blksz; y++)
      {
        int n = y*w+x;
        blk_angle[nblk]     += edgeangle[n]*edgeintensity[n];
        blk_intensity[nblk] += edgeintensity[n];
      }
    // normalise for by weighted sum and block size respectively
    blk_angle[nblk]     /= blk_intensity[nblk];
    blk_intensity[nblk] /= (ofield.blksz*ofield.blksz);
  }
  return;
}
void est2orient::edge_csvout()
{
  ofstream fout("C:\\Program Files\\Matlab7\\work\\"
    "oriented wavelets\\edgemap.csv",ios::binary);
  for (int n=0;n<h*w;n++)
  {
    fout << edgeintensity[n];
    if (n%w == w-1)
      fout << endl;
    else
      fout << ',';
  }
  fout.close();
  fout.open("C:\\Program Files\\Matlab7\\work\\"
    "oriented wavelets\\edgeangle.csv",ios::binary);
  for (int n=0;n<h*w;n++)
  {
    fout << edgeangle[n];
    if (n%w == w-1)
      fout << endl;
    else
      fout << ',';
  }
  fout.close();
  return;
}
double est2orient::fetch_residual(int orientblk, direction dir, char shift)
{
  const double halfpi = 1.5707963268;
  const double pi = 2*halfpi;
  // find angle (in radians) corresponding to chosen shift
  double shiftangle;
  if (dir == vertical)
    shiftangle = atan(ofield.oprec/((double)shift));
  else // oriented horizontal transform
    shiftangle = atan(((double)shift)/ofield.oprec);
  double anglediff = abs(blk_angle[orientblk] - shiftangle);
  if (anglediff > halfpi)
    anglediff = abs(anglediff - pi);
  return blk_intensity[orientblk]*anglediff;
}
double est2orient::fetch_loss(int n, direction dir, char shift)
{
  if (dir==both)
  {
    double vloss = fetch_loss(n, vertical, shift);
    double hloss = fetch_loss(n, horizontal, shift);
    return min(vloss,hloss);
  }
  const double alpha = 0.05; // continuity penalty weight
  const double beta  = 0.2; // switching penalty weight
  const double residual = fetch_residual(n,dir,shift);
  const int fieldh = h/ofield.blksz;
  const int fieldw = w/ofield.blksz;
  double cont_penalty   = 0;
  double switch_penalty = 0;
  if (dir==horizontal) // vertical shifts
  {
    // penalise against the horizontal neighbours
    if (n%fieldw != 0) // left edge of field
    {
      if (ofield.orientvec[n-1].hshift != 0)
        switch_penalty += blk_intensity[n-1];
      else cont_penalty += abs(shift - ofield.orientvec[n-1].vshift)
        *blk_intensity[n-1]/ofield.oprec;
    }
    if (n%fieldw != (fieldw-1)) // right edge of field
    {
      if (ofield.orientvec[n+1].hshift != 0)
        switch_penalty += blk_intensity[n+1];
      else cont_penalty += abs(shift - ofield.orientvec[n+1].vshift)
        *blk_intensity[n+1]/ofield.oprec;
    }
    // penalise against the vertical neighbours
    if (n >= fieldw) // skip if at top edge of field
    {
      if (ofield.orientvec[n-fieldw].hshift != 0)
        switch_penalty += blk_intensity[n-fieldw];
      else cont_penalty += abs(shift - ofield.orientvec[n-fieldw].vshift)
        *blk_intensity[n-fieldw]/ofield.oprec;
    }
    if (n < (fieldh-1)*fieldw) // bottom edge of field
    {
      if (ofield.orientvec[n+fieldw].hshift != 0)
        switch_penalty += blk_intensity[n+fieldw];
      else cont_penalty += abs(shift - ofield.orientvec[n+fieldw].vshift)
        *blk_intensity[n+fieldw]/ofield.oprec;
    }
  }
  else // (dir==vertical)
  {
    if (n%fieldw != 0) // left edge of field
    {
      if (ofield.orientvec[n-1].vshift != 0)
        switch_penalty += blk_intensity[n-1];
      else cont_penalty += abs(shift - ofield.orientvec[n-1].hshift)
        *blk_intensity[n-1]/ofield.oprec;
    }
    if (n%fieldw != (fieldw-1)) // right edge of field
    {
      if (ofield.orientvec[n+1].vshift != 0)
        switch_penalty += blk_intensity[n+1];
      else cont_penalty += abs(shift - ofield.orientvec[n+1].hshift)
        *blk_intensity[n+1]/ofield.oprec;
    }
    if (n >= fieldw) // skip if at top edge of field
    {
      if (ofield.orientvec[n-fieldw].vshift != 0)
        switch_penalty += blk_intensity[n-fieldw];
      else cont_penalty += abs(shift - ofield.orientvec[n-fieldw].hshift)
        *blk_intensity[n-fieldw]/ofield.oprec;
    }
    if (n < (fieldh-1)*fieldw) // bottom edge of field
    {
      if (ofield.orientvec[n+fieldw].vshift != 0)
        switch_penalty += blk_intensity[n+fieldw];
      else cont_penalty += abs(shift - ofield.orientvec[n+fieldw].hshift)
        *blk_intensity[n+fieldw]/ofield.oprec;
    }
  }
  cont_penalty   *= alpha;
  switch_penalty *= beta;
  return residual+cont_penalty+switch_penalty;
}
  //const int ybase = (n/fieldw)*ofield.blksz;
  //const int xbase = (n%fieldw)*ofield.blksz;
  //const int nbrsize = 4; // # of pixels beyond a 2x2 block
  //                       // eg 4 -> 10x10 pixel neighbourhood
  //if (dir==vertical)
  //{
  //  for (int y = -nbrsize; y<=1+nbrsize; y++)
  //    for (int x = -nbrsize; x<=1+nbrsize; x++)
  //    {
  //      cont_penalty += abs(shift - ofield.retrieve(y+ybase,x+xbase,vertical));
  //      if (ofield.retrieve(y+ybase,x+xbase,horizontal)!=0)
  //        switch_penalty += abs(shift);
  //    }
  //}
  //else // dir==horizontal
  //{
  //  char nbrshift;
  //  for (int y = -nbrsize; y<=1+nbrsize; y++)
  //    for (int x = -nbrsize; x<=1+nbrsize; x++)
  //    {
  //      cont_penalty += abs(shift - ofield.retrieve(y+ybase,x+xbase,horizontal));
  //      nbrshift = ofield.retrieve(y+ybase,x+xbase,vertical);
  //      if (nbrshift!=0)
  //        switch_penalty += abs(nbrshift);
  //    }
  //}
bool est2orient::test_loss(int n, direction dir, char shift, double &best)
{
  if (dir == both)
  {
    cerr << "test_loss direction argument can only be horizontal or vertical" << endl;
    exit(1);
  }
  const double alpha = 0.05;  // continuity penalty weight
  const double beta  = 0.2; // switching penalty weight
  const double residual = fetch_residual(n, dir, shift);
  if (residual >= best) // early exit
    return false;
  const int fieldh = h/ofield.blksz;
  const int fieldw = w/ofield.blksz;
  double cont_penalty   = 0;
  double switch_penalty = 0;
  if (dir==horizontal) // vertical shifts
  {
    // penalise against the horizontal neighbours
    if (n%fieldw != 0) // left edge of field
    {
      if (ofield.orientvec[n-1].hshift != 0)
        switch_penalty += blk_intensity[n-1];
      else cont_penalty += abs(shift - ofield.orientvec[n-1].vshift)
        *blk_intensity[n-1]/ofield.oprec;
    }
    if (n%fieldw != (fieldw-1)) // right edge of field
    {
      if (ofield.orientvec[n+1].hshift != 0)
        switch_penalty += blk_intensity[n+1];
      else cont_penalty += abs(shift - ofield.orientvec[n+1].vshift)
        *blk_intensity[n+1]/ofield.oprec;
    }
    // penalise against the vertical neighbours
    if (n >= fieldw) // skip if at top edge of field
    {
      if (ofield.orientvec[n-fieldw].hshift != 0)
        switch_penalty += blk_intensity[n-fieldw];
      else cont_penalty += abs(shift - ofield.orientvec[n-fieldw].vshift)
        *blk_intensity[n-fieldw]/ofield.oprec;
    }
    if (n < (fieldh-1)*fieldw) // bottom edge of field
    {
      if (ofield.orientvec[n+fieldw].hshift != 0)
        switch_penalty += blk_intensity[n+fieldw];
      else cont_penalty += abs(shift - ofield.orientvec[n+fieldw].vshift)
        *blk_intensity[n+fieldw]/ofield.oprec;
    }
  }
  else // (dir==vertical)
  {
    if (n%fieldw != 0) // left edge of field
    {
      if (ofield.orientvec[n-1].vshift != 0)
        switch_penalty += blk_intensity[n-1];
      else cont_penalty += abs(shift - ofield.orientvec[n-1].hshift)
        *blk_intensity[n-1]/ofield.oprec;
    }
    if (n%fieldw != (fieldw-1)) // right edge of field
    {
      if (ofield.orientvec[n+1].vshift != 0)
        switch_penalty += blk_intensity[n+1];
      else cont_penalty += abs(shift - ofield.orientvec[n+1].hshift)
        *blk_intensity[n+1]/ofield.oprec;
    }
    if (n >= fieldw) // skip if at top edge of field
    {
      if (ofield.orientvec[n-fieldw].vshift != 0)
        switch_penalty += blk_intensity[n-fieldw];
      else cont_penalty += abs(shift - ofield.orientvec[n-fieldw].hshift)
        *blk_intensity[n-fieldw]/ofield.oprec;
    }
    if (n < (fieldh-1)*fieldw) // bottom edge of field
    {
      if (ofield.orientvec[n+fieldw].vshift != 0)
        switch_penalty += blk_intensity[n+fieldw];
      else cont_penalty += abs(shift - ofield.orientvec[n+fieldw].hshift)
        *blk_intensity[n+fieldw]/ofield.oprec;
    }
  }
  cont_penalty   *= alpha;
  switch_penalty *= beta;
  double total_loss = residual + cont_penalty + switch_penalty;
  bool thresh_met = (total_loss < best);
  if (thresh_met)
    best = total_loss;
  return thresh_met;
}
// set shift closest to the edge angle 
void est2orient::set_shift(int n, bool forcedir, direction dir)
{
  const double qpi = 0.785398163; // pi/4
  const double halfpi = 1.5707963268;
  if (!forcedir)
  {
    if ((edgeangle[n]<qpi)&&(edgeangle[n]>-qpi))
      dir = horizontal;
    else dir = vertical;
  }
  if (dir == horizontal)
    ofield.orientvec[n].vshift = (char) round(edgeangle[n]/qpi*ofield.oprec);
  else
  {
    if (edgeangle[n]>0)
      ofield.orientvec[n].hshift = (char) round((halfpi-edgeangle[n])/qpi*ofield.oprec);
    else
      ofield.orientvec[n].hshift = (char) round((-halfpi-edgeangle[n])/qpi*ofield.oprec);
  }
  return;
}
void est2orient::choose_shift(int orientblk)
{
  double best = fetch_loss(orientblk,both,0);
  for (int sigma=-ofield.maxshift;sigma<=ofield.maxshift;sigma++)
  {
    if (test_loss(orientblk,vertical,sigma,best))
    {
      ofield.orientvec[orientblk].hshift = sigma;
      ofield.orientvec[orientblk].vshift = 0;
    }
    if (test_loss(orientblk,horizontal,sigma,best))
    {
      ofield.orientvec[orientblk].hshift = 0;
      ofield.orientvec[orientblk].vshift = sigma;
    }
  }
  return;
}
// seeding pass chooses shifts on the transform with the lowest
// zero-shift energy, but does not allow the other direction to be altered.
// subsequent passes allow either direction to be chosen, even while the
// threshold requirement is loosened
void est2orient::choose_orient()
{
  // gain of the edge detection filter
  const double gain = 2*(sqrt(2.0)+2);
  calc_raw_edges();
  calc_edge_blocks();
  double tarray[] = {40, 35, 30, 25, 20};
  //first pass
  cout << "pass 1" << "... ";
  for (int n=0;n<ofield.numblks;n++)
    if (blk_intensity[n] > gain*tarray[0])
      set_shift(n);
  //subsequent passes
  for (int i=1;i<5;i++)
  {
    double thresh = gain*tarray[i];
    cout << "pass " << (i+1) << "... ";
    for (int n=0;n<ofield.numblks;n++)
      if (blk_intensity[n] > thresh)
        choose_shift(n);
  }
  cout << endl;
  return;
}