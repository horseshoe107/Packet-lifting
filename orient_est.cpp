#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
extern double splines[];
extern const int splines_extent;
estorient::estorient(dwtnode *target):dwtnode(target->h,target->w,target->dwtbase)
{
  dwtlevel[vertical]=target->dwtlevel[vertical];
  dwtlevel[horizontal]=target->dwtlevel[horizontal];
  for (int n=0;n<h*w;n++)
    pixels[n] = target->pixels[n];
}
void estorient::init_orient(int setblksz, int setprec, int setmaxshift)
{
  ofield.init_orient(setblksz,setprec,setmaxshift);
  orientenergy.vJ = new double*[ofield.numblks];
  orientenergy.hJ = new double*[ofield.numblks];
  for (int n=0;n<ofield.numblks;n++)
  {
    orientenergy.vJ[n] = new double[2*ofield.maxshift+1];
    orientenergy.hJ[n] = new double[2*ofield.maxshift+1];
  }
  return;
}
void estorient::transpose()
{
  const int orienth = (h-1)/ofield.blksz+1;
  const int orientw = (w-1)/ofield.blksz+1;
  swap(orientenergy.hJ,orientenergy.vJ);
  if (orientenergy.hJ != NULL)
  {
    double **edest = new double*[ofield.numblks];
    for (int y=0;y<orienth;y++)
      for (int x=0;x<orientw;x++)
        edest[x*orienth+y]=orientenergy.hJ[y*orientw+x];
    delete[] orientenergy.hJ;
    orientenergy.hJ = edest;
  }
  if (orientenergy.vJ != NULL)
  {
    double **edest = new double*[ofield.numblks];
    for (int y=0;y<orienth;y++)
      for (int x=0;x<orientw;x++)
        edest[x*orienth+y]=orientenergy.vJ[y*orientw+x];
    delete[] orientenergy.vJ;
    orientenergy.vJ = edest;
  }
  dwtnode::transpose();
  return;
}
// Note that for images where the dimensions are not a multiple of blksz, the
// "incomplete" edge blocks will result in proportionately lower energy values
void estorient::calc_energies()
{
	const int orienth = (h-1)/ofield.blksz+1;
  const int orientw = (w-1)/ofield.blksz+1;
  int z, sker; // shift expressed in whole pixels, and fractional (via shift kernel LUT)
  bool ksig;   // flag indicating positive or negative shift
  for (int n=0;n<ofield.numblks;n++)
  {
    int xbase = (n%orientw)*ofield.blksz;
		int ybase = (n/orientw)*ofield.blksz;
    // iterate over all orientations
    for (int sigma=-ofield.maxshift;sigma<=ofield.maxshift;sigma++)
    { // initialise prediction energies to 0
      double accVenergy=0, accHenergy=0;
      double tmpV, tmpH;
      // iterate over each pixel in the block
      for (int y=ybase;(y<ybase+ofield.blksz)&&y<h;y++)
        for (int x=xbase;(x<xbase+ofield.blksz)&&x<w;x++)
        {
          tmpV = tmpH = pixels[y*w+x];
          // calculate the vertical high-pass energy
          kernel_selection(x,sigma,horizontal,z,sker,ksig);
          if (y==0) // top edge must be replicated
            tmpV -= filt(splines+sker,(y+1)*w+x,-z,splines_extent,horizontal,ksig);
          else if (y==h-1) // bottom edge must be replicated
            tmpV -= filt(splines+sker,(y-1)*w+x,z,splines_extent,horizontal,!ksig);
          else
            tmpV -= 0.5*(filt(splines+sker,(y-1)*w+x,z,splines_extent,horizontal,!ksig)
                       +filt(splines+sker,(y+1)*w+x,-z,splines_extent,horizontal,ksig) );
          accVenergy += tmpV*tmpV; // accumulate for all pixels in the block
          // calc horizontal high-pass energy
          kernel_selection(y,sigma,vertical,z,sker,ksig);
          if (x==0) // top edge must be replicated
            tmpH -= filt(splines+sker,y*w+x+1,-z,splines_extent,vertical,ksig);
          else if (x==w-1) // bottom edge must be replicated
            tmpH -= filt(splines+sker,y*w+x-1,z,splines_extent,vertical,!ksig);
          else
            tmpH -= 0.5*(filt(splines+sker,y*w+x-1,z,splines_extent,vertical,!ksig)
                       +filt(splines+sker,y*w+x+1,-z,splines_extent,vertical,ksig) );
          accHenergy += tmpH*tmpH; // accumulate for all pixels in the block
        }
      orientenergy.vJ[n][sigma+ofield.maxshift] = accVenergy;
      orientenergy.hJ[n][sigma+ofield.maxshift] = accHenergy;
    }
  }
  return;
}
// switching penalty is always dependent on the amount of vertical shift required
// (ie the shifting used for implementing the horizontal transform
double estorient::fetch_residual(int n, direction dir, char shift)
{
  const double alpha = 5.0;  // continuity penalty weight
  const double beta  = 50.0; // switching penalty weight
  const double energy =
    (dir==vertical)?orientenergy.vJ[n][ofield.maxshift+shift]
                   :orientenergy.hJ[n][ofield.maxshift+shift];
  const int fieldh = (h-1)/ofield.blksz+1;
  const int fieldw = (w-1)/ofield.blksz+1;
  double cont_penalty   = 0;
  double switch_penalty = 0;
  if (dir==horizontal) // vertical shifts
  {
    // penalise against the horizontal neighbours
    if (n%fieldw != 0) // skip if left edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n-1].vshift);
      if (ofield.orientvec[n-1].hshift != 0)
        switch_penalty += abs(shift); // abs(ofield.orientvec[n-1].hshift);
    }
    if (n%fieldw != (fieldw-1)) // skip if right edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n+1].vshift);
      if (ofield.orientvec[n+1].hshift != 0)
        switch_penalty += abs(shift); // abs(ofield.orientvec[n+1].hshift);
    }
    // penalise against the vertical neighbours
    if (n >= fieldw) // skip if at top edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n-fieldw].vshift);
      if (ofield.orientvec[n-fieldw].hshift != 0)
        switch_penalty += abs(shift); // abs(ofield.orientvec[n-fieldw].hshift);
    }
    if (n < (fieldh-1)*fieldw) // bottom edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n+fieldw].vshift);
      if (ofield.orientvec[n+fieldw].hshift != 0)
        switch_penalty += abs(shift); // abs(ofield.orientvec[n+fieldw].hshift);
    }
  }
  else // (dir==vertical)
  {
    // penalise against the vertical neighbours
    if (n >= fieldw) // skip if at top edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n-fieldw].hshift);
      switch_penalty += abs(ofield.orientvec[n-fieldw].vshift);
    }
    if (n < (fieldh-1)*fieldw) // bottom edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n+fieldw].hshift);
      switch_penalty += abs(ofield.orientvec[n+fieldw].vshift);
    }
    // penalise against the horizontal neighbours
    if (n%fieldw != 0) // left edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n-1].hshift);
      switch_penalty += abs(ofield.orientvec[n-1].vshift);
    }
    if (n%fieldw != (fieldw-1)) // right edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n+1].hshift);
      switch_penalty += abs(ofield.orientvec[n+1].vshift);
    }
  }
  cont_penalty   *= alpha;
  switch_penalty *= beta;
  return energy+cont_penalty+switch_penalty;
}
bool estorient::test_residual(int n, direction dir, char shift, double &best)
{
  const double alpha = 5.0;  // continuity penalty weight
  const double beta  = 50.0; // switching penalty weight
  const double energy =
    (dir==vertical)?orientenergy.vJ[n][ofield.maxshift+shift]
                   :orientenergy.hJ[n][ofield.maxshift+shift];
  if (energy >= best) // early exit
    return false;
  const int fieldh = (h-1)/ofield.blksz+1;
  const int fieldw = (w-1)/ofield.blksz+1;
  double cont_penalty   = 0;
  double switch_penalty = 0;
  if (dir==horizontal) // vertical shifts
  {
    // penalise against the horizontal neighbours
    if (n%fieldw != 0) // left edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n-1].vshift);
      if (ofield.orientvec[n-1].hshift != 0)
        switch_penalty += abs(shift); // abs(ofield.orientvec[n-1].hshift);
    }
    if (n%fieldw != (fieldw-1)) // right edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n+1].vshift);
      if (ofield.orientvec[n+1].hshift != 0)
        switch_penalty += abs(shift); // abs(ofield.orientvec[n+1].hshift);
    }
    // penalise against the vertical neighbours
    if (n >= fieldw) // skip if at top edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n-fieldw].vshift);
      if (ofield.orientvec[n-fieldw].hshift != 0)
        switch_penalty += abs(shift); // abs(ofield.orientvec[n-fieldw].hshift);
    }
    if (n < (fieldh-1)*fieldw) // bottom edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n+fieldw].vshift);
      if (ofield.orientvec[n+fieldw].hshift != 0)
        switch_penalty += abs(shift); // abs(ofield.orientvec[n+fieldw].hshift);
    }
  }
  else // (dir==vertical)
  {
    // penalise against the vertical neighbours
    if (n >= fieldw) // skip if at top edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n-fieldw].hshift);
      switch_penalty += abs(ofield.orientvec[n-fieldw].vshift);
    }
    if (n < (fieldh-1)*fieldw) // bottom edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n+fieldw].hshift);
      switch_penalty += abs(ofield.orientvec[n+fieldw].vshift);
    }
    // penalise against the horizontal neighbours
    if (n%fieldw != 0) // left edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n-1].hshift);
      switch_penalty += abs(ofield.orientvec[n-1].vshift);
    }
    if (n%fieldw != (fieldw-1)) // right edge of field
    {
      cont_penalty += abs(shift - ofield.orientvec[n+1].hshift);
      switch_penalty += abs(ofield.orientvec[n+1].vshift);
    }
  }
  cont_penalty   *= alpha;
  switch_penalty *= beta;
  double total_residual = energy + cont_penalty + switch_penalty;
  bool thresh_met = (total_residual < best);
  if (thresh_met)
    best = total_residual;
  return thresh_met;
}
void estorient::choose_shift(int orientblk, direction dir, double thresh)
{
  if (dir==both)
  {
    cerr << "Direction can only be vertical or horizontal!" << endl;
    exit(1);
  }
  char bestshift = 0;
  double zeroshift = fetch_residual(orientblk,dir,0);
  double best = thresh*zeroshift; // adjust benchmark by threshold
  for (int sigma=-ofield.maxshift;sigma<=ofield.maxshift;sigma++)
  {
    if (test_residual(orientblk,dir,sigma,best))
      bestshift = sigma;
  }
  if (dir==vertical)
    ofield.orientvec[orientblk].hshift = bestshift;
  else // horizontal transform
    ofield.orientvec[orientblk].vshift = bestshift;
  return;
}
void estorient::choose_orient()
{
  double tarray[] = {0.1, 0.4, 0.65, 0.8, 0.95};
  for (int i=0;i<5;i++)
  {
    cout << "pass " << (i+1) << "... ";
    double thresh = tarray[i];
    for (int n=0;n<ofield.numblks;n++)
    { // align the transform with the least 0-shift error
      if (orientenergy.vJ[n][ofield.maxshift] < orientenergy.hJ[n][ofield.maxshift])
      { // align vertical transform
        ofield.orientvec[n].vshift=0;
        choose_shift(n,vertical,thresh);
      }
      else // align horizontal transform
      {
        ofield.orientvec[n].hshift=0;
        choose_shift(n,horizontal,thresh);
      }
    }
  }
  cout << endl;
  return;
}
void estorient::legacy_choose_orient()
{
  for (int n=0;n<ofield.numblks;n++)
  {
    double vJ0 = orientenergy.vJ[n][ofield.maxshift];
    double hJ0 = orientenergy.hJ[n][ofield.maxshift];
    double vEdiff, hEdiff;
    double best = 0;
    direction tmpdim = vertical;
    char tmpshift = 0;
    for (int sigma=-ofield.maxshift;sigma<=ofield.maxshift;sigma++)
    {
      vEdiff = vJ0-orientenergy.vJ[n][sigma+ofield.maxshift];
      hEdiff = hJ0-orientenergy.hJ[n][sigma+ofield.maxshift];
      if (vEdiff > best)
      {
        best = vEdiff;
        tmpdim = vertical;
        tmpshift = sigma;
      }
      if (hEdiff > best)
      {
        best = hEdiff;
        tmpdim = horizontal;
        tmpshift = sigma;
      }
    }
    if (tmpdim == vertical)
    {
      if (best > 0.5*vJ0)
        ofield.orientvec[n].hshift = tmpshift;
      else ofield.orientvec[n].hshift = 0;
    }
    else //(tmpdim == horizontal)
    {
      if (best > 0.5*hJ0)
        ofield.orientvec[n].vshift = tmpshift;
      else ofield.orientvec[n].vshift = 0;
    }
  }
  return;
}