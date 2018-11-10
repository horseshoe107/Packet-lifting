#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
#include "orient.h"
extern double splines[];
extern const int splines_extent;
estorient::estorient(dwtnode *target):dwtnode(target->h,target->w,target->txbase)
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
  std::swap(orientenergy.hJ,orientenergy.vJ);
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
          if (x==0) // left edge must be replicated
            tmpH -= filt(splines+sker,y*w+x+1,-z,splines_extent,vertical,ksig);
          else if (x==w-1) // right edge must be replicated
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
    std::cerr << "Direction can only be vertical or horizontal!" << std::endl;
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
    std::cout << "pass " << (i+1) << "... ";
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
  std::cout << std::endl;
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
const float a = 2*log(2.0);
orienttree *root;
class quadtree_est
{
public:
  quadtree_est(){orientenergy[0]=NULL; orientenergy[1]=NULL;}
  ~quadtree_est();
  void init(int seth, int setw, int steminblksize, float setlambda, int maxshift);
  void firstpass_constructtree(int locy, int locx, int depth, orienttree *&node);
  void secondpass(orienttree *node);
  void thirdpass(direction pdir, char pshift, orienttree *node);
  void add_quadtree(int locy, int locx, int depth, orienttree *node, orientation *orientvec);
//private:
  int h, w;
  int size; // next power of 2 that is >= h and w
  int minblksize;
  float lambda;
  int maxshift, s;
  float **orientenergy[2]; // 0=vertical transform
} qtree;
void quadtree_est::init(int seth, int setw, int setminblksize, float setlambda, int setmaxshift)
{
  h = seth;
  w = setw;
  minblksize = setminblksize;
  lambda = setlambda;
  maxshift = setmaxshift;
  s = 2*maxshift+1; // number of shifts allowed per direction
  orientenergy[0] = new float*[h*w];
  orientenergy[1] = new float*[h*w];
  for (int n=0;n<(h*w);n++)
  {
    orientenergy[0][n] = new float[s];
    orientenergy[1][n] = new float[s];
  }
  for (size=1;size<std::max(h,w);size<<=1)
    ;
  return;
}
quadtree_est::~quadtree_est()
{
  for (int i=0;i<2;i++)
    if (orientenergy[i]!=NULL)
    {
      for (int n=0;n<(h*w);n++)
        if (orientenergy[i][n]!=NULL)
          delete[] orientenergy[i][n];
      delete[] orientenergy[i];
    }
}
void estorient2::calc_energies(int minblksize, float lambda)
{
  qtree.init(h,w,minblksize,lambda,this->ofield.maxshift);
  int z, sker; // shift expressed in whole pixels, and fractional (via shift kernel LUT)
  bool ksig;   // flag indicating positive or negative shift
  for (int y=0;y<h;y++)
    for (int x=0;x<w;x++)
    { // iterate over all orientations
      for (int sigma=-ofield.maxshift;sigma<=ofield.maxshift;sigma++)
      {
        double tmpV, tmpH;
        tmpV = tmpH = pixels[y*w+x];
        // calculate the vertical high-pass energy
        kernel_selection(x,sigma,horizontal,z,sker,ksig);
        if (y==0) // top edge must be replicated
          tmpV -= filt(splines+sker,(y+1)*w+x,-z,splines_extent,horizontal,ksig);
        else if (y==(h-1)) // bottom edge must be replicated
          tmpV -= filt(splines+sker,(y-1)*w+x,z,splines_extent,horizontal,!ksig);
        else
          tmpV -= 0.5*(filt(splines+sker,(y-1)*w+x,z,splines_extent,horizontal,!ksig)
                      +filt(splines+sker,(y+1)*w+x,-z,splines_extent,horizontal,ksig) );
        // calc horizontal high-pass energy
        kernel_selection(y,sigma,vertical,z,sker,ksig);
        if (x==0) // left edge must be replicated
          tmpH -= filt(splines+sker,y*w+x+1,-z,splines_extent,vertical,ksig);
        else if (x==(w-1)) // right edge must be replicated
          tmpH -= filt(splines+sker,y*w+x-1,z,splines_extent,vertical,!ksig);
        else
          tmpH -= 0.5*(filt(splines+sker,y*w+x-1,z,splines_extent,vertical,!ksig)
                      +filt(splines+sker,y*w+x+1,-z,splines_extent,vertical,ksig) );
        qtree.orientenergy[0][y*w+x][sigma+ofield.maxshift] = (float) tmpV*tmpV;
        qtree.orientenergy[1][y*w+x][sigma+ofield.maxshift] = (float) tmpH*tmpH;
      }
    }
  return;
}
void quadtree_est::firstpass_constructtree(int locy, int locx, int depth, orienttree *&node)
{ //determine the appropriate area of pixels
  int endy = std::min(h,locy+(size>>depth));
  int endx = std::min(w,locx+(size>>depth));
  int blocksize = std::max((endy-locy)*(endx-locx),0); // blocksize must be non-negative
  float b = lambda*blocksize/a; // knee of log-linear function
  if (node!=NULL)
    delete node;
  node = new orienttree(vertical,0,false);
  node->jprimecand = new jprimevec[2*s];
  for (int shift=0;shift<=2*maxshift;shift++)
  {
    float accE[2]={0.0,0.0};
    for (int y=locy;y<endy;y++)
      for (int x=locx;x<endx;x++)
      {
        accE[0] += orientenergy[0][y*w+x][shift];
        accE[1] += orientenergy[1][y*w+x][shift];
      }
    node->jprimecand[shift].jprime = (accE[0]<=b)? accE[0] : b*(1+log(accE[0]/b));
    node->jprimecand[shift+s].jprime = (accE[1]<=b)? accE[1] : b*(1+log(accE[1]/b));
    //node->jprimecand[shift].jprime = accE[0];
    //node->jprimecand[shift+s].jprime = accE[1];
  }
  if ((size>>depth) > minblksize) // iterate for the children
  {
    int d=depth+1;
    firstpass_constructtree(locy,locx,d,node->children[0]);
    firstpass_constructtree(locy,locx+(size>>d),d,node->children[1]);
    firstpass_constructtree(locy+(size>>d),locx,d,node->children[2]);
    firstpass_constructtree(locy+(size>>d),locx+(size>>d),d,node->children[3]);
  }
  else
    node->leaf = true;
  return;
}
int compute_child_Ls(direction childdir, char childshift, direction parentdir, int parentshift)
{
  // error checking
  if ((childdir==both)||(parentdir==both))
  {
    std::cerr << "direction supplied must be vertical or horizontal only" << std::endl;
    exit(1);
  }
  if (parentshift==0)
    if (childshift==0)
      return 2;
    else
      return abs(childshift)+4;
  if (childdir==parentdir) // code differentially
    childshift -= parentshift;
  if (childshift==0)
    return 3;
  else
    return abs(childshift)+4;
}
void quadtree_est::secondpass(orienttree *node)
{
  if (!(node->leaf))
  { // process children first
    secondpass(node->children[0]);
    secondpass(node->children[1]);
    secondpass(node->children[2]);
    secondpass(node->children[3]);
    for (int shift=0;shift<2*s;shift++) // calculate pruning for every possible shift
    {
      float childsum = 0;
      for (int i=0;i<4;i++)
        childsum += node->children[i]->Jvec[shift].J;
      if (childsum < node->jprimecand[shift].jprime)
      { // combined children's cost is lower so replace
        node->jprimecand[shift].jprime = childsum;
        node->jprimecand[shift].prune  = false;
      }
      else // prune children
        node->jprimecand[shift].prune = true;
    }
  }
  node->Jvec = new bestJvec[2*s];
  for (int pdir=vertical,pindex=0;pdir!=both;pdir++) // iterate over all possible parent shifts
    for (char pshift=-maxshift;pshift<=maxshift;pshift++,pindex++)
    {
      float tmpbestJ = node->jprimecand[maxshift].jprime
        + lambda*compute_child_Ls(vertical,0,(direction)pdir,pshift);
      direction tmpbestdir = vertical;
      char tmpbestshift = 0;
      // for a fixed parent shift, choose child shift that produces the min J
      for (int cdir=vertical,cindex=0; cdir!=both; cdir++)
        for (char cshift=-maxshift;cshift<=maxshift;cshift++,cindex++)
        {
          float currJ = node->jprimecand[cindex].jprime + lambda*
            compute_child_Ls((direction)cdir,cshift,(direction)pdir,pshift);
          if (currJ<tmpbestJ)
          {
            tmpbestJ     = currJ;
            tmpbestdir   = (direction) cdir;
            tmpbestshift = cshift;
          }
        }
      node->Jvec[pindex].J=tmpbestJ;
      node->Jvec[pindex].cdir=tmpbestdir;
      node->Jvec[pindex].cshift=tmpbestshift;
    }
  return;
}
void quadtree_est::thirdpass(direction pdir, char pshift, orienttree *node)
{
  if (pdir==both)
  {
    std::cerr << "Only vertical and horizontal directions allowed" << std::endl;
    exit(1);
  }
  // choose best shift for this node
  int index = pdir*s + pshift + maxshift;
  node->dir   = (node->Jvec[index].cdir == horizontal);
  node->shift = node->Jvec[index].cshift;
  if (!node->leaf) // if the node is already a leaf, finish
  { // choose whether to prune
    if (node->jprimecand[node->dir*s + node->shift + maxshift].prune)
    { // prune
      delete node->children[0];
      delete node->children[1];
      delete node->children[2];
      delete node->children[3];
      node->children[0]=NULL;
      node->children[1]=NULL;
      node->children[2]=NULL;
      node->children[3]=NULL;
      node->leaf = true;
      return;
    }
    else
    { // otherwise recurse
      thirdpass((direction)node->dir, node->shift, node->children[0]);
      thirdpass((direction)node->dir, node->shift, node->children[1]);
      thirdpass((direction)node->dir, node->shift, node->children[2]);
      thirdpass((direction)node->dir, node->shift, node->children[3]);
    }
  }
  return;
}
void estorient2::quadtree_estimate()
{
  qtree.firstpass_constructtree(0,0,0,root); // construct tree, calculate E + J
  std::cout << "cleared first pass" << std::endl;
  qtree.secondpass(root); // calculate best shift+J for every parent possibility
  std::cout << "cleared second pass" << std::endl;
  qtree.thirdpass(vertical, 0, root); // choose shifts and carry out prunings
  std::cout << "cleared third pass" << std::endl;
  codetree("encodedshiftfield.dat",root);
  std::cout << "cleared coding" << std::endl;
}
void quadtree_est::add_quadtree(int locy, int locx, int depth, orienttree *node, orientation *orientvec)
{
  if (node->leaf)
  {
    int endy = std::min(h,locy+(size>>depth));
    int endx = std::min(w,locx+(size>>depth));
    for (int y=locy;y<endy;y++) // locy+1 for more visibility in matlab, locy for the correct field
      for (int x=locx;x<endx;x++)
      { // there should be no overlap, but if there is always overwrite
        if (node->dir) // horizontal
        {
          orientvec[y*w+x].vshift = node->shift;
          orientvec[y*w+x].hshift = 0;
        }
        else
        {
          orientvec[y*w+x].hshift = node->shift;
          orientvec[y*w+x].vshift = 0;
        }
      }
  }
  else // (!node->leaf)
  { // recurse
    int d=depth+1;
    add_quadtree(locy,locx,d,node->children[0],orientvec);
    add_quadtree(locy,locx+(size>>d),d,node->children[1],orientvec);
    add_quadtree(locy+(size>>d),locx,d,node->children[2],orientvec);
    add_quadtree(locy+(size>>d),locx+(size>>d),d,node->children[3],orientvec);
  }
  return;
}
void estorient2::quadtree_flatten()
{
  ofield.init_orient(1,8,8,0,0);
  qtree.add_quadtree(0,0,0,root,ofield.orientvec);
}