#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
// clears existing orientation data and sets dimensions
// to match the image. call this whenever a new image has
// been loaded, and prior to calling init_orient
void orientationfield::clearfield(int seth, int setw)
{
  h = seth;
  w = setw;
  fieldtype = blockgrid;
  if (orientvec!=NULL)
  {
    delete[] orientvec;
    orientvec = NULL;
  }
  return;
}
int orientationfield::retrieve(int y, int x, direction dir)
{
  if (dir == both)
  {
    cerr << "Direction of the oriented transform can only "
      << "be vertical or horizontal" << endl;
    exit(1);
  }
  if (fieldtype == blockgrid)
    return block_retrieve(y,x,dir);
  else if (fieldtype == affinegrid)
    return affine_retrieve(y,x,dir);
  else // impossible error
  {
    cerr << "Orientation field type cannot be accessed" << endl;
    exit(1);
  }
}
int orientationfield::affine_retrieve(int y, int x, direction dir)
{
  const int o_w = (w-1)/blksz+1; // orientation field width
  const int o_h = (h-1)/blksz+1; // orientation field height
  int xnodeL = (x-blksz/2)/blksz; // location of nearest node to the left (lower node)
  int ynodeL = (y-blksz/2)/blksz;
  bool outofbounds = (xnodeL<0)||(xnodeL>(o_w-2));
  xnodeL = (xnodeL<0)?0:(xnodeL>(o_w-1))?o_w-1:xnodeL;
  int xnodeU = outofbounds?xnodeL:xnodeL+1;
  outofbounds = (ynodeL<0)||(ynodeL>(o_h-2));
  ynodeL = (ynodeL<0)?0:(ynodeL>o_h-1)?o_h-1:ynodeL;
  int ynodeU = outofbounds?ynodeL:ynodeL+1;
  double s1 = (x-xnodeL*blksz-0.5*blksz+0.5)/blksz; // normalised distance from the horizontal L node
  double s2 = (y-ynodeL*blksz-0.5*blksz+0.5)/blksz; // vertical normalised displacement
  int orientLL = (dir==vertical)?
    orientvec[ynodeL*o_w+xnodeL].hshift
    :orientvec[ynodeL*o_w+xnodeL].vshift;
  int orientLU = (dir==vertical)?
    orientvec[ynodeU*o_w+xnodeL].hshift
    :orientvec[ynodeU*o_w+xnodeL].vshift;
  int orientUL = (dir==vertical)?
    orientvec[ynodeL*o_w+xnodeU].hshift
    :orientvec[ynodeL*o_w+xnodeU].vshift;
  int orientUU = (dir==vertical)?
    orientvec[ynodeU*o_w+xnodeU].hshift
    :orientvec[ynodeU*o_w+xnodeU].vshift;
  if ((s1+s2) < 1) // top left triangle affine interpolation
    return (int) round(orientLL + s1*(orientUL-orientLL) + s2*(orientLU-orientLL));
  else // bottom right triangle affine interpolation
    return (int) round(orientUU + (1-s1)*(orientLU-orientUU) + (1-s2)*(orientUL-orientUU));
}
// this function should only be called when dwtnode::transpose
// has occured
void orientationfield::transpose()
{
  if (orientvec != NULL)
  {
    orientation *odest = new orientation[numblks];
    for (int y=0;y<=(h-1)/blksz;y++)
      for (int x=0;x<=(w-1)/blksz;x++)
        odest[x*((h-1)/blksz+1)+y]=orientvec[y*((w-1)/blksz+1)+x];
    for (int n=0;n<numblks;n++)
      swap(odest[n].hshift,odest[n].vshift);
    delete[] orientvec;
    orientvec = odest;
  }
  swap(h,w);
  return;
}
void dwtnode::switchkernel(shker &newk)
{
  k = &newk;
  if (ofield.oprec > k->lutprec)
  {
    cerr << "Interpolation LUT has insufficient precision for"
      << " current orientation field." << endl;
    exit(1);
  }
  if ((k->lutprec%ofield.oprec)!=0)
  {
    cerr << "Interpolation LUT must have a precision that is "
      << "divisible into the orientation field's shift precision."
      << endl;
    exit(2);
  }
  k->skip = k->lutprec/ofield.oprec;
  return;
}
// given the output destination pixel (n) and the desired shift
// (sigma) expressed in 1/oprec units, splits the shift into its
// integral (z) and fractional (LUT_index, shiftdirection) components.
//
// (LUT_index) selects a shift filter from k->lut. As the LUT
// contains filters implementing positive fractional shifts,
// (shiftdirection) indicates whether a negative fractional shift is
// actually desired. Usually, this is effected by reversing the filter
// coefficients, however inband shifting requires more sophisticated
// selection of (LUT_index)
void dwtnode::kernel_selection(int n, int sigma, int &z,
                      int &LUT_index, bool &shiftdirection)
{
  const int N = (k->veclen-1)/2; // centre of the kernel (each shift kernel has length 2N+1)
  const int vecskip = k->skip*k->veclen; // elements to skip per sub-pixel shift increment
  const int phaseskip = (k->lutprec/2+1)*k->veclen; // elements per polyphase set of filters
  const int shiftskip = k->p * phaseskip; // elements per integer shift increment
  z = divround(sigma,ofield.oprec); // find integer component of shift
  // select kernels for the fractional-pixel shifts
  shiftdirection = (sigma >= z*ofield.oprec);
  LUT_index = abs(sigma - z*ofield.oprec)*vecskip + N;
  if (k->inband) // if LUT acts directly on subband coefficients
  { // modify kernel for subband type of the destination pixel
    LUT_index += (n%k->p)*phaseskip;
    // modify kernel for subband pattern of the support
    // if applying a reverse shift on packet subbands, and centre of support
    if ((k->p==4)&&(shiftdirection==false)&&(mod(n+z,2)==1)) // support centre is a H (odd) pixel
      LUT_index += mod(z+2,k->p)*shiftskip;
    else
      LUT_index += mod(z,k->p)*shiftskip;
  }
  return;
}
void dwtnode::apply_oriented_LHlift(double a, direction dir)
{ // error checking
  if (k==NULL)
  {
    cerr << "No shift kernel LUT has been selected" << endl;
    exit(1);
  }
  if (dir==both)
  {
    cerr << "L-H lift; only horizontal or vertical allowed" << endl;
    exit(2);
  }
  const int s = 1<<dwtlevel[dir]; // stepsize
  const int t = 1<<operating_depth;
  const int last = (dir==vertical)?((h-1)/s)*s : ((w-1)/s)*s;
  const int N = (k->veclen-1)/2; // centre of the kernel (each shift kernel has length 2N+1)
  int sigma0, sigma1; // total shift in 1/oprec units
  int z0, z1;         // shift in whole pixels (rounded to nearest integer)
  int sker0, sker1;   // selectors used to access the shift kernel LUT
  bool ksig0, ksig1;  // flag indicating positive or negative shift
  if (dir == vertical)
    for (int y=s;y<h;y+=2*s) // iterating over the destination (H) rows
      for (int x=0;x<w;x+=t)
      { // find relative shift between (y-s,y) and (y,y+s) rows
        sigma0 = ofield.retrieve(y-1,x,vertical);
        sigma1 = ofield.retrieve(y,x,vertical);
        for (int n=1;n<s;n++)
        { // accumulate relative shifts (refer bk3p66)
          sigma0 += ofield.retrieve(y-1-n,x+divround(sigma0,ofield.oprec),vertical);
          sigma1 += ofield.retrieve(y+n,x+divround(-sigma1,ofield.oprec),vertical);
        }
        // determine which LUT filters will implement the shifts
        kernel_selection(x,-sigma0,z0,sker0,ksig0); // refer bk3p65
        kernel_selection(x,sigma1,z1,sker1,ksig1);
        if (y==last) // bottom edge must be replicated
          pixels[y*w+x] += 2*a*
              filt(&k->lut[sker0],(y-s)*w+x,-z0,N,horizontal,ksig0);
        else // NB: filter in the orthogonal direction to the transform
          pixels[y*w+x] += a*
            ( filt(&k->lut[sker0],(y-s)*w+x,-z0,N,horizontal,ksig0)
            + filt(&k->lut[sker1],(y+s)*w+x,-z1,N,horizontal,ksig1));
      }
  else // horizontal
    for (int y=0;y<h;y+=t)
      for (int x=s;x<w;x+=2*s) // iterating over the destination (H) cols
      { // find relative shift between (x-s,x) and (x,x+s) cols
        sigma0 = ofield.retrieve(y,x-1,horizontal);
        sigma1 = ofield.retrieve(y,x,horizontal);
        for (int n=1;n<s;n++)
        { // accumulate relative shifts of subsequent col pairs
          sigma0 += ofield.retrieve(divround(y*ofield.oprec+sigma0,ofield.oprec),x-1-n,horizontal);
          sigma1 += ofield.retrieve(divround(y*ofield.oprec-sigma1,ofield.oprec),x+n,horizontal);
        }
        kernel_selection(y,-sigma0,z0,sker0,ksig0); // choose LUT filters
        kernel_selection(y,sigma1,z1,sker1,ksig1);
        if (x==last) // bottom edge must be replicated
          pixels[y*w+x] += 2*a*
              filt(&k->lut[sker0],y*w+x-s,-z0,N,vertical,ksig0);
        else
          pixels[y*w+x] += a*
            ( filt(&k->lut[sker0],y*w+x-s,-z0,N,vertical,ksig0)
            + filt(&k->lut[sker1],y*w+x+s,-z1,N,vertical,ksig1));
      }
  return;
}
void dwtnode::apply_oriented_HLlift(double a, direction dir)
{
  if (k==NULL)
  {
    cerr << "No shift kernel LUT has been selected" << endl;
    exit(1);
  }
  if (dir==both)
  {
    cerr << "H-L lift; only horizontal or vertical allowed" << endl;
    exit(2);
  }
  const int s = 1<<dwtlevel[dir]; // stepsize
  const int t = 1<<operating_depth;
  const int last = (dir==vertical)?((h-1)/s)*s : ((w-1)/s)*s;
  const int N = (k->veclen-1)/2; // centre of the shift kernel
  int sigma0, sigma1; // total shift in 1/oprec units
  int z0, z1, sker0, sker1;
  bool ksig0, ksig1;  // flag indicating positive or negative shift
  if (dir == vertical)
    for (int y=0;y<h;y+=2*s) // lifting to the L rows
      for (int x=0;x<w;x+=t)
      { // NB: When y==0, the retrieve function might not return a
        // meaningful value. However, in this case sigma0 is not used
        sigma0 = ofield.retrieve(y-1,x,vertical);
        sigma1 = ofield.retrieve(y,x,vertical);
        for (int n=1;n<s;n++)
        { // accumulate relative shifts of subsequent row pairs
          sigma0 += ofield.retrieve(y-1-n,x+divround(sigma0,ofield.oprec),vertical);
          sigma1 += ofield.retrieve(y+n,x+divround(-sigma1,ofield.oprec),vertical);
        }
        kernel_selection(x,-sigma0,z0,sker0,ksig0);
        kernel_selection(x,sigma1,z1,sker1,ksig1);
        if (y==0) // top edge must be replicated
          pixels[y*w+x] += 2*a*
              filt(&k->lut[sker1],(y+s)*w+x,-z1,N,horizontal,ksig1);
        else if (y==last) // replicate bottom edge
          pixels[y*w+x] += 2*a*
              filt(&k->lut[sker0],(y-s)*w+x,-z0,N,horizontal,ksig0);
        else
          pixels[y*w+x] += a*
            ( filt(&k->lut[sker0],(y-s)*w+x,-z0,N,horizontal,ksig0)
            + filt(&k->lut[sker1],(y+s)*w+x,-z1,N,horizontal,ksig1));
      }
  else // horizontal
    for (int y=0;y<h;y+=t)
      for (int x=0;x<w;x+=2*s) // lifting to the L cols
      {
        sigma0 = ofield.retrieve(y,x-1,horizontal);
        sigma1 = ofield.retrieve(y,x,horizontal);
        for (int n=1;n<s;n++)
        { // accumulate relative shifts of subsequent col pairs
          sigma0 += ofield.retrieve(divround(y*ofield.oprec+sigma0,ofield.oprec),x-1-n,horizontal);
          sigma1 += ofield.retrieve(divround(y*ofield.oprec-sigma1,ofield.oprec),x+n,horizontal);
        }
        kernel_selection(y,-sigma0,z0,sker0,ksig0);
        kernel_selection(y,sigma1,z1,sker1,ksig1);
        if (x==0) // replicate left edge
          pixels[y*w+x] += 2*a*
              filt(&k->lut[sker1],y*w+x+s,-z1,N,vertical,ksig1);
        else if (x==last) // replicate right edge
          pixels[y*w+x] += 2*a*
              filt(&k->lut[sker0],y*w+x-s,-z0,N,vertical,ksig0);
        else
          pixels[y*w+x] += a*
            ( filt(&k->lut[sker0],y*w+x-s,-z0,N,vertical,ksig0)
            + filt(&k->lut[sker1],y*w+x+s,-z1,N,vertical,ksig1));
      }
  return;
}
void dwtnode::oriented_analysis(shker &shftkern, direction dir)
{
  if (dir==both)
  {
    cerr << "Cannot analyse in both directions unless two"
      << " shift kernels are supplied. Call oriented_analysis"
      << " (shker &shftbase, shker &shftpoly2) instead" << endl;
    exit(1);
  }
  switchkernel(shftkern);
  switch (dwtbase)
  {
  case w5x3:
    apply_oriented_LHlift(-0.5,dir);
    apply_oriented_HLlift(0.25,dir);
    apply_gain_factors(1,0.5,dir);
    break;
  case w9x7:
    apply_oriented_LHlift(-1.586134342,dir);
    apply_oriented_HLlift(-0.052980118,dir);
    apply_oriented_LHlift(0.882911075,dir);
    apply_oriented_HLlift(0.443506852,dir);
    apply_gain_factors(0.812893066,0.615087052,dir);
    break;
  default:
    cerr << "Unknown wavelet kernel used" << endl;
    exit(2);
  }
  dwtlevel[dir]++;
  return;
}
void dwtnode::oriented_synthesis(shker &shftkern, direction dir)
{
  if (dir==both)
  {
    cerr << "Cannot synthesise in both directions unless two"
      << " shift kernels are supplied. Call oriented_synthesis"
      << " (shker &shftbase, shker &shftpoly2) instead" << endl;
    exit(1);
  }
  switchkernel(shftkern);
  dwtlevel[dir]--;
  switch (dwtbase)
  {
  case w5x3:
    apply_gain_factors(1,2,dir);
    apply_oriented_HLlift(-0.25,dir);
    apply_oriented_LHlift(0.5,dir);
    break;
  case w9x7:
    apply_gain_factors(1.230174105,1.625786132,dir);
    apply_oriented_HLlift(-0.443506852,dir);
    apply_oriented_LHlift(-0.882911075,dir);
    apply_oriented_HLlift(0.052980118,dir);
    apply_oriented_LHlift(1.586134342,dir);
    break;
  default:
    cerr << "Unknown wavelet kernel used" << endl;
    exit(2);
  }
  return;
}
void dwtnode::oriented_analysis(shker &shftbase, shker &shftpoly2)
{
  oriented_analysis(shftbase,vertical);
  oriented_analysis(shftpoly2,horizontal);
  return;
}
void dwtnode::oriented_synthesis(shker &shftbase, shker &shftpoly2)
{
  oriented_synthesis(shftpoly2,horizontal);
  oriented_synthesis(shftbase,vertical);
  return;
}
void dwtnode::oriented_packet_analysis(shker &shftbase, shker &shftpoly4)
{
  oriented_analysis(shftbase,vertical);
  oriented_analysis(shftbase,vertical);
  oriented_analysis(shftpoly4,horizontal);
  oriented_analysis(shftpoly4,horizontal);
  return;
}
void dwtnode::oriented_packet_synthesis(shker &shftbase, shker &shftpoly4)
{
  switchkernel(shftpoly4);
  oriented_synthesis(shftpoly4,horizontal);
  oriented_synthesis(shftpoly4,horizontal);
  oriented_synthesis(shftbase,vertical);
  oriented_synthesis(shftbase,vertical);
  return;
}