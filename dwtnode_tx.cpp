#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
void dwtnode::extract_subband(int band)
{
  if ((dwtlevel[vertical]==0)||(dwtlevel[horizontal]==0))
  {
    cerr << "Subband cannot be extracted - image is still in"
      "baseband form in at least one direction." << endl;
    exit(1);
  }
  // calculate the offsets for the subband when it is in
  // interleaved form
  int xoff = (band%2);
  int yoff = (band/2);
  if (subbands[band]!=NULL)
    delete subbands[band];
  subbands[band] = new dwtnode((h+1-yoff)/2,(w+1-xoff)/2,dwtbase);
  // subband inherits dwtlevel from parent (default is 0)
  if (band==0) // LL
  {
    subbands[band]->dwtlevel[vertical] = this->dwtlevel[vertical]-1;
    subbands[band]->dwtlevel[horizontal] = this->dwtlevel[horizontal]-1;
		subbands[0]->ofield.inherit(this->ofield);
  }
  if (band==1) // HL
    subbands[band]->dwtlevel[vertical] = this->dwtlevel[vertical]-1;
  if (band==2) // LH
    subbands[band]->dwtlevel[horizontal] = this->dwtlevel[horizontal]-1;
  for (int y=0;y<subbands[band]->h;y++)
    for (int x=0;x<subbands[band]->w;x++)
      subbands[band]->pixels[y*subbands[band]->w+x]
        = pixels[(2*y+yoff)*w+2*x+xoff];
  return;
}
void dwtnode::insert_subband(int n, bool suppress_warnings)
{
  if (subbands[n]==NULL)
  {
    cerr << "Attempting to interleave from subband "
      << n << "; data does not exist!" << endl;
    exit(1);
  }
  // check for inconsistency of dwtlevel between parent and child
  if (!suppress_warnings)
  {
    bool warning=false;
    switch(n)
    {
    case 0:
      warning = (this->dwtlevel[vertical]!=subbands[n]->dwtlevel[vertical]+1)
        ||(this->dwtlevel[horizontal]!=subbands[n]->dwtlevel[horizontal]+1);
      break;
    case 1:
      warning = this->dwtlevel[vertical] != subbands[n]->dwtlevel[vertical]+1;
      break;
    case 2:
      warning = this->dwtlevel[horizontal]!=subbands[n]->dwtlevel[horizontal]+1;
      break;
    default: break;
    }
    if (warning)
      cerr << "Warning: subband " << n << "'s dwtlevel is incompatible"
        << " with parent - interleave may cause errors!" << endl;
  }
  int xoff = (n%2);
  int yoff = (n/2);
  for (int y=0;y<subbands[n]->h;y++)
    for (int x=0;x<subbands[n]->w;x++)
      pixels[(2*y+yoff)*w+2*x+xoff]
        = subbands[n]->pixels[y*subbands[n]->w+x];
  delete subbands[n];
  subbands[n] = NULL;
  return;
}
void dwtnode::interleave(bool suppress_warnings)
{
  for (int n=0;n<4;n++)
  {
    if (subbands[n] == NULL)
      continue;
    // check for inconsistency of dwtlevel between parent and child
    if (!suppress_warnings)
    {
      bool warning=false;
      switch(n)
      {
      case 0:
        warning = (this->dwtlevel[vertical]!=subbands[n]->dwtlevel[vertical]+1)
          ||(this->dwtlevel[horizontal]!=subbands[n]->dwtlevel[horizontal]+1);
        break;
      case 1:
        warning = this->dwtlevel[vertical] != subbands[n]->dwtlevel[vertical]+1;
        break;
      case 2:
        warning = this->dwtlevel[horizontal]!=subbands[n]->dwtlevel[horizontal]+1;
        break;
      default: break;
      }
      if (warning)
        cerr << "Warning: subband " << n << "'s dwtlevel is incompatible"
          << " with parent - interleave may cause errors!" << endl;
    }
    int xoff = (n%2);
    int yoff = (n/2);
    for (int y=0;y<subbands[n]->h;y++)
      for (int x=0;x<subbands[n]->w;x++)
        pixels[(2*y+yoff)*w+2*x+xoff]
          = subbands[n]->pixels[y*subbands[n]->w+x];
    delete subbands[n];
    subbands[n] = NULL;
  }
  return;
}
void dwtnode::transpose()
{
  ofield.transpose();
  if (pixels != NULL)
  {
    double *pixdest = new double[h*w];
    for (int y=0;y<h;y++)
      for (int x=0;x<w;x++)
        pixdest[x*h+y] = pixels[y*w+x];
    delete[] pixels;
    pixels = pixdest;
  }
  else
    cerr << "warning: no pixel data to transpose!" << endl;
  swap(h,w);
  for (int n=0;n<4;n++)
    if (subbands[n]!=NULL)
      subbands[n]->transpose();
  swap(subbands[1],subbands[2]);
  return;
}
// Calculates the output of the inner product of filter f and the image over the
// support centred at x,y+offset for vertical filtering or x+offset,y for
// horizontal filtering, where x and y are calculated from n = y*w + x
// Direction of filtering is specified by dir.
// The array holding the filter must have have 2N+1 taps, where f points to the
// centre coefficient. Even length filters can be zero-padded to meet this restriction.
// For normal operation of this function, the output is the inner product of ~f
// (the time-reversal of the filter) with the image, eg:
// out = sum_i sig[n+i]*f[-i]
// If the boolean argument forward is set false then the filter is reversed:
// out = sum_i sig[n+i]*f[i]
// This makes no difference when f is a zero phase filter, but may be useful for
// reversing the direction of shift filters.
// Pixel values beyond the boundary are extrapolated by constant extension. For
// an image that is not in baseband (if dwtlevel[dir] != 0), the last lowpass
// coefficient will be replicated, while the highpass coefficients will be set to 0.
double dwtnode::filt(double *f, int n, int offset, int N, direction dir, bool forward)
{
  if (dir == both)
  {
    cerr << "filt function only operates in one direction" << endl;
    exit(1);
  }
  double *sig = pixels+n;
  int skip = (dir == vertical)?w:1;
  int fsign = forward ? -1 : 1;
  // evaluate the indices for the first and last low-pass pixels in the column/row
  int first_n = (dir==vertical) ? n%w : (n/w)*w;
  int last_n = (dir==vertical)?
    first_n + (((h-1)>>dwtlevel[dir])<<dwtlevel[dir])*w
    : first_n + (((w-1)>>dwtlevel[dir])<<dwtlevel[dir]);
  double Lext = pixels[first_n]; // pixel values required for constant extension
  double Uext = pixels[last_n];
  // calculate how far n is from the edge of the image (0 means n is already on the edge)
  int Lbuf = (dir==vertical)? n/w : n%w; // also the index along the direction of filtering
  int Ubuf = (dir==vertical)? h-1-n/w : w-1-n%w;
  int L=N-offset, U=N+offset; // boundaries of filter support
  int s = 1<<dwtlevel[dir]; // step size when filtering over boundary edges
  int Umod = ((((Lbuf+Ubuf+s)>>dwtlevel[dir])<<dwtlevel[dir]) - Lbuf);
  double out=0;
  if (L > Lbuf) // extension of signal on left/top needed
  {
    for (int i=-Lbuf-s;i>=-L;i-=s)
      out += Lext * f[fsign*(i-offset)];
    L = Lbuf;
  }
  if (U > Ubuf) // extension on right/bot needed
  {
    for (int i=Umod;i<=U;i+=s)
      out += Uext * f[fsign*(i-offset)];
    U = Ubuf;
  }
  for (int i=-L;i<=U;i++)
    out += sig[i*skip]*f[fsign*(i-offset)];
  return out;
}
// Performs 3x3 averaging filter on the magnitude
// The 1-dimension filter coefficients are [1/4 1/2 1/4], and the support
// of the filter is centred on (x,y) and spaced by 2^dwtlevel
// Pixel values beyond the boundary are extrapolated by constant extension
double dwtnode::filt3x3abs(int y, int x)
{
  if ((y<0)||(y>=h)||(x<0)||(x>=w))
  {
    cerr << "filt3x3abs error: y,x outside boundaries of the image" << endl;
    exit(1);
  }
  const int xstep = 1 << dwtlevel[horizontal];
  const int ystep = 1 << dwtlevel[vertical];
  int yminus=((y-ystep)<0)    ?y:y-ystep;
  int yplus =((y+ystep)>(h-1))?y:y+ystep;
  int xminus=((x-xstep)<0)    ?x:x-xstep;
  int xplus =((x+xstep)>(w-1))?x:x+xstep;
  double out = abs(pixels[yminus*w+xminus]);
  out += abs(pixels[yminus*w+xplus]);
  out += abs(pixels[yplus*w+xminus]);
  out += abs(pixels[yplus*w+xplus]);
  out *= 0.5;
  out += abs(pixels[yminus*w+x]);
  out += abs(pixels[y*w+xminus]);
  out += abs(pixels[y*w+xplus]);
  out += abs(pixels[yplus*w+x]);
  out *= 0.5;
  out += abs(pixels[y*w+x]);
  out *= 0.25;
  return out;
}
void dwtnode::apply_LHlift(double a, direction dir)
{
  if (dir == both)
  {
    cerr << "apply_LHlift: only horizontal or vertical allowed" << endl;
    exit(1);
  }
  const int s = 1<<dwtlevel[dir]; // step size
  if (dir == vertical)
  {
    const int last = ((h-1)/s)*s;
    for (int y=s;y<h;y+=2*s)
      for (int x=0;x<w;x++)
        if (y==last) // bottom edge must be replicated
          pixels[y*w+x] += 2*a*pixels[(y-s)*w+x];
        else
          pixels[y*w+x] += a*
            (pixels[(y-s)*w+x] + pixels[(y+s)*w+x]);
  }
  else // horizontal
  {
    const int last = ((w-1)/s)*s;
    for (int y=0;y<h;y++)
      for (int x=s;x<w;x+=2*s)
        if (x==last) // right edge must be replicated
          pixels[y*w+x] += 2*a*pixels[y*w+(x-s)];
        else
          pixels[y*w+x] += a*
            (pixels[y*w+(x-s)] + pixels[y*w+(x+s)]);
  }
  return;
}
void dwtnode::apply_HLlift(double a, direction dir)
{
  if (dir == both)
  {
    cerr << "apply_HLlift: only horizontal or vertical allowed" << endl;
    exit(1);
  }
  const int s = 1<<dwtlevel[dir];
  if (dir == vertical)
  {
    const int last = ((h-1)/s)*s;
    for (int y=0;y<h;y+=2*s)
      for (int x=0;x<w;x++)
        if (y==0) // top edge must be replicated
          pixels[y*w+x] += 2*a*pixels[(y+s)*w+x];
        else if (y==last) // replicate bottom edge
          pixels[y*w+x] += 2*a*pixels[(y-s)*w+x];
        else
          pixels[y*w+x] += a*
            ( pixels[(y-s)*w+x] + pixels[(y+s)*w+x]);
  }
  else // horizontal
  {
    const int last = ((w-1)/s)*s;
    for (int y=0;y<h;y++)
      for (int x=0;x<w;x+=2*s)
        if (x==0) // left edge must be replicated
          pixels[y*w+x] += 2*a*pixels[y*w+(x+s)];
        else if (x==last) // replicate right edge
          pixels[y*w+x] += 2*a*pixels[y*w+(x-s)];
        else
          pixels[y*w+x] += a*
            ( pixels[y*w+(x-s)] + pixels[y*w+(x+s)]);
  }
  return;
}
// scale the low-pass and high-pass cosets by gains k0 and k1 respectively
void dwtnode::apply_gain_factors(double k0, double k1, direction dir)
{
  const int s = 1<<dwtlevel[dir];
  if (dir == vertical)
    for (int y=0;y<h;y+=s)
      if (y%(2*s)==0) //L row
        for (int x=0;x<w;x++)
          pixels[y*w+x] *= k0;
      else //H row
        for (int x=0;x<w;x++)
          pixels[y*w+x] *= k1;
  else // horizontal
    for (int x=0;x<w;x+=s)
      if (x%(2*s)==0) //L col
        for (int y=0;y<h;y++)
          pixels[y*w+x] *= k0;
      else //H col
        for (int y=0;y<h;y++)
          pixels[y*w+x] *= k1;
  return;
}
void dwtnode::analysis(direction dir)
{
  if (dir == both)
  {
    analysis(vertical);
    analysis(horizontal);
    return;
  }
  switch (dwtbase)
  {
  case w5x3:
    apply_LHlift(-0.5,dir);
    apply_HLlift(0.25,dir);
    apply_gain_factors(1,0.5,dir);
    break;
  case w9x7:
    apply_LHlift(-1.586134342,dir);
    apply_HLlift(-0.052980118,dir);
    apply_LHlift(0.882911075,dir);
    apply_HLlift(0.443506852,dir);
    apply_gain_factors(0.812893066,0.615087052,dir);
    break;
  default:
    cerr << "DWT kernel undefined!" << endl;
    exit(1);
    break;
  }
  dwtlevel[dir]++;
  return;
}
void dwtnode::synthesis(direction dir)
{
  if (dir == both)
  {
    synthesis(horizontal);
    synthesis(vertical);
    return;
  }
  dwtlevel[dir]--;
  switch (dwtbase)
  {
  case w5x3:
    apply_gain_factors(1,2,dir);
    apply_HLlift(-0.25,dir);
    apply_LHlift(0.5,dir);
    break;
  case w9x7:
    apply_gain_factors(1.230174105,1.625786132,dir);
    apply_HLlift(-0.443506852,dir);
    apply_LHlift(-0.882911075,dir);
    apply_HLlift(0.052980118,dir);
    apply_LHlift(1.586134342,dir);
    break;
  default:
    cerr << "DWT kernel undefined!" << endl;
    exit(1);
    break;
  }
  return;
}
void dwtnode::packet_analysis(direction dir)
{
  if (dir == both)
  {
    analysis(vertical);
    analysis(vertical);
    analysis(horizontal);
    analysis(horizontal);
  }
  else
  {
    analysis(dir);
    analysis(dir);
  }
  return;
}
void dwtnode::packet_synthesis(direction dir)
{
  if (dir == both)
  {
    synthesis(horizontal);
    synthesis(horizontal);
    synthesis(vertical);
    synthesis(vertical);
  }
  else
  {
    synthesis(dir);
    synthesis(dir);
  }
  return;
}