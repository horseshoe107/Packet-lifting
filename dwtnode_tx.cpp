#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
void dwtnode::extract_subband(int band)
{
  if ((dwtlevel[vertical]==0)||(dwtlevel[horizontal]==0))
  {
    cerr << "Subband cannot be extracted - image has not been "
      "analysed in at least one direction." << endl;
    exit(1);
  }
  // calculate the offsets for the subband when it is in
  // interleaved form
  int xoff = (band%2);
  int yoff = (band/2);
  subbands[band] = new dwtnode((h+1-yoff)/2,(w+1-xoff)/2,dwtbase);
  subbands[band]->packf = this->packf; // inherit packet lifting filters
  // subband inherits dwtlevel from parent (default is 0)
  if (band==0) // LL
  {
    subbands[band]->dwtlevel[vertical] = this->dwtlevel[vertical]-1;
    subbands[band]->dwtlevel[horizontal] = this->dwtlevel[horizontal]-1;
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
// Calculates the output of a filter with its support centred at pixel location
// n+z (or n+z*w in the case of vertical filtering). z can be set to 0, and
// is simply a convenience for implementing shift filters.
// Filtering is applied vertically or horizontally as denoted by dir. The filter
// must have have 2N+1 taps, where f points at the centre coefficient. For even
// length filters, the array should be zero-padded to meet this restriction.
// Replication of the end pixels is used for boundary extension.
double dwtnode::filt(double *f, int n, int z, int N, direction dir, bool forward)
{
  if (dir == both)
  {
    cerr << "filt function only operates in one direction" << endl;
    exit(1);
  }
  double out=0;
  double *sig = &pixels[n+((dir==vertical)?z*w:z)];
  int skip = (dir == vertical)?w:1;
  // constant boundary extension
  double Lext = (dir==vertical)?pixels[mod(n,w)]:pixels[divfloor(n,w)*w];
  double Uext = (dir==vertical)?pixels[mod(n,w) + (h-1)*w]:pixels[(divfloor(n,w)+1)*w - 1];
  if (false) // set true for zero boundary extension
    Lext = Uext = 0;
  int Lbuf = (dir==vertical)?divfloor(n,w):mod(n,w);
  int Ubuf = (dir==vertical)?h-1 - divfloor(n,w):w-1 - mod(n,w);
  Lbuf += z;
  Ubuf -= z;
  int L=N, U=N;
  if (forward) // normal filtering
  {
    if (N > Lbuf) // extension of signal required
    {
      for (int i=-N;i<-Lbuf;i++)
        out += Lext * f[-i];
      L = Lbuf;
    }
    if (N > Ubuf) // extension on right/bot needed
    {
      for (int i=Ubuf+1;i<=N;i++)
        out += Uext * f[-i];
      U = Ubuf;
    }
    for (int i=-L;i<=U;i++)
      out += sig[i*skip]*f[-i];
  }
  else // reverse filter coefficients (apply negative shift)
  {
    if (N > Lbuf)
    {
      for (int i=-N;i<-Lbuf;i++)
        out += Lext * f[i];
      L = Lbuf;
    }
    if (N > Ubuf)
    {
      for (int i=Ubuf+1;i<=N;i++)
        out += Uext * f[i];
      U = Ubuf;
    }
    for (int i=-L;i<=U;i++)
      out += sig[i*skip]*f[i];
  }
  return out;
}
void dwtnode::apply_LHlift(double a, direction dir)
{
  if (dir == both)
  {
    apply_LHlift(a,vertical);
    apply_LHlift(a,horizontal);
    return;
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
    apply_HLlift(a,vertical);
    apply_HLlift(a,horizontal);
    return;
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
  const int t = 1<<operating_depth;
  if (dir == vertical)
    for (int y=0;y<h;y+=s)
      if (y%(2*s)==0) //L row
        for (int x=0;x<w;x+=t)
          pixels[y*w+x] *= k0;
      else //H row
        for (int x=0;x<w;x+=t)
          pixels[y*w+x] *= k1;
  else // horizontal
    for (int x=0;x<w;x+=s)
      if (x%(2*s)==0) //L col
        for (int y=0;y<h;y+=t)
          pixels[y*w+x] *= k0;
      else //H col
        for (int y=0;y<h;y+=t)
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