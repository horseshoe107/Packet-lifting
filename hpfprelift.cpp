#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
#define H0_COEFFS {-0.125, 0.25, 0.75, 0.25, -0.125}
#define H0_EXTENT 2
#define G0_COEFFS {0.5, 1, 0.5}
#define G0_EXTENT 1
#define HPF_COEFFS {-2,0,4,3,-5,-19,38,-19,-5,3,4,0,-2} // divide by 64 later
#define HPF_EXTENT 6
static double h0_coeffs[] = H0_COEFFS;
static double g0_coeffs[] = G0_COEFFS;
static double *g0_filter = g0_coeffs + G0_EXTENT;
static double hpf_coeffs[] = HPF_COEFFS;
double * convolve_generate(double *f1, int n1, double *f2, int n2);
#define H0_HPF_EXTENT (H0_EXTENT+HPF_EXTENT)
static double conv_coeffs[2*H0_HPF_EXTENT+1];
static double *h0_hpf_filter =
  convolve_generate(h0_coeffs+H0_EXTENT,H0_EXTENT,hpf_coeffs+HPF_EXTENT,HPF_EXTENT);
double * convolve_generate(double *f1, int n1, double *f2, int n2)
{
  int n3 = n1+n2;
  for (int n=-n2;n<=n2;n++)
    f2[n] /= 64;
  double *f = conv_coeffs+n3;
  for (int n=-n3;n<=n3;n++)
    f[n]=0;
  for (int i=-n1;i<=n1;i++)
    for (int j=-n2;j<=n2;j++)
    {
      f[i+j] += f1[i]*f2[j];
    }
  return f;
}
// note, the LHlift step remains unchanged
void dwtnode::hpf_HLlift(double a, direction dir)
{
  if (k==NULL)
  {
    cerr << "No shift kernel LUT has been selected" << endl;
    exit(1);
  }
  if (dir!=vertical)
  {
    cerr << "experimental function; only vertical allowed for now" << endl;
    exit(2);
  }
  const int s = 1<<dwtlevel[dir]; // stepsize
  const int t = 1<<operating_depth;
  const int last = (dir==vertical)?((h-1)/s)*s : ((w-1)/s)*s;
  const int N = (k->veclen-1)/2; // centre of the shift kernel
  int sigma0, sigma1; // total shift in 1/oprec units
  int z0, z1, sker0, sker1;
  bool ksig0, ksig1;  // flag indicating positive or negative shift
  double *tempbuff = new double[(w>>operating_depth) + 2*H0_HPF_EXTENT];
  if (dir == vertical)
    for (int y=0;y<h;y+=2*s) // lifting to the L rows
    {
      // apply shift first, write output to temporary buffer
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
          tempbuff[(x>>operating_depth)+H0_HPF_EXTENT] = 2*a*
              filt(&k->lut[sker1],(y+s)*w+x,-z1,N,horizontal,ksig1);
        else if (y==last) // replicate bottom edge
          tempbuff[(x>>operating_depth)+H0_HPF_EXTENT] = 2*a*
              filt(&k->lut[sker0],(y-s)*w+x,-z0,N,horizontal,ksig0);
        else
          tempbuff[(x>>operating_depth)+H0_HPF_EXTENT] = a*
            ( filt(&k->lut[sker0],(y-s)*w+x,-z0,N,horizontal,ksig0)
            + filt(&k->lut[sker1],(y+s)*w+x,-z1,N,horizontal,ksig1));
      }
      // constant boundary extension
      for (int x=0;x<H0_HPF_EXTENT;x++)
      {
        tempbuff[x] = tempbuff[H0_HPF_EXTENT];
        tempbuff[(w>>operating_depth)+H0_HPF_EXTENT+x]
          = tempbuff[(w>>operating_depth)+H0_HPF_EXTENT-1];
      }
      for (int xsub=0;xsub<(w>>operating_depth);xsub++)
      {
        double v = 0;
        // apply bandpass + low-pass analysis composite filter
        for (int n=-H0_HPF_EXTENT;n<=H0_HPF_EXTENT;n++)
          v += tempbuff[xsub+H0_HPF_EXTENT+n]*h0_hpf_filter[-n];
        // apply synthesis filter and write output to the low-pass row
        for (int n=-G0_EXTENT;n<=G0_EXTENT;n++)
        {
          int x = (xsub+n) << operating_depth;
          if ((x>=0)&&(x<w))
            pixels[y*w+x] += v*g0_filter[n];
        }
      }
    }
  //else // horizontal
  //  for (int y=0;y<h;y+=t)
  //    for (int x=0;x<w;x+=2*s) // lifting to the L cols
  //    {
  //      sigma0 = ofield.retrieve(y,x-1,horizontal);
  //      sigma1 = ofield.retrieve(y,x,horizontal);
  //      for (int n=1;n<s;n++)
  //      { // accumulate relative shifts of subsequent col pairs
  //        sigma0 += ofield.retrieve(divround(y*ofield.oprec+sigma0,ofield.oprec),x-1-n,horizontal);
  //        sigma1 += ofield.retrieve(divround(y*ofield.oprec-sigma1,ofield.oprec),x+n,horizontal);
  //      }
  //      kernel_selection(y,-sigma0,z0,sker0,ksig0);
  //      kernel_selection(y,sigma1,z1,sker1,ksig1);
  //      if (x==0) // replicate left edge
  //        pixels[y*w+x] += 2*a*
  //            filt(&k->lut[sker1],y*w+x+s,-z1,N,vertical,ksig1);
  //      else if (x==last) // replicate right edge
  //        pixels[y*w+x] += 2*a*
  //            filt(&k->lut[sker0],y*w+x-s,-z0,N,vertical,ksig0);
  //      else
  //        pixels[y*w+x] += a*
  //          ( filt(&k->lut[sker0],y*w+x-s,-z0,N,vertical,ksig0)
  //          + filt(&k->lut[sker1],y*w+x+s,-z1,N,vertical,ksig1));
  //    }
  delete [] tempbuff;
  return;
}
void dwtnode::hpf_oriented_analysis(shker &shftkern, direction dir)
{
  switchkernel(shftkern);
  switch (dwtbase)
  {
  case w5x3:
    hpf_HLlift(-0.5,dir);
    apply_oriented_LHlift(-0.5,dir);
    apply_oriented_HLlift(0.25,dir);
    hpf_HLlift(-0.25,dir);
    apply_gain_factors(1,0.5,dir);
    cout << "stuff is happening" << endl;
    break;
  //case w9x7:
  //  apply_oriented_LHlift(-1.586134342,dir);
  //  apply_oriented_HLlift(-0.052980118,dir);
  //  apply_oriented_LHlift(0.882911075,dir);
  //  apply_oriented_HLlift(0.443506852,dir);
  //  apply_gain_factors(0.812893066,0.615087052,dir);
  //  break;
  default:
    cerr << "No wavelet kernels other than 5x3 permitted" << endl;
    exit(2);
  }
  dwtlevel[dir]++;
  return;
}
void dwtnode::hpf_oriented_synthesis(shker &shftkern, direction dir)
{
  switchkernel(shftkern);
  dwtlevel[dir]--;
  switch (dwtbase)
  {
  case w5x3:
    apply_gain_factors(1,2,dir);
    hpf_HLlift(0.25,dir);
    apply_oriented_HLlift(-0.25,dir);
    apply_oriented_LHlift(0.5,dir);
    hpf_HLlift(0.5,dir);
    break;
  default:
    cerr << "No wavelet kernels other than 5x3 permitted" << endl;
    exit(2);
  }
}