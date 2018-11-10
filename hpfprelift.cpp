#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
static class hpfprelift_init
{
public:
  hpfprelift_init();
} _do_it;
static double h0_coeffs[] = {-0.125, 0.25, 0.75, 0.25, -0.125};
#define H0_EXTENT 2
static double g0_coeffs[] = {0.5, 1, 0.5};
#define G0_EXTENT 1
static double g1_coeffs[] = {-0.25, -0.5, 1.5, -0.5, -0.25};
#define G1_EXTENT 2
static double hpf_coeffs[] = {-2,0,4,3,-5,-19,38,-19,-5,3,4,0,-2}; // divide by 64 later
#define HPF_EXTENT 6
static double spline_lut[] = {0, 0, 0, 1, 0, 0, 0,
-0.0057314933, 0.0213902240,-0.0798294029, 0.9674625206, 0.1198302579,-0.0315853847, 0.0084632783,
-0.0088233996, 0.0329293754,-0.1228941022, 0.8796881878, 0.2687588056,-0.0678352739, 0.0181764069,
-0.0097770136, 0.0364883116,-0.1361762329, 0.7514540266, 0.4325511590,-0.1018238761, 0.0272836254,
-0.0090991941, 0.0339586548,-0.1267354249, 0.5973263672, 0.5973263672,-0.1267354249, 0.0339586548};
#define SPLINE_EXTENT 3
#define LUT_SIZE 5
#define INBAND_EXTENT 8 // sets the limit on length of the inband filters
static const int hpf_h0_extent = H0_EXTENT+HPF_EXTENT;
static double *g0_filter = g0_coeffs + G0_EXTENT;
static double *g1_filter = g1_coeffs + G1_EXTENT;
static double *hpf_h0_filter;
static double *hpf_inband_lut;
double *initialise_hpf()
{
  double *filtcentre = hpf_coeffs+HPF_EXTENT;
  for (int n=-HPF_EXTENT;n<=HPF_EXTENT;n++)
    filtcentre[n] /= 64;
  return filtcentre;
}
//interleave using the even-indexed samples from f0 and odd-indexed samples from f1.
//if we are composing an inband filter from its polyphase components, f0 would
//normally be the low-pass response and f1 the high-pass.
void interleave_filters(double *f0, int k0, double *f1, int k1,
                        double *dest, int n, bool f0lowpass=true)
{
  for (int i=-n;i<=n;i++)
    dest[i] = 0;
  int n0 = (k0>=n) ? (n%2==0) ? n:(n-1) // n0 is even
    : (k0%2==0) ? k0:(k0-1);
  int n1 = (k1>=n) ? (n%2==0) ? (n-1):n // n1 is odd
    : (k1%2==0) ? k1-1:k1;
  double sum = 0;
  for (int i=-n0;i<=n0;i+=2)
    sum += dest[i] = f0[i];
  //cout << "Modifying central tap by " << sum << endl;
  dest[0] -= sum; // restore 0 dc
  sum = 0;
  for (int i=-n1;i<=n1;i+=2)
    sum += dest[i] = f1[i];
  //cout << "Modifying -1,1 taps by " << sum << endl;
  dest[-1] -= sum/2; // restore 0 dc
  dest[1]  -= sum/2;
  return;
}
hpfprelift_init::hpfprelift_init()
{
  double *hpf_filter = initialise_hpf();
  hpf_h0_filter = convolve(hpf_filter, HPF_EXTENT, h0_coeffs+H0_EXTENT, H0_EXTENT);
  const int g0_hpf_h0_extent = hpf_h0_extent + G0_EXTENT;
  const int g1_hpf_h0_extent = hpf_h0_extent + G1_EXTENT;
  double *g0_hpf_h0 = convolve(hpf_h0_filter, hpf_h0_extent, g0_filter, G0_EXTENT);
  double *g1_hpf_h0 = convolve(hpf_h0_filter, hpf_h0_extent, g1_filter, G1_EXTENT);
  const int inband_L_extent = g0_hpf_h0_extent + SPLINE_EXTENT;
  const int inband_H_extent = g1_hpf_h0_extent + SPLINE_EXTENT;
  double *inband_L;
  double *inband_H;
  const int veclength = 2*INBAND_EXTENT+1;
  hpf_inband_lut = new double[2*LUT_SIZE*veclength];
  for (int i=0;i<LUT_SIZE;i++)
  {
    inband_L = convolve(g0_hpf_h0,g0_hpf_h0_extent,
      spline_lut+SPLINE_EXTENT+i*(2*SPLINE_EXTENT+1),SPLINE_EXTENT);
    inband_H = convolve(g1_hpf_h0,g1_hpf_h0_extent,
      spline_lut+SPLINE_EXTENT+i*(2*SPLINE_EXTENT+1),SPLINE_EXTENT);
    interleave_filters(inband_L,inband_L_extent,inband_H,inband_H_extent,
      hpf_inband_lut+INBAND_EXTENT+i*veclength,INBAND_EXTENT);
    interleave_filters(inband_H,inband_H_extent,inband_L,inband_L_extent,
      hpf_inband_lut+INBAND_EXTENT+i*veclength+LUT_SIZE*veclength,INBAND_EXTENT);
    delete [] (inband_L - inband_L_extent);
    delete [] (inband_H - inband_H_extent);
  }
  //cleanup
  delete [] (g0_hpf_h0 - g0_hpf_h0_extent);
  delete [] (g1_hpf_h0 - g1_hpf_h0_extent);
}
void hpf_inband_lut_select(int sigma, int oprec, int &z, int &LUT_index, bool &direction)
{
  const int veclength = 2*INBAND_EXTENT+1;
  z = divround(sigma,oprec); // find integer compenent of shift
  direction = (sigma >= z*oprec);
  LUT_index = abs(sigma - z*oprec)*veclength + INBAND_EXTENT;
  if ((z%2)!=0)
    LUT_index += veclength*LUT_SIZE;
  return;
}
// note, the LHlift step remains unchanged
void dwtnode::hpf_HLlift(double a, direction dir)
{
  if (k==NULL)
  {
    cerr << "No shift kernel LUT has been selected" << endl;
    exit(1);
  }
  if (dir==both)
  {
    cerr << "Only vertical or horizontal directions allowed" << endl;
    exit(2);
  }
  const int s = 1<<dwtlevel[dir]; // stepsize
  const int t = 1<<operating_depth;
  const int last = (dir==vertical)?((h-1)/s)*s : ((w-1)/s)*s;
  const int N = (dir==vertical) ? (k->veclen-1)/2 // centre of the shift kernel
    : INBAND_EXTENT;
  int sigma0, sigma1; // total shift in 1/oprec units
  int z0, z1, lut0, lut1;
  bool ksig0, ksig1;  // flag indicating positive or negative shift
  double *tempbuff = new double[(w>>operating_depth) + 2*hpf_h0_extent];
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
        kernel_selection(x,-sigma0,z0,lut0,ksig0);
        kernel_selection(x,sigma1,z1,lut1,ksig1);
        if (y==0) // top edge must be replicated
          tempbuff[(x>>operating_depth)+hpf_h0_extent] = 2*a*
              filt(k->lut+lut1,(y+s)*w+x,-z1,N,horizontal,ksig1);
        else if (y==last) // replicate bottom edge
          tempbuff[(x>>operating_depth)+hpf_h0_extent] = 2*a*
              filt(k->lut+lut0,(y-s)*w+x,-z0,N,horizontal,ksig0);
        else
          tempbuff[(x>>operating_depth)+hpf_h0_extent] = a*
            ( filt(k->lut+lut0,(y-s)*w+x,-z0,N,horizontal,ksig0)
            + filt(k->lut+lut1,(y+s)*w+x,-z1,N,horizontal,ksig1));
      }
      // constant boundary extension
      for (int x=0;x<hpf_h0_extent;x++)
      {
        tempbuff[x] = tempbuff[hpf_h0_extent];
        tempbuff[(w>>operating_depth)+hpf_h0_extent+x]
          = tempbuff[(w>>operating_depth)+hpf_h0_extent-1];
      }
      for (int xsub=0;xsub<(w>>operating_depth);xsub+=2)
      { // find low-pass projection in the baseband by calculating
        // intermediate even (low-pass) samples
        double v = 0;
        // apply bandpass + low-pass analysis composite filter
        for (int n=-hpf_h0_extent;n<=hpf_h0_extent;n++)
          v += tempbuff[xsub+hpf_h0_extent+n]*hpf_h0_filter[-n];
        // apply synthesis filter and write output to the low-pass row
        for (int n=-G0_EXTENT;n<=G0_EXTENT;n++)
        {
          int x = (xsub+n) << operating_depth;
          if ((x>=0)&&(x<w))
            pixels[y*w+x] += v*g0_filter[n];
        }
      }
    }
  else // horizontal
    for (int x=0;x<w;x+=2*s) // lifting to the L cols
      for (int y=0;y<h;y+=2*t) // only lifting to the L (orthogonal) rows as well!
      {
        sigma0 = ofield.retrieve(y,x-1,horizontal);
        sigma1 = ofield.retrieve(y,x,horizontal);
        for (int n=1;n<s;n++)
        { // accumulate relative shifts of subsequent col pairs
          sigma0 += ofield.retrieve(divround(y*ofield.oprec+sigma0,ofield.oprec),x-1-n,horizontal);
          sigma1 += ofield.retrieve(divround(y*ofield.oprec-sigma1,ofield.oprec),x+n,horizontal);
        }
        hpf_inband_lut_select(-sigma0,ofield.oprec,z0,lut0,ksig0);
        hpf_inband_lut_select(sigma1,ofield.oprec,z1,lut1,ksig1);
        if (x==0) // replicate left edge
          pixels[y*w+x] += 2*a*
              filt(hpf_inband_lut+lut1,y*w+x+s,0,N,vertical,ksig1,true);
        else if (x==last) // replicate right edge
          pixels[y*w+x] += 2*a*
              filt(hpf_inband_lut+lut0,y*w+x-s,0,N,vertical,ksig0,true);
        else
          pixels[y*w+x] += a*
            ( filt(hpf_inband_lut+lut0,y*w+x-s,0,N,vertical,ksig0,true)
            + filt(hpf_inband_lut+lut1,y*w+x+s,0,N,vertical,ksig1,true));
      }
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
void dwtnode::hpf_oriented_analysis(shker &shftbase, shker &shftpoly2)
{
  hpf_oriented_analysis(shftbase, vertical);
  hpf_oriented_analysis(shftpoly2, horizontal);
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
void dwtnode::hpf_oriented_synthesis(shker &shftbase, shker &shftpoly2)
{
  hpf_oriented_synthesis(shftpoly2, horizontal);
  hpf_oriented_synthesis(shftbase, vertical);
}