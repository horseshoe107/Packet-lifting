#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
extern double splines[];
extern const int splines_extent;
extern const int lutprec;
// must be declared before hpf_prelift_init() runs
static const int lut_size = (lutprec/2+1); // number of vectors per integer-shift block
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
static double hpf_coeffs[] = {-2.0/64,0,4.0/64,3.0/64,-5.0/64,
  -19.0/64,38.0/64,-19.0/64,-5.0/64,3.0/64,4.0/64,0,-2.0/64};
#define HPF_EXTENT 6
static double *g0_filter = g0_coeffs + G0_EXTENT;
static double *g1_filter = g1_coeffs + G1_EXTENT;
static double *hpf_h0_filter;
static const int hpf_h0_extent = H0_EXTENT+HPF_EXTENT;
static double *hpf_inband_lut;
#define INBAND_EXTENT 8 // sets the limit on length of the inband filters
//interleave using the even-indexed samples from f0 and odd-indexed samples from f1.
//f0 and f1 should point to the centre of the filter, where the length of the filters
//are 2*k0+1 and 2*k1+2 respectively.
//the interleaved result is limited to length 2n+1. if k0<n or k1<n then the end
//coefficients will be trimmed, and both polyphase components will be adjusted to
//maintain dc gain of 0
void interleave_filters(double *f0, int k0, double *f1, int k1,
                        double *dest, int n)
{
  int i; double sum;
  for (i=-n;i<=n;i++)
    dest[i] = 0;
  int n0 = (k0>=n) ? (n%2==0) ? n:(n-1) // n0 is even
                   : (k0%2==0) ? k0:(k0-1);
  int n1 = (k1>=n) ? (n%2==0) ? (n-1):n // n1 is odd
                   : (k1%2==0) ? k1-1:k1;
  for (sum=0,i=-n0;i<=n0;i+=2)
    sum += dest[i] = f0[i];
  //cout << "Modifying central tap by " << sum << endl;
  dest[0] -= sum; // restore 0 dc
  for (sum=0,i=-n1;i<=n1;i+=2)
    sum += dest[i] = f1[i];
  //cout << "Modifying -1,1 taps by " << sum << endl;
  dest[-1] -= sum/2; // restore 0 dc
  dest[1]  -= sum/2;
  return;
}
hpfprelift_init::hpfprelift_init()
{
  double *hpf_filter = hpf_coeffs + HPF_EXTENT;
  hpf_h0_filter = convolve(hpf_filter, HPF_EXTENT, h0_coeffs+H0_EXTENT, H0_EXTENT);
	const int g0_hpf_h0_extent = hpf_h0_extent + G0_EXTENT;
  const int g1_hpf_h0_extent = hpf_h0_extent + G1_EXTENT;
  double *g0_hpf_h0 = convolve(hpf_h0_filter, hpf_h0_extent, g0_filter, G0_EXTENT);
  double *g1_hpf_h0 = convolve(hpf_h0_filter, hpf_h0_extent, g1_filter, G1_EXTENT);
  const int inband_L_extent = g0_hpf_h0_extent + splines_extent;
  const int inband_H_extent = g1_hpf_h0_extent + splines_extent;
  double *inband_L;
  double *inband_H;
  const int veclength = 2*INBAND_EXTENT+1;
  hpf_inband_lut = new double[2*lut_size*veclength];
  for (int i=0;i<lut_size;i++)
  {
    inband_L = convolve(g0_hpf_h0,g0_hpf_h0_extent,
      splines+splines_extent+i*(2*splines_extent+1),splines_extent);
    inband_H = convolve(g1_hpf_h0,g1_hpf_h0_extent,
      splines+splines_extent+i*(2*splines_extent+1),splines_extent);
    interleave_filters(inband_L,inband_L_extent,inband_H,inband_H_extent,
      hpf_inband_lut+INBAND_EXTENT+i*veclength,INBAND_EXTENT);
    interleave_filters(inband_H,inband_H_extent,inband_L,inband_L_extent,
      hpf_inband_lut+INBAND_EXTENT+i*veclength+lut_size*veclength,INBAND_EXTENT);
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
    LUT_index += lut_size*veclength;
  return;
}
double prelift_adaptivity_lookup(double sum, double diff)
{
  // calculated a std deviation of 0.4099 in logspace
  //1/exp(0.4099) = 0.6637
  //1/exp(4/3*0.4099) = 0.579
  //1/exp(5/3*0.4099) = 0.505
  //1/exp(2*0.4099) = 0.4405
  if (diff > 0.6637*sum)
    return 0.0;
  if (diff > 0.4405*sum)
    return 0.5;
  //if (diff > 0.579*sum)
  //  return 0.25;
  //if (diff > 0.505*sum)
  //  return 0.5;
  //if (diff > 0.4405*sum)
  //  return 0.75;
  return 1.0;
}
void dwtnode::hpf_HLlift(double a, direction dir, bool adapt)
{
  if (dir==both)
  {
    std::cerr << "hpf_HLlift: only horizontal or vertical allowed" << std::endl;
    exit(1);
  }
  if (dwtlevel[!dir] != ((dir==vertical)?0:1))
  {
    std::cerr << "hpf_HLlift cannot be used with this dwtlevel structure: vertical depth"
		  << dwtlevel[vertical] << ", horizontal depth " << dwtlevel[horizontal] << std::endl;
    exit(2);
  }
  const int s = 1<<dwtlevel[dir]; // stepsize
  const int last = (dir==vertical)?((h-1)/s)*s : ((w-1)/s)*s;
  const int N = (dir==vertical) ? splines_extent // centre of the shift kernel
    : INBAND_EXTENT;
  int sigma0, sigma1; // total shift in 1/oprec units
  int z0, z1, lut0, lut1, n;
  bool ksig0, ksig1;  // flag indicating positive or negative shift
  if (dir == vertical)
	{ // temporary row arrays for holding post-warping data
		dwtnode tmp0(1,w,this->txbase), tmp1(1,w,this->txbase);
		dwtnode sum((h+2*s-1)/(2*s),(w+1)/2,txbase,true), diff=sum;
    for (int y=0;y<h;y+=2*s) // lifting to the L rows
    { // apply warping first, write output to temporary buffers
      for (int x=0;x<w;x++)
      { // NB: When y==0, the retrieve function might not return a
        // meaningful value. However, in this case sigma0 is not used
        for (sigma0=0, sigma1=0, n=0;n<s;n++)
        { // accumulate relative shifts of subsequent row pairs
          sigma0 += ofield.retrieve(y-1-n,x+divround(sigma0,ofield.oprec),vertical);
          sigma1 += ofield.retrieve(y+n,x+divround(-sigma1,ofield.oprec),vertical);
        }
        kernel_selection(x,-sigma0,horizontal,z0,lut0,ksig0);
        kernel_selection(x,sigma1,horizontal,z1,lut1,ksig1);
        if (y==0) // top edge must be replicated
          tmp0.pixels[x] = tmp1.pixels[x] = a*
              filt(splines+lut1,(y+s)*w+x,-z1,N,horizontal,ksig1);
        else if (y==last) // replicate bottom edge
          tmp0.pixels[x] = tmp1.pixels[x] = a*
              filt(splines+lut0,(y-s)*w+x,-z0,N,horizontal,ksig0);
        else
        {
          tmp0.pixels[x] = a*filt(splines+lut0,(y-s)*w+x,-z0,N,horizontal,ksig0);
          tmp1.pixels[x] = a*filt(splines+lut1,(y+s)*w+x,-z1,N,horizontal,ksig1);
        }
      }
      for (int xsub=0;xsub<w;xsub+=2)
      { // find horizontal low-pass aliasing component by applying
        // highpass + low-pass analysis composite filter
        double v0 = tmp0.filt(hpf_h0_filter,xsub,0,hpf_h0_extent,horizontal);
        double v1 = tmp1.filt(hpf_h0_filter,xsub,0,hpf_h0_extent,horizontal);
        sum.pixels[(y/s/2)*sum.w+xsub/2] = v0+v1;
        diff.pixels[(y/s/2)*diff.w+xsub/2] = v0-v1;
      }
		}
		for (int y=0;y<h;y+=2*s)
			for (int xsub=0;xsub<w;xsub+=2)
			{ // if the orthogonal direction is oriented (nonzero shift) then we can assume it
        // doesn't have aliasing, so we skip the hpf_prelift step
        if (ofield.retrieve(y,xsub,!dir)==0)
        { // filter v0 and v1 appropriately
          double sum_n = sum.filt3x3abs(y/s/2,xsub/2);
				  double diff_n = diff.filt3x3abs(y/s/2,xsub/2);
          double b = sum.pixels[y/s/2*sum.w+xsub/2];
          if (adapt)
            b *= prelift_adaptivity_lookup(sum_n,diff_n);
          for (n=-G0_EXTENT;n<=G0_EXTENT;n++) // input-based filtering
          { // apply synthesis filter and write output to the low-pass row
            int x = xsub+n;
            if ((x>=0)&&(x<w))
              pixels[y*w+x] += b*g0_filter[n];
          }
        }
			}
	}
  else // horizontal
  {
    dwtnode sum((h+1)/2,(w+2*s-1)/(2*s),txbase,true), diff=sum;
    for (int x=0;x<w;x+=2*s) // lifting to the L cols
      for (int y=0;y<h;y+=2) // only lifting to the L (orthogonal) rows as well!
      {
        for (sigma0=0, sigma1=0, n=0;n<s;n++)
        { // accumulate relative shifts of subsequent col pairs
          sigma0 += ofield.retrieve(divround(y*ofield.oprec+sigma0,ofield.oprec),x-1-n,horizontal);
          sigma1 += ofield.retrieve(divround(y*ofield.oprec-sigma1,ofield.oprec),x+n,horizontal);
        }
        hpf_inband_lut_select(-sigma0,ofield.oprec,z0,lut0,ksig0);
        hpf_inband_lut_select(sigma1,ofield.oprec,z1,lut1,ksig1);
        double v0, v1;
        if (x==0) // replicate left edge
          v0=v1= a*filt(hpf_inband_lut+lut1,y*w+x+s,0,N,vertical,ksig1);
        else if (x==last) // replicate right edge
          v0=v1= a*filt(hpf_inband_lut+lut0,y*w+x-s,0,N,vertical,ksig0);
        else
        {
          v0 = a*filt(hpf_inband_lut+lut0,y*w+x-s,0,N,vertical,ksig0);
          v1 = a*filt(hpf_inband_lut+lut1,y*w+x+s,0,N,vertical,ksig1);
        }
        sum.pixels[(y/2)*sum.w+x/s/2] = v0+v1;
        diff.pixels[(y/2)*diff.w+x/s/2] = v0-v1;
      }
    for (int x=0;x<w;x+=2*s)
      for (int y=0;y<h;y+=2)
        if (ofield.retrieve(y,x,!dir)==0) // skip step if orthogonal direction is oriented
        { // filter v0 and v1
          double sum_n = sum.filt3x3abs(y/2,x/s/2);
				  double diff_n = diff.filt3x3abs(y/2,x/s/2);
          double b = sum.pixels[y/2*sum.w+x/s/2];
          if (adapt)
            b *= prelift_adaptivity_lookup(sum_n,diff_n);
          pixels[y*w+x] += b;
        }
  }
  return;
}
double update_adaptivity_lookup(double alias, double lowE)
{
  double aliasE = abs(alias);
  //if (aliasE < 0.6752*lowE)
  if (aliasE < 0.5*lowE)
    return 0;
  return 1;
}
void dwtnode::hpf_update_HLlift(double a, direction dir, bool adapt)
{
  if (dir==both)
  {
    std::cerr << "hpf_update_HLlift: only horizontal or vertical allowed" << std::endl;
    exit(1);
  }
  if (dwtlevel[!dir] != ((dir==vertical)?0:1))
  {
    std::cerr << "hpf_HLlift cannot be used with this dwtlevel structure: vertical depth"
		  << dwtlevel[vertical] << ", horizontal depth " << dwtlevel[horizontal] << std::endl;
    exit(2);
  }
  const int s = 1<<dwtlevel[dir]; // stepsize
  const int last = (dir==vertical)?((h-1)/s)*s : ((w-1)/s)*s;
  const int N = (dir==vertical) ? splines_extent // centre of the shift kernel
    : INBAND_EXTENT;
  int sigma0, sigma1; // total shift in 1/oprec units
  int z0, z1, lut0, lut1, n;
  bool ksig0, ksig1;  // flag indicating positive or negative shift
  if (dir == vertical)
	{ // temporary row arrays for holding post-warping data
		dwtnode tmp0(1,w,this->txbase), tmp1(1,w,this->txbase);
		dwtnode v0((h+2*s-1)/(2*s),(w+1)/2,txbase,true), v1=v0;
    dwtnode low0((h+2*s-1)/(2*s),(w+1)/2,txbase,true), low1=low0;
    for (int y=0;y<h;y+=2*s) // lifting to the L rows
    { // apply warping first, write output to temporary buffers
      for (int x=0;x<w;x++)
      { // NB: When y==0, the retrieve function might not return a
        // meaningful value. However, in this case sigma0 is not used
        for (sigma0=0, sigma1=0, n=0;n<s;n++)
        { // accumulate relative shifts of subsequent row pairs
          sigma0 += ofield.retrieve(y-1-n,x+divround(sigma0,ofield.oprec),vertical);
          sigma1 += ofield.retrieve(y+n,x+divround(-sigma1,ofield.oprec),vertical);
        }
        kernel_selection(x,-sigma0,horizontal,z0,lut0,ksig0);
        kernel_selection(x,sigma1,horizontal,z1,lut1,ksig1);
        if (y==0) // top edge must be replicated
          tmp0.pixels[x] = tmp1.pixels[x] = a*
              filt(splines+lut1,(y+s)*w+x,-z1,N,horizontal,ksig1);
        else if (y==last) // replicate bottom edge
          tmp0.pixels[x] = tmp1.pixels[x] = a*
              filt(splines+lut0,(y-s)*w+x,-z0,N,horizontal,ksig0);
        else
        {
          tmp0.pixels[x] = a*filt(splines+lut0,(y-s)*w+x,-z0,N,horizontal,ksig0);
          tmp1.pixels[x] = a*filt(splines+lut1,(y+s)*w+x,-z1,N,horizontal,ksig1);
        }
      }
      for (int xsub=0;xsub<w;xsub+=2)
      { // find horizontal low-pass aliasing component by applying
        // highpass + low-pass analysis composite filter
        v0.pixels[(y/s/2)*v0.w+xsub/2] = tmp0.filt(hpf_h0_filter,xsub,0,hpf_h0_extent,horizontal);
        v1.pixels[(y/s/2)*v1.w+xsub/2] = tmp1.filt(hpf_h0_filter,xsub,0,hpf_h0_extent,horizontal);
        low0.pixels[(y/s/2)*v0.w+xsub/2] = tmp0.filt(h0_coeffs+H0_EXTENT,xsub,0,H0_EXTENT,horizontal);
        low1.pixels[(y/s/2)*v0.w+xsub/2] = tmp1.filt(h0_coeffs+H0_EXTENT,xsub,0,H0_EXTENT,horizontal);
      }
		}
		for (int y=0;y<h;y+=2*s)
			for (int xsub=0;xsub<w;xsub+=2)
			{
        if (ofield.retrieve(y,xsub,!dir)==0) // if orthogonal direction is not oriented
        {
          double b0 = v0.pixels[y/s/2*v0.w+xsub/2];
          double b1 = v1.pixels[y/s/2*v1.w+xsub/2];
          if (adapt)
          {
            double low0_filt = low0.filt3x3abs(y/s/2,xsub/2);
				    double low1_filt = low1.filt3x3abs(y/s/2,xsub/2);
            b0 *= update_adaptivity_lookup(b0, low0_filt);
            b1 *= update_adaptivity_lookup(b1, low1_filt);
          }
          for (n=-G0_EXTENT;n<=G0_EXTENT;n++) // input-based filtering
          { // apply synthesis filter and write output to the low-pass row
            int x = xsub+n;
            if ((x>=0)&&(x<w))
              pixels[y*w+x] += (b0+b1)*g0_filter[n];
          }
        }
			}
	}
  else // horizontal
  { 
    dwtnode v0((h+1)/2,(w+2*s-1)/(2*s),txbase,true), v1=v0;
    dwtnode low0((h+1)/2,(w+2*s-1)/(2*s),txbase,true), low1=low0;
    for (int x=0;x<w;x+=2*s) // lifting to the L cols
      for (int y=0;y<h;y+=2) // only lifting to the L (orthogonal) rows as well!
      {
        for (sigma0=0, sigma1=0, n=0;n<s;n++)
        { // accumulate relative shifts of subsequent col pairs
          sigma0 += ofield.retrieve(divround(y*ofield.oprec+sigma0,ofield.oprec),x-1-n,horizontal);
          sigma1 += ofield.retrieve(divround(y*ofield.oprec-sigma1,ofield.oprec),x+n,horizontal);
        }
        hpf_inband_lut_select(-sigma0,ofield.oprec,z0,lut0,ksig0);
        hpf_inband_lut_select(sigma1,ofield.oprec,z1,lut1,ksig1);
        if (x==0) // replicate left edge
          v0.pixels[(y/2)*v0.w+x/s/2]=v1.pixels[(y/2)*v1.w+x/s/2]
              = a*filt(hpf_inband_lut+lut1,y*w+x+s,0,N,vertical,ksig1);
        else if (x==last) // replicate right edge
          v0.pixels[(y/2)*v0.w+x/s/2]=v1.pixels[(y/2)*v1.w+x/s/2]
              = a*filt(hpf_inband_lut+lut0,y*w+x-s,0,N,vertical,ksig0);
        else
        {
          v0.pixels[(y/2)*v0.w+x/s/2] = a*filt(hpf_inband_lut+lut0,y*w+x-s,0,N,vertical,ksig0);
          v1.pixels[(y/2)*v1.w+x/s/2] = a*filt(hpf_inband_lut+lut1,y*w+x+s,0,N,vertical,ksig1);
        }
      }
    dwtlevel[horizontal]++; // temporarily increase the level for filt3x3abs
    for (int x=0;x<w;x+=2*s)
      for (int y=0;y<h;y+=2)
      {
        if (ofield.retrieve(y,x,!dir)==0) // if orthogonal direction is not oriented
        {
          double b0 = v0.pixels[y/2*v0.w+x/s/2];
          double b1 = v1.pixels[y/2*v1.w+x/s/2];
          if (adapt)
          { // filt3x3abs will only access the odd columns (x=s,3s,5s,7s etc)
            double low0_filt = (x==0)? this->filt3x3abs(y,x+s):this->filt3x3abs(y,x-s);
            double low1_filt = ((x+s)>=w)? low0_filt : this->filt3x3abs(y,x+s);
            b0 *= update_adaptivity_lookup(b0, low0_filt);
            b1 *= update_adaptivity_lookup(b1, low1_filt);
          }
          pixels[y*w+x] += b0+b1;
        }
      }
    dwtlevel[horizontal]--;
  }
  return;
}
void dwtnode::hpf_oriented_analysis(direction dir, bool adapt)
{
	if (dir==both)
	{
		hpf_oriented_analysis(vertical,adapt);
		hpf_oriented_analysis(horizontal,adapt);
		return;
	}
  switch (txbase)
  {
  case w5x3:
    hpf_HLlift(-0.5,dir,adapt);
    apply_oriented_LHlift(-0.5,dir);
    hpf_update_HLlift(-0.25,dir,adapt);
    apply_oriented_HLlift(0.25,dir);
    apply_gain_factors(1,0.5,dir);
    break;
  default:
    std::cerr << "No wavelet kernels other than 5x3 permitted" << std::endl;
    exit(2);
  }
  dwtlevel[dir]++;
  return;
}
void dwtnode::hpf_oriented_synthesis(direction dir, bool adapt)
{
	if (dir==both)
	{
    hpf_oriented_synthesis(horizontal,adapt);
		hpf_oriented_synthesis(vertical,adapt);
		return;
	}
  dwtlevel[dir]--;
  switch (txbase)
  {
  case w5x3:
    apply_gain_factors(1,2,dir);
    apply_oriented_HLlift(-0.25,dir);
    hpf_update_HLlift(0.25,dir,adapt);
    apply_oriented_LHlift(0.5,dir);
    hpf_HLlift(0.5,dir,adapt);
    break;
  default:
    std::cerr << "No wavelet kernels other than 5x3 permitted" << std::endl;
    exit(2);
  }
}