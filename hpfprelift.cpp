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
static double hpf_coeffs[] = {-2,0,4,3,-5,-19,38,-19,-5,3,4,0,-2}; // divide by 64 later
#define HPF_EXTENT 6
static double *g0_filter = g0_coeffs + G0_EXTENT;
static double *g1_filter = g1_coeffs + G1_EXTENT;
static double *hpf_h0_filter;
static const int hpf_h0_extent = H0_EXTENT+HPF_EXTENT;
static double *hpf_inband_lut;
#define INBAND_EXTENT 8 // sets the limit on length of the inband filters
double *initialise_hpf()
{
  double *filtcentre = hpf_coeffs+HPF_EXTENT;
  for (int n=-HPF_EXTENT;n<=HPF_EXTENT;n++)
    filtcentre[n] /= 64;
  return filtcentre;
}
//interleave using the even-indexed samples from f0 and odd-indexed samples from f1.
//f0 and f1 should point to the centre of the filter, where the extent of the filters
//are k0 and k1 respectively.
//the interleaved result has an extent of n. if k0<n or k1<n then the end coefficients
//will be trimmed, and both polyphase components will be adjusted to maintain dc gain
//of 0
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
  double *hpf_filter = initialise_hpf();
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
void dwtnode::hpf_HLlift(double a, direction dir, bool adapt)
{
  if (dir==both)
  {
    cerr << "Only vertical or horizontal directions allowed" << endl;
    exit(1);
  }
  if (dwtlevel[!dir] != ((dir==vertical)?0:1))
  {
    cerr << "hpf_HLlift cannot be used with this dwtlevel structure: vertical depth"
		  << dwtlevel[vertical] << ", horizontal depth " << dwtlevel[horizontal] << endl;
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
    bool *adaptswitch = new bool[h*w];
		dwtnode tmp0(1,w,this->dwtbase), tmp1(1,w,this->dwtbase);
		dwtnode v0((h+1)/2,(w+1)/2,dwtbase,true), v1((h+1)/2,(w+1)/2,dwtbase,true);
    for (int y=0;y<h;y+=2*s) // lifting to the L rows
    {
      // apply shift first, write output to temporary buffer
      for (int x=0;x<w;x++)
      { // NB: When y==0, the retrieve function might not return a
        // meaningful value. However, in this case sigma0 is not used
        for (sigma0=0, sigma1=0, n=0;n<s;n++)
        { // accumulate relative shifts of subsequent row pairs
          sigma0 += ofield.retrieve(y-1-n,x+divround(sigma0,ofield.oprec),vertical);
          sigma1 += ofield.retrieve(y+n,x+divround(-sigma1,ofield.oprec),vertical);
        }
        //if (adapt&&(sigma0==0)&&(sigma1==0)) // adaptivity check
        //  adaptswitch[y*w+x]=false;
        //else adaptswitch[y*w+x]=true;
        kernel_selection(x,-sigma0,dir,z0,lut0,ksig0);
        kernel_selection(x,sigma1,dir,z1,lut1,ksig1);
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
        v0.pixels[(y/2)*v0.w+xsub/2] = tmp0.filt(hpf_h0_filter,xsub,0,hpf_h0_extent,horizontal);
        v1.pixels[(y/2)*v1.w+xsub/2] = tmp1.filt(hpf_h0_filter,xsub,0,hpf_h0_extent,horizontal);
      }
		}
		//v0.rawlwrite("barbv1.rawl");
		//v1.rawlwrite("barbv2.rawl");
		for (int y=0;y<h;y+=2*s)
			for (int xsub=0;xsub<w;xsub+=2)
			{ // filter v0 and v1 appropriately
				double v0_n = v0.pixels[y/2*v0.w+xsub/2]; // or use v0.filt
				double v1_n = v1.pixels[y/2*v1.w+xsub/2];
        for (n=-G0_EXTENT;n<=G0_EXTENT;n++) // input-based filtering
        { // apply synthesis filter and write output to the low-pass row
          int x = xsub+n;
          if ((x>=0)&&(x<w)) //&&adaptswitch[y*w+x] or && another boolean adaptivity function
            pixels[y*w+x] += (v0_n+v1_n)*g0_filter[n]; // * scaling adaptivity function
        }
			}
	}
  else // horizontal
    for (int x=0;x<w;x+=2*s) // lifting to the L cols
      for (int y=0;y<h;y+=2) // only lifting to the L (orthogonal) rows as well!
      {
        for (sigma0=0, sigma1=0, n=0;n<s;n++)
        { // accumulate relative shifts of subsequent col pairs
          sigma0 += ofield.retrieve(divround(y*ofield.oprec+sigma0,ofield.oprec),x-1-n,horizontal);
          sigma1 += ofield.retrieve(divround(y*ofield.oprec-sigma1,ofield.oprec),x+n,horizontal);
        }
        // lazy adaptivity
        //if (adapt&&(sigma0==0)&&(sigma1==0))
        //  continue;
        hpf_inband_lut_select(-sigma0,ofield.oprec,z0,lut0,ksig0);
        hpf_inband_lut_select(sigma1,ofield.oprec,z1,lut1,ksig1);
        if (x==0) // replicate left edge
          pixels[y*w+x] += 2*a*
              filt(hpf_inband_lut+lut1,y*w+x+s,0,N,vertical,ksig1);
        else if (x==last) // replicate right edge
          pixels[y*w+x] += 2*a*
              filt(hpf_inband_lut+lut0,y*w+x-s,0,N,vertical,ksig0);
        else
          pixels[y*w+x] += a*
            ( filt(hpf_inband_lut+lut0,y*w+x-s,0,N,vertical,ksig0)
            + filt(hpf_inband_lut+lut1,y*w+x+s,0,N,vertical,ksig1));
      }
  return;
}
void dwtnode::hpf_oriented_analysis(direction dir, bool adapt)
{
	if (dir==both)
	{
		hpf_oriented_analysis(vertical,adapt);
#ifdef HORIZONTALREADY
		hpf_oriented_analysis(horizontal,adapt);
#else
    oriented_analysis(horizontal);
#endif
		return;
	}
  switch (dwtbase)
  {
  case w5x3:
    hpf_HLlift(-0.5,dir,adapt);
    apply_oriented_LHlift(-0.5,dir);
    apply_oriented_HLlift(0.25,dir);
    //hpf_HLlift(-0.25,dir,adapt);
    apply_gain_factors(1,0.5,dir);
    break;
  default:
    cerr << "No wavelet kernels other than 5x3 permitted" << endl;
    exit(2);
  }
  dwtlevel[dir]++;
  return;
}
void dwtnode::hpf_oriented_synthesis(direction dir, bool adapt)
{
	if (dir==both)
	{
#ifdef HORIZONTALREADY
    hpf_oriented_synthesis(horizontal,adapt);
#else
    oriented_synthesis(horizontal);
#endif
		hpf_oriented_synthesis(vertical,adapt);
		return;
	}
  dwtlevel[dir]--;
  switch (dwtbase)
  {
  case w5x3:
    apply_gain_factors(1,2,dir);
    //hpf_HLlift(0.25,dir,adapt);
    apply_oriented_HLlift(-0.25,dir);
    apply_oriented_LHlift(0.5,dir);
    hpf_HLlift(0.5,dir,adapt);
    break;
  default:
    cerr << "No wavelet kernels other than 5x3 permitted" << endl;
    exit(2);
  }
}