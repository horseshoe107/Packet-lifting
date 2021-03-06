#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
static class packlift_filters
{
public:
  // constructor defined in antialias.cpp
  packlift_filters(char *);
  ~packlift_filters()
  {
    if (ht_coeff != NULL) delete[] ht_coeff;
    if (hc_coeff != NULL) delete[] hc_coeff;
    if (hp_coeff != NULL) delete[] hp_coeff;
  }
  int htN, hcN, hpN, offset;
  double *ht_coeff, *hc_coeff, *hp_coeff;
} f("filters\\icip_aa_filters.dat");
// Initialise antialiasing filters from external file.
// Note that although htN and hcN indicate the true filter length,
// zero padding will occur if necessary to ensure the stored
// coefficients are always odd length support. This allows
// compatibility with the dwtnode::filt() function.
// ht is zero padded on the left side, while hc (and hp, if
// support is added) is padded on the right side. Appropriate
// boundary selection (see packet_transfer and packet_cancel
// functions) must take into account this asymmetry.
packlift_filters::packlift_filters(char *fname)
{
  ht_coeff=NULL, hc_coeff=NULL, hp_coeff=NULL;
  std::ifstream filtin(fname,std::ios::binary);
  if (!filtin.good())
  {
    std::cerr << "Access of file " << fname << " unsuccessful." << std::endl;
    exit(1);
  }
  filtin >> htN;
  int n, h_even = 1-htN%2; // +1 if htN is even, 0 otherwise
  ht_coeff = new double[htN+h_even];
  for (n=0; n<htN; n++)
    filtin >> ht_coeff[n];
  if (htN%2 == 0)  ht_coeff[htN]=0; // pad right side with 0
  filtin >> hcN;
  h_even = 1-hcN%2;
  hc_coeff = new double[hcN+h_even];
  hc_coeff[0] = 0; // pad left side with 0
  for (n=h_even; n<hcN+h_even; n++)
    filtin >> hc_coeff[n];
  if (!filtin.good())
  {
    std::cerr << "Could not read in filters: " << fname
      << " file may be wrong format" << std::endl;
    exit(2);
  }
  if ((htN+hcN)%2 != 0)
  {
    std::cerr << "The transfer and cancellation filters must "
      << "both have even, or odd number of samples - "
      << "end-to-end filter must be linear phase." << std::endl;
    exit(3);
  }
  filtin >> hpN;
  hp_coeff = new double[hpN];
  if (!filtin.good())
    return; // assume hp not specified
  for (n=0;n<hpN;n++)
    filtin >> hp_coeff[n];
  int htmaxindex=0, hcmaxindex=0;
  double maxcoeff=ht_coeff[0];
  for (n=1;n<htN;n++)
    if (maxcoeff<ht_coeff[n])
    {
      maxcoeff = ht_coeff[n];
      htmaxindex = n;
    }
  maxcoeff = hc_coeff[0];
  for (n=1;n<hcN;n++)
    if (maxcoeff<hc_coeff[n])
    {
      maxcoeff = hc_coeff[n];
      hcmaxindex = n;
    }
  offset = htmaxindex - htN/2;
  offset += hcN/2 - hcmaxindex;
  offset >>= 1;
  //offset = -1; // set this to emulate old results
}
void packet_transfer(dwtnode &donor, dwtnode &receiver,
          bool analysis, direction dir)
{
  if (dir==both)
  {
    std::cerr << "packet_transfer: only horizontal or vertical allowed" << std::endl;
    exit(1);
  }
  int N = f.htN/2;
	int Nleft = ((f.htN%2)==0)? N-1:N;
  int acc_shft = (dir==horizontal)?-f.offset:-f.offset*receiver.w;
  int sign = (analysis)?+1:-1;
  int ystart, yend, xstart, xend;
  // Set boundaries for filtering. Note, for odd order filters,
  // the left/top margin is smaller than the right/bottom. This
  // is consistent with the chosen storage of ht_coeff when htN
  // is even, where the right-most element is zero padded. The
  // filter support will overextend to the left/top, but boundary
  // extended samples will multiply by the zero padded filter
  // coefficient (support is the array's reverse)
  if (dir==horizontal)
  {
    ystart = 0;
    yend = donor.h;
    xstart = Nleft;
    xend = donor.w-N;
  }
  else // vertical
  {
    ystart = Nleft;
    yend = donor.h-N;
    xstart = 0;
    xend = donor.w;
  }
  for (int y=ystart;y<yend;y++)
    for (int x=xstart;x<xend;x++)
      receiver.pixels[acc_shft+y*receiver.w+x] += sign*
        donor.filt(&f.ht_coeff[N],y*donor.w+x,0,N,dir,true);
  return;
}
void packet_cancel(dwtnode &donor, dwtnode &receiver,
          bool analysis, direction dir)
{
  if (dir==both)
  {
    std::cerr << "packet_cancel: only horizontal or vertical allowed" << std::endl;
    exit(1);
  }
  int N = f.hcN/2;
  int don_shft = (dir==horizontal)?f.offset:f.offset*receiver.w;
  int B = (f.htN+f.hcN-2)/2;
  int sign = (analysis)?+1:-1;
  int ystart, yend, xstart, xend;
  // define boundaries
  if (dir==horizontal)
  {
    ystart = 0;
    yend = donor.h;
    xstart = B-f.offset;
    xend = donor.w-B-f.offset;
  }
  else // vertical
  {
    ystart = B-f.offset;
    yend = donor.h-B-f.offset;
    xstart = 0;
    xend = donor.w;
  }
  for (int y=ystart;y<yend;y++)
    for (int x=xstart;x<xend;x++)
      donor.pixels[don_shft+y*donor.w+x] -= sign*
        receiver.filt(&f.hc_coeff[N],y*receiver.w+x,0,N,dir,true);
  return;
}
// writes the convolution of filters ht and hc to the array hconv
int conv(packlift_filters &f, double *hconv)
{
  int convN = f.htN + f.hcN - 1;
  //left side of ht is padded when htN is even
  int htpad = 1-f.htN%2;
  for (int n=0;n<convN;n++)
  {
    hconv[n] = 0;
    for (int m=0;m<=n;m++)
    {
      // if beyond the support of ht or hc, the product must be 0
      if (((n-m)>=f.htN)||(m>=f.hcN))
        continue;
      hconv[n] += f.ht_coeff[n-m+htpad]*f.hc_coeff[m];
    }
  }
  return convN/2;
}
double alphaTlookup(double hE,dwtnode &donorE,int y,int x,direction dir)
{
  const double T = 1.0;
  const double beta = 0.1;
	const double Tdonor = 2.0;
  double refE = beta*donorE.pixels[y*donorE.w+x] + T;
  double alpha;
  if (hE > 4*refE)
    alpha = 1.0;
  else if (hE > 2*refE)
    alpha = 0.75;
  else if (hE > 1.33333333*refE)
    alpha = 0.5;
  else if (hE > refE)
    alpha = 0.25;
  else
    alpha = 0;
	// attenuate alpha if surrounding donor energy is low
	// we only want to do packet lifting in large regions of aliasing
	// to avoid creating ringing artifacts
  int step = (dir==vertical)?donorE.w:1;
  for (int n=-4*step;n<5*step;n+=step)
  {
    if (donorE.pixels[y*donorE.w+x+n]<Tdonor)
    {
      if (alpha>0.5)
        alpha -= 0.25;
      else
        alpha *= 0.5;
    }
  }
  return alpha;
}
void packet_transfer_adaptive(dwtnode &donor, dwtnode &receiver,
          bool analysis, direction dir)
{
  if (dir==both)
  {
    std::cerr << "packet_transfer_adaptive: only horizontal or vertical allowed" << std::endl;
    exit(1);
  }
  int N = f.htN/2;                   // For odd order filters the
	int Nleft = ((f.htN%2)==0)? N-1:N; // left margin is shortened
  int acc_shft = (dir==horizontal)?-f.offset:-f.offset*receiver.w;
  int sign = (analysis)?+1:-1;
  int ystart, yend, xstart, xend, y, x;
  double *hee = new double[f.htN+f.hcN-1];
  int heeN = conv(f,hee);
  // Set boundaries for filtering
  if (dir==horizontal)
  {
    ystart = 0;
    yend = donor.h;
    xstart = Nleft;
    xend = donor.w-N;
  }
  else // vertical
  {
    ystart = Nleft;
    yend = donor.h-N;
    xstart = 0;
    xend = donor.w;
  }
  dwtnode hee_output(donor.h,donor.w,none,true);
  dwtnode htE(donor.h,donor.w,none,true);
  dwtnode donorE(donor.h,donor.w,none,true);
  for (y=0;y<donor.h;y++)
    for (x=0;x<donor.w;x++)
    {
      hee_output.pixels[y*donor.w+x] = donor.filt(&hee[heeN],y*donor.w+x,0,heeN,dir,true);
      donorE.pixels[y*donorE.w+x] = donor.filt3x3abs(y,x);
    }
  for (y=ystart+1;y<yend-1;y++)
    for (x=xstart+1;x<xend-1;x++)
    {
      double alpha, htE;
      htE = hee_output.filt3x3abs(y,x);
      alpha = alphaTlookup(htE,donorE,y,x,dir);
      receiver.pixels[acc_shft+y*receiver.w+x] += sign*alpha*
        donor.filt(&f.ht_coeff[N],y*donor.w+x,0,N,dir,true);
    }
  return;
}
// Carries out 2 lifting steps between a pair of subbands, with filters
// as specified by f.ht_coeff and f.hc_coeff
void packswap(dwtnode &donor, dwtnode &receiver,
          bool analysis, bool adaptive, direction dir)
{
  if (analysis)
  {
    if (adaptive)
      packet_transfer_adaptive(donor, receiver, analysis, dir);
    else
      packet_transfer(donor, receiver, analysis, dir);
    packet_cancel(donor, receiver, analysis, dir);
  }
  else
  {
    packet_cancel(donor, receiver, analysis, dir);
    if (adaptive)
      packet_transfer_adaptive(donor, receiver, analysis, dir);
    else
      packet_transfer(donor, receiver, analysis, dir);
  }
  return;
}
// Performs packet lifting between the appropriate packet subbands
// The parent band must be analysed to dwtlevel 1, and the child bands
// (LL1, HL1 and LH1) must be extracted and analysed to dwtlevel 1 in
// the direction of the packet lifting. The dwtlevels of these subbands
// will be restored to this structure before returning.
void dwtnode::packlift(direction dim, bool analysis, bool adaptive)
{
  bool horzsecond = false;
  subbands[0]->extract_subband(3); // the LL-HH packet is always used
  if ((dim == horizontal)||(dim == both))
  {
    if ((subbands[0]->dwtlevel[horizontal]!=1)||(subbands[1]->dwtlevel[horizontal]!=1))
    {
      std::cerr << "Warning: horizontal dwt levels are inconsistent; packet "
        << "lifting may be invalid" << std::endl;
    }
    // extract subband packets needed for horizontal antialiasing
    subbands[0]->extract_subband(1); // LL-HL
    subbands[1]->extract_subband(0); // HL-LL
    subbands[1]->extract_subband(2); // HL-LH
    if (!analysis) // if synthesis, undo horizontal transform first
    {
      packswap(*(subbands[0]->subbands[1]),*(subbands[1]->subbands[0]),
        analysis,adaptive,horizontal);
      packswap(*(subbands[0]->subbands[3]),*(subbands[1]->subbands[2]),
        analysis,adaptive,horizontal);
      subbands[1]->interleave();
    }
    else horzsecond = true;
  }
  if ((dim == vertical)||(dim == both))
  {
    if ((subbands[0]->dwtlevel[vertical]!=1)||(subbands[2]->dwtlevel[vertical]!=1))
    {
      std::cerr << "Warning: vertical dwt levels are inconsistent; packet "
        << "lifting may be invalid" << std::endl;
    }
    subbands[0]->extract_subband(2); // LL-LH
    subbands[2]->extract_subband(0); // LH-LL
    subbands[2]->extract_subband(1); // LH-HL
    packswap(*(subbands[0]->subbands[2]),*(subbands[2]->subbands[0]),
      analysis,adaptive,vertical);
    packswap(*(subbands[0]->subbands[3]),*(subbands[2]->subbands[1]),
      analysis,adaptive,vertical);
    subbands[2]->interleave();
  }
  if (horzsecond) // if analysis, horizontal transform comes second
  {
    packswap(*(subbands[0]->subbands[1]),*(subbands[1]->subbands[0]),
      analysis,adaptive,horizontal);
    packswap(*(subbands[0]->subbands[3]),*(subbands[1]->subbands[2]),
      analysis,adaptive,horizontal);
    subbands[1]->interleave();
  }
  subbands[0]->interleave();
  return;
}
void dwtnode::packlift_analysis(direction dim, bool adapt)
{
  analysis(dim);
  extract_subband(0);
  extract_subband(1);
  extract_subband(2);
  subbands[0]->analysis(both);
  subbands[1]->analysis(both);
  subbands[2]->analysis(both);
  packlift(dim,true,adapt);
  subbands[0]->synthesis(both);
  subbands[1]->synthesis(both);
  subbands[2]->synthesis(both);
  interleave();
  return;
}
void dwtnode::packlift_synthesis(direction dim, bool adapt)
{
  extract_subband(0);
  extract_subband(1);
  extract_subband(2);
  subbands[0]->analysis(both);
  subbands[1]->analysis(both);
  subbands[2]->analysis(both);
  packlift(dim,false,adapt);
  subbands[0]->synthesis(both);
  subbands[1]->synthesis(both);
  subbands[2]->synthesis(both);
  interleave();
  synthesis(both);
  return;
}