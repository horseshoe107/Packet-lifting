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
} f("sideinf\\icip_aa_filters.dat");
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
  ifstream filtin(fname,ios::binary);
  if (!filtin.good())
  {
    cerr << "Access of file " << fname << " unsuccessful." << endl;
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
    cerr << "Could not read in filters: " << fname
      << " file may be wrong format" << endl;
    exit(2);
  }
  if ((htN+hcN)%2 != 0)
  {
    cerr << "The transfer and cancellation filters must "
      << "both have even, or odd number of samples - "
      << "end-to-end filter must be linear phase." << endl;
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
int globalcount; // keeps a count of which packet pair is being swapped
void packet_transfer(dwtnode &donor, dwtnode &receiver,
          bool analysis, direction dir)
{
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
void average3x3abs(dwtnode &band, dwtnode &out)
{
	if ((band.h!=out.h)||(band.w!=out.w))
	{
		cerr << "Bands must have identical dimensions" << endl;
	}
	for (int y=1;y<band.h-1;y++)
		for (int x=1;x<band.w-1;x++)
		{
			double tmp = abs(band.pixels[(y-1)*band.w+(x-1)]);
			tmp += abs(band.pixels[(y-1)*band.w+(x+1)]);
			tmp += abs(band.pixels[(y+1)*band.w+(x-1)]);
			tmp += abs(band.pixels[(y+1)*band.w+(x+1)]);
			tmp *= 0.5;
			tmp += abs(band.pixels[(y-1)*band.w+x]);
			tmp += abs(band.pixels[y*band.w+(x-1)]);
			tmp += abs(band.pixels[y*band.w+(x+1)]);
			tmp += abs(band.pixels[(y+1)*band.w+x]);
			tmp *= 0.5;
			tmp += abs(band.pixels[y*band.w+x]);
			tmp *= 0.25;
			out.pixels[y*out.w+x] = tmp;
		}
}
double alphaTlookup(dwtnode &hE,dwtnode &donorE,int y,int x,direction dir)
{
  const double T = 1.0;
  const double beta = 0.1;
	const double Tdonor = 2.0;
  double hEp = hE.pixels[y*hE.w+x]; // aliased energy
  double refE = beta*donorE.pixels[y*donorE.w+x] + T; //
  double alpha;
  if (hEp < refE)
    alpha = 0;
  else if (0.75*hEp < refE)
    alpha = 0.25;
  else if (0.5*hEp < refE)
    alpha = 0.5;
  else if (0.25*hEp < refE)
    alpha = 0.75;
  else
    alpha = 1.0;
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
  int N = f.htN/2;                   // For odd order filters the
	int Nleft = ((f.htN%2)==0)? N-1:N; // left margin is shortened
  int acc_shft = (dir==horizontal)?-f.offset:-f.offset*receiver.w;
  int sign = (analysis)?+1:-1;
  int ystart, yend, xstart, xend, y, x;
  double *hee = new double[f.htN+f.hcN-1];
  int heeN = conv(f,hee);
  fstream fio("sideinf\\alpha_transfer.dat",ios::in|ios::out|ios::binary|ios::app);
  if (!analysis)
  {
    char c;
    for (int n=0;n<globalcount;n++)
      while((c=fio.get())!='\n') ; // skip lines
  }
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
  dwtnode hee_output(donor.h,donor.w,disabled,true);
  dwtnode htE(donor.h,donor.w,disabled,true);
  dwtnode donorE(donor.h,donor.w,disabled,true);
  dwtnode averaged_output(donor.h,donor.w,disabled,true);
  dwtnode alpha_map(donor.h,donor.w,disabled,true);
  for (y=ystart;y<yend;y++)
    for (x=xstart;x<xend;x++)
    {
      hee_output.pixels[y*donor.w+x] = donor.filt(&hee[heeN],y*donor.w+x,0,heeN,dir,true);
    }
	average3x3abs(donor,donorE);
	average3x3abs(hee_output,htE);
  for (y=ystart+1;y<yend-1;y++)
    for (x=xstart+1;x<xend-1;x++)
    {
      double alpha;
      alpha = alphaTlookup(htE,donorE,y,x,dir);
      alpha_map.pixels[y*alpha_map.w+x] = alpha*128;
      htE.pixels[y*htE.w+x] *= 8;
      donorE.pixels[y*donorE.w+x] *= 8;
      //if (analysis)
      //  fio << alpha << " ";
      //else
      //  fio >> alpha;
      receiver.pixels[acc_shft+y*receiver.w+x] += sign*alpha*
        donor.filt(&f.ht_coeff[N],y*donor.w+x,0,N,dir,true);
    }
  fio << endl;
  fio.close();
  //if ((globalcount==0)&&(analysis))
  //{
  //  htE.pgmwrite("tmp\\htE.pgm");
  //  donorE.pgmwrite("tmp\\donorE.pgm");
  //  alpha_map.pgmwrite("tmp\\alphaT.pgm");
  //}
  //if ((globalcount==2)&&(analysis))
  //{
  //  alpha_map.pgmwrite("tmp\\alphaT2.pgm");
  //}
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
  globalcount++;
  if (globalcount==4)
    globalcount=0;
  return;
}
// Perform packet lifting between the appropriate subbands
// It is assumed that the required child bands - eg. LL1, HL1
// and LH1 - have already been extracted, and additionally that
// the child bands have already been analysed into an interleaved
// state. Deeper child subband extraction is performed by the function.
void dwtnode::packlift(direction dim, bool analysis, bool adaptive)
{
  bool horzsecond = false;
  if (analysis)
    globalcount=0;
  else
    globalcount=2;
  subbands[0]->extract_subband(3); // the LL-HH packet is always used
  if ((dim == horizontal)||(dim == both))
  {
    if ((subbands[0]->dwtlevel[horizontal]!=1)
      ||(subbands[1]->dwtlevel[horizontal]!=1))
    {
      cerr << "Warning: horizontal dwt levels are inconsistent; packet "
        << "lifting may be invalid" << endl;
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
    if ((subbands[0]->dwtlevel[vertical]!=1)
      ||(subbands[2]->dwtlevel[vertical]!=1))
    {
      cerr << "Warning: vertical dwt levels are inconsistent; packet "
        << "lifting may be invalid" << endl;
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