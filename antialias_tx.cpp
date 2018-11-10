#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
int globalcount; // keeps a count of which packet pair is being swapped
void packet_transfer(dwtnode &donor, dwtnode &receiver,
          bool analysis, direction dir)
{
  packlift_filters *f = donor.packf;
  int N = f->htN/2;
  int acc_shft = (dir==horizontal)?-f->offset:-f->offset*receiver.w;
  int N_even = 1-f->htN%2; // 1 if htN is even, 0 otherwise
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
    xstart = N-N_even;
    xend = donor.w-N;
  }
  else // vertical
  {
    ystart = N-N_even;
    yend = donor.h-N;
    xstart = 0;
    xend = donor.w;
  }
  for (int y=ystart;y<yend;y++)
    for (int x=xstart;x<xend;x++)
      receiver.pixels[acc_shft+y*receiver.w+x] += sign*
        donor.filt(&f->ht_coeff[N],y*donor.w+x,0,N,dir,true);
  return;
}
void packet_cancel(dwtnode &donor, dwtnode &receiver,
          bool analysis, direction dir)
{
  packlift_filters *f = donor.packf;
  int N = f->hcN/2;
  int don_shft = (dir==horizontal)?f->offset:f->offset*receiver.w;
  int B = (f->htN+f->hcN-2)/2;
  int sign = (analysis)?+1:-1;
  int ystart, yend, xstart, xend;
  // define boundaries
  if (dir==horizontal)
  {
    ystart = 0;
    yend = donor.h;
    xstart = B-f->offset;
    xend = donor.w-B-f->offset;
  }
  else // vertical
  {
    ystart = B-f->offset;
    yend = donor.h-B-f->offset;
    xstart = 0;
    xend = donor.w;
  }
  for (int y=ystart;y<yend;y++)
    for (int x=xstart;x<xend;x++)
      donor.pixels[don_shft+y*donor.w+x] -= sign*
        receiver.filt(&f->hc_coeff[N],y*receiver.w+x,0,N,dir,true);
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
// computes the local average energy in a packet
double average3x3energy(dwtnode &band, int ycentre, int xcentre)
{
  double total=0;
  int y=ycentre, x=xcentre;
  total = band.pixels[(y-1)*band.w+(x-1)]*band.pixels[(y-1)*band.w+(x-1)];
  total += band.pixels[(y-1)*band.w+(x+1)]*band.pixels[(y-1)*band.w+(x+1)];
  total += band.pixels[(y+1)*band.w+(x-1)]*band.pixels[(y+1)*band.w+(x-1)];
  total += band.pixels[(y+1)*band.w+(x+1)]*band.pixels[(y+1)*band.w+(x+1)];
  total *= 0.5;
  total += band.pixels[(y-1)*band.w+x]*band.pixels[(y-1)*band.w+x];
  total += band.pixels[y*band.w+(x-1)]*band.pixels[y*band.w+(x-1)];
  total += band.pixels[y*band.w+(x+1)]*band.pixels[y*band.w+(x+1)];
  total += band.pixels[(y+1)*band.w+x]*band.pixels[(y+1)*band.w+x];
  total *= 0.5;
  total += band.pixels[y*band.w+x]*band.pixels[y*band.w+x];
  total *= 0.25;
  // 3x3 weighted smoothing filter of the packet energy
  return total;
}
double average3x3abs(dwtnode &band, int ycentre, int xcentre)
{
  double total=0;
  int y=ycentre, x=xcentre;
  total = abs(band.pixels[(y-1)*band.w+(x-1)]);
  total += abs(band.pixels[(y-1)*band.w+(x+1)]);
  total += abs(band.pixels[(y+1)*band.w+(x-1)]);
  total += abs(band.pixels[(y+1)*band.w+(x+1)]);
  total *= 0.5;
  total += abs(band.pixels[(y-1)*band.w+x]);
  total += abs(band.pixels[y*band.w+(x-1)]);
  total += abs(band.pixels[y*band.w+(x+1)]);
  total += abs(band.pixels[(y+1)*band.w+x]);
  total *= 0.5;
  total += abs(band.pixels[y*band.w+x]);
  total *= 0.25;
  // 3x3 weighted smoothing filter of the packet energy
  return total;
}
double average2x3abs(dwtnode &band, int y, int x, direction dir)
{ // shift forwards by 1/2 pixels in the given direction
  double total=0;
  if (dir==vertical)
  {
    total = abs(band.pixels[(y-1)*band.w+(x-1)]);
    total += abs(band.pixels[(y-1)*band.w+(x+1)]);
    total += abs(band.pixels[y*band.w+(x-1)]);
    total += abs(band.pixels[y*band.w+(x+1)]);
    total *= 0.5;
    total += abs(band.pixels[(y-1)*band.w+x]);
    total += abs(band.pixels[y*band.w+x]);
    total *= 0.25;
  }
  else
  {
    total = abs(band.pixels[(y-1)*band.w+(x-1)]);
    total += abs(band.pixels[(y-1)*band.w+x]);
    total += abs(band.pixels[(y+1)*band.w+(x-1)]);
    total += abs(band.pixels[(y+1)*band.w+x]);
    total *= 0.5;
    total += abs(band.pixels[y*band.w+(x-1)]);
    total += abs(band.pixels[y*band.w+x]);
    total *= 0.25;
  }
  return total;
}
// computes the local average energy, using an unweighted NxN window
double averageenergy(dwtnode &packet, int ycentre, int xcentre, int N)
{
  double total=0;
  int y,x;
  if (N%2==0)
  {
    cerr << "Window size must be odd" << endl;
    exit(1);
  }
  int n=N>>1;
  for (int y1=ycentre-n;y1<=ycentre+n;y1++)
    for (int x1=xcentre-n;x1<=xcentre+n;x1++)
    { // boundary replication should never actually be required
      if (y1<0)
        y=0;
      else if (y1>=packet.h)
        y=packet.h-1;
      else
        y=y1;
      if (x1<0)
        x=0;
      else if (x1>=packet.w)
        x=packet.w-1;
      else
        x=x1;
      total += packet.pixels[y*packet.w+x]*packet.pixels[y*packet.w+x];
    }
  return total;
}
double alphlookup(double sideE, double hE)
{
  const double beta=1.0;
  const double T = 1000.0;
  double refE = T + beta*hE;
  if (refE<sideE) // alpha = sideE/hE, but only if hE is sufficiently large
    return 1.0;
  if (0.75*refE<sideE)
    return 0.75;
  if (0.5*refE<sideE)
    return 0.5;
  if (0.25*refE<sideE)
    return 0.25;
  return 0.0;
}
//double alphaTlookup(double hE, double donorE)
//{ // note T scales with the dc gain of the averaging window
//  const double T = 60;
//  // decrease T/beta to lower the threshold of aliasing permitted,
//  // thus increasing alpha on average
//  const double beta = 0.02;
//  double refE = T + beta*donorE;
//  // test against tiered thresholds
//  {
//  if (hE < refE)
//    return 0.0;
//  if (0.765625*hE < refE)
//    return 0.125;
//  if (0.5625*hE < refE)
//    return 0.25;
//  if (0.390625*hE < refE)
//    return 0.375;
//  if (0.25*hE < refE)
//    return 0.5;
//  //if (0.140625*hE < refE)
//  //  return 0.625;
//  if (0.0625*hE < refE)
//    return 0.75;
//  }
//  return 1.0;
//}
double alphaTlookup(dwtnode &hE,dwtnode &donorE,int y,int x,direction dir)
{
  const double T = 1.0;
  const double beta = 0.1;
  double hEp = hE.pixels[y*hE.w+x];
  double refE = beta*donorE.pixels[y*donorE.w+x] + T;
  //if (hE > 1.5*donorE+T) // prevent ringing
  //  return 0;
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
  int step = (dir==vertical)?donorE.w:1;
  for (int n=-4*step;n<5*step;n+=step)
  {
    if (donorE.pixels[y*donorE.w+x+n]<2.0)
    {
      if (alpha>0.5)
        alpha -= 0.25;
      else
        alpha *= 0.5;
    }
  }
  return alpha;
}
double alphaClookup(double hcE, double accE)
{
  double T = 0.0;
  //const double beta = 2.0;
  double refE = accE*accE*accE;
  if (2*hcE>refE)
    return 0;
  //if (4*hcE>refE)
  //  return 0.5;
  //if (hcE+T>0.8*accE)
  //  return 0.3;
  //if (hcE+T>0.6*accE)
  //  return 0.6;
  return 1;
}
void packet_transfer_adaptive(dwtnode &donor, dwtnode &receiver,
          dwtnode &side, bool analysis, direction dir)
{
  packlift_filters *f = donor.packf;
  int N = f->htN/2;
  int acc_shft = (dir==horizontal)?-f->offset:-f->offset*receiver.w;
  int sign = (analysis)?+1:-1;
  int ystart, yend, xstart, xend, y, x;
  double *hee = new double[f->htN+f->hcN-1];
  int heeN = conv(*f,hee);
  fstream fio("sideinf\\alpha_transfer.dat",ios::in|ios::out|ios::binary|ios::app);
  if (!analysis)
  {
    char c;
    for (int n=0;n<globalcount;n++)
      while((c=fio.get())!='\n') ; // skip lines
  }
  // Set boundaries for filtering. Note, for odd order filters
  // the left/top margin is shorter than the right/bottom
  if (dir==horizontal)
  {
    ystart = 0;
    yend = donor.h;
    xstart = N-1+f->htN%2;
    xend = donor.w-N;
  }
  else // vertical
  {
    ystart = N-1+f->htN%2;
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
  for (y=1;y<donor.h-1;y++)
    for (x=1;x<donor.w-1;x++)
    {
      htE.pixels[y*htE.w+x] = average3x3abs(hee_output,y,x);
      donorE.pixels[y*donor.w+x] = average3x3abs(donor,y,x);
    }
  for (y=ystart+1;y<yend-1;y++)
    for (x=xstart+1;x<xend-1;x++)
    {
      double alpha;
      //htE.pixels[y*htE.w+x] = averageenergy(hee_output,y,x,9);
      //donorE.pixels[y*donorE.w+x] = averageenergy(donor,y,x,9);      
      alpha = alphaTlookup(htE,donorE,y,x,dir);
      //if ((globalcount==0)&&(x>121)&&(x<140)&&(y>21)&&(y<42))
      //  alpha=0;
      alpha_map.pixels[y*alpha_map.w+x] = alpha*128;
      htE.pixels[y*htE.w+x] *= 8;
      donorE.pixels[y*donorE.w+x] *= 8;
      //if (analysis)
      //  fio << alpha << " ";
      //else
      //  fio >> alpha;
      receiver.pixels[acc_shft+y*receiver.w+x] += sign*alpha*
        donor.filt(&f->ht_coeff[N],y*donor.w+x,0,N,dir,true);
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
void packet_cancel_adaptive(dwtnode &donor, dwtnode &receiver,
          bool analysis, direction dir)
{
  packlift_filters *f = donor.packf;
  int N = f->hcN/2;
  int don_shft = (dir==horizontal)?f->offset:f->offset*receiver.w;
  int B = (f->htN+f->hcN-2)/2;
  int sign = (analysis)?+1:-1;
  int ystart, yend, xstart, xend;
  //fstream fileio("sideinf\\alpha_cancel.dat",ios::in|ios::out|ios::binary|ios::app);
  //if (!analysis)
  //{
  //  char c;
  //  for (int n=0;n<globalcount;n++)
  //    while((c=fileio.get())!='\n') ; // skip lines
  //}
  // define boundaries
  if (dir==horizontal)
  {
    ystart = 0;
    yend = donor.h;
    xstart = B-f->offset;
    xend = donor.w-B-f->offset;
  }
  else // vertical
  {
    ystart = B-f->offset;
    yend = donor.h-B-f->offset;
    xstart = 0;
    xend = donor.w;
  }
  dwtnode hc_output(receiver.h,receiver.w,disabled,true);
  dwtnode accE(receiver.h,receiver.w,disabled,true);
  dwtnode hcE(receiver.h,receiver.w,disabled,true);
  dwtnode alpha_map(receiver.h,receiver.w,disabled,true);
  for (int y=ystart;y<yend;y++)
    for (int x=xstart;x<xend;x++)
    {
      hc_output.pixels[don_shft+y*receiver.w+x] =
        receiver.filt(&f->hc_coeff[N],y*receiver.w+x,0,N,dir,true);
    }
  for (int y=ystart+1;y<yend-1;y++)
    for (int x=xstart+1;x<xend-1;x++)
    {
      double alpha;
      //hcE = average3x3abs(hc_output,y,x);
      hcE.pixels[y*hcE.w+x] = average2x3abs(hc_output,y,x,dir);
      accE.pixels[y*accE.w+x] = average3x3abs(receiver,y,x);
      //hcE.pixels[y*hcE.w+x] /= accE.pixels[y*accE.w+x];
      //if (globalcount==2)
      alpha = 1;//alphaClookup(hcE.pixels[y*hcE.w+x],accE.pixels[y*accE.w+x]);
      //alpha_map.pixels[y*alpha_map.w+x] = alpha*128;
      //if (analysis)
      //  fileio << alpha << " ";
      //else
      //  fileio >> alpha;
      donor.pixels[don_shft+y*donor.w+x] -= sign*alpha*
        receiver.filt(&f->hc_coeff[N],y*receiver.w+x,0,N,dir,true);
    }
  //fileio << endl;
  //fileio.close();
  if ((globalcount==0)&&(analysis))
  {
    //hcE.pgmwrite("city_hcE.pgm");
    //accE.pgmwrite("city_accE.pgm");
    //alpha_map.pgmwrite("alpha.pgm");
  }
  return;
}
// Carries out 2 lifting steps between a pair of subbands, with filters
// as specified by f.ht_coeff and f.hc_coeff
void packswap(dwtnode &donor, dwtnode &receiver, dwtnode &side, 
          bool analysis, bool adaptive, direction dir)
{
  if (analysis)
  {
    if (adaptive)
      packet_transfer_adaptive(donor, receiver, side, analysis, dir);
    else
      packet_transfer(donor, receiver, analysis, dir);
    packet_cancel(donor, receiver, analysis, dir);
  }
  else
  {
    packet_cancel(donor, receiver, analysis, dir);
    if (adaptive)
      packet_transfer_adaptive(donor, receiver, side, analysis, dir);
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
  if (packf==NULL)
  {
    cerr << "Packet lifting filters must be loaded!" << endl;
    exit(1);
  }
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
    subbands[1]->extract_subband(1);
    subbands[1]->extract_subband(3);
    if (!analysis) // if synthesis, undo horizontal transform first
    {
      packswap(*(subbands[0]->subbands[1]),*(subbands[1]->subbands[0]),
        *(subbands[1]->subbands[1]),analysis,adaptive,horizontal);
      packswap(*(subbands[0]->subbands[3]),*(subbands[1]->subbands[2]),
        *(subbands[1]->subbands[3]),analysis,adaptive,horizontal);
      subbands[1]->interleave();
    }
    else horzsecond = true;
  }
  if ((dim == vertical)||(dim == both))
  {
    if ((subbands[0]->dwtlevel[vertical]!=1)
      ||(subbands[1]->dwtlevel[vertical]!=1))
    {
      cerr << "Warning: vertical dwt levels are inconsistent; packet "
        << "lifting may be invalid" << endl;
    }
    subbands[0]->extract_subband(2); // LL-LH
    subbands[2]->extract_subband(0); // LH-LL
    subbands[2]->extract_subband(1); // LH-HL
    subbands[2]->extract_subband(2);
    subbands[2]->extract_subband(3);
    packswap(*(subbands[0]->subbands[2]),*(subbands[2]->subbands[0]),
      *(subbands[2]->subbands[2]),analysis,adaptive,vertical);
    packswap(*(subbands[0]->subbands[3]),*(subbands[2]->subbands[1]),
      *(subbands[2]->subbands[3]),analysis,adaptive,vertical);
    subbands[2]->interleave();
  }
  if (horzsecond) // if analysis, horizontal transform comes second
  {
    packswap(*(subbands[0]->subbands[1]),*(subbands[1]->subbands[0]),
      *(subbands[1]->subbands[1]),analysis,adaptive,horizontal);
    packswap(*(subbands[0]->subbands[3]),*(subbands[1]->subbands[2]),
      *(subbands[1]->subbands[3]),analysis,adaptive,horizontal);
    subbands[1]->interleave();
  }
  subbands[0]->interleave();
  return;
}