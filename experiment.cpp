#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
extern double splines[];
extern const int splines_extent;
extern double inbandker_lut[];
extern const int inbandker_extent;
extern double packetker_lut[];
extern const int packetker_extent;
double mse(dwtnode &a, dwtnode &b)
{
  if ((a.h>>a.dwtlevel[vertical]!=b.h>>b.dwtlevel[vertical])
    ||(a.w>>a.dwtlevel[horizontal]!=b.w>>b.dwtlevel[horizontal]))
  {
    std::cerr << "Images must have identical h&w dimensions." << std::endl;
    exit(1);
  }
  double apix, bpix, acc=0;
  int adepth=0,bdepth=0;
  dwtnode *aptr, *bptr;
  for (aptr=&a; (adepth<a.dwtlevel[vertical]) && (aptr->subbands[0]!=nullptr); adepth++)
    aptr = aptr->subbands[0];
  for (bptr=&b; (bdepth<b.dwtlevel[vertical]) && (bptr->subbands[0]!=nullptr); bdepth++)
    bptr = bptr->subbands[0];
  int aystep=1<<(a.dwtlevel[vertical]-adepth);
  int axstep=1<<(a.dwtlevel[horizontal]-adepth);
  int bystep=1<<(b.dwtlevel[vertical]-bdepth);
  int bxstep=1<<(b.dwtlevel[horizontal]-bdepth);
  for (int ay=0,by=0;ay<aptr->h;ay+=aystep,by+=bystep)
    for (int ax=0,bx=0;ax<aptr->w;ax+=axstep,bx+=bxstep)
    {
      apix = aptr->pixels[ay*(aptr->w)+ax];
      bpix = bptr->pixels[by*(bptr->w)+bx];
      apix = (apix<0)?0:(apix>255)?255:apix; // ensure data is
      bpix = (bpix<0)?0:(bpix>255)?255:bpix; // [0,255] bounded
      acc += (apix-bpix)*(apix-bpix);
    }
  return acc/(a.h>>a.dwtlevel[vertical])/(a.w>>a.dwtlevel[horizontal]);
}
// applies a constant shift of (sigma/lutprec) pixels in the horizontal
// direction to the entire image
void dwtnode::shift(int sigma, direction dir)
{
  int z, sker;
  bool ksig;
  double *lut = (dwtlevel[dir]==0)?splines:
    (dwtlevel[dir]==1)?inbandker_lut:packetker_lut;
  const int N = (dwtlevel[dir]==0)?splines_extent:
    (dwtlevel[dir]==1)?inbandker_extent:packetker_extent;
  double *pixdest = new double[h*w];
  for (int x=0;x<w;x++)
  {
    kernel_selection(x,sigma,dir,z,sker,ksig);
    for (int y=0;y<h;y++)
    { // if LUT acts on subband coefficients, modify kernel selection
      pixdest[y*w+x] = filt(lut+sker,y*w+x,-z,N,dir,ksig);
    }
  }
  delete[] pixels;
  pixels = pixdest;
  return;
}
// shrink h, w dimensions of image and orientation field
void dwtnode::shrink(int depth)
{
  for (int i=0;i<depth;i++)
  {
    h = (h+1)/2;
    w = (w+1)/2;
  }
  ofield.h=h;
  ofield.w=w;
  if (ofield.blksz%(1<<depth)!=0)
    std::cerr << "Warning: orientation field is invalid" << std::endl;
  ofield.blksz >>= depth;
}
void dwtnode::halveimage()
{
  ofield.h=h=(h+1)/2;
  ofield.w=w=(w+1)/2;
  if (ofield.blksz==1)
  {
    std::cerr << "Warning: orientation field incorrect now" << std::endl;
  }
  else if (ofield.blksz%2!=0) // blksz needs to be multiple of 2
  {
    std::cerr << "Subsampling of shift field is not defined for "
      "block size " << ofield.blksz << std::endl;
    exit(1);
  }
  else
    ofield.blksz >>= 1;
}
// call the appropriate batch file of kdu_compress and kdu_expand commands
// also writes to fout the tested resolution and dwt filters used
void dwtnode::call_batch(testmode mode, char *Cdecomp, bool halfres, std::ofstream &fout)
{
  char *dwtmode_strings[]={" w5x3"," w9x7"," no dwt"};
  std::stringstream batch;
  if (halfres)
  {
    batch << "halfres.bat";
    fout << " half resolution";
    for (;(*Cdecomp!=',')&&(*Cdecomp!='\0');Cdecomp++)
      ; // skip the first level of Cdecomp
    Cdecomp++;
  }
  else if (mode==pyramid_test)
    batch << "pyramid.bat";
  else
    batch << "out.bat";
  fout <<dwtmode_strings[txbase]<<" Cdecomp:"<< Cdecomp << std::endl;
  batch <<" "<<h<<" "<<w<<" "<<txbase<<" \""<<Cdecomp<<"\"";
  system(batch.str().c_str());
}
void dwtnode::call_batch(testmode mode, char *Cdecomp, int depth, int layer, std::ofstream &fout)
{
  char *kdudwt_strings[]={" w5x3"," w9x7"};
  int dwt_mode = (txbase==w5x3)?0:1;
  std::stringstream batch;
  // skip levels from Cdecomp according to layer
  for (int i=0;i<layer;i++)
  { 
    for (;(*Cdecomp!=',')&&(*Cdecomp!='\0');Cdecomp++)
        ;
    if (*Cdecomp==',')
      Cdecomp++;
  }
  int layerh = h>>layer;
  int layerw = w>>layer;
  if (mode==pyramid_test)
  { // need to select batch file depending on number of files to encode
    batch << "pyramid_" << (depth-layer) << ".bat"; // select based on layer and depth
    dwt_mode = 0; // select which wavelet kernel is used by kdu_compress
  }
  else
  {
    batch << "out.bat";
  }
  fout << "depth " << depth << " layer " << layer << kdudwt_strings[dwt_mode]
    << " Cdecomp=" << Cdecomp << std::endl;
  batch <<" "<<layerh<<" "<<layerw<<" "<<dwt_mode<<" \""<<Cdecomp<<"\"";
  system(batch.str().c_str());
}
// wrappers for analysis-encode-decode-synthesis experimentation
void dwtnode::rawl_encode(bool halfres, bool adapt)
{
  if (halfres)
  {
    analysis(both);
    rawlwrite("tmp\\out.rawl");
  }
  else rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::rawl_encode(int depth, int layer, bool adapt)
{
  for (int i=0;i<layer;i++)
    analysis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::rawl_decode(char *bitrate, bool halfres, bool adapt)
{
  std::string fname = "tmp\\out";
  fname += bitrate;
  fname += ".rawl";
  rawlread((char *)fname.c_str());
}
void dwtnode::rawl_decode(char *bitrate, int depth, int layer, bool adapt)
{
  std::string fname = "tmp\\out";
  fname += bitrate;
  fname += ".rawl";
  rawlread((char *)fname.c_str());
  dwtlevel[vertical]=dwtlevel[horizontal]=depth-layer;
  for (int i=depth-1;i>=layer;i--)
    synthesis(both);
}
void dwtnode::packlift_encode(bool halfres, bool adapt)
{
  analysis(both);
  extract_subband(0);
  extract_subband(1);
  extract_subband(2);
  subbands[0]->analysis(both);
  subbands[1]->analysis(both);
  subbands[2]->analysis(both);
  packlift(both,true,adapt);
  subbands[0]->synthesis(both);
  if (halfres)
  {
    subbands[0]->rawlwrite("tmp\\out.rawl");
    return;
  }
  subbands[1]->synthesis(both);
  subbands[2]->synthesis(both);
  interleave();
  synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::packlift_encode(int depth, int layer, bool adapt)
{
  //packlift_analysis(both,adapt);
  //if (!halfres)
  //  synthesis(both);
  //rawlwrite("tmp\\out.rawl");
  //return;
}
void dwtnode::packlift_decode(char *bitrate, bool halfres, bool adapt)
{
  rawl_decode(bitrate,halfres,adapt);
  if (halfres) return;
  analysis(both);
  extract_subband(0);
  extract_subband(1);
  extract_subband(2);
  subbands[0]->analysis(both);
  subbands[1]->analysis(both);
  subbands[2]->analysis(both);
  packlift(both,false,adapt);
  subbands[0]->synthesis(both);
  subbands[1]->synthesis(both);
  subbands[2]->synthesis(both);
  interleave();
  synthesis(both);
  return;
}
void dwtnode::packlift_decode(char *bitrate, int depth, int layer, bool adapt)
{
 // rawl_decode(bitrate,halfres,adapt);
 // if (!halfres)
 // {
 //   analysis(both);
 //   packlift_synthesis(both,adapt);
 // }
	//return;
}
void dwtnode::packlift_2layer_encode(bool quartres, bool adapt)
{
  packlift_analysis(both,adapt);
  extract_subband(0);
  subbands[0]->packlift_analysis(both,adapt);
  subbands[0]->synthesis(both);
  if (quartres)
  {
    subbands[0]->rawlwrite("tmp\\out.rawl");
    return;
  }
  interleave();
  synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::packlift_2layer_decode(char *bitrate, bool quartres, bool adapt)
{
  rawl_decode(bitrate,quartres,adapt);
  if (quartres) return;
  analysis(both);
  extract_subband(0);
  subbands[0]->analysis(both);
  subbands[0]->packlift_synthesis(both,adapt);
  interleave();
  packlift_synthesis(both,adapt);
  return;
}
void dwtnode::orient_encode(bool halfres, bool adapt)
{
  oriented_analysis(both);
  if (!halfres)
    synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::orient_encode(int depth, int layer, bool adapt)
{
}
void dwtnode::orient_decode(char *bitrate, bool halfres, bool adapt)
{
  rawl_decode(bitrate,halfres,adapt);
  if (!halfres)
  {
    analysis(both);
    oriented_synthesis(both);
  }
  return;
}
void dwtnode::orient_decode(char *bitrate, int depth, int layer, bool adapt)
{
}
void dwtnode::orient2_packet_encode(bool halfres, bool adapt)
{
  oriented_packet_analysis(both);
  if (halfres)
  {
		extract_subband(0);
    subbands[0]->synthesis(both); // repacked for kdu encoding
    subbands[0]->rawlwrite("tmp\\out.rawl");
		return;
  }
  packet_synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::orient2_packet_decode(char *bitrate, bool halfres, bool adapt)
{
  rawl_decode(bitrate,halfres,adapt);
  if (halfres)
  { // unpack then resynthesise
    analysis(both);
    oriented_synthesis(both);
    return;
  }
  packet_analysis(both);
  oriented_packet_synthesis(both);
  return;
}
void dwtnode::orient_2layer_encode(bool halfres, bool adapt)
{
  oriented_analysis(both);
  extract_subband(0); // ofield automatically inherited
  subbands[0]->oriented_analysis(both); // excluding detail subbands reduces aliasing
  subbands[0]->synthesis(both);
  if (halfres)
  {
    subbands[0]->rawlwrite("tmp\\out.rawl");
    return;
  }
  interleave();
  synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::orient_2layer_decode(char *bitrate, bool halfres, bool adapt)
{
  rawl_decode(bitrate,halfres,adapt);
  analysis(both);
  if (halfres)
  {
    oriented_synthesis(both);
    return;
  }
  extract_subband(0);
  subbands[0]->analysis(both);
  subbands[0]->oriented_synthesis(both);
  interleave();
  oriented_synthesis(both);
  return;
}
void dwtnode::packlift_orient2_encode(bool halfres, bool adapt)
{
  oriented_analysis(both);
  extract_subband(0);
  extract_subband(1);
  extract_subband(2);
  subbands[0]->analysis(both);
  subbands[1]->analysis(both);
  subbands[2]->analysis(both);
  packlift(both,true,adapt);
  subbands[0]->synthesis(both);
  subbands[0]->oriented_analysis(both);
  subbands[0]->synthesis(both);
  if (halfres)
  {
    subbands[0]->rawlwrite("tmp\\out.rawl");
    return;
  }
  subbands[1]->synthesis(both);
  subbands[2]->synthesis(both);
  interleave();
  synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::packlift_orient2_decode(char *bitrate, bool halfres, bool adapt)
{
  rawl_decode(bitrate,halfres,adapt);
  analysis(both);
  if (halfres)
  {
    oriented_synthesis(both);
    return;
  }
  extract_subband(0);
  extract_subband(1);
  extract_subband(2);
  subbands[0]->analysis(both);
  subbands[0]->oriented_synthesis(both);
  subbands[0]->analysis(both);
  subbands[1]->analysis(both);
  subbands[2]->analysis(both);
  packlift(both,false,adapt);
  subbands[0]->synthesis(both);
  subbands[1]->synthesis(both);
  subbands[2]->synthesis(both);
  interleave();
  oriented_synthesis(both);
  return;
}
void dwtnode::hpfprelift_encode(bool halfres, bool adapt)
{
	hpf_oriented_analysis(both,adapt);
	if (!halfres)
		synthesis(both);
	rawlwrite("tmp\\out.rawl");
	return;
}
void dwtnode::hpfprelift_encode(int depth, int layer, bool adapt)
{
  dwtnode *curr = this;
  for (int i=0;i<depth;i++,curr=curr->subbands[0])
  {
    curr->hpf_oriented_analysis(both,adapt);
    curr->extract_subband(0);
  }
  for (int i=depth-1;i<=layer;i--)
  {

  }
}
void dwtnode::hpfprelift_decode(char *bitrate, bool halfres, bool adapt)
{
  rawl_decode(bitrate,halfres,adapt);
  if (!halfres)
  {
    analysis(both);
    hpf_oriented_synthesis(both,adapt);
  }
	return;
}
void dwtnode::hpfprelift_decode(char *bitrate, int depth, int layer, bool adapt)
{
}
void dwtnode::hpfprelift_2layer_encode(bool halfres, bool adapt)
{
  hpf_oriented_analysis(both,adapt);
  extract_subband(0); // ofield automatically inherited
  subbands[0]->hpf_oriented_analysis(both,adapt);
  //subbands[0]->oriented_analysis(both);
  subbands[0]->synthesis(both);
  if (halfres)
  {
    subbands[0]->rawlwrite("tmp\\out.rawl");
    return;
  }
  interleave();
  synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::hpfprelift_2layer_decode(char *bitrate, bool halfres, bool adapt)
{
  rawl_decode(bitrate,halfres,adapt);
  analysis(both);
  if (halfres)
  {
    hpf_oriented_synthesis(both,adapt);
    return;
  }
  extract_subband(0);
  subbands[0]->analysis(both);
  subbands[0]->hpf_oriented_synthesis(both,adapt);
  //subbands[0]->oriented_synthesis(both);
  interleave();
  hpf_oriented_synthesis(both,adapt);
  return;
}