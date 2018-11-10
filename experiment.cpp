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
    cerr << "Images must have identical h&w dimensions." << endl;
    exit(1);
  }
  double apix, bpix, acc=0;
  int aystep=1<<a.dwtlevel[vertical];
  int axstep=1<<a.dwtlevel[horizontal];
  int bystep=1<<b.dwtlevel[vertical];
  int bxstep=1<<b.dwtlevel[horizontal];
  for (int ay=0,by=0;ay<a.h;ay+=aystep,by+=bystep)
    for (int ax=0,bx=0;ax<a.w;ax+=axstep,bx+=bxstep)
    {
      apix = a.pixels[ay*(a.w)+ax];
      bpix = b.pixels[by*(b.w)+bx];
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
void dwtnode::halveimage()
{
  ofield.h=h=(h+1)/2;
  ofield.w=w=(w+1)/2;
  if (ofield.blksz==1)
  {
    cerr << "Warning: orientation field incorrect now" << endl;
  }
  else if (ofield.blksz%2!=0) // blksz needs to be multiple of 2
  {
    cerr << "Subsampling of shift field is not defined for ";
    cerr << "block size " << ofield.blksz << endl;
    exit(1);
  }
  else
    ofield.blksz >>= 1;
}
// call the appropriate batch file of kdu_compress and kdu_expand commands
// also writes to fout the tested resolution and dwt filters used
void dwtnode::call_batch(testmode mode, char *Cdecomp, bool halfres, ofstream &fout)
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
  else if ((mode==pyramid3x2)||(mode==pyramid2x3))
    batch << "pyramid.bat";
  else
    batch << "out.bat";
  fout <<dwtmode_strings[dwtbase]<<" Cdecomp:"<< Cdecomp << endl;
  batch <<" "<<h<<" "<<w<<" "<<dwtbase<<" \""<<Cdecomp<<"\"";
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
void dwtnode::rawl_decode(char *bitrate, bool halfres, bool adapt)
{
  string fname = "tmp\\out";
  fname += bitrate;
  fname += ".rawl";
  rawlread((char *)fname.c_str());
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
void dwtnode::orient_encode(bool halfres, bool adapt)
{
  oriented_analysis(both);
  if (!halfres)
    synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::orient_decode(char *bitrate, bool halfres, bool adapt)
{
  rawl_decode(bitrate,halfres,adapt);
  if (halfres)
    return;
  analysis(both);
  oriented_synthesis(both);
  return;
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
void dwtnode::orient2_encode(bool halfres, bool adapt)
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
void dwtnode::orient2_decode(char *bitrate, bool halfres, bool adapt)
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
void dwtnode::hpfprelift2_encode(bool halfres, bool adapt)
{
  hpf_oriented_analysis(both,adapt);
  extract_subband(0); // ofield automatically inherited
  //subbands[0]->hpf_oriented_analysis(both,adapt);
  subbands[0]->oriented_analysis(both);
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
void dwtnode::hpfprelift2_decode(char *bitrate, bool halfres, bool adapt)
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
  //subbands[0]->hpf_oriented_synthesis(both,adapt);
  subbands[0]->oriented_synthesis(both);
  interleave();
  hpf_oriented_synthesis(both,adapt);
  return;
}