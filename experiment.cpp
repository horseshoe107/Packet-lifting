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
void dwtnode::halveimage(bool doit)
{
  if (doit)
  {
    this->h>>=1;
    this->w>>=1;
  }
}
// call the appropriate batch file of kdu_compress and kdu_expand commands
// also writes to fout the tested resolution and dwt filters used
void dwtnode::call_batch(testmode mode, bool halfres,ofstream &fout)
{
  char *dwtmode_strings[]={"w5x3","w9x7","nonstandard dwt"};
  std::stringstream batch;
  if (halfres)
  {
    batch << "halfres.bat";
    fout << " half resolution";
  }
  else if ((mode==pyramid3x2)||(mode==pyramid2x3))
    batch << "pyramid.bat";
  else
    batch << "out.bat";
  fout << " " << dwtmode_strings[dwtbase] << endl;
  batch <<" "<<h<<" "<<w<<" "<<dwtbase;
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
void dwtnode::rawl_decode(char *bitrate, bool adapt)
{
  string fname = "tmp\\out";
  fname += bitrate;
  fname += ".rawl";
  rawlread((char *)fname.c_str());
}
void dwtnode::antialias_encode(bool halfres, bool adapt)
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
void dwtnode::antialias_decode(char *bitrate, bool adapt)
{
  string fname = "tmp\\out";
  fname += bitrate;
  fname += ".rawl";
  rawlread((char *)fname.c_str(),h,w);
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
  oriented_packet_analysis(both);
  if (halfres)
  {
		extract_subband(0);
		subbands[0]->oriented_synthesis(both);
    subbands[0]->rawlwrite("tmp\\out.rawl");
		delete subbands[0];
		subbands[0] = NULL;
		return;
  }
  packet_synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::orient_decode(char *bitrate, bool adapt)
{
  string fname = "tmp\\out";
  fname += bitrate;
  fname += ".rawl";
  rawlread((char *)fname.c_str(),h,w);
  packet_analysis(both);
  oriented_packet_synthesis(both);
  return;
}
void dwtnode::aa_orient_encode(bool out, bool adapt)
{
  oriented_packet_analysis(both);
  extract_subband(0);
  extract_subband(1);
  extract_subband(2);
  subbands[1]->analysis(horizontal);
  subbands[2]->analysis(vertical);
  packlift(both,true,adapt);
  if (out)
  {
    subbands[1]->synthesis(horizontal);
    subbands[2]->synthesis(vertical);
    interleave();
    oriented_synthesis(both);
    pgmwrite("tmp\\aaor_LL1.pgm");
    oriented_analysis(both);
    extract_subband(0);
    extract_subband(1);
    extract_subband(2);
    subbands[1]->analysis(horizontal);
    subbands[2]->analysis(vertical);
  }
  subbands[1]->synthesis(horizontal);
  subbands[2]->synthesis(vertical);
  interleave();
  packet_synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::aa_orient_decode(char *bitrate, bool adapt)
{
  string fname = "tmp\\out";
  fname += bitrate;
  fname += ".rawl";
  rawlread((char *)fname.c_str(),h,w);
  packet_analysis(both);
  extract_subband(0);
  extract_subband(1);
  extract_subband(2);
  subbands[1]->analysis(horizontal);
  subbands[2]->analysis(vertical);
  packlift(both,false,adapt);
  subbands[1]->synthesis(horizontal);
  subbands[2]->synthesis(vertical);
  interleave();
  oriented_packet_synthesis(both);
  return;
}
void dwtnode::orient2_encode(bool out, bool est)
{
  oriented_packet_analysis(both);
  extract_subband(0);
  if (out)
  {
    subbands[0]->oriented_synthesis(both);
    subbands[0]->pgmwrite("tmp\\or2_LL1.pgm");
    subbands[0]->oriented_analysis(both);
  }
  // only valid when the field is inherited
  // the other approach is to do an oriented synthesis back to LL1,
  // then packet synthesis down using a new orientation field
  subbands[0]->oriented_analysis(both);
  subbands[0]->synthesis(both); // synthesise back up to LL1
  interleave();
  packet_synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::orient2_decode(char *fname, bool est)
{
  rawlread(fname,h,w);
  packet_analysis(both);
  extract_subband(0);
  subbands[0]->analysis(both);
  subbands[0]->oriented_synthesis(both);
  interleave();
  oriented_packet_synthesis(both);
  return;
}
void dwtnode::aa_orient2_encode(bool out, bool est)
{
  oriented_packet_analysis(both);
  extract_subband(0);
  extract_subband(1);
  extract_subband(2);
  subbands[1]->analysis(horizontal);
  subbands[2]->analysis(vertical);
  packlift(both,true);
  subbands[1]->synthesis(horizontal);
  subbands[2]->synthesis(vertical);
  if (est)
  {
    estorient cpy(subbands[0]);
    cpy.init_orient(4,8,16);
    cpy.calc_energies();
    cpy.choose_orient();
    cpy.ofield.orientwrite("tmp\\LLaaorient.dat");
    subbands[0]->ofield.init_orient("tmp\\LLaaorient.dat");
  }
  subbands[0]->synthesis(both);
  if (out)
    subbands[0]->pgmwrite("tmp\\aaor2_LL.pgm");
  subbands[0]->oriented_packet_analysis(both);
  subbands[0]->synthesis(both);
  interleave();
  packet_synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::aa_orient2_decode(char *fname, bool est)
{
  rawlread(fname,h,w);
  analysis(both);
  extract_subband(0);
  if (est)
    subbands[0]->ofield.init_orient("tmp\\LLaaorient.dat");
  subbands[0]->packet_analysis(both);
  subbands[0]->oriented_packet_synthesis(both);
  interleave();
  oriented_analysis(both);
  extract_subband(0); // undo packet lifting
  extract_subband(1);
  extract_subband(2);
  subbands[1]->analysis(both);
  subbands[2]->analysis(both);
  packlift(both,false);
  subbands[1]->synthesis(horizontal);
  subbands[2]->synthesis(vertical);
  interleave();
  oriented_packet_synthesis(both);
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
void dwtnode::hpfprelift_decode(char *bitrate, bool adapt)
{
  string fname = "tmp\\out";
  fname += bitrate;
  fname += ".rawl";
  rawlread((char *)fname.c_str(),h,w);
	analysis(both);
	hpf_oriented_synthesis(both,adapt);
	return;
}