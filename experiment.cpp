#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
extern double splines[];
extern int splines_extent;
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
void orientationfield::orient_csvout()
{
  const int fieldw = w/blksz;
  ofstream fout("C:\\Program Files\\Matlab7\\work\\"
    "oriented wavelets\\orient.csv",ios::binary);
  fout << h << ',' << w << ',' << blksz << endl;
  for (int n=0;n<numblks;n++)
  {
    fout << (int) orientvec[n].hshift;
    if (n%fieldw == fieldw-1)
      fout << endl;
    else
      fout << ',';
  }
  for (int n=0;n<numblks;n++)
  {
    fout << (int) orientvec[n].vshift;
    if (n%fieldw == fieldw-1)
      fout << endl;
    else
      fout << ',';
  }
  fout.close();
  return;
}
// applies a constant shift of (sigma/lutprec) pixels in the horizontal
// direction to the entire image
void dwtnode::shift(int sigma)
{
  int z, sker;
  bool ksig;
  double *pixdest = new double[h*w];
  for (int y=0;y<h;y++)
    for (int x=0;x<w;x++)
    { // if LUT acts on subband coefficients, modify kernel selection
      kernel_selection(x,sigma,horizontal,z,sker,ksig);
      pixdest[y*w+x] = filt(splines+sker,y*w+x,-z,splines_extent,horizontal,ksig);
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
  //Legacy code chunk; this implementation reconstructs LL1
  //without access to LH1 and HL1. Substantial artifacts appear
  //due to this imperfection
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
  //in.oriented_packet_analysis(shftbase,shftpoly4);
  //in.extract_subband(0);
  //in.subbands[0]->oriented_synthesis(shftbase,shftpoly2);
  //in.subbands[0]->pgmwrite("barb_rightleg.pgm");

  //oriented_packet_analysis(shftbase,shftpoly4);
  //extract_subband(1);
  //subbands[1]->synthesis(vertical); // ordinary synthesis of HL1
  //extract_subband(2);
  //subbands[2]->synthesis(horizontal); // ordinary synthesis of LH1
  //this->operating_depth = 1; // only write to LL1 pixels
  //oriented_synthesis(shftpoly2,shftpoly4);
  //this->operating_depth = 0;
  //extract_subband(0);
  //if (out)
  //  subbands[0]->pgmwrite("tmp\\or2_LL1.pgm");
  //if (est)
  //  ;
  //subbands[0]->oriented_packet_analysis(shftbase,shftpoly4);
  //subbands[0]->packet_synthesis(both);
  //interleave();
  //synthesis(both);
  rawlwrite("tmp\\out.rawl");
  return;
}
void dwtnode::orient2_decode(char *fname, bool est)
{
  rawlread(fname,h,w);
  //Legacy code chunk
  packet_analysis(both);
  extract_subband(0);
  subbands[0]->analysis(both);
  subbands[0]->oriented_synthesis(both);
  interleave();
  oriented_packet_synthesis(both);
  //analysis(both);
  //extract_subband(0);
  //if (est)
  //  subbands[0]->ofield.init_orient("tmp\\LLorient.dat");
  //subbands[0]->packet_analysis(both);
  //subbands[0]->oriented_packet_synthesis(shftbase,shftpoly4);
  //extract_subband(1);
  //subbands[1]->analysis(vertical);
  //extract_subband(2);
  //subbands[2]->analysis(horizontal);
  //interleave(true); // suppress warnings about incompatible dwt levels
  //this->operating_depth = 1;
  //oriented_analysis(shftpoly2,shftpoly4);
  //this->operating_depth = 0;
  //oriented_packet_synthesis(shftbase,shftpoly4);
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
    cpy.ofield.orient_csvout();
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
	return;
}
void dwtnode::hpfprelift_decode(char *bitrate, bool adapt)
{
	return;
}