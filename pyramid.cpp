#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
static double mpegB[] = {2,0,-4,-3,5,19,26,19,5,-3,-4,0,2}; // 13-tap downsampler
static const int mpegN = 6;
static double h264[] = {1,0,-5,0,20,32,20,0,-5,0,1}; // 11-tap upsampler
static const int h264N = 5;
void dwtnode::upsample_lift(bool analysis)
{
  int sign = analysis?+1:-1;
  dwtnode *tmp = new dwtnode(h,w,txbase,true);
  dwtnode *tmp2 = new dwtnode(h,w,txbase,true);
  // copy coarse band into tmp to "upsample" with zeros
  for (int y=0;y<subbands[0]->h;y++)
    for (int x=0;x<subbands[0]->w;x++)
      tmp->pixels[2*y*w+2*x] = subbands[0]->pixels[y*subbands[0]->w+x];
  for (int y=0;y<h;y++) // vertical filter
    for (int x=0;x<w;x++)
      tmp2->pixels[y*w+x] = tmp->filt(&h264[h264N],y*w+x,0,h264N,vertical,true)/32;
  for (int y=0;y<h;y++) // horizontal filter
    for (int x=0;x<w;x++)
      pixels[y*w+x] -= sign*tmp2->filt(&h264[h264N],y*w+x,0,h264N,horizontal,true)/32;
  delete[] tmp;
  delete[] tmp2;
  return;
}
void dwtnode::downsample_lift(bool analysis)
{
  double sign = (analysis)?+1:-1;
  dwtnode *tmp = new dwtnode((h+1)/2,w,txbase,true);
  // vertical filtering first; skip every 2nd row
  for (int y=0;y<subbands[0]->h;y++)
    for (int x=0;x<w;x++)
      tmp->pixels[y*w+x] = filt(&mpegB[mpegN],2*y*w+x,0,mpegN,vertical,true)/64;
  // horizontal filter; skip every 2nd row and column
  for (int y=0;y<subbands[0]->h;y++)
    for (int x=0;x<subbands[0]->w;x++)
      subbands[0]->pixels[y*subbands[0]->w+x]
        += sign*tmp->filt(&mpegB[mpegN],y*w+2*x,0,mpegN,horizontal,true)/64;
  delete[] tmp;
  return;
}
void dwtnode::pyramid_encode(bool halfres, bool adapt)
{
  this->subbands[0] = new dwtnode((h+1)/2,(w+1)/2,txbase,true);
  downsample_lift(true);
  upsample_lift(true);
  if (halfres)
    subbands[0]->rawlwrite("tmp\\out.rawl");
  else
  {
    rawlwrite("tmp\\diff.rawl");
    subbands[0]->rawlwrite("tmp\\coarse.rawl");
  }
  return;
}
void dwtnode::pyramid_encode(int depth, int layer, bool adapt)
{
  switch (txbase)
  {
  case pyramid:
  case pyramid2x3:
  case pyramid3x2:
    break;
  default:
    cerr << "Non-pyramid transform defined!" << endl;
    exit(1);
  }
  dwtnode *curr = this;
  for (int i=0;i<depth;i++,curr=curr->subbands[0])
  {
    if (curr->subbands[0]!=nullptr)
      delete subbands[0];
    curr->subbands[0] = new dwtnode((curr->h+1)/2,(curr->w+1)/2,txbase,true);
    curr->analysis(both);
    string d_fname = "tmp\\diff";
    d_fname += i + ".rawl";
    if (i >= layer) // write out all layers needed to reconstruct
      curr->rawlwrite(d_fname.c_str());
  }
  curr->rawlwrite("tmp\\coarse.rawl");
  return;
}
void dwtnode::pyramid_decode(char *bitrate, bool halfres, bool adapt)
{
  if (halfres)
  {
    rawl_decode(bitrate,halfres,adapt);
    return;
  }
  string d_fname = "tmp\\diff";
  d_fname = d_fname + bitrate + ".rawl";
  rawlread(d_fname.c_str());
  string c_fname = "tmp\\coarse";
  c_fname = c_fname + bitrate + ".rawl";
  this->subbands[0] = new dwtnode(c_fname.c_str(),(h+1)/2,(w+1)/2,txbase);
  upsample_lift(false);
  return;
}
void dwtnode::pyramid_decode(char *bitrate, int depth, int layer, bool adapt)
{
  switch (txbase)
  {
  case pyramid:
  case pyramid2x3:
  case pyramid3x2:
    break;
  default:
    cerr << "Non-pyramid transform defined!" << endl;
    exit(1);
  }
  dwtnode *curr = this;
  for (int i=layer;i<depth;i++,curr=curr->subbands[0])
  {
    string d_fname = "tmp\\diff";
    d_fname += i;
    d_fname = d_fname + bitrate + ".rawl";
    rawlread(d_fname.c_str());
    curr->rawlread(d_fname.c_str());
    if (curr->subbands[0]!=nullptr)
      delete subbands[0];
    curr->subbands[0] = new dwtnode((curr->h+1)/2,(curr->w+1)/2,txbase,true);
  }
  string c_fname = "tmp\\coarse";
  c_fname = c_fname + bitrate + ".rawl";
  curr->rawlread(c_fname.c_str());
  for (int i=depth-1;i>=layer;i--)
  {
    curr=this;
    for (int n=0;n<i;n++)
      curr = curr->subbands[0];
    curr->synthesis(both);
  }
}
void dwtnode::lp3x2_halfres()
{
  this->subbands[0] = new dwtnode((h+1)/2,(w+1)/2,txbase,true);
  downsample_lift(true);
  upsample_lift(true);
  downsample_lift(true);
  delete [] this->pixels;
  pixels = subbands[0]->pixels;
  subbands[0]->pixels = NULL;
  delete subbands[0];
  subbands[0] = NULL;
  this->h=(h+1)/2;
  this->w=(w+1)/2;
}
void dwtnode::lp3x2_encode(bool halfres, bool adapt)
{
  this->subbands[0] = new dwtnode((h+1)/2,(w+1)/2,txbase,true);
  downsample_lift(true);
  upsample_lift(true);
  downsample_lift(true);
  if (halfres)
    subbands[0]->rawlwrite("tmp\\out.rawl");
  else
  {
    rawlwrite("tmp\\diff.rawl");
    subbands[0]->rawlwrite("tmp\\coarse.rawl");
  }
  return;
}
void dwtnode::lp3x2_decode(char *bitrate, bool halfres, bool adapt)
{
  if (halfres)
  {
    rawl_decode(bitrate,halfres,adapt);
    return;
  }
  string d_fname = "tmp\\diff";
  d_fname = d_fname + bitrate + ".rawl";
  rawlread((char *)d_fname.c_str());
  string c_fname = "tmp\\coarse";
  c_fname = c_fname + bitrate + ".rawl";
  this->subbands[0] = new dwtnode(c_fname.c_str(),(h+1)/2,(w+1)/2,txbase);
  downsample_lift(false);
  upsample_lift(false);
  return;
}
void dwtnode::lp2x3_halfres()
{
  this->subbands[0] = new dwtnode((h+1)/2,(w+1)/2,txbase,true);
  downsample_lift(true);
  delete [] this->pixels;
  pixels = subbands[0]->pixels;
  subbands[0]->pixels = NULL;
  delete subbands[0];
  subbands[0] = NULL;
  this->h=(h+1)/2;
  this->w=(w+1)/2;
}
void dwtnode::lp2x3_encode(bool halfres, bool adapt)
{
  this->subbands[0] = new dwtnode((h+1)/2,(w+1)/2,txbase,true);
  downsample_lift(true);
  upsample_lift(true);
  if (halfres)
    subbands[0]->rawlwrite("tmp\\out.rawl");
  else
  {
    rawlwrite("tmp\\diff.rawl");
    subbands[0]->rawlwrite("tmp\\coarse.rawl");
  }
  return;
}
void dwtnode::lp2x3_decode(char *bitrate, bool halfres, bool adapt)
{
  if (halfres)
  {
    rawl_decode(bitrate,halfres,adapt);
    return;
  }
  string d_fname = "tmp\\diff";
  d_fname = d_fname + bitrate + ".rawl";
  rawlread((char *)d_fname.c_str());
  string c_fname = "tmp\\coarse";
  c_fname = c_fname + bitrate + ".rawl";
  this->subbands[0] = new dwtnode(c_fname.c_str(),(h+1)/2,(w+1)/2,txbase);
  upsample_lift(false);
  downsample_lift(false);
  upsample_lift(false);
  return;
}