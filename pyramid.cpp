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
  dwtnode *tmp = new dwtnode(h,w,disabled,true);
  dwtnode *tmp2 = new dwtnode(h,w,disabled,true);
  // copy coarse band into tmp to "upsample" with zeros
  for (int y=0;y<subbands[0]->h;y++)
    for (int x=0;x<subbands[0]->w;x++)
      tmp->pixels[2*y*w+2*x] = subbands[0]->pixels[y*subbands[0]->w+x];
  // vertical filter
  for (int y=0;y<h;y++)
    for (int x=0;x<w;x++)
      tmp2->pixels[y*w+x] = tmp->filt(&h264[h264N],y*w+x,0,h264N,vertical,true)/32;
  for (int y=0;y<h;y++)
    for (int x=0;x<w;x++)
      pixels[y*w+x] -= sign*tmp2->filt(&h264[h264N],y*w+x,0,h264N,horizontal,true)/32;
  delete[] tmp;
  delete[] tmp2;
  return;
}
void dwtnode::downsample_lift(bool analysis)
{
  double sign = (analysis)?+1:-1;
  dwtnode *tmp = new dwtnode((h+1)/2,w,disabled,true);
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
void dwtnode::lp3x2_halfres()
{
  dwtbase = disabled;
  this->subbands[0] = new dwtnode((h+1)/2,(w+1)/2,disabled,true);
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
  dwtbase = disabled;
  this->subbands[0] = new dwtnode((h+1)/2,(w+1)/2,disabled,true);
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
void dwtnode::lp3x2_decode(char *bitrate, bool adapt)
{
  dwtbase = disabled;
  string c_fname = "tmp\\coarse";
  c_fname = c_fname + bitrate + ".rawl";
  subbands[0]->rawlread((char *)c_fname.c_str());
  string d_fname = "tmp\\diff";
  d_fname = d_fname + bitrate + ".rawl";
  rawlread((char *)d_fname.c_str());
  downsample_lift(false);
  upsample_lift(false);
  return;
}
void dwtnode::lp2x3_halfres()
{
  dwtbase = disabled;
  this->subbands[0] = new dwtnode((h+1)/2,(w+1)/2,disabled,true);
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
  dwtbase = disabled;
  this->subbands[0] = new dwtnode((h+1)/2,(w+1)/2,disabled,true);
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
void dwtnode::lp2x3_decode(char *bitrate, bool adapt)
{
  dwtbase = disabled;
  string c_fname = "tmp\\coarse";
  c_fname = c_fname + bitrate + ".rawl";
  subbands[0]->rawlread((char *)c_fname.c_str());
  string d_fname = "tmp\\diff";
  d_fname = d_fname + bitrate + ".rawl";
  rawlread((char *)d_fname.c_str());
  upsample_lift(false);
  downsample_lift(false);
  upsample_lift(false);
  return;
}