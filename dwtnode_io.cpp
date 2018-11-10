#include "stdafx.h"
#include <ctype.h>
#include "base.h"
#include "dwtnode.h"
dwtnode::dwtnode(int hset, int wset, dwttype type, bool initzero)
{
  dwtlevel[vertical]=0, dwtlevel[horizontal]=0;
  dwtbase=type, packf=NULL;
  h=hset; w=wset;
  pixels = new double[h*w];
  for (int n=0;n<4;n++)
    subbands[n] = NULL;
  if (initzero)
    for (int n=0;n<h*w;n++)
      pixels[n]=0;
}
// initialise with data from a file; currently only works for pgm
dwtnode::dwtnode(char *fname, dwttype type)
{
  dwtlevel[vertical]=0, dwtlevel[horizontal]=0;
  dwtbase=type, pixels=NULL, packf=NULL;
  if (strcmp(fname+strlen(fname)-4,".pgm") == 0)
    pgmread(fname);
  else
  {
    cerr << "File format not supported" << endl;
    exit(1);
  }
  for (int n=0;n<4;n++)
    subbands[n] = NULL;
}
dwtnode::dwtnode(char *fname, int hset, int wset, dwttype type, int expi)
{
  dwtlevel[vertical]=0, dwtlevel[horizontal]=0;
  dwtbase=type, pixels=NULL, packf=NULL;
  if (strcmp(fname+strlen(fname)-5,".rawl") == 0)
    rawlread(fname,hset,wset,expi);
  else
  {
    cerr << "File format not supported" << endl;
    exit(1);
  }
  for (int n=0;n<4;n++)
    subbands[n] = NULL;
}
void dwtnode::copy(dwtnode &target)
{
  this->h=target.h;
  this->w=target.w;
  dwtlevel[vertical]=target.dwtlevel[vertical];
  dwtlevel[horizontal]=target.dwtlevel[horizontal];
  dwtbase=target.dwtbase;
  if (pixels!=NULL)
    delete [] pixels;
  pixels = new double[h*w];
  for (int n=0;n<h*w;n++)
    pixels[n]=target.pixels[n];
  this->packf=target.packf;
  this->ofield.copy(target.ofield);
}
void chk_pgm_comment(ifstream &in, bool eatspace)
{
  char c;
  while (true)
  {
    c = in.peek();
    if (eatspace && isspace(c)) // eat whitespace
    {
      in.ignore();
      continue;
    }
    else if (c=='#') // comment
    { // ignore until '\n' encountered
      while (in.get(c))
        if (c=='\n')
          break;
      continue;
    }
    else
      return; // useful character; quit function
  }
}
void dwtnode::pgmread(char *fname)
{
  char c[2], tmp;
  int maxval;
  ifstream pgmin(fname,ios::binary);
  if (!pgmin.good())
  {
    cerr << "Access of file " << fname << " unsuccessful." << endl;
    exit(1);
  }
  pgmin >> c[0] >> c[1]; // read header
  if ((c[0]!='P')||(c[1]!='5'))
  { // "P5" is the recognised header code for pgm files
    cerr << "File does not have a valid PGM header." << endl;
    exit(1);
  }
  chk_pgm_comment(pgmin, true);
  pgmin >> w;
  chk_pgm_comment(pgmin, true);
  pgmin >> h;
  chk_pgm_comment(pgmin, true);
  pgmin >> maxval;
  if (maxval != 255)
    cerr << "Warning: PGM is not an 8-bit file" << endl;
  chk_pgm_comment(pgmin, false);
  if (pixels!=NULL) // clear any previous data in aoim object
    delete[] pixels;
  pixels = new double[h*w]; // assign space for pgm input
  pgmin.ignore(); // skip last whitespace character
  dwtlevel[vertical]=0; dwtlevel[horizontal]=0;
  for (int n=0;n<h*w;n++) // pixel values are in row-major form
  {
    pgmin.read(&tmp,1); // read image data
    if (tmp<0)
      pixels[n] = (double) tmp+256;
    else pixels[n] = (double) tmp;
  }
  pgmin.close();
  if ((h!=ofield.h)||(w!=ofield.w))
    ofield.clearfield(h,w);
  return;
}
void dwtnode::rawlread(char *fname, int hset, int wset, int expi)
{
  short t;
  h=hset; w=wset;
  dwtlevel[vertical]=0; dwtlevel[horizontal]=0;
  const double expf = 1<<expi;
  ifstream rawin(fname,ios::binary);
  if (!rawin.good())
  {
    cerr << "Access of file " << fname << " unsuccessful." << endl;
    exit(1);
  }
  if (pixels!=NULL)
    delete[] pixels;
  pixels = new double[h*w];
  for (int n=0;n<h*w;n++)
  {
    rawin.read((char *) &t,sizeof(short));
    pixels[n] = t/expf; // floating point division
  }
  rawin.close();
  if ((h!=ofield.h)||(w!=ofield.w))
    ofield.clearfield(h,w);
  return;
}
void dwtnode::rawlread(char *fname, int expi)
{
  rawlread(fname,h,w,expi);
  return;
}
bool dwtnode::yuvstreamread(ifstream &yuvin)
{
  dwtlevel[vertical]=0; dwtlevel[horizontal]=0;
  char tmp;
  if (pixels!=NULL)
    delete[] pixels;
  pixels = new double[h*w];
  for (int n=0;n<h*w;n++)
  {
    yuvin.read(&tmp,1);
    if (tmp<0)
      pixels[n] = (double) tmp+256;
    else pixels[n] = (double) tmp;
  }
  yuvin.ignore(h*w/2); // ignore colour components
  if ((h!=ofield.h)||(w!=ofield.w))
    ofield.clearfield(h,w);
  return yuvin.good();
}
void dwtnode::pgmwrite(char *fname)
{
  const int verstep = 1<<dwtlevel[vertical];
  const int horstep = 1<<dwtlevel[horizontal];
  bool warning = false;
  ofstream pgmout(fname,ios::binary);
  if (!pgmout.good())
  {
    cerr << "Access of file " << fname << " unsuccessful." << endl;
    exit(1);
  }
  // write compliant header
  pgmout << "P5" << '\n' << (w/horstep) << ' '
    << (h/verstep) << '\n' << "255" << endl;
  for (int y=0;y<h;y+=verstep)
    for (int x=0;x<w;x+=horstep)
    {
      if (pixels[y*w+x]<0)
      {
        warning = true;
        pgmout << (char) 0;
      }
      else if (pixels[y*w+x]>255)
      {
        warning = true;
        pgmout << (char) 255;
      }
      else
        pgmout << (char) floor(pixels[y*w+x]+0.5);
    }
  if (warning)
    cerr << "Warning: (pgmwrite) some pixels exceeded [0,255]" << endl;
  pgmout.close();
  return;
}
// write 16-bit values out in little-endian format. pixel
// values are prescaled up by 2^expi, with default expi=6
void dwtnode::rawlwrite(char *fname, int expi, bool allbands)
{
  const int verstep = allbands ? 1 : 1<<dwtlevel[vertical];
  const int horstep = allbands ? 1 : 1<<dwtlevel[horizontal];
  const int expf = 1<<expi; // expand by 2^6 = 64
  const int maxval = 1<<15;
  bool warning = false;
  ofstream rawout(fname,ios::binary);
  if (!rawout.good())
  {
    cerr << "Access of file " << fname << " unsuccessful." << endl;
    exit(1);
  }
  int j;
  short k;
  for (int y=0;y<h;y+=verstep)
    for (int x=0;x<w;x+=horstep)
    {
      j = (int) round(pixels[y*w+x]*expf);
      if (j<-maxval)
      {
        k = -maxval;
        warning = true;
      }
      else if (j>maxval-1)
      {
        k = maxval-1;
        warning = true;
      }
      else
        k = (short) j;
      rawout.write((char *) &k,sizeof(short));
    }
  if (warning)
    cerr << "Warning: (rawlwrite) some pixels exceeded ["
      << -maxval << "," << (maxval-1) << "] range" << endl;
  rawout.close();
  return;
}
void dwtnode::csvwrite(char *fname)
{
  const int verstep = 1<<dwtlevel[vertical];
  const int horstep = 1<<dwtlevel[horizontal];
  ofstream csvout(fname,ios::binary);
  if (!csvout.good())
  {
    cerr << "Access of file " << fname << " unsuccessful." << endl;
    exit(1);
  }
  for (int y=0;y<h;y+=verstep)
    for (int x=0;x<w;x+=horstep)
    {
      csvout << pixels[y*w+x];
      if (x==(w-1))
        csvout << endl;
      else csvout << ", ";
    }
  csvout.close();
  return;
}
void dwtnode::yuvstreamwrite(ofstream &yuvout)
{
  const int verstep = 1<<dwtlevel[vertical];
  const int horstep = 1<<dwtlevel[horizontal];
  bool warning = false;
  for (int y=0;y<h;y+=verstep)
    for (int x=0;x<w;x+=horstep)
    {
      if (pixels[y*w+x]<0)
      {
        warning = true;
        yuvout << (char) 0;
      }
      else if (pixels[y*w+x]>255)
      {
        warning = true;
        yuvout << (char) 255;
      }
      else
        yuvout << (char) floor(pixels[y*w+x]+0.5);
    }
  // write zeros to U&V components; note each colour component
  // has 1/4 the samples of the Y component
  for (int n=0;n<h*w/2/verstep/horstep;n++)
    yuvout << (char) 128;
  if (warning)
    cerr << "Warning: image exceeded pixel dynamic range" << endl;
  return;
}