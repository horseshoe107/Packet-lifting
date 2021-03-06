#include "stdafx.h"
#include <ctype.h>
#include "base.h"
#include "dwtnode.h"
dwtnode::dwtnode(int hset, int wset, txtype type, bool initzero)
{
  ofield.h=h=hset, ofield.w=w=wset, txbase=type;
  dwtlevel[vertical]=0, dwtlevel[horizontal]=0;
  pixels = new double[h*w];
  if (initzero)
    for (int n=0;n<h*w;n++)
      pixels[n]=0;
  for (int n=0;n<4;n++)
    subbands[n] = NULL;
}
// initialise with data from a file; currently only works for pgm
dwtnode::dwtnode(const char *fname, txtype type)
{
  dwtlevel[vertical]=0, dwtlevel[horizontal]=0;
  txbase=type, pixels=NULL;
  for (int n=0;n<4;n++)
    subbands[n] = NULL;
  if (strcmp(fname+strlen(fname)-4,".pgm") == 0)
    pgmread(fname);
  else
  {
    std::cerr << "File format not supported" << std::endl;
    exit(1);
  }
}
dwtnode::dwtnode(const char *fname, int hset, int wset, txtype type, int expi)
{
  dwtlevel[vertical]=0, dwtlevel[horizontal]=0;
  txbase=type, pixels=NULL;
  for (int n=0;n<4;n++)
    subbands[n] = NULL;
  if (strcmp(fname+strlen(fname)-5,".rawl") == 0)
    rawlread(fname,hset,wset,expi);
  else
  {
    std::cerr << "File format not supported" << std::endl;
    exit(1);
  }
}
dwtnode::dwtnode(const dwtnode &target)
{
  this->h=target.h;
  this->w=target.w;
  dwtlevel[vertical]=target.dwtlevel[vertical];
  dwtlevel[horizontal]=target.dwtlevel[horizontal];
  txbase=target.txbase;
  pixels = new double[h*w];
  for (int n=0;n<h*w;n++)
    pixels[n]=target.pixels[n];
  for (int n=0;n<4;n++)
  {
    if (target.subbands[n]==NULL)
      this->subbands[n]=NULL;
    else
      this->subbands[n]= new dwtnode(*target.subbands[n]);
  }
	ofield.copy(target.ofield);
}
dwtnode& dwtnode::operator=(const dwtnode &target)
{
  if (this == &target)
    return *this;
  this->h=target.h;
  this->w=target.w;
  dwtlevel[vertical]=target.dwtlevel[vertical];
  dwtlevel[horizontal]=target.dwtlevel[horizontal];
  txbase=target.txbase;
  if (pixels!=NULL)
    delete[] pixels;
  pixels = new double[h*w];
  for (int n=0;n<h*w;n++)
    pixels[n]=target.pixels[n];
  for (int n=0;n<4;n++)
  {
    if (this->subbands[n]!=NULL)
      delete this->subbands[n];
    if (target.subbands[n]==NULL)
      this->subbands[n]=NULL;
    else
      this->subbands[n]= new dwtnode(*target.subbands[n]);
  }
	ofield.copy(target.ofield);
  return *this;
}
void chk_pgm_comment(std::ifstream &in, bool eatspace)
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
void dwtnode::pgmread(const char *fname)
{
  char c[2], tmp;
  int maxval;
  std::ifstream pgmin(fname,std::ios::binary);
  if (!pgmin.good())
  {
    std::cerr << "Access of file " << fname << " unsuccessful." << std::endl;
    exit(1);
  }
  pgmin >> c[0] >> c[1]; // read header
  if ((c[0]!='P')||(c[1]!='5'))
  { // "P5" is the recognised header code for pgm files
    std::cerr << "File does not have a valid PGM header." << std::endl;
    exit(1);
  }
  chk_pgm_comment(pgmin, true);
  pgmin >> w;
  chk_pgm_comment(pgmin, true);
  pgmin >> h;
  chk_pgm_comment(pgmin, true);
  pgmin >> maxval;
  if (maxval != 255)
    std::cerr << "Warning: PGM is not an 8-bit file" << std::endl;
  chk_pgm_comment(pgmin, false);
  dwtlevel[vertical]=0; dwtlevel[horizontal]=0;
  if (pixels!=NULL)
    delete[] pixels;
  pixels = new double[h*w];
  pgmin.ignore(); // skip last whitespace character
  for (int n=0;n<h*w;n++) // pixel values are in row-major form
  {
    pgmin.read(&tmp,1); // read image data
    if (tmp<0) // data is stored as unsigned but will be read as signed
      pixels[n] = (double) tmp+256;
    else pixels[n] = (double) tmp;
  }
  pgmin.close();
  for (int n=0;n<4;n++)
    if (subbands[n]!=NULL)
    {
      delete subbands[n];
      subbands[n]=NULL;
    }
  if ((h!=ofield.h)||(w!=ofield.w))
    ofield.clearfield(h,w);
  return;
}
void dwtnode::rawlread(const char *fname, int hset, int wset, int expi)
{
  short t;
  h=hset; w=wset;
  dwtlevel[vertical]=0; dwtlevel[horizontal]=0;
  const double expf = 1<<expi;
  std::ifstream rawin(fname,std::ios::binary);
  if (!rawin.good())
  {
    std::cerr << "Access of file " << fname << " unsuccessful." << std::endl;
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
  for (int n=0;n<4;n++)
    if (subbands[n]!=NULL)
    {
      delete subbands[n];
      subbands[n]=NULL;
    }
  if ((h!=ofield.h)||(w!=ofield.w))
    ofield.clearfield(h,w);
  return;
}
void dwtnode::rawlread(const char *fname, int expi)
{  rawlread(fname,h,w,expi);  }
// stream object must be binary: ifstream fin(fname,ios::binary);
bool dwtnode::yuvstreamread(std::ifstream &yuvin)
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
  for (int n=0;n<4;n++)
    if (subbands[n]!=NULL)
    {
      delete subbands[n];
      subbands[n]=NULL;
    }
  if ((h!=ofield.h)||(w!=ofield.w))
    ofield.clearfield(h,w);
  return yuvin.good();
}
bool dwtnode::uyvystreamread(std::ifstream &uyvyin)
{
  dwtlevel[vertical]=0; dwtlevel[horizontal]=0;
  char tmp;
  if (pixels!=NULL)
    delete[] pixels;
  pixels = new double[h*w];
  for (int n=0;n<h*w;n++)
  {
    uyvyin.ignore(1); // ignore colour components
    uyvyin.read(&tmp,1);
    if (tmp<0)
      pixels[n] = (double) tmp+256;
    else pixels[n] = (double) tmp;
  }
  for (int n=0;n<4;n++)
    if (subbands[n]!=NULL)
    {
      delete subbands[n];
      subbands[n]=NULL;
    }
  if ((h!=ofield.h)||(w!=ofield.w))
    ofield.clearfield(h,w);
  return uyvyin.good();
}
void dwtnode::pgmwrite(const char *fname)
{
  const int verstep = 1<<dwtlevel[vertical];
  const int horstep = 1<<dwtlevel[horizontal];
  bool warning = false;
  std::ofstream pgmout(fname,std::ios::binary);
  if (!pgmout.good())
  {
    std::cerr << "Access of file " << fname << " unsuccessful." << std::endl;
    exit(1);
  }
  // write compliant header
  pgmout << "P5" << '\n' << (w/horstep) << ' '
    << (h/verstep) << '\n' << "255" << std::endl;
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
    std::cerr << "Warning: (pgmwrite) some pixels exceeded [0,255]" << std::endl;
  pgmout.close();
  return;
}
// write 16-bit values out in little-endian format. pixel
// values are prescaled up by 2^expi, with default expi=6
void dwtnode::rawlwrite(const char *fname, int expi, bool allbands)
{
  const int verstep = allbands ? 1 : 1<<dwtlevel[vertical];
  const int horstep = allbands ? 1 : 1<<dwtlevel[horizontal];
  const int expf = 1<<expi; // expand by 2^6 = 64
  const int maxval = 1<<15;
  bool warning = false;
  std::ofstream rawout(fname,std::ios::binary);
  if (!rawout.good())
  {
    std::cerr << "Access of file " << fname << " unsuccessful." << std::endl;
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
    std::cerr << "Warning: (rawlwrite) some pixels exceeded ["
      << -maxval << "," << (maxval-1) << "] range" << std::endl;
  rawout.close();
  return;
}
void dwtnode::csvwrite(const char *fname)
{
  const int verstep = 1<<dwtlevel[vertical];
  const int horstep = 1<<dwtlevel[horizontal];
  std::ofstream csvout(fname,std::ios::binary);
  if (!csvout.good())
  {
    std::cerr << "Access of file " << fname << " unsuccessful." << std::endl;
    exit(1);
  }
  for (int y=0;y<h;y+=verstep)
    for (int x=0;x<w;x+=horstep)
    {
      csvout << pixels[y*w+x];
      if (x==(w-1))
        csvout << std::endl;
      else csvout << ", ";
    }
  csvout.close();
  return;
}
void dwtnode::yuvstreamwrite(std::ofstream &yuvout)
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
    std::cerr << "Warning: image exceeded pixel dynamic range" << std::endl;
  return;
}