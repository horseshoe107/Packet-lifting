#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
// Initialise the orientation field. This should be called prior
// to orientation estimation, or can be used to set a simple field
// with a constant orientation
void orientationfield::init_orient
(int setblksz, int setprec, int setmaxshift, int defaulthorzshft, int defaultvertshft)
{
  blksz = setblksz;
  oprec = setprec;
  maxshift = setmaxshift;
  numblks = ((h-1)/blksz+1)*((w-1)/blksz+1);
  if (orientvec!=NULL)
    delete[] orientvec;
  orientvec = new orientation[numblks];
  for (int n=0;n<numblks;n++)
  {
    orientvec[n].hshift = defaulthorzshft;
    orientvec[n].vshift = defaultvertshft;
  }
  return;
}
// Load orientation field from an external file
void orientationfield::init_orient(char *fname)
{
  int tmp, hchk, wchk;
  std::ifstream fin(fname,std::ios::binary);
  if (!fin.good())
  {
    std::cerr << "Access of file " << fname << " unsuccessful." << std::endl;
    exit(1);
  }
  fin >> wchk >> hchk >> blksz >> oprec >> maxshift;
  if ((h!=hchk)||(w!=wchk))
  {
    std::cerr << "Orientation field must have dimensions matching the image" << std::endl;
    exit(1);
  }
  numblks=((h-1)/blksz+1)*((w-1)/blksz+1);
  if (orientvec!=NULL)
    delete[] orientvec;
  orientvec = new orientation[numblks];
  for (int n=0;n<numblks;n++)
  {
    fin >> tmp;
    orientvec[n].hshift = tmp;
    fin.ignore(10,',');
    fin >> tmp;
    orientvec[n].vshift = tmp;
  }
  return;
}
void orientationfield::clearfield(int seth, int setw)
{
  std::cerr << "Warning: resetting orientation field" << std::endl;
  h = seth;
  w = setw;
  if (orientvec!=NULL)
  {
    delete[] orientvec;
    orientvec = NULL;
  }
  return;
}
void orientationfield::copy(const orientationfield &target)
{
  this->h=target.h;
  this->w=target.w;
  blksz=target.blksz;
  numblks=target.numblks;
  oprec=target.oprec;
  maxshift=target.maxshift;
  fieldtype=target.fieldtype;
  if (orientvec!=NULL)
    delete [] orientvec;
  if (target.orientvec!=NULL)
  {
    orientvec = new orientation[numblks];
    for (int n=0;n<numblks;n++)
      orientvec[n]=target.orientvec[n];
  }
  else orientvec = NULL;
}
void orientationfield::inherit(orientationfield &parent)
{
  h=(parent.h+1)/2; w=(parent.w+1)/2;
  this->oprec     = parent.oprec;
  this->maxshift  = parent.maxshift;
  this->fieldtype = parent.fieldtype;
  if (parent.blksz==1)
  {
    std::cerr << "Warning: subsampling orientation field" << std::endl;
    this->blksz = 1;
    this->numblks = h*w;
    if (orientvec != NULL)
      delete[] orientvec;
    this->orientvec = new orientation[numblks];
    for (int y=0;y<h;y++)
      for (int x=0;x<w;x++)
      {
        orientation acc=parent.orientvec[2*y*parent.w+2*x];
        if ((2*y+1<parent.h)&&(2*x+1<parent.w))
        {
          acc += parent.orientvec[2*y*parent.w+2*x+1]
                + parent.orientvec[(2*y+1)*parent.w+2*x]
                + parent.orientvec[(2*y+1)*parent.w+2*x+1];
          acc >>= 2;
        }
        else if (2*y+1<parent.h)
        {
          acc += parent.orientvec[(2*y+1)*parent.w+2*x];
          acc >>= 1;
        }
        else if (2*x+1<parent.w)
        {
          acc += parent.orientvec[2*y*parent.w+2*x+1];
          acc >>= 1;
        }
        orientvec[y*w+x] = acc;
      }
  }
  else if (parent.blksz%2!=0) // blksz needs to be multiple of 2
  {
    std::cerr << "Subsampling of shift field is not defined for "
      << "block size " << parent.blksz << std::endl;
    exit(1);
  }
  else
  {
  this->blksz = parent.blksz/2;
  this->numblks   = parent.numblks;
  if (orientvec != NULL)
    delete[] orientvec;
  this->orientvec = new orientation[numblks];
  for (int n=0;n<numblks;n++)
    orientvec[n] = parent.orientvec[n];
  }
  return;
}
void orientationfield::orientwrite(char *fname)
{
	const int blkwidth = (w-1)/blksz + 1;
  std::ofstream fout(fname,std::ios::binary);
  if (!fout.good())
  {
    std::cerr << "Access of file " << fname << " unsuccessful." << std::endl;
    exit(1);
  }
  // write out header
  fout << w << ' ' << h << ' ' << blksz << ' '
    << oprec << ' ' << maxshift;
  for (int n=0;n<numblks;n++)
  {
		if (n%blkwidth==0)
			fout << std::endl;
    fout << (int)orientvec[n].hshift << ','
        << (int)orientvec[n].vshift << ' ';
  }
  return;
}