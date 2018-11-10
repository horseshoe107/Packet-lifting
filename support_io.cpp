#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
// Initialise the orientation field with preset shifts. This
// should be called prior to orientation estimation.
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
  for (int n=0;n<numblks;n++) {
    orientvec[n].hshift = defaulthorzshft;
    orientvec[n].vshift = defaultvertshft;
  }
  return;
}
// Load orientation field from an external file
void orientationfield::init_orient(char *fname)
{
  int tmp, hchk, wchk;
  ifstream fin(fname,ios::binary);
  if (!fin.good())
  {
    cerr << "Access of file " << fname << " unsuccessful." << endl;
    exit(1);
  }
  fin >> wchk >> hchk >> blksz >> oprec >> maxshift;
  if ((h!=hchk)||(w!=wchk))
  {
    cerr << "Orientation field must have dimensions matching the image" << endl;
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
  orientvec = new orientation[numblks];
  for (int n=0;n<numblks;n++)
    orientvec[n]=target.orientvec[n];
}
void orientationfield::inherit(orientationfield parent)
{
  if (parent.blksz >=2)
    this->blksz = parent.blksz/2;
  else // fix this later
  {
    cerr << "Cannot inherit from parent; shift field is too fine" << endl;
    exit(1);
  }
  this->fieldtype = parent.fieldtype;
  this->maxshift  = parent.maxshift;
  this->numblks   = parent.numblks;
  this->oprec     = parent.oprec;
  if (orientvec != NULL)
    delete[] orientvec;
  this->orientvec = new orientation[numblks];
  for (int n=0;n<numblks;n++)
    orientvec[n] = parent.orientvec[n];
  return;
}
void orientationfield::orientwrite(char *fname)
{
	const int blkwidth = (w-1)/blksz + 1;
  ofstream fout(fname,ios::binary);
  if (!fout.good())
  {
    cerr << "Access of file " << fname << " unsuccessful." << endl;
    exit(1);
  }
  // write out header
  fout << w << ' ' << h << ' ' << blksz << ' '
    << oprec << ' ' << maxshift;
  for (int n=0;n<numblks;n++)
  {
		if (n%blkwidth==0)
			fout << endl;
    fout << (int)orientvec[n].hshift << ','
        << (int)orientvec[n].vshift << ' ';
  }
  return;
}