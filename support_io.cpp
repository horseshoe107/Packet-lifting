#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
// Initialise antialiasing filters from external file.
// Note that although htN and hcN indicate the true filter length,
// zero padding will occur if necessary to ensure the stored
// coefficients are always odd length support. This allows
// compatibility with the dwtnode::filt() function.
// ht is zero padded on the left side, while hc (and hp, if
// support is added) is padded on the right side. Appropriate
// boundary selection (see packet_transfer and packet_cancel
// functions) must take into account this asymmetry.
packlift_filters::packlift_filters(char *fname)
{
  ht_coeff=NULL, hc_coeff=NULL, hp_coeff=NULL;
  ifstream filtin(fname,ios::binary);
  if (!filtin.good())
  {
    cerr << "Access of file " << fname << " unsuccessful." << endl;
    exit(1);
  }
  filtin >> htN;
  int n, h_even = 1-htN%2; // +1 if htN is even, 0 otherwise
  ht_coeff = new double[htN+h_even];
  for (n=0; n<htN; n++)
    filtin >> ht_coeff[n];
  if (htN%2 == 0)  ht_coeff[htN]=0; // pad right side with 0
  filtin >> hcN;
  h_even = 1-hcN%2;
  hc_coeff = new double[hcN+h_even];
  hc_coeff[0] = 0; // pad left side with 0
  for (n=h_even; n<hcN+h_even; n++)
    filtin >> hc_coeff[n];
  if (!filtin.good())
  {
    cerr << "Could not read in filters: " << fname
      << " file may be wrong format" << endl;
    exit(2);
  }
  if ((htN+hcN)%2 != 0)
  {
    cerr << "The transfer and cancellation filters must "
      << "both have even, or odd number of samples - "
      << "end-to-end filter must be linear phase." << endl;
    exit(3);
  }
  filtin >> hpN;
  hp_coeff = new double[hpN];
  if (!filtin.good())
    return; // assume hp not specified
  for (n=0;n<hpN;n++)
    filtin >> hp_coeff[n];
  int htmaxindex=0, hcmaxindex=0;
  double maxcoeff=ht_coeff[0];
  for (n=1;n<htN;n++)
    if (maxcoeff<ht_coeff[n])
    {
      maxcoeff = ht_coeff[n];
      htmaxindex = n;
    }
  maxcoeff = hc_coeff[0];
  for (n=1;n<hcN;n++)
    if (maxcoeff<hc_coeff[n])
    {
      maxcoeff = hc_coeff[n];
      hcmaxindex = n;
    }
  offset = htmaxindex - htN/2;
  offset += hcN/2 - hcmaxindex;
  offset >>= 1;
  //offset = -1; // set this to emulate old results
}
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
void orientationfield::copy(orientationfield &target)
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
  ofstream fout(fname,ios::binary);
  if (!fout.good())
  {
    cerr << "Access of file " << fname << " unsuccessful." << endl;
    exit(1);
  }
  // write out header
  fout << w << ' ' << h << ' ' << blksz << ' '
    << oprec << ' ' << maxshift << endl;
  for (int n=0;n<numblks;n++)
  {
    fout << (int)orientvec[n].hshift << ','
        << (int)orientvec[n].vshift << ' ' << endl;
  }
  return;
}
void dwtnode::load_packfilts(char *fname)
{
  if (packf!=NULL)
    delete packf;
  packf = new packlift_filters(fname);
  return;
}