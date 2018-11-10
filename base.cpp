#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
double log_lut[8]; // approximations are rounded down
static class base_init
{
public:
  base_init()
  { // initialise log_lut
    for (float i=0;i<8;i++)
      log_lut[(int)i] = log(1+i/8);
  }
} _do_it;
// a<1: returns a, a>1: returns ln(a)+1
double gamma(double a)
{
  if (a<0)
  {
    cerr << "log function cannot be called on a negative number" << endl;
    exit(1);
  }
  if (a<1) // linear section; return unchanged
    return a;
  const double ln2 = log(2.0);
  int i=0,j=0;
  while (a>UINT_MAX){a/=2; i++;}
  unsigned int b=(unsigned int)a;
  while (b!=1){b>>=1; j++;}
  if (j>=3) // when a>=2^3
  {
    b = (unsigned int)a;
    b >>= (j-3);
  }
  else // when 2^3>a>1
  {
    for (int n=j;n<3;n++)
      a *= 2;
    b = (unsigned int)a;
  }
  b &= 0x7;
  return (i+j)*ln2+log_lut[b]+1;
}