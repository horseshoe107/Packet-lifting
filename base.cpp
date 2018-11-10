#include "stdafx.h"
#include "base.h"
#include "dwtnode.h"
double * convolve(double *f1, int n1, double *f2, int n2)
{
  int n = n1+n2;
  double *conv_coeffs = new double[2*n+1];
  double *f = conv_coeffs+n;
  for (int i=-n;i<=n;i++)
    f[i]=0;
  for (int i=-n1;i<=n1;i++)
    for (int j=-n2;j<=n2;j++)
      f[i+j] += f1[i]*f2[j];
  return f;
}
direction operator!(direction dir)
{
  switch (dir)
  {
    case vertical:
      return horizontal;
    case horizontal:
      return vertical;
    default:
      cerr << "Inverse of this direction has no meaning" << endl;
      exit(1);
  }
}