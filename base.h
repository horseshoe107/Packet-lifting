#ifndef _BASE
#define _BASE
// integer division with fractions rounded to the nearest
// integer, rather than (usual) flooring towards 0.
inline int divround(int num, int den)
{
  if ((num>=0)^(den>=0)) // if only one is true
    return (num - den/2)/den;
  else // if both are positive/negative
    return (num + den/2)/den;
}
inline double round(double x){
  if (x>0) return floor(x+0.5);
  else     return ceil(x-0.5);
}
// integer division with fractions rounded towards -inf
inline int divfloor(int num, int den)
{
  if ((num<0)&&(den>0))
    return (num - den+1)/den;
  else if ((num>0)&&(den<0))
    return (num - den-1)/den;
  else
    return num/den;
}
// true modulus operator (always returns a positive value)
inline int mod(int num, int den)
{
  int a = num%den;
  while (a<0)
    a += den;
  return a;
}
double * convolve(double *, int, double *, int);
#endif