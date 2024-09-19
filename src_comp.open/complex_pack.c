#include <math.h>
#include "complex_pack.h"

/* simplify some complex number operations...! *
 * it is using gsl complex routine, 
 * but only to simplify the written expression... */

/* simple arithmatic */
gsl_complex cadd(a, b)
gsl_complex a, b;
{
  gsl_complex c;
  c = gsl_complex_add(a, b);
  return c;
}

gsl_complex csub(a, b) /* c = a-b */
gsl_complex a, b;
{
  gsl_complex c;
  c = gsl_complex_sub(a, b);
  return c;
}

gsl_complex cmul(a, b)
gsl_complex a, b;
{
  gsl_complex c;
  c = gsl_complex_mul(a, b);
  return c; 
}

gsl_complex cdiv(a, b) /* c = a/b */
gsl_complex a, b;
{
  gsl_complex c;
  c = gsl_complex_div(a, b);
  return c;
}

gsl_complex amul(a, x)
double a; 
gsl_complex x;
{
  gsl_complex c; 
  c = gsl_complex_mul_real(x, a);
  return c;
}

gsl_complex adiv(a, x)
double a;
gsl_complex x;
{
  gsl_complex c;
  c = gsl_complex_div_real(x, a);
  return c;
}

double c_abs(a)
gsl_complex a;
{
  double aabs=0;
  aabs = gsl_complex_abs(a);
  return aabs;
}
