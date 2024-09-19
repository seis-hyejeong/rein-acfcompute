#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fft_pack.h"

/* create base trace */
/* calculate number of nfft and in the main part, we should create zero series of nfft length. */
int trace_base(double dt, double t_length)
{
  /* time = t_pre + t_length + t_pre (additional tail) */
  int nlength, nfft;

  nfft = 1;
  nlength = (int) (t_length)/dt;
  while( nfft < nlength){ nfft *=2; }
  
  return nfft;
}

/* create waveform in frequency domain */
void trace_fft(double *xtrace, int nfft, gsl_complex *xfft)
{
  int i;
  double * xtrace_double = (double *)alloca(sizeof(double)*nfft*2);

  for(i=0; i < nfft; i++)
  {
    *(xtrace_double +2*i+0) = *(xtrace+i);
    *(xtrace_double +2*i+1) =0.;
  }

  gsl_fft_complex_radix2_forward(xtrace_double, 1, nfft);

  for(i=0; i < nfft; i++)
  {
    GSL_SET_COMPLEX(xfft+i, *(xtrace_double+2*i+0), *(xtrace_double+2*i+1));
  }

  return;
}

/* do inverse fourier transform: make into time series */
/* output: xtrace */
void trace_inv_fft(gsl_complex *xfft, int nfft, double *xtrace)
{
  double *xtrace_double;
  int i;
  xtrace_double = (double *)alloca(sizeof(double)*2*nfft);

  for(i=0; i < nfft; i++)
  {
    *(xtrace_double+2*i+0) = GSL_REAL(*(xfft+i));
    *(xtrace_double+2*i+1) = GSL_IMAG(*(xfft+i));
  }
  
  gsl_fft_complex_radix2_inverse(xtrace_double, 1, nfft);
  
  for(i=0; i < nfft; i++)
  {
    *(xtrace+i) = *(xtrace_double+2*i);
  }

  return;
}

/* multiply amplitude to waveform */
/* overwriting the trace with multiplying "amp" */
void trace_mul_amp(gsl_complex *xfft, double amp, int nfft)
{
  int i;
  for(i=0; i < nfft; i++)
  {
    *(xfft+i) = gsl_complex_mul_real(*(xfft+i), amp);
  }

  return;
}

/* give time shift to waveform */
/* it is shifting to the right
 * and multiply the given amplitude coeff. (if not, set amp=1
 * we overwrite the sequence...
 */
void trace_shift_amp(double tshift, double dt, int nfft, double amp, gsl_complex *xfft)
{
  int i;
  double df = 2*M_PI/(nfft*dt);
  
  for(i=0; i < nfft/2; i++)
  {
    *(xfft+i) = gsl_complex_mul(*(xfft+i), gsl_complex_polar(1.0, (-1)*df*i*tshift));
    *(xfft+i) = gsl_complex_mul_real(*(xfft+i), amp);
  }
  for(i= nfft/2; i < nfft; i++)
  {
    *(xfft+i) = gsl_complex_mul(*(xfft+i), gsl_complex_polar(1.0, (-1)*df*(i-nfft)*tshift));
    *(xfft+i) = gsl_complex_mul_real(*(xfft+i), amp);
  }

  return;
}

/* add two signals in frequency domain */
void trace_2add(int nfft, gsl_complex *x1, gsl_complex *x2)
{ // x2 = x1+x2
  int i;
  gsl_complex xtmp;

  for(i=0; i < nfft; i++)
  {
    xtmp = gsl_complex_add(*(x2+i), *(x1+i));
    *(x2+i) = xtmp;
  }

  return;
}

/* initialize trace */
void trace_initialize(int npts, gsl_complex *xfft)
{
  int i;
  for(i=0; i < npts; i++)
  {
    GSL_SET_COMPLEX(xfft+i, 0, 0);
  }

  return;
}

/* copy trace */
void trace_copy(int npts, gsl_complex *xorig, gsl_complex *xcopy)
{
  int i;
  for(i=0; i <npts; i++)
  {
    *(xcopy+i) = *(xorig+i);
  }
  return;
}
