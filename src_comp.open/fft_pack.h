#ifndef _fft_pack_h_
#define _fft_pack_h_

#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_fft.h"
#include "gsl/gsl_fft_real.h"
#include "gsl/gsl_fft_complex.h"

int trace_base(double, double);
void trace_fft(double *, int, gsl_complex *);
void trace_inv_fft(gsl_complex *, int, double *);
void trace_mul_amp(gsl_complex *, double, int);
void trace_shift_amp(double, double, int, double, gsl_complex *);
void trace_2add(int, gsl_complex *, gsl_complex *);
void trace_initialize(int, gsl_complex *);
void trace_copy(int, gsl_complex *, gsl_complex *);

#endif
