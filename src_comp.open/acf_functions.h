#ifndef _acf_functions_h_
#define _acf_functions_h_

#include <stdio.h>
#include <stdlib.h>
#include "complex_pack.h"
#include "math.h"

/* double gaussian(double awidth, int nfft, double deltat, gsl_complex *xarr1, gsl_complex *xarr2, gsl_complex *xarr3) */
double gaussian(double, int, double, gsl_complex *, gsl_complex *, gsl_complex *);
/* void autocorr(int nfft, gsl_complex *xarr) */
void autocorr(int, gsl_complex *);
/* void autocorr_smooth(int nfft, double deltat, double smooth, gsl_complex *xarr) */
void autocorr_smooth(int, double, double, gsl_complex *);
/* void ratio_spectral(int nfft, gsl_complex *numer, gsl_complex *denom, gsl_complex *ratio) */
void ratio_spectral(int, gsl_complex *, gsl_complex *, gsl_complex *);
/* void ratio_simplediv(int nfft, gsl_complex *numer, gsl_complex *denom, gsl_complex *ratio) */
void ratio_simplediv(int, gsl_complex *, gsl_complex *, gsl_complex *);

/* void hanning(int npts, float *window) */
void hanning(int, float *);
/* float area_func(int npts, float func[]) */
float area_func(int, float []);
/* void autocorr_smooth_func(int nfft, int nint_half, float *window, gsl_complex *xarr) */
void autocorr_smooth_func(int, int, float *, gsl_complex *);
/* following the definition of Tauzin 2019 employing the effect of noise.
void autocorr_smooth_t2019_noise(int nfft, double deltat, double smooth, gsl_complex *xarr, double snr)
*/
void autocorr_smooth_t2019_noise(int, double, double, gsl_complex *, double);
/* noise employed to ACF without whitening 
void autocorr_noise(int nfft, gsl_complex *xarr, double snr)
*/
void autocorr_noise(int, gsl_complex *, double);

#endif
