#include "acf_functions.h" 

double gaussian(double awidth, int nfft, double deltat, gsl_complex *xarr1, gsl_complex *xarr2, gsl_complex *xarr3)
{
  int i;
  double df, gaussamp, normalizer = 0;

  df = 2*M_PI/(nfft*deltat);
  for(i=0; i < nfft/2; i++)
  {
    gaussamp = exp(-pow(df*i,2)/(4*awidth*awidth)); 
    xarr1[i] = amul(gaussamp, xarr1[i]);
    xarr2[i] = amul(gaussamp, xarr2[i]);
    xarr3[i] = amul(gaussamp, xarr3[i]);
    normalizer += gaussamp;
  }
  for(i=nfft/2; i < nfft; i++)
  {
    gaussamp = exp(-pow(df*(i-nfft),2)/(4*awidth*awidth));
    xarr1[i] = amul(gaussamp, xarr1[i]);
    xarr2[i] = amul(gaussamp, xarr2[i]);
    xarr3[i] = amul(gaussamp, xarr3[i]);
    normalizer += gaussamp;
  }

  /* invert the normalizer */
  normalizer = nfft/normalizer;

  /* return the normalizing value */
  normalizer = 1./normalizer;

  return normalizer;
}

void autocorr(int nfft, gsl_complex *xarr)
{
  int i;
  double tmpdouble;
  gsl_complex xtmp; //, xconj;

  tmpdouble = 0.;
  for(i=0; i < nfft; i++)
  {
    tmpdouble = gsl_complex_abs2(xarr[i]);
    GSL_SET_COMPLEX(&xtmp, tmpdouble, 0.);

    xarr[i] = xtmp; 
  }

  return;
}

/* autocorrelation with the white noise: consider noise only when snr > 0 */
void autocorr_noise(int nfft, gsl_complex *xarr, double snr)
{
  int i;
  double tmpdouble, nsquare = 0;
 
  tmpdouble = 0.; /* we make summation of the square amplitudes of the xarr */
  for(i=0; i < nfft; i++)
  {
    tmpdouble += gsl_complex_abs2(xarr[i]);
  }
  nsquare = tmpdouble/(nfft);
  if(snr > 0)
  {
    tmpdouble = 1./(snr*snr -1);
    nsquare *= tmpdouble;
  }else{ nsquare = 0.; }

  /* now add the nsquare to the xarr */
  for(i=0; i < nfft; i++)
  {
    tmpdouble = gsl_complex_abs2(xarr[i]);
    tmpdouble += nsquare; /* add noise */

    GSL_SET_COMPLEX(xarr+i, tmpdouble, 0.);
  }
  
  return;
}

/* autocorrelation with smoothing -- denominator is smoothed. */
void autocorr_smooth(int nfft, double deltat, double smooth, gsl_complex *xarr)
{
  int i, j, nint, nelements;
  gsl_complex xtmp, xconj;
  double df, xxtmp, *smoothing;

  smoothing = calloc(nfft, sizeof(double));
  df = 1./(nfft*deltat);
  nint = floor(smooth/df/2);
  /*fprintf(stderr, "smoothing: %.4f %d\n", smooth, nint); */
  /* get smoothing weights */
  for(i=0; i < nint; i++)
  {
    xxtmp = 0;
    nelements = nint+i+1;
    for(j=0; j <=i+nint; j++)
    {
      xxtmp += gsl_complex_abs(xarr[j]);
    }
    xxtmp = xxtmp/nelements;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i=nint; i< nfft/2-nint; i++)
  {
    xxtmp = 0;
    nelements = 2*nint+1;
    for(j = i-nint; j <= i+nint; j++)
    {
      xxtmp += gsl_complex_abs(xarr[j]);
    }
    xxtmp = xxtmp/nelements;
    smoothing[i] = xxtmp*xxtmp; 
  }
  for(i = nfft/2-nint; i <= nfft/2; i++)
  {
    xxtmp = 0;
    nelements = nint + (nfft/2-i)+1;
    for(j=i - nint; j <= nfft/2; j++)
    {
      xxtmp += gsl_complex_abs(xarr[j]);
    }
    xxtmp = xxtmp/nelements;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i= nfft/2+1; i < nfft/2+ nint; i++)
  {
    xxtmp = 0;
    nelements = nint +i - nfft/2 +1; 
    for(j = nfft/2; j <= nint+i; j++)
    {
      xxtmp += gsl_complex_abs(xarr[j]);
    }
    xxtmp = xxtmp/nelements;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i = nfft/2+ nint; i < nfft - nint; i++)
  {
    xxtmp = 0;
    nelements = 2*nint +1;
    for(j= i - nint; j <= i+nint; j++)
    {
      xxtmp += gsl_complex_abs(xarr[j]);
    }
    xxtmp = xxtmp/nelements;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i = nfft-nint; i < nfft; i++)
  {
    xxtmp = 0;
    nelements = nfft -i + nint ;
    for(j = i - nint; j < nfft; j++)
    {
      xxtmp += gsl_complex_abs(xarr[j]);
    }
    xxtmp = xxtmp/nelements;
    smoothing[i] = xxtmp*xxtmp;
  }
 
  for(i=0; i < nfft; i++)
  {
    xconj = gsl_complex_conjugate(xarr[i]);
    xtmp = cmul(xconj, xarr[i]);
    /* now smoothing!! */
    xarr[i] = adiv(smoothing[i], xtmp);
  }

  free(smoothing);
  return;
}

void ratio_spectral(int nfft, gsl_complex *numer, gsl_complex *denom, gsl_complex *ratio)
{
  int i;
  gsl_complex tmpconj, finalnumer, finaldenom;

  for(i=0; i < nfft; i++)
  {
    tmpconj = gsl_complex_conjugate(denom[i]);
    finalnumer = cmul(numer[i], tmpconj); /* H* Z' */
    finaldenom = cmul(denom[i], tmpconj); /* Z* Z' */
 
    ratio[i] = cdiv(finalnumer, finaldenom);
  }

  return;
}

void ratio_simplediv(int nfft, gsl_complex *numer, gsl_complex *denom, gsl_complex *ratio)
{ /* well, the ratio will be double (bc double divided by double), but for the sake of uniform format, we do complex */
  int i;
  double tmpcomp, tmpnumer, tmpdenom;

  for(i=0; i < nfft; i++)
  {
    tmpnumer = gsl_complex_abs(numer[i]);
    tmpdenom = gsl_complex_abs(denom[i]);
    tmpcomp = tmpnumer/tmpdenom;
    GSL_SET_COMPLEX(ratio+i, tmpcomp, 0);
  }

  return;
}

/* the hanning window with a given length */
void hanning(int npts, float *window)
{
  int i, tmpnpts;
  float x = 0;

  tmpnpts = npts -1;
  for(i=0; i < npts; i++)
  {
    x = (float) i/tmpnpts;
    *(window+i) = 0.5*(1-cos(M_PI*2*x)); 
  }
 
  return;
}

/* calculate area of the function */
float area_func(int npts, float func[])
{
  float area = 0;
  int i = 0;
  for(i = 0; i < npts; i++)
  {
    area += func[i];
  }

  return area;
}

/* smoothing using certain time window */
/* The length of the function is 2*nint_half+1 and nint_half = floor(smooth/def/2) where df = 1/(nfft*deltat) */
void autocorr_smooth_func(int nfft, int nint_half, float *window, gsl_complex *xarr)
{
  int i, j;
  gsl_complex xtmp, xconj;
  double xxtmp, *smoothing;
  float area = 0;

  smoothing = calloc(nfft, sizeof(double)); 
  for(i=0; i < nint_half; i++)
  {
    xxtmp = 0;
    area = area_func(nint_half+i+1, window + nint_half -i);
    /*fprintf(stderr, "area: %.4f\n", area); */
    for(j=0; j <=i+nint_half; j++)
    {
      xxtmp += ( gsl_complex_abs(xarr[j])*window[nint_half -i+j] );
    }
    xxtmp = xxtmp/area;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i=nint_half; i< nfft/2-nint_half; i++)
  {
    xxtmp = 0;
    area = area_func(2*nint_half +1, window);
    for(j = i-nint_half; j <= i+nint_half; j++)
    {
      xxtmp += ( gsl_complex_abs(xarr[j])*window[j-i+nint_half] );
    }
    xxtmp = xxtmp/area;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i = nfft/2-nint_half; i <= nfft/2; i++)
  {
    xxtmp = 0;
    area = area_func(nint_half + (nfft/2-i)+1, window + nint_half - (nfft/2-i));
    for(j=i - nint_half; j <= nfft/2; j++)
    {
      xxtmp += ( gsl_complex_abs(xarr[j])*window[j-i + nint_half] );
    }
    xxtmp = xxtmp/area;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i= nfft/2+1; i < nfft/2+ nint_half; i++)
  {
    xxtmp = 0;
    area = area_func(nint_half +i - nfft/2 +1, window + nint_half - i + nfft/2);
    for(j = nfft/2; j <= nint_half+i; j++)
    {
      xxtmp += ( gsl_complex_abs(xarr[j])*window[ j - i + nint_half] );
    }
    xxtmp = xxtmp/area;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i = nfft/2+ nint_half; i < nfft - nint_half; i++)
  {
    xxtmp = 0;
    area = area_func(2*nint_half + 1, window);
    for(j= i - nint_half; j <= i+nint_half; j++)
    {
      xxtmp += ( gsl_complex_abs(xarr[j])*window[j - i + nint_half ] );
    }
    xxtmp = xxtmp/area;
    smoothing[i] = xxtmp*xxtmp;
  }
  for(i = nfft-nint_half; i < nfft; i++)
  {
    xxtmp = 0;
    area = area_func(nfft -i + nint_half, window + nint_half +1 + i - nfft);
    for(j = i - nint_half; j < nfft; j++)
    {
      xxtmp += ( gsl_complex_abs(xarr[j])*window[j - i + nint_half] );
    }
    xxtmp = xxtmp/area;
    smoothing[i] = xxtmp*xxtmp;
  }

  for(i=0; i < nfft; i++)
  {
    xconj = gsl_complex_conjugate(xarr[i]);
    xtmp = cmul(xconj, xarr[i]);
    /* now smoothing!! */
    xarr[i] = adiv(smoothing[i], xtmp);
  }

  free(smoothing);
  return;
}

/* calculate the ACF with whitening definition of Tauzin 2019 */
/* autocorrelation with smoothing -- denominator is smoothed. 
   Also with a SNR. When it's -1, we don't consider that. 
   If it's a positive number, we consider it.
*/
void autocorr_smooth_t2019_noise(int nfft, double deltat, double smooth, gsl_complex *xarr, double snr)
{
  int i, j, nint, nelements;
  gsl_complex xtmp; //, xconj;
  double df, xxtmp, *smoothing, *abs2array, nsquare;

  smoothing = calloc(nfft, sizeof(double));
  abs2array = calloc(nfft, sizeof(double));
  if(smoothing == NULL || abs2array == NULL){ fprintf(stderr, "wrong!!!!\n"); exit(-1); } /* terminate right away */
  for(i=0; i < nfft; i++)
  {
    abs2array[i] = gsl_complex_abs2(xarr[i]);
  }  // write from here.

  /* get the nsquare to add to the abs2 array */
  nsquare = 0.;
  if(snr > 0)
  {
    for(i=0; i < nfft; i++){ nsquare += abs2array[i];  }
    nsquare /= nfft;
    nsquare *= 1./(snr*snr -1);

    for(i=0; i < nfft; i++){ abs2array[i] += nsquare; } /* add the noise amplitude.  */
  }

  df = 1./(nfft*deltat);
  nint = floor(smooth/df/2);
  /*fprintf(stderr, "smoothing: %.4f %d\n", smooth, nint); */
  /* get smoothing weights */
  for(i=0; i < nint; i++)
  {
    xxtmp = 0;
    nelements = nint+i+1;
    for(j=0; j <=i+nint; j++)
    {
      xxtmp += abs2array[j]; 
    }
    xxtmp = xxtmp/nelements;
    smoothing[i] = xxtmp;
  }
  for(i=nint; i< nfft/2-nint; i++)
  {
    xxtmp = 0;
    nelements = 2*nint+1;
    for(j = i-nint; j <= i+nint; j++)
    {
      xxtmp += abs2array[j];
    }
    xxtmp = xxtmp/nelements;
    smoothing[i] = xxtmp;
  }
  for(i = nfft/2-nint; i <= nfft/2; i++)
  {
    xxtmp = 0;
    nelements = nint + (nfft/2-i)+1;
    for(j=i - nint; j <= nfft/2; j++)
    {
      xxtmp += abs2array[j];
    }
    xxtmp = xxtmp/nelements;
    smoothing[i] = xxtmp;
  }
  for(i= nfft/2+1; i < nfft/2+ nint; i++)
  {
    xxtmp = 0;
    nelements = nint +i - nfft/2 +1;
    for(j = nfft/2; j <= nint+i; j++)
    {
      xxtmp += abs2array[j];
    }
    xxtmp = xxtmp/nelements;
    smoothing[i] = xxtmp;
  }
  for(i = nfft/2+ nint; i < nfft - nint; i++)
  {
    xxtmp = 0;
    nelements = 2*nint +1;
    for(j= i - nint; j <= i+nint; j++)
    {
      xxtmp += abs2array[j];
    }
    xxtmp = xxtmp/nelements;
    smoothing[i] = xxtmp;
  }
  for(i = nfft-nint; i < nfft; i++)
  {
    xxtmp = 0;
    nelements = nfft -i + nint ;
    for(j = i - nint; j < nfft; j++)
    {
      xxtmp += abs2array[j];
    }
    xxtmp = xxtmp/nelements;
    smoothing[i] = xxtmp;
  }

  for(i=0; i < nfft; i++)
  {
    GSL_SET_COMPLEX(&xtmp, abs2array[i], 0.);
    /* now smoothing!! */
    xarr[i] = gsl_complex_div_real(xtmp, smoothing[i]);
  }

  free(smoothing);
  free(abs2array);
  return;
}


