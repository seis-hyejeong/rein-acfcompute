#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chris_only.h"
#include "calculation_w.h"
#include "readwrite_w.h"
#include "complex_pack.h"
#include "matrix_pack.h"
#include "acf_functions_compprint.h"
#include "fft_pack.h"
#include "sac_hk.h"

char fname_param[strlng], fname_model[strlng], fname_ray[strlng], fname_out[strlng], fname_r[strlng], fname_t[strlng], fname_z[strlng], fname_tmp[strlng];
char tmp_phline[maxchar];
FILE *file_trace, *file_raygeom, *file_model;

/* prototype of the subroutines */
int errmsg(char **);

int main(int argc, char **argv)
{
  int iverbose, iinput, isoflag[maxlay], iwater;
  int i, j, k, ii, jj, errn, nray, nlayer, nfft, nreceiverlayer;
  double snr, wsmooth, dw, dt, tlength, t0, gaussf, receiverD, t_first=0;
  double m_rho[maxlay], m_thick[maxlay], m_chris[81*maxlay], rot_chris[81*maxlay], hor_rotm[9];
  double i_baz[maxray], i_rayp[maxray], q0[maxray*3];
  double *dxfft0;
  gsl_complex tmpE[36*maxlay], tmpinvE[36*maxlay], Jmatrix[6*6], Jmatrix_w[8*8], A[36*maxlay], surftmp[6], surftmp2[6]; //*x1_response, *x2_response, *x3_response, surftmp[3];
  gsl_complex tmpE_w[4], tmpinvE_w[4], A_w[4]; 
  double q3sol[6*maxlay], dvec[3*6], q0tmp[3], water_alpha;
  float *xtrace, *numer1, *denom1, *numer2, *denom2, *numer3, *denom3;
  SACHEAD hdrsac;

  /* initialize all non-array variables as possible */
  iwater = 0; water_alpha = 0; 
  receiverD = 0; nreceiverlayer = 0; /* the index that we are above of */
  dt=0.; tlength=0.; t0=0.; gaussf=-1;
  iverbose=0; iinput=0;
  nray=0; nlayer=0; errn=0;  
  wsmooth = -1;
  snr = -1;

  sprintf(fname_tmp, "FILE NAME");
  if(argc==2)
  {
    switch(argv[1][0])
    {
      case 'h':
        errmsg_general(argv[0]);
        break;
                
      case 'm':
      case 'M':
        errmsg_model(fname_tmp);
        break;
                
      case 'r':
      case 'R':
      case 'G':
      case 'g':
        errmsg_ray(fname_tmp);
        break;
      }
      exit(-1);
  }
    
  if(argc < 10)
  {
    errmsg(argv);  exit(-1);
  }
    
  i=1;
  while(i < argc)
  {
    if(argv[i][0]=='-')
    {
      switch(argv[i][1])
      {
        case 'm':
          i++;
          strcpy(fname_model, *(argv+i));
          nlayer = read_model3(fname_model, m_thick, m_rho, isoflag, m_chris, &iwater, &water_alpha);
          /* if the iwater ==1, the water velocity should not be zero */
	  if(iwater==1 && water_alpha <1.E-4){ fprintf(stderr, "the water P velocity should not be zero. Now: %f\n", water_alpha); exit(-1); }
          break;
                    
        case 'r':
          i++;
          strcpy(fname_ray, *(argv+i));
          nray = read_raygeom(fname_ray, i_rayp, i_baz);
          break;
                    
        case 'p':
          i++;
          strcpy(fname_param, *(argv+i));
          read_params(fname_param, &iinput, &dt, &tlength, &t0, &gaussf);
          break;
                
	case 'd':
          i++;
	  sscanf(&argv[i][0], "%lf", &receiverD);
 	  break;
 
        case 'o':
          i++;
          strcpy(fname_out, *(argv+i));
          break;

        case 'w':
          i++;
          sscanf(&argv[i][0], "%lf", &wsmooth);
          break;
                    
        case 's':
          i++;
          sscanf(&argv[i][0], "%lf", &snr);
          break;

        case 'v':
          iverbose=1;
          break;
                    
        default:
          errn=1;
          fprintf(stderr, "you put in un-recognizable option\n");
          break;
      }
    }
    else{ fprintf(stderr, "there is problem in your alogorithm reading inputs\n"); exit(-1); }
    i++;
  }
  if(errn!=0){ fprintf(stderr, " due to error, we cannot continue...\n"); exit(-1);}
     
  /* print out read files */
  if(iverbose==1)
  {
    /* parameter file */
    fprintf(stderr, "\n paramter setups\n\
input phase: %d (1 is P 2 is SV)\n\
sample time intv (dt): %1.3f\n\
length of signal (default): %5.2f\n\
gaussian filter (in frequency domain..): %2.3f\n\
coordinate space: RTZ without freedom of selection\n\
\n", iinput, dt, tlength, gaussf);
    /* velocity */
    fprintf(stderr, "model saved from file: %s directly as christoffel matrix\n\
and the layer number is %d\n", fname_model, nlayer);
 
    /* ray paramter file */
    fprintf(stderr, "\n total %d number of rays made to propagate and the receiver locates at %.2f \n", nray, receiverD); 

    for(i=0; i < nray; i++)
    {
      fprintf(stderr, "%02d ray: %10.4e baz: %3.0f\n", i+1, *(i_rayp+i), *(i_baz+i));
    }
  }

  /* print that we are considering noise when snr > 0  and print error when snr <= 1 --> impossible. */
  if(snr > 1 ){ fprintf(stderr, "we are considering noise with the SNR %.6f\n", snr);  }
  else if(snr > 0){ fprintf(stderr, "SNR <=1 is impossible. TERMINATING.\n"); exit(-1);  }
  
  /* determine whether it is a water layer or not */

  /* get the receiver index */
  dw = receiverD; /* use dw tempoerarily */
  while( nreceiverlayer < nlayer && dw - m_thick[nreceiverlayer] > -1.E-5)
  {
    dw -= m_thick[nreceiverlayer]; 
    nreceiverlayer++;
  }
  if(iwater ==1 && nreceiverlayer <1){ fprintf(stderr, "the receiver should not exist on the water layer\n"); exit(-1);}
  if(nreceiverlayer == nlayer){ fprintf(stderr, "the receiver should not exist too deep. \n"); exit(-1); }
  if(dw < -1.E-3){ fprintf(stderr, "the receiver at %f is not at the boundary of the layers?\n", receiverD); exit(-1); }

  /* make input waveform */
  nfft = trace_base(dt, t0+tlength);
  dxfft0 = (double *)calloc(sizeof(double *), nfft);
  gsl_complex xfft0[nfft], x1_response[ nfft], x2_response[nfft], x3_response[nfft];
  xtrace = (float *)calloc(sizeof(float *), nfft*3);
  t0 = (nfft*dt)/2; /* it should be symmetrical */

  denom1 = calloc(sizeof(float), nfft);
  numer1 = calloc(sizeof(float), nfft);
  denom2 = calloc(sizeof(float), nfft);
  numer2 = calloc(sizeof(float), nfft);
  denom3 = calloc(sizeof(float), nfft);
  numer3 = calloc(sizeof(float), nfft);

  /* make the source trace */
  for(i=0; i < nfft; i++)
  {
    GSL_SET_COMPLEX(xfft0+i, 1., 0.);
  }

  /* make list of frequencies ( angular frequency ) */
  dw = (M_PI*2/dt)/nfft;

  /* save some information to the sac header */
  hdrsac = sac_null;
  hdrsac.npts = nfft;
  hdrsac.t1 = 0; hdrsac.a=0;
  hdrsac.b = -t0;
  hdrsac.delta = (float) dt;
  if(iinput==1){ sprintf(hdrsac.ka, "P");}
  else if(iinput==2){ sprintf(hdrsac.ka, "S"); }
  hdrsac.iftype = ITIME;
  hdrsac.iztype = IB;
  hdrsac.baz=-12345;
  hdrsac.user0=-12345;
  hdrsac.cmpaz=-12345;
  hdrsac.cmpinc=-12345;

  /* onto main calculation */
  for(i=0; i < nray; i++)
  {
    hdrsac.b = -t0;
    hdrsac.delta = (float) dt;
    hdrsac.npts = nfft;
    if(iverbose==1){ fprintf(stderr, "working on %d th ray\n", i+1); }
    errn=0;
    /* set horizontal slownesses (into 2D) */
    *(q0+i*3+0) = (-1.)*cos(*(i_baz+i)*deg2rad)**(i_rayp+i);
    *(q0+i*3+1) = (-1.)*sin(*(i_baz+i)*deg2rad)**(i_rayp+i);

    mod_elastic(nlayer, q0+3*i, m_chris, hor_rotm, rot_chris); /* rotate all tensors and also get matrix used for rotation */
    /* note that horizontal ray parameters do not change within the flat-layered system */
    mat_rot(hor_rotm, q0+3*i, q0tmp); /* in order to prevent some sign difference.. (not sure) */
    t_first=0.;

    /* copy xfft0 into three components */
    trace_copy(nfft, xfft0, x1_response);
    trace_copy(nfft, xfft0, x2_response);
    trace_copy(nfft, xfft0, x3_response);

    /* first calculate the ones that are independent to frequency: except for the water layer */
    if(iwater==1)
    { 
      for(j=0; j < 6; j++){ q3sol[j] = 0.; }
      for(j=0; j < 18; j++){ dvec[j] = 0.;} 
      q3sol[0] = -sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] ); q3sol[3] = sqrt( 1./water_alpha/water_alpha - q0tmp[0]*q0tmp[0] ); 
      dvec[3*0+0] = q0tmp[0]*water_alpha; dvec[3*0+2] = q3sol[0]*water_alpha;
      dvec[3*3+0] = q0tmp[0]*water_alpha; dvec[3*3+2] = q3sol[3]*water_alpha; 
     
      /* get the Ematrix 6 x 6 */
      E_matrix(rot_chris, *(q0tmp+0), q3sol, dvec, tmpE);
      /* extract to 2 x 2 matrix and get Einv */
      tmpE_w[2*0+0] = tmpE[6*2+0]; tmpE_w[2*1+0] = tmpE[6*5+0]; /* for n = 1 */
      tmpE_w[2*0+1] = tmpE[6*2+3]; tmpE_w[2*1+1] = tmpE[6*5+3]; /* for n = 4 */
      Einv_w_matrix(tmpE_w, tmpinvE_w);
    } 

    for(j=iwater; j < nlayer; j++) /* get matrix A for each layers */
    {
      errn = qsols2(*(q0tmp+0), rot_chris+81*j, *(m_rho+j), q3sol+6*j);
      if(errn!=0){ fprintf(stderr, "complex root coming out!!!\n"); goto M_00001;}
      errn = cal_oscillvec_eigen(rot_chris+81*j, *(m_rho+j), *(q0tmp+0), q3sol+6*j, dvec);
      if(errn !=0){ fprintf(stderr, "problem in calculating oscillation vector!!!\n"); goto M_00001; }
      if(j!=nlayer-1 && j >= nreceiverlayer){      t_first += fabs(*(q3sol+6*j+iinput-1)**(m_thick+j)); }
      E_matrix(rot_chris+81*j, *(q0tmp+0), q3sol+6*j, dvec, tmpE+j*36);
      Einv_matrix(tmpE+j*36, tmpinvE+j*36);
    }

    /* get matrix E for n-th layer (j=nlayer-1) */
    /* now frequency dependent part.. 
     * calculations will be quite large due to calculation with respect to
     * frequency vale... */
    for(k=0; k < nfft; k++)
    {
      double tmpw = 0;
      if(k<nfft/2+1){ tmpw = dw*k; } /* index order following gsl_fft */
      else{ tmpw = (k-nfft)*dw; }

      for(j=iwater; j < nlayer-1; j++) /* get matrix A for each layers */
      {
        A_matrix(tmpE+j*36, tmpinvE+j*36, q3sol+6*j, tmpw, *(m_thick+j), A+36*j); 
      }
      if(iwater==1)
      { /* use the Amatrix function specially made for water layer case */
        A_w_matrix(tmpE, tmpinvE_w, q3sol, tmpw, m_thick[0], A_w);
      }
 
      /* now calculate matrix J */
      Jw_matrix(nlayer-iwater,  tmpinvE+(nlayer-1)*36, A+36*iwater, Jmatrix);
      
      /* response at the top of the solid layer */
      if(iwater==0){  get_response(iinput, Jmatrix, surftmp); }
      else /* if iwater ==1 */
      {
        //fprintf(stderr, "I pass here\n"); 
        /* Jmatrix_w is a matrix to solve (u11, u12, u13, tau33, u03) */
        /* u1i are displacement at seafloor, tau33 is stress at seafloor, f(i) are coeff at the bottom of the model, u03 is at surface of water */
        /* first, initialize the Jmatrix_w */
        for(ii=0; ii <5; ii++)
        {
          for(jj=0; jj<5; jj++)
          {
            GSL_SET_COMPLEX(Jmatrix_w+5*ii+jj, 0, 0);
          }
        }
 
        /* second, make the matrix for iwater: first */
        for(ii=0; ii < 3; ii++) 
        {
          for(jj=0; jj<3; jj++){ Jmatrix_w[ii*5+jj] = Jmatrix[ii*6+jj]; }
          Jmatrix_w[ii*5+3] = Jmatrix[ii*6 + 5];
        }
      
        GSL_SET_COMPLEX(Jmatrix_w +3*5 + 2, -1, 0); /* at (4, 3) */
        GSL_SET_COMPLEX(Jmatrix_w +4*5 + 3, -1, 0); /* at (5, 4) */
        Jmatrix_w[3*5 + 4] = A_w[0]; /* put A11 to (4, 5) */
        Jmatrix_w[4*5 + 4] = A_w[2]; /* put A21 to (5, 5) */

        get_response_w(iinput, Jmatrix_w, surftmp);
      }

      /* now get response at the receiver if the receiver is located deeper than the layer */
      /* first, we need to get Jw2 */
      if(iwater < nreceiverlayer)
      { 
        Jw_matrix2(nreceiverlayer-iwater, A+36*iwater, Jmatrix);
        cmat_vec_mul(Jmatrix, surftmp, 6, 6, surftmp2); 
        /* copy back to the surftmp */
        for(j=0; j < 6; j++){ surftmp[j] = surftmp2[j]; }
      }

      *(x1_response+k) = gsl_complex_mul(*(x1_response+k), surftmp[0]);
      *(x2_response+k) = gsl_complex_mul(*(x2_response+k), surftmp[1]);
      *(x3_response+k) = gsl_complex_mul(*(x3_response+k), surftmp[2]);
    }

    /* we want to calculate the autocorrelation of the waveforms!! */
    
    if(wsmooth > 0)
    {
      // autocorr_smooth_t2019(int nfft, double deltat, double smooth, gsl_complex *xarr)
      autocorr_smooth_t2019_noise_comp(nfft, dt, wsmooth, x1_response, snr, denom1, numer1);
      autocorr_smooth_t2019_noise_comp(nfft, dt, wsmooth, x2_response, snr, denom2, numer2);
      autocorr_smooth_t2019_noise_comp(nfft, dt, wsmooth, x3_response, snr, denom3, numer3);
    }else /* when there is no smoothing */
    {
      autocorr_noise_comp(nfft, x1_response, snr, denom1, numer1);
      autocorr_noise_comp(nfft, x2_response, snr, denom2, numer2);
      autocorr_noise_comp(nfft, x3_response, snr, denom3, numer3);
    }
   
    /* give gaussian low pass */
    if(gaussf >0)
    {
      double gaussnorm = 0;
      gaussnorm = gaussian(gaussf, nfft, dt, x1_response, x2_response, x3_response);
      if(gaussnorm <= 1.E-5){ fprintf(stderr, "error in normalization\n"); exit(-1); }
//      else{ fprintf(stderr, "gaussf: %e\n", gaussnorm); }
    }

    /* shift the time so the peak goes to zero */
    //fprintf(stderr, "t_first is %f\n", t_first);
    trace_shift_amp(t0, dt, nfft, 1, x1_response);
    trace_shift_amp(t0, dt, nfft, 1, x2_response);
    trace_shift_amp(t0, dt, nfft, 1, x3_response);   

    /* inverse fft ----> save as sac file.*/
    trace_inv_fft(x1_response, nfft, dxfft0);
    double norm = 1./dxfft0[nfft/2];
    for(k =0; k < nfft; k++){ *(xtrace+k) = (float) *(dxfft0+k)*norm; }
    trace_inv_fft(x2_response, nfft, dxfft0);
    norm = 1./dxfft0[nfft/2];
    for(k=0; k < nfft; k++){ *(xtrace+nfft+k) = (float) *(dxfft0+k)*norm; }
    trace_inv_fft(x3_response, nfft, dxfft0);
    norm = 1./dxfft0[nfft/2];
    for(k=0; k < nfft; k++){ *(xtrace+nfft*2+k) = (float) *(dxfft0+k)*norm; }
    sprintf(fname_r, "%s.R.%03d.acf", fname_out, i); /* put the index in the name */ 
    sprintf(fname_t, "%s.T.%03d.acf", fname_out, i);
    sprintf(fname_z, "%s.Z.%03d.acf", fname_out, i);

    if(iverbose==1){ fprintf(stderr, "We will save as R %s T %s and Z %s\n", fname_r, fname_t, fname_z);}

    /* give the slowness as user6 */
    hdrsac.user6 = (float) i_rayp[i];
    /* R file */ hdrsac.cmpaz = (float) *(i_baz+i); hdrsac.cmpinc =90;
    wsac( fname_r, hdrsac, xtrace);
    /* T file */ hdrsac.cmpaz = (float) *(i_baz+i)+90;
    wsac( fname_t, hdrsac, xtrace+nfft);
    /* Z file */ hdrsac.cmpaz = 0; hdrsac.cmpinc=0;
    wsac( fname_z, hdrsac, xtrace+2*nfft);

    /* now write the denominator and the numerator */
    hdrsac.b = 0; hdrsac.npts = 0.5*nfft;
    hdrsac.delta = 1./(dt*nfft);

    sprintf(fname_r, "%s.R.%03d.denom", fname_out, i);
    sprintf(fname_t, "%s.T.%03d.denom", fname_out, i);
    sprintf(fname_z, "%s.Z.%03d.denom", fname_out, i);
    wsac(fname_r, hdrsac, denom1);
    wsac(fname_t, hdrsac, denom2);
    wsac(fname_z, hdrsac, denom3);

    /* now save the numerator */
    sprintf(fname_r, "%s.R.%03d.numer", fname_out, i);
    sprintf(fname_t, "%s.T.%03d.numer", fname_out, i);
    sprintf(fname_z, "%s.Z.%03d.numer", fname_out, i);
    wsac(fname_r, hdrsac, numer1);
    wsac(fname_t, hdrsac, numer2);
    wsac(fname_z, hdrsac, numer3);

    if(iverbose==1){ fprintf(stderr, "successfully wrote for baz %4.1f file as %s\n", *(i_baz+i), fname_z); }

    M_00001: ;
    if(errn!=0)
    {
      fprintf(stderr, "too large horizontal ray parameter...!\n");
    }
   /* convert it to time series : note that what caused error will become zero trace */
  }
      
  return 0;
}

/* subroutines */
int errmsg(char **input)
{
  fprintf(stderr, " ACF with smoothing \n\
  This code started as replacement of code anirec (Levin and Park) \n\
  THE TOP LAYER CAN BE A WATER LAYER, THE RECEIVER DP SHOULD BE DEEPER THAN THE WATER LAYER. \n\
  For calculating the Autocorrelation Function, this takes the definition of Tauzin. \n\
  File name of parameter files and such should be given and input and\n\
  output options or its name should be also set as option\n\n\
  [execution format] \n\
    %s [-d (receiverdp) -w (smoothing width) -s (snr to have noise)] -m (model file name) -p (parameter file) -r (ray geometry information file) -o (out filename) [-v(erbose)]\n\n\
  make sure that file contents are in exact format. \n\
  If you need information about the format of the file, please type 'h' \n\
  ENJOY! (February 2022 HK)\n\n", input[0]);

  return 0;
}
