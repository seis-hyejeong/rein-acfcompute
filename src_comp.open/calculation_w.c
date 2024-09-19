#include <math.h>

#include "calculation_w.h"
#include "matrix_pack.h"
#include "parameters.h"
 
/* This basically includes making matrices necessary for calculation */
/* notations follow Crampin 1977 */

/* prepare all christoffel tensors */
void mod_elastic(int nlayer, double qvec[], double i_tensor[], double *rotM, double *o_tensor)
{
  int i;
  double hangle=0;

  hangle = get_rothor_angle(qvec);
  get_2d_rotmat(hangle, rotM);

  for(i=0; i < nlayer; i++)
  {
    rot_tensor(rotM, i_tensor+i*tintv, o_tensor+i*tintv);
  }
  
  return;
}

/* get phase velocities for given ray paramter */
int qsols2(double rayp, double C[], double rho, double *q3sol)
{
  int errn =0; /* if return value errn is 1, it means this step is not succesful */
  double vq[3][2], F[9][3], poly[9], qsoltmp[6*2], tmp[3], tmp1[7], tmp2[7], det_tmp[3][7];
  int i, j, k, m, n;
  double tmpC=0;
  /* parameter space for solving polynomials */
  gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc (7);

  for(i=0; i <3; i++){ vq[i][0]=0; vq[i][1]=0; }
  vq[2][1] = 1;  vq[0][0] = rayp;

  /* obtain matrix F (equation (5) or Keith and Crampin 1977 a) */
  for(j=0; j<3; j++)
  {
    for(m=0; m <3; m++)
    { /* first, initialize the terms */
      F[j*3+m][0] =0; F[j*3+m][1] =0; F[j*3+m][2] =0;
      tmp[0]=0; tmp[1]=0; tmp[2]=0;
      
       for(k=0; k<3; k++)
       {
         for(n=0; n<3; n++)
         {
           tmpC = tensor_get(C,j+1, k+1, m+1, n+1);
           conv(vq[k], 2, vq[n], 2, tmp);
           amul_v(3, tmp, tmpC);
           sum_v(3, tmp, F[j*3+m]);
         }
       }
       if(j==m){ F[j*3+m][0] -= rho; }
    }
  }
  
  /* sub11 */
  conv(F[(2-1)*3+(2-1)], 3, F[(3-1)*3+(3-1)], 3, tmp1);
  conv(F[(2-1)*3+(3-1)], 3, F[(3-1)*3+(2-1)], 3, tmp2); amul_v(5, tmp2, -1);
  sum_v(5, tmp1, tmp2); // tmp2 = tmp1+tmp2;
  conv(F[(1-1)*3+(1-1)], 3, tmp2, 5, det_tmp[0]); /* x11*sub11 */

  /* sub12 */
  conv(F[(2-1)*3+(1-1)], 3, F[(3-1)*3+(3-1)], 3, tmp1); 
  conv(F[(2-1)*3+(3-1)], 3, F[(3-1)*3+(1-1)], 3, tmp2); amul_v(5, tmp2, -1);
  sum_v(5, tmp1, tmp2); // tmp2 = tmp1+tmp2;
  conv(F[(1-1)*3+(2-1)], 3, tmp2, 5, det_tmp[1]);

  /* sub13 */
  conv(F[(2-1)*3+(1-1)], 3, F[(3-1)*3+(2-1)], 3, tmp1); 
  conv(F[(2-1)*3+(2-1)], 3, F[(3-1)*3+(1-1)], 3, tmp2); amul_v(5, tmp2, -1);
  sum_v(5, tmp1, tmp2); 
  conv(F[(1-1)*3+(3-1)], 3, tmp2, 5, det_tmp[2]);

  /* det(F) = sub11 - sub12 + sub13 */
  amul_v(7, det_tmp[1], -1);
  sum_v(7, det_tmp[0], det_tmp[1]);
  sum_v(7, det_tmp[1], det_tmp[2]);
  for(i=0; i < 7; i++){ poly[i] = det_tmp[2][i]; }

  /* solve the polynomial */
  gsl_poly_complex_solve( poly, 7, w, qsoltmp); /* qsoltmp is complex solutions */
  gsl_poly_complex_workspace_free(w);

  /* handle obtained solutions */
  for(i=0; i <6; i++)
  {
    if(qsoltmp[2*i+1]>1E-4 || fabs( qsoltmp[2*i+0]) <1E-5){ errn =1;}
  }
  if(errn==1){ goto L_0001;}

  gsl_sort(qsoltmp, 2, 6); /* gives it in ascending order (increasing order) */

  /* first: P wave  --> goes to index 1, 4 */
  *(q3sol+0) = qsoltmp[2*2];
  *(q3sol+3) = qsoltmp[3*2];

  /* second: S1 wave --> goes to index 2, 5 */
  *(q3sol+1) = qsoltmp[1*2];
  *(q3sol+4) = qsoltmp[4*2];

  /* third: S2 wave  --> goes to index 3, 6 */
  *(q3sol+2) = qsoltmp[0*2];
  *(q3sol+5) = qsoltmp[5*2];

  L_0001: ;
  return errn;
}

/* for propagator part of matrices */ 

/* matrix E in Keith and Crampin 1977 */
/* output is matrix of complex numbers for convenience of overall calculation */
void E_matrix(double C[], double rayp, double qz[], double dvec[], gsl_complex *Emat)
{ 
  int j, n, m;
  double c, tmp=0;
  c = 1./rayp; /* the expression in Keith & Crampin 1977 use phase velocity! */

  for(n=0; n <6; n++)
  {
    for(j=0; j<3; j++){  GSL_SET_COMPLEX(Emat+j*6+n, dvec[3*n+j], 0.); }
  }

  for(n=0; n <6; n++)
  {
    for(j=0; j<3; j++)
    {  
      tmp=0;
      for(m=0; m <3; m++)
      {
        tmp += dvec[3*n+m]*(tensor_get(C, j+1, 3, m+1, 1)+ tensor_get(C, j+1, 3, m+1, 3)*c*qz[n]);
      }
      GSL_SET_COMPLEX(Emat+(j+3)*6+n, tmp, 0.); 
    }
  }

  return;
}

/* the water Einv_matrix (difference is the dimension of the matrix */
void Einv_w_matrix(gsl_complex Emat[], gsl_complex *Einvmat)
{
  int i, j, s;
  double tmp= 0;
 
  gsl_permutation * p = gsl_permutation_alloc(2);
  gsl_matrix * Ereal = gsl_matrix_calloc(2,2);
  gsl_matrix * invE = gsl_matrix_calloc(2,2);

  /* make it as real matrix to invert it ! */
  for(i=0; i < 2; i++)
  {
    for(j=0; j < 2; j++)
    {
      gsl_matrix_set(Ereal, i, j, GSL_REAL(Emat[2*i+j]));
    }
  }

  gsl_linalg_LU_decomp(Ereal, p, &s);
  gsl_linalg_LU_invert(Ereal, p, invE);

  for(i=0; i < 2; i++)
  {
    for(j=0; j < 2; j++)
    {
      tmp = gsl_matrix_get(invE, i, j); 
      GSL_SET_COMPLEX(Einvmat+2*i+j, tmp, 0);
    } 
  }

  gsl_matrix_free(Ereal); gsl_matrix_free(invE);
  gsl_permutation_free(p);

  return;
}
 
void Einv_matrix(gsl_complex Emat[], gsl_complex *Einvmat)
{    
  int i, j, s;
  double tmp=0;

  gsl_permutation * p = gsl_permutation_alloc(6);
  gsl_matrix * Ereal = gsl_matrix_calloc(6,6);
  gsl_matrix * invE = gsl_matrix_calloc(6,6);

  /* make it as real matrix to invert it! */
  for(i=0; i< 6; i++)
  {    
    for(j=0; j<6; j++)
    {
      gsl_matrix_set(Ereal, i, j, GSL_REAL(Emat[6*i+j]));
    }
  }

  gsl_linalg_LU_decomp(Ereal, p, &s);
  gsl_linalg_LU_invert(Ereal, p, invE);

  for(i=0; i <6; i++)
  {
    for(j=0; j<6; j++)
    {
      tmp = gsl_matrix_get(invE, i, j);
      GSL_SET_COMPLEX(Einvmat+6*i+j, tmp, 0);
    }
  }

  gsl_matrix_free(Ereal); gsl_matrix_free(invE);
  gsl_permutation_free(p);

  return;
}

/* get A matrix in water */
void A_w_matrix(gsl_complex Emat[], gsl_complex E_w_invmat[], double qz[], double w, double dm, gsl_complex *A)
{
  int j, n; 
  gsl_complex Dw[36], Dw_w[4], expw;
  
  for(n=0; n<6; n++)
  {
    expw = gsl_complex_polar(1., -w*dm*qz[n]);
    for(j=0; j <6; j++)
    {
      *(Dw+j*6+n) = cmul(expw, Emat[j*6+n]); /* D = expmat*E */ 
    }
  }

  /* get the Dw_w from the whole Dw */
  Dw_w[2*0+0] = Dw[6*2+0]; Dw_w[2*1+0] = Dw[6*5+0]; /* for n = 1 */
  Dw_w[2*0+1] = Dw[6*2+3]; Dw_w[2*1+1] = Dw[6*5+3]; /* for n = 4 */ 

  /* A = D*invE */
  cmat_mat_mul(Dw_w, E_w_invmat, 2, 2, 2, A);  

  return;
}

/* A = D*(invE) */ /* follow exactly Keith and Crampin 1970 */
void A_matrix(gsl_complex Emat[], gsl_complex Einvmat[], double qz[], double w, double dm, gsl_complex *A)
{
  int j, n;
  gsl_complex Dw[36], expw;

  for(n=0; n<6; n++)
  {
    expw = gsl_complex_polar(1.,-w*dm*qz[n]);
    for(j=0; j<6; j++)
    {
      *(Dw+j*6+n) = cmul(expw, Emat[j*6+n]); /* invD = iExpmat*invE */
    }
  }
  
  /* A = D*invE */
  cmat_mat_mul( Dw, Einvmat, 6, 6, 6, A);
  
  return;
}

double cal_poly(int order, double *coeffs, double x)
{ // it is highest order. not number of terms (for 4th order polynomial it is 4)
    int i;
    double xx, output=0;
    xx=1;

    for(i=0; i <=order; i++)
    {
        output += *(coeffs+i)*xx;
        xx = xx*x;
    }
    return output;
}

/***** MAIN propagator part!!!!!******************/
/* fn = inv(En)*A(n-1)*A(n-2)*...A1*v0 */
/* fn = J*v0
 * v0 = inv(J)*fn --> we will get matrix J
 */ 
/* fn = (f1 f2 f3 f4 f5 f6) and v0 = (u01, u02, u03, 0, 0, 0) */
void Jw_matrix(int nlayer, gsl_complex invEn[], gsl_complex A[], gsl_complex *J)
{
    int i,j ;
    gsl_complex Jtmp[36];
    
    /* initialize matrix J */
    for(i=0; i<6; i++)
    {
        for(j=0; j<6; j++)
        {
            if(i==j){ GSL_SET_COMPLEX(J+6*i+j, 1, 0);}
            else{ GSL_SET_COMPLEX(J+6*i+j, 0, 0); }
        }
    }
    
    for(i=0; i < nlayer-1; i++) /* add at left */
    {
        cmat_mat_mul(A+i*36, J, 6, 6, 6, Jtmp);
        for(j=0; j<36; j++){ J[j] = Jtmp[j]; }
    }
    
    /*multiply invEn at left */
    cmat_mat_mul(invEn, J, 6, 6, 6, Jtmp);
    for(j=0; j<36; j++){ J[j] = Jtmp[j]; }
    
    return;
}

void Jw_matrix2(int nlayer, gsl_complex A[], gsl_complex *J)
{
  int i, j;
  gsl_complex Jtmp[36];

  /* initialize matrix J */
  for(i=0; i<6; i++)
  {
    for(j=0; j < 6; j++)
    {
      if(i==j){ GSL_SET_COMPLEX(J+6*i+j, 1, 0); }
      else{ GSL_SET_COMPLEX(J+6*i+j, 0, 0); }
    }
  }

  for(i=0; i < nlayer; i++) /* add at left, we will multiply all nlayers */
  {
    cmat_mat_mul(A+i*36, J, 6, 6, 6, Jtmp);
    for(j=0; j<36; j++){ J[j] = Jtmp[j]; }
  }

  /* we do not need to add invE at the left*/
  return;
}

/* calculate for selected input */
/* we can have either P input or SV input */
/* solve matrix equation...! */
/* iinci =1 -> P wave incidence, iinci =2 -> S wave */
void get_response(int iinci, gsl_complex Jwmat[], gsl_complex *surf_u)
{ 
  int i,j, s;
  gsl_complex tmp;

  /* in = J11*(x) */
  gsl_matrix_complex *J11 = gsl_matrix_complex_calloc(3, 3);
  gsl_vector_complex *x = gsl_vector_complex_calloc(3);
  gsl_vector_complex *in = gsl_vector_complex_calloc(3); /* our target is to get x */
  gsl_permutation *p = gsl_permutation_alloc(3);
  
  /* initialize gsl complex array "in" */
  GSL_SET_COMPLEX(&tmp, 0, 0);
  for(i=0; i <3; i++){ gsl_vector_complex_set(in, i, tmp); }
  GSL_SET_COMPLEX(&tmp, 1., 0.);
  gsl_vector_complex_set(in, iinci-1, tmp);

  /* extract part of matrix */
  for(i=0; i<3; i++)
  {
    for(j=0; j<3; j++)
    {
      gsl_matrix_complex_set(J11, i, j, Jwmat[6*(i)+(j)]);
    }
  }

  gsl_linalg_complex_LU_decomp(J11, p, &s);
  gsl_linalg_complex_LU_solve(J11, p, in, x);

  for(i=0; i<3; i++)
  {
    tmp = gsl_vector_complex_get(x, i);
    *(surf_u+i) = tmp; 
  } 
  for(i=3; i < 6; i++)
  {
    GSL_SET_COMPLEX(surf_u+i, 0., 0.);
  }

  gsl_matrix_complex_free(J11);
  gsl_vector_complex_free(x);
  gsl_vector_complex_free(in);

  return;
}

/* get response in case of existing water layer */
void get_response_w(int iinci, gsl_complex Jwmat[], gsl_complex *surf_u)
{
  int i,j, s;
  gsl_complex tmp;

  /* in = J11*(x) */
  gsl_matrix_complex *J11 = gsl_matrix_complex_calloc(5, 5);
  gsl_vector_complex *x = gsl_vector_complex_calloc(5);
  gsl_vector_complex *in = gsl_vector_complex_calloc(5); /* our target is to get x */
  gsl_permutation *p = gsl_permutation_alloc(5);

  /* initialize gsl complex array "in" */
  GSL_SET_COMPLEX(&tmp, 0, 0);
  for(i=0; i <5; i++){ gsl_vector_complex_set(in, i, tmp); }
  GSL_SET_COMPLEX(&tmp, 1., 0.);
  gsl_vector_complex_set(in, iinci-1, tmp);

  /* copy the matrix to the one we will be using */
  for(i=0; i<5; i++)
  {
    for(j=0; j<5; j++)
    {
      gsl_matrix_complex_set(J11, i, j, Jwmat[5*(i)+(j)]);
    }
  }

  gsl_linalg_complex_LU_decomp(J11, p, &s);
  gsl_linalg_complex_LU_solve(J11, p, in, x);
  for(i=0; i<3; i++)
  {
    tmp = gsl_vector_complex_get(x, i);
    *(surf_u+i) = tmp; 
  }
  for(i=3; i < 5; i++)
  {
    GSL_SET_COMPLEX(surf_u+i, 0., 0.);
  }
  *(surf_u + 5) = gsl_vector_complex_get(x, 3); /* the tau33 */

  gsl_matrix_complex_free(J11);
  gsl_vector_complex_free(x);
  gsl_vector_complex_free(in);

  return;
}


/* get convolution of source signal and response to get output */
/* u1(R), u2(T), u3(Z) is all organized by frequency and signal is organized also in same frequency */
void conv_source(gsl_complex *u1, gsl_complex *u2, gsl_complex *u3, gsl_complex source[], int nfft)
{
  int i;
  gsl_complex tmp1, tmp2, tmp3;
  for(i=0; i < nfft; i++)
  {
    tmp1 = cmul(*(u1+i), source[i]); *(u1+i)=tmp1;
    tmp2 = cmul(*(u2+i), source[i]); *(u2+i)=tmp2;
    tmp3 = cmul(*(u3+i), source[i]); *(u3+i)=tmp3;
  }

  return;
}

void source(int nfft, double dt, double gaussf, gsl_complex *source)
{
  int i;
  double x[nfft];
  for(i=0; i < nfft; i++)
  {
      *(x+i) = exp( -pow(dt*(i-nfft/2),2)/(4*gaussf*gaussf));
  }
  trace_fft(x, nfft, source);
  trace_shift_amp(-(nfft/2)*dt, dt, nfft, 1, source);

  return;
}

/* rotate back to NEZ */
void output_invrot(int npts, double rotM [], double u1[], double u2[])
{
  int i;
  double tmpur[3] = {0, 0, 0};
  double tmpu[3] = {0, 0, 0};
  
  for(i=0; i < npts; i++)
  {
    tmpur[0] = u1[i]; tmpur[1] = u2[i];
    mat_rot_inv(rotM, tmpur, tmpu);
    u1[i] = tmpu[0]; u2[i] = tmpu[1]; 
  }
  
  return;
}

/* calculate P arrival time */
double p_tt(int nlayer, double qz_p[], double thickness[])
{
  int i;
  double tt=0;

  for(i=0; i < nlayer; i++)
  {
    tt += qz_p[i]*thickness[i];
  }

  return tt;
}
