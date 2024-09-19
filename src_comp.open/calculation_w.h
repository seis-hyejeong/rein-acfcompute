#ifndef _calculation_w_h_
#define _calculation_w_h_

#include "gsl/gsl_poly.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_linalg.h"

#include "parameters.h"
#include "matrix_pack.h"
#include "fft_pack.h"

void mod_elastic(int, double [], double [], double *, double *);
int qsols2(double , double [], double, double *);
void oscillvec(double [], double [], double [], double [], double [], double *, int);
void E_matrix(double [], double, double [], double [], gsl_complex *);
void Einv_matrix(gsl_complex [], gsl_complex *);
void Einv_w_matrix(gsl_complex [], gsl_complex *);
void A_matrix(gsl_complex [], gsl_complex [], double [], double, double, gsl_complex *);
void A_w_matrix(gsl_complex [], gsl_complex [], double [], double, double, gsl_complex *);
void cal_polymatrix(int, int, double *, double, double *);
double cal_poly(int, double *, double);
void Jw_matrix(int, gsl_complex [], gsl_complex [], gsl_complex *);
void Jw_matrix2(int, gsl_complex [], gsl_complex *);
void get_response(int, gsl_complex [], gsl_complex *);
void get_response_w(int, gsl_complex [], gsl_complex *);
void conv_source(gsl_complex *, gsl_complex *, gsl_complex *, gsl_complex [], int);
void source(int, double, double, gsl_complex *);
void output_invrot(int, double [], double [], double []);
double p_tt(int, double [], double []);

#endif
