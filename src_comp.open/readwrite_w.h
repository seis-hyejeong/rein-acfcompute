#ifndef _readwrite_w_h_
#define _readwrite_w_h_

#include "parameters.h"
#include "sac_hk.h"
int read_params(char [], int *, double *, double *,double *, double *);
int read_model3(char [], double *, double *, int *, double *, int *, double *);
int read_raygeom(char [], double *, double *);
void write_as_sac(int, SACHEAD, float, float, float, float);
void errmsg_general(char []);
void errmsg_param(char []);
void errmsg_model(char []);
void errmsg_ray(char []);

#endif
