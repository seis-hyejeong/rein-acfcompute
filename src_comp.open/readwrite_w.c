#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "matrix_pack.h"
#include "readwrite_w.h"
#include "parameters.h"

/* read information from parameter file */
int read_params(char fname[maxchar], int *iphase, double *dt, double *tl, double *t0, double *gaussianf)
{
  int errn=0;
  FILE *infile;
  char stmp[maxchar], tmp[100];
  infile = fopen(fname, "r");

  if(infile==NULL){ fprintf(stderr, "your file %s does not exist?\n", fname); exit(1); }

  while(fgets(stmp, maxchar, infile)!=NULL)
  {
    switch(stmp[0])
    {
      case '#':
        break;

      case 'I': // input phase
        sscanf(&stmp[1], "%s", tmp);
        *iphase = atoi(tmp);
        break;

      case 'T': // sample interval in sec and total time length
        sscanf(&stmp[1], "%lf/%lf", dt, tl);
        break;

      case 'O': // origin time signal (when to start your signal)
        sscanf(&stmp[1], "%s", tmp);
        *t0 = atof(tmp); // origin time
        if(*t0<5){ fprintf(stderr, "it should be larger than 5. please correct it and try again\n"); errn=1;}
        break;

      case 'F':
      case 'f': // asking for input peak's guassian width factor
        sscanf(&stmp[1], "%s", tmp);
        *gaussianf = atof(tmp);
        break;

      default:
        errn=1;
    }
  }

  if(errn==1)
  {
    errmsg_param(fname);
  }

  fclose(infile);
  return errn;
}

/* read input velocity model from file name */
int read_model3(char fname[maxchar], double *thick, double *rho, int *iso, double *chris, int *iwater, double *water_alpha)
{
  int i, nlayer=0;
  double tmpalpha, tmpbeta, tmpb, tmpc, tmpe, s_azi, s_tilt;
    
  FILE *infile;
  char stmp[maxchar], stmp2[maxchar];
  infile = fopen(fname, "r");

  if(infile==NULL){ fprintf(stderr, "your file %s does not exist?\n", fname); exit(1);}

  i=0; *iwater = 0;
  while(fgets(stmp, maxchar, infile)!=NULL)
  {
      tmpalpha=0.; tmpbeta=0.; tmpb=0; tmpc=0; tmpe=0; s_azi=0; s_tilt=0;

      if(stmp[0]=='#' || sscanf(stmp, "#%s", stmp2)==1 || sscanf(stmp, "# %s", stmp2)==1){}
      else if(sscanf(stmp, "%lf %lf %lf %lf %d %lf %lf %lf %lf %lf", thick+i, rho+i, &tmpalpha, &tmpbeta, iso+i, &tmpb, &tmpc, &tmpe, &s_azi, &s_tilt)!=10){ nlayer=-1; goto L_00001;}
      else
      {
        if( fabs(tmpbeta) < 1.E-4 && i ==0){ *iwater = 1; *water_alpha = tmpalpha; }
        set_elastic(*(rho+i), tmpalpha, tmpbeta, tmpb, tmpc, tmpe, s_azi, s_tilt, chris+81*i);
      
        i++;
      }
  }
  nlayer = i;

  if (i==0 || i==1){ errmsg_model(fname);}

L_00001: ;
  return nlayer;
}

/* read ray input file */
int read_raygeom(char fname[maxchar], double *rayp, double *baz)
{
  int i;
  char stmp[maxchar], stmp2[maxchar], tmp1[100], tmp2[100];
  FILE *infile;

  infile = fopen(fname, "r");
  i=0;
  while(fgets(stmp, maxchar, infile)!=NULL)
  {
    if(stmp[0]=='#' || sscanf(stmp, "#%s", stmp2)==1 || sscanf(stmp, "# %s", stmp2)==1){}
    else
    {
          sscanf(stmp, "%s %s", tmp1, tmp2);
          *(baz+i) = atof(tmp1);
          *(rayp+i) = atof(tmp2);

          i++;
    }
  }

  int nrays = i;

  return nrays;
}

/* write as sac file */
void initialize_sac(int npts, SACHEAD sachdr, float dt, float b, float t0, float rayp)
{
  sachdr = sac_null;
  sachdr.npts = npts;
  sachdr.t1 = t0; sachdr.a = t0;
  sachdr.b = b;
  sachdr.delta = dt;
  sprintf(sachdr.ka, "P");
  sachdr.iftype = ITIME;
  sachdr.iztype = IB;
  sachdr.user0 = rayp;

  return;
}

/* collection of all error messages */
/* general errormessage */
void errmsg_general(char progname[maxchar])
{
    char stmp[maxchar], fname[maxchar];
    sprintf(fname, "FILENAME");

    fprintf(stderr, " \n\
This program is:\n\
PAPS \n\
This code starts from the code that replaces the anirec (Levin and Park) \n\
File name of parameter files and such should be given and input and\n\
output options or its name should be also set as option\n\n\
[execution format] \n\
   %s -mod (model file name) -param (parameter file) -raygeom (ray geometry information file) -o (out filename) [-v(erbose)]\n\n\
make sure that file contents are in exact format. \n\n", progname);

    fprintf(stderr, "do you need information of the file format?\n\
            param/model/ray(geom) to see, 0 to exit: ");

    scanf("%s", stmp);
    switch(stmp[0])
    {
      case 'P':
      case 'p':
        errmsg_param(fname);
        break;

      case 'M':
      case 'm':
        errmsg_model(fname);
        break;

      case 'R':
      case 'r':
        errmsg_ray(fname);
        break;

      case '0':
        break;
    }

    return;
}


/*error parameter file */
void errmsg_param(char fname[maxchar])
{
    fprintf(stderr, " :::::::::::::::::::::::\n\
            [Parallel layered Anisotropic Propagator matrix based Synthetics]\n\
:::: %s :::: parameter file should be in following form :\n\
please write line by line.  \n\
# : can comment anything with # in front of a line. \n\
I(input phase): 1 for P wave input 2 for SV wave input \n\
T(dt/tlength) : sampling rate and time length separated by slash  \n\
O(origin time): time to start the signal (signal starts from t= 0s anyway). \n\
F(gaussian filter number): This is stdev in time domain, not frequency domain\n\
\n", fname);

    return;
}

/* collection of error messages */
void errmsg_model(char fname[])
{

  fprintf(stderr, "[Input Velocity model of TRAIL code]\n\
Your file %s is wrong \n\
each line should be in following format in following unit. \n\
::: thick[m] rho[kg/m^3] vp[m/s] vs[m/s] isoflag P_perturb(B) P_perturb(C) S_perturb(E) trend plunge[degree] \n\
::: adding comments with # is okay (not at the same line as the input line). \n\
::: velocity perturbation is 2theta, 4theta variance and they are in percentage unit in peak-to-peak calculation \n\
::: \n" , fname);

     return;
}

/* ray geometry file */
void errmsg_ray(char fname[maxchar])
{
    fprintf(stderr, " \n\
            [ input ray geometry file ] \n\
            Your lines in %s may be wrong. \n\
            each line should be in following format: \n\
            ::: baz [degree] ray parameter [s/m] NS displacement [m] EW displacement [m] \n\
            ::: Positive displacement is to N or E and Negative displacement refers to S or W. \n\n", fname);
    
    return;
}
