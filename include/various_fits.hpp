#ifndef various_fit_H
#define various_fit_H
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
double *polynimial_degree_n(int np1, double in);
double *poles_degree_n(int np1, double in);

double  **fit_polynomial(char **argv,const char* string, int n, int kmin,int ik1, char* name, double ***mass_jack_fit_GEVP,int  Njack , FILE *outfile );


double **fit_FX(char **argv, const char  *string,int ikt,int iks, double **kp_tot, double **FX,  struct header file_head ,int Njack ,int npar_fun,double *fit_function(int,double) );
double **interpolate_FX(char **argv, const char  *string,int ikt,int iks,double **kp_tot,double xmin,double xmax, double **FX,  struct header head ,int Njack ,int npar_fun,double *fit_function(int,double) );


#endif
