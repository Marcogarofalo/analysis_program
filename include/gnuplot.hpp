#ifndef gnuplot_H
#define gnuplot_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>



void plotting(int T, double **corr,int *tmin,int *tmax  , int *sep);
void plotting_fit_p(int T, double **corr,int tmin,int tmax, double *fit );
int plotting_fit(int T, double **corr,int tmin,int tmax, double *fit ,double *chi2);

void contraction_name(int ii,char *r);

void plotting_fit_pdf( char **option,const char *string,int L, int T, double beta, double k_sea,double mu_sea, double **corr,int tmin,int tmax, double *fit, const char *contraction ,struct kinematic kinematic_2pt); 
void plotting_fit_deg_n_pdf(char **option, const char *string,int npoints,int deg,double *x, double **y, double **fit,double *chi2 );
void plotting_fit_poles_deg_n_pdf(char **argv, const char *string,int npoints , int deg ,double *x, double **y, double **fit,double *chi2 );
void contraction_index(int *ii ,const char*r );

void plotting_fit_pdf_G(char **argv, const char *string ,double **corr,int tmin,int tmax, double *fit,const char *contraction , struct kinematic_G kinematic_2pt_G );


#endif
