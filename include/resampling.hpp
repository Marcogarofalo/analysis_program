#ifndef resampling_H
#define resampling_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>



void free_jack(int N,int var , int t, double ****in);
void write_jack(int N, double *jack, char *outname);
void write_jack_bin(int N, double *jack, char *outname);

double ***read_jack(int N, int var, int T);
void write_jack_corr(int N, int t,double **jack, char *outname);
//create_jack
//in[#conf.][#variable][#time_cordinate][#re or im]
//returns the jacknife configuration from the data ****in
//the last entry of [#conf] is the average
double ****create_jack(int  N, int var, int t, double ****in);
double *mean_jack(int N,int var,int t, int call, double ****jack, double function_jack(int,int,int,double ***) );
//mean_and_error_jack
//returns the mean and error from set of N  jacknife called *in  and the average stored in in[N]
double *mean_and_error_jack(int Np1, double *in);
double *mean_and_error_jack_biased(int Np1, double *in);
double *mean_and_error_jack_biased1(int Np1, double *in);

double *fake_jack(double mean,double error, int Njack,int seed);




/////////////boot
double ****create_boot(int  N, int Nboot, int var, int t, double ****in);
double *mean_and_error_boot(int Np1, double *in);
double *fake_boot(double mean,double error, int Njack,int seed);


/////////////////////
double *mean_and_error( const char *option , int Np1, double *in);  
double error_jackboot( const char *option , int Np1, double *in);  

const char  *smean_and_error( const char *option , int Np1, double *in);  
double ****create_resampling(const char *option, int  N, int var, int t, double ****in,int seed=123);
double *fake_sampling(const char *option,double mean,double error, int Njack,int seed);

double **covariance(const char *option , int Nobs, int Np1, double **in);
double **error_covariance(const char *option , int Nobs, int Np1, double **in);




double **fake_sampling_covariance(const char *option,double *mean, int Njack,int N, double **cov,int seed);


///////////////////////////////////////////////////////
double* malloc_copy_jackboot(int Np1,  double *a);
    
void sum_jackboot(int Np1,  double *r, double *a, double *b);
void sub_jackboot(int Np1,  double *r, double *a, double *b);
void mult_jackboot(int Np1,  double *r, double *a, double *b);
void div_jackboot(int Np1,  double *r, double *a, double *b);
void invert_jackboot(int Np1,  double *r, double *a);

void scalar_times_jackboot(int Np1,  double *r, double *a, double s);


#endif
 
