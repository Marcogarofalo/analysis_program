#ifndef read_H
#define read_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>


double ****calloc_corr(int N, int var, int t );
void copy_corr(int N, int var, int t, double ****in, double ****out);
void free_corr(int N, int var, int t, double ****out);
double ****binning(int N, int var, int t ,double ****data,int bin);
double ****binning_toNb(int N, int var, int t ,double ****data,int Nb);
double**** bin_intoN( double**** data, int nvar, int T, int Nconf_in, int Nb);
//data[#conf.][#variable][#time_cordinate][#re or im]
double ****read_datafile(int N, int var, int t, int bin );
void symmetrise_corr(int N, int var, int t ,double ****data);
void antisymmetrise_corr(int N, int var, int t ,double ****data);
void forward_derivative_corr(int N, int var, int t ,double ****data);
void symmetric_derivative_corr(int N, int var, int t ,double ****data);
void second_derivative_corr(int N, int var, int t ,double ****data);
void cbtc_corr(int N, int var, int t ,double ****data);
#endif    
