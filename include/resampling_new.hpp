#ifndef resampling_new_H
#define resampling_new_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <iostream>
#include "global.hpp"




class resampling_f {
public:
  int Njack;
  char option[NAMESIZE];
  int seed;
  // resampling_f() {
  //   printf("Error you need to initialise the resampling procedure");
  //   exit(2);
  // }
  // resampling_f(const char* name, int N);
  void free_jack(int var, int t, double**** in);
  void write_jack_bin(double* jack, char* outname);

  virtual double**** create(int  N, int var, int t, double**** in)=0 ;
  // double* mean_and_error(  double* in);
  virtual double comp_error(double* in)=0 ;
  virtual double* create_fake(double mean, double error, int seed)=0 ;
  // double** comp_covariance(int Nobs,  double** in);

  // double** fake_sampling_covariance( double* mean,  double** cov, int seed);
  // double* malloc_copy( double* a);
  void free_res(int var, int t, double**** in);
  void write_res_bin(int N, double* jack, char* outname);

};

class resampling_jack : public resampling_f {
public:
    resampling_jack();
    resampling_jack(int N);
    double**** create(int  N, int var, int t, double**** in);
    double comp_error(double* in);
    double* create_fake(double mean, double error, int seed);
};


class resampling_boot : public resampling_f {
public:   
    resampling_boot();
    resampling_boot(int N);
    resampling_boot(int N, int o_seed);
    double**** create(int  N, int var, int t, double**** in);
    double comp_error(double* in);
    double* create_fake(double mean, double error, int seed);

};



#endif