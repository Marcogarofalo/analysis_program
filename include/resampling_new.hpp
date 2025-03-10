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

#ifdef WITH_ARB
#include "arb.h"
#include "acb_calc.h"
#endif // WITH_ARB



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
  
  void write_jack_in_file(double* jack, const char* outname);
  void read_jack_from_file(double* jack, const char* name);



  virtual double**** create(int  N, int var, int t, double**** in) = 0;
  // double* mean_and_error(  double* in);
  virtual double comp_mean_unbias(double* in) = 0;
  virtual double comp_error(double* in) = 0;
  virtual double* create_fake(double mean, double error, int seed) = 0;
  // double** comp_covariance(int Nobs,  double** in);

  // double** fake_sampling_covariance( double* mean,  double** cov, int seed);
  // double* malloc_copy( double* a);
  void free_res(int var, int t, double**** in);
  void write_res_bin(int N, double* jack, char* outname);
  virtual double** comp_cov(int Nobs, double** in) = 0;
#ifdef WITH_ARB
  virtual void comp_cov_arb(arb_mat_t r, int Nobs, double** in, slong prec) = 0;
#endif // WITH_ARB

  double *create_copy(double *in);
  void copy(double *out, double *in);
  void add(double *out, double *in1, double *in2);
  void add(double *out, double *in1, double a);
  void sub(double *out, double *in1, double *in2);
  void mult(double *out, double *in1, double *in2);
  void mult(double *out, double *in1, double a);
  void div(double *out, double *in1, double *in2);
  void div(double *out, double *in1, double a);
  void linear_comb(double *out, double a, double *in1, double b, double *in2);

  void mult_error(double* out,  double* in, double d);
  void add_error_quadrature(double* out,  double* in, double d);
};

class resampling_jack : public resampling_f {
public:
  resampling_jack();
  resampling_jack(int N);
  double**** create(int  N, int var, int t, double**** in);
  double comp_mean_unbias(double* in);
  double comp_error(double* in);
  double* create_fake(double mean, double error, int seed);
  double** comp_cov(int Nobs, double** in);
#ifdef WITH_ARB
  void comp_cov_arb(arb_mat_t r, int Nobs, double** in, slong prec);
#endif // WITH_ARB

};


class resampling_boot : public resampling_f {
public:
  resampling_boot();
  resampling_boot(int N);
  resampling_boot(int N, int o_seed);
  double**** create(int  N, int var, int t, double**** in);
  double comp_mean_unbias(double* in);
  double comp_error(double* in);
  double* create_fake(double mean, double error, int seed);
  double** comp_cov(int Nobs, double** in);
#ifdef WITH_ARB
  void comp_cov_arb(arb_mat_t r, int Nobs, double** in, slong prec);
#endif // WITH_ARB

};

double** covariance_jack(int Nobs, int Np1, double** in);
double** covariance_boot(int Nobs, int Np1, double** in);

#endif