#ifndef indicaes_H
#define indices_H
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

 int index_n_twopt(int si,int ii,int ix0,int ik1,int ik2,int imom1,int imom2);
 int index_n_twopt_fit(int ik1,int ik2,int imom1,int imom2);
 int  index_ensemble_twopt_fit(struct header head,int ik1,int ik2,int imom1,int imom2);
 int index_n_twopt_G_fit(int ikt,int iks, int imom0, int imomt,int imoms);
int  index_ensemble_twopt_G_fit(struct header head,int ikt,int iks, int imom0, int imomt,int imoms);
 int index_n_threept(int si,int ii,int ix0,int ik1,int ik2,int ik3,int imom1,int imom2);
 int index_n_twoptgamma(int si,int ii,int ix0,int ikt,int iks,int imom0,int imomt,int imoms);
int index_n_minus_kappa(int ik);
 int index_n_minus_theta(int imom);
 
 
 int  index_mesonic_twopt(int si,int ii,int ix0,int ik1,int ik2,int imom1,int imom2);
#endif
 
