#ifndef FVE_H
#define FVE_H


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>


void  FVE(double w0, double w0a,double dl3b, double dl4b, double Bw,double fw,int Lsize, double mw ,double  dmpi2, double dfpi,double *KM, double *Kf);
void  FVE_K(double Bw,double fw,double  frac_Lw, double mlw, double msw ,double  dmpi2, double dfpi,double  dmK2, double dfK,double *KM, double *Kf);
void  FVE_Mpi(double afm,double dl1phys, double dl2phys,double dl3phys, double dl4phys,int Lsize,double  ampi, double afpi,int nloop, int  ischeme, double *KM, double *Kf);

double FVE_GL(double Lsize_w,double  aml, double af, double aB);
double FVE_GL_fast(double Lsize_w,double  aml, double af, double aB);
double FVE_GL_Mpi(double Lsize_w /* L/w0 */,double  xi,double f_PS);


#endif

