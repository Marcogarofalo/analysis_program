#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include <unistd.h>
#include <sys/time.h>
#include <fcntl.h>
#include "global.hpp"
 

#include "linear_fit.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"
double f(int n ,int Nvar,double *x,int Npar,double*P){
    
    return  P[0] * exp(x[0] * P[1]) ;
}

int main(){
    double **y=double_malloc_2(13,2);
    double **x=double_malloc_2(13,2);
    double **sigmax=double_malloc_2(13,2);
        
    x[0][0]=368.5;     y[0][0]=73.332;        sigmax[0][0]=0.66;    y[0][1]=1.5;
    x[1][0]=364.2;     y[1][0]=62.668;        sigmax[1][0]=0.66;    y[1][1]=1.0;
    x[2][0]=359.2;     y[2][0]=52.004;        sigmax[2][0]=0.66;    y[2][1]=0.8;
    x[3][0]=354.5;     y[3][0]=44.006;        sigmax[3][0]=0.66;    y[3][1]=0.7;
    x[4][0]=348.7;     y[4][0]=34.675;        sigmax[4][0]=0.66;    y[4][1]=1.2;
    x[5][0]=343.7;     y[5][0]=28.010;        sigmax[5][0]=0.66;    y[5][1]=1.6;
    x[6][0]=338.7;     y[6][0]=22.678;        sigmax[6][0]=0.66;    y[6][1]=1.2;
    x[7][0]=334.2;     y[7][0]=17.346;        sigmax[7][0]=0.66;    y[7][1]=1.5;
    x[8][0]=329.0;     y[8][0]=14.680;        sigmax[8][0]=0.66;    y[8][1]=1.6;
    x[9][0]=324.0;     y[9][0]=10.681;        sigmax[9][0]=0.66;    y[9][1]=1.2;
    x[10][0]=319.1;    y[10][0]=8.015;        sigmax[10][0]=0.66;   y[10][1]=0.8;
    x[11][0]=314.6;    y[11][0]=6.682;        sigmax[11][0]=0.66;   y[11][1]=1.0;
    x[12][0]=308.7;    y[12][0]=5.349;        sigmax[12][0]=0.66;   y[12][1]=1.5;
        
    
    int N=1;
    int *en=(int*) malloc(sizeof(int)*1);
    en[0]=13;
    int Nvar=1;
    int Npar=2;
    double *guess=(double*) malloc(sizeof(double)*2);
    guess[0]=1; guess[1]=0.001;
    //double  *non_linear_fit_Nf_sigmax(int N, int *ensemble ,double **x, double **sigmax, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess ){

    printf("only y errors\n");
    double *P=non_linear_fit_Nf( N, en ,x, y , Nvar,  Npar,  f ,guess ).P;
    for(int i=0;i<Npar;i++)
        printf("P%d=%g  \n",i,P[i]);
    int i,k,j;
    double **cov=covariance_non_linear_fit_Nf(N, en,x, y,P , Nvar,  Npar, f ); 
     for (i=0;i<Npar;i++){
            for (k=0;k<i;k++)
                    cov[i][k]/=sqrt(cov[i][i]*cov[k][k]);
            for (k=i+1;k<Npar;k++)
                    cov[i][k]/=sqrt(cov[i][i]*cov[k][k]);
    }
    for (i=0;i<Npar;i++){
        for (j=0;j<Npar;j++){
              printf("%g\t",  cov[i][j] );
        }
        if (i!=Npar) printf(" \n");
        else printf("\n");
    }
    
    free(P);
    printf("xy errors\n");
    P=non_linear_fit_Nf_sigmax( N, en ,x, sigmax, y , Nvar,  Npar,  f ,guess );
    for( i=0;i<Npar;i++)
        printf("P%d=%g  \n",i,P[i]);
    
    
    free(P);
    printf("xy errors with diagonal covariance \n");
    double **cov_yx1=double_malloc_2(en[0]*2,en[0]*2);
    for (i=0;i<en[0]*2;i++)
        for (j=0;j<en[0]*2;j++)
            cov_yx1[i][j]=0;
    for (i=0;i<en[0];i++)     
        cov_yx1[i][i]=1./(y[i][1]*y[i][1]);
    for (i=en[0];i<en[0]*2;i++)     
        cov_yx1[i][i]=1./(sigmax[i-en[0]][0]*sigmax[i-en[0]][0]);
        
    P=non_linear_fit_Nf_sigmax_covariance( N, en ,x, sigmax, y , Nvar,  Npar,  f ,guess,cov_yx1 );
    for(int i=0;i<Npar;i++)
        printf("P%d=%g  \n",i,P[i]);
    
    free(P);
    printf("xy errors   ,x iterative procedure\n");
    P=non_linear_fit_Nf_sigmax_iterative( N, en ,x, sigmax, y , Nvar,  Npar,  f ,guess );
    for(int i=0;i<Npar;i++)
        printf("P%d=%g  \n",i,P[i]);
    
    return 0;
}


