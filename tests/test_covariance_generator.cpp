#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "global.hpp"
#include "linear_fit.hpp"
#include "resampling.hpp"


int main(){
   double  **y,*x,**b,**r;
   int N=4,i,j,k;
   y=(double**) malloc(sizeof(double*)*N);
   r=(double**) malloc(sizeof(double*)*N);
  for (i=0;i<N;i++){
        y[i]=(double*) malloc(sizeof(double)*N);
        r[i]=(double*) calloc(N,sizeof(double));     
  }
  rand();
   for (i=0;i<N;i++){
        for (j=i;j<N;j++){
            y[i][j]=((double)rand())/RAND_MAX*3;
            y[j][i]=y[i][j];
            if (i==j) y[i][j]=fabs(y[i][j]+2);
        }
   }

   
   printf("MATRIX y\n");
   for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            printf("%g\t",y[i][j]);
        }
        printf("\n");
   }
   printf("cholesky_decomposition\n");
   b=cholesky_decomposition(y, N);
   for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            printf("%g\t",b[i][j]);
        }
        printf("\n");
   }
   for (i=0;i<N;i++){
        for (j=0;j<N;j++){
                for (k=0;k<N;k++)
                    r[i][j]+=b[i][k]*b[j][k];
        }
   }
      printf("L L^T-y\n");
   for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            printf("%g\t",r[i][j]-y[i][j]);
        }
        printf("\n");
   }
   printf("resampling array {0,10,20,..} with covariance y\n");
   double *mean=(double*) malloc(sizeof(double)*N);
   for (i=0;i<N;i++)
       mean[i]=i*10;
   int Njack=1e+4;
   double **jacknife=fake_sampling_covariance("jack",  mean, Njack ,N , y,1234);
   
   for (i=0;i<N;i++){
   double *res =mean_and_error_jack(Njack ,jacknife[i] );
   //double *res =mean_and_error("jack",Njack ,jacknife[0] ); //this does not work because return the last element of the jack series which is the mean
   
   printf("after resampling r[%d]=%f  err= %f  err^2=%f \n", i,res[0],res[1],res[1]*res[1]); free(res);
   }
  
   double **covJ=covariance("jack" , N,Njack, jacknife);
   printf("covariance of the jacknife \n");
   for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            printf("%g\t",covJ[i][j]);
        }
        printf("\n");
   }
   
}

