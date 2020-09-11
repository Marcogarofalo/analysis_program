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



int main(){
   double  **y,*x,**b,**r;
   int N=2,i,j,k;
   y=(double**) malloc(sizeof(double*)*N);
   r=(double**) malloc(sizeof(double*)*N);
  for (i=0;i<N;i++){
        y[i]=(double*) malloc(sizeof(double)*N);
        r[i]=(double*) calloc(N,sizeof(double));     
  }
   for (i=0;i<N;i++){
        for (j=i;j<N;j++){
            y[i][j]=(rand())*1.0/RAND_MAX;
            y[j][i]=y[i][j];
        if (i==j) y[i][j]+=1;
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
   
}

