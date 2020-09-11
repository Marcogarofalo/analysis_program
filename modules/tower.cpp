#define tower_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>



void free_2(int size ,double **p){
    int i;
    
    for(i=0;i<size;i++)
        free(p[i]);
    free(p);
}

void free_3( int size1,int size2, double ***p){
    int i,j;
    
    for(i=0;i<size1;i++){
        for(j=0;j<size2;j++)
            free(p[i][j]);
        free(p[i]);
    }
    free(p);
}

double  **double_malloc_2(int size1,int size2){
    int i;
    double **p;
  
    p=(double**) malloc(sizeof(double*)*size1);
    for(i=0;i<size1;i++){
        p[i]=(double*) malloc(sizeof(double)*size2);
    }
    return p;    
}
double  **double_calloc_2(int size1,int size2){
    int i;
    double **p;
  
    p=(double**) malloc(sizeof(double*)*size1);
    for(i=0;i<size1;i++){
        p[i]=(double*) calloc(size2,sizeof(double));
    }
    return p;    
}

double  ***double_malloc_3(int size1,int size2, int size3){
    int i,j;
    double ***p;
  
    p=(double***) malloc(sizeof(double**)*size1);
    for(i=0;i<size1;i++){
        p[i]=(double**) malloc(sizeof(double*)*size2);
        for(j=0;j<size2;j++){
            p[i][j]=(double*) malloc(sizeof(double)*size3);
        }
    }
    
   return p;    
}


double **swap_indices(int N,int Njack, double **in){
    double **out;
    out=double_malloc_2(Njack,N);
    
    for (int i =0;i<N;i++)
        for (int j =0;j<Njack;j++)
            out[j][i]=in[i][j];
    return out;
}
