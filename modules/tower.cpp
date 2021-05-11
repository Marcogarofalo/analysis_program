#define tower_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "tower.hpp"


void free_tower(int size ,void  **p){
    int i;
    
    for(i=0;i<size;i++)
        free(p[i]);
    free(p);
}

void free_2(int size ,double  **p){
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

void free_4( int size1,int size2,int size3, double ****p){
    int i,j,k;
    
    for(i=0;i<size1;i++){
        for(j=0;j<size2;j++){
            for(k=0;k<size3;k++){
                free(p[i][j][k]);
            }
            free(p[i][j]);
        }
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

int  **int_malloc_2(int size1,int size2){
    int i;
    int **p;
    
    p=(int**) malloc(sizeof(int*)*size1);
    for(i=0;i<size1;i++){
        p[i]=(int*) malloc(sizeof(int)*size2);
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

double  ****double_malloc_4(int size1,int size2, int size3,int size4){
    int i,j;
    double ****p;
    
    p=(double****) malloc(sizeof(double***)*size1);
    for(i=0;i<size1;i++){
        p[i]=(double***) malloc(sizeof(double**)*size2);
        for(j=0;j<size2;j++){
            p[i][j]=(double**) malloc(sizeof(double*)*size3);
            for(int k=0;k<size3;k++){
                p[i][j][k]=(double*) malloc(sizeof(double)*size4);
            }
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

template<typename T>
T   **malloc_2(int size1,int size2){
    int i;
    T **p;
    
    p=(T**) malloc(sizeof(T*)*size1);
    for(i=0;i<size1;i++){
        p[i]=(T*) malloc(sizeof(T)*size2);
    }
    return p;    
}
template<typename T>
T  **calloc_2(int size1,int size2){
    int i;
    T **p;
    
    p=(T**) malloc(sizeof(T*)*size1);
    for(i=0;i<size1;i++){
        p[i]=(T*) calloc(size2,sizeof(T));
    }
    return p;    
}

template<typename T>
T  ***malloc_3(int size1,int size2, int size3){
    int i,j;
    T ***p;
    
    p=(T***) malloc(sizeof(T**)*size1);
    for(i=0;i<size1;i++){
        p[i]=(T**) malloc(sizeof(T*)*size2);
        for(j=0;j<size2;j++){
            p[i][j]=(T*) malloc(sizeof(T)*size3);
        }
    }
    
    return p;    
}
