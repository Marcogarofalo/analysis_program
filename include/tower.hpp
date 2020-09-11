#ifndef tower_H
#define tower_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

void free_2(int size ,double **p);
void free_3( int size1,int size2, double ***p);

double  ***double_malloc_3(int size1,int size2, int size3);
double  **double_malloc_2(int size1,int size2);
double  **double_calloc_2(int size1,int size2);

double **swap_indices(int N,int Njack, double **in);

#endif
