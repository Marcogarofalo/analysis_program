#ifndef tower_H
#define tower_H

#include <stdlib.h>

void free_tower(int size ,void **p);
void free_2(int size ,double **p);
void free_3( int size1,int size2, double ***p);
void free_4( int size1,int size2,int size3, double ****p);

double  ***double_malloc_3(int size1,int size2, int size3);
double  **double_malloc_2(int size1,int size2);
int  **int_malloc_2(int size1,int size2);
double  **double_calloc_2(int size1,int size2);
double  ****double_malloc_4(int size1,int size2, int size3,int size4);
    
double **swap_indices(int N,int Njack, double **in);

template<typename T>  T  ****malloc_4(int size1,int size2, int size3, int size4);
template<typename T>  T  ***malloc_3(int size1,int size2, int size3);
template<typename T>  T  **malloc_2(int size1,int size2);
template<typename T>  T  **calloc_2(int size1,int size2);



#endif
