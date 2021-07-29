#define eigensystem_H
#include <iostream>
#include <complex>
#include <Eigen/Eigenvalues>  
#include <Eigen/Dense>  
#include "eigensystem.hpp"
#include <vector>

void complexEigenproblem(double **A, int N, double **eigenvalues, double **eigenvectors )
{
  int i,j;  
  Eigen::MatrixXcd a(N, N);
  typedef std::complex<double> C;
  for(i=0;i<N;i++){
      for(j=0;j<N;j++){
          a(i,j)= C(A[i+j*N][0],A[i+j*N][1]);   
      }
  }
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
  ces.compute(a);

  for(i=0;i<N;i++){
      eigenvalues[i][0]=real(ces.eigenvalues()(i)); 
      eigenvalues[i][1]=imag(ces.eigenvalues()(i));
  }
  for(i=0;i<N;i++){
      for(j=0;j<N;j++){
        eigenvectors[i+N*j][0]=real(ces.eigenvectors()(i,j)); 
        eigenvectors[i+N*j][1]=imag(ces.eigenvectors()(i,j));
      }
  }
  
    
}

/*
void real_generalysed_Eigenproblem(double ***A, double ***B, int N, double **eigenvalues, double ***eigenvectors )
{
  int i,j;  
  Eigen::MatrixXf a(N, N);
  Eigen::MatrixXf b(N, N);
  for(i=0;i<N;i++){
      for(j=0;j<N;j++){
          a(i,j)= A[i][j][0];   
          b(i,j)= B[i][j][0];   
      }
  }
  Eigen::GeneralizedEigenSolver<Eigen::MatrixXf> ges;
  ges.compute(a, b);
  for(i=0;i<N;i++){
      eigenvalues[i][0]=real(ges.eigenvalues()(i)); 
      eigenvalues[i][1]=imag(ges.eigenvalues()(i));
  }
  for(i=0;i<N;i++){
      for(j=0;j<N;j++){
        eigenvectors[i][j][0]=real(ges.eigenvectors()(i,j)); 
        eigenvectors[i][j][1]=imag(ges.eigenvectors()(i,j));
      }
  }

}
 */

// A utility function to swap two elements
void swap(int* a, int* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

/* This function takes last element as pivot, places
 *  the pivot element at its correct position in sorted
 *   array, and places all smaller (smaller than pivot)
 *  to left of pivot and all greater elements to right
 *  of pivot */
int partition (int *order, double *arr, int low, int high)
{
    double pivot = arr[order[high]];    // pivot
    int i = (low - 1);  // Index of smaller element
    
    for (int j = low; j <= high- 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (arr[order[j]] <= pivot)
        {
            i++;    // increment index of smaller element
            swap(&order[i], &order[j]);
        }
    }
    swap(&order[i + 1], &order[high]);
    return (i + 1);
}

/* The main function that implements QuickSort
 * arr[] --> Array to be sorted,
 * low  --> Starting index,
 * high  --> Ending index */
void quickSort(int *order, double *arr, int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
         *          at right place */
        int pi = partition(order,arr, low, high);
        
        // Separately sort elements before
        // partition and after partition
        quickSort(order, arr, low, pi - 1);
        quickSort(order, arr, pi + 1, high);
    }
}

void generalysed_Eigenproblem(double **A, double **B, int N, double ***eigenvalues, double ***eigenvectors )
{
  int i,j;
  Eigen::MatrixXcd a(N, N);
  Eigen::MatrixXcd b(N, N);
  typedef std::complex<double> C;
  for(i=0;i<N;i++){
      for(j=0;j<N;j++){
          a(i,j)= C(A[i+j*N][0],0*A[i+j*N][1]);   
          b(i,j)= C(B[i+j*N][0],0*B[i+j*N][1]);   
      }
  }  
  Eigen::MatrixXcd c(N,N);
  c  = b.inverse()*a;
  
//   std::cout<< "matrix c(t)\n" << a <<std::endl;
//   std::cout<< "matrix c(t0)\n" << b <<std::endl;
//   std::cout<< "matrix c(t0)^-1\n" << b.inverse() <<std::endl;
//   std::cout<< "matrix I\n" << b.inverse()*b <<std::endl;
//   std::cout<< "matrix I\n" << b*b.inverse() <<std::endl;
//   std::cout<< "matrix c(t0)^{-1}*c(t)\n" << c <<std::endl;

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
  ces.compute(c);

  Eigen::VectorXcd v;
  std::complex<double> norm;
  
  double *lambda=(double*) malloc(sizeof(double)*N);
  int *order=(int*) malloc(sizeof(int)*N);
  for(int i=0;i<N;i++){
      lambda[i]=real(ces.eigenvalues()(i));
      order[i]=i;
  }
  quickSort(order, lambda, 0,N-1);
  
  for(i=0;i<N;i++){
      int ii=order[N-1-i];
      (*eigenvalues)[i][0]=lambda[ii];
      (*eigenvalues)[i][1]=0;
      v=ces.eigenvectors().col(i);
      norm=v.transpose()*a*v;
      //std::cout<< norm<<std::endl;
      norm=sqrt(norm);
      for(j=0;j<N;j++){
          
          (*eigenvectors)[j+ii*N][0]=real(v(j));
          (*eigenvectors)[j+ii*N][1]=imag(v(j));
      }
  }
  free(lambda);free(order);
  
//   for(i=0;i<N;i++){
//       v=ces.eigenvectors().col(i);
//       norm=v.transpose()*a*v;
//       //std::cout<< norm<<std::endl;
//       norm=sqrt(norm);
//     //  v=-v/norm;
//       (*eigenvalues)[i][0]=real(ces.eigenvalues()(i));
//       (*eigenvalues)[i][1]=imag(ces.eigenvalues()(i));
//       for(j=0;j<N;j++){
//         (*eigenvectors)[j+i*N][0]=real(v(j));
//         (*eigenvectors)[j+i*N][1]=imag(v(j));
//       }
//       /* the eigenvalues are lambda_n ~ exp(-En(t-t0))+exp(-En(T-(t-t0)))
//        * I am interesting in E0 which is smaller the the others E0 < En
//        * So I want the biggest lambda lambda_0 > lambda_1 */
//       if ((*eigenvalues)[0][0] < (*eigenvalues)[i][0]){ 
//   	    (*eigenvalues)[i][0]=(*eigenvalues)[0][0];
//   	    (*eigenvalues)[i][1]=(*eigenvalues)[0][1];
//       	    for(j=0;j<N;j++){
//         	(*eigenvectors)[j+i*N][0]=(*eigenvectors)[j+0*N][0];
// 	        (*eigenvectors)[j+i*N][1]=(*eigenvectors)[j+0*N][0];
//             }
//       	    
// 	        (*eigenvalues)[0][0]=real(ces.eigenvalues()(i));
//             (*eigenvalues)[0][1]=imag(ces.eigenvalues()(i));
//       	    for(j=0;j<N;j++){
//         	(*eigenvectors)[j+0*N][0]=real(v(j));
// 	        (*eigenvectors)[j+0*N][1]=imag(v(j));
//             }
// 
//       }
//    }
      
      
  
}
