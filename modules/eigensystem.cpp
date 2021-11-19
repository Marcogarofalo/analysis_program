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
  
  
  
  
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
  ces.compute(c);
  
  Eigen::LLT<Eigen::MatrixXcd> lltOfA(c); // compute the Cholesky decomposition of A
  if(lltOfA.info() == Eigen::NumericalIssue){
       std::cout<<"non semi-positive definitie matrix!"<<std::endl;
       std::cout<< lltOfA.info() << std::endl;
       std::cout<< c<< std::endl;
       for(int i=0;i<N;i++)
           std::cout<< ces.eigenvalues()[i]<< std::endl;
  }    

  Eigen::VectorXcd v;
  std::complex<double> norm;
  
  double *lambda=(double*) malloc(sizeof(double)*N);
  int *order=(int*) malloc(sizeof(int)*N);
  for(int i=0;i<N;i++){
      lambda[i]=real(ces.eigenvalues()[i]);
      order[i]=i;
  }
  quickSort(order, lambda, 0,N-1);
  
  for(i=0;i<N;i++){
      int ii=order[N-1-i];
      (*eigenvalues)[i][0]=lambda[ii];
      (*eigenvalues)[i][1]=0;
      v=ces.eigenvectors().col(ii);
      
      Eigen::VectorXcd res=c*v-ces.eigenvalues()[ii]*v;
      double sum=0;
      for (int i =0; i<N;i++)
          sum+=real(res(i));
      if ( sum > 1e-6 ){
          printf("error eigenvalues\n");
          std::cout <<c*v<< std::endl;
          std::cout <<  "   l=" <<ces.eigenvalues()[ii]<< std::endl;
          std::cout <<ces.eigenvalues()[ii]*v<< std::endl;
          exit(1);
      }
//       Eigen::VectorXcd v1=ces.eigenvectors().col(order[(N-1-i+1)%N]);
//       Eigen::VectorXcd v2=ces.eigenvectors().col(order[(N-1-i+2)%N]);
//       if ( v1.adjoint()*(v) > 1e-6 ){
//             printf("error eigevectors not ortogonal\n");
//             
//       }
     /* if (N==3){
        std::cout<< lambda[order[(N-1-i+2)%N]] << "\n";  
        std::cout<< v2 << "\n\n";
        
        std::cout<< lambda[order[(N-1-i+1)%N]] << "\n";  
        std::cout<< v1 << "\n\n";
        std::cout<< lambda[order[ii]] << "\n";  
        std::cout<< v << "\n\n";
        std::cout<< v1.adjoint()*(v) << std::endl;
      }        */  
      int sing =1;
      double vmax=fabs(real(v(0)));
      for(j=1;j<N;j++){
        if (fabs(vmax) < fabs(real(v(j))) ) vmax=(real(v(j)));
      }
      if (vmax <0) sing=-1;
      for(j=0;j<N;j++){
          (*eigenvectors)[j+ii*N][0]=sing*real(v(j));
          (*eigenvectors)[j+ii*N][1]=sing*imag(v(j));
      }
  }
  free(lambda);free(order);
  
  
}





void complexeigenSolver(double **A, int N, double ***eigenvalues, double ***eigenvectors )
{
  int i,j;
  Eigen::MatrixXcd a(N, N);
  
  typedef std::complex<double> C;
  for(i=0;i<N;i++){
      for(j=0;j<N;j++){
          a(i,j)= C(A[i+j*N][0],0*A[i+j*N][1]);   
      }
  }  
  
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
  ces.compute(a);

  Eigen::VectorXcd v;
  std::complex<double> norm;
  
  double *lambda=(double*) malloc(sizeof(double)*N);
  int *order=(int*) malloc(sizeof(int)*N);
  for(int i=0;i<N;i++){
      lambda[i]=real(ces.eigenvalues()[i]);
      order[i]=i;
  }
  quickSort(order, lambda, 0,N-1);
  
  for(i=0;i<N;i++){
      int ii=order[N-1-i];
      (*eigenvalues)[i][0]=lambda[ii];
      (*eigenvalues)[i][1]=0;
      v=ces.eigenvectors().col(ii);
 
      for(j=0;j<N;j++){
          
          (*eigenvectors)[j+ii*N][0]=real(v(j));
          (*eigenvectors)[j+ii*N][1]=imag(v(j));
      }
  }
  free(lambda);free(order);
  
}
