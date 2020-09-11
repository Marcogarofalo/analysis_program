#define eigensystem_H
#include <iostream>
#include <complex>
#include <Eigen/Eigenvalues>  
#include <Eigen/Dense>  
#include "eigensystem.hpp"


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
void generalysed_Eigenproblem(double **A, double **B, int N, double **eigenvalues, double **eigenvectors )
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
  
/*  std::cout<< "matrix c(t)\n" << a <<std::endl;
  std::cout<< "matrix c(t0)\n" << b <<std::endl;
  std::cout<< "matrix c(t0)^-1\n" << b.inverse() <<std::endl;
  std::cout<< "matrix I\n" << b.inverse()*b <<std::endl;
  std::cout<< "matrix I\n" << b*b.inverse() <<std::endl;
  std::cout<< "matrix c(t0)^-1*c(t)\n" << c <<std::endl;
*/
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
  ces.compute(c);

  Eigen::VectorXcd v;
  std::complex<double> norm;
  for(i=0;i<N;i++){
      v=ces.eigenvectors().col(i);
      norm=v.transpose()*a*v;
      //std::cout<< norm<<std::endl;
      norm=sqrt(norm);
    //  v=-v/norm;
      eigenvalues[i][0]=real(ces.eigenvalues()(i));
      eigenvalues[i][1]=imag(ces.eigenvalues()(i));
      for(j=0;j<N;j++){
        eigenvectors[j+i*N][0]=real(v(j));
        eigenvectors[j+i*N][1]=imag(v(j));
      }
      /* the eigenvalues are lambda_n ~ exp(-En(t))+exp(-En(T-t))
       * I am interesting in E0 which is smaller the the others E0 < En
       * then i want the beggest lambda lambda_0 > lambda_1 */
      if (eigenvalues[0][0] < eigenvalues[i][0]){ 
  	    eigenvalues[i][0]=eigenvalues[0][0];
  	    eigenvalues[i][1]=eigenvalues[0][1];
      	    for(j=0;j<N;j++){
        	eigenvectors[j+i*N][0]=eigenvectors[j+0*N][0];
	        eigenvectors[j+i*N][1]=eigenvectors[j+0*N][0];
            }
      	    
	        eigenvalues[0][0]=real(ces.eigenvalues()(i));
            eigenvalues[0][1]=imag(ces.eigenvalues()(i));
      	    for(j=0;j<N;j++){
        	eigenvectors[j+0*N][0]=real(v(j));
	        eigenvectors[j+0*N][1]=imag(v(j));
            }

      }
  }
}
