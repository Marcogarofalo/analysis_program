#define eigensystem_C
#include <iostream>
#include <complex>
#include <Eigen/Eigenvalues>  
#include <Eigen/Dense>  
#include "eigensystem.hpp"
#include "sorting.hpp"
// #include <vector>

void complexEigenproblem(double** A, int N, double** eigenvalues, double** eigenvectors) {
    int i, j;
    Eigen::MatrixXcd a(N, N);
    typedef std::complex<double> C;
    for (i = 0;i < N;i++) {
        for (j = 0;j < N;j++) {
            a(i, j) = C(A[i + j * N][0], A[i + j * N][1]);
        }
    }
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(a);

    for (i = 0;i < N;i++) {
        eigenvalues[i][0] = real(ces.eigenvalues()(i));
        eigenvalues[i][1] = imag(ces.eigenvalues()(i));
    }
    for (i = 0;i < N;i++) {
        for (j = 0;j < N;j++) {
            eigenvectors[i + N * j][0] = real(ces.eigenvectors()(i, j));
            eigenvectors[i + N * j][1] = imag(ces.eigenvectors()(i, j));
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



int generalysed_Eigenproblem(double** A, double** B, int N, double*** eigenvalues, double*** eigenvectors, int verbosity) {
    int i, j;
    Eigen::MatrixXcd a(N, N);
    Eigen::MatrixXcd b(N, N);
    typedef std::complex<double> C;
    for (i = 0;i < N;i++) {
        for (j = 0;j < N;j++) {
            a(i, j) = C(A[i + j * N][0], A[i + j * N][1]);
            b(i, j) = C(B[i + j * N][0], B[i + j * N][1]);
        }
    }
    Eigen::MatrixXcd c(N, N);
    c = b.inverse() * a;

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    // Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ces(a,b);
    ces.compute(c);

    //   Eigen::LLT<Eigen::MatrixXcd> lltOfA(c); // compute the Cholesky decomposition of A
    //   if(lltOfA.info() == Eigen::NumericalIssue){
    //        std::cout<<"non semi-positive definitie matrix!"<<std::endl;
    //        std::cout<< lltOfA.info() << std::endl;
    //        std::cout<< c<< std::endl;
    //        for(int i=0;i<N;i++)
    //            std::cout<< ces.eigenvalues()[i]<< std::endl;
    //   }    
    int err = 0;
    if (verbosity >= 0) {
        for (int i = 0;i < N;i++) {
            if (imag(ces.eigenvalues()[i]) > 1e-8 || real(ces.eigenvalues()[i]) < 0) {
                err = 1;
            }
        }
        for (int i = 0;i < N;i++) {
            for (int j = 0;j < N;j++) {
                if (imag(ces.eigenvectors().col(i)[j]) > 1e-8) {
                    err = 2;
                }
            }
        }
        if (err >= 1) {
            std::cout << "non semi-positive definitie matrix!" << std::endl;
            std::cout << c << std::endl;
            for (int i = 0;i < N;i++) {
                std::cout << ces.eigenvalues()[i] << std::endl;
            }
        }
        if (err >= 2) {
            std::cout << "eigenvectors (columns) are coplex!" << std::endl;
            std::cout << ces.eigenvectors() << std::endl;
        }
        if (verbosity >= 3) {
            printf("the GEVP\n");
            printf("{");
            for (int i = 0;i < N;i++) {
                printf("{");
                for (int j = 0;j < N;j++)
                    printf("%.15f,\t", real(c(i, j)));
                printf("},\n");
            }
            printf("}\n");
            for (int i = 0;i < N;i++) {
                std::cout << ces.eigenvalues()[i] << std::endl;
            }
            std::cout << ces.eigenvectors() << std::endl;
        }
    }

    Eigen::VectorXcd v;

    double* lambda = (double*)malloc(sizeof(double) * N);
    int* order = (int*)malloc(sizeof(int) * N);
    for (int i = 0;i < N;i++) {
        lambda[i] = real(ces.eigenvalues()[i]);
        // lambda[i] = (ces.eigenvalues()[i]);
        order[i] = i;
    }
    quickSort(order, lambda, 0, N - 1);

    for (i = 0;i < N;i++) {
        int ii = order[N - 1 - i];
        (*eigenvalues)[i][0] = lambda[ii];
        (*eigenvalues)[i][1] = 0;
        v = ces.eigenvectors().col(ii);

        Eigen::VectorXcd res = c * v - ces.eigenvalues()[ii] * v;
        // Eigen::VectorXcd res = a * v - ces.eigenvalues()[ii] * b* v;
        double sum = res.norm();
        if (sum > 1e-6) {
            printf("error eigenvalues: deviation = %g\n",sum);
            std::cout << c * v << std::endl;
            std::cout << "   l=" << ces.eigenvalues()[ii] << std::endl;
            std::cout << ces.eigenvalues()[ii] * v << std::endl;
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
        int sing = 1;
        double vmax = real(v(0));
        for (j = 1;j < N;j++) {
            if (fabs(vmax) < fabs(real(v(j)))) vmax = (real(v(j)));
        }
        if (vmax < 0) sing = -1;
        for (j = 0;j < N;j++) {
            (*eigenvectors)[j + i * N][0] = sing * real(v(j));
            (*eigenvectors)[j + i * N][1] = sing * imag(v(j));
        }
    }
    free(lambda);free(order);
    return err;
}


void GEVP_real(double** A, double** B, int N, double*** eigenvalues, double*** eigenvectors, int verbosity) {
    int i, j;
    Eigen::MatrixXd a(N, N);
    Eigen::MatrixXd b(N, N);
    for (i = 0;i < N;i++) {
        for (j = 0;j < N;j++) {
            a(i, j) = A[i + j * N][0];
            b(i, j) = B[i + j * N][0];
        }
    }
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
    ges.compute(a, b);

    int err = 0;
    if (verbosity >= 0) {
        for (int i = 0;i < N;i++) {
            if (imag(ges.eigenvalues()[i]) > 1e-8 || real(ges.eigenvalues()[i]) < 0) {
                err = 1;
            }
        }
        for (int i = 0;i < N;i++) {
            for (int j = 0;j < N;j++) {
                if (imag(ges.eigenvectors().col(i)[j]) > 1e-8) {
                    err = 2;
                }
            }
        }
        if (err >= 1) {
            std::cout << "non semi-positive definitie matrix!" << std::endl;
            std::cout << a << std::endl;
            for (int i = 0;i < N;i++) {
                std::cout << ges.eigenvalues()[i] << std::endl;
            }
        }
        if (err >= 2) {
            std::cout << "eigenvectors (columns) are coplex!" << std::endl;
            std::cout << ges.eigenvectors() << std::endl;
        }

    }

    Eigen::VectorXcd v;
    // std::complex<double> norm;

    double* lambda = (double*)malloc(sizeof(double) * N);
    int* order = (int*)malloc(sizeof(int) * N);
    for (int i = 0;i < N;i++) {
        lambda[i] = real(ges.eigenvalues()[i]);
        order[i] = i;
    }
    quickSort(order, lambda, 0, N - 1);

    for (i = 0;i < N;i++) {
        int ii = order[N - 1 - i];
        (*eigenvalues)[i][0] = lambda[ii];
        (*eigenvalues)[i][1] = 0;
        v = ges.eigenvectors().col(ii);

        int sing = 1;
        double vmax = real(v(0));
        for (j = 1;j < N;j++) {
            if (fabs(vmax) < fabs(real(v(j)))) vmax = (real(v(j)));
        }
        if (vmax < 0) sing = -1;
        for (j = 0;j < N;j++) {
            (*eigenvectors)[j + i * N][0] = sing * real(v(j));
            (*eigenvectors)[j + i * N][1] = sing * imag(v(j));
        }
    }
    free(lambda);free(order);
}





void complexeigenSolver(double** A, int N, double*** eigenvalues, double*** eigenvectors) {
    int i, j;
    Eigen::MatrixXcd a(N, N);

    typedef std::complex<double> C;
    for (i = 0;i < N;i++) {
        for (j = 0;j < N;j++) {
            a(i, j) = C(A[i + j * N][0], 0 * A[i + j * N][1]);
        }
    }

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(a);

    Eigen::VectorXcd v;
    // std::complex<double> norm;

    double* lambda = (double*)malloc(sizeof(double) * N);
    int* order = (int*)malloc(sizeof(int) * N);
    for (int i = 0;i < N;i++) {
        lambda[i] = real(ces.eigenvalues()[i]);
        order[i] = i;
    }
    quickSort(order, lambda, 0, N - 1);

    for (i = 0;i < N;i++) {
        int ii = order[N - 1 - i];
        (*eigenvalues)[i][0] = lambda[ii];
        (*eigenvalues)[i][1] = 0;
        v = ces.eigenvectors().col(ii);

        for (j = 0;j < N;j++) {

            (*eigenvectors)[j + i * N][0] = real(v(j));
            (*eigenvectors)[j + i * N][1] = imag(v(j));
        }
    }
    free(lambda);free(order);

}
