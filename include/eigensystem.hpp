#ifndef eigensystem_H
#define eigensystem_H

#include <iostream>
#include <complex>
#include <Eigen/Eigenvalues>  
#include <Eigen/Dense>  



void complexEigenproblem(double **A, int N, double **eigenvalues, double **eigenvectors );
/*
void real_generalysed_Eigenproblem(double ***A, double ***B, int N, double **eigenvalues, double ***eigenvectors );
 */
int generalysed_Eigenproblem(double **A, double **B, int N, double ***eigenvalues, double ***eigenvectors,int verbosity=0 );
void GEVP_real(double **A, double **B, int N, double ***eigenvalues, double ***eigenvectors , int verbosity=0);

#endif
