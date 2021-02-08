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
void generalysed_Eigenproblem(double **A, double **B, int N, double ***eigenvalues, double ***eigenvectors );
#endif
