#ifndef pion_H
#define pion_H


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h> 


double Mw2_fw(int n, int Nvar, double *x,int Npar,double  *P);
double Mw2_fw_a0_minus(  double x,int Npar,double  *P);
double Mw2_fw_polynomial(int n, int Nvar, double *x,int Npar,double  *P);
double fw_over_Mw_of_Mw_minus_exp(double  input, double x,int Npar,double  *P);
double Mw2_over_fw2_chiral_FVE_a0_minus(double  input, double x,int Npar,double  *P);



double m_over_f_xi(double xi,int Npar,double *P);
double m_over_f(double x,int Npar,double *P);
double m_over_f_pol(double x,int Npar,double *P);
double **fit_Mpi_fw(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ );
double **fit_Mpi_fw_polynomial(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ );

double Mw2_fw_polynomial_6(int n, int Nvar, double *x,int Npar,double  *P);
double m_over_f_pol_6(double x,int Npar,double *P);
double **fit_Mpi_fw_polynomial_6(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ );

double Mw2_fw_chiral_c2(int n, int Nvar, double *x,int Npar,double  *P);

double Mw2_fw_chiral_c2_a0_minus(  double x,int Npar,double  *P);
double **fit_Mpi_fw_chiral_c2(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ );


double Mw2_fw_chiral_FVE(int n, int Nvar, double *x,int Npar,double  *P);
double **fit_Mpi_fw_chiral_FVE(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ );
double Mw2_fw_chiral_FVE_P1_P3(int n, int Nvar, double *x,int Npar,double  *P);
double **fit_Mpi_fw_chiral_FVE_P1_P3(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ );

double **fit_Mpi_fw_chiral_FVE_prior(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ );
double **fit_Mpi_fw_chiral_FVE_P40_prior(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ );

double Mw2_chiral_FVE_a0_minus( double, double x,int Npar,double  *P);
double Mw2_chiral_FVE_a0_minus_P1_P3( double, double x,int Npar,double  *P);

double fPSw_chiral_FVE(  double x,int Npar,double  *P);
double fPSw_chiral_FVE_P1_P3(  double x,int Npar,double  *P);
double **fit_Mpi_fw_chiral_FVE_flag(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ );

double fw_of_Mw2_chiral_FVE(int n, int Nvar, double *x,int Npar,double  *P);
double **fit_fw_of_Mw_chiral_FVE(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ );
double fw_of_Mw2_physical_point(double x,int Npar,double  *P);

double **fit_Mpi_fw_chiral_FVE_P40(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ );
double **fit_fw_of_Mw_chiral_FVE_P40(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ );


//Pion_clover.c
struct fit_result fit_Mpi_fw_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info);
struct fit_result fit_Mpi_fw_chiral_FVE_clover_treshold(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info, double threshold);
struct fit_result fit_Mpi_fwMpi4_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info);
struct fit_result fit_Mpi_fwMpi4_chiral_FVE_clover_threshold(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info, double threshold);

#endif
