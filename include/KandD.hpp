#ifndef  KandD_H
#define KandD_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

double one_line(int n, int Nvar, double *x,int Npar,double  *P);
double two_lines(int n, int Nvar, double *x,int Npar,double  *P);

double **fit_MK_double_chiral_FVE(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1);
double **fit_B0_double_chiral_FVE(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1);

double **fit_MK_double_chiral_FVE_P40(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1);

double **fit_MD_double_chiral(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack r1);
double **fit_MD_double_chiral_P30(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack r1);

double **fit_MDs_double_chiral(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack r1);




//from Meson mass
double **fit_fK_double_chiral_FVE_from_M(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1);
double **fit_fD_chiral_continuum_from_M(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1);
double **fit_fDs_chiral_continuum_from_M(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1);

double **w0_from_MDs(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1);



// K_chiral_fit.c
double **fit_MK_double_chiral_B0_FVE(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1);

double         **fit_Omegaw0_from_M(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1);


//KandD_clover
double ** fit_MK_fK_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info,  struct result_jack *r1, const char *prefix ,char **argv);
double **fit_fKoverfpi_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info , struct result_jack *r1, const char *prefix,char **argv);

double **fit_MK_Mpi_fK_fpi_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info , struct result_jack *r1, const char *prefix,char **argv);

double **fit_MK_fK_chiral_spline_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info , struct result_jack *r1, const char *prefix,char **argv);

double **fit_MK_Mpi_fK_fpi_chiral_spline_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info , struct result_jack *r1, const char *prefix,char **argv);



#endif
