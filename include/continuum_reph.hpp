#ifndef continuum_reph_H
#define continuum_reph_H
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

struct fit_result fit_FA_pion_generic(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info);
struct fit_result fit_FAxg_pion_generic(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info);
struct fit_result fit_FAV_5_pion(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ , int AV);

struct fit_result fit_FAV_K(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info);
struct fit_result fit_FAV_D(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info);
struct fit_result fit_FAV_Ds(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info);


struct fit_result fit_HA_pion_generic(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info);
struct fit_result fit_HAV_K(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info);
struct fit_result fit_HAV_D(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info);
struct fit_result fit_HAV_Ds(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info);


double   *compute_Rmur_auto_plateau(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,double *mass_rest,double *oPp_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym );
double   *H_over_H0_autoplateaux(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,double *mass_rest,double *oPp_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym );
double   *H_minus_H0_autoplateaux(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,double *mass_rest,double *oPp_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym );

double **subtract_non_Lorentz_and_int_Ds(char **argv, const char  *string,double **phys_point,double xmin,double xmax,  struct header *head ,int Njack ,int npar_fun,double *fit_function(int,double),struct fit_result fit_out, struct fit_type fit_info  , struct reph_jack *gJ,int e,const char *AV);
double **subtract_non_Lorentz_and_int_K(char **argv, const char  *string,double **phys_point,double xmin,double xmax,  struct header *head ,int Njack ,int npar_fun,double *fit_function(int,double),struct fit_result fit_out, struct fit_type fit_info  , struct reph_jack *gJ,int e,const char *AV);
double **subtract_non_Lorentz_and_int(char **argv, const char  *string,int ikt,int iks,double **kp_tot,double xmin,double xmax,  struct header *head ,int Njack ,int npar_fun,double *fit_function(int,double),struct fit_result fit_out, struct fit_type fit_info  , struct reph_jack *gJ,int e,const char *AV);
double derivative_xi(int ind,int Njack ,struct fit_result fit_out, struct fit_type fit_info, double *x);


double inter_2(double x1, double x2, double y1, double y2, double x);
double inter_4(double xs1, double xs2,double xc1, double xc2 , double ys1c1, double ys2c1, double ys1c2, double ys2c2, double xs, double xc);

struct fit_result fit_FAV_Kphys(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info);
struct fit_result fit_FAV_Dphys(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info);
struct fit_result fit_FAV_Dsphys(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info);


struct fit_result fit_FAV_pion_treshold_Mpi(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info,double r0_min, double r0_max);
struct fit_result fit_FAV_Kphys_treshold(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info, double r0_min, double r0_max);
struct fit_result fit_FAV_Dphys_treshold(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info, double r0_min, double r0_max);
struct fit_result fit_FAV_Dsphys_treshold(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info,double r0_min, double r0_max);
#endif
