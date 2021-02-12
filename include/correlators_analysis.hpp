#ifndef correlators_analysis_H
#define correlators_analysis_H

double M_eff_T(  int t, int T, double **in);
double M_eff_sinh_T(  int t, int T, double **in);
double M_eff_log(  int t, int T, double **in);
double identity(  int t, int T, double **in);
double shift_corr(  int t, int T, double **in);
double M_eff_log_shift(  int t, int T, double **in);
double constant_fit(int n, int Nvar, double *x,int Npar,double  *P);

double two_particle_energy(int t,int T , double **in);
double   *plateau_correlator_function(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,FILE **plateaux_masses,FILE *outfile,  int index , const char *description , double (*fun)(int ,int  , double ** ));

struct fit_result fit_function_to_corr(char **option ,struct kinematic kinematic_2pt ,  char* name, double ****conf_jack ,FILE **plateaux_masses,FILE *outfile,  int index, int re_im , const char *description , struct fit_type fit_info,  char * namefile_jack );

struct fit_result fit_fun_to_fun_of_corr(char **option ,struct kinematic kinematic_2pt ,  char* name, double ****conf_jack ,FILE **plateaux_masses,FILE *outfile,  double fun_of_corr(int,double****,int,struct fit_type ) , const char *description , struct fit_type fit_info,  char * namefile_jack );


#endif

