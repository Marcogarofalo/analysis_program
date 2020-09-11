#ifndef m_eff_H
#define m_eff_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

double *fit_plateaux_G(char **option,struct kinematic_G kinematic_2pt_G , char* name,const char *description,double **mt,double **r, int Njack,FILE *plateaux_masses,FILE *outfile);


extern double *constant_fit_m_eff(int M, double in);
extern double M_eff_t( int t,int L0, double ***in, int num_corr);


extern double *create_jmeff_from_jcorr(int Njack,int L0, int tmin ,int tmax, double ****in, int num_corr );

extern double ratio_i( int t, double ***in);
extern double ratio_0( int t, double ***in);


extern double   *compute_effective_mass(char **option,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,FILE **plateaux_masses,FILE *outfile ,int index, const char* description);
extern double   *compute_effective_mass_GEVP(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack, int var ,FILE *plateaux_masses,FILE *outfile );
extern double   *compute_f_PS_ls_ss(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,int Njack ,FILE *plateaux_masses,FILE *outfile );
double   *compute_f_PS_GEVP(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack,double *mass_jack_fit_k2k1, int Njack, int var ,FILE *plateaux_masses,FILE *outfile );


double Rmur(int t, double ***in,double mass, double oPp,struct kinematic_G kinematic_2pt_G,int index,int *sym);

extern double   *compute_Zf_PS_ll(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,int Njack ,FILE *plateaux_masses,FILE *outfile );
extern double   *compute_oPp_ll(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,int Njack ,FILE *plateaux_masses,FILE *outfile );
extern double   *compute_f_PS_ll(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1, double *oPp_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile );

extern double   *compute_Rmur(char **option ,struct kinematic_G kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1, double* mass_rest, double *oPp_PS_jack_fit, int Njack ,FILE *plateaux_masses,FILE *outfile , int index,int *sym);
double         *compute_CAmur(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack, int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym );
double         *compute_Rmur_from_meff(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest, int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym );

double   *H_over_H0(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest, double *oPp_PS_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym );
double   *H_minus_H0(char **option ,struct kinematic_G kinematic_2pt_G , char* name,double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest, double *oPp_PS_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym );
double   *H_minus_H0_HA(char **option ,struct kinematic_G kinematic_2pt_G , char* name,double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest, double *oPp_PS_jack_fit, double *Zf_PS_jack_fit ,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym );


extern double   *compute_effective_mass_out_max_twist(char **option,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,FILE *plateaux_masses,FILE *outfile ,int index);
extern double   *compute_Zf_PS_ll_out_max_twist(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,int Njack ,FILE *plateaux_masses,FILE *outfile );
double   *compute_f_PS_ls_ss_out_max_twist(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,int Njack ,FILE *plateaux_masses,FILE*plateaux_mpcac,FILE *outfile );

double   *compute_effective_mass_out_max_twist_K(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,FILE *plateaux_masses,FILE *outfile,  int index  );

#endif
