#ifndef FVE_H
#define FVE_H


double   *H_over_H0_vir(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym );

double   *H_minus_H0_HA_vir(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest, double *Zf_PS_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym );


#endif


