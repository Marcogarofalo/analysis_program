#ifndef correlators_analysis_H
#define correlators_analysis_H

double M_eff_T(  int t, int T, double **in);
double M_eff_sinh_T(  int t, int T, double **in);

double two_particle_energy(int t,int T , double **in);
double   *plateau_correlator_function(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,FILE **plateaux_masses,FILE *outfile,  int index , const char *description , double (*fun)(int ,int  , double ** ));


#endif

