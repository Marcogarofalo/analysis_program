#ifndef gamma_analysis_H
#define gamma_analysis_H


double  *analysis_gamma (  int var, int rep, int nconf,int flow, double *a , double *function(int var, int order,int flow ,double *ah));
void gamma_correlator(char** option, struct kinematic kinematic_2pt,
    char* name, double**** data, int Confs, const char* plateaux_masses, FILE* outfile,
    int index, const char* description);

double *identity_func_gamma(int var, int order,int flow ,double *ah);
double *mass_gamma(int var, int order,int flow ,double *ah);

void gamma_correlator_func(char** option, struct kinematic kinematic_2pt,
    char* name, double**** data, int Confs, const char* plateaux_masses, FILE* outfile,
    int index, const char* description,double* function(int var, int order, int flow, double* ah));

#endif
 
