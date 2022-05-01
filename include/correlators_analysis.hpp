#ifndef correlators_analysis_H
#define correlators_analysis_H
#include "resampling.hpp"
#include "global.hpp"
#include "non_linear_fit.hpp"
double laplacian_M_eff_T(int t, int T, double** in);
double der2corr_M_eff_T(int t, int T, double** in);
double M_eff_T(int t, int T, double** in);
double M_eff_sinh_T(int t, int T, double** in);
double M_eff_log(int t, int T, double** in);
double identity(int t, int T, double** in);
double shift_corr(int t, int T, double** in);
double M_eff_log_shift(int t, int T, double** in);
double shift_and_M_eff_sinh_T(int t, int T, double** in);
double log_corr(int t, int T, double** in);


double M_eff_t0_sinh_T(int t, int t0, int T, double** in);
double M_eff_sinh_T_ct_ctp(int t, int T, double ct, double ctp);


void check_correlatro_counter(int i);

double constant_fit(int n, int Nvar, double* x, int Npar, double* P);

double shift_and_M_eff_acosh(int t, int T, double** in);


double* plateau_correlator_function(char** option, struct kinematic kinematic_2pt, char* name, double**** conf_jack, int Njack, const char* plateaux_masses, FILE* outfile, int index, const char* description, double (*fun)(int, int, double**), FILE* file_jack, struct fit_type fit_info = default_fit_info);
// jackknife plateau_correlator_function(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,const char  *plateaux_masses,FILE *outfile,  int index , const char *description , double (*fun)(int ,int  , double ** ),  FILE * file_jack);


struct fit_result fit_function_to_corr(char** option, struct kinematic kinematic_2pt, char* name, double**** conf_jack, const char* plateaux_masses, FILE* outfile, int index, int re_im, const char* description, struct fit_type fit_info, FILE* file_jack);
/****************************************************************************
 *
 * option[1] = blind/read_plateaux/see
 * option[3] = path to the directory with subdirectories: jack , out
 * option[4] = resampling type
 * option[5] = do not set to "pdf"
 * option[6] = basename to look for in the plateaux.txt file
 * description = name to give to the fit, the one get by r
****************************************************************************/
struct fit_result fit_fun_to_fun_of_corr(char** option, struct kinematic kinematic_2pt, char* name, double**** conf_jack, const char* plateaux_masses, FILE* outfile, double fun_of_corr(int, double****, int, struct fit_type), const char* description, struct fit_type fit_info, FILE* file_jack);


void add_correlators(char** option, int& ncorr_conf_jack, double****& conf_jack, double** fun_of_corr(int, double****, int, struct fit_type), struct fit_type fit_info);
void add_one_correlators(char** option, int& ncorr_conf_jack, double****& conf_jack, struct fit_type fit_info, double** r);
void zero_corr(double* zeros, int Njack, FILE* jack_file);

void write_jack(double* corr, int Njack, FILE* jack_file);
#endif

