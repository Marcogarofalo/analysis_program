#ifndef  global_fit_KandD_H
#define global_fit_KandD_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

struct fit_result fit_fK_of_M(struct database_file_jack  *jack_files,  struct header *head ,int Njack, struct data_jack *gJ,int ***mass_index, struct result_jack *r1 ,struct fit_type fit_info);
struct fit_result global_fit_Omega_from_Mpi_MK(struct database_file_jack  *jack_files,  struct header *head ,int Njack, struct data_jack *gJ,int ***mass_index, struct result_jack *r1 ,struct fit_type fit_info);
#endif
