#ifndef m_eff1_H
#define m_eff1_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>


double   *compute_effective_mass1(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,  int index , struct observable *obs, const char *description);

#endif
