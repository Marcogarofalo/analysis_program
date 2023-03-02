#ifndef read_nissa_H
#define read_nissa_H

#include <fstream>
#include <string.h>

#include "tower.hpp"

void read_single_nissa(double** out, std::string contraction, std::string gamma, std::string namefile);
double*** return_all_nissa(std::string namefile, int Ncorr, int T, bool check = true);

#endif // !read_nissa_H

