#ifndef read_nissa_H
#define read_nissa_H

#include <fstream>
#include <string.h>
#include <vector>
#include "tower.hpp"

void read_single_nissa(double** out, std::string contraction, std::string gamma, std::string namefile);
double*** return_all_nissa(std::string namefile, int Ncorr, int T, bool check = true);
void read_all_nissa(double*** out, std::string namefile, int Ncorr, int T, bool check = true);

class struct_nissa_out_info {
public:
    int Ncorr;
    int Ncontr;
    int Ngamma;
    int T;
    std::vector<std::string> gamma;
    struct_nissa_out_info(std::string namefile);
    void print();
    std::vector<int> inidices_of_gamma(std::vector<std::string> mygammas);
};
#endif // !read_nissa_H

