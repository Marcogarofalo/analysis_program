#ifndef fit_all_H
#define fit_all_H


#include "global.hpp"

class generic_header {
public:
    int T;
    int L;
    double mu;
    std::vector<double> mus;
    std::vector<double> theta;

};

class data_all {

public:
    int Njack;
    int Nobs;
    double** jack;

    generic_header header;


};

struct fit_result fit_all_data(char** argv, std::vector<data_all> gjack,
    double lhs_fun(int, int, int, std::vector<data_all>, struct fit_type),
    struct fit_type fit_info, const char* label);

#endif

