#ifndef HLT_H
#define HLT_H
#include "global.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"
#include "non_linear_fit.hpp"

enum HLT_b{
    HLT_EXP_b=0
} ;

class  HLT_type {
public:
    int Tmax;
    int T;
    int Njack;
    int id;
    double E0;
    std::vector<double> Es = {};
    std::vector<double> lambdas = {};
    double **A;
    double *R;

    HLT_type(int tmax, int L0, double E0, int njack, HLT_b type_b, double alpha=0 );
    ~HLT_type();
    
    double* HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
    FILE* outfile, double fun_of_corr(int, double****, int, struct fit_type), const char* description, struct HLT_type HLT_info, FILE* file_jack);


};



#endif // !HLT_H
