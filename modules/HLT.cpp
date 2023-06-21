#define HLT_C
#include "HLT.hpp"
#include "global.hpp"
#include "tower.hpp"



HLT_type::HLT_type(int tmax, int L0, double E0, int njack, HLT_b type_b, double alpha ){
    Tmax=tmax;
    T=L0;
    Njack=njack;
    R=(double*) malloc(sizeof(double)*Tmax);
    A=malloc_2<double>(Tmax,Tmax);

    if(type_b==HLT_EXP_b){
        for (size_t t = 0; t < Tmax; t++)     {
            R[t]=1.0/(t+1.0)+1.0/(T-t-1.0);
            for (size_t r = 0; r < Tmax; r++)     {
                A[t][r]=exp(-(r+t+2-alpha)*E0) / (r+t+2-alpha);
                A[t][r]+=exp(-(T-r+t-alpha)*E0) / (T-r+t-alpha);
                A[t][r]+=exp(-(T+r-t-alpha)*E0) / (T+r-t-alpha);
                A[t][r]+=exp(-(2*T-r-t-2-alpha)*E0) / (2*T-r-t-2-alpha);
            }
        }
        

    }
}

HLT_type::~HLT_type(){
    free(R);
    for (size_t i = 0; i < Tmax; i++)   {
        free(A[i]);
    }
    free(A);
}


double * HLT_type::HLT_of_corr(char** option, double**** conf_jack, const char* plateaux_masses,
    FILE* outfile, double fun_of_corr(int, double****, int, struct fit_type), const char* description, struct HLT_type HLT_info, FILE* file_jack){

    double** r=malloc_2<double>(Tmax, Njack);
    for (int t = 0;t < Tmax;t++)
        for (int j = 0;j < Njack;j++)
            r[t][j] = conf_jack[j][id][t][0];


    double **cov = myres->comp_cov( Tmax,  r);

};
