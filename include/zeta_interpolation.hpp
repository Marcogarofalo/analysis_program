#ifndef zeta_interpolation_H
#define zeta_interpolation_H

#include <vector>
#include "global.hpp"
#include "header_phi4.hpp"



class  zeta_interpolation{
private:
    double ****grid;// [en][ index for dvec]  [ mass =0,1,2 mean-err,mean,mean+err] [ (kmax-kmin)/h  ]
    double ****kint;
    int etot;
    int ntot=5;
//     double h=1e-4;
    double h=1e-3;
    double maxk;
    double hg=0.01;
    double maxg=4;
    std::vector<int> L;
    int Nm=7;
//      int Nm=3;
    double  **m;//[en]  [ mass =0,1,2 mean-err,mean,mean+err] 
    std::vector< std::vector<int> > mom;
    double ****krange; // [en][ index for dvec]  [mass][ kmin,kmax, iters] 
    bool allocated=1;
public:
    zeta_interpolation(){};
    void Init(char *resampling, std::vector<int>  myen,  std::vector<cluster::IO_params> paramsj, std::vector<data_phi> gjack );
    double compute(double inL, int n, double mass,  double k );
    ~zeta_interpolation();
} ;



class  zeta_interpolation_qsqg{
private:
    double ***grid;// [ index for dvec]  [ gamma ] [ qsq  ]
    double ***qsqrange;// [ index for dvec]  [ gamma ] [ 3: qsqmin, qsqmax, Nqsq  ]
    int ntot=5;
    double hqsq=1e-3;
    double maxqsq=3;
    double hg=0.1;
    double maxg=3;
    double Ng;
    std::vector< std::vector<int> > mom;
    bool allocated=1;
public:
    zeta_interpolation_qsqg(){};
    void Init( );
    double compute(double qsq , int n, double gamma );
    ~zeta_interpolation_qsqg();
} ;

EXTERN zeta_interpolation zeta;
EXTERN zeta_interpolation_qsqg zeta_qsqg;

#endif

