#ifndef lhs_function_H
#define lhs_function_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "linear_fit.hpp"
#include "non_linear_fit.hpp"
#include "mutils.hpp"
#include "resampling.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "global.hpp"
#include "eigensystem.hpp"
#include "gamma_analysis.hpp"
#include "tower.hpp"

using namespace std;

// template<int tau>
double r_AWI(int j, double ****in,int t ,struct fit_type fit_info){
    int T=file_head.l0;
    //int tptau=(t+tau)%T;
    double num= in[j][8][t][0];//(in[j][0][(t+1)%T][0]-in[j][0][t][0])*in[j][7][tptau][0];
    double den= in[j][9][t][0];//(in[j][5][t][0])*in[j][7][tptau][0];
    return (-2*num/(4.*den)) ;
    //return num;
}



double m_PCAC(int j, double ****in,int t ,struct fit_type fit_info){
    int T=file_head.l0;
    double num= (in[j][0][(t+1)%T][0]-in[j][0][t][0]);
    double den= (in[j][1][t][0]);
    return (num/(4*den)) ;
    //return num;
}

#endif

