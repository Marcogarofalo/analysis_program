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


double ration_corr_min_half(int j, double ****in,int t ,struct fit_type fit_info){
    int inum=fit_info.corr_id[0];
    int iden=fit_info.corr_id[1];
    int T=file_head.l0;
    //int tptau=(t+tau)%T;
    double num= in[j][inum][t][0];//(in[j][0][(t+1)%T][0]-in[j][0][t][0])*in[j][7][tptau][0];
    double den= in[j][iden][t][0];//(in[j][5][t][0])*in[j][7][tptau][0];
    return (-2*num/(4.*den)) ;
    //return num;
}

// template<int tau>
double r_AWI(int j, double ****in,int t ,struct fit_type fit_info){
    int T=file_head.l0;
    //int tptau=(t+tau)%T;
    double num= in[j][8][t][0];//(in[j][0][(t+1)%T][0]-in[j][0][t][0])*in[j][7][tptau][0];
    double den= in[j][9][t][0];//(in[j][5][t][0])*in[j][7][tptau][0];
    return (-2*num/(4.*den)) ;
    //return num;
}
double r_AWI_loc(int j, double ****in,int t ,struct fit_type fit_info){
    int T=file_head.l0;
    double num= in[j][13][t][0];//(in[j][0][(t+1)%T][0]-in[j][0][t][0])*in[j][7][tptau][0];
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


double m_PCAC_loc(int j, double ****in,int t ,struct fit_type fit_info){
    int T=file_head.l0;
    double num= (in[j][11][(t+1)%T][0]-in[j][11][t][0]);
    double den= (in[j][1][t][0]);
    return (num/(4*den)) ;
    //return num;
}

double **num_awi_wallphi(int j, double ****in,int t,struct fit_type fit_info ){
    
    int id=fit_info.corr_id[0];
    int T=fit_info.T;
    error(fit_info.N!=1,1,"der2_der2_corr","works only with one corr");
    double **r=double_malloc_2(fit_info.N,2);
    r[0][0]=in[j][id][(t+1)%T][0]-in[j][id][t][0];
    r[0][1]=in[j][id][(t+1)%T][1]-in[j][id][t][1];
    return r; 
}

#endif

