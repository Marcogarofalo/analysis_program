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

#include "header_phi4.hpp"
using namespace std;


template<int ix,int iy,int iz>
double m_eff_of_sum(int j, double ****in,int t ,struct fit_type fit_info){
    
    double ct= (in[j][ix][t][0]+in[j][iy][t][0]+in[j][iz][t][0])/3. ;
    double ctp= (in[j][ix][t+1][0]+in[j][iy][t+1][0]+in[j][iz][t+1][0])/3. ;
    
    return M_eff_in_inp(t,ct,ctp);
}

template<int ix,int iy,int iz>
double sum_corr_directions_shift(int j, double ****in,int t ,struct fit_type fit_info){
    
    return (in[j][ix][t][0]+in[j][iy][t][0]+in[j][iz][t][0])/3. ;
    
}


double lhs_four_BH_0(int j, double ****in,int t ,struct fit_type fit_info){
    
    double r;
    int T=file_head.l0;
    double L3=(double)file_head.l1*file_head.l2*file_head.l3;
    
    r=in[j][8][t][0];
    r-=in[j][0][(T/2-t+T)%T][0] *in[j][0][T/8][0] ;
    r-=in[j][0][(T/2-T/8)][0] *in[j][0][t][0] ;
    r-=in[j][0][T/2][0]*in[j][0][(t-T/8+T)%T][0] ;
    
    r/=(in[j][0][T/2][0] * in[j][0][ (t-T/8+T)%T  ][0]  );
    
    // to_do:
    // I think that there is a 4 M factor missing
    r*=L3;
    r/=(t-T/8);
    
    
    return r;
    
}
inline double den(double m0, double m1  , int t,int T){
    return exp(-T*m0-t*m1+t*m0) +exp(-T*m1-t*m0+t*m1);
}

template<int id>
double lhs_E2_div_shift(int j, double ****in,int t ,struct fit_type fit_info){
    
    double r;
    int T=file_head.l0;
    
    double m0=fit_info.ext_P[0][j];
    double m1=fit_info.ext_P[1][j];
    
    
    double ct=in[j][id][t][0]/ den(m0,m1, t ,T)  ;
    ct-=in[j][id][t+1][0]/  den(m0,m1, t+1 ,T)  ;
    //ct*=den;
    
    double ctp1=in[j][id][t+1][0]/  den(m0,m1, t+1 ,T)  ;
    ctp1-=in[j][id][t+2][0]/  den(m0,m1, t+2 ,T) ;
    //ctp1*=den;
    
    r=M_eff_sinh_T_ct_ctp(t,T, ct, ctp1);
    //error(1==1,1,"lhs_E2_div_shift","t+1 deve essere diviso il suo denominatore");
    return r;
    
}

double lhs_four_BH_1(int j,double ****in,int t ,struct fit_type fit_info){
    
    double r;
    int T=file_head.l0;
    double L3=(double)file_head.l1*file_head.l2*file_head.l3;
    
    r=in[j][9][t][0];
    r-=in[j][1][(T/2-t+T)%T][0] *in[j][1][T/8][0] ;
    r-=in[j][1][(T/2-T/8)][0] *in[j][1][t][0] ;
    r-=in[j][1][T/2][0]*in[j][1][(t-T/8+T)%T][0] ;
    
    r/=(in[j][1][T/2][0] * in[j][1][ (t-T/8+T)%T  ][0]  );
    r*=L3;
    r/=(t-T/8);
    
    return r;
}

double lhs_four_BH(int j, double ****in,int t ,struct fit_type fit_info){
    
    double r;
    int T=file_head.l0;
    double L3=(double)file_head.l1*file_head.l2*file_head.l3;
    double disc;
    
    r=in[j][10][t][0];
    r/=( in[j][1][ (t-T/8+T)%T  ][0]  *in[j][0][T/2][0]  );
    
    disc=1.;
    
    r-=disc;
    r*=L3;
    r/=(t-T/8);
    
    return r;
}




double lhs_four_BH_0_s(int j, double ****in,int t ,struct fit_type fit_info ,int t1, int t2,int t4, int id){
    
    double r;
    int T=file_head.l0;
    double L3=(double)file_head.l1*file_head.l2*file_head.l3;
    t1=t1%T;
    t2=t2%T;
    t4=t4%T;
    
    r=in[j][id][t][0];
    r-=in[j][0][(t4-t+T)%T][0] *in[j][0][(t2-t1+T)%T][0] ;
    r-=in[j][0][(t4-t2+T)%T][0] *in[j][0][(t-t1+T)%T][0] ;
    r-=in[j][0][(t4-t1+T)%T][0]*in[j][0][(t-t2+T)%T][0] ;
    
    r/=(in[j][0][(t4-t1+T)%T][0] * in[j][0][ (t-t2+T)%T  ][0]  );
    
    r*=L3/2;
    r/=(t-t2);
    
    return r;
    
}


double lhs_four_BH_1_s(int j, double ****in,int t ,struct fit_type fit_info ,int t1, int t2,int t4, int id){
    
    double r;
    int T=file_head.l0;
    double L3=(double)file_head.l1*file_head.l2*file_head.l3;
    int t25=((T*2)/5);
    t1=t1%T;
    t2=t2%T;
    t4=t4%T;
    
    r=in[j][id][t][0];
    r-=in[j][1][(t4-t+T)%T][0] *in[j][1][(t2-t1+T)%T][0] ;
    r-=in[j][1][(t4-t2+T)%T][0] *in[j][1][(t-t1+T)%T][0] ;
    r-=in[j][1][(t4-t1+T)%T][0]*in[j][1][(t-t2+T)%T][0] ;
    
    r/=(in[j][1][(t4-t1+T)%T][0] * in[j][1][ (t-t2+T)%T  ][0]  );
    
    r*=L3/2;
    r/=(t-t2);
    
    return r;
    
}



double lhs_four_BH_s(int j, double ****in,int t ,struct fit_type fit_info ,int t1, int t2,int t4, int id){
    
    double r;
    int T=file_head.l0;
    double L3=(double)file_head.l1*file_head.l2*file_head.l3;
      
    r=in[j][id][t][0];
    r/=(in[j][0][(t4-t1+T)%T][0] * in[j][1][ (t-t2+T)%T  ][0]  );
    r-=1.;
    r*=((double) L3/2.);
    
    r/=(t-t2);
    
    return r;
    
}


double lhs_four_BH_10_s(int j, double ****in,int t ,struct fit_type fit_info ,int t1, int t2,int t4, int id){
    
    double r;
    int T=file_head.l0;
    double L3=(double)file_head.l1*file_head.l2*file_head.l3;
      
    r=in[j][id][t][0];
    r/=(in[j][1][(t4-t1+T)%T][0] * in[j][0][ (t-t2+T)%T  ][0]  );
    r-=1.;
    r*=((double) L3/2.);
    r/=(t-t2);
    
    return r;
    
}




double lhs_four_BH_0_03t16(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_0_s(j, in, t, fit_info, 0 ,3, 16, 15);
}
double lhs_four_BH_1_03t16(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_1_s(j, in, t, fit_info, 0 ,3, 16, 16);
}
double lhs_four_BH_03t16(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_s(j, in, t, fit_info, 0 ,3, 16, 17);
}


double lhs_four_BH_0_04t16(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_0_s(j, in, t, fit_info, 0 ,4, 16, 18);
}
double lhs_four_BH_1_04t16(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_1_s(j, in, t, fit_info, 0 ,4, 16, 19);
}
double lhs_four_BH_04t16(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_s(j, in, t, fit_info, 0 ,4, 16, 20);
}


double lhs_four_BH_0_03t20(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_0_s(j, in, t, fit_info, 0 ,3, 20, 21);
}
double lhs_four_BH_1_03t20(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_1_s(j, in, t, fit_info, 0 ,3, 20, 22);
}
double lhs_four_BH_03t20(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_s(j, in, t, fit_info, 0 ,3, 20, 23);
}


double lhs_four_BH_0_04t20(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_0_s(j, in, t, fit_info, 0 ,4, 20, 24);
}
double lhs_four_BH_1_04t20(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_1_s(j, in, t, fit_info, 0 ,4, 20, 25);
}
double lhs_four_BH_04t20(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_s(j, in, t, fit_info, 0 ,4, 20, 26);
}


double lhs_four_BH_0_05t20(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_0_s(j, in, t, fit_info, 0 ,5, 20, 27);
}
double lhs_four_BH_1_05t20(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_1_s(j, in, t, fit_info, 0 ,5, 20, 28);
}
double lhs_four_BH_05t20(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_s(j, in, t, fit_info, 0 ,5, 20, 29);
}


double lhs_four_BH_10_03t16(int j, double ****in,int t ,struct fit_type fit_info){
    return lhs_four_BH_10_s(j, in, t, fit_info, 0 ,3, 16, 30);
}




double lhs_four_BH_01_tx_tf(int j, double ****in,int t ,struct fit_type fit_info ,int tx, int tf,int comp0 ,int comp1,int id){
    
    double r;
    int T=file_head.l0;
    double L3=(double)file_head.l1*file_head.l2*file_head.l3;
    
    r=in[j][id][t][0];
    r/=(in[j][comp0][tf%T][0] * in[j][comp1][ (t-tx+T)%T  ][0]  );
    r-=1.;
    r*=((double) L3/2.);
    
    //r/=(t-tx);
    return r;
}


template<int delta ,int tx, int tf, int comp0 ,int comp1,int id >
double lhs_four_BH_01_tx_tf_shifetd(int j, double ****in,int t ,struct fit_type fit_info ){
    
    double r=-lhs_four_BH_01_tx_tf(j, in, t , fit_info,tx,tf,comp0,comp1,id );
    
    r+=lhs_four_BH_01_tx_tf(j, in, t+delta , fit_info,tx,tf,comp0,comp1,id );
    r/=delta;
    return r;
}


double lhs_four_BH_00_tx_tf(int j, double ****in,int t ,struct fit_type fit_info ,int tx, int tf,int comp0 ,int comp1,int id){
    
    double r;
    int T=file_head.l0;
    double L3=(double)file_head.l1*file_head.l2*file_head.l3;
    
    
    r=in[j][id][t][0];
    r-=in[j][comp0][tf%T][0] * in[j][comp1][ (t-tx+T)%T  ][0];
    r-=in[j][comp0][(tf-t+T)%T][0] * in[j][comp1][ tx%T  ][0];
    r-=in[j][comp0][(tf-tx)%T][0] * in[j][comp1][ t  ][0];
    
    r/=(in[j][comp0][tf%T][0] * in[j][comp1][ (t-tx+T)%T  ][0]  );
    r*=((double) L3/2.);
    
    //r/=(t-tx);
    return r;
}


template<int delta ,int tx, int tf, int comp0 ,int comp1,int id >
double lhs_four_BH_00_tx_tf_shifetd(int j, double ****in,int t ,struct fit_type fit_info ){
    
    double r=-lhs_four_BH_00_tx_tf(j, in, t , fit_info,tx,tf,comp0,comp1,id );
    r+=lhs_four_BH_00_tx_tf(j, in, t+delta , fit_info,tx,tf,comp0,comp1,id );
    r/=delta;
    return r;
}


double lhs_four_BH_no_sub(int j, double ****in,int t ,struct fit_type fit_info){
    
    double r;
    int T=file_head.l0;
    double L3=(double)file_head.l1*file_head.l2*file_head.l3;
    int t25=((T*2)/5);
    
    r=in[j][10][t][0];
    r/=( in[j][1][ (t-T/8+T)%T  ][0]  *in[j][0][T/2][0]  );
    
    r*=L3;
    
    return r;
}

double GEVP_shift_matrix(int j, double ****in,int t,struct fit_type fit_info ){
    double ct,ctp;
    int N=2;
    int T=file_head.l0;
    double **M=double_calloc_2(N*N,2);// [NxN] [reim ]
    double **Mt0=double_calloc_2(N*N,2);
    
    double **lambdat=double_malloc_2(N,2);// [N] [reim]
    double **lambdatp1=double_malloc_2(N,2);// [N] [reim]
    double **vec=double_malloc_2(N*N,2);
    int t0=3;
    double r;
    
    
    double s[4],s0[4];
    
    int i00=2, i11=3, i01=12;
    
    s[0]=in[j][i00][t][0]-in[j][i00][(t+1)%T][0]; // two0_to_two0
    s[1]=in[j][i01][t][0]-in[j][i01][(t+1)%T][0]; // two0_to_two1
    s[3]=in[j][i11][t][0]-in[j][i11][(t+1)%T][0]; // two1_to_two1
    
    //t
    M[0][0]=s[0];
    M[1][0]=s[1];
    M[3][0]=s[3];
    M[2][0]=M[1][0];
    
    s0[0]=in[j][i00][t0%T][0]-in[j][i00][(t0+1)%T][0]; // two0_to_two0
    s0[1]=in[j][i01][t0%T][0]-in[j][i01][(t0+1)%T][0]; // two0_to_two1
    s0[3]=in[j][i11][t0%T][0]-in[j][i11][(t0+1)%T][0]; // two1_to_two1
    
    Mt0[0][0]=s0[0];
    Mt0[1][0]=s0[1];
    Mt0[3][0]=s0[3];
    Mt0[2][0]=Mt0[1][0];
    
    generalysed_Eigenproblem(M,Mt0,N,&lambdat,&vec); 
    
    
    //t+1
    s[0]=in[j][i00][(t+1)%T][0]-in[j][i00][(t+2)%T][0]; // two0_to_two0
    s[1]=in[j][i01][(t+1)%T][0]-in[j][i01][(t+2)%T][0]; // two0_to_two1
    s[3]=in[j][i11][(t+1)%T][0]-in[j][i11][(t+2)%T][0]; // two1_to_two1
    
    //t
    M[0][0]=s[0];
    M[1][0]=s[1];
    M[3][0]=s[3];
    M[2][0]=M[1][0];
    
    
    generalysed_Eigenproblem(M,Mt0,2,&lambdatp1,&vec); 
    
    if((t-t0)>=0)
        r=M_eff_sinh_T_ct_ctp(t-t0,T, lambdat[0][0], lambdatp1[0][0]);
    else 
        r=M_eff_sinh_T_ct_ctp( t-t0, T,  lambdat[1][0], lambdatp1[1][0]);
    
    
    free_2(N*N,M);
    free_2(N*N,Mt0);
    free_2(N,lambdat);
    free_2(N,lambdatp1);
    free_2(N*N,vec);
    
    return r;
}


#endif

