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

template<int id, int ids>
double remove_exp_and_meff(int j, double ****in,int t ,struct fit_type fit_info){
    
    double E2=fit_info.ext_P[0][j];
    double A01=fit_info.ext_P[1][j];
    int T=file_head.l0;
    double cs=in[j][id][t][0]-A01*A01* in[j][ids][t][0]/2.0;
    double csp1=in[j][id][t+1][0]-A01*A01 * in[j][ids][t+1][0]/2.0;
    
    double mass=M_eff_in_inp(t, cs, csp1) ;
    
    return  mass;
 
}


template<int id>
double single_corr(int j, double ****in,int t ,struct fit_type fit_info){
    
    return in[j][id][t][0];
}

template<int id6, int id4, int id2>
double subtracted_phi3(int j, double ****in,int t ,struct fit_type fit_info){
    
    double phi2=fit_info.ext_P[0][j];
    return in[j][id6][t][0] - 2*in[j][id4][t][0]*phi2 + in[j][id2][t][0] * phi2*phi2;
}

template<int id4, int id2>
double subtracted_phi(int j, double ****in,int t ,struct fit_type fit_info){
    
    double phi2=fit_info.ext_P[0][j];
    return in[j][id4][t][0]-in[j][id2][t][0]*phi2;
}


template<int id, int idc>
double two_to_two_con(int j, double ****in,int t ,struct fit_type fit_info){
    int T=file_head.l0;
    double **ct=double_malloc_2(T, 1);
    for (int i =t; i< t+3;i++){
        ct[i][0]=in[j][id][i][0]- in[j][idc][i][0]*in[j][idc][i][0];
    }
    
    double r=shift_and_M_eff_sinh_T(  t,  T, ct);
    free_2(T,ct);
    return r;
}

template< int idc>
double one_to_one_sq(int j, double ****in,int t ,struct fit_type fit_info){
    int T=file_head.l0;
    double **ct=double_malloc_2(T, 1);
    for (int i =t; i< t+3;i++){
        ct[i][0]=in[j][idc][i][0]*in[j][idc][i][0];
    }
    
    double r=shift_and_M_eff_sinh_T(  t,  T, ct);
    free_2(T,ct);
    return r;
}

template<int id>
double me_oPp(int j, double ****in,int t ,struct fit_type fit_info){
    
    double E=fit_info.ext_P[0][j];
    double T=file_head.l0;
    double ct= in[j][id][t][0]/(exp(-E*t)+exp(-E*(T-t))  );
    ct*=2*E;
    
    return sqrt(ct);
}

double diff_meff(int j, double ****in,int t ,struct fit_type fit_info){
    int id1=fit_info.corr_id[0];
    int id2=fit_info.corr_id[1];
    if(fit_info.corr_id.size()!=2){printf("ratio_meff: corr_id must be a vector of the two masses you want to subtract\n");exit(1);}
    int T=file_head.l0;
    double ct= in[j][id1][t][0];
    double ctp= in[j][id1][(t+1)%T][0];
    double m1= M_eff_in_inp(t,ct,ctp);
     ct= in[j][id2][t][0];
     ctp= in[j][id2][(t+1)%T][0];
    double  m2= M_eff_in_inp(t,ct,ctp);

    return m1-m2;
}

double ratio_meff(int j, double ****in,int t ,struct fit_type fit_info){
    int id1=fit_info.corr_id[0];
    int id2=fit_info.corr_id[1];
    if(fit_info.corr_id.size()!=2){printf("ratio_meff: corr_id must be a vector of the two masses you want to subtract\n");exit(1);}
    
    int T=file_head.l0;
    double ct= in[j][id1][t][0];
    double ctp= in[j][id1][(t+1)%T][0];
    double m1= M_eff_in_inp(t,ct,ctp);
     ct= in[j][id2][t][0];
     ctp= in[j][id2][(t+1)%T][0];
    double m2= M_eff_in_inp(t,ct,ctp);

    return m1/m2;
}

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


template<int ix,int iy,int iz>
double sum_corr_weight_shift(int j, double ****in,int t ,struct fit_type fit_info){
    int T=file_head.l0;
    double m0=fit_info.ext_P[0][j];
    double m1=fit_info.ext_P[1][j];
    
    double ct=(in[j][ix][t][0]+in[j][iy][t][0]+in[j][iz][t][0])/3.;
    double ctp=(in[j][ix][t+1][0]+in[j][iy][t+1][0]+in[j][iz][t+1][0])/3.;
    
    ct=ct/( exp(- (m0+m1)*T/2.) * cosh((m0-m1)* (t-T/2.))  );
    ctp=ctp/( exp(- (m0+m1)*T/2.) * cosh((m0-m1)* (t+1.-T/2.))  );
    return ctp-ct;
}


template<int id,int tf>
double matrix_element_k3pi(int j, double ****in,int t ,struct fit_type fit_info){
    int T=file_head.l0;
    return (in[j][id][t][0])/sqrt(in[j][5][t][0]*in[j][1][(tf-t+T)%T][0] )  ;    
}
double matrix_element_k3pi_T_2(int j, double ****in,int t ,struct fit_type fit_info){
    int T=file_head.l0;
    return (in[j][119][t][0])/sqrt(in[j][5][t][0]*in[j][1][(T/2-t+T)%T][0] )  ;
//     return sqrt(in[j][5][t][0]*in[j][1][(T/2-t+T)%T][0] )  ;
    
}


template<int id,int tf>
double matrix_element_3pik(int j, double ****in,int t ,struct fit_type fit_info){
    int T=file_head.l0;
    return (in[j][id][t][0])/sqrt(in[j][1][t][0]*in[j][5][(tf-t+T)%T][0] )  ;    
}
double matrix_element_3pik_T_2(int j, double ****in,int t ,struct fit_type fit_info){
    int T=file_head.l0;
    return (in[j][123][t][0])/sqrt(in[j][1][t][0]*in[j][5][(T/2-t+T)%T][0] )  ;
    //     return sqrt(in[j][5][t][0]*in[j][1][(T/2-t+T)%T][0] )  ;
    
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
    
    
    generalysed_Eigenproblem(M,Mt0,N,&lambdatp1,&vec); 
    
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


void r_equal_value_or_vector(double &r, double **lambdat, double **vec, fit_type fit_info, int t, int t0 ){
    int n=fit_info.n;
    int N=fit_info.N; //if eigenvector N= component+ id_eigenvector * components
    if (fit_info.value_or_vector==0){
        if((t-t0)>=0)
            r=lambdat[n][0];
        else 
            r=lambdat[N-1-n][0];
    }       
    else if (fit_info.value_or_vector==1){
        
        if((t-t0)>=0){
            r=vec[n][0];
        }
        else{ 
            int comps=sqrt(N);
            int id=n/comps;
            int comp=n%comps;
            id = comp + (comps-1-id) *comps;
            r=vec[id][0];
        }
    }
    else{
        printf("you need to specify if you want the:\n\
                - eigenvalues (fit_info.value_or_vector=0) \n\
                - eigenvectors (fit_info.value_or_vector=1)");
        exit(1);
    }
}


/**********************************
 * you need to fill
 * std::vector<int> fit_info.corr_id
 * with the id of the correlators in the order:
 * [first row], [second row], ...
 * M_00, M_01, ..., M_0N, M_11, ...
**********************************/
double GEVP_matrix(int j, double ****in,int t,struct fit_type fit_info ){
    double ct,ctp;
    int N=fit_info.N;
    if (fit_info.value_or_vector==1){
        N=sqrt(fit_info.N);
        error(fit_info.N!=(N*N) ,1,"GEVP_matrix_p1",
          "when you want the eigenvector N must be the quare of the size of the matrix: fit_info.N=%d ", fit_info.N );
    }
    int ncorr=fit_info.corr_id.size();
    error(ncorr!=(N*N+N)/2 ,1,"GEVP_matrix",
          "you need to provide (N^2+N)/2 to populate the top triangular matrix NxN:\n  N=%d    ncorr=%d\n",N,ncorr  );
    
    int T=file_head.l0;
    double **M=double_calloc_2(N*N,2);// [NxN] [reim ]
    double **Mt0=double_calloc_2(N*N,2);
    
    double **lambdat=double_malloc_2(N,2);// [N] [reim]
    double **vec=double_malloc_2(N*N,2);
    int t0=fit_info.t0_GEVP%T;
    if (fit_info.GEVP_swap_t_t0){
        t0=t;
        t=fit_info.t0_GEVP%T;
    }
    else if (fit_info.GEVP_tpt0)
        t0= (fit_info.t0_GEVP+t)%T ; 
    double r;
    
    double *s=(double*) malloc(sizeof(double)*N);
    double *s0=(double*) malloc(sizeof(double)*N);
    
    
    
     //t
    int count=0;
    for (int i=0;i<N;i++){
        for (int k=i;k<N;k++){
            int corr_ik= fit_info.corr_id[count];
            int ik=i+k*N;
            int ki=k+i*N;
            //printf("%d  %g\n",ik,in[j][corr_ik][t][0]);
            M[ik][0]  = in[j][corr_ik][t][0];
            Mt0[ik][0]= in[j][corr_ik][t0][0];
            M[ki][0]=M[ik][0];
            Mt0[ki][0]=Mt0[ik][0];
            count++;
        }
        
    }
    if (t==0 && j==0 && fit_info.n==0){
        for (int i=0;i<N;i++){
            for (int j=0;j<N;j++)
                printf("%.5e\t",Mt0[i+j*N][0]);
            printf("\n");
        }
        
    }
//     error(!is_it_positive_lex_reim(Mt0, N) , 1, "GEVP_matrix:", "GEVP_matrix M(t0) not positive defined"  ) ;
    generalysed_Eigenproblem(M,Mt0,N,&lambdat,&vec); 
//     if (t==t0){
//         for (int i=0;i<N;i++)
//             printf("%g\t",lambdat[i][0]);
//         printf("\n");
//     }
    
//     int n=fit_info.n;
//     if((t-t0)>=0)
//         r=lambdat[n][0];
//     else 
//         r=lambdat[N-1-n][0];
            
    int n=fit_info.n;
    r_equal_value_or_vector(r,  lambdat, vec, fit_info,  t, t0);
    
    free_2(N*N,M);
    free_2(N*N,Mt0);
    free_2(N,lambdat);
//     free_2(N,lambdatp1);
    free_2(N*N,vec);
    
    return r;
}



/**********************************
 * you need to fill
 * std::vector<int> fit_info.corr_id
 * with the id of the correlators in the order:
 * [first row], [second row], ...
 * M_00, M_01, ..., M_0N, M_11, ... + M_00(py), M_00(pz), M_22(py), M_22(pz),
**********************************/
double GEVP_matrix_p1(int j, double ****in,int t,struct fit_type fit_info ){
    double ct,ctp;
    
    int N=fit_info.N;
    if (fit_info.value_or_vector==1){
        N=sqrt(fit_info.N);
        error(fit_info.N!=(N*N) ,1,"GEVP_matrix_p1",
          "when you want the eigenvector N must be the quare of the size of the matrix: fit_info.N=%d ", fit_info.N );
    }
    int ncorr=fit_info.corr_id.size()/3;
    error(ncorr!=(N*N+N)/2 ,1,"GEVP_matrix_p1",
          "you need to provide (N^2+N)/2 to populate the top triangular matrix NxN:\n  N=%d    ncorr=%d\n",N,ncorr  );
    
    int T=file_head.l0;
    double **M=double_calloc_2(N*N,2);// [NxN] [reim ]
    double **Mt0=double_calloc_2(N*N,2);
    
    double **lambdat=double_malloc_2(N,2);// [N] [reim]
    double **vec=double_malloc_2(N*N,2);
    
    
    int t0=fit_info.t0_GEVP%T;
    if (fit_info.GEVP_swap_t_t0){
        t0=t;
        t=fit_info.t0_GEVP%T;
    }
    else if (fit_info.GEVP_tpt0)
        t0= (fit_info.t0_GEVP+t)%T ; 
    
    double r;
    
    double *s=(double*) malloc(sizeof(double)*N);
    double *s0=(double*) malloc(sizeof(double)*N);
    
    
    
     //t
    int count=0;
    for (int dir=0; dir<3;dir++){
        for (int i=0;i<N;i++){
            for (int k=i;k<N;k++){
                int corr_ik= fit_info.corr_id[count];
                int ik=i+k*N;
                int ki=k+i*N;
                //printf("%d  %g\n",ik,in[j][corr_ik][t][0]);
                M[ik][0]  = in[j][corr_ik][t][0];
                Mt0[ik][0]= in[j][corr_ik][t0][0];
                M[ki][0]=M[ik][0];
                Mt0[ki][0]=Mt0[ik][0];
                count++;
            }
        }
    }
    auto v= fit_info.corr_id;
    M[0][0]/=2.0;
    Mt0[0][0]/=2.0;
    
    M[8][0]/=2.0;
    Mt0[8][0]/=2.0;
    
//     error(!is_it_positive_lex_reim(Mt0, N) , 1, "GEVP_matrix:", "GEVP_matrix M(t0) not positive defined"  ) ;
//   
    //printf("t=%d\n", t);
    generalysed_Eigenproblem(M,Mt0,N,&lambdat,&vec); 
//     if (t==t0){
//         for (int i=0;i<N;i++)
//             for (int j=0;j<N;j++)
//                 printf("%.15f\t",M[i+j*N][0]);
//         printf("\n");
//         for (int i=0;i<N;i++)
//             for (int j=0;j<N;j++)
//                 printf("%.15f\t",Mt0[i+j*N][0]);
//         printf("\n");
//         
//     for (int i=0;i<N;i++)
//         printf("%g\t",lambdat[i][0]);
//     printf("\n");
//     }
    int n=fit_info.n;
    r_equal_value_or_vector(r,  lambdat, vec, fit_info,  t, t0);
    
    free_2(N*N,M);
    free_2(N*N,Mt0);
    free_2(N,lambdat);
    free_2(N*N,vec);
    
    return r;
}



/**********************************
 * you need to fill
 * std::vector<int> fit_info.corr_id
 * with the id of the correlators in the order:
 * [first row], [second row], ...
 * M_00, M_01, ..., M_0N, M_11, ... + M_00(py), M_00(pz), M_22(py), M_22(pz),
**********************************/
double GEVP_matrix_4_p1(int j, double ****in,int t,struct fit_type fit_info ){
    double ct,ctp;
    int N=fit_info.N;
    int ncorr=fit_info.corr_id.size()/3;
//     error(fit_info.corr_id.size()!=12,1,"GEVP_matrix_p1" ," careful populating the GEVP, we to do manually the sum on the directions xyz");
    error(ncorr!=(N*N+N)/2 ,1,"GEVP_matrix_p1",
          "you need to provide (N^2+N)/2 to populate the top triangular matrix NxN:\n  N=%d    ncorr=%d\n",N,ncorr  );
    
    int T=file_head.l0;
    double **M=double_calloc_2(N*N,2);// [NxN] [reim ]
    double **Mt0=double_calloc_2(N*N,2);
    
    double **lambdat=double_malloc_2(N,2);// [N] [reim]
    double **vec=double_malloc_2(N*N,2);
    int t0=fit_info.t0_GEVP%T;
    if (fit_info.GEVP_swap_t_t0){
        t0=t;
        t=fit_info.t0_GEVP%T;
    }
    else if (fit_info.GEVP_tpt0)
        t0= (fit_info.t0_GEVP+t)%T ; 
     
    double r;
    
    double *s=(double*) malloc(sizeof(double)*N);
    double *s0=(double*) malloc(sizeof(double)*N);
    
    
    
     //t
    int count=0;
    for (int dir=0; dir<3;dir++){
        for (int i=0;i<N;i++){
            for (int k=i;k<N;k++){
                int corr_ik= fit_info.corr_id[count];
                int ik=i+k*N;
                int ki=k+i*N;
                //printf("%d  %g\n",ik,in[j][corr_ik][t][0]);
                M[ik][0]  = in[j][corr_ik][t][0];
                Mt0[ik][0]= in[j][corr_ik][t0][0];
                M[ki][0]=M[ik][0];
                Mt0[ki][0]=Mt0[ik][0];
                count++;
            }
        }
    }
    auto v= fit_info.corr_id;
    M[0][0]/=2.0;
    Mt0[0][0]/=2.0;
    
    M[10][0]/=2.0;
    Mt0[10][0]/=2.0;
    
    M[3][0]/=3.0;
    Mt0[3][0]/=3.0;
    
    M[7][0]/=3.0;
    Mt0[7][0]/=3.0;
    
    M[11][0]/=3.0;
    Mt0[11][0]/=3.0;
    
    M[15][0]/=3.0;
    Mt0[15][0]/=3.0;
    
    
//     error(!is_it_positive_lex_reim(Mt0, N) , 1, "GEVP_matrix:", "GEVP_matrix M(t0) not positive defined"  ) ;
//     
    generalysed_Eigenproblem(M,Mt0,N,&lambdat,&vec); 

    int n=fit_info.n;
    if((t-t0)>=0)
        r=lambdat[n][0];
    else 
        r=lambdat[N-1-n][0];
            
    
    free_2(N*N,M);
    free_2(N*N,Mt0);
    free_2(N,lambdat);
//     free_2(N,lambdatp1);
    free_2(N*N,vec);
    
    return r;
}


/**********************************
 * you need to fill
 * std::vector<int> fit_info.corr_id
 * with the id of the correlators in the order:
 * [first row], [second row], ...
 * M_00, M_01, ..., M_0N, M_11, ... + M_00(py), M_00(pz), M_22(py), M_22(pz),
**********************************/
double GEVP_matrix_p11(int j, double ****in,int t,struct fit_type fit_info ){
    double ct,ctp;
    int N=fit_info.N;
    int ncorr=fit_info.corr_id.size()/3;
//     error(fit_info.corr_id.size()!=12,1,"GEVP_matrix_p1" ," careful populating the GEVP, we to do manually the sum on the directions xyz");
    error(ncorr!=(N*N+N)/2 ,1,"GEVP_matrix_p1",
          "you need to provide (N^2+N)/2 to populate the top triangular matrix NxN:\n  N=%d    ncorr=%d\n",N,ncorr  );
    
    int T=file_head.l0;
    double **M=double_calloc_2(N*N,2);// [NxN] [reim ]
    double **Mt0=double_calloc_2(N*N,2);
    
    double **lambdat=double_malloc_2(N,2);// [N] [reim]
    double **vec=double_malloc_2(N*N,2);
    int t0=fit_info.t0_GEVP%T;
    if (fit_info.GEVP_swap_t_t0){
        t0=t;
        t=fit_info.t0_GEVP%T;
    }
    else if (fit_info.GEVP_tpt0)
        t0= (fit_info.t0_GEVP+t)%T ; 
    
    double r;
    
    double *s=(double*) malloc(sizeof(double)*N);
    double *s0=(double*) malloc(sizeof(double)*N);
    
    
    
     //t
    int count=0;
    for (int dir=0; dir<3;dir++){
        for (int i=0;i<N;i++){
            for (int k=i;k<N;k++){
                int corr_ik= fit_info.corr_id[count];
                int ik=i+k*N;
                int ki=k+i*N;
                //printf("%d  %g\n",ik,in[j][corr_ik][t][0]);
                M[ik][0]  = in[j][corr_ik][t][0];
                Mt0[ik][0]= in[j][corr_ik][t0][0];
                M[ki][0]=M[ik][0];
                Mt0[ki][0]=Mt0[ik][0];
                count++;
            }
        }
    }
    auto v= fit_info.corr_id;
   
    
//     error(!is_it_positive_lex_reim(Mt0, N) , 1, "GEVP_matrix:", "GEVP_matrix M(t0) not positive defined"  ) ;
    generalysed_Eigenproblem(M,Mt0,N,&lambdat,&vec); 
    int n=fit_info.n;
    if((t-t0)>=0)
        r=lambdat[n][0];
    else 
        r=lambdat[N-1-n][0];
            
    
    free_2(N*N,M);
    free_2(N*N,Mt0);
    free_2(N,lambdat);
    free_2(N*N,vec);
    
    return r;
}

#endif

