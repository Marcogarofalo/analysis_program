#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"
#include "global_reph.hpp"
#include "resampling.hpp"
#include "read.hpp"
//#include "m_eff.hpp"
#include "gnuplot.hpp"
//#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "rand.hpp"
#include "non_linear_fit.hpp"
//#include "pion.hpp"
//#include "KandD.hpp"
#include "indices.hpp"
#include "continuum_reph.hpp"
#include "fve.hpp"
#include "tower.hpp"
#include "mutils.hpp"

#define r0A 5.31
#define r0B 5.77
#define r0D 7.60

double FA_FV_ChPT(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3];
    
    double FA,FV;

    FA=P[0]*Mpiw0+P[1]*Mpiw0/(w0*w0);
    return FA;
}

double FA_FV_ChPT_NNLO(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3];
    
    double FA,FV;

    FA=P[0]*Mpiw0+ P[1]*pow(Mpiw0,3)+P[2]*pow(Mpiw0,3)*xG +P[3]*Mpiw0/(w0*w0);
    return FA;
}

double FA_FVxg_ChPT_NNLO(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3];
    
    double FA,FV;

    FA=P[0]*Mpiw0+ P[1]*pow(Mpiw0,3)+P[2]*pow(Mpiw0,3)*xG +P[3]*Mpiw0/(w0*w0);
    FA=FA*xG;
    return FA;
}


double FA_FV_ChPT_pole(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3];
    
    double FA,FV;

    FA= P[0]*Mpiw0 / (1. - P[1]*Mpiw0*Mpiw0*(1.-xG))  +P[2]*Mpiw0/(w0*w0);
    return FA;
}
double FA_FV_pole(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0))*Mpiw0 / (1. - (P[1]+P[3]/(w0*w0))*Mpiw0*Mpiw0*(1.-xG)) ;
    return FA;
}

double FA_FV_Pi_poly2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+ (P[1]+P[3]/(w0*w0))*xG*Mpiw0*Mpiw0)*Mpiw0;
    return FA;
}

double FA_FV_Pi_poly3(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+ (P[1]+P[3]/(w0*w0))*xG*Mpiw0*Mpiw0+  P[4]*xG*xG)*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly2_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+ (P[1]+P[3]/(w0*w0))*xG*Mpiw0*Mpiw0+P[4]*p2+P[5]*k2)*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly2_p2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+ (P[1]+P[3]/(w0*w0))*xG*Mpiw0*Mpiw0+P[4]*p2)*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly2_k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+ (P[1]+P[3]/(w0*w0))*xG*Mpiw0*Mpiw0+P[4]*k2)*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly2_p_k(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p=x[9], k=x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+ (P[1]+P[3]/(w0*w0))*xG*Mpiw0*Mpiw0+P[4]*p+P[5]*k)*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly2_M_Pi_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0) +P[4]*Mpiw0*Mpiw0+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*Mpiw0*Mpiw0+P[6]*p2+P[7]*k2)*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly2_M_Pi_p2k2_p4k4(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0) +P[4]*Mpiw0*Mpiw0+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*Mpiw0*Mpiw0+P[6]*p2+P[7]*k2+P[8]*p2*p2+P[9]*k2*k2)*Mpiw0;
    return FA;
}

double FA_FV_Pi_fpi(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) + (P[3]*xi +P[4]/(w0*w0))*(xG))*Mpiw0/fpiw0;
    return FA;
}

double FA_FV_Pi_fpi_aMpi(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0)+ P[5]*Mpiw0*Mpiw0/(w0*w0)  + (P[3]*xi +P[4]/(w0*w0)  )*(xG))*Mpiw0/fpiw0;
    return FA;
}
double FA_FV_Pi_fpi_aMpi_M4(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0)+ P[5]*Mpiw0*Mpiw0/(w0*w0) + P[6]*xi*xi + (P[3]*xi +P[4]/(w0*w0)  )*(xG))*Mpiw0/fpiw0;
    return FA;
}

double FA_FV_Pi_fpi_aMpi_M4_M4x(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0)+ P[5]*Mpiw0*Mpiw0/(w0*w0) + P[6]*xi*xi + (P[3]*xi +P[4]/(w0*w0) + P[7]*Mpiw0*Mpiw0/(w0*w0) + P[8]*xi*xi )*(xG))*Mpiw0/fpiw0;
    return FA;
}

double FA_FV_Pi_fpi_aMpi_M4_M4x_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0)+ P[5]*Mpiw0*Mpiw0/(w0*w0) + P[6]*xi*xi + (P[3]*xi +P[4]/(w0*w0) + P[7]*Mpiw0*Mpiw0/(w0*w0) + P[8]*xi*xi )*(xG))*Mpiw0/fpiw0+ P[9]*p2+P[10]*k2;
    return FA;
}

double FA_Pi_fpi_aMpi_M4_M4x_fix(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA;
    double PP[9];
    PP[0]=0.01434;
    PP[1]=0.01037;
    PP[2]=-0.17329;
    PP[3]=-0.01527;
    PP[4]=-0.01344;
    PP[5]=P[0];
    PP[6]=P[1];
    PP[7]=P[2];
    PP[8]=P[3];
    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA= (PP[0]+PP[1]*xi+PP[2]/(w0*w0)+ PP[5]*Mpiw0*Mpiw0/(w0*w0) + PP[6]*xi*xi + (PP[3]*xi +PP[4]/(w0*w0) + PP[7]*Mpiw0*Mpiw0/(w0*w0) + PP[8]*xi*xi )*(xG))*Mpiw0/fpiw0;
    return FA;
}
double FV_Pi_fpi_aMpi_M4_M4x_fix(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA;
    double PP[9];
    PP[0]=+0.02507;
    PP[1]=+0.14781;
    PP[2]=-0.17806;
    PP[3]=-0.08657;
    PP[4]=+0.00112;
    PP[5]=P[0];
    PP[6]=P[1];
    PP[7]=P[2];
    PP[8]=P[3];
    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA= (PP[0]+PP[1]*xi+PP[2]/(w0*w0)+ PP[5]*Mpiw0*Mpiw0/(w0*w0) + PP[6]*xi*xi + (PP[3]*xi +PP[4]/(w0*w0) + PP[7]*Mpiw0*Mpiw0/(w0*w0) + PP[8]*xi*xi )*(xG))*Mpiw0/fpiw0;
    return FA;
}

double FA_FV_Pi_fpi_aMpix(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0)+ P[5]*Mpiw0*Mpiw0/(w0*w0)  + (P[3]*xi +P[4]/(w0*w0) + P[6]*Mpiw0*Mpiw0/(w0*w0)  )*(xG))*Mpiw0/fpiw0;
    return FA;
}


double FA_FV_Pi_w0(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) + (P[3]*xi +P[4]/(w0*w0))*(xG))*Mpiw0;
    return FA;
}

double FA_FV_Pi_fpi_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) + (P[3]*xi +P[4]/(w0*w0))*(xG))*Mpiw0/fpiw0 + P[5]*p2+P[6]*k2;
    return FA;
}


double FA_FV_Pi_poly2_M_Pi_p2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0) +P[4]*Mpiw0*Mpiw0+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*Mpiw0*Mpiw0+P[6]*p2)*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly2_M_Pi_k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0) +P[4]*Mpiw0*Mpiw0+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*Mpiw0*Mpiw0+P[6]*k2)*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly2_M_Pi_p_k(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p=x[9], k=x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0) +P[4]*Mpiw0*Mpiw0+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*Mpiw0*Mpiw0+P[6]*p+P[7]*k)*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly2_M_Pi_p2k2_kp(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], kp=x[9]*x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0) +P[4]*Mpiw0*Mpiw0+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*Mpiw0*Mpiw0+P[6]*p2 + P[7]*k2 + P[8]*kp)*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly2_M_Pi_p2k2x(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], kp=x[9]*x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0) +P[4]*Mpiw0*Mpiw0+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*Mpiw0*Mpiw0+P[6]*p2 + P[7]*k2/xG )*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly2_M_Pi_p2kx(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], kp=x[9]*x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0) +P[4]*Mpiw0*Mpiw0+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*Mpiw0*Mpiw0+P[6]*p2 + P[7]*sqrt(k2)/xG )*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly2_M_Pi_p2k2x_kx(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10], kp=x[9]*x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0) +P[4]*Mpiw0*Mpiw0+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*Mpiw0*Mpiw0+P[6]*p2 + P[7]*k2/xG+P[8]*sqrt(k2)/xG )*Mpiw0;
    return FA;
}
double FA_FV_Pi_poly3_M_Pi_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0) +P[4]*Mpiw0*Mpiw0+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*Mpiw0*Mpiw0+P[6]*p2+P[7]*k2 +P[8]*xG*xG*Mpiw0*Mpiw0)*Mpiw0;
    return FA;
}
double fpiw0_fit(int n, int Nvar, double *x,int Npar,double  *P){
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3];
    double FA,H,fpw0;
    fpw0=P[4]*(1.+P[5]*Mpiw0*Mpiw0+P[6]/(w0*w0)+P[7]*Mpiw0*Mpiw0*Mpiw0*Mpiw0 );
    
    return fpw0;
}

double HA_pole(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3];
    double FA,H,fpw0;

    FA= FA_FV_pole(n,  Nvar, x,Npar,P);//(P[0]+P[2]/(w0*w0))*Mpiw0 / (1. - (P[1]+P[3]/(w0*w0))*Mpiw0*Mpiw0*(1.-xG)) ;
    fpw0=fpiw0_fit(n,  Nvar, x,Npar,P);//P[4]*(1.+P[5]*Mpiw0*Mpiw0+P[6]/(w0*w0) );
    H=(Mpiw0*xG/2.)  *FA + fpw0;
    return H;
}

double FA_FV_pole_xg_const(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0))*Mpiw0 / (1. - (P[1]+P[3]/(w0*w0))*Mpiw0*Mpiw0*(1.-xG)) +P[4]/(w0*w0*xG);
    return FA;
}

double FA_FV_K_pole(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0)*MKw0 / (1. - (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*MKw0*MKw0*(1.-xG)) ;

    return FA;
}


double FA_FV_K_fK(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) + (P[3]*xi +P[4]/(w0*w0))*(xG))*MKw0/fKw0;
    return FA;
}

double FA_FV_K_fK_Mpi(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    double xiK=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) + (P[3]*xiK +P[4]/(w0*w0))*(xG))*MKw0/fKw0;
    return FA;
}

double FA_FV_K_fK_Mpi_heavy(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    double xiK=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) + (P[3]+P[4]*xi +P[5]/(w0*w0))*(xG))*MKw0/fKw0;
    return FA;
}


double FA_FV_K_fK_Mpi_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    double xiK=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) + (P[3]*xiK +P[4]/(w0*w0))*(xG))*MKw0/fKw0+P[5]*p2+P[6]*k2;;
    return FA;
}

double FA_FV_K_fK_Mpi_heavy_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    double xiK=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) + (P[3]+P[4]*xi +P[5]/(w0*w0))*(xG))*MKw0/fKw0+P[6]*p2+P[7]*k2;;
    return FA;
}

double FA_FV_K_fK_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) + (P[3]*xi +P[4]/(w0*w0))*(xG))*MKw0/fKw0+P[5]*p2+P[6]*k2;
    return FA;
}

double FA_FV_K_fK_MK(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) +P[5]*MKw0*MKw0/(w0*w0) + (P[3]*xi +P[4]/(w0*w0)  +P[6]*MKw0*MKw0/(w0*w0) )*(xG))*MKw0/fKw0;
    return FA;
}

double FA_FV_K_fK_MK_M4(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) +P[5]*MKw0*MKw0/(w0*w0) +P[7]*xi*xi + (P[3]*xi +P[4]/(w0*w0)  +P[6]*MKw0*MKw0/(w0*w0) + P[8]*xi*xi )*(xG))*MKw0/fKw0;
    return FA;
}

double FA_FV_K_fK_MK_M4_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) +P[5]*MKw0*MKw0/(w0*w0) +P[7]*xi*xi + (P[3]*xi +P[4]/(w0*w0)  +P[6]*MKw0*MKw0/(w0*w0) + P[8]*xi*xi )*(xG))*MKw0/fKw0+P[9]*p2+P[10]*k2;
    return FA;
}

double FA_FV_K_fK_Mpi_Mpi4(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    double xiK=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);

    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) +P[5]*Mpiw0*Mpiw0/(w0*w0) +P[7]*xi*xi + (P[3]*xiK +P[4]/(w0*w0)  +P[6]*Mpiw0*Mpiw0/(w0*w0) + P[8]*xi*xi )*(xG))*MKw0/fKw0;
    return FA;
}

double FA_FV_K_fK_Mpi_Mpi4_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    double xiK=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);

    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) +P[5]*Mpiw0*Mpiw0/(w0*w0) +P[7]*xi*xi + (P[3]*xiK +P[4]/(w0*w0)  +P[6]*Mpiw0*Mpiw0/(w0*w0) + P[8]*xi*xi )*(xG))*MKw0/fKw0+P[9]*p2+P[10]*k2;
    return FA;
}

double FA_FV_K_fK_Mpi_Mpi4_heavy(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    double xiK=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);

    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) +P[6]*Mpiw0*Mpiw0/(w0*w0) +P[7]*xi*xi + (P[3]+P[4]*xi +P[5]/(w0*w0)  +P[8]*Mpiw0*Mpiw0/(w0*w0) + P[9]*xi*xi )*(xG))*MKw0/fKw0;
    return FA;
}

double FA_FV_K_fK_Mpi_Mpi4_heavy_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    double xiK=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) +P[6]*Mpiw0*Mpiw0/(w0*w0) +P[7]*xi*xi + (P[3]+P[4]*xi +P[5]/(w0*w0)  +P[8]*Mpiw0*Mpiw0/(w0*w0) + P[9]*xi*xi )*(xG))*MKw0/fKw0+P[10]*p2+P[11]*k2;

    return FA;
}


double FA_FV_K_fK_MK_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11], fKw0=x[12];
    
    double FA,FV;

    double xi=MKw0*MKw0/(16*M_PI*M_PI*fKw0*fKw0);
    FA= (P[0]+P[1]*xi+P[2]/(w0*w0) +P[5]*MKw0*MKw0/(w0*w0) + (P[3]*xi +P[4]/(w0*w0)  +P[6]*MKw0*MKw0/(w0*w0) )*(xG))*MKw0/fKw0+P[7]*p2+P[8]*k2;
    return FA;
}
double FA_FV_K_pole_M_Pi_K(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+ +P[6]*MKw0*MKw0+P[8]*MKw0*MKw0/(w0*w0))*MKw0 / (1. - (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0+P[9]*MKw0*MKw0/(w0*w0))*MKw0*MKw0*(1.-xG)) ;

    return FA;
}

double FA_FV_K_poly2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5];
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*MKw0*MKw0)*MKw0;
    return FA;
}
double FA_FV_K_poly3_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5], p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*MKw0*MKw0  + P[6]*p2+P[7]*k2+P[8]*xG*xG)*MKw0;
    return FA;
}



double FA_FV_K_poly2_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5] ,p2=x[9]*x[9] , k2=x[10]*x[10];
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*MKw0*MKw0+ P[6]*p2+P[7]*k2)*MKw0;
    return FA;
}
double FA_FV_K_poly2_M_Pi_K(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5] ,p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0+P[8]*MKw0*MKw0/(w0*w0))+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0+P[9]*MKw0*MKw0/(w0*w0))*xG*MKw0*MKw0)*MKw0;
    return FA;
}
double FA_FV_K_poly2_M_Pi_K_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5] ,p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0+P[8]*MKw0*MKw0/(w0*w0))+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0+P[9]*MKw0*MKw0/(w0*w0))*xG*MKw0*MKw0+P[10]*p2+P[11] *k2)*MKw0;
    return FA;
}
double FA_FV_K_poly3_M_Pi_K_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5] ,p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0+P[8]*MKw0*MKw0/(w0*w0))+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0+P[9]*MKw0*MKw0/(w0*w0))*xG*MKw0*MKw0+P[10]*p2+P[11] *k2+P[12]*xG*xG)*MKw0;
    return FA;
}

double fKw0_fit(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5];
    double FA,H,fpw0;

    
    fpw0=(P[6]+P[9]*MKw0*MKw0)*(1.+(P[7]+P[10]*MKw0*MKw0)*Mpiw0*Mpiw0+(P[8]+P[11]*MKw0*MKw0)/(w0*w0) );
    return fpw0;
}
double HA_K_pole(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5];
    double FA,H,fpw0;

    FA= FA_FV_K_pole(n,  Nvar, x,Npar,P);
    fpw0=fKw0_fit(n,  Nvar, x,Npar,P);
    H=(MKw0*xG/2.)  *FA + fpw0;
    return H;
}


double FA_FV_K_pole_xg_const(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5];
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0)*MKw0 / (1. - (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*MKw0*MKw0*(1.-xG)) +P[6]/(w0*w0*xG);

    return FA;
}

double FA_FV_flat(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2];
    
    double FA,FV;

    FA=P[0]+P[1]*mw+P[2]/(w0*w0);
    return FA;
}

double FA_FV_flat_x_m(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2];
    
    double FA,FV;

    FA=P[0]+P[1]/(w0*w0);
    return FA;
}

double FA_FV_K_ChPT(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5];
    
    double FA,FV;

    FA=(P[0]+P[2]*MKw0)*Mpiw0+(P[1]+P[3]*MKw0)*Mpiw0/(w0*w0);
    return FA;
}

double FA_FV_K_NNLO_ChPT(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5];
    
    double FA,FV;

    //FA=(P[0]+P[2]*MKw0)*Mpiw0+(P[1]+P[3]*MKw0)*Mpiw0/(w0*w0);
    FA=(P[0]+P[4]*MKw0)*Mpiw0+ (P[1]+P[5]*MKw0)*pow(Mpiw0,3)+(P[2]+P[6]*MKw0)*pow(Mpiw0,3)*xG +(P[3]+P[7]*MKw0)*Mpiw0/(w0*w0);
    return FA;
}


double FA_FV_K_chiral_continuum(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], ms=x[3];
    double FA,FV;

    FA=(P[0]+ms*P[4])+(P[1]+ms*P[5])*mw+(P[2]+ms*P[6])*xG+(P[3]+ms*P[7])/(w0*w0);
    return FA;
    
}
double FA_FV_D_pole(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7]  ;
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0)*MDw0 / (1. - (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*MDw0*MDw0*(1.-xG)) ;

    return FA;
}

double FA_FV_Dphys_poly2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=1+( P[3] +P[4]*xi+P[5] *MDw0*MDw0/(w0*w0))*(1-xG);
    FA*=P[0]*(1+P[1]*xi+P[2]*MDw0*MDw0/(w0*w0));
    return FA;
}

double FA_FV_Dphys_poly2_fixa(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=1+( P[2] +P[3]*xi)*(1-xG);
    FA*=P[0]*(1+P[1]*xi);
    return FA;
}

double FA_FV_Dphys_poly2_fixa_simply(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=P[0]+P[1]*xi+(P[2]+P[3]*xi)*xG;
    return FA;
}
double FA_FV_Dphys_pole_fixa_simply(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=(P[0]+P[1]*xi)/(1+(P[2]+P[3]*xi)*xG);
    return FA;
}

double FA_FV_Dphys_poly2_simply(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=P[0]+P[1]*xi+P[2]*xG+P[3]/(w0*w0);
    return FA;
}


double FA_FV_Dphys_pole_simply(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=(  (P[0]+P[1]*xi)/(1+P[2]*xG)    )+   P[3]/(w0*w0);
    return FA;
}

double FA_FV_Dphys_poly2_simply_Mx(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=P[0]+P[1]*xi+(P[2]+P[3]*xi)*xG+P[4]/(w0*w0);
    return FA;
}

double FA_FV_Dphys_pole_simply_Mx(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=(  (P[0]+P[1]*xi)/(1+(P[2]+P[3]*xi)*xG)    )+   P[4]/(w0*w0);
    return FA;
}

double FA_FV_Dphys_pole_simply_ax(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=(  (P[0]+P[1]*xi)/(1+P[2]*xG)    )+   P[3]/(w0*w0)+ P[4]*xG/(w0*w0);
    return FA;
}

double FA_FV_Dphys_pole_simply_ax_Mx(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=(  (P[0]+P[1]*xi)/(1+(P[2]+P[3]*xi)*xG)    )+   P[4]/(w0*w0)+ P[5]*xG/(w0*w0);
    return FA;
}

double FA_FV_Dphys_pole_simply_apole(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=(  (P[0]+P[1]*xi)/(1+(P[2]+P[4]/(w0*w0))*xG)    )+   P[3]/(w0*w0);
    return FA;
}

double FA_FV_Dphys_poly2_simply_ax(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=P[0]+P[1]*xi+P[2]*xG+P[3]/(w0*w0)+P[4]*xG/(w0*w0);
    return FA;
}

double FA_FV_Dphys_poly2_simply_ax_Mx(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=P[0]+P[1]*xi+(P[2]+P[3]*xi)*xG+P[4]/(w0*w0)+P[5]*xG/(w0*w0);
    return FA;
}

double FA_FV_Dphys_poly2_simply_ax_a4x(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=P[0]+P[1]*xi+P[2]*xG+P[3]/(w0*w0)+P[4]*xG/(w0*w0)+P[5]/(w0*w0*w0*w0)+P[6]*xG/(w0*w0*w0*w0);
    return FA;
}



double FA_FV_Dphys_poly2_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=1+( P[3] +P[4]*xi+P[5] *MDw0*MDw0/(w0*w0))*(1-xG);
    FA*=P[0]*(1+P[1]*xi+P[2]*MDw0*MDw0/(w0*w0));
    FA+=P[6]*p2 +P[7]*k2;
    return FA;
}
double FA_FV_Dphys_pole(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=P[0]*(1+P[1]*xi+P[2]*MDw0*MDw0/(w0*w0));
    FA/=1-( P[3] +P[4]*xi+P[5] *MDw0*MDw0/(w0*w0))*(1-xG);

    return FA;
}
double FA_FV_Dphys_pole_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=P[0]*(1+P[1]*xi+P[2]*MDw0*MDw0/(w0*w0));
    FA/=(1-( P[3] +P[4]*xi+P[5] *MDw0*MDw0/(w0*w0))*(1-xG));
    FA+=P[6]*p2 +P[7]*k2;

    return FA;
}


double FA_FV_D_pole_M_Pi_D(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7]  ;
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MDw0+P[8]*MDw0/(w0*w0))*MDw0 / (1. - (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MDw0+P[9]*MDw0/(w0*w0))*MDw0*MDw0*(1.-xG)) ;

    return FA;
}


double FA_FV_D_poly2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ;
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG )*MDw0;
    return FA;
}

double FA_FV_D_poly3(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ;
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG+   P[6]*xG*xG)*MDw0;
    return FA;
}
double FA_FV_D_poly2_p2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] , p2=x[9]*x[9];
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG+ P[6]*p2)*MDw0;
    return FA;
}

double FA_FV_D_poly2_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] , p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*MDw0*MDw0+ P[6]*p2+P[7]*k2)*MDw0;
    return FA;
}

double FA_FV_D_poly3_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] , p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*xG*MDw0*MDw0+ P[6]*p2+P[7]*k2  +P[8]*xG*xG)*MDw0;
    return FA;
}


double FA_FV_D_poly2_M_Pi_D(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] , p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0 +P[6]*MDw0+P[8]*MDw0/(w0*w0))+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MDw0+P[9]*MDw0/(w0*w0))*xG*MDw0*MDw0)*MDw0;
    return FA;
}

double FA_FV_D_poly2_M_Pi_D_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] , p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;
    FA=FA_FV_D_poly2_M_Pi_D(n,Nvar,x,Npar,P)+(P[10]*p2+P[11]*k2)*MDw0;
    //FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0 +P[6]*MDw0+P[8]*MDw0*MDw0/(w0*w0))+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MDw0+P[9]*MDw0*MDw0/(w0*w0))*xG*MDw0*MDw0+ P[10]*p2+P[11]*k2)*MDw0;
    return FA;
}
double FA_FV_D_poly3_M_Pi_D_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] , p2=x[9]*x[9], k2=x[10]*x[10];
    
    double FA,FV;
    FA=FA_FV_D_poly2_M_Pi_D_p2k2(n,Nvar,x,Npar,P)+(P[12]*xG*xG)*MDw0;
    //FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0 +P[6]*MDw0+P[8]*MDw0*MDw0/(w0*w0))+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MDw0+P[9]*MDw0*MDw0/(w0*w0))*xG*MDw0*MDw0+ P[10]*p2+P[11]*k2)*MDw0;
    return FA;
}

double fDw0_fit(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7]  ;
    double FA,H,fpw0;

    
    fpw0=(P[6]+P[9]*MDw0)*(1.+(P[7]+P[10]*MDw0)*Mpiw0*Mpiw0+(P[8]+P[11]*MDw0)/(w0*w0) );
    return fpw0;
}
double HA_D_pole(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7]  ;
    double FA,H,fpw0;

    FA= FA_FV_D_pole(n,  Nvar, x,Npar,P);
    fpw0=fDw0_fit(n,  Nvar, x,Npar,P);
    H=(MDw0*xG/2.)  *FA + fpw0;
    return H;
}

double FA_FV_D_pole_xg_const(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7]  ;
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0)*MDw0 / (1. - (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0)*MDw0*MDw0*(1.-xG)) +P[6]/(w0*w0*xG);

    return FA;
}


double FA_FV_Ds_pole(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ;
    
    double FA,FV;

    //FA=(P[0]+P[3]*MDw0+P[5]*MKw0)/(xG-1+P[1])+  (P[2]+P[4]*MDw0+P[6]*MKw0)/(w0*w0);
    FA= (P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0+P[8]*MDw0)*MDsw0 / (1. - (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0+P[9]*MDw0)*MDsw0*MDsw0*(1.-xG)) ;
    return FA;
}
double FA_FV_Ds_pole_M_Ds_K_Pi(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ;
    
    double FA,FV;

    //FA=(P[0]+P[3]*MDw0+P[5]*MKw0)/(xG-1+P[1])+  (P[2]+P[4]*MDw0+P[6]*MKw0)/(w0*w0);
    FA= (P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0)*MDsw0 / (1. - (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0)*MDsw0*MDsw0*(1.-xG)) ;
    return FA;
}


double fDsw0_fit(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ;
    double FA,H,fpw0;

    
    fpw0=(P[10]+P[13]*MDw0+P[16]*MKw0*MKw0)*(1.+(P[11]+P[14]*MDw0+P[17]*MKw0*MKw0)*Mpiw0*Mpiw0+(P[12]+P[15]*MDw0+P[18]*MKw0*MKw0)/(w0*w0) );
    return fpw0;
}
double HA_Ds_pole(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ;
    double FA,H,fpw0;

    FA= FA_FV_Ds_pole(n,  Nvar, x,Npar,P);
    fpw0=fDsw0_fit(n,  Nvar, x,Npar,P);
    H=(MDsw0*xG/2.)  *FA + fpw0;
    return H;
}

double FA_FV_Ds_pole_xg_const(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ;
    
    double FA,FV;

    FA= (P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0+P[8]*MDw0)*MDsw0 / (1. - (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0+P[9]*MDw0)*MDsw0*MDsw0*(1.-xG))+P[10]/(w0*w0*xG) ;
    return FA;
}
double FA_FV_Ds_line_xg_const(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ;
    
    double FA,FV;

    FA= (P[0]+P[1]/(w0*w0)+P[2]*Mpiw0*Mpiw0+P[3]*MKw0*MKw0+P[4]*MDw0+P[5]*xG)*MDsw0 +P[6]/(w0*w0*xG) ;
    return FA;
}
double FA_FV_Ds_poly2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ;
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0+P[8]*MDw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0+P[9]*MDw0)*xG*MDsw0*MDsw0)*MDsw0;
    return FA;
}

double FA_FV_Ds_poly2_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ,p2=x[9]*x[9] , k2=x[10]*x[10];;
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0+P[8]*MDw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0+P[9]*MDw0)*xG*MDsw0*MDsw0  +P[10]*p2+P[11]*k2 )*MDsw0;
    return FA;
}


double FA_FV_Ds_poly3_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ,p2=x[9]*x[9] , k2=x[10]*x[10];;
    
    double FA,FV;
    FA=FA_FV_Ds_poly2_p2k2( n,  Nvar, x, Npar, P)+ MDsw0*P[12]*xG*xG;
    
    return FA;
}


double FA_FV_Ds_poly2_aMpi_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ,p2=x[9]*x[9] , k2=x[10]*x[10];;
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0/(w0*w0)+P[6]*MKw0*MKw0+P[8]*MDw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0/(w0*w0)+P[7]*MKw0*MKw0+P[9]*MDw0)*xG*MDsw0*MDsw0  +P[10]*p2+P[11]*k2 )*MDsw0;
    return FA;
}

double FA_FV_Ds_poly2_p2k2_Mpi_overall(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ,p2=x[9]*x[9] , k2=x[10]*x[10];;
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0+P[8]*MDw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0+P[9]*MDw0)*xG*MDsw0*MDsw0  +P[10]*p2+P[11]*k2 )*Mpiw0;
    return FA;
}



double FA_FV_Ds_poly2_M_Ds_K_Pi(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ;
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0)*xG*MDsw0*MDsw0  )*MDsw0;
    return FA;
}
double FA_FV_Ds_poly3_M_Ds_K_Pi(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ;
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0)*xG*MDsw0*MDsw0   + P[8]*xG*xG)*MDsw0;
    return FA;
}
double FA_FV_Ds_poly2_M_Ds_K_Pi_p2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8], p2=x[9]*x[9] ;
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0)*xG*MDsw0*MDsw0 +P[8]*p2  )*MDsw0;
    return FA;
}
double FA_FV_Ds_poly2_M_Ds_K_Pi_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8], p2=x[9]*x[9] , k2=x[10]*x[10];
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*Mpiw0*Mpiw0+P[6]*MKw0*MKw0)+ (P[1]+P[3]/(w0*w0)+P[5]*Mpiw0*Mpiw0+P[7]*MKw0*MKw0)*xG*MDsw0*MDsw0 +P[8]*p2+P[9]*k2  )*MDsw0;
    return FA;
}
double FA_FV_Ds_poly3_M_Ds_K_Pi_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8], p2=x[9]*x[9] , k2=x[10]*x[10];
    
    double FA,FV;
    FA=FA_FV_Ds_poly2_p2k2( n,  Nvar, x, Npar, P)+ MDsw0*P[10]*xG*xG;
    
    return FA;
}
double FA_FV_Ds_poly2_M_Ds_K(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ;
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0)+P[4]*MKw0*MKw0)+ (P[1]+P[3]/(w0*w0)+P[5]*MKw0*MKw0)*xG*MDsw0*MDsw0)*MDsw0;
    return FA;
}

double FA_FV_Ds_poly2_M_Ds(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3],  msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ;
    
    double FA,FV;

    FA= ((P[0]+P[2]/(w0*w0))+ (P[1]+P[3]/(w0*w0))*xG*MDsw0*MDsw0)*MDsw0;
    return FA;
}



double FA_FV_Dsphys_poly2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=1+( P[3] +P[4]*xi+P[5] *MDsw0*MDsw0/(w0*w0))*(1-xG);
    FA*=P[0]*(1+P[1]*xi+P[2]*MDsw0*MDsw0/(w0*w0));
    return FA;
}




double FA_FV_Dsphys_poly2_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=1+( P[3] +P[4]*xi+P[5] *MDsw0*MDsw0/(w0*w0))*(1-xG);
    FA*=P[0]*(1+P[1]*xi+P[2]*MDsw0*MDsw0/(w0*w0));
    FA+=P[6]*p2 +P[7]*k2;
    return FA;
}
double FA_FV_Dsphys_pole(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7] ,MDsw0=x[8] ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=P[0]*(1+P[1]*xi+P[2]*MDsw0*MDsw0/(w0*w0));
    FA/=1-( P[3] +P[4]*xi+P[5] *MDsw0*MDsw0/(w0*w0))*(1-xG);

    return FA;
}
double FA_FV_Dsphys_pole_p2k2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double mw=x[0], w0=x[1], xG=x[2], Mpiw0=x[3], msw=x[4],    MKw0=x[5],  mcw0=x[6], MDw0=x[7],MDsw0=x[8]  ,  p2=x[9]*x[9], k2=x[10]*x[10], fpiw0=x[11];
    
    double FA,FV;

    double xi=Mpiw0*Mpiw0/(16*M_PI*M_PI*fpiw0*fpiw0);
    FA=P[0]*(1+P[1]*xi+P[2]*MDsw0*MDsw0/(w0*w0));
    FA/=1-( P[3] +P[4]*xi+P[5] *MDsw0*MDsw0/(w0*w0))*(1-xG);
    FA+=P[6]*p2 +P[7]*k2;

    return FA;
}


int ***mass_index;



double **fit_to_tif(int Npar,int jack_tot,double **fit)
{
    int j,i;
    double **tif=(double**) malloc(sizeof(double*)*jack_tot); 
    for (j=0;j<jack_tot;j++){
        tif[j]=(double*) malloc(sizeof(double)*Npar);
        for (i=0;i<Npar;i++)
            tif[j][i]=fit[i][j];
    }
    return tif;
}

void free_tif(int jack_tot,double **tif){
    int j;
    for (j=0;j<jack_tot;j++)
        free(tif[j]);
    free(tif);
}

 
void print_chiral_extrapolation(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, struct reph_jack *grephJ, struct header *head,  const char *AV,const char *namefile){
    
    int i,j,ix;
    char namefunA[NAMESIZE],namefunAA[NAMESIZE];
    char namefunB[NAMESIZE];
    char namefunD[NAMESIZE];
    char namefun[NAMESIZE];
    FILE *fA=NULL, *fB=NULL, *fD=NULL,*f=NULL;
    double xg;
   
for(ix=0;ix<24;ix++ ){
    xg=0.05*ix+0.00;
    mysprintf(namefunA,NAMESIZE,"%s/%s_chiral_extrapolation_xg%.2f_A.txt",argv[2],namefile,xg);
    mysprintf(namefunB,NAMESIZE,"%s/%s_chiral_extrapolation_xg%.2f_B.txt",argv[2],namefile,xg);
    mysprintf(namefunD,NAMESIZE,"%s/%s_chiral_extrapolation_xg%.2f_D.txt",argv[2],namefile,xg);
    mysprintf(namefun,NAMESIZE,"%s/%s_chiral_extrapolation_xg%.2f.txt",argv[2],namefile,xg);
    fA=open_file(namefunA,"w+");
    fB=open_file(namefunB,"w+");
    fD=open_file(namefunD,"w+");
    f=open_file(namefun,"w+");
    double h=0.002;
    double **tif=fit_to_tif(fit_info.Npar,jack_tot,fit_out.P);   
    double *r=(double*) malloc(sizeof(double)*jack_tot);
    double *m;
    double **x=(double**) malloc(sizeof(double*)*jack_tot);
    
    int i_mK=index_ensemble_twopt_fit(head[0],0,1,0,0);
    int i_mD=index_ensemble_twopt_fit(head[0],3,0,0,0);
    int i_mDs=index_ensemble_twopt_fit(head[0],3,1,0,0);
    for(j=0;j<jack_tot;j++){
        x[j]=(double*) malloc(sizeof(double)*fit_info.Nvar);
        x[j][0]=phys_point[j][0];//mlw
        x[j][1]=phys_point[j][1];//=1e+10;//r0
        x[j][2]=xg;//=1.;//xg
        x[j][3]=phys_point[jack_tot-1][3];//=result.MpiMeV[j]*result.w0MeV[j];
        x[j][4]=phys_point[j][4];//=result.msw[j];
        x[j][5]=phys_point[j][5];//=result.MKMeV[j]*result.w0MeV[j];
        x[j][6]=phys_point[j][6];//=result.mcw[j];
        x[j][7]=phys_point[j][7];//=result.MDMeV[j]*result.w0MeV[j];
        x[j][8]=phys_point[j][8];//=result.MDsMeV[j]*result.w0MeV[j];
        x[j][9]=0;
        x[j][10]=0;
        x[j][11]=phys_point[j][11];//=result.fpiMeV_exp[j]*result.w0MeV[j];
        x[j][12]=phys_point[j][12];//=result.fpiMeV_exp[j]*result.w0MeV[j];?? fk?
    }
    for (i=0;i<500;i++){
        
        for(j=0;j<jack_tot;j++){
              x[j][1]=r0A;//=1e+10;//r0
              x[j][3]=phys_point[jack_tot-1][3]+((double)i)*h;//=result.MpiMeV[j]*result.w0MeV[j];
              x[j][11]=phys_point[jack_tot-1][11]-(phys_point[jack_tot-1][11] -grephJ[10].Zf_PS[0][jack_tot-1] *grephJ[10].w0[jack_tot-1]  )*((double)i)/389.0;//=result.fpiMeV_exp[j]*result.w0MeV[j];

              r[j]=fit_info.function(fit_info.N,fit_info.Npar,x[j],fit_info.Npar,tif[j]);
        }
        m=mean_and_error_jack_biased(jack_tot,r);
        fprintf(fA,"%g  %g  %g  %g  %g\n",x[jack_tot-1][3],m[0],m[1],x[jack_tot-1][2],x[jack_tot-1][3]/x[jack_tot-1][11]);
        free(m);
        
        for(j=0;j<jack_tot;j++){
              x[j][1]=r0B;//=1e+10;//r0
              x[j][3]=phys_point[jack_tot-1][3]+((double)i)*h;//=result.MpiMeV[j]*result.w0MeV[j];
              x[j][11]=phys_point[jack_tot-1][11]-(phys_point[jack_tot-1][11] -grephJ[6].Zf_PS[0][jack_tot-1] *grephJ[6].w0[jack_tot-1]  )*((double)i)/389.0;//=result.fpiMeV_exp[j]*result.w0MeV[j];
//x[j][11]=phys_point[jack_tot-1][11]-((double)i)*h;//=result.fpiMeV_exp[j]*result.w0MeV[j];

              r[j]=fit_info.function(fit_info.N,fit_info.Npar,x[j],fit_info.Npar,tif[j]);
        }
        m=mean_and_error_jack_biased(jack_tot,r);
        fprintf(fB,"%g  %g  %g  %g  %g\n",x[jack_tot-1][3],m[0],m[1],x[jack_tot-1][2],x[jack_tot-1][3]/x[jack_tot-1][11]);
        free(m);
        
        for(j=0;j<jack_tot;j++){
              x[j][1]=r0D;//=1e+10;//r0
              x[j][3]=phys_point[jack_tot-1][3]+((double)i)*h;//=result.MpiMeV[j]*result.w0MeV[j];
              x[j][11]=phys_point[jack_tot-1][11]-(phys_point[jack_tot-1][11] -grephJ[2].Zf_PS[0][jack_tot-1] *grephJ[2].w0[jack_tot-1]  )*((double)i)/250.0;//=result.fpiMeV_exp[j]*result.w0MeV[j];
              r[j]=fit_info.function(fit_info.N,fit_info.Npar,x[j],fit_info.Npar,tif[j]);
        }
        m=mean_and_error_jack_biased(jack_tot,r);
        fprintf(fD,"%g  %g  %g  %g  %g\n",x[jack_tot-1][3],m[0],m[1],x[jack_tot-1][2],x[jack_tot-1][3]/x[jack_tot-1][11]);
        free(m);
        
        
        for(j=0;j<jack_tot;j++){
              x[j][1]=1e+10;//r0
              x[j][3]=phys_point[jack_tot-1][3]+((double)i)*h;//=result.MpiMeV[j]*result.w0MeV[j];
              x[j][11]=phys_point[jack_tot-1][11];
              r[j]=fit_info.function(fit_info.N,fit_info.Npar,x[j],fit_info.Npar,tif[j]);
        }
        m=mean_and_error_jack_biased(jack_tot,r);
        fprintf(f,"%g  %g  %g  %g  %g\n",x[jack_tot-1][3],m[0],m[1],x[jack_tot-1][2],x[jack_tot-1][3]/x[jack_tot-1][11]);
        free(m);

        
    }
    
    free(r);
    free_tif(jack_tot,tif);
    free_tif(jack_tot,x);
    fclose(fA);fclose(fB);fclose(fD);fclose(f);
}
    
}

void print_continuum_extrapolation(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, const char *AV,const char *namefile){
    
    int i,j;
    char name[NAMESIZE];
    FILE *fc=NULL;
    mysprintf(name,NAMESIZE,"%s/%s_continuum_extrapolation_275MeV_xg0.3.txt",argv[2],namefile);
    fc=open_file(name,"w+");
    double h=0.005;
    double **tif=fit_to_tif(fit_info.Npar,jack_tot,fit_out.P);   
    double *r=(double*) malloc(sizeof(double)*jack_tot);
    double *m;
    double **x=(double**) malloc(sizeof(double*)*jack_tot);
    for(j=0;j<jack_tot;j++){
        x[j]=(double*) malloc(sizeof(double)*fit_info.Nvar);
        x[j][0]=phys_point[j][0];//mlw
        x[j][1]=4.5;//=1e+10;//r0
        x[j][2]=0.3;//=1.;//xg
        x[j][3]=275*0.474/197.326963;//=result.MpiMeV[j]*result.w0MeV[j];275MeV
        x[j][4]=phys_point[j][4];//=result.msw[j];
        x[j][5]=phys_point[j][5];//=result.MKMeV[j]*result.w0MeV[j];
        x[j][6]=phys_point[j][6];//=result.mcw[j];
        x[j][7]=phys_point[j][7];//=result.MDMeV[j]*result.w0MeV[j];
        x[j][8]=phys_point[j][8];//=result.MDsMeV[j]*result.w0MeV[j];
        x[j][9]=0;
        x[j][10]=0;
        x[j][11]=phys_point[j][11];//=result.fpiMeV_exp[j]*result.w0MeV[j];       
        x[j][12]=phys_point[j][12];//=result.fKMeV_exp[j]*result.w0MeV[j];       
    }
    for (i=0;i<1500;i++){
        for(j=0;j<jack_tot;j++){
              x[j][1]=4.5+((double)i)*h;//=1e+10;//r0
              r[j]=fit_info.function(fit_info.N,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
        }
        m=mean_and_error_jack_biased(jack_tot,r);
        fprintf(fc,"%g  %g  %g  %g\n",1./(x[jack_tot-1][1]*x[jack_tot-1][1]),m[0],m[1],x[jack_tot-1][2]);
        free(m);
        
    }
    for (i=1500;i<2500;i++){
        for(j=0;j<jack_tot;j++){
              x[j][1]=4.5+2*pow((double) i,1.5)*h;//=1e+10;//r0
              r[j]=fit_info.function(fit_info.N,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
        }
        m=mean_and_error_jack_biased(jack_tot,r);
        fprintf(fc,"%g  %g  %g  %g\n",1./(x[jack_tot-1][1]*x[jack_tot-1][1]),m[0],m[1],x[jack_tot-1][2]);
        free(m);
        
    }
    
    free(r);
    free_tif(jack_tot,tif);
    free_tif(jack_tot,x);
    fclose(fc);
    
}
void print_chiral_continuum_fit(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, const char *AV,const char *namefile, struct header *head ,struct reph_jack *gJ){
    
    int i,j;
    char name[NAMESIZE];
    FILE *fc=NULL,*fcA=NULL,*fcB=NULL,*fcD=NULL;
    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum.txt",argv[2],namefile);
    fc=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_A.txt",argv[2],namefile);
    fcA=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_B.txt",argv[2],namefile);
    fcB=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_D.txt",argv[2],namefile);
    fcD=open_file(name,"w+");
    
    double h=0.01;
    double **tif=fit_to_tif(fit_info.Npar,jack_tot,fit_out.P);   
    double *r=(double*) malloc(sizeof(double)*jack_tot);
    double *m;
    double **x=(double**) malloc(sizeof(double*)*jack_tot);
    for(j=0;j<jack_tot;j++){
        x[j]=(double*) malloc(sizeof(double)*fit_info.Nvar);
        x[j][0]=phys_point[j][0];//mlw
        x[j][1]=1e+10;//r0
        x[j][2]=0;//=1.;//xg
        x[j][3]=phys_point[j][3];//=result.MpiMeV[j]*result.w0MeV[j];275MeV
        x[j][4]=phys_point[j][4];//=result.msw[j];
        x[j][5]=phys_point[j][5];//=result.MKMeV[j]*result.w0MeV[j];
        x[j][6]=phys_point[j][6];//=result.mcw[j];
        x[j][7]=phys_point[j][7];//=result.MDMeV[j]*result.w0MeV[j];
        x[j][8]=phys_point[j][8];//=result.MDsMeV[j]*result.w0MeV[j];
        x[j][9]=0;
        x[j][10]=0;
        x[j][11]=phys_point[j][11];//=result.fpiMeV_exp[j]*result.w0MeV[j];
        x[j][12]=phys_point[j][12];//=result.fKMeV_exp[j]*result.w0MeV[j];
        
    }
    for (i=0;i<500;i++){
        for(j=0;j<jack_tot;j++){
              x[j][2]=((double)i)*h;//=1e+10;//xG
              r[j]=fit_info.function(fit_info.N,fit_info.Npar,x[j],fit_info.Npar,tif[j]);
        }
        m=mean_and_error_jack_biased(jack_tot,r);
        fprintf(fc,"%g  %g  %g \n",x[jack_tot-1][2],m[0],m[1]);
        free(m);
        
    }
    for (i=0;i<500;i++){
        for(j=0;j<jack_tot;j++){
              x[j][1]=r0A;
              x[j][2]=((double)i)*h;//=1e+10;//xG
              r[j]=fit_info.function(fit_info.N,fit_info.Npar,x[j],fit_info.Npar,tif[j]);
        }
        m=mean_and_error_jack_biased(jack_tot,r);
        fprintf(fcA,"%g  %g  %g \n",x[jack_tot-1][2],m[0],m[1]);
        free(m);
        for(j=0;j<jack_tot;j++){
              x[j][1]=r0B;
              r[j]=fit_info.function(fit_info.N,fit_info.Npar,x[j],fit_info.Npar,tif[j]);
        }
        m=mean_and_error_jack_biased(jack_tot,r);
        fprintf(fcB,"%g  %g  %g \n",x[jack_tot-1][2],m[0],m[1]);
        free(m);
        for(j=0;j<jack_tot;j++){
              x[j][1]=r0D;
              r[j]=fit_info.function(fit_info.N,fit_info.Npar,x[j],fit_info.Npar,tif[j]);
        }
        m=mean_and_error_jack_biased(jack_tot,r);
        fprintf(fcD,"%g  %g  %g \n",x[jack_tot-1][2],m[0],m[1]);
        free(m);
    }
    fclose(fc);
    int ikt=3,iks=1;
   int     i_mDs=index_ensemble_twopt_fit(head[0],ikt,iks,0,0);
   int     i_mD=index_ensemble_twopt_fit(head[0],ikt,0,0,0);
   int     i_mK=index_ensemble_twopt_fit(head[0],0,iks,0,0);
   int     i_mp=index_ensemble_twopt_fit(head[0],0,0,0,0);

    for (int e=0;e<ensembles_reph;e++){
        if (e==0)         mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_D15.48.txt",argv[2],namefile);
        else if (e==1)    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_D20.48.txt",argv[2],namefile);
        else if (e==2)    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_D30.48.txt",argv[2],namefile);
        else if (e==3)    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_B25.32.txt",argv[2],namefile);
        else if (e==4)    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_B35.32.txt",argv[2],namefile);
        else if (e==5)    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_B55.32.txt",argv[2],namefile);
        else if (e==6)    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_B75.32.txt",argv[2],namefile);
        else if (e==7)    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_A30.32.txt",argv[2],namefile);
        else if (e==8)    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_A40.32.txt",argv[2],namefile);
        else if (e==9)    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_A60.24.txt",argv[2],namefile);
        else if (e==10)   mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_A80.24.txt",argv[2],namefile);
        fc=open_file(name,"w+");

        for(iks=0;iks<3;iks++){
        for(ikt=0;ikt<head[e].nk;ikt++){   
 
            i_mDs=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);
            i_mD=index_ensemble_twopt_fit(head[e],ikt,0,0,0);
            i_mK=index_ensemble_twopt_fit(head[e],0,iks,0,0);
            i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
            for(j=0;j<jack_tot;j++){
                    x[j][1]=gJ[e].w0[j];
                    x[j][3]=gJ[e].M_PS[i_mp][j]*gJ[e].w0[j];//=result.MpiMeV[j]*result.w0MeV[j];275MeV
                    x[j][5]=gJ[e].M_PS[i_mK][j]*gJ[e].w0[j];//=result.MKMeV[j]*result.w0MeV[j];
                    x[j][7]=gJ[e].M_PS[i_mD][j]*gJ[e].w0[j];//=result.MDMeV[j]*result.w0MeV[j];
                    x[j][8]=gJ[e].M_PS[i_mDs][j]*gJ[e].w0[j];//=result.MDsMeV[j]*result.w0MeV[j];
                    x[j][11]=gJ[e].Zf_PS[i_mp][j]*gJ[e].w0[j];
                    x[j][12]=gJ[e].Zf_PS[i_mK][j]*gJ[e].w0[j];
            }
            for (i=0;i<500;i++){
                for(j=0;j<jack_tot;j++){
                    x[j][2]=((double)i)*h;//=1e+10;//xG
                    r[j]=fit_info.function(fit_info.N,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
                }
                m=mean_and_error_jack_biased(jack_tot,r);
                fprintf(fc,"%g  %g  %g \n",x[jack_tot-1][2],m[0],m[1]);
                free(m);
            }
            fprintf(fc,"\n\n");
        }}
        fclose(fc); 
    }
    

    free(r);
    free_tif(jack_tot,tif);
    free_tif(jack_tot,x);
    fclose(fcA);    fclose(fcB);    fclose(fcD);
    
}


void print_chiral_continuum_fit_phys(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, const char *AV,const char *namefile, struct header *head ,struct reph_jack *gJ){
    
    int i,j;
    char name[NAMESIZE];
    FILE *fc=NULL,*fcA=NULL,*fcB=NULL,*fcD=NULL;
    
    
    double h=0.01;
    double **tif=fit_to_tif(fit_info.Npar,jack_tot,fit_out.P);   
    double *r=(double*) malloc(sizeof(double)*jack_tot);
    double *m,aMD;
    double **x=(double**) malloc(sizeof(double*)*jack_tot);
    for(j=0;j<jack_tot;j++){
        x[j]=(double*) malloc(sizeof(double)*fit_info.Nvar);
        x[j][0]=phys_point[j][0];//mlw
        x[j][1]=1e+10;//r0
        x[j][2]=0;//=1.;//xg
        x[j][3]=phys_point[j][3];//=result.MpiMeV[j]*result.w0MeV[j];275MeV
        x[j][4]=phys_point[j][4];//=result.msw[j];
        x[j][5]=phys_point[j][5];//=result.MKMeV[j]*result.w0MeV[j];
        x[j][6]=phys_point[j][6];//=result.mcw[j];
        x[j][7]=phys_point[j][7];//=result.MDMeV[j]*result.w0MeV[j];
        x[j][8]=phys_point[j][8];//=result.MDsMeV[j]*result.w0MeV[j];
        x[j][9]=0;
        x[j][10]=0;
        x[j][11]=phys_point[j][11];//=result.fpiMeV_exp[j]*result.w0MeV[j];
        x[j][12]=phys_point[j][12];//=result.fKMeV_exp[j]*result.w0MeV[j];
        
    }
    
   int ikt=3,iks=1;
   int     i_mDs=index_ensemble_twopt_fit(head[0],ikt,iks,0,0);
   int     i_mD=index_ensemble_twopt_fit(head[0],ikt,0,0,0);
   int     i_mK=index_ensemble_twopt_fit(head[0],0,iks,0,0);
   int     i_mp=index_ensemble_twopt_fit(head[0],0,0,0,0);
    int i_mDs1c1, i_mDs1c2, i_mDs2c1, i_mDs2c2, i_mK1, i_mD1;
    int iK1, iK2, iD1, iD2, iDd1c1, iDd2c1, iDd1c2, iDd2c2;

    
    
    for (int e=0;e<ensembles_reph;e++){
        if (e==0)         mysprintf(name,NAMESIZE,"%s/%s_phys_chiral_continuum_D15.48.txt",argv[2],namefile);
        else if (e==1)    mysprintf(name,NAMESIZE,"%s/%s_phys_chiral_continuum_D20.48.txt",argv[2],namefile);
        else if (e==2)    mysprintf(name,NAMESIZE,"%s/%s_phys_chiral_continuum_D30.48.txt",argv[2],namefile);
        else if (e==3)    mysprintf(name,NAMESIZE,"%s/%s_phys_chiral_continuum_B25.32.txt",argv[2],namefile);
        else if (e==4)    mysprintf(name,NAMESIZE,"%s/%s_phys_chiral_continuum_B35.32.txt",argv[2],namefile);
        else if (e==5)    mysprintf(name,NAMESIZE,"%s/%s_phys_chiral_continuum_B55.32.txt",argv[2],namefile);
        else if (e==6)    mysprintf(name,NAMESIZE,"%s/%s_phys_chiral_continuum_B75.32.txt",argv[2],namefile);
        else if (e==7)    mysprintf(name,NAMESIZE,"%s/%s_phys_chiral_continuum_A30.32.txt",argv[2],namefile);
        else if (e==8)    mysprintf(name,NAMESIZE,"%s/%s_phys_chiral_continuum_A40.32.txt",argv[2],namefile);
        else if (e==9)    mysprintf(name,NAMESIZE,"%s/%s_phys_chiral_continuum_A60.24.txt",argv[2],namefile);
        else if (e==10)   mysprintf(name,NAMESIZE,"%s/%s_phys_chiral_continuum_A80.24.txt",argv[2],namefile);
        fc=open_file(name,"w+");
        i_mDs1c1=index_ensemble_twopt_fit(head[e],3,1,0,0);
        i_mDs2c1=index_ensemble_twopt_fit(head[e],3,2,0,0);
        i_mDs1c2=index_ensemble_twopt_fit(head[e],4,1,0,0);
        i_mDs2c2=index_ensemble_twopt_fit(head[e],4,2,0,0);
        i_mD=index_ensemble_twopt_fit(head[e],3,0,0,0);
        i_mD1=index_ensemble_twopt_fit(head[e],4,0,0,0);
        i_mK=index_ensemble_twopt_fit(head[e],0,1,0,0);
        i_mK1=index_ensemble_twopt_fit(head[e],0,2,0,0);
        i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
        for(j=0;j<jack_tot;j++){
                    x[j][0]=head[e].k[head[e].nk+0]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                    x[j][1]=gJ[e].w0[j];//w0
                    x[j][3]=gJ[e].M_PS[i_mp][j]*gJ[e].w0[j];//x_gamma
                        
                    x[j][4]=head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j];//ms*w0
                    
                    aMD=inter_2(head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+2]*gJ[e].w0[j]/gJ[e].Zp[j],
                                gJ[e].M_PS[i_mK][j],                 gJ[e].M_PS[i_mK1][j], 
                                result.msw[j]);
                    x[j][5]=aMD*gJ[e].w0[j];//M_K r0    
                    
                    x[j][6]=head[e].k[head[e].nk+3]*gJ[e].w0[j]/gJ[e].Zp[j];//mc*w0
                    
                    aMD=inter_2(head[e].k[head[e].nk+3]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+4]*gJ[e].w0[j]/gJ[e].Zp[j],
                                gJ[e].M_PS[i_mD][j],                 gJ[e].M_PS[i_mD1][j], 
                                result.mcw[j]);
                    x[j][7]=aMD*gJ[e].w0[j];//M_K r0                  
                    
                    aMD=inter_4(head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+2]*gJ[e].w0[j]/gJ[e].Zp[j],
                                    head[e].k[head[e].nk+3]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+4]*gJ[e].w0[j]/gJ[e].Zp[j],
                                    gJ[e].M_PS[i_mDs1c1][j],                 gJ[e].M_PS[i_mDs2c1][j], 
                                    gJ[e].M_PS[i_mDs1c2][j],                 gJ[e].M_PS[i_mDs2c2][j], 
                                    result.msw[j], result.mcw[j]);
                    x[j][8]=aMD*gJ[e].w0[j];//M_D r0
                    
                    x[j][9]=0;
                    x[j][10]=0;
                    
                    x[j][11]=gJ[e].Zf_PS[i_mp][j] *  gJ[e].w0[j];
                    
                    x[j][12]=inter_2(head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+2]*gJ[e].w0[j]/gJ[e].Zp[j],
                                gJ[e].Zf_PS[i_mK][j]*  gJ[e].w0[j],                 gJ[e].Zf_PS[i_mK1][j]*  gJ[e].w0[j], 
                                result.msw[j]);
         }  
        
        
        
        for(iks=0;iks<3;iks++){
        for(ikt=0;ikt<head[e].nk;ikt++){   
 
             
            //i_mDs=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);
              for (i=0;i<500;i++){   //v_gamma loop
                for(j=0;j<jack_tot;j++){
                      
                    x[j][2]=((double)i)*h;//x_gamma

                    
                   /* if (fabs(x[j][1]-gJ[e].w0[j])>0.1) printf("here 1\n");
                    if (fabs(x[j][2]-((double)i)*h)>0.1) printf("here 2\n");
                    if (fabs(x[j][3]-gJ[e].M_PS[i_mp][j]*gJ[e].w0[j])>0.1) printf("here 3\n");
                    if (fabs(x[j][5]-gJ[e].M_PS[i_mK][j]*gJ[e].w0[j])>0.1) printf("here 4\n");
                    if (fabs(x[j][7]-gJ[e].M_PS[i_mD][j]*gJ[e].w0[j])>0.1) printf("here 5\n");
                    //if (fabs(x[j][8]-gJ[e].M_PS[i_mDs][j]*gJ[e].w0[j])>1) printf("here 6\n");
                    if (fabs(x[j][11]-gJ[e].Zf_PS[i_mp][j]*gJ[e].w0[j])>0.1) printf("here 7\n");
                    if (fabs(x[j][12]-gJ[e].Zf_PS[i_mK][j]*gJ[e].w0[j])>0.1) printf("here 8\n");*/
                    r[j]=fit_info.function(fit_info.N,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
                    // if (strcmp(namefile,"FA_H_K_fK")==0 && i==499 ) printf("xg= %g  %g \n" , x[j][2], r[j] );
                    //if (i==499  && j==0) printf( "fk=%g\n", x[j][12]);
                    //if (strcmp(namefile,"FA_H_K_fK")==0 && i==499 && fabs(r[j]+1.7 )>6) printf("error j=%d   r=%g   fk=%g\n",j,r[j],x[j][12]);
                }
                m=mean_and_error_jack_biased(jack_tot,r);
                fprintf(fc,"%g  %g  %g \n",x[jack_tot-1][2],m[0],m[1]);
                free(m);
            }
            fprintf(fc,"\n\n");
        }}
        fclose(fc); 
    }
    

    free(r);
    free_tif(jack_tot,tif);
    free_tif(jack_tot,x);
    
}



void reweighting_for_plot(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, const char *AV,const char *namefile,struct header *head , struct reph_jack *gJ ){

    int i,j;
    char name[NAMESIZE];
    FILE *fc=NULL;
    mysprintf(name,NAMESIZE,"%s/%s_reweighting_275MeV_xg0.3.txt",argv[2],namefile);
    fc=open_file(name,"w+");
    int i_mD,i_mDs,i_mK,i_mp,e,ikt,iks;
    double h=0.005;
    double **tif=fit_to_tif(fit_info.Npar,jack_tot,fit_out.P);   

    double *r=(double*) malloc(sizeof(double)*jack_tot);
    double *r1=(double*) malloc(sizeof(double)*jack_tot);
    double *r2=(double*) malloc(sizeof(double)*jack_tot);
    double *m,*m1,*m2;
    double **x=(double**) malloc(sizeof(double*)*jack_tot);
    
    double **interpolation_FX;

                            

    for(j=0;j<jack_tot;j++){
        x[j]=(double*) malloc(sizeof(double)*fit_info.Nvar);
        x[j][0]=phys_point[j][0];//mlw
        x[j][2]=0.3;//=1.;//xg
        x[j][3]=275*0.474/197.326963;//=result.MpiMeV[j]*result.w0MeV[j];275MeV
        x[j][4]=phys_point[j][4];//=result.msw[j];
        x[j][5]=phys_point[j][5];//=result.MKMeV[j]*result.w0MeV[j];
        x[j][6]=phys_point[j][6];//=result.mcw[j];
        x[j][7]=phys_point[j][7];//=result.MDMeV[j]*result.w0MeV[j];
        x[j][8]=phys_point[j][8];//=result.MDsMeV[j]*result.w0MeV[j];
        x[j][9]=0;
        x[j][10]=0;
        x[j][11]=phys_point[j][11];//=result.fpiMeV_exp[j]*result.w0MeV[j];        
        x[j][12]=phys_point[j][12];//=result.fKMeV_exp[j]*result.w0MeV[j];        
    }
    for (e=0;e<ensembles_reph;e++){
    for(iks=1;iks<2;iks++){
    for(ikt=3;ikt<4;ikt++){     

        i_mDs=index_ensemble_twopt_fit(head[e],3,1,0,0);
        i_mD=index_ensemble_twopt_fit(head[e],3,0,0,0);
        i_mK=index_ensemble_twopt_fit(head[e],0,1,0,0);
        i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
        if(strcmp(AV,"A")==0  ){
            //interpolation_FX=(double**) malloc(sizeof(double*)*1);
            //interpolation_FX[0]=(double*) malloc(sizeof(double)*1);
            interpolation_FX=interpolate_FX(argv, namefile, ikt, iks,gJ[e].xG,0.285,0.315,gJ[e].FA,  head[e] , jack_tot ,1/*npar_fun,double*/,polynimial_degree_n/*poles_degree_n*/ );
        }
        else if(strcmp(AV,"V")==0){
            interpolation_FX=interpolate_FX(argv, namefile, ikt, iks,gJ[e].xG,0.285,0.315,gJ[e].FV,  head[e] , jack_tot ,1/* npar_fun,double*/,polynimial_degree_n);
        }
        else if(strcmp(AV,"{A^H}")==0 || strcmp(AV,"{V^H}")==0  || strcmp(AV,"{V^HA}")==0){
            //interpolation_FX=interpolate_FX(argv, namefile, ikt, iks,gJ[e].xG,0.285,0.315,gJ[e].FA_from_H0,  head[e] , jack_tot ,1/* npar_fun,double*/,polynimial_degree_n);
            interpolation_FX=subtract_non_Lorentz_and_int(argv, namefile, ikt, iks,gJ[e].xG,0.27,0.33,  head , jack_tot ,1/*npar_fun,double*/,polynimial_degree_n/*poles_degree_n*/,fit_out, fit_info  , gJ, e,AV );
        }
        
        else error(0==0,1,"argument of function reweighting_for_plot AV is nor A neither V","");

        for(j=0;j<jack_tot;j++){
              x[j][1]=gJ[e].w0[j];//=1e+10;//r0
              x[j][3]=275*0.474/197.326963;//=result.MpiMeV[j]*result.w0MeV[j];275MeV
              x[j][5]=phys_point[j][5];//=result.MKMeV[j]*result.w0MeV[j];
              x[j][7]=phys_point[j][7];//=result.MDMeV[j]*result.w0MeV[j];
              x[j][8]=phys_point[j][8];//=result.MDsMeV[j]*result.w0MeV[j];
              r[j]=fit_info.function(fit_info.N,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
              
              x[j][3]=gJ[e].M_PS[i_mp][j]*gJ[e].w0[j];//M_\pi r0
              x[j][5]=gJ[e].M_PS[i_mK][j]*gJ[e].w0[j];//M_K r0
              x[j][7]=gJ[e].M_PS[i_mD][j]*gJ[e].w0[j];//M_D r0
              x[j][8]=gJ[e].M_PS[i_mDs][j]*gJ[e].w0[j];//M_Ds r0
              r[j]=r[j]/fit_info.function(fit_info.N,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
              
              if(strcmp(AV,"A")==0){
                            r1[j]=interpolation_FX[0][j]   *gJ[e].ZV[j];
              }
              else if(strcmp(AV,"V")==0){
                            r1[j]=interpolation_FX[0][j]   *gJ[e].ZA[j];
              }
              else if(strcmp(AV,"{A^H}")==0){
                            r1[j]=interpolation_FX[0][j] ;
              }
              else if(strcmp(AV,"{V^H}")==0){
                            r1[j]=interpolation_FX[0][j]  ;
              }
              else if(strcmp(AV,"{V^HA}")==0){
                            r1[j]=interpolation_FX[0][j]  ;
              }
              else error(0==0,1,"argument of function reweighting_for_plot AV is nor A neither V","");
              
              
              r2[j]=r[j]*r1[j];
        }
        m=mean_and_error_jack_biased(jack_tot,r);
        m1=mean_and_error_jack_biased(jack_tot, r1);
        m2=mean_and_error_jack_biased(jack_tot, r2);
        fprintf(fc,"%d  %d   %d  %g  %g  %g\t\t    %g  %g\t\t %g  %g\n",e,ikt,iks,x[jack_tot-1][1]/*r0*/,m2[0],m2[1],m1[0],m1[1],x[jack_tot-1][2]/*xg*/,x[jack_tot-1][3]/(0.474/197.326963)/*Mpi*/);
        free(m);free(m1);free(m2);
        free_tif(1,interpolation_FX);
    }}
    }
    free(r);free(r1);free(r2);
    free_tif(jack_tot,tif); 

   
    free_tif(jack_tot,x);

    fclose(fc);        
    
}

void interpolating_for_plot(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, const char *AV,const char *namefile,struct header *head , struct reph_jack *gJ ){

    int i,j;
    char name[NAMESIZE];
    FILE *fc=NULL,*fc_K=NULL,*fc_Ds;
    mysprintf(name,NAMESIZE,"%s/%s_Pi_interpolation_xg0.3.txt",argv[2],namefile);
    fc=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/%s_K_interpolation_xg0.3.txt",argv[2],namefile);
    fc_K=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/%s_Ds_interpolation_xg0.3.txt",argv[2],namefile);
    fc_Ds=open_file(name,"w+");
    int i_mD,i_mDs,i_mK,i_mp,e,ikt,iks;
    double h=0.005;
    double **tif=fit_to_tif(fit_info.Npar,jack_tot,fit_out.P);   

    double *r=(double*) malloc(sizeof(double)*jack_tot);
    double *r1=(double*) malloc(sizeof(double)*jack_tot);
    double *r2=(double*) malloc(sizeof(double)*jack_tot);
    double *m,*m1,*m2;
    double **x=(double**) malloc(sizeof(double*)*jack_tot);
    
    double **interpolation_FX;

    double ave_xg=0.3, dxg=1;
    double xg_min=0.27;
    double xg_max=0.33;

    for(j=0;j<jack_tot;j++){
        x[j]=(double*) malloc(sizeof(double)*fit_info.Nvar);
        x[j][0]=phys_point[j][0];//mlw
        x[j][2]=0.3;//=1.;//xg
        x[j][3]=275*0.474/197.326963;//=result.MpiMeV[j]*result.w0MeV[j];275MeV
        x[j][4]=phys_point[j][4];//=result.msw[j];
        x[j][5]=phys_point[j][5];//=result.MKMeV[j]*result.w0MeV[j];
        x[j][6]=phys_point[j][6];//=result.mcw[j];
        x[j][7]=phys_point[j][7];//=result.MDMeV[j]*result.w0MeV[j];
        x[j][8]=phys_point[j][8];//=result.MDsMeV[j]*result.w0MeV[j];
        x[j][9]=0;
        x[j][10]=0;
        x[j][11]=phys_point[j][11];//=result.fpiMeV_exp[j]*result.w0MeV[j];
        x[j][12]=phys_point[j][12];//=result.fKMeV_exp[j]*result.w0MeV[j];
    }
            fprintf(fc,"\n\n");
    for (e=0;e<ensembles_reph;e++){
        if(e==3 || e==7) {fprintf(fc,"\n\n");fprintf(fc_K,"\n\n");fprintf(fc_Ds,"\n\n");}
        i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
 
        interpolation_FX=subtract_non_Lorentz_and_int(argv, namefile, 0/*ikt*/, 0,gJ[e].xG,xg_min,xg_max,  head , jack_tot ,1/*npar_fun,double*/,polynimial_degree_n/*poles_degree_n*/,fit_out,  fit_info  , gJ, e,AV );
        m=mean_and_error_jack_biased(jack_tot,interpolation_FX[0]);
        fprintf(fc,"%g  %g   %g   %g\n",(xg_min+xg_max)/2., m[0],m[1],gJ[e].M_PS[i_mp][jack_tot-1]*gJ[e].w0[jack_tot-1]);
        free(m);
        free_tif(1,interpolation_FX);
        
        interpolation_FX=subtract_non_Lorentz_and_int_K(argv, namefile,phys_point,xg_min,xg_max,  head , jack_tot ,1/*npar_fun,double*/,polynimial_degree_n/*poles_degree_n*/,fit_out,  fit_info  , gJ, e,AV );
        m=mean_and_error_jack_biased(jack_tot,interpolation_FX[0]);
        //printf(      "%g  %g   %g   %g\n",(xg_min+xg_max)/2., m[0],m[1],gJ[e].M_PS[i_mp][jack_tot-1]*gJ[e].w0[jack_tot-1]);
        fprintf(fc_K,"%g  %g   %g   %g\n",(xg_min+xg_max)/2., m[0],m[1],gJ[e].M_PS[i_mp][jack_tot-1]*gJ[e].w0[jack_tot-1]);
        free(m);
        free_tif(1,interpolation_FX);
        
        interpolation_FX=subtract_non_Lorentz_and_int_Ds(argv, namefile,phys_point,0.28,0.32,  head , jack_tot ,1/*npar_fun,double*/,polynimial_degree_n/*poles_degree_n*/,fit_out,  fit_info  , gJ, e,AV );
        m=mean_and_error_jack_biased(jack_tot,interpolation_FX[0]);
        //printf(      "%g  %g   %g   %g\n",(xg_min+xg_max)/2., m[0],m[1],gJ[e].M_PS[i_mp][jack_tot-1]*gJ[e].w0[jack_tot-1]);
        fprintf(fc_Ds,"%g  %g   %g   %g\n",(xg_min+xg_max)/2., m[0],m[1],gJ[e].M_PS[i_mp][jack_tot-1]*gJ[e].w0[jack_tot-1]);
        free(m);
        free_tif(1,interpolation_FX);
        
        
    }        
      
    free(r);free(r1);free(r2);
    free_tif(jack_tot,tif); 

    free_tif(jack_tot,x);
    fclose(fc);  fclose(fc_K);    fclose(fc_Ds);    

}



void print_data_without_p2k2(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, const char *AV,const char *namefile,struct header *head , struct reph_jack *gJ ){

    int i,j;
    char name[NAMESIZE];
    FILE *fc=NULL;
    int i_mD,i_mDs,i_mK,i_mp,e,ikt,iks;
    double *m,cp2,ck2;
    int Neff=0;
    int imom0,imoms,imomt;
    double *xp=(double*) malloc(sizeof(double)*fit_info.Nvar); 
    double *F=(double*) malloc(sizeof(double)*jack_tot);
    int Njack=jack_tot;
    double xg_min=0.001;
            
    for (e=0;e<ensembles_reph;e++){
        if (e==0)         mysprintf(name,NAMESIZE,"%s/D15.48/out/%s_%s.txt",argv[4],namefile,AV);
        else if (e==1)    mysprintf(name,NAMESIZE,"%s/D20.48/out/%s_%s.txt",argv[4],namefile,AV);
        else if (e==2)    mysprintf(name,NAMESIZE,"%s/D30.48/out/%s_%s.txt",argv[4],namefile,AV);
        else if (e==3)    mysprintf(name,NAMESIZE,"%s/B25.32/out/%s_%s.txt",argv[4],namefile,AV);
        else if (e==4)    mysprintf(name,NAMESIZE,"%s/B35.32/out/%s_%s.txt",argv[4],namefile,AV);
        else if (e==5)    mysprintf(name,NAMESIZE,"%s/B55.32/out/%s_%s.txt",argv[4],namefile,AV);
        else if (e==6)    mysprintf(name,NAMESIZE,"%s/B75.32/out/%s_%s.txt",argv[4],namefile,AV);
        else if (e==7)    mysprintf(name,NAMESIZE,"%s/A30.32/out/%s_%s.txt",argv[4],namefile,AV);
        else if (e==8)    mysprintf(name,NAMESIZE,"%s/A40.32/out/%s_%s.txt",argv[4],namefile,AV);
        else if (e==9)    mysprintf(name,NAMESIZE,"%s/A60.24/out/%s_%s.txt",argv[4],namefile,AV);
        else if (e==10)    mysprintf(name,NAMESIZE,"%s/A80.24/out/%s_%s.txt",argv[4],namefile,AV);
        
        fc=open_file(name,"w+");
 
        for(iks=0;iks<3;iks++){
        for(ikt=0;ikt<head[e].nk;ikt++){   
            i_mDs=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);
            i_mD=index_ensemble_twopt_fit(head[e],ikt,0,0,0);
            i_mK=index_ensemble_twopt_fit(head[e],0,iks,0,0);
            i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
            for(imom0=0;imom0<head[e].nmoms;imom0++){       
            for(imomt=0;imomt<head[e].nmoms;imomt++){
            for(imoms=0;imoms<head[e].nmoms;imoms++){
                i=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);
                if (gJ[e].xG[i][Njack-1]>xg_min){
                for (j=0;j<Njack;j++){ 
                xp[0]=head[e].k[head[e].nk+ikt]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                xp[1]=gJ[e].w0[j];//w0
                xp[2]=gJ[e].xG[i][j];//x_gamma
                xp[3]=gJ[e].M_PS[i_mp][j]*gJ[e].w0[j];//x_gamma
                    
                xp[4]=head[e].k[head[e].nk+iks]*gJ[e].w0[j]/gJ[e].Zp[j];//ms*w0
                xp[5]=gJ[e].M_PS[i_mK][j]*gJ[e].w0[j];//M_K r0
                        
                xp[6]=head[e].k[head[e].nk+ikt]*gJ[e].w0[j]/gJ[e].Zp[j];//mc*w0
                xp[7]=gJ[e].M_PS[i_mD][j]*gJ[e].w0[j];//M_D r0
                xp[8]=gJ[e].M_PS[i_mDs][j]*gJ[e].w0[j];//M_Ds r0
                
                xp[9]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3]);
                xp[9]=2*sin(xp[9]/2.);
                
                xp[10]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3]);
                xp[10]=2*sin(xp[10]/2.);
                
                xp[11]=gJ[e].Zf_PS[i_mp][j] *  gJ[e].w0[j];
                xp[12]=gJ[e].Zf_PS[i_mK][j] *  gJ[e].w0[j];
                
                    cp2=0;//derivative_xi(9, j , fit_out,  fit_info, xp);
                    ck2=0;//derivative_xi(10, j , fit_out,  fit_info, xp);
                    if(strcmp(AV,"A")==0)  F[j]=gJ[e].FA[i][j] * gJ[e].ZV[j]; 
                    else if(strcmp(AV,"V")==0) F[j]=gJ[e].FV[i][j] * gJ[e].ZA[j];
                    else if(strcmp(AV,"{A^H}")==0) F[j]=gJ[e].FA_from_H0[i][j] ;
                    else if(strcmp(AV,"{V^H}")==0) F[j]=gJ[e].FV_from_H0[i][j] * gJ[e].ZA[j];
                    else if(strcmp(AV,"{V^HA}")==0) F[j]=gJ[e].FV_from_H0_HA[i][j] * gJ[e].ZAV[j];
                    else error(0==0,1,"argument of function  AV is nor A neither V","");
                    F[j]=F[j]-cp2*xp[9]*xp[9]-ck2*xp[10]*xp[10];
                }
                m=mean_and_error("jack",jack_tot,F);
                fprintf(fc,"%g  %g   %g   %g  \t %g  %g  %g  %g\n",xp[2],xp[2] ,m[0],m[1],xp[3],xp[5],xp[7],xp[8]);
                free(m);
                }
            }}}
            fprintf(fc,"\n\n");
        }}
        fclose(fc);

    }
    
    free(xp);free(F);
    
}

void print_data_without_p2k2_phys(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, const char *AV,const char *namefile,struct header *head , struct reph_jack *gJ ){

    int i,j;
    char name[NAMESIZE];
    FILE *fc=NULL;
    int i_mD,i_mDs,i_mK,i_mp,e,ikt,iks;
    double *m,cp2,ck2;
   
    int Neff=0;
    int imom0,imoms,imomt;
    double **xp=double_malloc_2(jack_tot,fit_info.Nvar);//(double*) malloc(sizeof(double)*fit_info.Nvar); 
    double *F=(double*) malloc(sizeof(double)*jack_tot);
    int Njack=jack_tot;
    double xg_min=0.001,aMD;


    int i_mDs1c1, i_mDs1c2, i_mDs2c1, i_mDs2c2, i_mK1, i_mD1;
    int iK1, iK2, iD1, iD2, iDs1c1, iDs2c1, iDs1c2, iDs2c2;
    for (e=0;e<ensembles_reph;e++){
        if (e==0)         mysprintf(name,NAMESIZE,"%s/D15.48/out/%s_phys_%s.txt",argv[4],namefile,AV);
        else if (e==1)    mysprintf(name,NAMESIZE,"%s/D20.48/out/%s_phys_%s.txt",argv[4],namefile,AV);
        else if (e==2)    mysprintf(name,NAMESIZE,"%s/D30.48/out/%s_phys_%s.txt",argv[4],namefile,AV);
        else if (e==3)    mysprintf(name,NAMESIZE,"%s/B25.32/out/%s_phys_%s.txt",argv[4],namefile,AV);
        else if (e==4)    mysprintf(name,NAMESIZE,"%s/B35.32/out/%s_phys_%s.txt",argv[4],namefile,AV);
        else if (e==5)    mysprintf(name,NAMESIZE,"%s/B55.32/out/%s_phys_%s.txt",argv[4],namefile,AV);
        else if (e==6)    mysprintf(name,NAMESIZE,"%s/B75.32/out/%s_phys_%s.txt",argv[4],namefile,AV);
        else if (e==7)    mysprintf(name,NAMESIZE,"%s/A30.32/out/%s_phys_%s.txt",argv[4],namefile,AV);
        else if (e==8)    mysprintf(name,NAMESIZE,"%s/A40.32/out/%s_phys_%s.txt",argv[4],namefile,AV);
        else if (e==9)    mysprintf(name,NAMESIZE,"%s/A60.24/out/%s_phys_%s.txt",argv[4],namefile,AV);
        else if (e==10)    mysprintf(name,NAMESIZE,"%s/A80.24/out/%s_phys_%s.txt",argv[4],namefile,AV);
        
        fc=open_file(name,"w+");
        fprintf(fc,"#xg  dxg   F   dF  Mpi  MK MD MDs  P^2  K^2 \n" );
        i_mDs1c1=index_ensemble_twopt_fit(head[e],3,1,0,0);
        i_mDs2c1=index_ensemble_twopt_fit(head[e],3,2,0,0);
        i_mDs1c2=index_ensemble_twopt_fit(head[e],4,1,0,0);
        i_mDs2c2=index_ensemble_twopt_fit(head[e],4,2,0,0);
        i_mD=index_ensemble_twopt_fit(head[e],3,0,0,0);
        i_mD1=index_ensemble_twopt_fit(head[e],4,0,0,0);
        i_mK=index_ensemble_twopt_fit(head[e],0,1,0,0);
        i_mK1=index_ensemble_twopt_fit(head[e],0,2,0,0);
        i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
        for (j=0;j<Njack;j++){
            xp[j][0]=head[e].k[head[e].nk+0]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
            xp[j][1]=gJ[e].w0[j];//w0
            xp[j][3]=gJ[e].M_PS[i_mp][j]*gJ[e].w0[j];//
            xp[j][4]=head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j];//ms*w0
            aMD=inter_2(head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+2]*gJ[e].w0[j]/gJ[e].Zp[j],
                            gJ[e].M_PS[i_mK][j],                 gJ[e].M_PS[i_mK1][j], 
                            result.msw[j]);
            xp[j][5]=aMD*gJ[e].w0[j];//M_K r0    
                
            xp[j][6]=head[e].k[head[e].nk+3]*gJ[e].w0[j]/gJ[e].Zp[j];//mc*w0
            aMD=inter_2(head[e].k[head[e].nk+3]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+4]*gJ[e].w0[j]/gJ[e].Zp[j],
                            gJ[e].M_PS[i_mD][j],                 gJ[e].M_PS[i_mD1][j], 
                            result.mcw[j]);
            
            if (j==Njack-1 && e==0   || j==Njack-1 && e==3 || j==Njack-1 && e==7   ){
                printf("%s MD\n",name);
            if (e==0)
                printf("res=%f  %d\n",aMD/0.0619*197,e  );
            if (e==3)
                printf("res=%f  %d\n",aMD/0.0815*197,e  );
            if (e==7)
                printf("res=%f  %d\n",aMD/0.0885*197,e  );
            }
            xp[j][7]=aMD*gJ[e].w0[j];//M_K r0                  xp[8]=gJ[e].M_PS[i_mDs][j]*gJ[e].w0[j];//M_Ds r0
            aMD=inter_4(head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+2]*gJ[e].w0[j]/gJ[e].Zp[j],
                                head[e].k[head[e].nk+3]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+4]*gJ[e].w0[j]/gJ[e].Zp[j],
                                gJ[e].M_PS[i_mDs1c1][j],                 gJ[e].M_PS[i_mDs2c1][j], 
                                gJ[e].M_PS[i_mDs1c2][j],                 gJ[e].M_PS[i_mDs2c2][j], 
                                result.msw[j], result.mcw[j]);
            if (j==Njack-1 && e==0   || j==Njack-1 && e==3 || j==Njack-1 && e==7   ){
                printf("%s\n",name);
            printf("mus1=%f    mus2=%f   muc1=%f   muc2=%f  \n",head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+2]*gJ[e].w0[j]/gJ[e].Zp[j],
                                head[e].k[head[e].nk+3]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+4]*gJ[e].w0[j]/gJ[e].Zp[j]);
            printf("inter to     %f    %f\n",result.msw[j], result.mcw[j]);
            printf("Ms1c1=%f   Ms2c1=%f  Ms1c2=%f   Ms2c2=%f\n",gJ[e].M_PS[i_mDs1c1][j],   gJ[e].M_PS[i_mDs2c1][j],     gJ[e].M_PS[i_mDs1c2][j],     gJ[e].M_PS[i_mDs2c2][j]);
            printf("res=%f \n",aMD );
            if (e==0)
                printf("res=%f  %d\n",aMD/0.0619*197,e  );
            if (e==3)
                printf("res=%f  %d\n",aMD/0.0815*197,e  );
            if (e==7)
                printf("res=%f  %d\n",aMD/0.0885*197,e  );
            
            
            
            }
            xp[j][8]=aMD*gJ[e].w0[j];//M_D r0        
            xp[j][11]=gJ[e].Zf_PS[i_mp][j] *  gJ[e].w0[j];
               
            xp[j][12]=inter_2(head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+2]*gJ[e].w0[j]/gJ[e].Zp[j],
                            gJ[e].Zf_PS[i_mK][j]*  gJ[e].w0[j],                 gJ[e].Zf_PS[i_mK1][j]*  gJ[e].w0[j], 
                            result.msw[j]);
            
            
        
        }
        
        for(iks=0;iks<3;iks++){
        for(ikt=0;ikt<head[e].nk;ikt++){   
            i_mDs=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);
            for(imom0=0;imom0<head[e].nmoms;imom0++){       
            for(imomt=0;imomt<head[e].nmoms;imomt++){
            for(imoms=0;imoms<head[e].nmoms;imoms++){
                i=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);
                iK1=index_ensemble_twopt_G_fit(head[e],0,1,imom0,imomt,imoms);
                iK2=index_ensemble_twopt_G_fit(head[e],0,2,imom0,imomt,imoms);
                iD1=index_ensemble_twopt_G_fit(head[e],3,0,imom0,imomt,imoms);
                iD2=index_ensemble_twopt_G_fit(head[e],4,0,imom0,imomt,imoms);
                iDs1c1=index_ensemble_twopt_G_fit(head[e],3,1,imom0,imomt,imoms);
                iDs2c1=index_ensemble_twopt_G_fit(head[e],3,2,imom0,imomt,imoms);
                iDs1c2=index_ensemble_twopt_G_fit(head[e],4,1,imom0,imomt,imoms);
                iDs2c2=index_ensemble_twopt_G_fit(head[e],4,2,imom0,imomt,imoms);
                if (gJ[e].xG[i][Njack-1]>xg_min){
                for (j=0;j<Njack;j++){
                
                xp[j][2]=gJ[e].xG[i][j];//x_gamma
   
                
                
                double p=2.*sin(  (2.*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3])/2.);
                double k=2.*sin(  (2.*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3])/2.);

                if(iks==0   && ikt ==0)   
                    aMD=gJ[e].M_PS[i_mp][j];
                else if (iks==1 && ikt ==0  || iks==2 && ikt ==0 )
                    aMD=xp[j][5]/gJ[e].w0[j];
                else if (iks==0 && ikt ==3  || iks==0 && ikt ==4 )
                    aMD=xp[j][7]/gJ[e].w0[j];
                else if (iks==1 && ikt ==3  || iks==1 && ikt ==4 || iks==2 && ikt ==3  || iks==2 && ikt ==4 )    
                    aMD=xp[j][8]/gJ[e].w0[j];//M_D r0
                    
                double E_p=aMD*aMD;
                E_p+=p*p;
                E_p=sqrt(E_p);
                double E_k=2.*asinh(   sqrt(k*k)/2.   );
                double kdp=E_p*E_k- p*k; 
                xp[j][2]=2*kdp/(aMD*aMD);//  /M_D^2
                
                xp[j][9]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3]);
                xp[j][9]=2*sin(xp[j][9]/2.);
                
                xp[j][10]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3]);
                xp[j][10]=2*sin(xp[j][10]/2.);
                
                 
                    cp2=0;//derivative_xi(9, j , fit_out,  fit_info, xp[j]);
                    ck2=0;//derivative_xi(10, j , fit_out,  fit_info, xp[j]);
                    if(strcmp(AV,"A")==0)  F[j]=gJ[e].FA[i][j] * gJ[e].ZV[j]; 
                    else if(strcmp(AV,"V")==0) F[j]=gJ[e].FV[i][j] * gJ[e].ZA[j];
                    else if(strcmp(AV,"{A^H}")==0){
                        F[j]=gJ[e].FA_from_H0[i][j] ;
                         if(iks==0   && ikt ==0)
                             F[j]=gJ[e].FA_from_H0[i][j] ;
                         else if (iks==1 && ikt ==0  || iks==2 && ikt ==0 )
                             F[j]=inter_2(head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+2]*gJ[e].w0[j]/gJ[e].Zp[j],
                                        gJ[e].FA_from_H0[iK1][j],                 gJ[e].FA_from_H0[iK2][j], 
                                         result.msw[j]);
                         else if (iks==0 && ikt ==3  || iks==0 && ikt ==4 )
                             F[j]=inter_2(head[e].k[head[e].nk+3]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+4]*gJ[e].w0[j]/gJ[e].Zp[j],
                                        gJ[e].FA_from_H0[iD1][j],                 gJ[e].FA_from_H0[iD2][j], 
                                        result.mcw[j]);    
                         else if (iks==1 && ikt ==3  || iks==1 && ikt ==4 || iks==2 && ikt ==3  || iks==2 && ikt ==4 )    
                             F[j]=inter_4(head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+2]*gJ[e].w0[j]/gJ[e].Zp[j],
                                        head[e].k[head[e].nk+3]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+4]*gJ[e].w0[j]/gJ[e].Zp[j],
                                        gJ[e].FA_from_H0[iDs1c1][j],                 gJ[e].FA_from_H0[iDs2c1][j], 
                                        gJ[e].FA_from_H0[iDs1c2][j],                 gJ[e].FA_from_H0[iDs2c2][j], 
                                            result.msw[j], result.mcw[j]);
                        
                    }
                    else if(strcmp(AV,"{V^H}")==0) F[j]=gJ[e].FV_from_H0[i][j] * gJ[e].ZA[j];
                    else if(strcmp(AV,"{V^HA}")==0){
                        if(iks==0   && ikt ==0)
                             F[j]=gJ[e].FV_from_H0_HA[i][j];
                         else if (iks==1 && ikt ==0  || iks==2 && ikt ==0 )
                             F[j]=inter_2(head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+2]*gJ[e].w0[j]/gJ[e].Zp[j],
                                        gJ[e].FV_from_H0_HA[iK1][j],                 gJ[e].FV_from_H0_HA[iK2][j], 
                                         result.msw[j]);
                         else if (iks==0 && ikt ==3  || iks==0 && ikt ==4 )
                             F[j]=inter_2(head[e].k[head[e].nk+3]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+4]*gJ[e].w0[j]/gJ[e].Zp[j],
                                        gJ[e].FV_from_H0_HA[iD1][j],                 gJ[e].FV_from_H0_HA[iD2][j], 
                                        result.mcw[j]);    
                         else if (iks==1 && ikt ==3  || iks==1 && ikt ==4 || iks==2 && ikt ==3  || iks==2 && ikt ==4 )    
                             F[j]=inter_4(head[e].k[head[e].nk+1]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+2]*gJ[e].w0[j]/gJ[e].Zp[j],
                                          head[e].k[head[e].nk+3]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+4]*gJ[e].w0[j]/gJ[e].Zp[j],
                                        gJ[e].FV_from_H0_HA[iDs1c1][j],                 gJ[e].FV_from_H0_HA[iDs2c1][j], 
                                        gJ[e].FV_from_H0_HA[iDs1c2][j],                 gJ[e].FV_from_H0_HA[iDs2c2][j], 
                                            result.msw[j], result.mcw[j]);
                         F[j]*=gJ[e].ZAV[j];
                        
                    }
                    else error(0==0,1,"argument of function  AV is nor A neither V","");
                    F[j]=F[j]-cp2*xp[j][9]*xp[j][9]-ck2*xp[j][10]*xp[j][10];
                }
                m=mean_and_error("jack",jack_tot,F);
                //fprintf(fc,"%g  %g   %g   %g  \t %g  %g  %g  %g\n",xp[Njack-1][2],xp[Njack-1][2] ,m[0],m[1],xp[Njack-1][3],xp[Njack-1][5],xp[Njack-1][7],xp[Njack-1][8]);
                fprintf(fc,"%g  %g   %g   %g     \t",xp[Njack-1][2],xp[Njack-1][2] ,m[0],m[1] );
                for (j=0;j<Njack;j++)
                    F[j]=xp[j][3];   //Mpi w0
                m=mean_and_error("jack",jack_tot,F);
                fprintf(fc,"%g  %g ",m[0],m[1]);
                for (j=0;j<Njack;j++)
                    F[j]=xp[j][5];   //MK w0
                m=mean_and_error("jack",jack_tot,F);
                fprintf(fc,"%g  %g ",m[0],m[1]);
                for (j=0;j<Njack;j++)
                    F[j]=xp[j][7];   //MD w0
                m=mean_and_error("jack",jack_tot,F);
                fprintf(fc,"%g  %g ",m[0],m[1]);
                for (j=0;j<Njack;j++)
                    F[j]=xp[j][8];   //MDs w0
                m=mean_and_error("jack",jack_tot,F);
                fprintf(fc,"%g  %g \t",m[0],m[1]);
                
                fprintf(fc,"%g  %g\n",xp[Njack-1][9],xp[Njack-1][10] );
                free(m);
                }
            }}}
            fprintf(fc,"\n\n");
        }}
        fclose(fc);

    }
    
    free_2(Njack,xp);free(F);
    
}


double  **taylor_expand_fit_fun(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, const char *AV ,int order, double xg_val){

    int i,j;
    char name[NAMESIZE];
    FILE *fc=NULL;
    int i_mD,i_mDs,i_mK,i_mp,e,ikt,iks;
    double *m,cp2,ck2;
    int Neff=0;
    int imom0,imoms,imomt;
    double **F=(double**) malloc(sizeof(double*)*(order+1));
    int Njack=jack_tot;
    double xg_min=0.001;
    double **x=(double**) malloc(sizeof(double*)*jack_tot);
    double h=1e-5;   
    double *d;
    double **tif=fit_to_tif(fit_info.Npar,jack_tot,fit_out.P);    
    for(j=0;j<jack_tot;j++)
            x[j]=(double*) malloc(sizeof(double)*fit_info.Nvar);
    
    
    for (i=0;i<=order;i++){
        F[i]=(double*) malloc(sizeof(double)*(jack_tot));
        for(j=0;j<jack_tot;j++){
            x[j][0]=phys_point[j][0];//mlw
            x[j][1]=phys_point[j][1];//r0
            x[j][2]=xg_val;//=1.;//xg
            x[j][3]=phys_point[j][3];//=result.MpiMeV[j]*result.w0MeV[j];275MeV
            x[j][4]=phys_point[j][4];//=result.msw[j];
            x[j][5]=phys_point[j][5];//=result.MKMeV[j]*result.w0MeV[j];
            x[j][6]=phys_point[j][6];//=result.mcw[j];
            x[j][7]=phys_point[j][7];//=result.MDMeV[j]*result.w0MeV[j];
            x[j][8]=phys_point[j][8];//=result.MDsMeV[j]*result.w0MeV[j];
            x[j][9]=0;
            x[j][10]=0;
            x[j][11]=phys_point[j][11];//=result.fpiMeV_exp[j]*result.w0MeV[j];
            x[j][12]=phys_point[j][12];//=result.fKMeV_exp[j]*result.w0MeV[j];
            d=derN_fun_Nf_var_h(1, fit_info.Nvar, x[j], fit_info.Npar, tif[j], fit_info.function,  h, i);
            F[i][j]=d[2];
            free(d);
        }
        m=mean_and_error("jack",jack_tot,F[i]);
        printf("d%d=%g   %g\t",i,m[0],m[1]);   
        free(m);
    
    }
    printf("\n");
    free_tif(jack_tot,tif);
    free_tif(jack_tot,x);
    return F;
}

void  print_fit_info(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, struct reph_jack *grephJ, struct header *head , const char *AV,const char *namefile){
    int i,j,k;
    
    double **fit=fit_out.P,*chi2m;
    
    double **tif,*tmp,**tmp2,*fk;
    double **Ci=(double**) malloc(sizeof(double*)*fit_info.Npar);
    char nametex[NAMESIZE],namegp[NAMESIZE];
    mysprintf(nametex,NAMESIZE,"%s/%s.tex",argv[2],namefile);
    mysprintf(namegp,NAMESIZE,"%s/%s.gp",argv[2],namefile);
    FILE *ftex=NULL, *fgp=NULL;
    ftex=open_file(nametex,"w+");
    fgp=open_file(namegp,"w+");
    const char *s;
    tmp=(double*)  malloc(sizeof(double) *jack_tot);
    tmp2=(double**) malloc(sizeof(double*)*fit_info.Npar);
    for (i=0;i<fit_info.Npar;i++)
        tmp2[i]=(double*)  malloc(sizeof(double) *jack_tot);

    tif=fit_to_tif(fit_info.Npar,jack_tot,fit);    
    
    chi2m=mean_and_error_jack_biased(jack_tot,fit_out.chi2);
    
    for(j=0;j<jack_tot;j++){
        fk=der_fun_Nf_h(fit_info.N,  fit_info.Nvar, phys_point[j], fit_info.Npar,tif[j],  fit_info.function,  0.00001);
        for (i=0;i<fit_info.Npar;i++)
            tmp2[i][j]=fk[i]*fit[i][j];
        free(fk);
    }
    for (i=0;i<fit_info.Npar;i++){
        Ci[i]=mean_and_error_jack_biased(jack_tot,fit[i]);
    }
    fprintf(ftex,"\\begin{align}\n");
    fprintf(ftex,"& \\chi^2/d.o.f.= %+.5f \\pm \t%.2g \\\\ \n",chi2m[0],chi2m[1]);
    for (i=0;i<fit_info.Npar;i++){
        fprintf(ftex,"& P_{%d}= %+.5f \\pm \t%.2g   \\\\ \n",i,Ci[i][0],Ci[i][1]);
        //s=smean_and_error("jack",jack_tot,fit[i]);
        //fprintf(ftex,"P_{%d}= %s  &\\quad  ",i,s);
        //free((void*)s);
        //s=smean_and_error("jack",jack_tot,tmp2[i]);
        //fprintf(ftex,"& phys\\quad %s\\\\ \n",s);
        //free((void*)s);
    }
    fprintf(ftex,"\\end{align}\n");
    
    for (i=0;i<fit_info.Npar;i++){
        fprintf(fgp,"P%d= %+.5f;\t",i,Ci[i][0]);
    }
    fprintf(fgp,"\n");
    for (i=0;i<fit_info.Npar;i++){       
        free(Ci[i]);
    }
/*    
    if (strcmp(AV,"H")!=0){
        for(j=0;j<jack_tot;j++){
            if (strcmp(AV,"A")==0  || strcmp(AV,"V")==0  )
                phys_point[j][2]=1.;
            else if (strcmp(AV,"{A^H}")==0  || strcmp(AV,"{V^H}")==0 || strcmp(AV,"{V^HA}")==0  )    
                phys_point[j][2]=0.;
            else 
                error(1==1,1,"print_fit_info","AV is not A neither V, AV=%s",AV);
                
        }
        for(j=0;j<jack_tot;j++)
            tmp[j]=fit_info.function(fit_info.N,fit_info.Npar,phys_point[j],fit_info.Npar,tif[j]);

        Ci[0]=mean_and_error_jack_biased(jack_tot,tmp);
        fprintf(ftex,"\\begin{equation}\n    F_%s(a=0, M_P=M_\\pi,x_\\gamma=%.0f)= %f \\pm %.2g\n\\end{equation}\n", AV,phys_point[jack_tot-1][2],Ci[0][0],Ci[0][1]);        
        free(Ci[0]);

            
        for(j=0;j<jack_tot;j++){
            tmp[j]=fit[0][j];
            if (fit_info.Npar>=5)
                tmp[j]+=fit[4][j]*phys_point[j][3]*phys_point[j][3]; //P_4*Mpiw0^2
            if (fit_info.Npar>=10){
                tmp[j]+=fit[6][j]*phys_point[j][5]*phys_point[j][5]; //P_6*MKw0^2
                tmp[j]+=fit[8][j]*phys_point[j][7]; //P_8*MDw0
            }
        }
        Ci[0]=mean_and_error_jack_biased(jack_tot,tmp);
        if (fit_info.Npar<6) fprintf(ftex,"\\begin{equation}\n    C_0= P_0 =%f \\pm %.2g\n\\end{equation}\n",Ci[0][0],Ci[0][1]);
        if (fit_info.Npar>=6 && fit_info.Npar<10 ) fprintf(ftex,"\\begin{equation}\n    C_0= P_0+P_4 M_\\pi^2r_0^2  =%f \\pm %.2g\n\\end{equation}\n",Ci[0][0],Ci[0][1]);
        if (fit_info.Npar>=10)fprintf(ftex,"\\begin{equation}\n    C_0= P_0+P_4 M_\\pi^2r_0^2  +P_6 M_K^2r_0^2+P_8 M_Dr_0 =%f \\pm %.2g\n\\end{equation}\n",Ci[0][0],Ci[0][1]);
        free(Ci[0]);
        for(j=0;j<jack_tot;j++){
            tmp[j]=fit[1][j];
            if (fit_info.Npar>=6)
                tmp[j]+=fit[5][j]*phys_point[j][3]*phys_point[j][3]; //P_4*Mpiw0^2
            if (fit_info.Npar>=10){
                tmp[j]+=fit[7][j]*phys_point[j][5]*phys_point[j][5]; //P_4*MKw0^2
                tmp[j]+=fit[9][j]*phys_point[j][7]; //P_4*MDw0
            }
        }
        Ci[0]=mean_and_error_jack_biased(jack_tot,tmp);
        if (fit_info.Npar<5) fprintf(ftex,"\\begin{equation}\n    C_1= P_1 =%f \\pm %.2g\n\\end{equation}\n",Ci[0][0],Ci[0][1]);
        if (fit_info.Npar>=6 && fit_info.Npar<10 )  fprintf(ftex,"\\begin{equation}\n    C_1= P_1+P_5 M_\\pi^2r_0^2  =%f \\pm %.2g\n\\end{equation}\n",Ci[0][0],Ci[0][1]);
        if (fit_info.Npar>=10)fprintf(ftex,"\\begin{equation}\n    C_1= P_1+P_5 M_\\pi^2r_0^2  +P_7 M_K^2r_0^2+P_9 M_Dr_0 =%f \\pm %.2g\n\\end{equation}\n",Ci[0][0],Ci[0][1]);
        free(Ci[0]);
    }
    if (strcmp(AV,"H")==0){
        for(j=0;j<jack_tot;j++)
            tmp[j]=fit_info.f1(fit_info.N,fit_info.Npar,phys_point[j],fit_info.Npar,tif[j]);
        Ci[0]=mean_and_error_jack_biased(jack_tot,tmp);
        fprintf(ftex,"\\begin{equation}\n    F_A(a=0, M_P=M_\\pi,x_\\gamma=1)= %g \\pm %.2g\n\\end{equation}\n",Ci[0][0],Ci[0][1]);        
        free(Ci[0]);
        for(j=0;j<jack_tot;j++)
            tmp[j]=fit_info.f2(fit_info.N,fit_info.Npar,phys_point[j],fit_info.Npar,tif[j])/result.w0MeV[j];
        Ci[0]=mean_and_error_jack_biased(jack_tot,tmp); 
        fprintf(ftex,"\\begin{equation}\n    f_p(a=0, M_P=M_\\pi)= (%g \\pm %.2g)\\, MeV\n\\end{equation}\n",Ci[0][0],Ci[0][1]);
        free(Ci[0]);
    }
     */
     
     print_chiral_extrapolation(argv, jack_tot,  fit_out,   fit_info,  phys_point, grephJ, head ,   AV,  namefile);
     print_continuum_extrapolation(argv, jack_tot,  fit_out,   fit_info,  phys_point,   AV,  namefile);
     print_chiral_continuum_fit(argv, jack_tot,  fit_out,   fit_info,  phys_point,   AV,  namefile,  head ,   grephJ);
     print_chiral_continuum_fit_phys(argv, jack_tot,  fit_out,   fit_info,  phys_point,   AV,  namefile,  head ,   grephJ);
     //interpolating_for_plot(argv, jack_tot,  fit_out,   fit_info,  phys_point,    AV,  namefile,  head ,   grephJ );
     print_data_without_p2k2(argv, jack_tot,  fit_out,   fit_info,  phys_point,   AV, namefile,  head ,   grephJ );
     print_data_without_p2k2_phys(argv, jack_tot,  fit_out,   fit_info,  phys_point,   AV, namefile,  head ,   grephJ );

    
    //for (i=0;i<fit_info.Npar;i++)        
    //    free(Ci[i]);
    //free(Ci);
    fprintf(ftex,"{\\tiny\\begin{gather}\n C=\\begin{pmatrix}\n");
    for(j=0;j<jack_tot;j++){
        for (i=0;i<fit_info.Npar;i++){
            for (k=0;k<i;k++)
                    jack_tot,fit_out.C[i][k][j]/=sqrt(fit_out.C[i][i][j]*fit_out.C[k][k][j]);
            for (k=i+1;k<fit_info.Npar;k++)
                    jack_tot,fit_out.C[i][k][j]/=sqrt(fit_out.C[i][i][j]*fit_out.C[k][k][j]);
        }
    }

     for (i=0;i<fit_info.Npar;i++){
        for (j=0;j<fit_info.Npar;j++){
            
            s=smean_and_error("jack",jack_tot,fit_out.C[i][j]);
            if (j==0)  fprintf(ftex,"%s",  s );
            else       fprintf(ftex,"& %s", s );
            free((void*)s);
        }
        if (i!=fit_info.Npar) fprintf(ftex,"\\\\ \n");
        else fprintf(ftex,"\n");
    }
    fprintf(ftex,"\\end{pmatrix}\n\\end{gather}}\n");
    
    int order=2;
    double **taylor=taylor_expand_fit_fun( argv, jack_tot,   fit_out,   fit_info,  phys_point,   AV, order,0);
    double **cov=covariance("jack",order+1,jack_tot,taylor);
    for(i=0;i<=order;i++){
        Ci[0]=mean_and_error_jack_biased(jack_tot,taylor[i]); 
        fprintf(ftex,"\\begin{equation}\n \\frac{1}{%d!} \\partial_{x_\\gamma}^{%d}  F_%s\\bigg|_{x_\\gamma=0}= %f \\pm %.2g\n\\end{equation}\n", i,i,AV,Ci[0][0],Ci[0][1]);        
        free(Ci[0]);
        
        free(taylor[i]);
    }
    free(taylor);
    fprintf(ftex,"{\\tiny\\begin{gather}\n C=\\begin{pmatrix}\n");
     for (i=0;i<order+1;i++){
        for (j=0;j<order+1;j++){
            if (j==0)  fprintf(ftex,"%g",  cov[i][j] );
            else       fprintf(ftex,"& %g", cov[i][j] );
        }
        if (i!=order+1) fprintf(ftex,"\\\\ \n");
        else fprintf(ftex,"\n");
    }
    fprintf(ftex,"\\end{pmatrix}\n\\end{gather}}\n");
    free_tif(order+1,cov);
    
    for (i=0;i<fit_info.Npar;i++){
        free(fit_out.P[i]);
        for (j=0;j<fit_info.Npar;j++){
            free(fit_out.C[i][j]);
        }
        free(fit_out.C[i]);
        free(tmp2[i]);
    }     
    free(tmp2);
    free(fit_out.chi2);
    free(fit_out.P);
    free(fit_out.C);
    
    free_tif(jack_tot,tif);
    //free_tif(fit_info.Npar,fit);
    fclose(ftex);fclose(fgp);free(chi2m);free(tmp); free(Ci);
    
}

double FV_slope_1_xg(int n, int Nvar, double *x,int Npar,double  *P){
    return P[0]*(1+P[1]*(1-x[0]));
}
double FV_slope_0_xg(int n, int Nvar, double *x,int Npar,double  *P){
    return P[0]+P[1]*(x[0]);
}
double FV_pole(int n, int Nvar, double *x,int Npar,double  *P){
    return P[0]/(1+P[1]*x[0]);
}

void  compute_FApmFV(char **argv, double **phys_point, struct fit_all fit_chi2_good, const char *meson, const char *option){
    int i,j,order=3,ord,k;
    int ord1=1;
    int Njack=fit_chi2_good.out[0].Njack;
    int Nvar=fit_chi2_good.info[0].Nvar;
    int Nfits=fit_chi2_good.Nfits/2;
    double **taylor,**taylor1;
    double **tay,CD_AV;
    double **tmp=(double**)  malloc((order+1)*sizeof(double*));
    double **tmpe=(double**) malloc((order+1)*sizeof(double*));
    double *C;
    double *ave=(double*)  calloc((order+1),sizeof(double));
    double *sigma=(double*)  calloc((order+1),sizeof(double));
    double ***tmpc=(double***) malloc(sizeof(double**)* Nfits);
    double **cov=(double**) malloc(sizeof(double*)* (order+1));
    
    
    error(strcmp(option,"pm")!=0  && strcmp(option,"_correlated_")!=0 && strcmp(option,"_correlated_pole_")!=0 , 1, "compute_FApmFV", "option = %s\n  while ammitted options are : pm  _correlated_   _correlated_pole_",option);
    
    
    for(ord=0;ord<order+1;ord++){
        tmp[ord]=(double*) malloc(fit_chi2_good.Nfits*sizeof(double));
        tmpe[ord]=(double*) malloc(fit_chi2_good.Nfits*sizeof(double));
        cov[ord]=(double*) calloc((order+1),sizeof(double));
    }
    char name[NAMESIZE];
    mysprintf(name,NAMESIZE,"%s/systematics_FA%sFV_%s.tex",argv[2],option,meson);
    FILE *f=NULL;
    f=open_file(name,"w+");
    printf("result including  systematics in %s, combining %d mesures\n",meson,Nfits);
    error(fit_chi2_good.Nfits%2!=0,1,"error in compute_FApmFV","fit_chi2_good.Nfits should be a multiple of 2 , also check the order");
    for (i=0;i<Nfits;i++){
        taylor=taylor_expand_fit_fun( argv, fit_chi2_good.out[i].Njack,   fit_chi2_good.out[i],   fit_chi2_good.info[i],  phys_point,   meson, ord1,0);//taylor of FA
        k=i+Nfits;//index to FV
        taylor1=taylor_expand_fit_fun( argv, fit_chi2_good.out[k].Njack,   fit_chi2_good.out[k],   fit_chi2_good.info[k],  phys_point,   meson, ord1,0);//taylor of FV
        tay=double_malloc_2(order+1,Njack);
        for (j=0;j<fit_chi2_good.out[i].Njack;j++){   
            if (strcmp(option,"pm")==0){
                tay[3][j]=(taylor[0][j]+taylor[1][j])-(taylor1[0][j]+taylor1[1][j]);//FA+FV in xg=1
                tay[2][j]=(taylor[0][j]+taylor[1][j])+(taylor1[0][j]+taylor1[1][j]);//FA-FV in xg=1
                tay[1][j]=taylor[0][j]-taylor1[0][j];//FA-FV in xg=0
                tay[0][j]=taylor[0][j]+taylor1[0][j];//FA+FV in xg=0
            }
            else if (strcmp(option,"_correlated_")==0){
                tay[3][j]=taylor1[1][j];//D_V
                tay[2][j]=taylor[1][j];//D_A
                tay[1][j]=taylor1[0][j];//C_V
                tay[0][j]=taylor[0][j];//C_A
            }
            else if (strcmp(option,"_correlated_pole_")==0){
                tay[3][j]=-taylor1[1][j]/taylor1[0][j];//D_V
                tay[2][j]=-taylor[1][j]/taylor[0][j];//D_A
                tay[1][j]=taylor1[0][j];//C_V
                tay[0][j]=taylor[0][j];//C_A
            }
            
        }
        tmpc[i]=covariance("jack",order+1,jack_tot,tay);

        for(ord=0;ord<order+1;ord++){
            C=mean_and_error_jack_biased(jack_tot,tay[ord]);
            tmp[ord][i]=C[0];
            tmpe[ord][i]=C[1]*C[1];
            free(C);
        }
        free_tif(ord1+1,taylor);
        free_tif(ord1+1,taylor1);
        free_tif(order+1,tay);

        
    }
    printf("average\n");
    for(ord=0;ord<order+1;ord++){
        for (i=0;i<Nfits;i++){
            ave[ord]+=tmp[ord][i];
            sigma[ord]+=tmpe[ord][i];
            for(k=0;k<order+1;k++)
                cov[ord][k]+=tmpc[i][ord][k];
        }
    }
    for(ord=0;ord<order+1;ord++)
        ave[ord]/=(Nfits);
    
    for(ord=0;ord<order+1;ord++){
        for (i=0;i<Nfits;i++){
            sigma[ord]+=(tmp[ord][i]-ave[ord])*(tmp[ord][i]-ave[ord]);
            
        }
        for(k=0;k<order+1;k++)
            for (i=0;i<Nfits;i++)
                cov[ord][k]+=(tmp[ord][i]-ave[ord])*(tmp[k][i]-ave[k]);
            
        sigma[ord]/=(Nfits);
        sigma[ord]=sqrt(sigma[ord]);
        
        printf("d%d=%g  +- %g\n",ord,ave[ord],sigma[ord]);
        if (strcmp(option,"pm")==0){
            if (ord==0)
                fprintf(f,"\\begin{equation}\n F_A+F_V|_{x_\\gamma=0}= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);
            else if (ord==1)
                fprintf(f,"\\begin{equation}\n F_A-F_V|_{x_\\gamma=0}= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);
            else if (ord==2)
                fprintf(f,"\\begin{equation}\n F_A+F_V|_{x_\\gamma=1}=(F_A+\\partial_{x} F_A)+(F_V+\\partial_{x} F_V)|_{x_\\gamma=0}= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);
            else if (ord==3)
                fprintf(f,"\\begin{equation}\n F_A-F_V|_{x_\\gamma=1}=(F_A+\\partial_{x} F_A)-(F_V+\\partial_{x} F_V)|_{x_\\gamma=0}= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);
        }
        else if (strcmp(option,"_correlated_")==0){
            if (ord==0)
                fprintf(f,"\\begin{equation}\n C_A= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);
            else if (ord==1)
                fprintf(f,"\\begin{equation}\n C_V= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);
            else if (ord==2)
                fprintf(f,"\\begin{equation}\n D_A= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);
            else if (ord==3)
                fprintf(f,"\\begin{equation}\n D_V= %f \\pm %.2g\n\\end{equation}\n",ave[ord],sigma[ord]);
            
        }
        else if (strcmp(option,"_correlated_pole_")==0){
            if (ord==0)
                fprintf(f,"\\begin{equation}\n \\tilde C_A= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);
            else if (ord==1)
                fprintf(f,"\\begin{equation}\n \\tilde C_V= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);
            else if (ord==2)
                fprintf(f,"\\begin{equation}\n \\tilde D_A= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);
            else if (ord==3)
                fprintf(f,"\\begin{equation}\n \\tilde D_V= %f \\pm %.2g\n\\end{equation}\n",ave[ord],sigma[ord]);
            
        }
        free(tmpe[ord]);
    }
    
    for(ord=0;ord<order+1;ord++){
        for(i=ord+1;i<order+1;i++)
            cov[ord][i]/=sqrt( cov[ord][ord]*cov[i][i]  );
        for(i=0;i<ord;i++)
            cov[ord][i]/=sqrt( cov[ord][ord]*cov[i][i]  );
    }
    for(ord=0;ord<order+1;ord++)
        cov[ord][ord]=1.;
    
    fprintf(f,"{\\tiny\\begin{gather}\n C=\\begin{pmatrix}\n");
    for(ord=0;ord<(order+1);ord++){
        free(tmp[ord]);
        for(i=0;i<(order+1);i++){
            //cov[ord][i]/=Nfits;
            if (i==0)  fprintf(f,"%.3f",  cov[ord][i] );
            else       fprintf(f,"& %.3f", cov[ord][i] );
            //fprintf(f,"%g\t",cov[ord][i]);
        }
        if (ord!=order+1) fprintf(f,"\\\\ \n");
        else fprintf(f,"\n");
    }
    fprintf(f,"\\end{pmatrix}\n\\end{gather}}\n");
    free_2(order+1,cov);
    fclose(f);
    free(tmp);free(tmpe);free(sigma);free(ave);
    free_3(Nfits,order+1,tmpc);
    if (strcmp(option,"_correlated_")==0 || strcmp(option,"_correlated_pole_")==0 ){
        for (k=0;k<Nfits*2;k++){
            for (i=0;i<fit_chi2_good.info[k].Npar;i++){
                free(fit_chi2_good.out[k].P[i]);
                for (j=0;j<fit_chi2_good.info[k].Npar;j++){
                    free(fit_chi2_good.out[k].C[i][j]);
                }
                free(fit_chi2_good.out[k].C[i]);
            }     
            free(fit_chi2_good.out[k].chi2);
            free(fit_chi2_good.out[k].P);
            free(fit_chi2_good.out[k].C);
        }
        free(fit_chi2_good.info);
        free(fit_chi2_good.out);
    }
    
}

void compute_systematics(char **argv, double **phys_point, struct fit_all fit_chi2_good, const char *AV){
    int i,j,order=4,ord,k;
    int ord1=1;
    int Njack=fit_chi2_good.out[0].Njack;
    int Nvar=fit_chi2_good.info[0].Nvar;
    int Nfits=fit_chi2_good.Nfits;
    double **taylor,**tay=double_malloc_2(5,Njack);
    double **tmp=(double**)  malloc((order+1)*sizeof(double*));
    double **tmpe=(double**) malloc((order+1)*sizeof(double*));
    double *C;
    double *ave=(double*)  calloc((order+1),sizeof(double));
    double *sigma=(double*)  calloc((order+1),sizeof(double));
    double ***tmpc=(double***) malloc(sizeof(double**)* Nfits);
    double **cov=(double**) malloc(sizeof(double*)* (order+1));
    for(ord=0;ord<order+1;ord++){
        tmp[ord]=(double*) malloc(fit_chi2_good.Nfits*sizeof(double));
        tmpe[ord]=(double*) malloc(fit_chi2_good.Nfits*sizeof(double));
        cov[ord]=(double*) calloc((order+1),sizeof(double));
    }
    char name[NAMESIZE];
    mysprintf(name,NAMESIZE,"%s/systematics_coefficients_%s.tex",argv[2],AV);
    FILE *f=NULL;
    f=open_file(name,"w+");
    printf("result including  systematics in %s, combining %d mesures\n",AV,Nfits);

    for (i=0;i<fit_chi2_good.Nfits;i++){
        tay=taylor_expand_fit_fun( argv, fit_chi2_good.out[i].Njack,   fit_chi2_good.out[i],   fit_chi2_good.info[i],  phys_point,   AV, ord1/*order*/,0);
        taylor=double_malloc_2(order+1,Njack);
        for (j=0;j<fit_chi2_good.out[i].Njack;j++){        //set taylor[2] as the coefficient of the pole
            taylor[0][j]=tay[0][j];
            taylor[1][j]=tay[1][j];
            taylor[2][j]=-tay[1][j]/tay[0][j];
            taylor[3][j]=-tay[1][j]/(tay[0][j]+tay[1][j]);
            taylor[4][j]=(tay[0][j]+tay[1][j]);
        }
        tmpc[i]=covariance("jack",order+1,jack_tot,taylor);
        for(ord=0;ord<order+1;ord++){
            C=mean_and_error_jack_biased(jack_tot,taylor[ord]);
            tmp[ord][i]=C[0];
            tmpe[ord][i]=C[1]*C[1];
            free(C);
        }
        free_tif(order+1,taylor);free_2(ord1+1,tay);
        
    }
    printf("average\n");
    for(ord=0;ord<order+1;ord++){
        for (i=0;i<fit_chi2_good.Nfits;i++){
            ave[ord]+=tmp[ord][i];
            sigma[ord]+=tmpe[ord][i];
            for(k=0;k<order+1;k++)
                cov[ord][k]+=tmpc[i][ord][k];
        }
    }
    for(ord=0;ord<order+1;ord++)
        ave[ord]/=fit_chi2_good.Nfits;
    
    for(ord=0;ord<order+1;ord++){
        for (i=0;i<fit_chi2_good.Nfits;i++){
            sigma[ord]+=(tmp[ord][i]-ave[ord])*(tmp[ord][i]-ave[ord]);
            
        }
        for(k=0;k<order+1;k++)
            for (i=0;i<fit_chi2_good.Nfits;i++)
                cov[ord][k]+=(tmp[ord][i]-ave[ord])*(tmp[k][i]-ave[k]);
            
        sigma[ord]/=fit_chi2_good.Nfits;
        sigma[ord]=sqrt(sigma[ord]);
        
        printf("d%d=%g  +- %g\n",ord,ave[ord],sigma[ord]);
        if (ord<2)
            fprintf(f,"\\begin{equation}\n \\frac{1}{%d!} \\partial_{x_\\gamma}^{%d}  F\\bigg|_{x_\\gamma=0}= %f \\pm %.2g\n\\end{equation}\n\t\t",ord,ord,ave[ord],sigma[ord]);
        else if (ord==2)
            fprintf(f,"\\begin{equation}\n  \\tilde D=-\\frac{\\partial_{x_\\gamma}  F}{F}\\bigg|_{x_\\gamma=0}= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);    
        else if (ord==3)
            fprintf(f,"\\begin{equation}\n  a=-\\frac{\\partial_{x_\\gamma}  F}{F+\\partial_{x_\\gamma}  F}\\bigg|_{x_\\gamma=0}= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);    
        else if (ord==4)
            fprintf(f,"\\begin{equation}\n  F(x_\\gamma=1)=F+\\partial_{x_\\gamma}  F\\bigg|_{x_\\gamma=0}= %f \\pm %.2g\n\\end{equation}\n\t\t",ave[ord],sigma[ord]);    
        
        
        free(tmpe[ord]);
    }
    
    for(ord=0;ord<order+1;ord++){
        for(i=ord+1;i<order+1;i++)
            cov[ord][i]/=sqrt( cov[ord][ord]*cov[i][i]  );
        for(i=0;i<ord;i++)
            cov[ord][i]/=sqrt( cov[ord][ord]*cov[i][i]  );
    }
    for(ord=0;ord<order+1;ord++)
        cov[ord][ord]=1.;
    
    fprintf(f,"{\\tiny\\begin{gather}\n C=\\begin{pmatrix}\n");
    for(ord=0;ord<(order+1);ord++){
        free(tmp[ord]);
        for(i=0;i<(order+1);i++){
            //cov[ord][i]/=Nfits;
            if (i==0)  fprintf(f,"%.3f",  cov[ord][i] );
            else       fprintf(f,"& %.3f", cov[ord][i] );
            //fprintf(f,"%g\t",cov[ord][i]);
        }
        if (ord!=order+1) fprintf(f,"\\\\ \n");
        else fprintf(f,"\n");
    }
    fprintf(f,"\\end{pmatrix}\n\\end{gather}}\n");
    free_2(order+1,cov);
    fclose(f);
    free(tmp);free(tmpe);free(sigma);free(ave);
    free_3(Nfits,order+1,tmpc);
     
    
    //Syntetic data//
    FILE *f1=NULL;
    mysprintf(name,NAMESIZE,"%s/systematics_data_%s.txt",argv[2],AV);
    f1=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/systematics_data_%s.tex",argv[2],AV);
    f=open_file(name,"w+");
    
    int Nxg=11;
    double **x=(double **) malloc(sizeof(double*)*fit_chi2_good.out[0].Njack);
    double **tif;
    double **tmp1=double_malloc_2(Nxg,Njack); 
    tmp=(double**)  malloc((Nxg)*sizeof(double*));
    tmpe=(double**) malloc((Nxg)*sizeof(double*));
    ave=(double*)  calloc((Nxg),sizeof(double));
    sigma=(double*)  calloc((Nxg),sizeof(double));
    tmpc=double_malloc_3(Nfits,Nxg,Nxg);
    cov=double_calloc_2(Nxg,Nxg);
    double **corr=double_calloc_2(Nxg,Nxg);
    
    for(ord=0;ord<Nxg;ord++){
        tmp[ord]=(double*) malloc(fit_chi2_good.Nfits*sizeof(double));
        tmpe[ord]=(double*) malloc(fit_chi2_good.Nfits*sizeof(double));
    }
    for (i=0;i<Nfits;i++){
        tif=fit_to_tif(fit_chi2_good.info[i].Npar, Njack, fit_chi2_good.out[i].P);
        for(k=0;k<Nxg;k++){
            for(j=0;j<Njack;j++){
                x[j]=(double*) malloc(sizeof(double)*Nvar);
                x[j][0]=phys_point[j][0];//mlw
                x[j][1]=phys_point[j][1];//r0/a 1e+10
                x[j][2]=0.1*((double) k);//=1.;//xg
                x[j][3]=phys_point[j][3];//=result.MpiMeV[j]*result.w0MeV[j];
                x[j][4]=phys_point[j][4];//=result.msw[j];
                x[j][5]=phys_point[j][5];//=result.MKMeV[j]*result.w0MeV[j];
                x[j][6]=phys_point[j][6];//=result.mcw[j];
                x[j][7]=phys_point[j][7];//=result.MDMeV[j]*result.w0MeV[j];
                x[j][8]=phys_point[j][8];//=result.MDsMeV[j]*result.w0MeV[j];
                x[j][9]=0;
                x[j][10]=0;
                x[j][11]=phys_point[j][11];//=result.fpiMeV_exp[j]*result.w0MeV[j];
                x[j][12]=phys_point[j][12];//=result.fKMeV_exp[j]*result.w0MeV[j];
                 
                tmp1[k][j]=fit_chi2_good.info[i].function(fit_chi2_good.info[i].N,Nvar,x[j],fit_chi2_good.info[i].Npar,tif[j]);

                free(x[j]);
            }
            //C=mean_and_error("jack",Njack,tmp1);
            C=mean_and_error_jack_biased(Njack,tmp1[k]);
            if(k==0) printf("%g\n",C[0]);
            tmp[k][i]=C[0];
            tmpe[k][i]=C[1]*C[1];
            free(C);
        }
        tmpc[i]=covariance("jack",Nxg,jack_tot,tmp1);

        free_tif(Njack,tif);
    }
    free(x);
    printf("average\n");
    for(ord=0;ord<Nxg;ord++){
        for (i=0;i<fit_chi2_good.Nfits;i++){
                ave[ord]+=tmp[ord][i];
                sigma[ord]+=tmpe[ord][i];
                for(k=0;k<Nxg;k++)
                    cov[ord][k]+=tmpc[i][ord][k];

        }
    }
    for(ord=0;ord<Nxg;ord++)
        ave[ord]/=fit_chi2_good.Nfits;
    
    fprintf(f,"\\begin{eqnarray} && \\qquad \\qquad \\qquad \\qquad   \\qquad \\qquad \\qquad \\qquad \\qquad   F_%c \\quad {\\rm Correlation} \\quad  {\\rm  Matrix} \\nonumber \\\\ &&\\footnotesize{\n",AV[1]);
    
    fprintf(f,"\\begin{tabular}{|c|c|c|}\n");
    fprintf(f,"\\hline\n");
    fprintf(f," $x_\\gamma$  & $F_%c $  & $\\Delta_{F_%c} $\\\\ \\hline\n",AV[1],AV[1]);
    for(ord=0;ord<Nxg;ord++){    
        for (i=0;i<fit_chi2_good.Nfits;i++){
                sigma[ord]+=(tmp[ord][i]-ave[ord])*(tmp[ord][i]-ave[ord]);
            
        }
        for(k=0;k<Nxg;k++)
            for (i=0;i<fit_chi2_good.Nfits;i++)
                    cov[ord][k]+=(tmp[ord][i]-ave[ord])*(tmp[k][i]-ave[k]);
                    
        sigma[ord]/=fit_chi2_good.Nfits;
        sigma[ord]=sqrt(sigma[ord]);
        
        
            
        printf("%.2g   %g   %g\n",0.1*ord,ave[ord],sigma[ord]);
        fprintf(f1,"%.2g   %g   %g\n",0.1*ord,ave[ord],sigma[ord]);
        fprintf(f,"%.2g &  %g &  %g \\\\ \\hline \n",0.1*ord,ave[ord],sigma[ord]);
        free(tmpe[ord]);
        
    }
    for(ord=0;ord<Nxg;ord++)
        for(k=0;k<Nxg;k++)
            cov[ord][k]/=fit_chi2_good.Nfits;
    fprintf(f,"\\end{tabular}\n");
    for(ord=0;ord<Nxg;ord++){
        for(i=ord+1;i<Nxg;i++)
            corr[ord][i]=cov[ord][i]/sqrt( cov[ord][ord]*cov[i][i]  );
        for(i=0;i<ord;i++)
            corr[ord][i]=cov[ord][i]/sqrt( cov[ord][ord]*cov[i][i]  );
    }
    for(ord=0;ord<Nxg;ord++)
        corr[ord][ord]=1.;
    
    fprintf(f,"\\quad \\left(\\begin{tabular}{");
    for(ord=0;ord<Nxg-1;ord++) 
        fprintf(f,"c|");
    fprintf(f,"c}");
    
    for(ord=0;ord<Nxg;ord++){
        free(tmp[ord]);
        for(i=0;i<Nxg;i++){
            //cov[ord][i]/=Nfits;
//            fprintf(f,"%g\t",cov[ord][i]);
            if (i==0)  fprintf(f,"%.3f",  corr[ord][i] );
            else       fprintf(f,"& %.3f", corr[ord][i] );

        }
        fprintf(f,"\\\\ \n");
    }
    fprintf(f,"\\end{tabular}\\right) \\nonumber\n  }\\end{eqnarray}\n");
    
    
    free(tmp);free(tmpe);
    fclose(f);  fclose(f1);
    free_2(order+1,corr);
    
    //fit of the synthetic data///////////////////////////////////////////////////////////////////////////////////////////////////
    mysprintf(name,NAMESIZE,"%s/fit_synthetic_data_%s.tex",argv[2],AV);
    f=open_file(name,"w+");
    
    double ***y=double_malloc_3(Njack,Nxg,2);      
    double **yj;
    double **xg=double_malloc_2(Nxg,1);
    double **fit=double_malloc_2(2,Njack);
    double *fit1;
    int en[1];
    en[0]=Nxg;
    
    int yn=is_it_positive( cov,  Nxg);
    while(yn==1){
        printf("covariance matrix not positive defined adding 0.0001*cov[0][0]*I \n");
        for(i=0;i<Nxg;i++)
             cov[i][i]+=cov[0][0]*1e-12;
        /*for(i=0;i<Nxg;i++)
             for(j=0;j<Nxg;j++)
                 cov[i][j]=0;
        for(i=0;i<Nxg;i++) 
            cov[i][i]=sigma[i]*sigma[i];
        */
        
        yn=is_it_positive( cov,  Nxg);  
        printf("now the matrix is positive defined.  %d\n",yn);
    }
     for(i=0;i<Nxg;i++){
             for(j=0;j<Nxg;j++)
                 printf("%g\t",cov[i][j]);
             printf("\n");
     }
    printf("\n");
    yj=fake_sampling_covariance("jack",ave,  Njack, Nxg, cov,33);
    double **covtmp=covariance("jack",Nxg,jack_tot,yj);
    for(i=0;i<Nxg;i++){
             for(j=0;j<Nxg;j++)
                 printf("%g\t",covtmp[i][j]);
             printf("\n");
     }
    for(i=0;i<Nxg;i++){
        printf("xg=%g\n",i*0.1);
    C=mean_and_error_jack_biased(Njack,yj[i]);
    printf("test cov=%f  +- %f   \n",C[0],C[1]);
    free(C);
    
    }
    /*
    if (yn==0){
        yj=fake_sampling_covariance("jack",ave,  Njack, Nxg, cov,33);
    }
    else if (yn==1){
        printf("covariance matrix not positive defined adding 0.0001*I \n");
        yj=double_calloc_2(Njack,Nxg);
        for(i=0;i<Nxg;i++){
          double *yjt=fake_sampling("jack",ave[i],sigma[i],Njack,33);  
          for (j=0;j<Njack;j++)
                yj[j][i]=yjt[j];
          free(yjt);
        }
    }*/
        
    for(i=0;i<Nxg;i++){
        xg[i][0]=0.1*((double) i);
        for (j=0;j<Njack;j++){
        y[j][i][0]=yj[i][j];
        y[j][i][1]=sqrt(cov[i][i]);
        }
    }
        
   
    double guess[2];
    guess[0]=1;guess[1]=0.1;
    double *chi2j=(double*) malloc(sizeof(double)*Njack);
    
    for (j=0;j<Njack;j++){
        fit1=non_linear_fit_Nf(1, en, xg , y[j] , 1,  2, FV_slope_1_xg ,guess);   
        chi2j[j]=compute_chi_non_linear_Nf(1, en,xg, y[j],fit1 ,1/* Nvar*/, 2/* Npar*/, FV_slope_1_xg  )/(Nxg-2);
        fit[0][j]=fit1[0];
        fit[1][j]=fit1[1];
        free(fit1);            
    }
    fprintf(f,"\\begin{gather} \n  F_%c(x_\\gamma)=F_%c(x_\\gamma=1)(1+\\lambda_%c(1-x_\\gamma))\n \\end{gather}  ",AV[1],AV[1],AV[1]);
    C=mean_and_error_jack_biased(Njack,chi2j);
    fprintf(f,"\\begin{gather}\n   \\chi^2/d.o.f.= %g \\pm %.2g\n \\\\ \n",C[0],C[1]);
    free(C);
    C=mean_and_error_jack_biased(Njack,fit[0]);
    fprintf(f,"  F_%c(x_\\gamma=1)= %f \\pm %.2g\n \\\\ \n",AV[1],C[0],C[1]);
    free(C);
    
    C=mean_and_error_jack_biased(Njack,fit[1]);
    fprintf(f,"  \\lambda_%c= %f \\pm %.2g\n\\end{gather}\n",AV[1],C[0],C[1]);
    free(C);

    for (j=0;j<Njack;j++){
        fit1=non_linear_fit_Nf(1, en, xg , y[j] , 1,  2, FV_slope_0_xg ,guess);
        chi2j[j]=compute_chi_non_linear_Nf(1, en,xg, y[j],fit1 ,1/* Nvar*/, 2/* Npar*/, FV_slope_0_xg  )/(Nxg-2);
        fit[0][j]=fit1[0];
        fit[1][j]=fit1[1];
        free(fit1);            
    }
    fprintf(f,"\\begin{gather} \n  F_%c(x_\\gamma)=C_%c+D_%c x_\\gamma\n \\end{gather}  ",AV[1],AV[1],AV[1]);
    C=mean_and_error_jack_biased(Njack,chi2j);
    fprintf(f,"\\begin{gather}\n   \\chi^2/d.o.f.= %g \\pm %.2g\n \\\\ \n",C[0],C[1]);
    free(C);
    C=mean_and_error_jack_biased(Njack,fit[0]);
    fprintf(f,"   C_%c= %f \\pm %.2g\n \\\\ \n",AV[1],C[0],C[1]);
    free(C);
    
    C=mean_and_error_jack_biased(Njack,fit[1]);
    fprintf(f,"  D_%c= %f \\pm %.2g\n\\end{gather}\n\t\t",AV[1],C[0],C[1]);
    free(C);

    for (j=0;j<Njack;j++){
        fit1=non_linear_fit_Nf(1, en, xg , y[j] , 1,  2, FV_pole ,guess);   
        chi2j[j]=compute_chi_non_linear_Nf(1, en,xg, y[j],fit1 ,1/* Nvar*/, 2/* Npar*/, FV_pole  )/(Nxg-2);
        fit[0][j]=fit1[0];
        fit[1][j]=fit1[1];
        free(fit1);            
    }
    fprintf(f,"\\begin{gather} \n  F_%c(x_\\gamma)=\\frac{\\tilde{C}_%c}{1+ \\tilde{D}_%c x_\\gamma}\n \\end{gather}  ",AV[1],AV[1],AV[1]);
    C=mean_and_error_jack_biased(Njack,chi2j);
    fprintf(f,"\\begin{gather}\n   \\chi^2/d.o.f.= %g \\pm %.2g\n \\\\ \n",C[0],C[1]);
    free(C);
    C=mean_and_error_jack_biased(Njack,fit[0]);
    fprintf(f,"   \\tilde{C}_%c= %f \\pm %.2g\n \\\\ \n",AV[1],C[0],C[1]);
    free(C);
    
    C=mean_and_error_jack_biased(Njack,fit[1]);
    fprintf(f,"  \\tilde{D}_%c= %f \\pm %.2g\n\\end{gather}\n\t\t",AV[1],C[0],C[1]);
    free(C);

    
    //Free//
    free(chi2j);
    free_2(Nxg,tmp1);free_2(2,fit);free_2(Nxg,yj);  free_2(Nxg,xg);free_3(Njack,Nxg,y);
    fclose(f);free(ave);free(sigma);   free_2(2,cov);

      

   
    for (k=0;k<Nfits;k++){
        for (i=0;i<fit_chi2_good.info[k].Npar;i++){
            free(fit_chi2_good.out[k].P[i]);
            for (j=0;j<fit_chi2_good.info[k].Npar;j++){
                free(fit_chi2_good.out[k].C[i][j]);
            }
            free(fit_chi2_good.out[k].C[i]);
        }     
        free(fit_chi2_good.out[k].chi2);
        free(fit_chi2_good.out[k].P);
        free(fit_chi2_good.out[k].C);
    }
   free(fit_chi2_good.info);
   free(fit_chi2_good.out);
}

/*void compute_systematics_jack(char **argv, double **phys_point, struct fit_all fit_chi2_good, const char *AV){
    int i,j,order=2,ord,k;
    int Njack=fit_chi2_good.out[0].Njack;
    int Nvar=fit_chi2_good.info[0].Nvar;
    int Nfits=fit_chi2_good.Nfits;
    double **taylor;
    double **tmp=(double**)  malloc((order+1)*sizeof(double*));
    double **tmpe=(double**) malloc((order+1)*sizeof(double*));
    double *C;
    double *ave=(double*)  calloc((order+1),sizeof(double));
    double *sigma=(double*)  calloc((order+1),sizeof(double));
    double **cov;
    for(ord=0;ord<order+1;ord++){
        tmp[ord]=(double*) calloc(Njack,sizeof(double));
        tmpe[ord]=(double*) malloc(fit_chi2_good.Nfits*sizeof(double));
    }
    char name[NAMESIZE];
    mysprintf(name,NAMESIZE,"%s/systematics_coefficients_%s.txt",argv[2],AV);
    FILE *f=NULL;
    f=open_file(name,"w+");
    printf("result including  systematics in %s, combining %d mesures\n",AV,Nfits);

    for (i=0;i<fit_chi2_good.Nfits;i++){
        taylor=taylor_expand_fit_fun( argv, fit_chi2_good.out[i].Njack,   fit_chi2_good.out[i],   fit_chi2_good.info[i],  phys_point,   AV, order);
        //cov[i]=covariance("jack",order+1,jack_tot,taylor);
        for(ord=0;ord<order+1;ord++){
            for(j=0;j<Njack;j++)
                
                tmp[ord][j]+=taylor[ord][j];
        }
        free_tif(order+1,taylor);
        
    }
    printf("average\n");
    for(ord=0;ord<order+1;ord++){
       for(j=0;j<Njack;j++)
           tmp[ord][j]/=Nfits;
       C=mean_and_error_jack_biased(jack_tot,tmp[ord]); 
       printf("d%d=%g  +- %g\n",ord,C[0],C[1]);
       fprintf(f,"%.15g\t%.15g\t\t",C[0],C[1]);
       free(C);
    }
    fprintf(f,"\n\n");
    cov=covariance("jack",order+1,Njack,tmp);
    for(ord=0;ord<order+1;ord++){
        for(i=0;i<order+1;i++){
            fprintf(f,"%g\t",cov[ord][i]);
        }
        fprintf(f,"\n");
        free(tmp[ord]);
    }
    fclose(f);free(tmp);
    

    mysprintf(name,NAMESIZE,"%s/systematics_data_%s.txt",argv[2],AV);
    f=open_file(name,"w+");
    int Nxg=11;
    double **x=(double **) malloc(sizeof(double*)*fit_chi2_good.out[0].Njack);
    double **tif;
    double *tmp1=(double*)  malloc((Njack)*sizeof(double));
    tmp=(double**)  malloc((Nxg)*sizeof(double*));
    tmpe=(double**) malloc((Nxg)*sizeof(double*));
    ave=(double*)  calloc((Nxg),sizeof(double));;
    sigma=(double*)  calloc((Nxg),sizeof(double));;
    for(ord=0;ord<Nxg;ord++){
        tmp[ord]=(double*) calloc(Njack,sizeof(double));
        tmpe[ord]=(double*) malloc(fit_chi2_good.Nfits*sizeof(double));
    }
    for (i=0;i<Nfits;i++){
        tif=fit_to_tif(fit_chi2_good.info[i].Npar, Njack, fit_chi2_good.out[i].P);
        for(k=0;k<Nxg;k++){
            for(j=0;j<Njack;j++){
                x[j]=(double*) malloc(sizeof(double)*Nvar);
                x[j][0]=phys_point[j][0];//mlw
                x[j][1]=phys_point[j][1];//r0/a 1e+10
                x[j][2]=0.1*((double) k);//=1.;//xg
                x[j][3]=phys_point[j][3];//=result.MpiMeV[j]*result.w0MeV[j];
                x[j][4]=phys_point[j][4];//=result.msw[j];
                x[j][5]=phys_point[j][5];//=result.MKMeV[j]*result.w0MeV[j];
                x[j][6]=phys_point[j][6];//=result.mcw[j];
                x[j][7]=phys_point[j][7];//=result.MDMeV[j]*result.w0MeV[j];
                x[j][8]=phys_point[j][8];//=result.MDsMeV[j]*result.w0MeV[j];
                x[j][9]=0;
                x[j][10]=0;
                
                tmp[k][j]+=fit_chi2_good.info[i].function(fit_chi2_good.info[i].N,Nvar,x[j],fit_chi2_good.info[i].Npar,tif[j]);

                free(x[j]);
            }
           
        }
        free_tif(Njack,tif);
    }
    
    printf("average\n");
    for(ord=0;ord<Nxg;ord++){
       for(j=0;j<Njack;j++)
           tmp[ord][j]/=Nfits;
       
       C=mean_and_error_jack_biased(jack_tot,tmp[ord]); 
       printf("%d    %g   %g\n",ord,C[0],C[1]);
       fprintf(f,"%d\t%.15g\t%.15g\n",(0.1)*ord,C[0],C[1]);
       free(C);
    }
    fprintf(f,"\n\n");
    cov=covariance("jack",Nxg,Njack,tmp);
    for(ord=0;ord<Nxg;ord++){
        for(i=0;i<Nxg;i++){
            fprintf(f,"%g\t",cov[ord][i]);
        }
        fprintf(f,"\n");
        free(tmp[ord]);
    }
    fclose(f);free(tmp);
    
    
    
    for (k=0;k<Nfits;k++){
        for (i=0;i<fit_chi2_good.info[k].Npar;i++){
            free(fit_chi2_good.out[k].P[i]);
            for (j=0;j<fit_chi2_good.info[k].Npar;j++){
                free(fit_chi2_good.out[k].C[i][j]);
            }
            free(fit_chi2_good.out[k].C[i]);
        }     
        free(fit_chi2_good.out[k].chi2);
        free(fit_chi2_good.out[k].P);
        free(fit_chi2_good.out[k].C);
    }
   free(fit_chi2_good.info);
   free(fit_chi2_good.out);
}
*/




static void  read_file_head_jack(FILE *stream,struct header *head)
{
    int i,nk_old,nmoms_old;
    if(head->allocated==1){
        nk_old=head->nk;
        nmoms_old=head->nmoms;
    }
    
    fread(&(head->twist),sizeof(int),1,stream);
    fread(&(head->nf),sizeof(int),1,stream);
    fread(&(head->nsrc),sizeof(int),1,stream);
    fread(&(head->l0),sizeof(int),1,stream);
    fread(&(head->l1),sizeof(int),1,stream);
    fread(&(head->l2),sizeof(int),1,stream);
    fread(&(head->l3),sizeof(int),1,stream);
    fread(&(head->nk),sizeof(int),1,stream);
    
    fread(&(head->nmoms),sizeof(int),1,stream);
    if(head->allocated==1)
        error(nk_old!=head->nk || nmoms_old!=head->nmoms,1,"read_file_head_jack", " file head nk or nmoms has changed");
    
    fread(&(head->beta),sizeof(double),1,stream);
    fread(&(head->ksea),sizeof(double),1,stream);
    fread(&(head->musea),sizeof(double),1,stream);
    fread(&(head->csw),sizeof(double),1,stream);
   
    if(head->allocated==0){
        head->k=(double*) malloc(sizeof(double)*2*head->nk);
        (*head).mom=(double**) malloc(head->nmoms*sizeof(double*));
        for(i=0;i<(head->nmoms);i++) {
                (*head).mom[i]=(double*) malloc(sizeof(double)*4);

        }
        head->allocated=1;
    }
    for(i=0;i<2*head->nk;++i)
    	fread(&(head->k[i]),sizeof(double),1,stream);
    
    for(i=0;i<head->nmoms;i++) {
        fread(&(head->mom[i][0]),sizeof(double),1,stream);
        fread(&(head->mom[i][1]),sizeof(double),1,stream);
        fread(&(head->mom[i][2]),sizeof(double),1,stream);
        fread(&(head->mom[i][3]),sizeof(double),1,stream);

    }
    
    /*
    head->k=(double*) malloc(sizeof(double)*2*head->nk);
    for(i=0;i<2*head->nk;++i)
    	fread(&(head->k[i]),sizeof(double),1,stream);
    
    head->mom=(double**) malloc(sizeof(double*)*head->nmoms);
    for(i=0;i<head->nmoms;i++) {
    	head->mom[i]=(double*) malloc(sizeof(double)*4);
        fread(&(head->mom[i][0]),sizeof(double),1,stream);
        fread(&(head->mom[i][1]),sizeof(double),1,stream);
        fread(&(head->mom[i][2]),sizeof(double),1,stream);
        fread(&(head->mom[i][3]),sizeof(double),1,stream);

    }
    */
    
    
}

int setup_reading_single_jack( struct header *head, FILE **f, const char *name){
    int N;
    
    *f=open_file(name,"r");
    read_file_head_jack(*f,head);
    fread(&N,sizeof(int),1,*f);
    
    return N;
}

void setup_reading_list_jack( struct  database_file_reph_jack *jack_files, struct header *head){
     int N,N1;
     
      N=setup_reading_single_jack(head,&(jack_files->f_FA_from_H0)  ,jack_files->FA_from_H0 );
      
      N1=setup_reading_single_jack(head,&(jack_files->f_FV_from_H0)       ,jack_files->FV_from_H0 );
      error(N!=N1,1,"setup_reading_list_jack","jacknifes have not the same number");
      N1=setup_reading_single_jack(head,&(jack_files->f_FV_from_H0_HA)       ,jack_files->FV_from_H0_HA );
      error(N!=N1,1,"setup_reading_list_jack","jacknifes have not the same number");
      N1=setup_reading_single_jack(head,&(jack_files->f_Zf_PS)       ,jack_files->Zf_PS );
      error(N!=N1,1,"setup_reading_list_jack","jacknifes have not the same number");
      N1=setup_reading_single_jack(head,&(jack_files->f_FAp)         ,jack_files->FAp );
      error(N!=N1,1,"setup_reading_list_jack","jacknifes have not the same number");
      N1=setup_reading_single_jack(head,&(jack_files->f_FA)          ,jack_files->FA );
      error(N!=N1,1,"setup_reading_list_jack","jacknifes have not the same number");
      N1=setup_reading_single_jack(head,&(jack_files->f_FV)          ,jack_files->FV );
      error(N!=N1,1,"setup_reading_list_jack","jacknifes have not the same number");
      N1=setup_reading_single_jack(head,&(jack_files->f_xG)          ,jack_files->xG );
      error(N!=N1,1,"setup_reading_list_jack","jacknifes have not the same number");
      N1=setup_reading_single_jack(head,&(jack_files->f_M_PS)        ,jack_files->M_PS );
      error(N!=N1,1,"setup_reading_list_jack","jacknifes have not the same number");
      
      jack_files->Njack=N;
     /*jack_files->f_FA_from_H0=fopen(jack_files->FA_from_H0,"r");
     error(jack_files->f_FA_from_H0==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->FA_from_H0);
     read_file_head_jack(jack_files->f_FA_from_H0,head);
     fread(&(jack_files->Njack),sizeof(int),1,jack_files->f_FA_from_H0);
     
     
     jack_files->f_Zf_PS=fopen(jack_files->Zf_PS,"r");
     error(jack_files->f_Zf_PS==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->Zf_PS);
     read_file_head_jack(jack_files->f_Zf_PS,head);
     fread(&N,sizeof(int),1,jack_files->f_Zf_PS);

     jack_files->f_FAp=fopen(jack_files->FAp,"r");
     error(jack_files->f_FAp==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->FAp);
     read_file_head_jack(jack_files->f_FAp,head);
     fread(&N,sizeof(int),1,jack_files->f_FAp);
     
     
     jack_files->f_FA=fopen(jack_files->FA,"r");
     error(jack_files->f_FA==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->FA);
     read_file_head_jack(jack_files->f_FA,head);
     fread(&N,sizeof(int),1,jack_files->f_FA);
     
     
     jack_files->f_FV=fopen(jack_files->FV,"r");
     error(jack_files->f_FV==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->FV);
     read_file_head_jack(jack_files->f_FV,head);
     fread(&N,sizeof(int),1,jack_files->f_FV);
     
     jack_files->f_xG=fopen(jack_files->xG,"r");
     error(jack_files->f_xG==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->xG);
     read_file_head_jack(jack_files->f_xG,head);
     fread(&N,sizeof(int),1,jack_files->f_xG);
     
     
     jack_files->f_M_PS=fopen(jack_files->M_PS,"r");
     error(jack_files->f_M_PS==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->M_PS);
     read_file_head_jack(jack_files->f_M_PS,head);
     fread(&N,sizeof(int),1,jack_files->f_M_PS);
     
     error(jack_files->Njack!=N,1,"setup_reading_single_jack", " files \n %s has %d elements \n %s has %d elements\n ",
         jack_files->FAp, jack_files->Njack,jack_files->FV,N); */
     
}

void  setup_reading_jack(char **option,struct  database_file_reph_jack *jack_files, struct header *head,const char  *name)  {
     
    if (strcmp(option[1],"manual")==0){
        mysprintf(jack_files->FA_from_H0,NAMESIZE,"%s/FA_from_H0_jack",name);
        mysprintf(jack_files->FV_from_H0,NAMESIZE,"%s/FV_from_H0_jack",name);
        mysprintf(jack_files->FV_from_H0_HA,NAMESIZE,"%s/FV_from_H0_HA_jack",name);
        mysprintf(jack_files->FAp,NAMESIZE,"%s/FAp_jack",name);
        mysprintf(jack_files->FA,NAMESIZE,"%s/FA_jack",name);
        mysprintf(jack_files->FV,NAMESIZE,"%s/FV_jack",name);  
        mysprintf(jack_files->xG,NAMESIZE,"%s/xG_jack",name);  
        mysprintf(jack_files->M_PS,NAMESIZE,"%s/M_{PS}_jack",name);  
        mysprintf(jack_files->Zf_PS,NAMESIZE,"%s/Zf_{PS}_jack",name);  
    }
    else if (strcmp(option[1],"auto")==0){
        mysprintf(jack_files->FA_from_H0,NAMESIZE,"%s/FA_from_H0_autoplateaux_jack",name);
        mysprintf(jack_files->FV_from_H0,NAMESIZE,"%s/FV_from_H0_autoplateaux_jack",name);
        mysprintf(jack_files->FV_from_H0_HA,NAMESIZE,"%s/FV_from_H0_HA_jack",name);
        mysprintf(jack_files->FAp,NAMESIZE,"%s/FAp_jack",name);
        mysprintf(jack_files->FA,NAMESIZE,"%s/FA_autoplateaux_jack",name);
        mysprintf(jack_files->FV,NAMESIZE,"%s/FV_autoplateaux_jack",name);  
        mysprintf(jack_files->xG,NAMESIZE,"%s/xG_jack",name);  
        mysprintf(jack_files->M_PS,NAMESIZE,"%s/M_{PS}_jack",name);  
        mysprintf(jack_files->Zf_PS,NAMESIZE,"%s/Zf_{PS}_jack",name);  
    }
    else 
        error(0==0,1,"setup_reading_jack"," argv[1] is not manual neither auto");
    
    setup_reading_list_jack(jack_files, head);
    
}



void  files_declarations(char **option, struct database_file_reph_jack **jack_files,struct header **head){
  
    (*head)=(struct header*) malloc(sizeof(struct header)*ensembles_reph);
    (*jack_files)=(struct database_file_reph_jack *) malloc (sizeof(struct database_file_reph_jack )*ensembles_reph);
    int i;
    for (i=0;i<ensembles_reph;i++)
        (*head)[i].allocated=0;
    char name[NAMESIZE];
    
    mysprintf(name,NAMESIZE,"%s/D15.48/jackknife",option[4]);
    setup_reading_jack( option,&((*jack_files)[0]),&((*head)[0]),name);  
    (*jack_files)[0].a=0.0619;
    if (ensembles_reph>1){
        mysprintf(name,NAMESIZE,"%s/D20.48/jackknife",option[4]);
        setup_reading_jack(option, &((*jack_files)[1]),&((*head)[1]),name);  
        (*jack_files)[1].a=0.0619;
    }
    if (ensembles_reph>2){
        mysprintf(name,NAMESIZE,"%s/D30.48/jackknife",option[4]);        
        setup_reading_jack(option, &((*jack_files)[2]),&((*head)[2]),name);  
        (*jack_files)[2].a=0.0619;
    }
    if (ensembles_reph>3){
        mysprintf(name,NAMESIZE,"%s/B25.32/jackknife",option[4]);
        setup_reading_jack(option, &((*jack_files)[3]),&((*head)[3]),name);  
        (*jack_files)[3].a=0.0815;
    }
    if (ensembles_reph>4){
        mysprintf(name,NAMESIZE,"%s/B35.32/jackknife",option[4]);
        setup_reading_jack(option, &((*jack_files)[4]),&((*head)[4]),name);  
        (*jack_files)[4].a=0.0815;
    }
    if (ensembles_reph>5){
        mysprintf(name,NAMESIZE,"%s/B55.32/jackknife",option[4]);
        setup_reading_jack(option, &((*jack_files)[5]),&((*head)[5]),name);  
        (*jack_files)[5].a=0.0815;
    }
    if (ensembles_reph>6){
        mysprintf(name,NAMESIZE,"%s/B75.32/jackknife",option[4]);
        setup_reading_jack(option, &((*jack_files)[6]),&((*head)[6]),name);  
        (*jack_files)[6].a=0.0815;
    }
    if (ensembles_reph>7){
        mysprintf(name,NAMESIZE,"%s/A30.32/jackknife",option[4]);
        setup_reading_jack(option, &((*jack_files)[7]),&((*head)[7]),name);  
        (*jack_files)[7].a=0.0885;
    }
    if (ensembles_reph>8){
        mysprintf(name,NAMESIZE,"%s/A40.32/jackknife",option[4]);
        setup_reading_jack(option, &((*jack_files)[8]),&((*head)[8]),name);  
        (*jack_files)[8].a=0.0885;
    }
    if (ensembles_reph>9){
        mysprintf(name,NAMESIZE,"%s/A60.24/jackknife",option[4]);
        setup_reading_jack(option, &((*jack_files)[9]),&((*head)[9]),name);  
        (*jack_files)[9].a=0.0885;
    }
    if (ensembles_reph>10){
        mysprintf(name,NAMESIZE,"%s/A80.24/jackknife",option[4]);
        setup_reading_jack(option, &((*jack_files)[10]),&((*head)[10]),name);  
        (*jack_files)[10].a=0.0885;
    }
   
}


double *create_fake_jack_Zp(int ensemble, struct database_file_reph_jack *jack_files)
{
    double *tmp;
    
    
     if (fabs(jack_files[ensemble].a-0.0619)<0.000001){
         //tmp=fake_jack(0.516,0.002,jack_files[ensemble].Njack);//M1 2Gev 
         tmp=fake_jack(0.545,0.002,jack_files[ensemble].Njack,123);//M2 2Gev 
     }
     if (fabs(jack_files[ensemble].a-0.0815)<0.000001){
        //tmp=fake_jack(0.509,0.004,jack_files[ensemble].Njack);//M1 2Gev  
        tmp=fake_jack(0.546,0.002,jack_files[ensemble].Njack,1234);//M2 2Gev  
    }
     if (fabs(jack_files[ensemble].a-0.0885)<0.000001){
        //tmp=fake_jack(0.529,0.007,jack_files[ensemble].Njack);//M1 2Gev  
        tmp=fake_jack(0.574,0.004,jack_files[ensemble].Njack,12345);//M2 2Gev  
     }
     
     return tmp;
}
double *create_fake_jack_ZA(int ensemble, struct database_file_reph_jack *jack_files)
{
    double *tmp;
    
    
     if (fabs(jack_files[ensemble].a-0.0619)<0.000001){
         //tmp=fake_jack(0.762,0.004,jack_files[ensemble].Njack);//M1 2Gev 
         tmp=fake_jack(0.752,0.002,jack_files[ensemble].Njack,123);//M2 2Gev 
     }
     if (fabs(jack_files[ensemble].a-0.0815)<0.000001){
        //tmp=fake_jack(0.737,0.005,jack_files[ensemble].Njack);//M1 2Gev  
        tmp=fake_jack(0.714,0.002,jack_files[ensemble].Njack,1234);//M2 2Gev  
    }
     if (fabs(jack_files[ensemble].a-0.0885)<0.000001){
        //tmp=fake_jack(0.531,0.008,jack_files[ensemble].Njack);//M1 2Gev  
        tmp=fake_jack(0.703,0.002,jack_files[ensemble].Njack,12345);//M2 2Gev  
     }
     
     return tmp;
}

double *create_fake_jack_ZV(int ensemble, struct database_file_reph_jack *jack_files)
{
    double *tmp;
    
    
     if (fabs(jack_files[ensemble].a-0.0619)<0.000001){
         tmp=fake_jack(0.6531,0.0002,jack_files[ensemble].Njack,123);//WTI
     }
     if (fabs(jack_files[ensemble].a-0.0815)<0.000001){
        tmp=fake_jack(0.6095,0.0003,jack_files[ensemble].Njack,1234);//WTI
    }
     if (fabs(jack_files[ensemble].a-0.0885)<0.000001){
        tmp=fake_jack(0.5920,0.0004,jack_files[ensemble].Njack,12345);//WTI
     }
     
     return tmp;
}

double *create_fake_jack_w0(int ensemble, struct database_file_reph_jack *jack_files)
{
    double *tmp;
    
     if (fabs(jack_files[ensemble].a-0.0619)<0.000001){
         tmp=fake_jack(7.60,0.008,jack_files[ensemble].Njack,123);
     }
     if (fabs(jack_files[ensemble].a-0.0815)<0.000001){
        tmp=fake_jack(5.77,0.006,jack_files[ensemble].Njack,1234);
    }
     if (fabs(jack_files[ensemble].a-0.0885)<0.000001){
        tmp=fake_jack(5.31,0.008,jack_files[ensemble].Njack,12345);
     }
     
     return tmp;
}


void  read_files_jack( struct database_file_reph_jack *jack_files, struct header *head,int ***mass_index, struct reph_jack **dataJ){
      
      int i,ikt,iks;
      int iG,iG_max;
      int im,im_max;
      int imom0,imoms,imomt;
      int e; 
      
      (*dataJ)=(struct reph_jack *) malloc(sizeof(struct reph_jack)*ensembles_reph);

      for(e=0;e<ensembles_reph;e++){
            iG_max=index_ensemble_twopt_G_fit(head[e],head[e].nk-1,  head[e].nk-1,  head[e].nmoms-1,  head[e].nmoms-1,  head[e].nmoms-1); 
            im_max=index_ensemble_twopt_fit(head[e],  head[e].nk-1,head[e].nk-1,head[e].nmoms-1,head[e].nmoms-1);      

            (*dataJ)[e].FA_from_H0=(double**) malloc(sizeof(double*)*(iG_max+1));
            (*dataJ)[e].FV_from_H0=(double**) malloc(sizeof(double*)*(iG_max+1));
            (*dataJ)[e].FV_from_H0_HA=(double**) malloc(sizeof(double*)*(iG_max+1));
            (*dataJ)[e].FAp=(double**) malloc(sizeof(double*)*(iG_max+1));
            (*dataJ)[e].FA=(double**) malloc(sizeof(double*)*(iG_max+1));
            (*dataJ)[e].M_PS=(double**) malloc(sizeof(double*)*(im_max+1));
            (*dataJ)[e].Zf_PS=(double**) malloc(sizeof(double*)*(im_max+1));
            (*dataJ)[e].FV=(double**) malloc(sizeof(double*)*(iG_max+1));
            (*dataJ)[e].xG=(double**) malloc(sizeof(double*)*(iG_max+1));

            for(iks=0;iks<3;iks++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ikt=0;ikt<head[e].nk;ikt++){
            for(imoms=0;imoms<head[e].nmoms;imoms++){
            for(imomt=0;imomt<head[e].nmoms;imomt++){
            for(imom0=0;imom0<head[e].nmoms;imom0++){       

                iG=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);   
                im=index_ensemble_twopt_fit(head[e], ikt,iks,imom0,imoms);      
                (*dataJ)[e].FA_from_H0[iG]=(double*) malloc(sizeof(double)*jack_files[e].Njack);
                (*dataJ)[e].FV_from_H0[iG]=(double*) malloc(sizeof(double)*jack_files[e].Njack);
                (*dataJ)[e].FV_from_H0_HA[iG]=(double*) malloc(sizeof(double)*jack_files[e].Njack);
	      	    (*dataJ)[e].FAp[iG]=(double*) malloc(sizeof(double)*jack_files[e].Njack);
                (*dataJ)[e].FA[iG]=(double*) malloc(sizeof(double)*jack_files[e].Njack);
                (*dataJ)[e].FV[iG]=(double*) malloc(sizeof(double)*jack_files[e].Njack);
                (*dataJ)[e].xG[iG]=(double*) malloc(sizeof(double)*jack_files[e].Njack);
                
                fread((*dataJ)[e].FA_from_H0[iG],   sizeof(double),   jack_files[e].Njack,   jack_files[e].f_FA_from_H0 );
                fread((*dataJ)[e].FV_from_H0[iG],   sizeof(double),   jack_files[e].Njack,   jack_files[e].f_FV_from_H0 );
                fread((*dataJ)[e].FV_from_H0_HA[iG],   sizeof(double),   jack_files[e].Njack,   jack_files[e].f_FV_from_H0_HA );
                fread((*dataJ)[e].FAp[iG],   sizeof(double),   jack_files[e].Njack,   jack_files[e].f_FAp );
                fread((*dataJ)[e].FA[iG],   sizeof(double),   jack_files[e].Njack,   jack_files[e].f_FA );
                fread((*dataJ)[e].FV[iG],        sizeof(double),   jack_files[e].Njack,  jack_files[e].f_FV );
                fread((*dataJ)[e].xG[iG],        sizeof(double),   jack_files[e].Njack,  jack_files[e].f_xG );
                if (imomt==0){
                    (*dataJ)[e].M_PS[im]=(double*)  malloc(sizeof(double)*jack_files[e].Njack);
                    (*dataJ)[e].Zf_PS[im]=(double*) malloc(sizeof(double)*jack_files[e].Njack);
                    
                    fread((*dataJ)[e].M_PS[im],   sizeof(double),   jack_files[e].Njack,   jack_files[e].f_M_PS );
                    fread((*dataJ)[e].Zf_PS[im],   sizeof(double),   jack_files[e].Njack,   jack_files[e].f_Zf_PS );
                }
                
            }}}}}
             //(*dataJ)[e].Zp=create_fake_jack_Zp(e, jack_files);
             //(*dataJ)[e].w0=create_fake_jack_w0(e, jack_files);
             //(*dataJ)[e].ZA=create_fake_jack_ZA(e, jack_files);
             //(*dataJ)[e].ZV=create_fake_jack_ZV(e, jack_files);
            
      }
    
}


void correct_FVE_pion(struct database_file_reph_jack *jack_files, struct header *head,int jack_tot, struct reph_jack **gJ){
    int e,im,iks,ikt,imoms,imom0,iG,iE,j ,Lsize;
    double mw,w0,Mpi,fpi, KM,Kf;
    
    for (e=0;e<ensembles_reph;e++){
          im=index_ensemble_twopt_fit(head[e], 0,0,0,0);      
          for (j=0;j<jack_tot;j++){
                mw=head[e].k[head[e].nk]/(*gJ)[e].Zp[j];
                w0=(*gJ)[e].w0[j];
                Mpi=(*gJ)[e].M_PS[im][j];
                fpi=(*gJ)[e].Zf_PS[im][j];
                Lsize=head[e].l0;
                
                FVE(  2.249  /*r_0 in GeV-1*/, w0, 3.772  /*l3b*/, 4.911 /*l4b*/, 2.508*2.249 /*Bw0*/,     0.1209 *2.249, Lsize, mw ,  Mpi*Mpi,  fpi, &KM, &Kf);
                //FVE(  2.249  /*r_0 in GeV-1*/, 5.23 , 3.772  /*l3b*/, 4.911 /*l4b*/, 2.508* 2.249,  0.1209 * 2.249, 24, 0.0040*5.23/0.510  , 0.5810,  0.3370, &KM, &Kf);
                for(iks=0;iks<1;iks++){     //for(ik1=0;ik1<=ik2;ik1++){
                for(ikt=0;ikt<1;ikt++){
                for(imoms=0;imoms<head[e].nmoms;imoms++){
                for(imom0=0;imom0<head[e].nmoms;imom0++){       
                    iE=index_ensemble_twopt_fit(head[e], ikt,iks,imom0,imoms);      
                    
                            printf("KM2 =%.7g  Kf=%g\n",KM*KM,Kf);
                            (*gJ)[e].M_PS[iE][j]*=KM;
                            (*gJ)[e].Zf_PS[iE][j]*=Kf;
                    
                    
                }}}}
          }
    }

}

struct reph_jack *create_generalised_jack( struct database_file_reph_jack *jack_files, struct header *head,int *jack_tot,int ***mass_index, struct reph_jack **dJ){
      int j,e,e1,iks,ikt,counter;
      double ***M_PS_GEVP_jack_tot;
      struct reph_jack *gJ;
      int iG,iG_max;
      int imom0,imoms,imomt;
      int im,im_max;

       
   
      gJ=(struct reph_jack *) malloc (sizeof(struct reph_jack )*ensembles_reph);
      
      *jack_tot=0;
      for(e=0;e<ensembles_reph;e++){
          *jack_tot+=jack_files[e].Njack;
      }
      *jack_tot=*jack_tot-ensembles_reph+1;
      
      for(e=0;e<ensembles_reph;e++){
          iG_max=index_ensemble_twopt_G_fit(head[e],head[e].nk-1,  head[e].nk-1,  head[e].nmoms-1,  head[e].nmoms-1,  head[e].nmoms-1); 
          im_max=index_ensemble_twopt_fit(head[e],  head[e].nk-1,  head[e].nk-1,head[e].nmoms-1,head[e].nmoms-1);      

          gJ[e].FA_from_H0=(double**) malloc(sizeof(double*)*(iG_max+1));
          gJ[e].FV_from_H0=(double**) malloc(sizeof(double*)*(iG_max+1));
          gJ[e].FV_from_H0_HA=(double**) malloc(sizeof(double*)*(iG_max+1));
          gJ[e].FAp=(double**) malloc(sizeof(double*)*(iG_max+1));
          gJ[e].FA=(double**) malloc(sizeof(double*)*(iG_max+1));
          gJ[e].M_PS=(double**) malloc(sizeof(double*)*(im_max+1));
          gJ[e].Zf_PS=(double**) malloc(sizeof(double*)*(im_max+1));
          gJ[e].FV=(double**) malloc(sizeof(double*)*(iG_max+1));
          gJ[e].xG=(double**) malloc(sizeof(double*)*(iG_max+1));
         
          gJ[e].Zp=(double*) malloc(sizeof(double)*(*jack_tot));//jack_tot is a pointer here
          gJ[e].w0=(double*) malloc(sizeof(double)*(*jack_tot));//jack_tot is a pointer here
          gJ[e].ZA=(double*) malloc(sizeof(double)*(*jack_tot));//jack_tot is a pointer here
          gJ[e].ZV=(double*) malloc(sizeof(double)*(*jack_tot));//jack_tot is a pointer here
          gJ[e].ZAV=(double*) malloc(sizeof(double)*(*jack_tot));//jack_tot is a pointer here
          for(iks=0;iks<3;iks++){     //for(ik1=0;ik1<=ik2;ik1++){
          for(ikt=0;ikt<head[e].nk;ikt++){
          for(imoms=0;imoms<head[e].nmoms;imoms++){
          for(imomt=0;imomt<head[e].nmoms;imomt++){
          for(imom0=0;imom0<head[e].nmoms;imom0++){       
                 iG=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);   

                 gJ[e].FA_from_H0[iG]=(double*) calloc(*jack_tot,sizeof(double));
                 gJ[e].FV_from_H0[iG]=(double*) calloc(*jack_tot,sizeof(double));
                 gJ[e].FV_from_H0_HA[iG]=(double*) calloc(*jack_tot,sizeof(double));
                 gJ[e].FAp[iG]=(double*)        calloc(*jack_tot,sizeof(double));
                 gJ[e].FA[iG]=(double*)         calloc(*jack_tot,sizeof(double));
                 gJ[e].FV[iG]=(double*)         calloc(*jack_tot,sizeof(double));
                 gJ[e].xG[iG]=(double*)         calloc(*jack_tot,sizeof(double));
                  if(e==3 && iG==0)
                  if (e==3 && iG==0) printf("allocating e=%d   iG=%d   jack_tot=%d   , jack_files[e].Njack-1=%d\n", e,iG,*jack_tot,jack_files[e].Njack-1);
                 if (imomt==0){
                    im=index_ensemble_twopt_fit(head[e], ikt,iks,imom0,imoms);
                                         error(im>im_max,1,"error ","im>im_max");
                    gJ[e].M_PS[im]=(double*) calloc(*jack_tot,sizeof(double));
                    gJ[e].Zf_PS[im]=(double*) calloc(*jack_tot,sizeof(double));
                 }
                 
                 counter=0;
                 for(e1=0;e1<ensembles_reph;e1++){
                      for(j=0;j<(jack_files[e1].Njack-1);j++){
                           if (e==e1){
                                gJ[e].FA_from_H0[iG][j+counter]=(*dJ)[e].FA_from_H0[  iG ][j];
                                gJ[e].FV_from_H0[iG][j+counter]=(*dJ)[e].FV_from_H0[  iG ][j];
                                gJ[e].FV_from_H0_HA[iG][j+counter]=(*dJ)[e].FV_from_H0_HA[  iG ][j];
                                gJ[e].FAp[iG][j+counter]=(*dJ)[e].FAp[  iG ][j];
                                gJ[e].FA[iG][j+counter]=(*dJ)[e].FA[  iG ][j];
                                gJ[e].FV[iG][j+counter]=(*dJ)[e].FV[  iG ][j];
                                gJ[e].xG[iG][j+counter]=(*dJ)[e].xG[  iG ][j];
                                if (imomt==0){
                                    gJ[e].M_PS[im][j+counter]=(*dJ)[e].M_PS[  im ][j];
                                    gJ[e].Zf_PS[im][j+counter]=(*dJ)[e].Zf_PS[  im ][j];
                                }
                            
                           }
                           else{ 
                                gJ[e].FA_from_H0[iG][j+counter]=(*dJ)[e].FA_from_H0[iG][ jack_files[e].Njack-1 ];
                                gJ[e].FV_from_H0[iG][j+counter]=(*dJ)[e].FV_from_H0[iG][ jack_files[e].Njack-1 ];
                                gJ[e].FV_from_H0_HA[iG][j+counter]=(*dJ)[e].FV_from_H0_HA[iG][ jack_files[e].Njack-1 ];
                                gJ[e].FAp[iG][j+counter]=(*dJ)[e].FAp[iG][ jack_files[e].Njack-1 ];
                                gJ[e].FA[iG][j+counter]=(*dJ)[e].FA[iG][ jack_files[e].Njack-1 ];
                                gJ[e].FV[iG][j+counter]=(*dJ)[e].FV[iG][ jack_files[e].Njack-1 ];   
                                gJ[e].xG[iG][j+counter]=(*dJ)[e].xG[iG][ jack_files[e].Njack-1 ];
                                if (imomt==0){
                                    gJ[e].M_PS[im][j+counter]=(*dJ)[e].M_PS[im][ jack_files[e].Njack-1 ];
                                    gJ[e].Zf_PS[im][j+counter]=(*dJ)[e].Zf_PS[im][ jack_files[e].Njack-1 ];
                                }
                            
                           }
                      }
                      counter+=jack_files[e1].Njack-1;
                 }
                 gJ[e].FA_from_H0[iG][*jack_tot-1]=(*dJ)[e].FA_from_H0[iG][jack_files[e].Njack-1];
                 gJ[e].FV_from_H0[iG][*jack_tot-1]=(*dJ)[e].FV_from_H0[iG][jack_files[e].Njack-1];
                 gJ[e].FV_from_H0_HA[iG][*jack_tot-1]=(*dJ)[e].FV_from_H0_HA[iG][jack_files[e].Njack-1];
                 gJ[e].FAp[iG][*jack_tot-1]=(*dJ)[e].FAp[iG][jack_files[e].Njack-1];
                 gJ[e].FA[iG][*jack_tot-1]=(*dJ)[e].FA[iG][jack_files[e].Njack-1];
                 gJ[e].FV[iG][*jack_tot-1]=(*dJ)[e].FV[iG][jack_files[e].Njack-1];
                 gJ[e].xG[iG][*jack_tot-1]=(*dJ)[e].xG[iG][jack_files[e].Njack-1]; 
                 if (imomt==0){
                    gJ[e].M_PS[im][*jack_tot-1]=(*dJ)[e].M_PS[im][jack_files[e].Njack-1];
                    gJ[e].Zf_PS[im][*jack_tot-1]=(*dJ)[e].Zf_PS[im][jack_files[e].Njack-1];
                 }
          }}}}}

      }
      
     double *w0A,*w0B,*w0D, *ZpA,*ZpB,*ZpD,*ZVA,*ZVB,*ZVD,*ZAA,*ZAB,*ZAD,*ZAVA,*ZAVB,*ZAVD;
      
     w0D=fake_jack(r0D,0.008,*jack_tot,123);
     w0B=fake_jack(r0B,0.006,*jack_tot,1234);
     w0A=fake_jack(r0A,0.008,*jack_tot,12345);
     
     
        //ZpD=fake_jack(0.516,0.002,*jack_tot,321);//M1 2Gev 
     ZpD=fake_jack(0.545,0.002,*jack_tot,321);//M2 2Gev 
        //ZpB=fake_jack(0.509,0.004,*jack_tot,3214);//M1 2Gev  
     ZpB=fake_jack(0.546,0.002,*jack_tot,3214);//M2 2Gev  
        //ZpA=fake_jack(0.529,0.007,*jack_tot,32145);//M1 2Gev  
     ZpA=fake_jack(0.574,0.004,*jack_tot,32145);//M2 2Gev  
     
    
         //tmp=fake_jack(0.762,0.004,*jack_tot);//M1 2Gev 
     ZAD=fake_jack(0.752,0.002,*jack_tot,213);//M2 2Gev 
        //tmp=fake_jack(0.737,0.005,*jack_tot);//M1 2Gev  
     ZAB=fake_jack(0.714,0.002,*jack_tot,2134);//M2 2Gev  
            //tmp=fake_jack(0.531,0.008,*jack_tot);//M1 2Gev  
     ZAA=fake_jack(0.703,0.002,*jack_tot,21345);//M2 2Gev  
    
    
     
     ZVD=fake_jack(0.6531,0.0002,*jack_tot,312);//WTI
     ZVB=fake_jack(0.6095,0.0003,*jack_tot,3124);//WTI
     ZVA=fake_jack(0.5920,0.0004,*jack_tot,3125);//WTI
     
     ZAVD=fake_jack(1.1434,0.0021,*jack_tot,111);//M2
     ZAVB=fake_jack(1.1638,0.0026,*jack_tot,222);//M2
     ZAVA=fake_jack(1.1560,0.0032,*jack_tot,333);//M2
     
     for(j=0;j<(*jack_tot);j++)
     for(e=0;e<ensembles_reph;e++){
     if (e>=0 && e<3){
        gJ[e].w0[j]=w0D[j];
        gJ[e].Zp[j]=ZpD[j];
        gJ[e].ZA[j]=ZAD[j];
        gJ[e].ZV[j]=ZVD[j];
        gJ[e].ZAV[j]=ZAVD[j];
      }
      else if (e>=3 && e<7){
        gJ[e].w0[j]=w0B[j];
        gJ[e].Zp[j]=ZpB[j];
        gJ[e].ZA[j]=ZAB[j];
        gJ[e].ZV[j]=ZVB[j];
        gJ[e].ZAV[j]=ZAVB[j];
      }
      else if (e>=7){
        gJ[e].w0[j]=w0A[j];
        gJ[e].Zp[j]=ZpA[j];
        gJ[e].ZA[j]=ZAA[j];
        gJ[e].ZV[j]=ZVA[j];
        gJ[e].ZAV[j]=ZAVA[j];
      }
     }
      int seedw0=412;
      result.w0fm=fake_jack(0.474,0.014,*jack_tot,seedw0);
      result.w0MeV=fake_jack(0.474/197.326963,0.014/197.326963 ,*jack_tot,seedw0);
      result.mlw=fake_jack(3.70*0.474/197.326963,0.17*0.474/197.326963 ,*jack_tot,4123);
      result.msw=fake_jack(99.6*0.474/197.326963,4.3*0.474/197.326963 ,*jack_tot,41235);
      result.mcw=fake_jack(1176*0.474/197.326963,39*0.474/197.326963 ,*jack_tot,521);

      result.MpiMeV=fake_jack(v_MpiMeV,err_MpiMeV,*jack_tot,5213);
      result.MKMeV=fake_jack(v_MKMeV,err_MKMeV,*jack_tot,5215);
      result.MDMeV=fake_jack(v_MDMeV,err_MDMeV,*jack_tot,621);
      result.MDsMeV=fake_jack(v_MDsMeV,err_MDsMeV,*jack_tot,6213);
      
      result.fpiMeV_exp=fake_jack(v_fpiMeV_exp,err_fpiMeV_exp,*jack_tot,62134);
      result.fKMeV_exp=fake_jack(v_fKMeV_exp,err_fKMeV_exp,*jack_tot,62135);
      
      for(e=0;e<ensembles_reph;e++){
          for(iks=0;iks<3;iks++){     //for(ik1=0;ik1<=ik2;ik1++){
          for(ikt=0;ikt<head[e].nk;ikt++){
          for(imoms=0;imoms<head[e].nmoms;imoms++){
          for(imomt=0;imomt<head[e].nmoms;imomt++){
          for(imom0=0;imom0<head[e].nmoms;imom0++){       
               iG=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);   
               im=index_ensemble_twopt_fit(head[e], ikt,iks,imom0,imoms);      
               free((*dJ)[e].FA_from_H0[iG]);
               free((*dJ)[e].FV_from_H0[iG]);
               free((*dJ)[e].FV_from_H0_HA[iG]);
               free((*dJ)[e].FAp[iG]);
               free((*dJ)[e].FA[iG]);
               free((*dJ)[e].FV[iG]);
               free((*dJ)[e].xG[iG]);
               if (imomt==0){
                    free((*dJ)[e].M_PS[ im ]);
                    free((*dJ)[e].Zf_PS[ im ]);
               }
          }}}}}
          free((*dJ)[e].FAp); free((*dJ)[e].FV);  free((*dJ)[e].FA_from_H0); free((*dJ)[e].FA); free((*dJ)[e].xG); free((*dJ)[e].M_PS); free((*dJ)[e].Zf_PS); 
          free((*dJ)[e].FV_from_H0);  free((*dJ)[e].FV_from_H0_HA);
      }
      free(*dJ);            
      
     free(w0A);free(w0B);free(w0D); free(ZpA);free(ZpB);free(ZpD);free(ZVA);free(ZVB);free(ZVD);free(ZAA);free(ZAB);free(ZAD);free(ZAVA);free(ZAVB);free(ZAVD);    
     
     //correct_FVE_pion(jack_files, head, *jack_tot, &gJ);
     return gJ;
}


void free_data( struct database_file_reph_jack **jack_files, struct header **head,int *jack_tot, struct reph_jack **grephJ){
      
    int e,i;
    int imoms,imomt,imom0,iks,ikt;
    int iG,im;
    
    for(e=0;e<ensembles_reph;e++){
          for(iks=0;iks<3;iks++){     //for(ik1=0;ik1<=ik2;ik1++){
          for(ikt=0;ikt<(*head)[e].nk;ikt++){
          for(imoms=0;imoms<(*head)[e].nmoms;imoms++){
          for(imomt=0;imomt<(*head)[e].nmoms;imomt++){
          for(imom0=0;imom0<(*head)[e].nmoms;imom0++){       
               iG=index_ensemble_twopt_G_fit((*head)[e],ikt,iks,imom0,imomt,imoms);   
               im=index_ensemble_twopt_fit((*head)[e], ikt,iks,imom0,imoms);      
               free((*grephJ)[e].FA_from_H0[ iG ]);
               free((*grephJ)[e].FV_from_H0[ iG ]);
               free((*grephJ)[e].FV_from_H0_HA[ iG ]);
               free((*grephJ)[e].FAp[ iG ]);
               free((*grephJ)[e].FA[ iG ]);
               free((*grephJ)[e].FV[ iG ]);
               free((*grephJ)[e].xG[ iG ]);
               if (imomt==0){
                    free((*grephJ)[e].M_PS[ im ]);
                    free((*grephJ)[e].Zf_PS[ im ]);
               }
          }}}}}
          free((*grephJ)[e].FAp); free((*grephJ)[e].FV);  free((*grephJ)[e].FA_from_H0); free((*grephJ)[e].FA); free((*grephJ)[e].xG); free((*grephJ)[e].M_PS); free((*grephJ)[e].Zf_PS);
          free((*grephJ)[e].w0); free((*grephJ)[e].Zp);  free((*grephJ)[e].ZV); free((*grephJ)[e].ZA);   free((*grephJ)[e].FV_from_H0); free((*grephJ)[e].FV_from_H0_HA);
          free((*grephJ)[e].ZAV);
    }
    free((*grephJ));   
     ///////////////////////free header///////////////
    for (e=0;e<ensembles_reph;e++){
        free((*head)[e].k);
        for (i=0;i<(*head)[e].nmoms;i++)
            free((*head)[e].mom[i]);
        free((*head)[e].mom);
    }
    free((*head));
    for (e=0;e<ensembles_reph;e++){

             fclose((*jack_files)[e].f_FA_from_H0);
             fclose((*jack_files)[e].f_FV_from_H0);
             fclose((*jack_files)[e].f_FV_from_H0_HA);
             fclose((*jack_files)[e].f_FAp);
             fclose((*jack_files)[e].f_FA);
             fclose((*jack_files)[e].f_FV);
             fclose((*jack_files)[e].f_xG);
             fclose((*jack_files)[e].f_M_PS);
             fclose((*jack_files)[e].f_Zf_PS);
    }
    free((*jack_files)); 
  
}
double **init_phys_point(int jack_tot){
    int j;
    double **phys_point=(double**) malloc(sizeof(double*)*jack_tot);

    for(j=0;j<jack_tot;j++){
        phys_point[j]=(double*) malloc(sizeof(double)*13);
        phys_point[j][0]=result.mlw[j];
        phys_point[j][1]=1e+10;
        phys_point[j][2]=0.;
        phys_point[j][3]=result.MpiMeV[j]*result.w0MeV[j];
        phys_point[j][4]=result.msw[j];
        phys_point[j][5]=result.MKMeV[j]*result.w0MeV[j];
        phys_point[j][6]=result.mcw[j];
        phys_point[j][7]=result.MDMeV[j]*result.w0MeV[j];
        phys_point[j][8]=result.MDsMeV[j]*result.w0MeV[j];
        phys_point[j][9]=0;
        phys_point[j][10]=0;
        phys_point[j][11]=result.fpiMeV_exp[j]*result.w0MeV[j];
        phys_point[j][12]=result.fKMeV_exp[j]*result.w0MeV[j];
    }
    printf("physical continuum point:\n");
    printf("Mpiw0= %g  \t Mpi=%g MeV\n",result.MpiMeV[jack_tot-1]*result.w0MeV[jack_tot-1],result.MpiMeV[jack_tot-1] );
    printf("MKw0= %g  \t MD=%g MeV\n",result.MKMeV[jack_tot-1]*result.w0MeV[jack_tot-1],result.MKMeV[jack_tot-1] );
    printf("MDw0= %g  \t MK=%g MeV\n",result.MDMeV[jack_tot-1]*result.w0MeV[jack_tot-1],result.MDMeV[jack_tot-1] );
    printf("MDsw0= %g  \t MDs=%g MeV\n",result.MDsMeV[jack_tot-1]*result.w0MeV[jack_tot-1],result.MDsMeV[jack_tot-1] );
    printf("r0=%g fm   \t   r0=%g MeV\n",result.w0fm[jack_tot-1],result.w0MeV[jack_tot-1]);

    return phys_point;
}
void    free_results(){
     
     free(result.w0fm);
     free(result.w0MeV);
     free(result.mlw);
     free(result.msw);
     free(result.mcw);
     free(result.MpiMeV);
     free(result.MKMeV);
     free(result.MDMeV);
     free(result.MDsMeV);
}
int main(int argc, char **argv){
    
    int i,j,e,Ns;
    struct header  *head;
    struct database_file_reph_jack  *jack_files;
    //double ***M_PS_GEVP_jack,***f_PS_jack;
    //double ***M_PS_GEVP_jack_tot,***f_PS_jack_tot;
    struct  reph_jack *grephJ;
    double *tmp,**fit;
    struct  reph_jack *rephJ;
    
    error(argc<5,1,"main ",
         "choose manual or auto plateaux \n usage: ./fit_all_reph   manual/auto   output_path_1  -I inpath_1       output_path_2  -I inpath_2  ...");

    error(strcmp(argv[1],"manual")!=0 && strcmp(argv[1],"auto")!=0 ,2,"main ",
         "choose manual or auto plateaux \n usage: ./fit_all_beta   manual/auto");
    error(strcmp(argv[3],"-I")!=0,1,"main ", "choose manual or auto plateaux \n usage: ./fit_all_reph   manual/auto   output_path  -I inpa.hpp");
    
    error(((argc+1)%3)!=0,1,"main","\n usage: ./fit_all_reph   manual/auto   output_path_1  -I inpath_1       output_path_2  -I inpath_2  ...");
    int Nsets=((argc+1)/3)-1;
    printf("analysing %d set(s) of data\n",Nsets);
    
    int im,ikt,iks;
    int imomt,imoms,imom0,iG;

    struct fit_type fit_info;
    struct fit_result  fit_out;
    struct fit_all  fit_chi2_good ,fit_FA_FV, fit_FA_FV_pole;
    fit_chi2_good.Nfits=0;
    fit_info.Nvar=13;
    fit_FA_FV.Nfits=0;
    fit_FA_FV_pole.Nfits=0;

    char *argvNs[5];
    double **phys_point;
    
    files_declarations(argv ,&jack_files,&head);
    read_files_jack(jack_files,head,mass_index,&rephJ);
    grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
    
    
    printf("ensambles\n");
    for (e=0;e<ensembles_reph;e++){
        printf("L%dT%d  %s beta=%.3f  aM_PS=%f    aM_PS=%f MeV ",head[e].l1,head[e].l0,jack_files[e].FA_from_H0,head[e].beta,grephJ[e].M_PS[0][jack_tot-1] ,grephJ[e].M_PS[0][jack_tot-1]*grephJ[e].w0[jack_tot-1]/result.w0MeV[jack_tot-1]);
        printf("r0=%f  m_lr_0=%f mu=",grephJ[e].w0[jack_tot-1], grephJ[e].w0[jack_tot-1]*head[e].k[head[e].nk]/grephJ[e].Zp[jack_tot-1] );
        for(ikt=0;ikt<head[e].nk;ikt++)
            printf("%.5f   ", head[e].k[head[e].nk+ikt] );
        printf("\n");
        printf("Njack=%d \n",jack_files[e].Njack);
    }
    printf("\n");
    
    phys_point=init_phys_point(jack_tot);  
    error(jack_tot!=139,3,"main","The program can not combine the systematics of set of different jack_tot lenght");
    free_data( &jack_files, &head,&jack_tot, &grephJ);

  /*
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

    files_declarations(argvNs ,&jack_files,&head);
    read_files_jack(jack_files,head,mass_index,&rephJ);
    
    grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
    //M_PS_GEVP_jack_tot=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, rephJ.M_PS_GEVP_jack);//every time jack_tot is initialised
    //f_PS_jack_tot=create_generalised_jack( jack_files, head, &jack_tot,mass_index, rephJ.f_PS_jack);
   
    //printf("ensambles\n");
    //for (e=0;e<ensembles_reph;e++){
        //printf("L%dT%d  %s beta=%.3f  M_PS=%f  ",head[e].l1,head[e].l0,jack_files[e].FA_from_H0,head[e].beta,grephJ[e].M_PS[0][jack_tot-1] );
        //printf("r0=%f  m_lr_0=%f mu=",grephJ[e].w0[jack_tot-1], grephJ[e].w0[jack_tot-1]*head[e].k[head[e].nk]/grephJ[e].Zp[jack_tot-1] );
        //for(ikt=0;ikt<head[e].nk;ikt++)
            //printf("%.5f   ", head[e].k[head[e].nk+ikt] );
        //printf("\n");
        //printf("Njack=%d \n",jack_files[e].Njack);
    //}
    //printf("\n");
    
    
    
    
    
    
    ////////////////////////FA
    
   // result.fw=(double*)  malloc(sizeof(double) *jack_tot);
 
    
      
    printf("///////////////////////////////////////ChPT FA///////////////////////\n");
    
    printf("\n\ntrying generic fit\n");
    fit_info.Npar=2;
    fit_info.N=1;
    fit_info.function=FA_FV_ChPT;
    
    
  
    //fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
    //reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"A","FA_pion_ChPT",head , grephJ );       
    //print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_pion_ChPT");


    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);

    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_pion_ChPT");

   
    if(2+2==5){
                printf("\n\n///////////////////////////////////////ChPT FA  NNLO///////////////////////\n");
                printf("ChPT NNLO\n");
                fit_info.Npar=4;
                fit_info.N=1;
                fit_info.function=FA_FV_ChPT_NNLO;
                
                fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
                print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_pion_ChPT_NNLO");

                
                printf("\n\n///////////////////////////////////////ChPT xg*FA  NNLO///////////////////////\n");
                printf("ChPT NNLO   xg*F\n");
                fit_info.Npar=4;
                fit_info.N=1;
                fit_info.function=FA_FVxg_ChPT_NNLO;
                
                fit_out=fit_FAxg_pion_generic(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
                print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FAxg_ChPT_NNLO");


                printf("\n\n///////////////////////////////////////ChPT FA  pole///////////////////////\n");
                printf("ChPT pole\n");
                fit_info.Npar=3;
                fit_info.N=1;
                fit_info.function=FA_FV_ChPT_pole;
                
                    
                fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
                print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_ChPT_pole");
    }
 

    printf("\n\n///////////////////////////////////////ChPT FA pole1  ///////////////////////\n");
    printf("ChPT pole\n");
    fit_info.Npar=4;
    fit_info.N=1;
    fit_info.function=FA_FV_pole;
    
    
    //fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
    //print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_pole");
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_pole");
 
//    printf("\n\n///////////////////////////////////////ChPT HA pole1  ///////////////////////\n");
//    printf("ChPT pole\n");
//    fit_info.Npar=8;
//    fit_info.N=1;
//    fit_info.function=HA_pole;
//    fit_info.f1=FA_FV_pole;
//    fit_info.f2=fpiw0_fit;
    
    
//   fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,2 ,fit_info);
//    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "H","HA_pole");

    if(2+2==5){
            printf("\n\n///////////////////////////////////////ChPT FA pole1 +xg ///////////////////////\n");
            printf("ChPT pole\n");
            fit_info.Npar=5;
            fit_info.N=1;
            fit_info.function=FA_FV_pole_xg_const;
            
            
            fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_pole_xg_const");
    }
    printf("\n\n/////////////////////////////////////// FA_from_H0 poly2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=4;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly2;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly2");
    
    printf("\n\n/////////////////////////////////////// FA_from_H0 poly3  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=5;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly3;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly3");
   
    printf("\n\n/////////////////////////////////////// FA_from_H0 poly2 p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=6;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly2_p2k2;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly2_p2k2");
    
    printf("\n\n/////////////////////////////////////// FA_from_H0 poly2 p2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=5;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly2_p2;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly2_p2");
    
    printf("\n\n/////////////////////////////////////// FA_from_H0 poly2 k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=5;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly2_k2;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly2_k2");
    
    printf("\n\n/////////////////////////////////////// FA_from_H0 poly2 p_k  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=6;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly2_p_k;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly2_p_k");
  
  
    printf("\n\n/////////////////////////////////////// FA_from_H0 poly2 M_Pi p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=8;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly2_M_Pi_p2k2;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly2_M_Pi_p2k2");
    
    printf("\n\n/////////////////////////////////////// FA_from_H0 poly2 M_Pi p2k2 kp ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=9;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly2_M_Pi_p2k2_kp;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly2_M_Pi_p2k2_kp");
    if (2+2==5){
                    printf("\n\n/////////////////////////////////////// FA_from_H0 poly2 M_Pi p2k2x kx ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=9;
                    fit_info.N=1;
                    fit_info.function=FA_FV_Pi_poly2_M_Pi_p2k2x_kx;
                    
                    
                    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly2_M_Pi_p2k2x_kx");
                    printf("\n\n/////////////////////////////////////// FA_from_H0 poly2 M_Pi p2k2x  ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=8;
                    fit_info.N=1;
                    fit_info.function=FA_FV_Pi_poly2_M_Pi_p2k2x;
                    
                    
                    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly2_M_Pi_p2k2x");
                    printf("\n\n/////////////////////////////////////// FA_from_H0 poly2 M_Pi p2kx  ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=8;
                    fit_info.N=1;
                    fit_info.function=FA_FV_Pi_poly2_M_Pi_p2kx;
                    
                    
                    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly2_M_Pi_p2kx");
        
    }
    printf("\n\n/////////////////////////////////////// FA_from_H0 poly3 M_Pi p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=9;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly3_M_Pi_p2k2;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly3_M_Pi_p2k2");
  
    printf("\n\n/////////////////////////////////////// FA_from_H0 poly2 M_Pi p2k2 p4k4 ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=10;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly2_M_Pi_p2k2_p4k4;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_poly2_M_Pi_p2k2_p4k4");
    
          ///////////////////////free///////////////
    free_data( &jack_files, &head,&jack_tot, &grephJ);

    }
    
      
    
 
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FA_Pion");
    fit_chi2_good.Nfits=0;
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

    files_declarations(argvNs ,&jack_files,&head);
    read_files_jack(jack_files,head,mass_index,&rephJ);
    grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);

    
    
    
    
    printf("\n\n/////////////////////////////////////// FV Pion///////////////////////\n");
    if(2+2==5){
                        printf("\n\n///////////////////////////////////////ChPT FV///////////////////////\n");
                        printf("trying generic fit\n");
                        fit_info.Npar=2;
                        fit_info.N=1;
                        fit_info.function=FA_FV_ChPT;
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_ChPT");
                        
                        printf("\n\n///////////////////////////////////////ChPT NNLO FV///////////////////////\n");
                        printf("ChPT NNLO\n");
                        fit_info.Npar=4;
                        fit_info.N=1;
                        fit_info.function=FA_FV_ChPT_NNLO;
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_ChPT_NNLO");
                    
                        printf("\n\n///////////////////////////////////////ChPT FV  pole///////////////////////\n");
                        printf("ChPT pole\n");
                        fit_info.Npar=3;
                        fit_info.N=1;
                        fit_info.function=FA_FV_ChPT_pole;
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_ChPT_pole");

                        printf("\n\n///////////////////////////////////////ChPT FV pole1  ///////////////////////\n");
                        printf("ChPT pole\n");
                        fit_info.Npar=4;
                        fit_info.N=1;
                        fit_info.function=FA_FV_pole;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_pole");
                        
                        printf("\n\n///////////////////////////////////////ChPT FV_from_H0 pole1  ///////////////////////\n");
                        printf("ChPT pole\n");
                        fit_info.Npar=4;
                        fit_info.N=1;
                        fit_info.function=FA_FV_pole;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_pole");
                        printf("\n\n///////////////////////////////////////ChPT FV pole1 +xg const  ///////////////////////\n");
                        printf("ChPT pole\n");
                        fit_info.Npar=5;
                        fit_info.N=1;
                        fit_info.function=FA_FV_pole_xg_const;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_pole_xg_const");

                        
                        printf("\n\n/////////////////////////////////////// FV poly2  ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=4;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly2;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_Pi_poly2");
                    
                        printf("\n\n/////////////////////////////////////// FV_from_H0 poly2  ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=4;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly2;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly2");
                        
                        printf("\n\n/////////////////////////////////////// FV_from_H0 poly3  ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=5;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly3;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly3");
                        
                        printf("\n\n/////////////////////////////////////// FV_from_H0 poly2 p2k2 ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=6;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly2_p2k2;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly2_p2k2");

                        printf("\n\n/////////////////////////////////////// FV_from_H0 poly2 M_Pi p2k2 ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=8;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly2_M_Pi_p2k2;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly2_M_Pi_p2k2");

                        
                        printf("\n\n/////////////////////////////////////// FV_from_H0 poly2 M_Pi p2 ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=7;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly2_M_Pi_p2;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly2_M_Pi_p2");
                        
                        printf("\n\n/////////////////////////////////////// FV_from_H0 poly2 M_Pi k2 ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=7;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly2_M_Pi_k2;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly2_M_Pi_k2");
                        
                        printf("\n\n/////////////////////////////////////// FV_from_H0 poly2 M_Pi p_k ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=8;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly2_M_Pi_p_k;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly2_M_Pi_p_k");



                        printf("\n\n/////////////////////////////////////// FV_from_H0 poly2 M_Pi p2k2 kp ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=9;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly2_M_Pi_p2k2_kp;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly2_M_Pi_p2k2_kp");
                        printf("\n\n/////////////////////////////////////// FV_from_H0 poly2 M_Pi p2k2 kx ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=9;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly2_M_Pi_p2k2x_kx;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly2_M_Pi_p2k2x_kx");
                        printf("\n\n/////////////////////////////////////// FV_from_H0 poly2 M_Pi p2k2x ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=8;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly2_M_Pi_p2k2x;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly2_M_Pi_p2k2x");
                        printf("\n\n/////////////////////////////////////// FV_from_H0 poly2 M_Pi p2kx ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=8;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly2_M_Pi_p2kx;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly2_M_Pi_p2kx");

                        printf("\n\n/////////////////////////////////////// FV_from_H0 poly2 M_Pi p2k2 p4k4 ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=10;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly2_M_Pi_p2k2_p4k4;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly2_M_Pi_p2k2_p4k4");

                    
                    printf("\n\n/////////////////////////////////////// FV_from_H0 poly3 M_Pi p2k2 ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=9;
                        fit_info.N=1;
                        fit_info.function=FA_FV_Pi_poly3_M_Pi_p2k2;
                        
                        
                        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Pi_poly3_M_Pi_p2k2");
    }
   
    printf("\n\n///////////////////////////////////////FV_from_H0_HA pole1  ///////////////////////\n");
    printf("ChPT pole\n");
    fit_info.Npar=4;
    fit_info.N=1;
    fit_info.function=FA_FV_pole;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_pole");
      
    printf("\n\n/////////////////////////////////////// FV_from_H0_HA poly2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=4;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly2;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_poly2");
    
    printf("\n\n/////////////////////////////////////// FV_from_H0_HA poly3  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=5;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly3;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_poly3");
    
    printf("\n\n/////////////////////////////////////// FV_from_H0_HA poly2 p2k2 ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=6;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly2_p2k2;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_poly2_p2k2");

    printf("\n\n/////////////////////////////////////// FV_from_H0_HA poly2 M_Pi p2k2 ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=8;
    fit_info.N=1;
    fit_info.function=FA_FV_Pi_poly2_M_Pi_p2k2;
    
    
    fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_poly2_M_Pi_p2k2");


    free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
 
 
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FV_Pion");
    fit_chi2_good.Nfits=0;
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

    files_declarations(argvNs ,&jack_files,&head);
    read_files_jack(jack_files,head,mass_index,&rephJ);
    grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);

    printf("///////////////////////////////////////FA K///////////////////////\n///////////////////////////////////////FA K///////////////////////\n///////////////////////////////////////FA K///////////////////////\n///////////////////////////////////////FA K///////////////////////\n///////////////////////////////////////FA K///////////////////////\n///////////////////////////////////////FA K///////////////////////\n");
    if(2+2==5){
                        printf("\n\n FA for K meson\n");
                        fit_info.Npar=4;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_ChPT;
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_K_ChPT");
                        
                        printf("\n\n//////////////////////////////////////\n FA for K meson\n////////////////////////////\n");
                        fit_info.Npar=8;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_NNLO_ChPT;
                        
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_K_NNLO_ChPT");
                        
                        printf("\n\n//////////////////////////////////////\n FA for K pole\n////////////////////////////\n");

                        fit_info.Npar=6;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_pole;
                        
                        
                        //fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
                        //print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_K_pole");
                        
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_pole");
                    
                        printf("\n\n//////////////////////////////////////\n HA for K pole\n////////////////////////////\n");
                    //    printf("ChPT pole\n");
                    //    fit_info.Npar=12;
                    //    fit_info.N=2;
                    //    fit_info.function=HA_K_pole;
                    //    fit_info.f1=FA_FV_K_pole;
                    //    fit_info.f2=fKw0_fit;
                        
                    
                    //    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,2 ,fit_info);
                    //    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "H","HA_K_pole");

                        printf("\n\n//////////////////////////////////////\n FA for K pole xg const\n////////////////////////////\n");

                        fit_info.Npar=7;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_pole_xg_const;
                        
                        
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_K_pole_xg_const");
    }
    printf("\n\n//////////////////////////////////////\n FA for K pole\n////////////////////////////\n");

    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=FA_FV_K_pole;
    
    
      
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_pole");
  
    printf("\n\n//////////////////////////////////////\n FA for K pole M_Pi_K\n////////////////////////////\n");

    fit_info.Npar=10;
    fit_info.N=2;
    fit_info.function=FA_FV_K_pole_M_Pi_K;
    
    
      
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_pole_M_Pi_K");
  
    printf("\n\n/////////////////////////////////////// FA_from_H0 K poly2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=FA_FV_K_poly2;
    
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_poly2");
    
    printf("\n\n/////////////////////////////////////// FA_from_H0 K poly3 p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=9;
    fit_info.N=2;
    fit_info.function=FA_FV_K_poly3_p2k2;
    
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_poly3_p2k2");
    
    printf("\n\n/////////////////////////////////////// FA_from_H0 K poly2 p2k2 ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=8;
    fit_info.N=2;
    fit_info.function=FA_FV_K_poly2_p2k2;
    
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_poly2_p2k2");
    
    printf("\n\n/////////////////////////////////////// FA_from_H0 K M_Pi_K poly2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=10;
    fit_info.N=2;
    fit_info.function=FA_FV_K_poly2_M_Pi_K;
    
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_poly2_M_Pi_K");
      

    printf("\n\n/////////////////////////////////////// FA_from_H0 K M_Pi_K poly2 p2k2 ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=12;
    fit_info.N=2;
    fit_info.function=FA_FV_K_poly2_M_Pi_K_p2k2;
    
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_poly2_M_Pi_K_p2k2");
      
    printf("\n\n/////////////////////////////////////// FA_from_H0 K poly3 p2 k2 ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=13;
    fit_info.N=2;
    fit_info.function=FA_FV_K_poly3_M_Pi_K_p2k2;
    
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_poly3_M_Pi_K_p2k2");
    
    
    free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
 
 
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FA_K");
    fit_chi2_good.Nfits=0;
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

    files_declarations(argvNs ,&jack_files,&head);
    read_files_jack(jack_files,head,mass_index,&rephJ);
    grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);

    
    printf("///////////////////////////////////////FV K///////////////////////\n///////////////////////////////////////FV K//////////////////\n");
    if (2+2==5){
                        printf("\n\n FV for K meson\n");
                        fit_info.Npar=4;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_ChPT;
                        
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_K_ChPT");
                        
                        
                        printf("\n\n//////////////////////////////////////\n FV for K meson\n////////////////////////////\n");
                        fit_info.Npar=8;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_NNLO_ChPT;
                        
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_K_NNLO_ChPT");

                        printf("\n\n//////////////////////////////////////\n FV for K pole\n////////////////////////////\n");

                        fit_info.Npar=6;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_pole;
                        
                        
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_K_pole");
                    
                        printf("\n\n//////////////////////////////////////\n FV_from_H0 for K pole\n////////////////////////////\n");

                        fit_info.Npar=6;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_pole;
                        
                        
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_K_pole");
                    
                        
                        //printf("\n\n//////////////////////////////////////\n FV for K pole xg const\n////////////////////////////\n");

                        //fit_info.Npar=7;
                        //fit_info.N=2;
                        //fit_info.function=FA_FV_K_pole_xg_const;
                        
                        
                        //fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                        //print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_K_pole_xg_const");
                        
                        printf("\n\n/////////////////////////////////////// FV K poly2  ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=6;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_poly2;
                        
                        
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_K_poly2");
                        
                        printf("\n\n/////////////////////////////////////// FV_from_H0 K poly2  ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=6;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_poly2;
                        
                        
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_K_poly2");
                        
                        //printf("\n\n/////////////////////////////////////// FV_from_H0 K poly3  ///////////////////////\n");
                        //printf("Polynomial\n");
                        //fit_info.Npar=7;
                        //fit_info.N=2;
                        //fit_info.function=FA_FV_K_poly3;
                        
                        
                        //fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        //print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_K_poly3");
                    
                        printf("\n\n/////////////////////////////////////// FV from H0 K poly2 p2k2  ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=8;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_poly2_p2k2;
                        
                        
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_K_poly2_p2k2");
                    
                        printf("\n\n/////////////////////////////////////// FV from H0 K poly2 M_Pi_K p2k2  ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=12;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_poly2_M_Pi_K_p2k2;
                        
                        
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_K_poly2_M_Pi_K_p2k2");
                    
                        printf("\n\n/////////////////////////////////////// FV from H0 K poly3 M_Pi_K p2k2  ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=13;
                        fit_info.N=2;
                        fit_info.function=FA_FV_K_poly3_M_Pi_K_p2k2;
                        
                        
                        fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_K_poly3_M_Pi_K_p2k2");
    }
    
    printf("\n\n//////////////////////////////////////\n FV_from_H0_HA for K pole\n////////////////////////////\n");

    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=FA_FV_K_pole;
    
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_pole");
   
    printf("\n\n//////////////////////////////////////\n FV_from_H0_HA for K pole M_Pi_K\n////////////////////////////\n");

    fit_info.Npar=10;
    fit_info.N=2;
    fit_info.function=FA_FV_K_pole_M_Pi_K;
    
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_pole_M_Pi_K");
   
    
     printf("\n\n/////////////////////////////////////// FV_from_H0_HA K poly2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=FA_FV_K_poly2;
    
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_poly2");
    
    printf("\n\n/////////////////////////////////////// FV from H0 HA K poly2 p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=8;
    fit_info.N=2;
    fit_info.function=FA_FV_K_poly2_p2k2;
    
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_poly2_p2k2");
  
    printf("\n\n/////////////////////////////////////// FV from H0 HA K poly3 p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=9;
    fit_info.N=2;
    fit_info.function=FA_FV_K_poly3_p2k2;
    
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_poly3_p2k2");
  
    
    printf("\n\n/////////////////////////////////////// FV from H0 HA K poly2 M_Pi_K   ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=10;
    fit_info.N=2;
    fit_info.function=FA_FV_K_poly2_M_Pi_K;
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_poly2_M_Pi_K");
    
    printf("\n\n/////////////////////////////////////// FV from H0 HA K poly2 M_Pi_K p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=12;
    fit_info.N=2;
    fit_info.function=FA_FV_K_poly2_M_Pi_K_p2k2;
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_poly2_M_Pi_K_p2k2");
    
    printf("\n\n/////////////////////////////////////// FV from H0 HA K poly3 M_Pi_K p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=13;
    fit_info.N=2;
    fit_info.function=FA_FV_K_poly3_M_Pi_K_p2k2;
    
    
    fit_out=fit_FAV_K(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_poly3_M_Pi_K_p2k2");

    
    
     free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
 
 
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FV_K");
    fit_chi2_good.Nfits=0;
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

    files_declarations(argvNs ,&jack_files,&head);
    read_files_jack(jack_files,head,mass_index,&rephJ);
    grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
     printf("///////////////////////////////////////FA D///////////////////////\n///////////////////////////////////////FA D///////////////////////\n///////////////////////////////////////FA D///////////////////////\n///////////////////////////////////////FA D///////////////////////\n///////////////////////////////////////FA D///////////////////////\n///////////////////////////////////////FA D///////////////////////\n");
  
    printf("\n\n FA for D meson\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=FA_FV_D_pole;
    
 //   fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
 //   print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_D_pole");
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_D_pole");
    
    
    printf("\n\n FA H for D meson pole M_Pi_D\n");
    fit_info.Npar=10;
    fit_info.N=2;
    fit_info.function=FA_FV_D_pole_M_Pi_D;
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_D_pole_M_Pi_D");
   
  //  printf("\n\n////////////\n HA for D pole\n/////////////\n");
//    printf("ChPT pole\n");
//    fit_info.Npar=12;
//    fit_info.N=2;
//    fit_info.function=HA_D_pole;
//    fit_info.f1=FA_FV_D_pole;
//    fit_info.f2=fDw0_fit;
        
//    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,2 ,fit_info);
//    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "H","HA_D_pole");
    
//    printf("\n\n FA for D meson +xg const\n");
//    fit_info.Npar=7;
//    fit_info.N=2;
//    fit_info.function=FA_FV_D_pole_xg_const;
    
//    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
//    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_D_pole_xg_const");

    printf("\n\n/////////////////////////////////////// FA_from_H0 D poly2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=FA_FV_D_poly2;
    
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_D_poly2");
  
    if(2+2==5){
                    printf("\n\n/////////////////////////////////////// FA_from_H0 D poly3  ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=7;
                    fit_info.N=2;
                    fit_info.function=FA_FV_D_poly3;
                    
                    
                    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_D_poly3");
    }
   
    printf("\n\n/////////////////////////////////////// FA_from_H0 D poly2 p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=8;
    fit_info.N=2;
    fit_info.function=FA_FV_D_poly2_p2k2;
    
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_D_poly2_p2k2");

    printf("\n\n/////////////////////////////////////// FA_from_H0 D poly3 p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=9;
    fit_info.N=2;
    fit_info.function=FA_FV_D_poly3_p2k2;
    
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_D_poly3_p2k2");

    printf("\n\n/////////////////////////////////////// FA_from_H0 D  M_Pi_D poly2   ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=10;
    fit_info.N=2;
    fit_info.function=FA_FV_D_poly2_M_Pi_D;
    
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_D_poly2_M_Pi_D");

    
    printf("\n\n/////////////////////////////////////// FA_from_H0 D  M_Pi_D poly2 p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=12;
    fit_info.N=2;
    fit_info.function=FA_FV_D_poly2_M_Pi_D_p2k2;
    
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_D_poly2_M_Pi_D_p2k2");

    printf("\n\n/////////////////////////////////////// FA_from_H0 D  M_Pi_D poly3 p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=13;
    fit_info.N=2;
    fit_info.function=FA_FV_D_poly3_M_Pi_D_p2k2;
    
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_D_poly3_M_Pi_D_p2k2");
    
     free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
 
 
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FA_D");
    fit_chi2_good.Nfits=0;
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

    files_declarations(argvNs ,&jack_files,&head);
    read_files_jack(jack_files,head,mass_index,&rephJ);
    grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
     printf("///////////////////////////////////////FV D///////////////////////\n///////////////////////////////////////FV D///////////////////////\n///////////////////////////////////////FV D///////////////////////\n///////////////////////////////////////FV D///////////////////////\n///////////////////////////////////////FV D///////////////////////\n///////////////////////////////////////FV D///////////////////////\n");
    if(2+2==5){
                            printf("\n\n FV for D meson\n");
                            fit_info.Npar=6;
                            fit_info.N=2;
                            fit_info.function=FA_FV_D_pole;
                            
                            fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_D_pole");
                            
                            printf("\n\n FV_from_H0 for D meson\n");
                            fit_info.Npar=6;
                            fit_info.N=2;
                            fit_info.function=FA_FV_D_pole;
                            
                            fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_D_pole");
                        
                            printf("\n\n FV for D meson +xg const\n");
                            fit_info.Npar=7;
                            fit_info.N=2;
                            fit_info.function=FA_FV_D_pole_xg_const;
                            
                            fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_D_pole_xg_const");
                            
                            printf("\n\n/////////////////////////////////////// FV D poly2  ///////////////////////\n");
                            printf("Polynomial\n");
                            fit_info.Npar=6;
                            fit_info.N=2;
                            fit_info.function=FA_FV_D_poly2;
                            
                            
                            fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_D_poly2");
                        
                            printf("\n\n/////////////////////////////////////// FV_from_H0 D poly2  ///////////////////////\n");
                            printf("Polynomial\n");
                            fit_info.Npar=6;
                            fit_info.N=2;
                            fit_info.function=FA_FV_D_poly2;
                            
                            
                            fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_D_poly2");
                            
                            printf("\n\n/////////////////////////////////////// FV_from_H0 D poly3  ///////////////////////\n");
                            printf("Polynomial\n");
                            fit_info.Npar=7;
                            fit_info.N=2;
                            fit_info.function=FA_FV_D_poly3;
                            
                            
                            fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_D_poly3");

                            printf("\n\n/////////////////////////////////////// FV_from_H0 D poly2 p2 ///////////////////////\n");
                            printf("Polynomial\n");
                            fit_info.Npar=7;
                            fit_info.N=2;
                            fit_info.function=FA_FV_D_poly2_p2;
                            
                            
                            fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_D_poly2_p2");

                            printf("\n\n/////////////////////////////////////// FV_from_H0 D poly2 p2k2 ///////////////////////\n");
                            printf("Polynomial\n");
                            fit_info.Npar=8;
                            fit_info.N=2;
                            fit_info.function=FA_FV_D_poly2_p2k2;
                            
                            
                            fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_D_poly2_p2k2");
                            
                            printf("\n\n/////////////////////////////////////// FV_from_H0 D poly2 M_Pi_D p2k2  ///////////////////////\n");
                            printf("Polynomial\n");
                            fit_info.Npar=12;
                            fit_info.N=2;
                            fit_info.function=FA_FV_D_poly2_M_Pi_D_p2k2;
                            
                            
                            fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_D_poly2_M_Pi_D_p2k2");
    }
    
    printf("\n\n FV_from_H0 for D meson pole\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=FA_FV_D_pole;
                            
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_D_pole");
    
    printf("\n\n FV_from_H0 for D meson pole M_Pi_D\n");
    fit_info.Npar=10;
    fit_info.N=2;
    fit_info.function=FA_FV_D_pole_M_Pi_D;
                            
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_D_pole_M_Pi_D");
    
     printf("\n\n/////////////////////////////////////// FV_from_H0 HA D poly2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=FA_FV_D_poly2;
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_D_poly2");
    printf("\n\n/////////////////////////////////////// FV_from_H0 HA D poly2 p2k2 ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=8;
    fit_info.N=2;
    fit_info.function=FA_FV_D_poly2_p2k2;
    
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_D_poly2_p2k2");
    printf("\n\n/////////////////////////////////////// FV_from_H0 HA D poly3 p2k2 ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=9;
    fit_info.N=2;
    fit_info.function=FA_FV_D_poly3_p2k2;
    
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_D_poly3_p2k2");
    printf("\n\n/////////////////////////////////////// FV_from_H0 HA D poly2 M_Pi_D   ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=10;
    fit_info.N=2;
    fit_info.function=FA_FV_D_poly2_M_Pi_D;
    
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_D_poly2_M_Pi_D");
    
    printf("\n\n/////////////////////////////////////// FV_from_H0 HA D poly2 M_Pi_D p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=12;
    fit_info.N=2;
    fit_info.function=FA_FV_D_poly2_M_Pi_D_p2k2;
    
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_D_poly2_M_Pi_D_p2k2");

    printf("\n\n/////////////////////////////////////// FV_from_H0 HA D poly3 M_Pi_D p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=13;
    fit_info.N=2;
    fit_info.function=FA_FV_D_poly3_M_Pi_D_p2k2;
    
    
    fit_out=fit_FAV_D(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_D_poly3_M_Pi_D_p2k2");


    free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
 
 
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FV_D");
    fit_chi2_good.Nfits=0;
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

    files_declarations(argvNs ,&jack_files,&head);
    read_files_jack(jack_files,head,mass_index,&rephJ);
    grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
       printf("///////////////////////////////////////FA Ds///////////////////////\n///////////////////////////////////////FA Ds///////////////////////\n///////////////////////////////////////FA Ds///////////////////////\n///////////////////////////////////////FA Ds///////////////////////\n///////////////////////////////////////FA Ds///////////////////////\n///////////////////////////////////////FA Ds///////////////////////\n");
   
   
    printf("\n\n FA for Ds meson\n");
    fit_info.Npar=10;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_pole;
    
//    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
//    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_Ds_pole");
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_pole");
    
    printf("\n\n FA for Ds meson pole M_Ds_K_Pi\n");
    fit_info.Npar=8;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_pole_M_Ds_K_Pi;
    
//    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);
//    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_Ds_pole");
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_pole_M_Ds_K_Pi");
    if(2+2==5){
                            printf("\n\n////////////\n HA for Ds pole\n/////////////\n");
                            printf("ChPT pole\n");
                            fit_info.Npar=19;
                            fit_info.N=4;
                            fit_info.function=HA_Ds_pole;
                            fit_info.f1=FA_FV_Ds_pole;
                            fit_info.f2=fDsw0_fit;
                            
                            
                            fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,2 ,fit_info);
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "H","HA_Ds_pole");
                            
                            printf("\n\n FA for Ds meson + xg const\n");
                            fit_info.Npar=11;
                            fit_info.N=4;
                            fit_info.function=FA_FV_Ds_pole_xg_const;

                            fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);  
                            reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"A","FA_Ds_pole_xg_const",head , grephJ );
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_Ds_pole_xg_const");

                            printf("\n\n FA for Ds meson line + xg const\n");
                            fit_info.Npar=7;
                            fit_info.N=4;
                            fit_info.function=FA_FV_Ds_line_xg_const;

                            fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,0 ,fit_info);  
                            reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"A","FA_Ds_line_xg_const",head , grephJ );
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "A","FA_Ds_line_xg_const");
    }
    printf("\n\n/////////////////////////////////////// FA_from_H0 Ds poly2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=10;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_poly2;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_poly2");
  
    printf("\n\n/////////////////////////////////////// FA_from_H0 Ds poly2 p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=12;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_poly2_p2k2;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{A^H}","FA_H_Ds_poly2_p2k2",head , grephJ );
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_poly2_p2k2");
    
    printf("\n\n/////////////////////////////////////// FA_from_H0 Ds poly3 p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=13;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_poly3_p2k2;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{A^H}","FA_H_Ds_poly3_p2k2",head , grephJ );
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_poly3_p2k2");
    if(2+2==5){
                            printf("\n\n/////////////////////////////////////// FA_from_H0 Ds aMpi poly2 p2k2  ///////////////////////\n");
                            printf("Polynomial\n");
                            fit_info.Npar=12;
                            fit_info.N=4;
                            fit_info.function=FA_FV_Ds_poly2_aMpi_p2k2;
                            
                            
                            fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
                            reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{A^H}","FA_H_Ds_poly2_aMpi_p2k2",head , grephJ );
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_poly2_aMpi_p2k2");

                        
                            printf("\n\n/////////////////////////////////////// FA_from_H0 Ds poly2 p2k2 Mpi_overall///////////////////////\n");
                            printf("Polynomial\n");
                            fit_info.Npar=12;
                            fit_info.N=4;
                            fit_info.function=FA_FV_Ds_poly2_p2k2_Mpi_overall;
                            
                            
                            fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
                            print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_poly2_p2k2_Mpi_overall");
    }
    printf("\n\n/////////////////////////////////////// FA_from_H0 Ds poly2 M_Ds_K_Pi ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=8;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_poly2_M_Ds_K_Pi;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{A^H}","FA_H_Ds_poly2_M_Ds_K_Pi",head , grephJ );
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_poly2_M_Ds_K_Pi");
   
    printf("\n\n/////////////////////////////////////// FA_from_H0 Ds poly2 M_Ds_K_Pi_p2k2 ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=10;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_poly2_M_Ds_K_Pi_p2k2;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{A^H}","FA_H_Ds_poly2_M_Ds_K_Pi_p2k2",head , grephJ );
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_poly2_M_Ds_K_Pi_p2k2");

    printf("\n\n/////////////////////////////////////// FA_from_H0 Ds poly3 M_Ds_K_Pi p2k2 ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=11;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_poly3_M_Ds_K_Pi_p2k2;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{A^H}","FA_H_Ds_poly3_M_Ds_K_Pi_p2k2",head , grephJ );
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_poly3_M_Ds_K_Pi_p2k2");
    if(2+2==5){
                        printf("\n\n/////////////////////////////////////// FA_from_H0 Ds poly2 M_Ds_K_Pi p2 ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=9;
                        fit_info.N=4;
                        fit_info.function=FA_FV_Ds_poly2_M_Ds_K_Pi_p2;
                        
                        
                        fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
                        reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{A^H}","FA_H_Ds_poly2_M_Ds_K_Pi_p2",head , grephJ );
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_poly2_M_Ds_K_Pi_p2");
                        
                        printf("\n\n/////////////////////////////////////// FA_from_H0 Ds poly2 M_Ds_K ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=6;
                        fit_info.N=4;
                        fit_info.function=FA_FV_Ds_poly2_M_Ds_K;
                        
                        
                        fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_poly2_M_Ds_K");
                    
                        printf("\n\n/////////////////////////////////////// FA_from_H0 Ds poly2 M_Ds ///////////////////////\n");
                        printf("Polynomial\n");
                        fit_info.Npar=4;
                        fit_info.N=4;
                        fit_info.function=FA_FV_Ds_poly2_M_Ds;
                        
                        
                        fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
                        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Ds_poly2_M_Ds");
                        }
    free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
 
 
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FA_Ds");
    fit_chi2_good.Nfits=0;
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

    files_declarations(argvNs ,&jack_files,&head);
    read_files_jack(jack_files,head,mass_index,&rephJ);
    grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
    
    
    printf("///////////////////////////////////////FV Ds///////////////////////\n//////////////////////////////////////////////////////////////\n//////////////////////////////////////////////////////////////\n//////////////////////////////////////////////////////////////\n//////////////////////////////////////////////////////////////\n//////////////////////////////////////////////////////////////\n");
    if(2+2==5){
                    printf("\n\n FV for Ds meson\n");
                    fit_info.Npar=10;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_pole;
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_Ds_pole");
                
                    printf("\n\n FV_from_H0 for Ds meson\n");
                    fit_info.Npar=10;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_pole;
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Ds_pole");
                
                    
                    printf("\n\n FV for Ds meson + xg const\n");
                    fit_info.Npar=11;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_pole_xg_const;
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"V","FV_Ds_pole_xg_const",head , grephJ );       
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_Ds_pole_xg_const");
                    
                    
                    printf("\n\n/////////////////////////////////////// FV Ds poly2  ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=10;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_poly2;
                    
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_Ds_poly2");
                
                    printf("\n\n/////////////////////////////////////// FV Ds poly2  ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=10;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_poly2;
                    
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,1 ,fit_info);
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "V","FV_Ds_poly2");
                
                    printf("\n\n/////////////////////////////////////// FV_from_H0 Ds poly2  ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=10;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_poly2;
                    
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Ds_poly2");
                    
                    printf("\n\n/////////////////////////////////////// FV_from_H0 Ds poly2 p2k2  ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=12;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_poly2_p2k2;
                    
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{V^H}","FV_H_Ds_poly2_p2k2",head , grephJ );
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Ds_poly2_p2k2");
                    printf("\n\n/////////////////////////////////////// FV_from_H0 Ds poly2 aMpi p2k2  ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=12;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_poly2_aMpi_p2k2;
                    
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{V^H}","FV_H_Ds_poly2_aMpi_p2k2",head , grephJ );
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Ds_poly2_aMpi_p2k2");

                    printf("\n\n/////////////////////////////////////// FV_from_H0 Ds poly2 p2k2 Mpi_overall ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=12;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_poly2_p2k2_Mpi_overall;
                    
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Ds_poly2_p2k2_Mpi_overall");
                
                
                    printf("\n\n/////////////////////////////////////// FV_from_H0 Ds poly2 M_Ds_K_Pi ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=8;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_poly2_M_Ds_K_Pi;
                    
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{V^H}","FV_H_Ds_poly2_M_Ds_K_Pi",head , grephJ );
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Ds_poly2_M_Ds_K_Pi");
                
                    printf("\n\n/////////////////////////////////////// FV_from_H0 Ds poly2 M_Ds_K_Pi p2k2 ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=10;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_poly2_M_Ds_K_Pi_p2k2;
                    
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{V^H}","FV_H_Ds_poly2_M_Ds_K_Pi_p2k2",head , grephJ );
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Ds_poly2_M_Ds_K_Pi_p2k2");
                
                    
                    printf("\n\n/////////////////////////////////////// FV_from_H0 Ds poly3 M_Ds_K_Pi ///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=9;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_poly3_M_Ds_K_Pi;
                    
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{V^H}","FV_H_Ds_poly3_M_Ds_K_Pi",head , grephJ );
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Ds_poly3_M_Ds_K_Pi");
                
                    printf("\n\n/////////////////////////////////////// FV_from_H0 Ds poly2 M_Ds_K_Pi p2///////////////////////\n");
                    printf("Polynomial\n");
                    fit_info.Npar=9;
                    fit_info.N=4;
                    fit_info.function=FA_FV_Ds_poly2_M_Ds_K_Pi_p2;
                    
                    
                    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,4 ,fit_info);
                    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{V^H}","FV_H_Ds_poly2_M_Ds_K_Pi_p2",head , grephJ );
                    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^H}","FV_H_Ds_poly2_M_Ds_K_Pi_p2");
    }
    printf("\n\n/////////////////////////////////////// FV_from_H0 HA Ds pole  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=10;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_pole;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Ds_pole");
    
    printf("\n\n/////////////////////////////////////// FV_from_H0 HA Ds pole M_Ds_K_Pi ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=8;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_pole_M_Ds_K_Pi;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Ds_pole_M_Ds_K_Pi");
    
 
    printf("\n\n/////////////////////////////////////// FV_from_H0 HA Ds poly2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=10;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_poly2;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Ds_poly2");
    printf("\n\n/////////////////////////////////////// FV_from_H0 HA Ds poly2 p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=12;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_poly2_p2k2;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{V^HA}","FV_HA_Ds_poly2_p2k2",head , grephJ );
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Ds_poly2_p2k2");

    printf("\n\n/////////////////////////////////////// FV_from_H0 HA Ds poly3 p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=13;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_poly3_p2k2;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{V^HA}","FV_HA_Ds_poly3_p2k2",head , grephJ );
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Ds_poly3_p2k2");

    printf("\n\n/////////////////////////////////////// FV_from_H0 HA Ds poly2 M_Ds_K_Pi ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=8;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_poly2_M_Ds_K_Pi;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Ds_poly2_M_Ds_K_Pi");
    printf("\n\n/////////////////////////////////////// FV_from_H0 HA Ds poly2 M_Ds_K_Pi p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=10;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_poly2_M_Ds_K_Pi_p2k2;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{V^HA}","FV_HA_Ds_poly2_M_Ds_K_Pi_p2k2",head , grephJ );
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Ds_poly2_M_Ds_K_Pi_p2k2");

    printf("\n\n/////////////////////////////////////// FV_from_H0 HA Ds poly3 M_Ds_K_Pi p2k2  ///////////////////////\n");
    printf("Polynomial\n");
    fit_info.Npar=11;
    fit_info.N=4;
    fit_info.function=FA_FV_Ds_poly3_M_Ds_K_Pi_p2k2;
    
    
    fit_out=fit_FAV_Ds(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
    fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    reweighting_for_plot( argvNs,jack_tot,  fit_out,  fit_info, phys_point,"{V^HA}","FV_HA_Ds_poly3_M_Ds_K_Pi_p2k2",head , grephJ );
    print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Ds_poly3_M_Ds_K_Pi_p2k2");

 
         ///////////////////////free///////////////
    free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
 
 
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FV_Ds");*/

printf("///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    intermopation at ms and mc   bar{MS} 2 GeV /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n");  
  
    fit_chi2_good.Nfits=0;
    fit_FA_FV.Nfits=0;

    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
        
      /*  printf("\n\n///////////////////////////////////////ChPT FA Pi_fpi  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_fpi");
        */
        
        
        printf("\n\n///////////////////////////////////////ChPT FA Pi_fpi   treshold Mpi ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi;
        
        
        
        fit_out=fit_FAV_pion_treshold_Mpi(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,0,1e+5);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_fpi_treshold_Mpi");
    
        printf("\n\n///////////////////////////////////////ChPT FA Pi_fpi p2k2   treshold Mpi ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_p2k2;
        
        
        
        fit_out=fit_FAV_pion_treshold_Mpi(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,0,1e+5);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_fpi_p2k2_treshold_Mpi");
       
       /* printf("\n\n///////////////////////////////////////ChPT FA Pi_fpi_aMpi  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_aMpi;
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_fpi_aMpi");

        printf("\n\n///////////////////////////////////////ChPT FA Pi_fpi_aMpi_M4  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_aMpi_M4;
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_fpi_aMpi_M4");
*/
        printf("\n\n///////////////////////////////////////ChPT FA Pi_fpi_aMpi_M4_M4x  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=9;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_aMpi_M4_M4x;
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_fpi_aMpi_M4_M4x");
        
        
         printf("\n\n///////////////////////////////////////ChPT FA Pi_fpi_aMpi_M4_M4x_p2k2  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=11;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_aMpi_M4_M4x_p2k2;
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_fpi_aMpi_M4_M4x_p2k2");

/*        printf("\n\n///////////////////////////////////////ChPT FA Pi_fpi_aMpi_M4_M4x  fix ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_Pi_fpi_aMpi_M4_M4x_fix;
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_fpi_aMpi_M4_M4x_fix");

        printf("\n\n///////////////////////////////////////ChPT FA Pi_fpi_aMpix  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_aMpix;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_fpi_aMpix");


        printf("\n\n///////////////////////////////////////ChPT FA Pi_w0  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_w0;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_w0");

        printf("\n\n///////////////////////////////////////ChPT FA Pi_fpi_p2k2  ///////////////////////\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_p2k2;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Pi_fpi_p2k2");
  */      
        printf("\n\n FA for Pi meson poly2 only simply A\n");
        fit_info.Npar=3;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_pion_treshold_Mpi(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-0.1,r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Piphys_poly2_simply_A");

        printf("\n\n FA for Pi meson poly2 only simply B\n");
        fit_info.Npar=3;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_pion_treshold_Mpi(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0B-0.1,r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Piphys_poly2_simply_B");
        
        printf("\n\n FA for Pi meson poly2 only simply D\n");
        fit_info.Npar=3;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_pion_treshold_Mpi(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0D-0.1,r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Piphys_poly2_simply_D");

     
    
        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FA_Piphys");
    fit_chi2_good.Nfits=0;
    
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
        
        printf("\n\n///////////////////////////////////////ChPT FV Pi_fpi  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_fpi");
        
        
        printf("\n\n///////////////////////////////////////ChPT FV FA_H_Pi_fpi_treshold_Mpi  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi;
        
        
        
        fit_out=fit_FAV_pion_treshold_Mpi(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,0,1e+5);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_fpi_treshold_Mpi");
        
        printf("\n\n///////////////////////////////////////ChPT FV FA_H_Pi_fpi_p2k2_treshold_Mpi  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_p2k2;
        
        
        
        fit_out=fit_FAV_pion_treshold_Mpi(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,0,1e+5);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_fpi_p2k2_treshold_Mpi");
      
        
        
     /*   
        printf("\n\n///////////////////////////////////////ChPT FV Pi_fpi_aMpi  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_aMpi;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_fpi_aMpi");
       
        printf("\n\n///////////////////////////////////////ChPT FV Pi_fpi_aMpi_M4  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_aMpi_M4;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_fpi_aMpi_M4");
       */
        printf("\n\n///////////////////////////////////////ChPT FV Pi_fpi_aMpi_M4_M4x  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=9;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_aMpi_M4_M4x;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_fpi_aMpi_M4_M4x");
       
        
        
        printf("\n\n///////////////////////////////////////ChPT FV Pi_fpi_aMpi_M4_M4x_p2k2  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=11;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_aMpi_M4_M4x_p2k2;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_fpi_aMpi_M4_M4x_p2k2");
       
       /* 
        printf("\n\n///////////////////////////////////////ChPT FV Pi_fpi_aMpi_M4_M4x fix  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FV_Pi_fpi_aMpi_M4_M4x_fix;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_fpi_aMpi_M4_M4x_fix");
       
        
        printf("\n\n///////////////////////////////////////ChPT FV Pi_fpi_aMpix  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_aMpix;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_fpi_aMpix");
       
        
        printf("\n\n///////////////////////////////////////ChPT FV Pi_w0  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_w0;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_w0");
        
        printf("\n\n///////////////////////////////////////ChPT FV Pi_fpi_p2k2  ///////////////////////\n");
        printf("ChPT pole\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_Pi_fpi_p2k2;
        
        
        
        fit_out=fit_FA_pion_generic(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Pi_fpi_p2k2");
    */
        
               printf("\n\n FV for Pi poly2 simply A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_pion_treshold_Mpi(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0A-0.1, r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Piphys_poly2_simply_A");

        printf("\n\n FV for Pi poly2 simply B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_pion_treshold_Mpi(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0B-0.1, r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Piphys_poly2_simply_B");

        printf("\n\n FV for Pi poly2  simply D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_pion_treshold_Mpi(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0D-0.1, r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Piphys_poly2_simply_D");

    
        
        
        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FV_Piphys");
    fit_chi2_good.Nfits=0;
    compute_FApmFV(argvNs, phys_point,  fit_FA_FV, "Pi","pm");
    compute_FApmFV(argvNs, phys_point,  fit_FA_FV, "Pi","_correlated_");
    fit_FA_FV.Nfits=0;

printf("\n\n///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n K\n/////////////////////////////////////////////////////////////////////////////////////////////\n");
   
     
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
   /*    
        printf("\n\n///////////////////////////////////////ChPT FA K_fK  ///////////////////////\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK;
        
        
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK");
        
        printf("\n\n///////////////////////////////////////ChPT FA FA_K_fK_MK  ///////////////////////\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_MK;
        
        
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_MK");
    
        printf("\n\n///////////////////////////////////////ChPT FA FA_K_fK_MK_M4  ///////////////////////\n");
        fit_info.Npar=9;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_MK_M4;
        
        
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_MK_M4");
       
        printf("\n\n///////////////////////////////////////ChPT FA FA_K_fK_p2k2  ///////////////////////\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_p2k2;
        
        
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);        
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_p2k2");
        
        
        printf("\n\n///////////////////////////////////////ChPT FA FA_K_fK_MK_M4_p2k2  ///////////////////////\n");
        fit_info.Npar=11;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_MK_M4_p2k2;
        
        
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_MK_M4_p2k2");
   
    
        printf("\n\n///////////////////////////////////////ChPT FA FA_K_fK_MK_p2k2  ///////////////////////\n");
        fit_info.Npar=9;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_MK_p2k2;
        
        
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_MK_p2k2");
    */
    
       /* printf("\n\n///////////////////////////////////////ChPT FA K_fK_Mpi  ///////////////////////\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi;
        
        
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_Mpi");
      
        printf("\n\n///////////////////////////////////////ChPT FA FA_K_fK_Mpi_Mpi4  ///////////////////////\n");
        fit_info.Npar=9;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_Mpi4;
        
        
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_Mpi_Mpi4");
       
        
        
         printf("\n\n///////////////////////////////////////ChPT FA K_fK_Mpi_p2k2  ///////////////////////\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_p2k2;
        
        
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_Mpi_p2k2");
      
        printf("\n\n///////////////////////////////////////ChPT FA FA_K_fK_Mpi_Mpi4_p2k2  ///////////////////////\n");
        fit_info.Npar=11;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_Mpi4_p2k2;
        
        
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_Mpi_Mpi4_p2k2");
       */
        printf("\n\n///////////////////////////////////////ChPT FA K_fK_Mpi_heavy  ///////////////////////\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_heavy;
        
        
        
        fit_out=fit_FAV_Kphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,0,1e+5);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_Mpi_heavy");
      
        printf("\n\n///////////////////////////////////////ChPT FA FA_K_fK_Mpi_Mpi4_heavy  ///////////////////////\n");
        fit_info.Npar=10;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_Mpi4_heavy;
        
        
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_Mpi_Mpi4_heavy");
       
        
        
         printf("\n\n///////////////////////////////////////ChPT FA K_fK_Mpi_heavy_p2k2  ///////////////////////\n");
        fit_info.Npar=8;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_heavy_p2k2;
        
        
        
        fit_out=fit_FAV_Kphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,0,1e+5);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_Mpi_heavy_p2k2");
      
        printf("\n\n///////////////////////////////////////ChPT FA FA_K_fK_Mpi_Mpi4_heavy_p2k2  ///////////////////////\n");
        fit_info.Npar=12;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_Mpi4_heavy_p2k2;
        
        
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_K_fK_Mpi_Mpi4_heavy_p2k2");
       
        
        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FA_Kphys");
    fit_chi2_good.Nfits=0;
    
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
    /*    
        printf("\n\n///////////////////////////////////////ChPT FV K_fK  ///////////////////////\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK;
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);        
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK");
        printf("\n\n///////////////////////////////////////ChPT FV K_fK_MK  ///////////////////////\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_MK;
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_MK");
        
        printf("\n\n///////////////////////////////////////ChPT FV K_fK_MK_M4  ///////////////////////\n");
        fit_info.Npar=9;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_MK_M4;
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_MK_M4");
        
        
        printf("\n\n///////////////////////////////////////ChPT FV K_fK_p2k2  ///////////////////////\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_p2k2;
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_p2k2");
        
         printf("\n\n///////////////////////////////////////ChPT FV K_fK_MK_M4_p2k2  ///////////////////////\n");
        fit_info.Npar=11;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_MK_M4_p2k2;
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_MK_M4_p2k2");
        
        
        printf("\n\n///////////////////////////////////////ChPT FV K_fK_MK_p2k2  ///////////////////////\n");
        fit_info.Npar=9;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_MK_p2k2;
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_MK_p2k2");
        */
    /*
         printf("\n\n///////////////////////////////////////ChPT FV K_fK_Mpi  ///////////////////////\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi;
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_Mpi");
    
        printf("\n\n///////////////////////////////////////ChPT FV K_fK_Mpi_Mpi4  ///////////////////////\n");
        fit_info.Npar=9;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_Mpi4;
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_Mpi_Mpi4");
        
        printf("\n\n///////////////////////////////////////ChPT FV K_fK_Mpi_p2k2  ///////////////////////\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_p2k2;
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_Mpi_p2k2");
    
        printf("\n\n///////////////////////////////////////ChPT FV K_fK_Mpi_Mpi4  ///////////////////////\n");
        fit_info.Npar=11;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_Mpi4_p2k2;
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_Mpi_Mpi4_p2k2");
        
        
        */
        printf("\n\n///////////////////////////////////////ChPT FV K_fK_Mpi_heavy  ///////////////////////\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_heavy;
        
        fit_out=fit_FAV_Kphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,0,1e+5);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_Mpi_heavy");
    
        printf("\n\n///////////////////////////////////////ChPT FV K_fK_Mpi_Mpi4_heavy  ///////////////////////\n");
        fit_info.Npar=10;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_Mpi4_heavy;
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_Mpi_Mpi4_heavy");
        
        printf("\n\n///////////////////////////////////////ChPT FV K_fK_Mpi_heavy_p2k2  ///////////////////////\n");
        fit_info.Npar=8;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_heavy_p2k2;
        
        fit_out=fit_FAV_Kphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,0,1e+5);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_Mpi_heavy_p2k2");
    
        printf("\n\n///////////////////////////////////////ChPT FV K_fK_Mpi_Mpi4_heavy  ///////////////////////\n");
        fit_info.Npar=12;
        fit_info.N=1;
        fit_info.function=FA_FV_K_fK_Mpi_Mpi4_heavy_p2k2;
        
        fit_out=fit_FAV_Kphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_K_fK_Mpi_Mpi4_heavy_p2k2");
        
    
        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FV_Kphys");
    fit_chi2_good.Nfits=0;
    compute_FApmFV(argvNs, phys_point,  fit_FA_FV, "K","pm");
    compute_FApmFV(argvNs, phys_point,  fit_FA_FV, "K","_correlated_");
    fit_FA_FV.Nfits=0;
   

printf("\n\n///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n D \n/////////////////////////////////////////////////////////////////////////////////////////////\n");
       
    
   
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
/*
        printf("\n\n FA for D meson poly2\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2");
*/

/*      
 *      printf("\n\n FA for D meson poly2 only A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-0.1,r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2_A");

        printf("\n\n FA for D meson poly2 only B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0B-0.1,r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2_B");
        
        printf("\n\n FA for D meson poly2 only D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0D-0.1,r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2_D");
        
        */
        
        printf("\n\n FA for D meson poly2 only simply A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-0.1,r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2_simply_A");

        printf("\n\n FA for D meson poly2 only simply B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0B-0.1,r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2_simply_B");
        
        printf("\n\n FA for D meson poly2 only simply D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0D-0.1,r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2_simply_D");
        
        printf("\n\n FA for D meson poly2  simply Mx\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_Mx;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2_simply_Mx");
        
       /*
        printf("\n\n FA for D meson poly2  simply ax\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2_simply_ax");
        */
        
        printf("\n\n FA for D meson poly2  simply ax Mx\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax_Mx;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2_simply_ax_Mx");
      
        
        printf("\n\n FA for D meson pole only simply A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_fixa_simply;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-0.1,r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_pole_simply_A");

        printf("\n\n FA for D meson pole only simply B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_fixa_simply;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0B-0.1,r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_pole_simply_B");
        
        printf("\n\n FA for D meson poly2 only simply D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_fixa_simply;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0D-0.1,r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_pole_simply_D");
 
        printf("\n\n FA for D meson pole simply Mx\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_Mx;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-10,r0A+10);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_pole_simply_Mx");

        printf("\n\n FA for D meson pole simply ax Mx\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_ax_Mx;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-10,r0A+10);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_pole_simply_ax_Mx");

        /*
        
        printf("\n\n FA for D meson poly2_p2k2\n");
        fit_info.Npar=8;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_p2k2;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2_p2k2");
        
        printf("\n\n FA for D meson pole\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_pole");
        
        printf("\n\n FA for D meson pole_p2k2\n");
        fit_info.Npar=8;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_p2k2;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_pole_p2k2");
       */
        
        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FA_Dphys");
    fit_chi2_good.Nfits=0;
    

    
    
   
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);

        printf("\n\n FA for D meson poly2  simply Mx\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_Mx;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2_simply_Mx");
        
        
        printf("\n\n FA for D meson poly2  simply ax Mx\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax_Mx;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_poly2_simply_ax_Mx");
      
    
        
        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FA_Dphys_poly");
    fit_chi2_good.Nfits=0;

    
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);

        
        printf("\n\n FA for D meson pole simply Mx\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_Mx;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-10,r0A+10);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV_pole=save_fit(fit_FA_FV_pole,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_pole_simply_Mx");

        printf("\n\n FA for D meson pole simply ax Mx\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_ax_Mx;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-10,r0A+10);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV_pole=save_fit(fit_FA_FV_pole,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dphys_pole_simply_ax_Mx");

        
        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FA_Dphys_pole");
    fit_chi2_good.Nfits=0;
    
    
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
/*
        printf("\n\n FV for D poly2\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2");

        printf("\n\n FV for D poly2 A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0A-0.1, r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2_A");

        printf("\n\n FV for D poly2 B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0B-0.1, r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2_B");

        printf("\n\n FV for D poly2 D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0D-0.1, r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2_D");
  */       
         printf("\n\n FV for D poly2 simply A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0A-0.1, r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2_simply_A");

        printf("\n\n FV for D poly2 simply B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0B-0.1, r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2_simply_B");

        printf("\n\n FV for D poly2  simply D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0D-0.1, r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2_simply_D");
         
         printf("\n\n FV for D pole simply A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_fixa_simply;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0A-0.1, r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_pole_simply_A");

        printf("\n\n FV for D pole simply B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_fixa_simply;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0B-0.1, r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_pole_simply_B");

        printf("\n\n FV for D pole  simply D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_fixa_simply;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0D-0.1, r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_pole_simply_D");
       
         printf("\n\n FV for D pole simply Mx \n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_Mx;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0A-10.1, r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_pole_simply_Mx");

         printf("\n\n FV for D pole simply ax Mx \n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_ax_Mx;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0A-10.1, r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_pole_simply_ax_Mx");

        
        printf("\n\n FV for D poly2 simply Mx \n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_Mx;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2_simply_Mx");
       
  /*      
        printf("\n\n FV for D poly2 simply ax\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2_simply_ax");
    */  
        printf("\n\n FV for D poly2 simply ax Mx\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax_Mx;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2_simply_ax_Mx");
      
        
      
/*        
        printf("\n\n FV for D poly2_p2k2\n");
        fit_info.Npar=8;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_p2k2;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2_p2k2");
        
        printf("\n\n FV for D pole\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_pole");
        
        printf("\n\n FV for D pole_p2k2\n");
        fit_info.Npar=8;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_p2k2;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_pole_p2k2");
        */
        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FV_Dphys");
    fit_chi2_good.Nfits=0;
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);

        
        printf("\n\n FV for D poly2 simply Mx \n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_Mx;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2_simply_Mx");
       
 
        printf("\n\n FV for D poly2 simply ax Mx\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax_Mx;
        
        fit_out=fit_FAV_Dphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_poly2_simply_ax_Mx");
      
        
        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FV_Dphys_poly");
    fit_chi2_good.Nfits=0;
    compute_FApmFV(argvNs, phys_point,  fit_FA_FV, "D_poly","pm");
    compute_FApmFV(argvNs, phys_point,  fit_FA_FV, "D_poly","_correlated_");
    fit_FA_FV.Nfits=0;
    
     for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);

         printf("\n\n FV for D pole simply Mx \n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_Mx;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0A-10.1, r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV_pole=save_fit(fit_FA_FV_pole,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_pole_simply_Mx");

         printf("\n\n FV for D pole simply ax Mx \n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_ax_Mx;
        
        fit_out=fit_FAV_Dphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info, r0A-10.1, r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV_pole=save_fit(fit_FA_FV_pole,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dphys_pole_simply_ax_Mx");

        
        
        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FV_Dphys_pole");
    fit_chi2_good.Nfits=0;
    compute_FApmFV(argvNs, phys_point,  fit_FA_FV_pole, "D_pole","_correlated_pole_");
    fit_FA_FV_pole.Nfits=0;
     
  
printf("\n\n///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n Ds \n/////////////////////////////////////////////////////////////////////////////////////////////\n");
       
    
    
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
/*
        printf("\n\n FA for Ds meson poly2\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dsphys_poly2;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2");
        
        
        printf("\n\n FA for Ds meson poly2 only A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-0.1,r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_A");
        
        printf("\n\n FA for Ds meson poly2 only B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0B-0.1,r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_B");
        
        printf("\n\n FA for Ds meson poly2 only D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0D-0.1,r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_D");
 */     
         printf("\n\n FA for Ds meson poly2 only simply A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply ;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-0.1,r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_simply_A");
        
        printf("\n\n FA for Ds meson poly2 only simply  B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply ;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0B-0.1,r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_simply_B");
        
        printf("\n\n FA for Ds meson poly2 only simply D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0D-0.1,r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_simply_D");
      
        printf("\n\n FA for Ds meson poly2  simply Mx \n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_Mx ;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_simply_Mx");
       /* 
        printf("\n\n FA for Ds meson poly2  simply  ax\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax ;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_simply_ax");
        */
        printf("\n\n FA for Ds meson poly2  simply  ax M4\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax_Mx ;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_simply_ax_Mx");
      
        
      /*  
        printf("\n\n FA for Ds meson poly2  simply  ax a4x\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax_a4x ;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_simply_ax_a4x");
        */
        
         printf("\n\n FA for Ds meson pole only simply A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_fixa_simply ;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-0.1,r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_pole_simply_A");
        
        printf("\n\n FA for Ds meson pole only simply  B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_fixa_simply ;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0B-0.1,r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_pole_simply_B");
        
        printf("\n\n FA for Ds meson pole only simply D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_fixa_simply;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0D-0.1,r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_pole_simply_D");
      
         printf("\n\n FA for Ds meson pole  simply Mx\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_Mx ;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_pole_simply_Mx");
        /*
         printf("\n\n FA for Ds meson pole  simply ax \n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_ax ;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_pole_simply_ax");
        */
         printf("\n\n FA for Ds meson pole  simply ax Mx \n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_ax_Mx ;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_pole_simply_ax_Mx");
     /*   
        printf("\n\n FA for Ds meson pole  simply  apole\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_apole ;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_pole_simply_apole");
        */
/*             
        printf("\n\n FA for Ds meson poly2_p2k2\n");
        fit_info.Npar=8;
        fit_info.N=1;
        fit_info.function=FA_FV_Dsphys_poly2_p2k2;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_p2k2");
        
        printf("\n\n FA for Ds meson pole\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dsphys_pole;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_pole");
        
        printf("\n\n FA for Ds meson pole_p2k2\n");
        fit_info.Npar=8;
        fit_info.N=1;
        fit_info.function=FA_FV_Dsphys_pole_p2k2;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_pole_p2k2");
     */  
        
        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FA_Dsphys");
    fit_chi2_good.Nfits=0;
    
       
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);

      
        printf("\n\n FA for Ds meson poly2  simply Mx \n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_Mx ;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_simply_Mx");
     
        printf("\n\n FA for Ds meson poly2  simply  ax M4\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax_Mx ;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,3 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_poly2_simply_ax_Mx");
      

        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FA_Dsphys_poly");
    fit_chi2_good.Nfits=0;
   
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);

      
        printf("\n\n FA for Ds meson pole  simply Mx\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_Mx ;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV_pole=save_fit(fit_FA_FV_pole,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_pole_simply_Mx");

        printf("\n\n FA for Ds meson pole  simply ax Mx \n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_ax_Mx ;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,3 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV_pole=save_fit(fit_FA_FV_pole,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{A^H}","FA_H_Dsphys_pole_simply_ax_Mx");

        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FA_Dsphys_pole");
    fit_chi2_good.Nfits=0;
    
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
/*
        printf("\n\n FV for Ds poly2\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dsphys_poly2;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly2");
        
        
        printf("\n\n FV for Ds meson poly2 only A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0A-0.1,r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly2_A");
        
        printf("\n\n FV for Ds meson poly2 only B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0B-0.1,r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly2_B");
        
        printf("\n\n FV for Ds meson poly2 only D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0D-0.1,r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly2_D");
*/
        printf("\n\n FV for Ds meson poly2 only simply A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0A-0.1,r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly2_simply_A");
        
        printf("\n\n FV for Ds meson poly2 only simply B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0B-0.1,r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly2_simply_B");
        
        printf("\n\n FV for Ds meson poly2 only simply D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_fixa_simply;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0D-0.1,r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly2_simply_D");
         
        printf("\n\n FV for Ds meson poly2  simply  Mx\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_Mx;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly2_simply_Mx");
      
        /*
        printf("\n\n FV for Ds meson poly2  simply ax \n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly2_simply_ax");
      */
        
        printf("\n\n FV for Ds meson poly2  simply ax Mx\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax_Mx;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly2_simply_ax_Mx");
      
     /*   
        printf("\n\n FV for Ds meson poly2  simply ax a4x\n");
        fit_info.Npar=7;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax_a4x;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly2_simply_ax_a4x");
      */
     
        printf("\n\n FV for Ds meson pole only simply A\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_fixa_simply;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0A-0.1,r0A+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_pole_simply_A");
        
        printf("\n\n FV for Ds meson pole only simply B\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_fixa_simply;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0B-0.1,r0B+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_pole_simply_B");
        
        printf("\n\n FV for Ds meson pole only simply D\n");
        fit_info.Npar=4;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_fixa_simply;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0D-0.1,r0D+0.1);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_pole_simply_D");
   
      
        printf("\n\n FV for Ds meson pole  simply Mx \n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_Mx;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_pole_simply_Mx");
       /* 
        printf("\n\n FV for Ds meson pole  simply ax\n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_ax;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_pole_simply_ax");
        */
        printf("\n\n FV for Ds meson pole  simply ax Mx\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_ax_Mx;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_pole_simply_ax_Mx");
        
 /*        
        printf("\n\n FV for Ds meson pole  simply apole \n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_apole;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_pole_simply_apole");
        
       
        printf("\n\n FV for Ds poly2_p2k2\n");
        fit_info.Npar=8;
        fit_info.N=1;
        fit_info.function=FA_FV_Dsphys_poly2_p2k2;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly2_p2k2");
        
        printf("\n\n FV for Ds pole\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dsphys_pole;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_pole");
        
        printf("\n\n FV for Ds pole_p2k2\n");
        fit_info.Npar=8;
        fit_info.N=1;
        fit_info.function=FA_FV_Dsphys_pole_p2k2;
        
        fit_out=fit_FAV_Dsphys(jack_files,  head ,jack_tot, grephJ,5 ,fit_info);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_pole_p2k2");
 */       
        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FV_Dsphys");
    fit_chi2_good.Nfits=0;
     
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
   
        printf("\n\n FV for Ds meson poly  simply Mx \n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_Mx;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly_simply_Mx");

        printf("\n\n FV for Ds meson poly  simply ax Mx\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_poly2_simply_ax_Mx;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV=save_fit(fit_FA_FV,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_poly_simply_ax_Mx");
        

        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FV_Dsphys_poly");
    fit_chi2_good.Nfits=0;
    compute_FApmFV(argvNs, phys_point,  fit_FA_FV, "Ds_poly","pm");
    compute_FApmFV(argvNs, phys_point,  fit_FA_FV, "Ds_poly","_correlated_");
    fit_FA_FV.Nfits=0;
    
    
    for(Ns=0;Ns<Nsets;Ns++){
    
        argvNs[0]=argv[0];
        argvNs[1]=argv[1];
        argvNs[2]=argv[2+Ns*3];
        argvNs[3]=argv[3+Ns*3];
        argvNs[4]=argv[4+Ns*3];

        files_declarations(argvNs ,&jack_files,&head);
        read_files_jack(jack_files,head,mass_index,&rephJ);
        grephJ=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, &rephJ);
   
        printf("\n\n FV for Ds meson pole  simply Mx \n");
        fit_info.Npar=5;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_Mx;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV_pole=save_fit(fit_FA_FV_pole,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_pole_simply_Mx");

        printf("\n\n FV for Ds meson pole  simply ax Mx\n");
        fit_info.Npar=6;
        fit_info.N=1;
        fit_info.function=FA_FV_Dphys_pole_simply_ax_Mx;
        
        fit_out=fit_FAV_Dsphys_treshold(jack_files,  head ,jack_tot, grephJ,5 ,fit_info,r0A-10.1,r0A+10.1);
        fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
        fit_FA_FV_pole=save_fit(fit_FA_FV_pole,fit_info,fit_out);
        print_fit_info( argvNs,jack_tot,  fit_out,  fit_info, phys_point, grephJ, head, "{V^HA}","FV_HA_Dsphys_pole_simply_ax_Mx");
        

        free_data( &jack_files, &head,&jack_tot, &grephJ);
    }
    
    compute_systematics(argvNs, phys_point,  fit_chi2_good, "FV_Dsphys_pole");
    fit_chi2_good.Nfits=0;
    compute_FApmFV(argvNs, phys_point,  fit_FA_FV_pole, "Ds_pole","_correlated_pole_");
    fit_FA_FV_pole.Nfits=0;
    
    free_tif(jack_tot,phys_point);
    free_results();

    return 0;
    
}
