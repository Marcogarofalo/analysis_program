#define pion_FVE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"

#include "resampling.hpp"
#include "read.hpp"
//#include "m_eff.hpp"
#include "gnuplot.hpp"
//#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "rand.hpp"
#include "non_linear_fit.hpp"
#include "pion.hpp"
#include "fve.hpp"
#include "mutils.hpp"

#include <unistd.h>

#include <omp.h>

double flagl3b=3.53;
double flagl4b=4.73;


double Mw2_chiral_FVE_a0_minus(double  input, double x,int Npar,double  *P){
    
    int  n=0;
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b;
    if (Npar==6)
         Bw=P[0], fw=P[1], l3b=P[2], P2=P[3], l4b=P[4];
    else if (Npar==4)//flag 2019
         Bw=P[0], fw=P[1], l3b=flagl3b, P2=P[2], l4b=flagl4b;
    double mw=x;
    double KM=1,Kf=1;
    
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
   
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    if (n==0){
        
        Mw2=1+xi*log(xi)+P1*xi;
        Mw2*=2*Bw*mw;
        
        Mw2-=input;
        
    }
    if (n==1){
        error(1==1,3,"Mw2_chiral_FVE_a0_minus", "the function was called with n=1, in this way it returns f_PS");
        Mw2=1-2*xi*log(xi)+P3*xi;
        Mw2*=fw;
        
    }
     return Mw2;
}

double Mw2_over_fw2_chiral_FVE_a0_minus(double  input, double x,int Npar,double  *P){
    
    int  n=0;
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b;
    if (Npar>=6 )
         Bw=P[0], fw=P[1], l3b=P[2], P2=P[3], l4b=P[4];
    else if (Npar==4)//flag 2019
         Bw=P[0], fw=P[1], l3b=flagl3b, P2=P[2], l4b=flagl4b;
    double mw=x;
    double KM=1,Kf=1;
    
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
   
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    double fw2;
        
        Mw2=1+xi*log(xi)+P1*xi;
        Mw2*=2*Bw*mw;
        
        //Mw2-=input;
        
        fw2=(1-2*xi*log(xi)+P3*xi)*fw;
        fw2*=fw2;
        
        Mw2=Mw2/fw2  -input;
        
     return Mw2;
}


double Mw2_chiral_FVE_a0_minus_P1_P3(double  input, double x,int Npar,double  *P){
    
    int  n=0;
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b, P1,P3;
    if (Npar==6)
         Bw=P[0], fw=P[1], P1=P[2], P2=P[3], P3=P[4];
    else if (Npar==4)//flag 2019
         Bw=P[0], fw=P[1], l3b=flagl3b, P2=P[2], l4b=flagl4b;
    double mw=x;
    double KM=1,Kf=1;
    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    if (n==0){
        
        Mw2=1+xi*log(xi)+P1*xi;
        Mw2*=2*Bw*mw;
        
        Mw2-=input;
        
    }
    if (n==1){
        error(1==1,3,"Mw2_chiral_FVE_a0_minus", "the function was called with n=1, in this way it returns f_PS");
        Mw2=1-2*xi*log(xi)+P3*xi;
        Mw2*=fw;
        
    }
     return Mw2;
}

double fPSw_chiral_FVE(  double x,int Npar,double  *P){
    
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b;
    if (Npar>=6)
         Bw=P[0], fw=P[1], l3b=P[2], P2=P[3], l4b=P[4];
    else if (Npar==4)//flag 2019
         Bw=P[0], fw=P[1], l3b=flagl3b, P2=P[2], l4b=flagl4b;
    double mw=x;
    double KM=1,Kf=1;
    
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    Mw2=1-2*xi*log(xi)+P3*xi;
    Mw2*=fw;
        
    return Mw2;
}
double fPSw_chiral_FVE_P1_P3(  double x,int Npar,double  *P){
    
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b,P1,P3;
    if (Npar==6)
         Bw=P[0], fw=P[1], P1=P[2], P2=P[3], P3=P[4];
    else if (Npar==4)//flag 2019
         Bw=P[0], fw=P[1], l3b=flagl3b, P2=P[2], l4b=flagl4b;
    double mw=x;
    double KM=1,Kf=1;
    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    Mw2=1-2*xi*log(xi)+P3*xi;
    Mw2*=fw;
        
    return Mw2;
}

double Mw2_fw_chiral_FVE(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b, P4,ZP;
    if (Npar==6)
         Bw=P[0], fw=P[1], l3b=P[2], P2=P[3], l4b=P[4], P4=P[5];
    else if (Npar==4)//flag 2019
         Bw=P[0], fw=P[1], l3b=flagl3b, P2=P[2], l4b=flagl4b, P4=P[3];
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
    int Lsize=(int(x[4]));
    double KM,Kf;
    
   
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    
    FVE( v_w0GeV, w0,l3b, l4b, Bw, fw, Lsize, mw ,  dmpi2,  dfpi, &KM, &Kf);
    
    double M,L,Rf,RM;
    M=sqrt(dmpi2)*w0/v_w0MeV;
    L=Lsize*v_w0fm/w0;
    Rf=exp(4.58982-0.0138032*M-2.013*L);
    RM=exp(3.81729-0.0130342*M-2.1714*L);
    Kf=1-Rf;
    KM=RM+1;
    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    if (n==0){
        
        Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2;
        Mw2*=2*Bw*mw*KM*KM;
        
    }
    if (n==1){
        
        Mw2=1-2*xi*log(xi)+P3*xi+ (1./(w0*w0))*P4;
        Mw2*=fw*Kf;
        
    }
     return Mw2;
    
}


double Mw2_fw_chiral_FVE_P1_P3(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b, P4,ZP,P1,P3;
    if (Npar==6)
         Bw=P[0], fw=P[1], P1=P[2], P2=P[3], P3=P[4], P4=P[5];
    else if (Npar==4)//flag 2019
         Bw=P[0], fw=P[1], l3b=flagl3b, P2=P[2], l4b=flagl4b, P4=P[3], ZP=P[4];
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
    int Lsize=(int(x[4]));
    double KM,Kf;
    
   
    //double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    //double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    
    //FVE( v_w0GeV, w0,l3b, l4b, Bw, fw, Lsize, mw ,  dmpi2,  dfpi, &KM, &Kf);
    
    double M,L,Rf,RM;
    M=sqrt(dmpi2)*w0/v_w0MeV;
    L=Lsize*v_w0fm/w0;
    Rf=exp(4.58982-0.0138032*M-2.013*L);
    RM=exp(3.81729-0.0130342*M-2.1714*L);
    Kf=1-Rf;
    KM=RM+1;

    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    if (n==0){
        
        Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2;
        Mw2*=2*Bw*mw*KM*KM;
        
    }
    if (n==1){
        
        Mw2=1-2*xi*log(xi)+P3*xi+ (1./(w0*w0))*P4;
        Mw2*=fw*Kf;
        
    }
    
     return Mw2;
    
}

double Mw2_fw_chiral_FVE_prior(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b, P4,ZPA,ZPB;
    if (Npar==8)
         Bw=P[0], fw=P[1], l3b=P[2], P2=P[3], l4b=P[4], P4=P[5], ZPA=P[6], ZPB=P[7];
    else if (Npar==7)//flag 2019
         Bw=P[0], fw=P[1], l3b=P[2], P2=P[3], l4b=P[4], P4=0, ZPA=P[5], ZPB=P[6];
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
    int Lsize=(int(x[4]));
    double KM,Kf;
    
   
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    if(n==0 || n==2 )
        mw=mw/ZPA;
    if(n==1 || n==3)
        mw=mw/ZPB;
    FVE( v_w0GeV,w0,l3b, l4b, Bw, fw, Lsize, mw ,  dmpi2,  dfpi, &KM, &Kf);

    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    if (n==0 ||n==1){
        Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2;
        Mw2*=2*Bw*mw*KM*KM;
    }
    if (n==2 || n==3){
        Mw2=1-2*xi*log(xi)+P3*xi+ (1./(w0*w0))*P4;
        Mw2*=fw*Kf;
    }
    if (n==4 ){
        Mw2=ZPA;
    }
    if (n==5 ){
        Mw2=ZPB;//printf("ZPB=%f\n",ZPB);
    }
    return Mw2;
    
}

double Mw2_fw_chiral_FVE_P40(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b, P4=0;
    if (Npar==5)
         Bw=P[0], fw=P[1], l3b=P[2], P2=P[3], l4b=P[4] ;
    else if (Npar==3)//flag 2019
         Bw=P[0], fw=P[1], l3b=flagl3b, P2=P[2], l4b=flagl4b ;
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
    int Lsize=(int(x[4]));
    double KM,Kf;
    
   
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    
    FVE( v_w0GeV, w0,l3b, l4b, Bw, fw, Lsize, mw ,  dmpi2,  dfpi, &KM, &Kf);

    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    if (n==0){
        
        Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2;
        Mw2*=2*Bw*mw*KM*KM;
        
    }
    if (n==1){
        
        Mw2=1-2*xi*log(xi)+P3*xi+ (1./(w0*w0))*P4;
        Mw2*=fw*Kf;
        
    }
     return Mw2;
}


double fw_of_Mw2_physical_point(double x,int Npar,double  *P){
    
    double r,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b, P4,c3;
    if (Npar==4)
         fw=P[0], l3b=P[1], l4b=P[2], P4=P[3];
    else if (Npar==3)//flag 2019
        fw=P[0], l3b=flagl3b, l4b=P[1];
    double  Mw2=x;
    double KM,Kf;
    
   
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    
    //FVE( v_w0GeV, w0,l3b, l4b, 0.5, fw, Lsize, dmpi2*w0*w0 ,  dmpi2,  dfpi, &KM, &Kf);

    xi=(Mw2)/(16.*pi*pi*fw*fw);
   
        
    r=1-2*xi*log(xi)+P3*xi;
    r*=fw;
        
    return r;
    
}

double fw_of_Mw2_chiral_FVE(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
//    double Bw, fw, l3b, P2, l4b, P4;
    double    fw=P[0], l3b=flagl3b, l4b=P[1], P4=P[2], c2w4=0, c3=0;
    if (Npar>3) c2w4=P[3];
    double  w0=x[0], dmpi2=x[1], dfpi=x[2];
    int Lsize=(int(x[3]));
    double KM,Kf;
    
   
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    
    FVE( v_w0GeV, w0,l3b, l4b, 0.5, fw, Lsize, dmpi2*w0*w0 ,  dmpi2,  dfpi, &KM, &Kf);
   
    double M,L,Rf,RM;
    M=sqrt(dmpi2)*w0/v_w0MeV;
    L=Lsize*v_w0fm/w0;
    Rf=exp(4.58982-0.0138032*M-2.013*L);
    RM=exp(3.81729-0.0130342*M-2.1714*L);
    Kf=1-Rf;
    KM=RM+1;
//  //FVE_K(2.550*w0GeV ,0.1214*w0GeV ,6.1185*2.282/w0GeV/*6.1185*/, 0.030572*w0GeV/2.282/*mlw*/, 0.23441*w0GeV/2.282 /*msw*/ ,0.41728  , 0.34120,  1.75986  , 0.40402, &KM, &Kf);
//printf("K:  KM=%.10f    Kf=%.10f\n", KM,Kf);

    xi=(dmpi2*w0*w0/(KM*KM))/(16.*pi*pi*fw*fw);
    xi=(dmpi2*w0*w0)/(16.*pi*pi*fw*fw);
   
        
//    Mw2=1-2*xi*log(xi)+P3*xi+ (1./(w0*w0))*(P4-4.*c2w4*log(xi)/pow(4*pi*fw,4))+c3*xi*xi;
    Mw2=1-2*xi*log(xi)+P3*xi+ (1./(w0*w0))*(P4);

    Mw2*=fw*Kf;
        
    return Mw2;
}

double fw_over_Mw_of_Mw_minus_exp(double  input, double x,int Npar,double  *P){//(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi,f_m;
    double pi=3.141592653589793;
//    double Bw, fw, l3b, P2, l4b, P4;
    double    fw=P[0], l3b=flagl3b, l4b=P[1], P4=P[2], c2w4=0, c3=0;
    if (Npar>3) c2w4=P[3];
    //double  w0=x[0], dmpi2=x[1], dfpi=x[2];
    double Mw=x;
    double KM=1.,Kf=1.;
    
   
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
        
    xi=(Mw*Mw)/(16.*pi*pi*fw*fw);
   
        
//    Mw2=1-2*xi*log(xi)+P3*xi+ (1./(w0*w0))*(P4-4.*c2w4*log(xi)/pow(4*pi*fw,4))+c3*xi*xi;
    f_m=1-2*xi*log(xi)+P3*xi;
    f_m*=fw/Mw;
    f_m-=input;
    
    return f_m;
}

double fw_of_Mw2_chiral_FVE_P40(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b, P4;
    
        fw=P[0], l3b=flagl3b, l4b=P[1], P4=0;
    double  w0=x[0], dmpi2=x[1], dfpi=x[2];
    int Lsize=(int(x[3]));
    double KM,Kf;
    
   
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    
    FVE( v_w0GeV, w0,l3b, l4b, 0.5, fw, Lsize, dmpi2*w0*w0 ,  dmpi2,  dfpi, &KM, &Kf);

    xi=(dmpi2*w0*w0/(KM*KM))/(16.*pi*pi*fw*fw);
   
        
    Mw2=1-2*xi*log(xi)+P3*xi+ (1./(w0*w0))*P4;
    Mw2*=fw*Kf;
        
    return Mw2;
    
}

double **fit_Mpi_fw_chiral_FVE_P1_P3(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit; 
   int i,j,e,im;  
   int Npar=6;
   int Nvar=5;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=0;
   int n,count,N=2;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=ensembles;
   en[1]=ensembles;
   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
    guess[0]=2.05478;
    guess[1]=0.113021;
    guess[2]=2.54451;
    guess[3]=0.410103;
    guess[4]=4.61247;
    guess[5]=-0.272884;
   
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double**) malloc(sizeof(double*)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }
   

   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j];
                        rm[j]*=rm[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
                
             //   printf("%g     %g      %g   %g\n",x[e+count][1],x[e+count][0],fit[e+count][0],fit[e+count][1]);
        }
        count+=en[n];
   }

  
   #pragma omp parallel for  private(tmp,i,count,n,e,im,x)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2)
   for (j=0;j<Njack;j++){
        count=0;
        x=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                x[e+count][1]=gJ[e].w0[j];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][j]*gJ[e].M_PS_jack[im][j];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][j];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
            }
            count+=en[n];
        }

        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, Mw2_fw_chiral_FVE_P1_P3,guess );
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, Mw2_fw_chiral_FVE_P1_P3  );
        if (j==Njack-1) printf("chi2/dof (Nb-1)=%g\n",chi2[j]/(en_tot-Npar));
              
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }         
        
        free(tmp);
        for (e=0;e<en_tot;e++)
            free(x[e]);   
        free(x);
   }  
   
   
   
// make_plots_MPi_fPi(en,);
   printf("w0/a[fm]     m*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err    KM2/Kf\n"); 
   double KM,Kf,K;
   count=0;
   x=(double**) malloc(sizeof(double*)*(en_tot));

   for (n=0;n<N;n++){
        printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=gJ[e].w0[Njack-1];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][Njack-1]*gJ[e].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
                
                FVE( v_w0GeV, gJ[e].w0[Njack-1] ,r[2][Njack-1], r[4][Njack-1], r[0][Njack-1], r[1][Njack-1],head[e].l1, x[e+count][0] ,  x[e+count][2], x[e+count][3], &KM, &Kf);
                
                double M,L,Rf,RM;
                M=sqrt(x[e+count][2])*gJ[e].w0[Njack-1]/v_w0MeV;
                L=x[e+count][4]*v_w0fm/gJ[e].w0[Njack-1];
                Rf=exp(4.58982-0.0138032*M-2.013*L);
                RM=exp(3.81729-0.0130342*M-2.1714*L);
                Kf=1-Rf;
                KM=RM+1;
                if (n==0)
                    K=KM*KM;
                else if (n==1)
                    K=Kf;
                printf("%g     %g      %g   %g     %g \n",x[e+count][1],x[e+count][0],fit[e+count][0]/K,fit[e+count][1]/K,K);
                
        }
        count+=en[n];
    }

   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
   printf("$\\chi^2/dof=%f+-%f$\n",chi2m[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
   free(rm);free(chi2m);
   
  
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);   free(fit[e]);
  
   }
   
   free(fit);  

   free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   return r;
    
} 




double **fit_Mpi_fw_chiral_FVE(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit; 
   int i,j,e,im;  
   int Npar=6;
   int Nvar=5;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=0;
   int n,count,N=2;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=ensembles;
   en[1]=ensembles;
   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
    guess[0]=2.05478;
    guess[1]=0.113021;
    guess[2]=2.54451;
    guess[3]=0.410103;
    guess[4]=4.61247;
    guess[5]=-0.272884;
   
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double**) malloc(sizeof(double*)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }
   

   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j];
                        rm[j]*=rm[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
                
             //   printf("%g     %g      %g   %g\n",x[e+count][1],x[e+count][0],fit[e+count][0],fit[e+count][1]);
        }
        count+=en[n];
   }

  
   #pragma omp parallel for  private(tmp,i,count,n,e,im,x)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2)
   for (j=0;j<Njack;j++){
        count=0;
        x=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                x[e+count][1]=gJ[e].w0[j];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][j]*gJ[e].M_PS_jack[im][j];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][j];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
            }
            count+=en[n];
        }

        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, Mw2_fw_chiral_FVE,guess );
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, Mw2_fw_chiral_FVE  );
              
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }         
        
        free(tmp);
        for (e=0;e<en_tot;e++)
            free(x[e]);   
        free(x);
   }  
   
   
   
// make_plots_MPi_fPi(en,);
   printf("w0/a[fm]     m*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err    KM2/Kf\n"); 
   double KM,Kf,K;
   count=0;
   x=(double**) malloc(sizeof(double*)*(en_tot));

   for (n=0;n<N;n++){
        printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=gJ[e].w0[Njack-1];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][Njack-1]*gJ[e].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
                
                FVE( v_w0GeV, gJ[e].w0[Njack-1] ,r[2][Njack-1], r[4][Njack-1], r[0][Njack-1], r[1][Njack-1],head[e].l1, x[e+count][0] ,  x[e+count][2], x[e+count][3], &KM, &Kf);
                
                double M,L,Rf,RM;
                M=sqrt(x[e+count][2])*gJ[e].w0[Njack-1]/v_w0MeV;
                L=x[e+count][4]*v_w0fm/gJ[e].w0[Njack-1];
                Rf=exp(4.58982-0.0138032*M-2.013*L);
                RM=exp(3.81729-0.0130342*M-2.1714*L);
                Kf=1-Rf;
                KM=RM+1;
                if (n==0)
                    K=KM*KM;
                else if (n==1)
                    K=Kf;
                printf("%g     %g      %g   %g     %g \n",x[e+count][1],x[e+count][0],fit[e+count][0]/K,fit[e+count][1]/K,K);
        }
        count+=en[n];
    }

   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
   printf("$\\chi^2/dof=%f+-%f$\n",chi2m[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
   free(rm);free(chi2m);
   
  
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);   free(fit[e]);
  
   }
   
   free(fit);  

   free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   return r;
    
} 















double **fit_Mpi_fw_chiral_FVE_prior(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit; 
   int i,j,e,im;  
   int Npar=8;
   int Nvar=5;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=0;
   int n,count,N=6;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=3;   //M_PS A
   en[1]=2;   // M_PS_B
   en[2]=3;   //f_PS A
   en[3]=2;    //f_PS B
   
   en[4]=1;    //Z_P A
   en[5]=1;    // Z_P B
   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   int ii;
   
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
    guess[0]=2.05478;
    guess[1]=0.113021;
    guess[2]=2.54451;
    guess[3]=0.410103;
    guess[4]=4.61247;
    guess[5]=-0.272884;
    guess[6]=gJ[0].Zp[Njack-1]+0.0001;
    guess[7]=gJ[3].Zp[Njack-1]+0.0001;
    
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double**) malloc(sizeof(double*)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }
   
   

   count=0;
   for (n=0;n<N;n++){
       // printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                ii=index_a[n][e];
                im=mass_index[ii][ik2][ik1];
                if(n==0 ||n==1 ){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[ii].M_PS_jack[im][j]   *  gJ[ii].w0[j];
                        rm[j]*=rm[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==2 ||n==3 ){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[ii].f_PS_jack[im][j]   *  gJ[ii].w0[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==4 ||n==5 ){
                    for (j=0;j<Njack;j++)
                        rm[j]=gJ[ii].Zp[j];//0 is one ensemble of the lattice spacing A
                    
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                    for (j=0;j<Njack;j++)
                        rm[j]=gJ[0].Zp[Njack-1];
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[ii].k[head[ii].nk+ik2]*gJ[ii].w0[Njack-1];//mu*w0
                x[e+count][1]=gJ[ii].w0[Njack-1];//w0
                x[e+count][2]=gJ[ii].M_PS_jack[im][Njack-1]*gJ[ii].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][3]=gJ[ii].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][4]=double(head[ii].l1);//L
             //   printf("%g     %g      %g   %g\n",x[e+count][1],x[e+count][0],fit[e+count][0],fit[e+count][1]);
        }
        count+=en[n];
   }

    
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, Mw2_fw_chiral_FVE_prior,guess );

        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, Mw2_fw_chiral_FVE_prior  );
              //printf("tmp[6]=%g\n",tmp[6]);
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }                
        free(tmp);

   }  
   
// make_plots_MPi_fPi(en,);
   printf("w0/a[fm]     m*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err    KM2/Kf\n"); 
   double KM,Kf,K;
   count=0;
   for (n=0;n<N;n++){
        printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                                
               // x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[Njack-1];///gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=gJ[e].w0[Njack-1];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][Njack-1]*gJ[e].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
                
                FVE( v_w0GeV, gJ[e].w0[Njack-1] ,r[2][Njack-1], r[4][Njack-1], r[0][Njack-1], r[1][Njack-1],head[e].l1, x[e+count][0] ,  x[e+count][2], x[e+count][3], &KM, &Kf);
                if (n==0)
                    K=KM*KM;
                else if (n==1)
                    K=KM*KM;
                else if (n==2)
                    K=Kf;
                else if (n==3||n==4)
                    K=1;
                printf("%g     %g      %g   %g     %g \n",x[e+count][1],x[e+count][0],fit[e+count][0]/K,fit[e+count][1]/K,K);
        }
        count+=en[n];
    }

   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
   printf("$\\chi^2=%f+-%f$\n",chi2m[0],chi2m[1]);
   free(rm);free(chi2m);
   
  
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);   free(fit[e]);
  
   }
   
   free(fit);  

   free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);

   return r;
    
} 


double **fit_Mpi_fw_chiral_FVE_P40_prior(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit; 
   int i,j,e,im;  
   int Npar=7;
   int Nvar=5;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=0;
   int n,count,N=6;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=3;   //M_PS A
   en[1]=2;   // M_PS_B
   en[2]=3;   //f_PS A
   en[3]=2;    //f_PS B
   
   en[4]=1;    //Z_P A
   en[5]=1;    // Z_P B
   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   int ii;
   
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
    guess[0]=2.05478;
    guess[1]=0.113021;
    guess[2]=2.54451;
    guess[3]=0.410103;
    guess[4]=4.61247;
    //guess[5]=-0.272884;
    guess[5]=gJ[0].Zp[Njack-1]+0.0001;
    guess[6]=gJ[3].Zp[Njack-1]+0.0001;
    
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double**) malloc(sizeof(double*)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }
   
   

   count=0;
   for (n=0;n<N;n++){
       // printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                ii=index_a[n][e];
                im=mass_index[ii][ik2][ik1];
                if(n==0 ||n==1 ){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[ii].M_PS_jack[im][j]   *  gJ[ii].w0[j];
                        rm[j]*=rm[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==2 ||n==3 ){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[ii].f_PS_jack[im][j]   *  gJ[ii].w0[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==4 ||n==5 ){
                    for (j=0;j<Njack;j++)
                        rm[j]=gJ[ii].Zp[j];//0 is one ensemble of the lattice spacing A
                    
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                    for (j=0;j<Njack;j++)
                        rm[j]=gJ[0].Zp[Njack-1];
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[ii].k[head[ii].nk+ik2]*gJ[ii].w0[Njack-1];//mu*w0
                x[e+count][1]=gJ[ii].w0[Njack-1];//w0
                x[e+count][2]=gJ[ii].M_PS_jack[im][Njack-1]*gJ[ii].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][3]=gJ[ii].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][4]=double(head[ii].l1);//L
             //   printf("%g     %g      %g   %g\n",x[e+count][1],x[e+count][0],fit[e+count][0],fit[e+count][1]);
        }
        count+=en[n];
   }

   
   #pragma omp parallel for  private(tmp,i)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)     
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, Mw2_fw_chiral_FVE_prior,guess );

        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, Mw2_fw_chiral_FVE_prior  );
              //printf("tmp[6]=%g\n",tmp[6]);
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }                
        free(tmp);

   }  
   
// make_plots_MPi_fPi(en,);
   printf("w0/a[fm]     m*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err    KM2/Kf\n"); 
   double KM,Kf,K;
   count=0;
   for (n=0;n<N;n++){
        printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                                
               // x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[Njack-1];///gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=gJ[e].w0[Njack-1];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][Njack-1]*gJ[e].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
                
                FVE( v_w0GeV, gJ[e].w0[Njack-1] ,r[2][Njack-1], r[4][Njack-1], r[0][Njack-1], r[1][Njack-1],head[e].l1, x[e+count][0] ,  x[e+count][2], x[e+count][3], &KM, &Kf);
                if (n==0)
                    K=KM*KM;
                else if (n==1)
                    K=KM*KM;
                else if (n==2)
                    K=Kf;
                else if (n==3||n==4)
                    K=1;
                printf("%g     %g      %g   %g     %g \n",x[e+count][1],x[e+count][0],fit[e+count][0]/K,fit[e+count][1]/K,K);
        }
        count+=en[n];
    }

   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
   printf("$\\chi^2=%f+-%f$\n",chi2m[0],chi2m[1]);
   free(rm);free(chi2m);
   
  
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);   free(fit[e]);
  
   }
   
   free(fit);  

   free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);

   return r;
    
} 


double **fit_Mpi_fw_chiral_FVE_flag(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit; 
   int i,j,e,im;  
   int Npar=4;
   int Nvar=5;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=0;
   int n,count,N=2;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=ensembles;
   en[1]=ensembles;
   int en_tot=0;

   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
    guess[0]=2.05478;
    guess[1]=0.113021;
    guess[2]=0.410103;
    guess[3]=-0.272884;
    
   
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double**) malloc(sizeof(double*)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }
   

   count=0;
   for (n=0;n<N;n++){
       // printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j];
                        rm[j]*=rm[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=gJ[e].w0[Njack-1];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][Njack-1]*gJ[e].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
             //   printf("%g     %g      %g   %g\n",x[e+count][1],x[e+count][0],fit[e+count][0],fit[e+count][1]);
        }
        count+=en[n];
   }

   
   #pragma omp parallel for  private(tmp,i)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)    
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, Mw2_fw_chiral_FVE,guess );

        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, Mw2_fw_chiral_FVE  );
              
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }                
        free(tmp);

   } 
   
// make_plots_MPi_fPi(en,);
   printf("w0/a[fm]     m*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err     KM2/Kf\n"); 
   double KM,Kf,K;
   count=0;
   for (n=0;n<N;n++){
        printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                                
               // x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=gJ[e].w0[Njack-1];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][Njack-1]*gJ[e].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
                
                FVE( v_w0GeV, gJ[e].w0[Njack-1], flagl3b, flagl4b, r[0][Njack-1], r[1][Njack-1],head[e].l1, x[e+count][0] ,  x[e+count][2], x[e+count][3], &KM, &Kf);
                if (n==0)
                    K=KM*KM;
                else if (n==1)
                    K=Kf;
                printf("%g     %g      %g   %g     %g \n",x[e+count][1],x[e+count][0],fit[e+count][0]/K,fit[e+count][1]/K,K);
        }
        count+=en[n];
    }

   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
   printf("$\\chi^2/dof=%f+-%f$\n",chi2m[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
   free(rm);free(chi2m);
   
  
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);   free(fit[e]);
  
   }
   
   free(fit);  

   free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   return r;
    
} 



double **fit_fw_of_Mw_chiral_FVE(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit; 
   int i,j,e,im;  
   int Npar=3;
   int Nvar=4;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=0;
   int n,count,N=1;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=ensembles;
   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
    guess[0]=0.106845;
    guess[1]=-4.44555;
    guess[2]=-0.0875071;
    if (Npar>3)guess[3]=-0.028469;
    if (Npar>4) guess[4]=1;
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double**) malloc(sizeof(double*)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }

   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                      /*  if(e==0)rm[j]/=(1-0.009);
                        if(e==1)rm[j]/=(1-0.0145);
                        if(e==2)rm[j]/=(1-0.0056);*/
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                   
                }
                
                 //x[e+count]=(double*) malloc(sizeof(double)*Nvar);               
           
        }
        count+=en[n];
   }

   
    
   #pragma omp parallel for  private(tmp,i,n,e,x,count,im)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2)
   for (j=0;j<Njack;j++){
       x=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);   
                x[e+count][0]=gJ[e].w0[j];//w0
                x[e+count][1]=gJ[e].M_PS_jack[im][j]*gJ[e].M_PS_jack[im][j];//MPS^2
                x[e+count][2]=gJ[e].f_PS_jack[im][j];//f_PS
                x[e+count][3]=double(head[e].l1);//f_PS
            }
            count+=en[n];
        }  
 
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fw_of_Mw2_chiral_FVE,guess );

        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fw_of_Mw2_chiral_FVE  );
             
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }                
        free(tmp);
        for (e=0;e<en_tot;e++)
            free(x[e]);
        free(x);

   }   
   x=(double**) malloc(sizeof(double*)*(en_tot));

// make_plots_MPi_fPi(en,);
   printf("w0/a[fm]     (M_PS*w0/KM)^2       f_PS*w0/Kf    err K\n"); 
   double KM,Kf,K;
   double w0; 
   count=0;
   for (n=0;n<N;n++){
        printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=gJ[e].w0[Njack-1];//w0
                w0=gJ[e].w0[Njack-1];
                x[e+count][1]=gJ[e].M_PS_jack[im][Njack-1]*gJ[e].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][2]=gJ[e].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][3]=double(head[e].l1);//f_PS
                
                FVE( v_w0GeV, gJ[e].w0[Njack-1] ,flagl3b, r[1][Njack-1], 0.5, r[0][Njack-1],head[e].l1, x[e+count][1]*x[e+count][0]*x[e+count][0] ,  x[e+count][1], x[e+count][2], &KM, &Kf);
                
                double M,L,Rf,RM;
                M=sqrt(x[e+count][1])*w0/v_w0MeV;
                L=(double(head[e].l1))*v_w0fm/w0;
                Rf=exp(4.58982-0.0138032*M-2.013*L);
                RM=exp(3.81729-0.0130342*M-2.1714*L);
                Kf=1-Rf;
                KM=RM+1; 
                    
                K=Kf;

                //FVE( v_w0GeV, gJ[e].w0[Njack-1] ,flagl3b, r[1][Njack-1], 2.11586, r[0][Njack-1],head[e].l1, head[e].k[head[e].nk+ik2]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1] ,  x[e+count][1], x[e+count][2], &KM, &Kf);

                printf("%g     %g      %g   %g   %g   \n",x[e+count][0],x[e+count][1]*w0*w0/(KM*KM),fit[e+count][0]/K,fit[e+count][1]/K,K);
        }
        count+=en[n];
    }

   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
   printf("$\\chi^2/dof=%f+-%f$\n",chi2m[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
   free(rm);free(chi2m);
   
  
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);   free(fit[e]);
  
   }
   
   free(fit);  

   free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   return r;
    
} 

double **fit_Mpi_fw_chiral_FVE_P40(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit; 
   int i,j,e,im;  
   int Npar=5;
   int Nvar=5;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=0;
   int n,count,N=2;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=ensembles;
   en[1]=ensembles;
   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
    guess[0]=2.05478;
    guess[1]=0.113021;
    guess[2]=2.54451;
    guess[3]=0.410103;
    guess[4]=4.61247;
   

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double**) malloc(sizeof(double*)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }
   

   count=0;
   for (n=0;n<N;n++){
       // printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j];
                        rm[j]*=rm[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
             //   printf("%g     %g      %g   %g\n",x[e+count][1],x[e+count][0],fit[e+count][0],fit[e+count][1]);
        }
        count+=en[n];
   }

   
   
   #pragma omp parallel for  private(tmp,i,n,e,x,count,im)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2,en_tot)
   for (j=0;j<Njack;j++){
        x=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);   
                x[e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                x[e+count][1]=gJ[e].w0[j];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][j]*gJ[e].M_PS_jack[im][j];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][j];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
            }
            count+=en[n];
        }
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, Mw2_fw_chiral_FVE_P40,guess );

        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, Mw2_fw_chiral_FVE_P40  );
              
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }                
        free(tmp);
        for (e=0;e<en_tot;e++)
            free(x[e]);
        free(x);

   } 

   
   x=(double**) malloc(sizeof(double*)*(en_tot));
// make_plots_MPi_fPi(en,);
   printf("w0/a[fm]     m*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err    KM2/Kf\n"); 
   double KM,Kf,K;
   count=0;
   for (n=0;n<N;n++){
        printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=gJ[e].w0[Njack-1];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][Njack-1]*gJ[e].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
                
                FVE( v_w0GeV, gJ[e].w0[Njack-1] ,r[2][Njack-1], r[4][Njack-1], r[0][Njack-1], r[1][Njack-1],head[e].l1, x[e+count][0] ,  x[e+count][2], x[e+count][3], &KM, &Kf);
                if (n==0)
                    K=KM*KM;
                else if (n==1)
                    K=Kf;
                printf("%g     %g      %g   %g     %g \n",x[e+count][1],x[e+count][0],fit[e+count][0]/K,fit[e+count][1]/K,K);
        }
        count+=en[n];
    }

   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
   printf("$\\chi^2=%f+-%f$\n",chi2m[0],chi2m[1]);
   free(rm);free(chi2m);
   
  
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);   free(fit[e]);
  
   }
   
   free(fit);  

   free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   return r;
    
} 



double **fit_fw_of_Mw_chiral_FVE_P40(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit; 
   int i,j,e,im;  
   int Npar=2;
   int Nvar=4;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=0;
   int n,count,N=1;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=ensembles;
   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
    guess[0]=0.110791;
    guess[1]=flagl4b;

   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double**) malloc(sizeof(double*)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){    
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
               // printf("%d  %d  %d\n",j,i,count);
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }


   count=0;
   for (n=0;n<N;n++){
       // printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
               
             //   printf("%g     %g      %g   %g\n",x[e+count][1],x[e+count][0],fit[e+count][0],fit[e+count][1]);
        }
        count+=en[n];
   }

   
   #pragma omp parallel for  private(tmp,i,n,e,x,count,im)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2)
   for (j=0;j<Njack;j++){
       x=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);   
                x[e+count][0]=gJ[e].w0[j];//w0
                x[e+count][1]=gJ[e].M_PS_jack[im][j]*gJ[e].M_PS_jack[im][j];//MPS^2
                x[e+count][2]=gJ[e].f_PS_jack[im][j];//f_PS
                x[e+count][3]=double(head[e].l1);//f_PS
            }
            count+=en[n];
        }
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fw_of_Mw2_chiral_FVE_P40,guess );

        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fw_of_Mw2_chiral_FVE_P40  );
              
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }                
        free(tmp);
         for (e=0;e<en_tot;e++)
            free(x[e]);
        free(x);

   }  
   
// make_plots_MPi_fPi(en,);
   printf("w0/a[fm]     (M_PS*w0/KM)^2       f_PS*w0/Kf    err  K\n"); 
   double KM,Kf,K;
   double w0; 
   count=0;
   for (n=0;n<N;n++){
        printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=gJ[e].w0[Njack-1];//w0
                w0=gJ[e].w0[Njack-1];
                x[e+count][1]=gJ[e].M_PS_jack[im][Njack-1]*gJ[e].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][2]=gJ[e].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][3]=double(head[e].l1);//f_PS
                
                FVE( v_w0GeV, gJ[e].w0[Njack-1] ,flagl3b, r[1][Njack-1], 0.5, r[0][Njack-1],head[e].l1, x[e+count][1]*x[e+count][0]*x[e+count][0] ,  x[e+count][1], x[e+count][2], &KM, &Kf);
                    K=Kf;
                printf("%g     %g      %g   %g   %g \n",x[e+count][0],x[e+count][1]*w0*w0/(KM*KM),y[Njack-1][e+count][0]/K,y[Njack-1][e+count][1]/K,K);
        }
        count+=en[n];
    }

   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
   printf("$\\chi^2=%f+-%f$\n",chi2m[0],chi2m[1]);
   free(rm);free(chi2m);
   
  
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);   free(fit[e]);
  
   }
   
   free(fit);  

   free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   return r;
    
} 
