#define CONTROL

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
#include "KandD.hpp"
#include "global_fit_KandD.hpp"
#include "tower.hpp"
#include "mutils.hpp"

#include "fve.hpp"


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


int ***mass_index;

double ***init_omega_jacob(){
    int e,ik;
    double ***r;
    
    r=double_malloc_3(ensembles,3, 2);
    
    if(ensembles>8){
        //B25.32
        /*r[8][0][0]=0.669; r[8][0][1]=0.003;  
        r[8][1][0]=0.698; r[8][1][1]=0.002;    
        r[8][2][0]=0.727; r[8][2][1]=0.002;  */
        r[8][0][0]=0.6763; r[8][0][1]=0.0020;  
        r[8][1][0]=0.7042; r[8][1][1]=0.0016;    
        r[8][2][0]=0.7318; r[8][2][1]=0.0014;  

    }if(ensembles>7){
        //C60
        /*r[7][0][0]=0.558; r[7][0][1]=0.007;  
        r[7][1][0]=0.586; r[7][1][1]=0.005;    
        r[7][2][0]=0.609; r[7][2][1]=0.005;  
        */
        r[7][0][0]=0.5765; r[7][0][1]=0.0025;  
        r[7][1][0]=0.5924; r[7][1][1]=0.0022;    
        r[7][2][0]=0.6081; r[7][2][1]=0.0020;
        
    }if(ensembles>6){
        //B072
        /*r[6][0][0]=0.682; r[6][0][1]=0.006;  
        r[6][1][0]=0.701; r[6][1][1]=0.005;  
        r[6][2][0]=0.719; r[6][2][1]=0.004;  */
        r[6][0][0]=0.6835; r[6][0][1]=0.0018;  
        r[6][1][0]=0.7031; r[6][1][1]=0.0016;  
        r[6][2][0]=0.7224; r[6][2][1]=0.0014;


    }if(ensembles>5){
        //B14.64
       /* r[5][0][0]=0.666; r[5][0][1]=0.002;  
        r[5][1][0]=0.695; r[5][1][1]=0.001;  
        r[5][2][0]=0.722; r[5][2][1]=0.001;  */
        r[5][0][0]=0.6606; r[5][0][1]=0.0034;  
        r[5][1][0]=0.6903; r[5][1][1]=0.0027;  
        r[5][2][0]=0.7194; r[5][2][1]=0.0022; 
    }if(ensembles>4){
        //B25.48
        /*r[4][0][0]=0.667; r[4][0][1]=0.002;  
        r[4][1][0]=0.696; r[4][1][1]=0.002;  
        r[4][2][0]=0.725; r[4][2][1]=0.001;  */
        r[4][0][0]=0.6652; r[4][0][1]=0.0042;  
        r[4][1][0]=0.6933; r[4][1][1]=0.0032;  
        r[4][2][0]=0.7218; r[4][2][1]=0.0027;  

    }if(ensembles>3){
        //A12
        /*r[3][0][0]=0.776; r[3][0][1]=0.005;  
        r[3][1][0]=0.811; r[3][1][1]=0.004;  
        r[3][2][0]=0.844; r[3][2][1]=0.004;  */
        r[3][0][0]=0.7843; r[3][0][1]=0.0041;  
        r[3][1][0]=0.8204; r[3][1][1]=0.0034;  
        r[3][2][0]=0.8564; r[3][2][1]=0.0029; 


    }if(ensembles>2){
        //A30
        /*r[2][0][0]=0.789; r[2][0][1]=0.003;  
        r[2][1][0]=0.823; r[2][1][1]=0.002;  
        r[2][2][0]=0.860; r[2][2][1]=0.002; */
        r[2][0][0]=0.7889; r[2][0][1]=0.0022;  
        r[2][1][0]=0.8255; r[2][1][1]=0.0018;  
        r[2][2][0]=0.8613; r[2][2][1]=0.0015;


    }if(ensembles>1){
        //A40
        /*r[1][0][0]=0.803; r[1][0][1]=0.003;  
        r[1][1][0]=0.840; r[1][1][1]=0.002;  
        r[1][2][0]=0.873; r[1][2][1]=0.002;*/
        r[1][0][0]=0.8064; r[1][0][1]=0.0028;  
        r[1][1][0]=0.8401; r[1][1][1]=0.0022;  
        r[1][2][0]=0.8739; r[1][2][1]=0.0019; 
        
    }if(ensembles>0){
        //A53
        /*r[0][0][0]=0.811; r[0][0][1]=0.003;  
        r[0][1][0]=0.843; r[0][1][1]=0.002;  
        r[0][2][0]=0.878; r[0][2][1]=0.002;  */
        r[0][0][0]=0.8132; r[0][0][1]=0.0019;  
        r[0][1][0]=0.8464; r[0][1][1]=0.0015;  
        r[0][2][0]=0.8799; r[0][2][1]=0.0013;

    }
    return r;
    
}




double fit_Fpi_and_Mpi(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b, P4,ZP;
         Bw=P[0], fw=P[1], l3b=P[2], P2=P[3], l4b=P[4], P4=P[5];
   
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
    int Lsize=(int(x[4]));
    double KM,Kf;
    
   
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    
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
    if (n==2){
         Mw2=KM;
    }
    if (n==3){
         Mw2=Kf;
    }
     return Mw2;
    
}




double fit_Fpi_and_Mpi_GL(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b, P4,ZP;
         Bw=P[0], fw=P[1], l3b=P[2], P2=P[3], l4b=P[4], P4=P[5];
   
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
    int Lsize=(int(x[4]));
    double L_w=(x[4]) /w0;
    double KM,Kf;
   
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    
    double Delta=FVE_GL_fast( L_w, mw, fw, Bw);
    
    
    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    if (n==0){
        
        Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2;
        Mw2*=2*Bw*mw*(1-0.25 *Delta)*(1-0.25 *Delta);
        
    }
    if (n==1){
        
        Mw2=fw*(1-2*xi*log(xi)+P3*xi)+(1./(w0*w0))*P4;
        Mw2*=(1+Delta);
        
    }
    if (n==2){
        Mw2=(1-0.25 *Delta);  //KM   M(inf)=M(L)/KM
        
    }
    if (n==3){
        Mw2=(1+Delta);  //Kf   M(inf)=M(L)/Kf
        
    }
    
     return Mw2;
    
}



double fit_FpiMpi4_and_Mpi2_GL(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b, P4,ZP;
         Bw=P[0], fw=P[1], l3b=P[2], P2=P[3], l4b=P[4], P4=P[5] ;
   
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
    int Lsize=(int(x[4]));
    double L_w=(x[4]) /w0;
    double KM,Kf;
   
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    
    double Delta=FVE_GL_fast( L_w, mw, fw, Bw);
    
    
    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    if (n==0){
        
        Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2;
        Mw2*=2*Bw*mw*(1-0.25 *Delta)*(1-0.25 *Delta);
        
    }
    if (n==1){// fit fot w0^5fpi Mpi^4
        double ML=sqrt(2*Bw*mw)*L_w;
        Mw2= (4.*pi)* (4.*pi)* (4.*pi)* (4.*pi);
        Mw2*= fw * fw *fw * fw * fw *xi*xi;
        Mw2*= (1.+ 2.*(l4b-l3b) *xi + xi*xi *P[6] +(1./(w0*w0))*P4 );
        Mw2*=(1.+ P[7]*xi*xi *exp(-ML)/ pow(ML,3./2.));
        
    }
    if (n==2){
        Mw2=(1-0.25 *Delta);  //KM   M(inf)=M(L)/KM
        
    }
    if (n==3){
        Mw2=1;  //KfMpi4   M(inf)=M(L)/Kf
        
    }
    
     return Mw2;
    
}


double fit_FpiMpi4_and_Mpi2_linear(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw, fw, l3b, P2, l4b, P4,ZP;
         Bw=P[0], fw=P[1], l3b=P[2], P2=P[3], l4b=P[4], P4=P[5] ;
   
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
    int Lsize=(int(x[4]));
    double L_w=(x[4]) /w0;
    double KM,Kf;
   
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    
    double Delta=FVE_GL_fast( L_w, mw, fw, Bw);
    
    
    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    if (n==0){
        
        Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2;
        Mw2*=2*Bw*mw*(1-0.25 *Delta)*(1-0.25 *Delta);
        
    }
    if (n==1){// fit fot w0^5fpi Mpi^4
        double ML=sqrt(2*Bw*mw)*L_w;
        Mw2= (4.*pi)* (4.*pi)* (4.*pi)* (4.*pi);
        Mw2*= fw * fw *fw * fw * fw *xi*xi;
        Mw2*= (1.+ 2.*(l4b-l3b) *xi + xi*xi *0. +(1./(w0*w0))*P4 );
        Mw2*=(1.+ 0.*xi*xi *exp(-ML)/ pow(ML,3./2.));
        
    }
    if (n==2){
        Mw2=(1-0.25 *Delta);  //KM   M(inf)=M(L)/KM
        
    }
    if (n==3){
        Mw2=1;  //KfMpi4   M(inf)=M(L)/Kf
        
    }
    
     return Mw2;
    
}


double fit_FK_and_MK_GL(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double   l3b, P2, l4b, P4,ZP;
        
   
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3],  msw=x[5],   MK2w2=x[6] , Bw=x[8],  fw=x[9];
    int Lsize=(int(x[4]));
    double L_w=(x[4]) /w0;
    double KM,Kf;
   
    double P1=-l3b-2*log( v_Mpiw0 /(4*pi*fw));
    double P3=2*l4b+4*log(   v_Mpiw0/(4*pi*fw) );
    
    
    double Delta=FVE_GL_fast( L_w, mw, fw, Bw);
    Delta=Delta*(3./8.); // to go from Delta_pi to Delta_K
    
    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    if (n==0){
        
        Mw2=1+P[1]*mw+P[2]*mw*mw+ P[3]/(w0*w0);
        Mw2*=P[0]*(mw+msw);//*(1-0.25 *Delta)*(1-0.25 *Delta);
        
    }
    if (n==1){
        
        Mw2=P[4]*(1-(3./4.)*xi*log(xi)+P[5]*xi+(1./(w0*w0))*P[6]);
        Mw2*=(1+Delta);
        
    }
    if (n==2){
        Mw2=1.;//(1-0.25 *Delta);  //KM   M(inf)=M(L)/KM
        
    }
    if (n==3){
        Mw2=(1+Delta);  //Kf   M(inf)=M(L)/Kf
        
    }
    
     return Mw2;
    
}



double fit_FKoverFpi_GL(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    
    double     P0=P[0], P1=P[1], P2=P[2], fw=P[3];
   
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3], Bw=x[8];
    int Lsize=(int(x[4]));
    double L_w=(x[4]) /w0;
    double KM,Kf;
   
    
    
    double Delta=FVE_GL_fast( L_w, mw, fw, Bw);
    double DeltaK=Delta*(3./8.);
    
    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    
    if (n==0){
        
        Mw2 = P0 *(1 + (5./4.)* xi* log(xi) + P1 *xi + P2/(w0*w0)  ) ;
        Mw2*=(1+DeltaK)/(1+Delta);
    }
    if (n==1){
        
        Mw2=(1+DeltaK)/(1+Delta);
        
    }
    
    
     return Mw2;
    
}

double fit_MK_Mpi_FK_Fpi_GL(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    
    double     P0=P[0], P1=P[1], P2=P[2],  P3=P[3] ,P4=P[4], P5=P[5] ,P6=P[6], fw=P[7];
   
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3],  msw=x[5],   MK2w2=x[6] , Bw=x[8] ;

    int Lsize=(int(x[4]));
    double L_w=(x[4]) /w0;
    double KM,Kf;
    
    double Delta=FVE_GL_fast( L_w, mw, fw, Bw);
    double DeltaK=Delta*(3./8.);
    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
     if (n==0){
        Mw2 = 0.5*(1. + msw/mw) *P0* (1. - xi* log(xi) + P1* xi + P2 *xi*xi + P3/(w0*w0) );
        Mw2*=1./((1-0.25 *Delta)*(1-0.25 *Delta));
    }
    if (n==1){
        Mw2 = P4 *(1 + (5./4.)* xi* log(xi) + P5 *xi + P6/(w0*w0)  ) ;
        Mw2*=(1+DeltaK)/(1+Delta);
    }
    if (n==2){
        Mw2=1./((1-0.25 *Delta));
    }
    if (n==3){
        Mw2=(1+DeltaK)/(1+Delta);
    }
     return Mw2;
    
}

double fit_MK_Mpi_FK_Fpi_GL_fix_f(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    
    double     P0=P[0], P1=P[1], P2=P[2],  P3=P[3] ,P4=P[4], P5=P[5] ,P6=P[6];
   
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3],  msw=x[5],   MK2w2=x[6] , Bw=x[8],fw=x[9] ;

    int Lsize=(int(x[4]));
    double L_w=(x[4]) /w0;
    double KM,Kf;
   
    
    
    double Delta=FVE_GL_fast( L_w, mw, fw, Bw);
    double DeltaK=Delta*(3./8.);
    
    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
     if (n==0){
        
        Mw2 = 0.5*(1. + msw/mw) *P0* (1. - xi* log(xi) + P1* xi + P2 *xi*xi + P3/(w0*w0) );
        Mw2*=1./((1-0.25 *Delta)*(1-0.25 *Delta));
    }
    if (n==1){
        
        Mw2 = P4 *(1 + (5./4.)* xi* log(xi) + P5 *xi + P6/(w0*w0)  ) ;
        Mw2*=(1+DeltaK)/(1+Delta);
    }
    if (n==2){
        
        Mw2=1./((1-0.25 *Delta));
        
    }
    if (n==3){
        
        Mw2=(1+DeltaK)/(1+Delta);
        
    }
    
    
     return Mw2;
    
}


int ***init_mass_index_ave_r(struct header *head)
{
     int k1, k2,i;
     int nk,e;
     int ***mass_index;
     
     mass_index=(int ***)  malloc(sizeof(int**)*ensembles);
     for (e=0;e<ensembles;e++){
        nk=head[e].nk;
        mass_index[e]=(int**) malloc(sizeof(int*)*nk);
        for (k2=0;k2<nk;k2++){
            mass_index[e][k2]=(int*) malloc(sizeof(int)*(k2+1));
        }

        i=0;
        for (k2=0;k2<nk;k2++)
        for (k1=0;k1<=k2;k1++)
        {
            mass_index[e][k2][k1]=i;
            i++;
        }

     }
     return mass_index;
}

double global_fK_from_M(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi,r;
    double pi=3.141592653589793;
    
    double Mpiw=x[0], MKw=x[1], w0=x[2], Mpi2=x[3], fpi=x[4], frac_Lw=x[7],  Bw=x[8];
    double fw=x[9],  MK2=x[5], fK=x[6];
    
    //double    P0_w=Bw, P1_w=P[0], P3ww=P[1];
    //double   Pf1w=P[0],  Pf2w=P[1],  Pf4www=P[2];
    
    double KM=1.,Kf=1.0;
    
   //FVE_K( Bw, fw, frac_Lw,  Mpiw*Mpiw/(2.*Bw)/*mlw*/, MKw*MKw/Bw-Mpiw*Mpiw/(2.*Bw) /*msw*/ ,Mpi2,  fpi,MK2, fK,&KM, &Kf);
        
    //    fKw=Pf1w*( 1.- (3./2.)* xi*log(xi)+Pf2w*xi+(1/(w0*w0))*Pf4www)*Kf;
    
    xi=Mpiw*Mpiw/(16*pi*pi*fw*fw);
    double P1=P[0]+P[3]*MKw*MKw;
    double P2=P[1]+P[4]*MKw*MKw;
    double P4=P[2]+P[5]*MKw*MKw;
    
    r=P1*(1.- (3./2.)* xi*log(xi)+ P2*xi +  P4 *(1/(w0*w0))    )*Kf;
    
    return r;
    
}

double  global_Omega_Mpi_MK(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi,r;
    double pi=3.141592653589793;
    
    double Mpiw=x[0], MKw=x[1], w0=x[2];
    
    r=P[0]+P[1]*Mpiw*Mpiw+P[2]*MKw*MKw+P[3]/(w0*w0);
    
    return r;
    
}
double  global_Omega_MK(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi,r;
    double pi=3.141592653589793;
    
    double Mpiw=x[0], MKw=x[1], w0=x[2];
    
    r=P[0]+P[1]*MKw*MKw+P[2]/(w0*w0);
    
    return r;
    
}
double  global_Omega_propMK(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi,r;
    double pi=3.141592653589793;
    
    double Mpiw=x[0], MKw=x[1], w0=x[2];
    
    r=(P[0])*MKw*MKw+P[1]/(w0*w0);
    
    return r;
    
}
/*
static void  read_file_head_jack(FILE *stream,struct header *head)
{
    int i;
    
    fread(&(head->twist),sizeof(int),1,stream);
    fread(&(head->nf),sizeof(int),1,stream);
    fread(&(head->nsrc),sizeof(int),1,stream);
    fread(&(head->l0),sizeof(int),1,stream);
    fread(&(head->l1),sizeof(int),1,stream);
    fread(&(head->l2),sizeof(int),1,stream);
    fread(&(head->l3),sizeof(int),1,stream);
    fread(&(head->nk),sizeof(int),1,stream);
    fread(&(head->nmoms),sizeof(int),1,stream);
    
    fread(&(head->beta),sizeof(double),1,stream);
    fread(&(head->ksea),sizeof(double),1,stream);
    fread(&(head->musea),sizeof(double),1,stream);
    fread(&(head->csw),sizeof(double),1,stream);
   
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
}*/

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
    
}


int setup_reading_single_jack( struct header *head, FILE **f, const char *name){
    int N;
    
    *f=open_file(name,"r");
    read_file_head_jack(*f,head);
    fread(&N,sizeof(int),1,*f);
    
    return N;
}/*
void setup_reading_single_jack( struct  database_file_jack *jack_files, struct header *head){
     int N;
     
     jack_files->f_M_PS=fopen(jack_files->M_PS,"r");
     error(jack_files->f_M_PS==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->M_PS);
     read_file_head_jack(jack_files->f_M_PS,head);
     fread(&(jack_files->Njack),sizeof(int),1,jack_files->f_M_PS);
     
     jack_files->f_M_PS_GEVP=fopen(jack_files->M_PS_GEVP,"r");
     error(jack_files->f_M_PS_GEVP==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->M_PS_GEVP);
     read_file_head_jack(jack_files->f_M_PS_GEVP,head);
     fread(&(jack_files->Njack),sizeof(int),1,jack_files->f_M_PS_GEVP);
     
     jack_files->f_f_PS=fopen(jack_files->f_PS,"r");
     error(jack_files->f_f_PS==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->f_PS);
     read_file_head_jack(jack_files->f_f_PS,head);
     fread(&N,sizeof(int),1,jack_files->f_f_PS);
     
     jack_files->f_f_PS_ls_ss=fopen(jack_files->f_PS_ls_ss,"r");
     error(jack_files->f_f_PS_ls_ss==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->f_PS_ls_ss);
     read_file_head_jack(jack_files->f_f_PS_ls_ss,head);
     fread(&N,sizeof(int),1,jack_files->f_f_PS_ls_ss);
     
     error(jack_files->Njack!=N,1,"setup_reading_single_jack", " files \n %s has %d elements \n %s has %d elements\n ",
         jack_files->M_PS_GEVP, jack_files->Njack,jack_files->f_PS,N); 
     
}*/

void setup_reading_list_jack( struct  database_file_jack *jack_files, struct header *head){
     int N,N1;
     
      N=setup_reading_single_jack(head,&(jack_files->f_M_PS)  ,jack_files->M_PS );
      
      N1=setup_reading_single_jack(head,&(jack_files->f_M_PS_GEVP)       ,jack_files->M_PS_GEVP );
      error(N!=N1,1,"setup_reading_list_jack","jacknifes  in %s have not the same number",jack_files->M_PS_GEVP  );
      N1=setup_reading_single_jack(head,&(jack_files->f_f_PS)       ,jack_files->f_PS );
      error(N!=N1,1,"setup_reading_list_jack","jacknifes in %s have not the same number", jack_files->f_PS );
      N1=setup_reading_single_jack(head,&(jack_files->f_f_PS_ls_ss)       ,jack_files->f_PS_ls_ss );
      error(N!=N1,1,"setup_reading_list_jack","jacknifes in %s have not the same number", jack_files->f_PS_ls_ss);
      
      
      jack_files->Njack=N;
    
     
}

void  setup_reading_jack(char **argv,struct  database_file_jack *jack_files, struct header *head,const char  *name)  {

    
        mysprintf(jack_files->M_PS,NAMESIZE,"%s/M_{PS}_%s",name,argv[1]);
        mysprintf(jack_files->f_PS,NAMESIZE,"%s/Zf_{PS}_%s",name,argv[1]);  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  I rename for simplicity  Zf_PS->f_PS !!!!!!!!!!!!!!!!!!!!!!!!!1
        

        mysprintf(jack_files->M_PS_GEVP,NAMESIZE,"%s/M_{PS}^{GEVP}_%s",name,argv[1]);
        mysprintf(jack_files->f_PS_ls_ss,NAMESIZE,"%s/f_{PS}_ls_ss_%s",name,argv[1]);   
        
        setup_reading_list_jack(jack_files, head);

}



void  files_declarations(char **argv,struct database_file_jack **jack_files,struct header **head){
  
    (*jack_files)=(struct database_file_jack *) malloc (sizeof(struct database_file_jack )*ensembles);    
    if( strcmp(argv[1],"jack")==0){
                mysprintf((*jack_files)[0].sampling,NAMESIZE,"jack");
    }
    if( strcmp(argv[1],"boot")==0){
                mysprintf((*jack_files)[0].sampling,NAMESIZE,"boot");
    }
    int i;
    (*head)=(struct header*) malloc(sizeof(struct header)*ensembles_reph);
    for (i=0;i<ensembles_reph;i++)
        (*head)[i].allocated=0;
    
    setup_reading_jack( argv,&((*jack_files)[0]),&((*head)[0]),"/media/marco/Elements/work/PRACE/beta1.726/cA211ab.53.24/analysis/main/jackknife");  
    (*jack_files)[0].a=0.096;
    if (ensembles>1){
        setup_reading_jack(argv, &((*jack_files)[1]),&((*head)[1]),"/media/marco/Elements/work/PRACE/beta1.726/cA211ab.40.24/analysis/main/jackknife");  
        (*jack_files)[1].a=0.096;
    }
    if (ensembles>2){
        setup_reading_jack(argv, &((*jack_files)[2]),&((*head)[2]),"/media/marco/Elements/work/PRACE/beta1.726/cA211ab.30.32/analysis/main/jackknife");  
        (*jack_files)[2].a=0.096;
    }
    if (ensembles>3){
     //   setup_reading_jack(argv, &jack_files[3],&head[3],"../../beta1.726/cA211ab.12.48/analysis/main/jackknife");  
        setup_reading_jack(argv, &((*jack_files)[3]),&((*head)[3]),"/media/marco/Elements/work/PRACE/beta1.726/cA211ab.12.48_no_rew/analysis/main/jackknife");  
        (*jack_files)[3].a=0.096;
    }
    if (ensembles>4){
        setup_reading_jack(argv, &((*jack_files)[4]),&((*head)[4]),"/media/marco/Elements/work/PRACE/beta1.778/cB211ab.25.48/analysis/main/jackknife");  
        (*jack_files)[4].a=0.081;
    }
    if (ensembles>5){
        setup_reading_jack(argv,&((*jack_files)[5]),&((*head)[5]),"/media/marco/Elements/work/PRACE/beta1.778/cB211ab.14.64/analysis/main/jackknife");  
        (*jack_files)[5].a=0.070;
    }
    if (ensembles>6){
        setup_reading_jack(argv,&((*jack_files)[6]),&((*head)[6]),"/media/marco/Elements/work/PRACE/beta1.778/cB211ab.072.64/analysis/main/jackknife");  
        (*jack_files)[6].a=0.081;
    }
    
    if (ensembles>7){
        setup_reading_jack(argv, &((*jack_files)[7]),&((*head)[7]),"/media/marco/Elements/work/PRACE/beta1.836/cC211ab.06.80/analysis/main/jackknife");  
        (*jack_files)[7].a=0.070;
    }
    if (ensembles>8){
        setup_reading_jack(argv, &((*jack_files)[8]),&((*head)[8]),"/media/marco/Elements/work/PRACE/beta1.778/cB211ab.25.32/analysis/main/jackknife");  
        (*jack_files)[8].a=0.081;
    }
    
    
    
    //take the header form ensemble 4
   /* if (ensembles>5){
        jack_files[5].a=0.096;jack_files[5].Njack=100;     
        fseek(jack_files[4].f_M_PS_GEVP,0,SEEK_SET);
        read_file_head_jack(jack_files[4].f_M_PS_GEVP,&head[5]);
        head[5].k[head[5].nk]=0.0012;
        head[5].musea=0.0012;
        head[5].l0=96;head[5].l1=48;head[5].l2=48;head[5].l3=48;
    }*/
    
 /*   if (ensembles>6){
        jack_files[6].a=0.070;jack_files[6].Njack=100;
        fseek(jack_files[4].f_M_PS_GEVP,0,SEEK_SET);
        read_file_head_jack(jack_files[4].f_M_PS_GEVP,&head[6]);
        head[6].k[head[6].nk]=0.0006;
        head[5].musea=0.0006;
        head[6].l0=160;head[6].l1=80;head[6].l2=80;head[6].l3=80;
    }
  
    if (ensembles>6){
        fseek(jack_files[4].f_M_PS_GEVP,sizeof(int),SEEK_CUR);
    }*/
}
/*
void  read_files_jack( struct database_file_jack *jack_files, struct header *head,int ***mass_index, double ***M_PS_GEVP_jack, double ***f_PS_jack){
      
      int i,ik1,ik2;
      
      for(i=0;i<ensembles;i++){
            M_PS_GEVP_jack[i]=(double**) malloc(sizeof(double*)*head[i].nk*head[i].nk);
            f_PS_jack[i]=(double**) malloc(sizeof(double*)*head[i].nk*head[i].nk);
            for(ik1=0;ik1<2;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[i].nk;ik2++){
   
                M_PS_GEVP_jack[i][mass_index[i][ik2][ik1]]=(double*) malloc(sizeof(double)*jack_files[i].Njack);
                f_PS_jack[i][mass_index[i][ik2][ik1]]=(double*) malloc(sizeof(double)*jack_files[i].Njack);

                fread(M_PS_GEVP_jack[i][mass_index[i][ik2][ik1]],   sizeof(double),   jack_files[i].Njack,   jack_files[i].f_M_PS_GEVP );
                fread(f_PS_jack[i][mass_index[i][ik2][ik1]],        sizeof(double),   jack_files[i].Njack,  jack_files[i].f_f_PS );
            }
            }
            
      }
    
}
*/


void  read_files_jack( struct database_file_jack *jack_files, struct header *head,int ***mass_index, struct data_jack *dataJ){
      
      int i,ik1,ik2;
      
      for(i=0;i<ensembles;i++){
            
          
            dataJ[i].M_PS_jack=(double**) malloc(sizeof(double*)*head[i].nk*head[i].nk);
            dataJ[i].f_PS_jack=(double**) malloc(sizeof(double*)*head[i].nk*head[i].nk);

            dataJ[i].M_PS_GEVP_jack=(double**) malloc(sizeof(double*)*head[i].nk*head[i].nk);
            dataJ[i].f_PS_ls_ss_jack=(double**) malloc(sizeof(double*)*head[i].nk*head[i].nk);


            for(ik1=0;ik1<4;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[i].nk;ik2++){
               
                dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=(double*) malloc(sizeof(double)*jack_files[i].Njack);
                dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=(double*) malloc(sizeof(double)*jack_files[i].Njack);
                dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=(double*) malloc(sizeof(double)*jack_files[i].Njack);
                dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=(double*) malloc(sizeof(double)*jack_files[i].Njack);
                
                
                fread(dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]],   sizeof(double),   jack_files[i].Njack,   jack_files[i].f_M_PS );
                fread(dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]],        sizeof(double),   jack_files[i].Njack,  jack_files[i].f_f_PS );
                fread(dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]],   sizeof(double),   jack_files[i].Njack,   jack_files[i].f_M_PS_GEVP );
                fread(dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]],        sizeof(double),   jack_files[i].Njack,  jack_files[i].f_f_PS_ls_ss );
                /*
                if (ik2==0 && ik1==0){
                
               
                 if (i==0){
                  // dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.166249064969,4.35388233333e-4,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0713600049579,1.84703451511e-4,jack_files[i].Njack);
                   
                   dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.166249064969,4.35388233333e-4,jack_files[i].Njack);
                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0713600049579,1.84703451511e-4,jack_files[i].Njack);
                  
                }
               else if (i==1){
                   //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.145705225353, 1.19948022076e-3,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0678279826442,2.25713230671e-4,jack_files[i].Njack);
                   
                   dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.145705225353, 1.19948022076e-3,jack_files[i].Njack);
                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0678279826442,2.25713230671e-4,jack_files[i].Njack);
                }
                else if (i==2){
                  //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.125775584839, 2.94951334969e-4,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0666157700528, 1.34015838616e-4,jack_files[i].Njack);
                     dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.125775584839, 2.94951334969e-4,jack_files[i].Njack);
                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0666157700528, 1.34015838616e-4,jack_files[i].Njack);
                    
                }
                else if (i==3){
                   //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0799,0.0002,jack_files[i].Njack);
                   dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0799,0.0002,jack_files[i].Njack);

                
                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0622,0.0005,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0622,0.0005,jack_files[i].Njack);

                }
               else if (i==4){
                   //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.10429304904, 9.99583163124e-5,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0576087558838, 7.3009468564e-5,jack_files[i].Njack);
                   
                    dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.10429304904, 9.99583163124e-5,jack_files[i].Njack);
                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0576087558838, 7.3009468564e-5,jack_files[i].Njack);
                }
                
                else if (i==5){
                   //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0567528, 0.0000691,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0526446, 0.0000852,jack_files[i].Njack);
                   
                   dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0567528, 0.0000691,jack_files[i].Njack);
                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0526446, 0.0000852,jack_files[i].Njack);
                }
                else if (i==6){
                   //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0799,0.0002,jack_files[i].Njack);
                   dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,4.723748e-02,8.6e-5,jack_files[i].Njack);

                   //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0782,0.0002,jack_files[i].Njack);//lower by supposed FSE of pion mass splitting
                  //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0799,0.02,jack_files[i].Njack);//larger error to eliminate it from the fit

                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0622,0.0005,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,4.496516e-2,1.3653e-4,jack_files[i].Njack);

                }
                
                }*/
            }
            }
      }
 
}


void free_data( struct database_file_jack **jack_files, struct header **head,int *jack_tot, struct data_jack **gJ){
      
    int e,i,j;
    int imoms,imomt,imom0,iks,ikt;
    int iG,im,ik1,ik2;
         

    for(e=0;e<ensembles;e++){
          free((*gJ)[e].KM);
          free((*gJ)[e].Kf);
          free((*gJ)[e].w0);
          free((*gJ)[e].Zp);
          
          for(ik1=0;ik1<4;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
           for(ik2=ik1;ik2<(*head[e]).nk;ik2++){   
                 im=mass_index[e][ik2][ik1];
                 free((*gJ[e]).M_PS_jack[im]);     
                 free((*gJ[e]).f_PS_jack[im]);     
                 free((*gJ[e]).M_PS_GEVP_jack[im]);     
                 free((*gJ[e]).f_PS_ls_ss_jack[im]);   
           }}
                 free((*gJ[e]).M_PS_jack);     
                 free((*gJ[e]).f_PS_jack);     
                 free((*gJ[e]).M_PS_GEVP_jack);     
                 free((*gJ[e]).f_PS_ls_ss_jack);   
           
    }
    
    free((*gJ));   
     ///////////////////////free header///////////////
    for (e=0;e<ensembles_reph;e++){
        free((*head)[e].k);
        for (i=0;i<(*head)[e].nmoms;i++)
            free((*head)[e].mom[i]);
        free((*head)[e].mom);
    }
    free((*head));
    for (e=0;e<ensembles_reph;e++){

             fclose((*jack_files)[e].f_M_PS);
             fclose((*jack_files)[e].f_f_PS);
             fclose((*jack_files)[e].f_M_PS_GEVP);
             fclose((*jack_files)[e].f_f_PS_ls_ss);
    }
    free((*jack_files)); 
  
}
//P[0]=B,  P[1]=f, P[2]=P1, P[3]=P2
//x[0]=m_l, x[1]=w0,  x[2]=KM
double *fun_Mw2_k( int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double *Mw2_k=(double*) calloc(Npar,sizeof(double));
    double Bw=P[0], fw=P[1], P1=P[2], P2=P[3];
    double ml=x[0], w0=x[1], KM=x[2];
    
    xi=2*Bw*ml*w0/(16.*pi*pi*fw*fw);
    
    Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2;
    Mw2*=2*Bw*w0*ml*KM;
    
    Mw2_k[0]=Mw2/Bw+ 2*Bw*w0*ml*KM* ( log(xi)+ 1+ P1  )*(xi/Bw);
    Mw2_k[1]= 2*Bw*w0*ml*KM* ( log(xi)+ 1+ P1  )*(-2*xi/fw);
    Mw2_k[2]=2*Bw*w0*ml*KM*(xi);
    Mw2_k[3]=2*Bw*w0*ml*KM*(1/(w0*w0));
    
    return Mw2_k;
    
}




void mud(double xi,int Npar,double *P, double *f, double *df){
    
    double h=0.001;
  
    
    *f=m_over_f_xi( xi, Npar,P);
    
    xi=xi-2.*h;
    *df=m_over_f_xi( xi, Npar,P);
    
    xi=xi+h;
    *df-=8*m_over_f_xi( xi, Npar,P);
    
    xi=xi+2*h;
    *df+=8*m_over_f_xi( xi, Npar,P);
    
    xi=xi+h;
    *df-=m_over_f_xi( xi, Npar,P);
    
    xi=xi-2*h;
    *df/=(12.*h);
    
}



void  create_fake_distribution(const char *jackboot,double **w0A,double **w0B,double **w0C,double **ZpA,double **ZpB,double **ZpC,int jack_tot, const char *scaletype, const char *Mtype){
    
      if (strcmp(scaletype,"w0")==0){
          //(*w0A)=fake_sampling(jackboot,1.8381,  0.0037,jack_tot,123);  // petros-from jacob 31/07/20 fit   in M_PS^2/f_PS^2
          //(*w0B)=fake_sampling(jackboot,2.1316 ,0.0024 ,jack_tot,1234);// petros-from jacob 31/07/20 fit  M_PS^2/f_PS^2
          //(*w0C)=fake_sampling(jackboot,2.5039, 0.0017,jack_tot,12345);//  petros-from jacob 31/07/20 fit  M_PS^2/f_PS^2
          (*w0A)=fake_sampling(jackboot,1.8353,  0.0035,jack_tot,123);  // petros-from jacob 31/07/20 fit   in M_PS^2/f_PS^2
          (*w0B)=fake_sampling(jackboot,2.1300 ,0.0017 ,jack_tot,1234);// petros-from jacob 31/07/20 fit  M_PS^2/f_PS^2
          (*w0C)=fake_sampling(jackboot,2.5039, 0.0017,jack_tot,12345);//  petros-from jacob 31/07/20 fit  M_PS^2/f_PS^2
          
      }
      else if(strcmp(scaletype,"t0/w0")==0){
          (*w0A)=fake_sampling(jackboot,1.33582,0.001038,jack_tot,123); //t/w0
          (*w0B)=fake_sampling(jackboot,1.52764 ,0.000342,jack_tot,1234);// t/w0
          (*w0C)=fake_sampling(jackboot,1.77671,0.00048,jack_tot,12345);
     }
     
     if(strcmp(Mtype,"M1")==0){
          /*(*ZpA)=fake_sampling(jackboot, 0.474,0.002,jack_tot,321);//M1 2Gev  Enrico-Matteo  19/05/2020  constant
          (*ZpB)=fake_sampling(jackboot, 0.482,0.0030,jack_tot,3214);//M1 2Gev Enrico-Matteo  19/05/2020  constant
          (*ZpC)=fake_sampling(jackboot, 0.493,0.0030,jack_tot,32145);//M1 2Gev  Enrico-Matteo  19/05/2020  constant  */
          //(*ZpA)=fake_sampling(jackboot,0.471,0.002,jack_tot,321);//M2 2Gev  Enrico-Matteo  28/05/2020  constant
          //(*ZpB)=fake_sampling(jackboot,0.479,0.002,jack_tot,3214);//M2 2Gev Enrico-Matteo  28/05/2020  constant
         //(*ZpC)=fake_sampling(jackboot,0.4915,0.002,jack_tot,32145);//M2 2Gev  Enrico-Matteo  28/05/2020  constant
         (*ZpA)=fake_sampling(jackboot,0.474,0.002,jack_tot,321);//M2 2Gev  Enrico-Matteo  28/05/2020  constant
         (*ZpB)=fake_sampling(jackboot,0.479,0.003,jack_tot,3214);//M2 2Gev Enrico-Matteo  28/05/2020  constant
         (*ZpC)=fake_sampling(jackboot,0.489,0.002,jack_tot,32145);//M2 2Gev  Enrico-Matteo  28/05/2020  constant
         
         
     }
     if(strcmp(Mtype,"M1a")==0){
         (*ZpA)=fake_sampling(jackboot,0.4748,0.0024,jack_tot,321);//M2 2Gev  petros report Oct/2020  
         (*ZpB)=fake_sampling(jackboot,0.4784,0.0034,jack_tot,3214);//M2 2Gev petros report Oct/2020
         (*ZpC)=fake_sampling(jackboot,0.4871,0.0026,jack_tot,32145);//M2 2Gev  petros report Oct/2020
     }
     if(strcmp(Mtype,"M1b")==0){
         (*ZpA)=fake_sampling(jackboot,0.4748,0.0026,jack_tot,321);//M2 2Gev  petros report Oct/2020  
         (*ZpB)=fake_sampling(jackboot,0.4775,0.0027,jack_tot,3214);//M2 2Gev petros report Oct/2020
         (*ZpC)=fake_sampling(jackboot,0.4873,0.0027,jack_tot,32145);//M2 2Gev  petros report Oct/2020
     }
     else if(strcmp(Mtype,"M2a")==0){ //(9)??
         double cp=1.03496; 
         /*(*ZpA)=fake_sampling(jackboot,cp* 0.4744,cp*0.0016,jack_tot,321);//M2 2Gev  Enrico-Matteo  19/05/2020  constant
         (*ZpB)=fake_sampling(jackboot,cp* 0.4777,cp*0.0020,jack_tot,3214);//M2 2Gev Enrico-Matteo  19/05/2020  constant
         (*ZpC)=fake_sampling(jackboot,cp* 0.4821,cp*0.0017,jack_tot,32145);//M2 2Gev  Enrico-Matteo  19/05/2020  constant  */
         //(*ZpA)=fake_sampling(jackboot,0.500,0.001,jack_tot,321);//M2 2Gev  Enrico-Matteo  28/05/2020  constant
         //(*ZpB)=fake_sampling(jackboot,0.501,0.002,jack_tot,3214);//M2 2Gev Enrico-Matteo  28/05/2020  constant
         //(*ZpC)=fake_sampling(jackboot,0.503,0.002,jack_tot,32145);//M2 2Gev  Enrico-Matteo  28/05/2020  constant
         //(*ZpA)=fake_sampling(jackboot,0.496,0.002,jack_tot,321);//M2 2Gev  Enrico-Matteo  28/05/2020  constant
         //(*ZpB)=fake_sampling(jackboot,0.495,0.002,jack_tot,3214);//M2 2Gev Enrico-Matteo  28/05/2020  constant
         //(*ZpC)=fake_sampling(jackboot,0.499,0.002,jack_tot,32145);//M2 2Gev  Enrico-Matteo  28/05/2020  constant
         (*ZpA)=fake_sampling(jackboot,0.5050,0.0024,jack_tot,321);//M2 2Gev  petros report Oct/2020  
         (*ZpB)=fake_sampling(jackboot,0.5013,0.0026,jack_tot,3214);//M2 2Gev petros report Oct/2020
         (*ZpC)=fake_sampling(jackboot,0.5024,0.0023,jack_tot,32145);//M2 2Gev  petros report Oct/2020
     }
     else if(strcmp(Mtype,"M2b")==0){//(9)??
         double cp= 0.972164;
         /*(*ZpA)=fake_sampling(jackboot,cp* 0.5140,cp*0.0015,jack_tot,321);//M2 2Gev  Enrico-Matteo  19/05/2020  constant
         (*ZpB)=fake_sampling(jackboot,cp* 0.514,cp*0.0020,jack_tot,3214);//M2 2Gev Enrico-Matteo  19/05/2020  constant
         (*ZpC)=fake_sampling(jackboot,cp* 0.5170,cp*0.0016,jack_tot,32145);//M2 2Gev  Enrico-Matteo  19/05/2020  constant
         (*ZpA)=fake_sampling(jackboot,0.508,0.002,jack_tot,321);//M2 2Gev  Enrico-Matteo  19/05/2020  constant
         (*ZpB)=fake_sampling(jackboot,0.505,0.002,jack_tot,3214);//M2 2Gev Enrico-Matteo  19/05/2020  constant
         (*ZpC)=fake_sampling(jackboot,0.505,0.002,jack_tot,32145);//M2 2Gev  Enrico-Matteo  19/05/2020  constant*/
         //(*ZpA)=fake_sampling(jackboot,0.508,0.001,jack_tot,321);//M2 2Gev  Enrico-Matteo  28/05/2020  constant
         //(*ZpB)=fake_sampling(jackboot,0.5054,0.002,jack_tot,3214);//M2 2Gev Enrico-Matteo  28/05/2020  constant
         //(*ZpC)=fake_sampling(jackboot,0.506,0.002,jack_tot,32145);//M2 2Gev  Enrico-Matteo  28/05/2020  constant
         //(*ZpA)=fake_sampling(jackboot,0.502,0.001,jack_tot,321);//M2 2Gev  Enrico-Matteo  28/05/2020  constant
         //(*ZpB)=fake_sampling(jackboot,0.499,0.002,jack_tot,3214);//M2 2Gev Enrico-Matteo  28/05/2020  constant
         //(*ZpC)=fake_sampling(jackboot,0.501,0.002,jack_tot,32145);//M2 2Gev  Enrico-Matteo  28/05/2020  constant
         (*ZpA)=fake_sampling(jackboot,0.5119,0.0024,jack_tot,321);//M2 2Gev  petros report Oct/2020  
         (*ZpB)=fake_sampling(jackboot,0.5075,0.0023,jack_tot,3214);//M2 2Gev petros report Oct/2020
         (*ZpC)=fake_sampling(jackboot,0.5060,0.0024,jack_tot,32145);//M2 2Gev  petros report Oct/2020
     }
      //w0A=fake_sampling(jackboot,1.8346689, 0.005178046,*jack_tot);
      //w0A=fake_sampling(jackboot,1.83005,3.48173757101e-3,*jack_tot,rand());// MG fit in M_PS^2
      //(*w0A)=fake_sampling(jackboot,1.8355,0.0042,jack_tot,123); // PD fit
      
      
      //w0B=fake_sampling(jackboot,2.1330729,0.00468807,*jack_tot,rand());// MG fit in M_PS^2
      //(*w0B)=fake_sampling(jackboot,2.1347,0.0047,jack_tot,1234);// MG fit in M_PS^2
      //(*w0B)=fake_sampling(jackboot,2.12650,0.0023,jack_tot,1234);// MG fit 26/11/2119  in M_PS^2/f_PS^2
      //(
      
      //(*w0C)=fake_sampling(jackboot,2.49879971,0.0034,jack_tot,12345);// MG fit in M_PS^2
      //(*w0C)=fake_sampling(jackboot,1.77671,0.00048,jack_tot,12345);//  t/w0

      //  ZpA=fake_sampling(jackboot,0.459,0.005,*jack_tot);//M1 2Gev  uncorrected
      //  ZpA=fake_sampling(jackboot,0.527,0.004,*jack_tot);//M2 2Gev  uncorrected
      //  ZpA=fake_sampling(jackboot, 0.485,0.005,*jack_tot);//M1 2Gev  29/3/2019  a^2g*^2
      //ZpA=fake_sampling(jackboot,0.502,0.004,*jack_tot);//M2 2Gev  29/3/2019  a^2g*^2
      //  ZpA=fake_sampling(jackboot, 0.482,0.005,*jack_tot);//M1 2Gev  29/3/2019  a^inf g_0^2
      //  ZpA=fake_sampling(jackboot, 0.530,0.004,*jack_tot);//M2 2Gev  29/3/2019  a^inf g_0^2
      //ZpA=fake_sampling(jackboot,0.471,0.005,*jack_tot);//M1 2Gev  fiorenza  
      //ZpA=fake_sampling(jackboot,0.491,0.004,*jack_tot);//M2a 2Gev  fiorenza  
      //ZpA=fake_sampling(jackboot, 0.508,0.003,*jack_tot,rand());//M2b 2Gev  fiorenza  
      // (*ZpA)=fake_sampling(jackboot, 0.4770,0.0045,jack_tot,321);//M1 2Gev  Petros  2/9/2019
       //(*ZpA)=fake_sampling(jackboot, 0.4628,0.0052,jack_tot,321);//M1 2Gev  Petros  15/12/2019  quadratic
       //(*ZpA)=fake_sampling(jackboot, 0.478,0.004,jack_tot,321);//M1 2Gev  Petros  15/12/2019   constant

      
      //  ZpB=fake_sampling(jackboot,0.471,0.007,*jack_tot);//M1 2Gev  uncorrected
      //  ZpB=fake_sampling(jackboot,0.502,0.005,*jack_tot);//M2 2Gev  uncorrected
      //  ZpB=fake_sampling(jackboot,0.484,0.007,*jack_tot);//M1 2Gev  29/3/2019  a^2g*^2
      //ZpB=fake_sampling(jackboot,0.491,0.005,*jack_tot);//M2 2Gev  29/3/2019  a^2g*^2
      //  ZpB=fake_sampling( jackboot,0.485,0.007,*jack_tot);//M1 2Gev  29/3/2019  a^inf g_0^2
      //  ZpB=fake_sampling( jackboot,0.509,0.005,*jack_tot);//M2 2Gev  29/3/2019  a^inf g_0^2
     // ZpB=fake_sampling(jackboot,0.476,0.008,*jack_tot);//M1 2Gev  fiorenza
      //ZpB=fake_sampling(jackboot,0.486,0.005,*jack_tot);//M2a 2Gev  fiorenza
      //ZpB=fake_sampling(jackboot,0.500,0.003,*jack_tot,rand());//M2b 2Gev  fiorenza
      //(*ZpB)=fake_sampling(jackboot, 0.4860,0.0070,jack_tot,3214);//M1 2Gev  Petros  2/9/2019
      //(*ZpB)=fake_sampling(jackboot, 0.4780,0.0070,jack_tot,3214);//M1 2Gev  Petros  15/12/2019 quadratic 
      //(*ZpB)=fake_sampling(jackboot, 0.487,0.0040,jack_tot,3214);//M1 2Gev  Petros  15/12/2019  constant
      

     //  ZpC=fake_sampling(jackboot,0.492333,0.003,*jack_tot);// extrapolated from M2 2Gev  16/2/2019  linear
      //ZpC=fake_sampling(jackboot,0.497,0.003,*jack_tot);// extrapolated from M2 2Gev  16/2/2019  RF
      //(*ZpC)=fake_sampling(jackboot, 0.4860,0.0030,jack_tot,32145);//M1 2Gev  Petros  15/12/2019 quadratic
      //(*ZpC)=fake_sampling(jackboot, 0.484,0.0060,jack_tot,32145);//M1 2Gev  Petros  15/12/2019 constant
      
       int seedw0=213;
      result.w0fm=fake_sampling(jackboot,v_w0fm,err_w0fm,jack_tot,seedw0);
      result.w0MeV=fake_sampling(jackboot,v_w0MeV,err_w0fm/197.326963 ,jack_tot,seedw0);
      result.MpiMeV=fake_sampling(jackboot,v_MpiMeV,err_MpiMeV,jack_tot,2134);
      result.MKMeV=fake_sampling(jackboot,v_MKMeV,err_MKMeV,jack_tot,21345);
      result.MDMeV=fake_sampling(jackboot,v_MDMeV,err_MDMeV,jack_tot,321);
      result.MDsMeV=fake_sampling(jackboot,v_MDsMeV,err_MDsMeV,jack_tot,3124);
      result.fpiMeV_exp=fake_sampling(jackboot,v_fpiMeV_exp,err_fpiMeV_exp,jack_tot,31245);
      result.MOmegaMeV=fake_sampling(jackboot,v_MOmegaMeV,err_MOmegaMeV,jack_tot,111);
      
      result.mlw=(double*) malloc(sizeof(double)*jack_tot);
      result.Bw=(double*) malloc(sizeof(double)*jack_tot);
      result.fw=(double*) malloc(sizeof(double)*jack_tot);
      //result.w0fm= (double*) malloc(sizeof(double)*jack_tot);
      //result.w0MeV= (double*) malloc(sizeof(double)*jack_tot);

    
}



void init_Z( struct database_file_jack *jack_files, struct header *head,int jack_tot, struct data_jack **gJ, const char *scaletype, const char *Mtype){
      int j;
      double *w0A,*w0B,*w0C, *ZpA,*ZpB,*ZpC;
      create_fake_distribution(jack_files[0].sampling, &w0A, &w0B, &w0C, &ZpA, &ZpB, &ZpC,jack_tot,scaletype,Mtype);
      for(j=0;j<jack_tot;j++){
      if (ensembles>0){
        (*gJ)[0].w0[j]=w0A[j];
        (*gJ)[0].Zp[j]=ZpA[j];
      }
      if (ensembles>1){
        (*gJ)[1].w0[j]=w0A[j];
        (*gJ)[1].Zp[j]=ZpA[j];        
      }
      if (ensembles>2){
        (*gJ)[2].w0[j]=w0A[j];
        (*gJ)[2].Zp[j]=ZpA[j];        
      }
      if (ensembles>3){
        (*gJ)[3].w0[j]=w0A[j];
        (*gJ)[3].Zp[j]=ZpA[j];
      }
      if (ensembles>4){
        (*gJ)[4].w0[j]=w0B[j];
        (*gJ)[4].Zp[j]=ZpB[j];
      }
      if (ensembles>5){
        (*gJ)[5].w0[j]=w0B[j];
        (*gJ)[5].Zp[j]=ZpB[j];
      }
      if (ensembles>6){
        (*gJ)[6].w0[j]=w0B[j];
        (*gJ)[6].Zp[j]=ZpB[j];
      }
      
      if (ensembles>7){
        (*gJ)[7].w0[j]=w0C[j];
        (*gJ)[7].Zp[j]=ZpC[j];
      }
      if (ensembles>8){
        (*gJ)[8].w0[j]=w0B[j];
        (*gJ)[8].Zp[j]=ZpB[j];
      } 
          
    }
    free(w0A);free(w0B);free(w0C); free(ZpA);free(ZpB);free(ZpC);
    
}



struct data_jack *create_generalised_boot( struct database_file_jack *jack_files, struct header *head,int *jack_tot,int ***mass_index, struct data_jack *dJ){
      int j,e,e1,ik1,ik2,counter;
      int im;
      double ***M_PS_GEVP_jack_tot;
      struct data_jack *gJ;
      double ***omega;
      omega=init_omega_jacob();
      
      gJ=(struct data_jack *) malloc (sizeof(struct data_jack )*ensembles);
      
      *jack_tot=0;
      for(e=0;e<ensembles-1;e++){
          error(jack_files[e].Njack!=jack_files[e+1].Njack,1,"create_generalised_boot","bootstrap of the file %d has different number",e+1 );
      }
      *jack_tot=jack_files[0].Njack;

      for(e=0;e<ensembles;e++){
          gJ[e].M_PS_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          gJ[e].f_PS_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          
          gJ[e].M_PS_GEVP_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          gJ[e].f_PS_ls_ss_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          
          gJ[e].KM=(double*) malloc(sizeof(double)*head[e].nk*head[e].nk);
          gJ[e].Kf=(double*) malloc(sizeof(double)*head[e].nk*head[e].nk);
          gJ[e].w0=(double*) calloc(*jack_tot,sizeof(double));
          gJ[e].Zp=(double*) calloc(*jack_tot,sizeof(double));
            for(ik1=0;ik1<4;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[e].nk;ik2++){   
                 im=mass_index[e][ik2][ik1];
                 gJ[e].M_PS_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 gJ[e].f_PS_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 
                 gJ[e].M_PS_GEVP_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 gJ[e].f_PS_ls_ss_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 
                 
                 gJ[e].KM[im]=dJ[e].KM[im];
                 gJ[e].Kf[im]=dJ[e].Kf[im];
                 
                 
                 for(j=0;j<(*jack_tot);j++){
                           
                           gJ[e].M_PS_jack[im][j]=dJ[e].M_PS_jack[  im ][j];
                           gJ[e].f_PS_jack[im][j]=dJ[e].f_PS_jack[  im ][j];
                           
                           gJ[e].M_PS_GEVP_jack[im][j]=dJ[e].M_PS_GEVP_jack[  im ][j];
                           gJ[e].f_PS_ls_ss_jack[im][j]=dJ[e].f_PS_ls_ss_jack[  im ][j];
                           
                 }
                 
            }}
            gJ[e].M_Omega_jack=(double**) malloc(sizeof(double*)*3);
            gJ[e].M_Omega_jack[0]= fake_sampling(jack_files[0].sampling, omega[e][0][0],omega[e][0][1],*jack_tot,e);        
            gJ[e].M_Omega_jack[1]= fake_sampling(jack_files[0].sampling, omega[e][1][0],omega[e][1][1],*jack_tot,e);        
            gJ[e].M_Omega_jack[2]= fake_sampling(jack_files[0].sampling, omega[e][2][0],omega[e][2][1],*jack_tot,e);        
            
      }

      
/*
      double *w0A,*w0B,*w0C, *ZpA,*ZpB,*ZpC;
      create_fake_distribution(jack_files[0].sampling, &w0A, &w0B, &w0C, &ZpA, &ZpB, &ZpC,*jack_tot);
      
      
      if (ensembles>0){
        gJ[0].w0=w0A;
        gJ[0].Zp=ZpA;
      }
      if (ensembles>1){
        gJ[1].w0=w0A;
        gJ[1].Zp=ZpA;        
      }
      if (ensembles>2){
        gJ[2].w0=w0A;
        gJ[2].Zp=ZpA;        
      }
      if (ensembles>3){
        gJ[3].w0=w0A;
        gJ[3].Zp=ZpA;
      }
      if (ensembles>4){
        gJ[4].w0=w0B;
        gJ[4].Zp=ZpB;
      }
      if (ensembles>5){
        gJ[5].w0=w0B;
        gJ[5].Zp=ZpB;
      }
      if (ensembles>6){
        gJ[6].w0=w0B;
        gJ[6].Zp=ZpB;
      }
      if (ensembles>7){
        gJ[7].w0=w0C;
        gJ[7].Zp=ZpC;
      }
      if (ensembles>8){
        gJ[8].w0=w0B;
        gJ[8].Zp=ZpB;
      }
      */
     // free(w0A); free(w0B); free(w0C); free(ZpA); free(ZpB); free(ZpC);
      
      
      for(e=0;e<ensembles;e++){
            for(ik1=0;ik1<2;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[e].nk;ik2++){
               free(dJ[e].M_PS_jack[ mass_index[e][ik2][ik1] ]);
               free(dJ[e].f_PS_jack[ mass_index[e][ik2][ik1] ]);
               free(dJ[e].M_PS_GEVP_jack[ mass_index[e][ik2][ik1] ]);
               free(dJ[e].f_PS_ls_ss_jack[ mass_index[e][ik2][ik1] ]);

            }
            }
            free(dJ[e].KM);free(dJ[e].Kf);
            free(dJ[e].M_PS_jack); free(dJ[e].M_PS_GEVP_jack); free(dJ[e].f_PS_jack); free(dJ[e].f_PS_ls_ss_jack);
            free(omega[e][0]);free(omega[e][1]);free(omega[e][2]);
      }
      free(dJ);          free(omega);  
      
     return gJ;
}
struct data_jack *create_generalised_jack( struct database_file_jack *jack_files, struct header *head,int *jack_tot,int ***mass_index, struct data_jack *dJ){
      int j,e,e1,ik1,ik2,counter;
      int im;
      double ***M_PS_GEVP_jack_tot;
      struct data_jack *gJ;
      double ***omega;
      omega=init_omega_jacob();
 
      gJ=(struct data_jack *) malloc (sizeof(struct data_jack )*ensembles);
      
      *jack_tot=0;
      for(e=0;e<ensembles;e++){
          *jack_tot+=jack_files[e].Njack;
      }
      *jack_tot=*jack_tot-ensembles+1;

      for(e=0;e<ensembles;e++){
          gJ[e].M_PS_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          gJ[e].f_PS_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          
          gJ[e].M_PS_GEVP_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          gJ[e].f_PS_ls_ss_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          
          
          gJ[e].KM=(double*) malloc(sizeof(double)*head[e].nk*head[e].nk);
          gJ[e].Kf=(double*) malloc(sizeof(double)*head[e].nk*head[e].nk);
          gJ[e].w0=(double*) malloc((*jack_tot)*sizeof(double));
          gJ[e].Zp=(double*) malloc((*jack_tot)*sizeof(double));
            for(ik1=0;ik1<4;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[e].nk;ik2++){   
                 im=mass_index[e][ik2][ik1];
                 gJ[e].M_PS_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 gJ[e].f_PS_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 
                 gJ[e].M_PS_GEVP_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 gJ[e].f_PS_ls_ss_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 
                 gJ[e].KM[im]=dJ[e].KM[im];
                 gJ[e].Kf[im]=dJ[e].Kf[im];
                 
                 counter=0;
                 for(e1=0;e1<ensembles;e1++){
                      for(j=0;j<(jack_files[e1].Njack-1);j++){
                           if (e==e1){
                           gJ[e].M_PS_jack[im][j+counter]=dJ[e].M_PS_jack[  im ][j];
                           gJ[e].f_PS_jack[im][j+counter]=dJ[e].f_PS_jack[  im ][j];
                           
                           gJ[e].M_PS_GEVP_jack[im][j+counter]=dJ[e].M_PS_GEVP_jack[  im ][j];
                           gJ[e].f_PS_ls_ss_jack[im][j+counter]=dJ[e].f_PS_ls_ss_jack[  im ][j];
                           }
                           else{ 
                           gJ[e].M_PS_jack[im][j+counter]=dJ[e].M_PS_jack[im][ jack_files[e].Njack-1 ];
                           gJ[e].f_PS_jack[im][j+counter]=dJ[e].f_PS_jack[im][ jack_files[e].Njack-1 ];   
                           
                           gJ[e].M_PS_GEVP_jack[im][j+counter]=dJ[e].M_PS_GEVP_jack[im][ jack_files[e].Njack-1 ];
                           gJ[e].f_PS_ls_ss_jack[im][j+counter]=dJ[e].f_PS_ls_ss_jack[im][ jack_files[e].Njack-1 ];    
                           }
                      }
                      counter+=jack_files[e1].Njack-1;
                 }
                 gJ[e].M_PS_jack[im][*jack_tot-1]=dJ[e].M_PS_jack[im][jack_files[e].Njack-1];
                 gJ[e].f_PS_jack[im][*jack_tot-1]=dJ[e].f_PS_jack[im][jack_files[e].Njack-1];
                 gJ[e].M_PS_GEVP_jack[im][*jack_tot-1]=dJ[e].M_PS_GEVP_jack[im][jack_files[e].Njack-1];
                 gJ[e].f_PS_ls_ss_jack[im][*jack_tot-1]=dJ[e].f_PS_ls_ss_jack[im][jack_files[e].Njack-1];
                 
            }}

            counter=0;
            gJ[e].M_Omega_jack=(double**) malloc(sizeof(double*)*3);
            gJ[e].M_Omega_jack[0]= fake_sampling(jack_files[0].sampling, omega[e][0][0],omega[e][0][1],*jack_tot,e);        
            gJ[e].M_Omega_jack[1]= fake_sampling(jack_files[0].sampling, omega[e][1][0],omega[e][1][1],*jack_tot,e);        
            gJ[e].M_Omega_jack[2]= fake_sampling(jack_files[0].sampling, omega[e][2][0],omega[e][2][1],*jack_tot,e);        

            /*gJ[e].Zp=(double*) malloc(sizeof(double)*(*jack_tot));//jack_tot is a pointer here
            gJ[e].w0=(double*) malloc(sizeof(double)*(*jack_tot));//jack_tot is a pointer here
            for(e1=0;e1<ensembles;e1++){
                for(j=0;j<(jack_files[e1].Njack-1);j++){
                    if (e==e1){
                        gJ[e].Zp[j+counter]=dJ[e].Zp[j];
                        gJ[e].w0[j+counter]=dJ[e].w0[j];

                    }
                    else{ 
                        gJ[e].Zp[j+counter]=dJ[e].Zp[ jack_files[e].Njack-1 ];//jack_tot is a pointer here
                        gJ[e].w0[j+counter]=dJ[e].w0[ jack_files[e].Njack-1 ];
                           
                    }
                }
                counter+=jack_files[e1].Njack-1;
            } 
            gJ[e].Zp[*jack_tot-1]=dJ[e].Zp[ jack_files[e].Njack-1 ];
            gJ[e].w0[*jack_tot-1]=dJ[e].w0[ jack_files[e].Njack-1 ];*/
            
      }
      /*double *w0A,*w0B,*w0C, *ZpA,*ZpB,*ZpC;
      create_fake_distribution(jack_files[0].sampling, &w0A, &w0B, &w0C, &ZpA, &ZpB, &ZpC,*jack_tot);
      
      if (ensembles>0){
        gJ[0].w0=w0A;
        gJ[0].Zp=ZpA;
      }
      if (ensembles>1){
        gJ[1].w0=w0A;
        gJ[1].Zp=ZpA;        
      }
      if (ensembles>2){
        gJ[2].w0=w0A;
        gJ[2].Zp=ZpA;        
      }
      if (ensembles>3){
        gJ[3].w0=w0A;
        gJ[3].Zp=ZpA;
      }
      if (ensembles>4){
        gJ[4].w0=w0B;
        gJ[4].Zp=ZpB;
      }
      if (ensembles>5){
        gJ[5].w0=w0B;
        gJ[5].Zp=ZpB;
      }
      if (ensembles>6){
        gJ[6].w0=w0B;
        gJ[6].Zp=ZpB;
      }
      
      if (ensembles>7){
        gJ[7].w0=w0C;
        gJ[7].Zp=ZpC;
      }
      if (ensembles>8){
        gJ[8].w0=w0B;
        gJ[8].Zp=ZpB;
      } */
      for(e=0;e<ensembles;e++){
            for(ik1=0;ik1<2;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[e].nk;ik2++){
               free(dJ[e].M_PS_jack[ mass_index[e][ik2][ik1] ]);
               free(dJ[e].f_PS_jack[ mass_index[e][ik2][ik1] ]);
               free(dJ[e].M_PS_GEVP_jack[ mass_index[e][ik2][ik1] ]);
               free(dJ[e].f_PS_ls_ss_jack[ mass_index[e][ik2][ik1] ]);

            }
            }
            free(dJ[e].KM);free(dJ[e].Kf);
            free(dJ[e].M_PS_jack); free(dJ[e].M_PS_GEVP_jack); free(dJ[e].f_PS_jack); free(dJ[e].f_PS_ls_ss_jack);
            free(omega[e][0]);free(omega[e][1]);free(omega[e][2]);
      }
      free(dJ);            free(omega);
      
     return gJ;
}

//table 3 of arXiv:hep-lat/0503014  fit: log R= A+B*M+C*L // R=(M_L-M_inf)/M_inf //  K*M_inf=M_L
void KM_FSE(struct database_file_jack *jack_files, struct header *head, struct data_jack *dJ){
    double RM,L,M;
    int ik1,ik2,im,e;
    
    for (e=0;e<ensembles;e++){
            dJ[e].KM=(double*) malloc(sizeof(double)*head[e].nk*head[e].nk);
            for(ik1=0;ik1<2;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[e].nk;ik2++){   
                im=mass_index[e][ik2][ik1];
                L=jack_files[e].a*head[e].l1;
                M=dJ[e].M_PS_GEVP_jack[im][jack_files[e].Njack-1]*197.326963/jack_files[e].a;
                
                RM=exp(3.81729-0.0130342*M-2.1714*L);
                dJ[e].KM[im]=RM+1;
                if (ik1==0 && ik2==0)printf("%f   %f   im=%d\n",dJ[e].M_PS_GEVP_jack[im][jack_files[e].Njack-1],dJ[e].KM[im], im );
            }}
            
    }
    
}

void Kf_FSE(struct database_file_jack *jack_files, struct header *head, struct data_jack *dJ){
    double Rf,L,M;
    int ik1,ik2,im,e;
    
    for (e=0;e<ensembles;e++){
            dJ[e].Kf=(double*) malloc(sizeof(double)*head[e].nk*head[e].nk);
            for(ik1=0;ik1<2;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[e].nk;ik2++){   
                im=mass_index[e][ik2][ik1];
                L=jack_files[e].a*head[e].l1;
                M=dJ[e].M_PS_GEVP_jack[im][jack_files[e].Njack-1]*197.326963/jack_files[e].a;
                
                Rf=exp(4.58982-0.0138032*M-2.013*L);
                dJ[e].Kf[im]=-Rf+1;
            }}
    }
    
}

void print_chiral_continuum_fit(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, const char *AV,const char *namefile, struct header *head ,struct data_jack *gJ){
    
    int i,j;
    char name[NAMESIZE];
    FILE *fc=NULL,*fcA=NULL,*fcB=NULL,*fcC=NULL;
    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum.txt",argv[2],namefile);
    fc=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_A.txt",argv[2],namefile);
    fcA=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_B.txt",argv[2],namefile);
    fcB=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_C.txt",argv[2],namefile);
    fcC=open_file(name,"w+");
    int N=fit_info.N,n;
    
    double h=0.0003;
    double **tif=fit_to_tif(fit_info.Npar,jack_tot,fit_out.P);   
    double *r=(double*) malloc(sizeof(double)*jack_tot);
    double *m;
    double **x=(double**) malloc(sizeof(double*)*jack_tot);
    for(j=0;j<jack_tot;j++){
        x[j]=(double*) malloc(sizeof(double)*fit_info.Nvar);
        //mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
        x[j][0]=0;//mlw
        x[j][1]=1e+5;//r0
        x[j][2]=1e+8;//Mpi2
        x[j][3]=1e+8;//fpi
        x[j][4]=1e+12;//L such that L/w=1e+4
        
    }
    for (i=0;i<100;i++){
        fprintf(fc,"%g \t",((double)i)*h);
        for (n=0;n<N;n++){
            for(j=0;j<jack_tot;j++){
                x[j][0]=((double)i)*h;//=1e+10;//xG
                r[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
            }
            m=mean_and_error_jack_biased(jack_tot,r);
            fprintf(fc," %g  %g \t",m[0],m[1]);
            free(m);
        }
        fprintf(fc,"\n");
    }
    for (i=2;i<100;i++){
        fprintf(fcA,"%g \t",((double)i)*h);
        for (n=0;n<N;n++){
            for(j=0;j<jack_tot;j++){
                x[j][1]=gJ[0].w0[j];
                x[j][0]=((double)i)*h;//=1e+10;//xG
                r[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
            }
            m=mean_and_error_jack_biased(jack_tot,r);
            fprintf(fcA," %g  %g \t",m[0],m[1]);
            free(m);
        }
        fprintf(fcA,"%g \n",gJ[0].w0[jack_tot-1]);

        fprintf(fcB,"%g \t",((double)i)*h);
        for (n=0;n<N;n++){
            for(j=0;j<jack_tot;j++){
                x[j][1]=gJ[4].w0[j];
                x[j][0]=((double)i)*h;//=1e+10;//xG
                r[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
            }
            m=mean_and_error_jack_biased(jack_tot,r);
            fprintf(fcB," %g  %g \t",m[0],m[1]);
            free(m);
        }
        fprintf(fcB,"\n");

       fprintf(fcC,"%g \t",((double)i)*h);
        for (n=0;n<N;n++){
            for(j=0;j<jack_tot;j++){
                x[j][1]=gJ[7].w0[j];
                x[j][2]=((double)i)*h;//=1e+10;//xG
                r[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
            }
            m=mean_and_error_jack_biased(jack_tot,r);
            fprintf(fcC," %g  %g \t",m[0],m[1]);
            free(m);
        }
        fprintf(fcC,"\n");

    }
    fclose(fc);
   

    free(r);
    free_2(jack_tot,tif);
    free_2(jack_tot,x);
    fclose(fcA);    fclose(fcB);    fclose(fcC);
    
}



void  print_fit_info(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, struct result_jack &r1, struct data_jack *grephJ, struct header *head , const char *AV,const char *namefile){
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
    
    printf("chi2=%g\n",fit_out.chi2[jack_tot-1] );
    chi2m=mean_and_error( argv[1],jack_tot,fit_out.chi2);
    //for(j=0;j<jack_tot;j++)
    //       error(fabs(fit_out.chi2[j]-chi2m[0])/chi2m[1]>3.0/sqrt(jack_tot-2), 44,"print fit info","chi2 of jacknife=%d  is %g    while the average is %g +- %g",j,fit_out.chi2[j],chi2m[0],chi2m[1] );
      
    /*
    for(j=0;j<jack_tot;j++){
        fk=der_fun_Nf_h(fit_info.N,  fit_info.Nvar, phys_point[j], fit_info.Npar,tif[j],  fit_info.function,  0.00001);
        for (i=0;i<fit_info.Npar;i++)
            tmp2[i][j]=fk[i]*fit[i][j];
        free(fk);
    }*/
    for (i=0;i<fit_info.Npar;i++){
        Ci[i]=mean_and_error( argv[1],jack_tot,fit[i]);
    }
    fprintf(ftex,"\\begin{align}\n");
    fprintf(ftex,"& \\chi^2/d.o.f.= %+.5f \\pm \t%.2g \\\\ \n",chi2m[0],chi2m[1]);
    for (i=0;i<fit_info.Npar;i++){
         if (strcmp(namefile,"fit_Mpi_Fpi")==0 ||  strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M1")==0  ||  strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M2a")==0 ||  strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M2b")==0 ||strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M1a")==0 || strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M1b")==0  ||
           strcmp(namefile,"fit_FpiMpi4_GL_w0_M1a")==0        ){
             if(i==0)     fprintf(ftex,"& Bw_{0}= %+.5f \\pm \t%.2g   \\\\ \n",Ci[i][0],Ci[i][1]);
             else if(i==1)     fprintf(ftex,"& fw_{0}= %+.5f \\pm \t%.2g   \\\\ \n",Ci[i][0],Ci[i][1]);
             else if(i==2)     fprintf(ftex,"& \\bar{\\ell_3}= %+.5f \\pm \t%.2g   \\\\ \n",Ci[i][0],Ci[i][1]);
             else if(i==3)     fprintf(ftex,"& P_2= %+.5f \\pm \t%.2g   \\\\ \n",Ci[i][0],Ci[i][1]);
             else if(i==4)     fprintf(ftex,"& \\bar{\\ell_4}= %+.5f \\pm \t%.2g   \\\\ \n",Ci[i][0],Ci[i][1]);
             else if(i==5)     fprintf(ftex,"& P_4= %+.5f \\pm \t%.2g   \\\\ \n",Ci[i][0],Ci[i][1]);
             else
                 fprintf(ftex,"& P_{%d}= %+.5f \\pm \t%.2g   \\\\ \n",i,Ci[i][0],Ci[i][1]);

        }
        else
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
   /////////////////////////////////////////////////////////// compute m_ud
    double in;
    double *w0_estimate=(double*) malloc(sizeof(double)*jack_tot);
    double *xi=(double*) malloc(sizeof(double)*jack_tot);
    double *tmp3=(double*) malloc(sizeof(double)*jack_tot);
    double *w0_MILC_MeV=fake_sampling(argv[1],v_w0MeV,err_w0fm/197.326963 ,jack_tot,123);
    
    double **x=(double**) malloc(sizeof(double*)*jack_tot);
    for(j=0;j<jack_tot;j++){
        x[j]=(double*) malloc(sizeof(double)*fit_info.Nvar);
        //mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
        x[j][0]=0;//mlw
        x[j][1]=1e+6;//r0
        x[j][2]=phys_point[j][2];//Mpi2
        x[j][3]=phys_point[j][3];//fpi
        x[j][4]=1e+10;//L such that L/w=1e+4
        
    }
    
    for (j=0;j<jack_tot;j++){
        
       
       in=r1.MpiMeV[j]*r1.MpiMeV[j]/(r1.fpiMeV_exp[j]*r1.fpiMeV_exp[j]);
       xi[j]=rtbis(Mw2_over_fw2_chiral_FVE_a0_minus,in,fit_info.Npar,tif[j], 0.0001, 0.01, 1e-10);//gives mw
       phys_point[j][0]=xi[j];
       r1.fpiw[j]=fPSw_chiral_FVE(  xi[j],fit_info.Npar,tif[j]);
       w0_estimate[j]=r1.fpiw[j]/(v_fpiMeV_exp/197.326963);
       if(strcmp(namefile,"fit_FpiMpi4_GL_w0_M1a")==0 ){
         /*  xi[j]=rtbis_func_eq_input(fit_info.function ,
                                     0,//double n
                                     fit_info.Nvar, x,fit_info.Npar, tif[j], 
                                     0,//ivar
                                     in, 0.0001, 0.01, 1e-10);//gives mw  */
           
           x[j][0]=xi[j];
           r1.fpiw[j]=fit_info.function(1,fit_info.Nvar,x[j],fit_info.Npar,tif[j])  / pow(fit_info.function(0,fit_info.Nvar,x[j],fit_info.Npar,tif[j]),2);
           w0_estimate[j]=r1.fpiw[j]/(v_fpiMeV_exp/197.326963);
           
           w0_estimate[j]=0.171;//first guess
           double res=1;
           while (  res >1e-6){
                in=r1.MpiMeV[j]*r1.MpiMeV[j]*r1.MpiMeV[j]*r1.MpiMeV[j]*r1.fpiMeV_exp[j];
                in*=pow(w0_estimate[j]/197.326963,5);//input w0^5 fpi Mpi^4  // we need w0 in MeV-1
                
                xi[j]=rtbis_func_eq_input(fit_info.function ,
                                            1,//double n
                                            fit_info.Nvar, x[j],fit_info.Npar, tif[j], 
                                            0,//ivar
                                            in, 0.0001, 0.01, 1e-10); // find mw such imposing in= w0^5 fpi Mpi^4 
                x[j][0]=xi[j]; 
                double tmp_Mw2=fit_info.function(0,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);// compute Mpi^2w0 at mw
                double w0_tmp=sqrt(tmp_Mw2)/(r1.MpiMeV[j]/197.326963);
                res=w0_estimate[j]-w0_tmp;
                
                w0_estimate[j]=w0_estimate[j] +res/2.;
                r1.fpiw[j]=fit_info.function(1,fit_info.Nvar,x[j],fit_info.Npar,tif[j])  / pow(fit_info.function(0,fit_info.Nvar,x[j],fit_info.Npar,tif[j]),2);
                //printf("guess w0=%f   res=%f\n" ,w0_estimate[j],res);
           }
       }
       
       r1.w0fm[j]= w0_estimate[j] ;
       r1.w0MeV[j]= w0_estimate[j]/197.326963 ;
       r1.mlw[j]=xi[j];
       xi[j]=xi[j]/(r1.w0MeV[j]);
       r1.Bw[j]=fit_out.P[0][j];
       r1.fw[j]=fit_out.P[1][j];
       
       
       tmp3[j]=r1.fpiMeV_exp[j]/(r1.fw[j] /r1.w0MeV[j]);
       
       
    }
    double **C1=(double**) malloc(sizeof(double*)*2);
    
    fprintf(ftex,"Imposing $M_\\pi =%.2f \\pm %.2g$ and $f_\\pi=%.4f\\pm %.2g$ MeV\n",v_MpiMeV,err_MpiMeV,v_fpiMeV_exp,err_fpiMeV_exp);
    printf("Imposing $M_\\pi =%.2f \\pm %.2g$ and $f_\\pi=%.4f\\pm %.2g$ MeV\n",v_MpiMeV,err_MpiMeV,v_fpiMeV_exp,err_fpiMeV_exp);
    
    C1[0]=mean_and_error(argv[1],jack_tot,xi);
    fprintf(ftex,"\\begin{gather}\n   m_{ud}=(%g\\pm%.2g) MeV  \\\\ \n",C1[0][0],C1[0][1]);
    printf("\\begin{gather}\n   m_{ud}=(%g\\pm%.2g) MeV  \\\\ \n",C1[0][0],C1[0][1]);
    

    C1[1]=mean_and_error(argv[1],jack_tot,w0_estimate);
    fprintf(ftex,"w_0=(%g\\pm%.2g) fm (\\mbox{from }\\, f_\\pi)   \\\\ \n",C1[1][0],C1[1][1]);
    printf("w_0=(%g\\pm%.2g) fm (\\mbox{from }\\, f_\\pi)   \\\\ \n",C1[1][0],C1[1][1]);
    
    free(C1[1]);
    C1[1]=mean_and_error(argv[1],jack_tot,tmp3);
    fprintf(ftex,"f_\\pi/f=(%g\\pm%.2g)    \n\\end{gather}\n",C1[1][0],C1[1][1]);
    printf("f_\\pi/f=(%g\\pm%.2g)    \n\\end{gather}\n",C1[1][0],C1[1][1]);
    
    free(tmp3);

    free_2(2,C1);
    free(w0_estimate);
    free(xi);
    ////////////////////////print fit
    
   print_chiral_continuum_fit(argv, jack_tot,  fit_out,   fit_info,  phys_point,  AV, namefile,  head , grephJ);
 
  /*  int order=2;
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
    */
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
    free_2(jack_tot,x);
    free_2(jack_tot,tif);
    //free_tif(fit_info.Npar,fit);
    fclose(ftex);fclose(fgp);free(chi2m);free(tmp); free(Ci);
    
}


double **init_phys_point(int jack_tot){
    int j;
    double **phys_point=(double**) malloc(sizeof(double*)*jack_tot);

    for(j=0;j<jack_tot;j++){
        phys_point[j]=(double*) malloc(sizeof(double)*5);
        //double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
        phys_point[j][0]=0;
        phys_point[j][1]=0;
        phys_point[j][2]=result.MpiMeV[j];
        phys_point[j][3]=result.fpiMeV_exp[j];
       
    }
    /*printf("physical continuum point:\n");
    printf("Mpiw0= %g  \t Mpi=%g MeV\n",result.MpiMeV[jack_tot-1]*result.w0MeV[jack_tot-1],result.MpiMeV[jack_tot-1] );
    printf("MKw0= %g  \t MD=%g MeV\n",result.MKMeV[jack_tot-1]*result.w0MeV[jack_tot-1],result.MKMeV[jack_tot-1] );
    printf("MDw0= %g  \t MK=%g MeV\n",result.MDMeV[jack_tot-1]*result.w0MeV[jack_tot-1],result.MDMeV[jack_tot-1] );
    printf("MDsw0= %g  \t MDs=%g MeV\n",result.MDsMeV[jack_tot-1]*result.w0MeV[jack_tot-1],result.MDsMeV[jack_tot-1] );
    printf("r0=%g fm   \t   r0=%g MeV\n",result.w0fm[jack_tot-1],result.w0MeV[jack_tot-1]);
*/
    return phys_point;
}

int main(int argc, char **argv){
    
    int i,j;
    struct header  *head;
    struct database_file_jack  *jack_files;
    //double ***M_PS_GEVP_jack,***f_PS_jack;
    double ***M_PS_GEVP_jack_tot,***f_PS_jack_tot;
   
    double *tmp,**fit,*tmp1;
    
   error(argc!=3,1,"main ",
         "usage: ./fit_all_beta   jack/boot output_folder_analysis");

   error(strcmp(argv[1],"jack")!=0 && strcmp(argv[1],"boot")!=0 ,2,"main ",
         "choose jack or boot \n usage: ./fit_all_beta   jack/boot");
 
   // head=(struct header*) malloc(sizeof(struct header)*ensembles);
    
    

//    M_PS_GEVP_jack=(double***) malloc(sizeof(double**)*ensembles);
 //   f_PS_jack=(double***) malloc(sizeof(double**)*ensembles);
//    alloca_data_jack(dataJ);
    dataJ=(struct data_jack *) malloc(sizeof(struct data_jack)*ensembles);

    files_declarations(argv,&jack_files,&head);
    mass_index=init_mass_index_ave_r(head);
   // read_files_jack(jack_files,head,mass_index,M_PS_GEVP_jack,f_PS_jack);
    read_files_jack(jack_files,head,mass_index,dataJ);
    KM_FSE(jack_files, head, dataJ);
    Kf_FSE(jack_files, head, dataJ);
    
    int atimesObs,e;
    atimesObs=2*(2+1);
    //lattice_spacings*(N )
   //N=obs+ prior
   //obs=2, i.e. M and f
    index_a=(int**) malloc(sizeof(int*)*(atimesObs)); 
    for(i=0;i<atimesObs;i++)
       index_a[i]=(int*) malloc(sizeof(int)*ensembles); 
   
   for (e=0;e<ensembles;e++){
        index_a[0][e]=e;  
        index_a[1][e]=e+3;
        index_a[2][e]=e;
        index_a[3][e]=e+3;
        index_a[4][e]=0;
        index_a[5][e]=3;
   }

    
    
   printf("ensembles\n");
    for (i=0;i<ensembles;i++){
        
        printf("L%dT%d N=%d  musea=%.5f  KM=%f  Kf=%f ",head[i].l1,head[i].l0, jack_files[i].Njack-1 , head[i].musea, dataJ[i].KM[0], dataJ[i].Kf[0]); 
        tmp=mean_and_error(argv[1],jack_files[i].Njack,   dataJ[i].M_PS_jack[0] );
        printf("M_PS=%g  +-  %g  ",  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_files[i].Njack,   dataJ[i].f_PS_jack[0] );
        printf("  f_PS=%g  +-  %g  \t ",  tmp[0],tmp[1] );
        free(tmp);
        
        tmp1=(double*) malloc(sizeof(double)* jack_files[i].Njack);
        for(j=0;j<jack_files[i].Njack;j++)
            tmp1[j]= dataJ[i].M_PS_jack[0][j]*dataJ[i].M_PS_jack[0][j]/( dataJ[i].f_PS_jack[0][j]* dataJ[i].f_PS_jack[0][j]);
        tmp=mean_and_error(argv[1],jack_files[i].Njack,   tmp1 );
        printf("  M_PS^2/f_PS^2=%g  +-  %g  \n ",  tmp[0],tmp[1] );
        free(tmp);free(tmp1);
        
    }
 
    if( strcmp(argv[1],"jack")==0){
                gjack=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, dataJ);
                mysprintf(jack_files[0].sampling,NAMESIZE,"jack");
    }
    if( strcmp(argv[1],"boot")==0){
                gjack=create_generalised_boot( jack_files, head, &jack_tot ,mass_index, dataJ);
                mysprintf(jack_files[0].sampling,NAMESIZE,"boot");
    }
    init_Z( jack_files, head, jack_tot, &gjack, "w0","M1");

    
    int im,ik1,ik2;
   for(ik1=0;ik1<1;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
   for(ik2=ik1;ik2<4;ik2++){
    printf("ensambles after generalised jack ik2=%d  ik1=%d\n",ik2,ik1);

    for (i=0;i<ensembles;i++){
        im=mass_index[i][ik2][ik1];
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].M_PS_jack[im] );
        printf("L%dT%d  mu2=%.5f mu1=%.5f  M_PS=%g  +-  %g  ",head[i].l1,head[i].l0, head[i].k[head[i].nk+ik2], head[i].k[head[i].nk+ik1] ,  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].f_PS_jack[im] );
        printf("  f_PS=%g  +-  %g   \t",  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].Zp );
        printf("  Zp=%g  +-  %g    %s\n",  tmp[0],tmp[1]  ,jack_files[i].f_PS);
        free(tmp);
        
    }
    printf("\n");
    }}
 /*   for(ik1=0;ik1<1;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
   for(ik2=4;ik2<7;ik2++){
    printf("ensambles after generalised jack ik2=%d  ik1=%d\n",ik2,ik1);
    for (i=0;i<ensembles;i++){
        im=mass_index[i][ik2][ik1];
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].M_PS_GEVP_jack[im] );
        printf("L%dT%d  mu2=%.5f mu1=%.5f  M_PS=%g  +-  %g  ",head[i].l1,head[i].l0, head[i].k[head[i].nk+ik2], head[i].k[head[i].nk+ik1] ,  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].f_PS_ls_ss_jack[im] );
        printf("  f_PS=%g  +-  %g   \t",  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].Zp );
        printf("  Zp=%g  +-  %g    %s\n",  tmp[0],tmp[1]  ,jack_files[i].f_PS);
        free(tmp);
        
    }
     printf("\n");
    }}
    for(ik1=1;ik1<2;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
   for(ik2=4;ik2<7;ik2++){
    printf("ensambles after generalised jack ik2=%d  ik1=%d\n",ik2,ik1);
    for (i=0;i<ensembles;i++){
        im=mass_index[i][ik2][ik1];
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].M_PS_GEVP_jack[im] );
        printf("L%dT%d  mu2=%.5f mu1=%.5f  M_PS=%g  +-  %g  ",head[i].l1,head[i].l0, head[i].k[head[i].nk+ik2], head[i].k[head[i].nk+ik1] ,  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].f_PS_ls_ss_jack[im] );
        printf("  f_PS=%g  +-  %g   \t",  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].Zp );
        printf("  Zp=%g  +-  %g    %s\n",  tmp[0],tmp[1]  ,jack_files[i].f_PS);
        free(tmp);
        
    }
    printf("\n");
    }}
   */ 
  
    double *xi=(double*) malloc(sizeof(double)*jack_tot);
    double *fw_phys=(double*) malloc(sizeof(double)*jack_tot);
    
    double *Mw2=(double*) malloc(sizeof(double)*jack_tot);
    //double *B_point=(double*) malloc(sizeof(double)*2);
        double **Ci;
    double *w0_estimate=(double*) malloc(sizeof(double)*jack_tot);
    result.fpiw=(double*) malloc(sizeof(double)*jack_tot);

    struct fit_type fit_info;
    struct fit_result  fit_out;
    fit_info.Nvar=13;
    
    double    **phys_point=init_phys_point(jack_tot);  
    
    
    printf("\n\nOMEGA MK\n");

    Ci=(double**) malloc(sizeof(double*)*2);
    
    fit=fit_Omegaw0_from_M(jack_files, head , jack_tot, mass_index, gjack,  &result );
    
    tmp1=(double*) malloc(sizeof(double)*jack_tot);
    for (j=0;j<jack_tot;j++){
        tmp1[j]=fit[0][j]*197.326963;
    }

    Ci[0]=mean_and_error(argv[1],jack_tot,tmp1);
    printf("Imposing $M_\\K =%.2f \\pm %.2g$ , $M_\\pi=%.4f\\pm %2g$ fm  and M_\\Omega=%g \\pm %g\n",v_MKMeV,err_MKMeV,v_MpiMeV,err_MpiMeV,v_MOmegaMeV,err_MOmegaMeV);
    printf("   w_0=(%g\\pm%.2g) fm \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1]);
     
    
    for (i=0;i<1;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);    
    free(tmp1);
/*
    printf("\n\n///////////////////////////////////////Pion of m_l ///////////////////////\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=fit_Fpi_and_Mpi;
        
        
        
    fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info);
    //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, gjack, head, "pion","fit_Mpi_Fpi");
  */  
    double in;
    double **tif;
    double **C1;
   double *tmp3=(double*) malloc(sizeof(double)*jack_tot);

    printf("\n\n///////////////////////////////////////Pion Mpi^2 and fpiMpi^4  GL  w0 M1a ///////////////////////\n");
    fit_info.Npar=8;
    fit_info.N=2;
    fit_info.function=fit_FpiMpi4_and_Mpi2_GL;
        
    init_Z( jack_files, head, jack_tot, &gjack, "w0","M1a");

        
    fit_out=fit_Mpi_fwMpi4_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info);
 
    //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion","fit_FpiMpi4_GL_w0_M1a");
 
    
    printf("\n\n///////////////////////////////////////Pion Mpi^2 and fpiMpi^4  linear  w0 M1a ///////////////////////\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=fit_FpiMpi4_and_Mpi2_linear;
        
    init_Z( jack_files, head, jack_tot, &gjack, "w0","M1a");

    double threshold_Mpiw=0.20 ;//0.164 ~19MeV 
    fit_out=fit_Mpi_fwMpi4_chiral_FVE_clover_threshold(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info, threshold_Mpiw);
 
    //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion","fit_FpiMpi4_GL_w0_M1a_threshold");
 
   
   
   
    printf("\n\n///////////////////////////////////////Pion of m_l GL   w0 M1a ///////////////////////\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=fit_Fpi_and_Mpi_GL;
        
    init_Z( jack_files, head, jack_tot, &gjack, "w0","M1a");

        
    fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info);
 
    //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion","fit_Mpi_Fpi_GL_w0_M1a");
    
    threshold_Mpiw=0.20;
    fit_out=fit_Mpi_fw_chiral_FVE_clover_treshold(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info, threshold_Mpiw);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion","fit_Mpi_Fpi_GL_w0_M1a");
    
    
    printf("\n\n///////////////////////////////////////Pion of m_l GL   w0 M1b ///////////////////////\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=fit_Fpi_and_Mpi_GL;
        
    init_Z( jack_files, head, jack_tot, &gjack, "w0","M1b");

        
    fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info);
 
    //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion","fit_Mpi_Fpi_GL_w0_M1b");
    
    
    printf("\n\n///////////////////////////////////////Pion of m_l GL   w0 M2a ///////////////////////\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=fit_Fpi_and_Mpi_GL;
        
         init_Z( jack_files, head, jack_tot, &gjack, "w0","M2a");
   
        
     fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info);
    
    //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result,gjack, head, "pion","fit_Mpi_Fpi_GL_w0_M2a");
      printf("\n\n///////////////////////////////////////Pion of m_l GL   w0 M2b ///////////////////////\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=fit_Fpi_and_Mpi_GL;
        
         init_Z( jack_files, head, jack_tot, &gjack, "w0","M2b");

    tmp3=(double*) malloc(sizeof(double)*jack_tot);   
     fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info);

     //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion","fit_Mpi_Fpi_GL_w0_M2b");
    
    
    printf("\n\n///////////////////////////////////////K of m_s GL   w0 M2b ///////////////////////\n");

    fit_info.Npar=7;
    fit_info.N=2;
    fit_info.Nvar=10;
    fit_info.function=fit_FK_and_MK_GL;
    
     Ci=(double**) malloc(sizeof(double*)*2);
    result.fkw=(double*) malloc(sizeof(double)*jack_tot);
    result.msw=(double*) malloc(sizeof(double)*jack_tot);
    
    result.fk_fpi=(double*) malloc(sizeof(double)*jack_tot);
    result.ms_mud=(double*) malloc(sizeof(double)*jack_tot);
    
    fit=fit_MK_fK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, "FK_and_MK_GL_M2b" ,argv);
    //fit_MK_double_chiral_FVE_P40(jack_files, head , jack_tot, mass_index, gjack,  &result );
    for (j=0;j<jack_tot;j++){
        result.msw[j]=fit[0][j];
        result.fkw[j]=fit[1][j];
        
        result.ms_mud[j]=result.msw[j]/result.mlw[j];
       // printf("%f\n",result.ms_mud[j]);
        result.fk_fpi[j]=result.fkw[j]/result.fpiw[j];
        fit[0][j]=fit[0][j]/result.w0MeV[j];
        fit[1][j]=fit[1][j]/result.w0MeV[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);

    
    
    printf("Imposing $M_K =%.2f \\pm %.2g$ and $w_0=%.4f$ fm\n",v_MKMeV,err_MKMeV,result.w0fm[jack_tot-1]);
    printf("\\begin{gather}\n   m_{s}=(%g\\pm%.2g) MeV   \\\\ \n",Ci[0][0],Ci[0][1]);
    printf("   f_{K}=(%g\\pm%.2g) MeV  \n\\end{gather} \n",Ci[1][0],Ci[1][1]);
    
     for (i=0;i<2;i++)
        free(Ci[i]);
     
    Ci[0]=mean_and_error(argv[1],jack_tot,result.ms_mud);
    Ci[1]=mean_and_error(argv[1],jack_tot,result.fk_fpi);
    
    printf("\\begin{gather}\n m_{s}/m_{ub}=(%g\\pm%.2g),\\quad \t  f_{K}/f_{\\pi}=(%g\\pm%.2g)   \n\\end{gather} \n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1]);
    
    
    
    
    for (i=0;i<2;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    
 /*   printf("\n\n///////////////////////////////////////fKoverfpi   w0 M2b ///////////////////////\n");

    fit_info.Npar=4;
    fit_info.N=1;
    fit_info.function=fit_FKoverFpi_GL;
    
     Ci=(double**) malloc(sizeof(double*)*2);
       
    fit=fit_fKoverfpi_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, "fKoverfpi_GL_M2b" ,argv);
    //fit_MK_double_chiral_FVE_P40(jack_files, head , jack_tot, mass_index, gjack,  &result );
   
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);

    
    
    printf("\\begin{gather}\n   f_{k}/f_\\pi=(%g\\pm%.2g) MeV    \n",Ci[0][0],Ci[0][1]);
    printf("     \n\\end{gather} \n");

    
    
    for (i=0;i<2;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);*/
 
 
     printf("\n\n///////////////////////////////////////MKoverMpi fKoverfpi   w0 M2b ///////////////////////\n");

    fit_info.Npar=8;
    fit_info.N=2;
    fit_info.function=fit_MK_Mpi_FK_Fpi_GL;
    
     Ci=(double**) malloc(sizeof(double*)*2);
       
    fit=fit_MK_Mpi_fK_fpi_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, "MK_Mpi_fK_fpi_GL_M2b" ,argv);
    //fit_MK_double_chiral_FVE_P40(jack_files, head , jack_tot, mass_index, gjack,  &result );
    for (j=0;j<jack_tot;j++){
        result.msw[j]=fit[0][j];        
        result.ms_mud[j]=result.msw[j]/result.mlw[j];
        result.fk_fpi[j]=fit[1][j];
        
        fit[0][j]=fit[0][j]/result.w0MeV[j];
        fit[1][j]=fit[1][j]/result.w0MeV[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);

    
    
    printf("Imposing $M_K =%.2f \\pm %.2g$ and $w_0=%.4f$ fm\n",v_MKMeV,err_MKMeV,result.w0fm[jack_tot-1]);
    printf("\\begin{gather}\n   m_{s}=(%g\\pm%.2g) MeV   \\\\ \n",Ci[0][0],Ci[0][1]);
    //printf("\\end{gather} \n");
    
     for (i=0;i<2;i++)
        free(Ci[i]);
     
    Ci[0]=mean_and_error(argv[1],jack_tot,result.ms_mud);
    Ci[1]=mean_and_error(argv[1],jack_tot,result.fk_fpi);
    
    printf(" m_{s}/m_{ub}=(%g\\pm%.2g),\\\\ \n  f_{K}/f_{\\pi}=(%g\\pm%.2g)   \n\\end{gather} \n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1]);
    
    
    
    for (i=0;i<2;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    
   /* 
         printf("\n\n///////////////////////////////////////MKoverMpi fKoverfpi   w0 M2b  fix f from pion///////////////////////\n");

    fit_info.Npar=7;
    fit_info.N=2;
    fit_info.function=fit_MK_Mpi_FK_Fpi_GL_fix_f;
    
     Ci=(double**) malloc(sizeof(double*)*2);
       
    fit=fit_MK_Mpi_fK_fpi_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, "MK_Mpi_fK_fpi_GL_fix_f_M2b" ,argv);
    //fit_MK_double_chiral_FVE_P40(jack_files, head , jack_tot, mass_index, gjack,  &result );
    for (j=0;j<jack_tot;j++){
        result.msw[j]=fit[0][j];        
        result.ms_mud[j]=result.msw[j]/result.mlw[j];
        result.fk_fpi[j]=fit[1][j];
        
        fit[0][j]=fit[0][j]/result.w0MeV[j];
        fit[1][j]=fit[1][j]/result.w0MeV[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);

    
    
    printf("Imposing $M_K =%.2f \\pm %.2g$ and $w_0=%.4f$ fm\n",v_MKMeV,err_MKMeV,result.w0fm[jack_tot-1]);
    printf("\\begin{gather}\n   m_{s}=(%g\\pm%.2g) MeV   \\\\ \n",Ci[0][0],Ci[0][1]);
    //printf("\\end{gather} \n");
    
     for (i=0;i<2;i++)
        free(Ci[i]);
     
    Ci[0]=mean_and_error(argv[1],jack_tot,result.ms_mud);
    Ci[1]=mean_and_error(argv[1],jack_tot,result.fk_fpi);
    
    printf(" m_{s}/m_{ub}=(%g\\pm%.2g), \\\\ \n f_{K}/f_{\\pi}=(%g\\pm%.2g)   \n\\end{gather} \n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1]);
    
    
    
    free(tif);
    for (i=0;i<2;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    */
  /*  
    printf("\n\n///////////////////////////////////////K of m_s GL   w0 M2b spline ///////////////////////\n");

    fit_info.Npar=7;
    fit_info.N=2;
    fit_info.function=fit_FK_and_MK_GL;
    
     Ci=(double**) malloc(sizeof(double*)*2);
    result.fkw=(double*) malloc(sizeof(double)*jack_tot);
    result.msw=(double*) malloc(sizeof(double)*jack_tot);
    
    result.fk_fpi=(double*) malloc(sizeof(double)*jack_tot);
    result.ms_mud=(double*) malloc(sizeof(double)*jack_tot);
    
    fit=fit_MK_fK_chiral_spline_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, "FK_and_MK_spline_GL_M2b" ,argv);
    //fit_MK_double_chiral_FVE_P40(jack_files, head , jack_tot, mass_index, gjack,  &result );
    for (j=0;j<jack_tot;j++){
        result.msw[j]=fit[0][j];
        result.fkw[j]=fit[1][j];
        
        result.ms_mud[j]=result.msw[j]/result.mlw[j];
       // printf("%f\n",result.ms_mud[j]);
        result.fk_fpi[j]=result.fkw[j]/result.fpiw[j];
        fit[0][j]=fit[0][j]/result.w0MeV[j];
        fit[1][j]=fit[1][j]/result.w0MeV[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);

    
    
    printf("Imposing $M_K =%.2f \\pm %.2g$ and $w_0=%.4f$ fm\n",v_MKMeV,err_MKMeV,result.w0fm[jack_tot-1]);
    printf("\\begin{gather}\n   m_{s}=(%g\\pm%.2g) MeV   \\\\ \n",Ci[0][0],Ci[0][1]);
    printf("   f_{K}=(%g\\pm%.2g) MeV \\\\ \n",Ci[1][0],Ci[1][1]);
    
     for (i=0;i<2;i++)
        free(Ci[i]);
     
    Ci[0]=mean_and_error(argv[1],jack_tot,result.ms_mud);
    Ci[1]=mean_and_error(argv[1],jack_tot,result.fk_fpi);
    
    printf(" m_{s}/m_{ub}=(%g\\pm%.2g),\\\\ \n  f_{K}/f_{\\pi}=(%g\\pm%.2g)   \n\\end{gather} \n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1]);
    
    
    
    for (i=0;i<2;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    
    
    printf("\n\n///////////////////////////////////////MKoverMpi fKoverfpi   w0 M2b spline ///////////////////////\n");

    fit_info.Npar=8;
    fit_info.N=2;
    fit_info.function=fit_MK_Mpi_FK_Fpi_GL;
    
     Ci=(double**) malloc(sizeof(double*)*2);
       
    fit=fit_MK_Mpi_fK_fpi_chiral_spline_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, "MK_Mpi_fK_fpi_GL_M2b" ,argv);
    //fit_MK_double_chiral_FVE_P40(jack_files, head , jack_tot, mass_index, gjack,  &result );
    for (j=0;j<jack_tot;j++){
        result.msw[j]=fit[0][j];        
        result.ms_mud[j]=result.msw[j]/result.mlw[j];
        result.fk_fpi[j]=fit[1][j];
        
        fit[0][j]=fit[0][j]/result.w0MeV[j];
        fit[1][j]=fit[1][j]/result.w0MeV[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);

    
    
    printf("Imposing $M_K =%.2f \\pm %.2g$ and $w_0=%.4f$ fm\n",v_MKMeV,err_MKMeV,result.w0fm[jack_tot-1]);
    printf("\\begin{gather}\n   m_{s}=(%g\\pm%.2g) MeV   \\\\ \n",Ci[0][0],Ci[0][1]);
    //printf("\\end{gather} \n");
    
     for (i=0;i<2;i++)
        free(Ci[i]);
     
    Ci[0]=mean_and_error(argv[1],jack_tot,result.ms_mud);
    Ci[1]=mean_and_error(argv[1],jack_tot,result.fk_fpi);
    
    printf(" m_{s}/m_{ub}=(%g\\pm%.2g),\\\\ \n  f_{K}/f_{\\pi}=(%g\\pm%.2g)   \n\\end{gather} \n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1]);
    
    
    
    for (i=0;i<2;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    */
    
    
    
    return 0;
}
