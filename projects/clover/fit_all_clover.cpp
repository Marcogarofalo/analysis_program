#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <memory>

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
#include <vector>



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
        r[8][0][0]=0.6699; r[8][0][1]=0.0027;  
        r[8][1][0]=0.6990; r[8][1][1]=0.0022;    
        r[8][2][0]=0.7274; r[8][2][1]=0.0018;  

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
        r[4][0][0]=0.6699; r[4][0][1]=0.0027;  
        r[4][1][0]=0.6990; r[4][1][1]=0.0022;  
        r[4][2][0]=0.7274; r[4][2][1]=0.0018;  

    }if(ensembles>3){
        //A12
        /*r[3][0][0]=0.776; r[3][0][1]=0.005;  
        r[3][1][0]=0.811; r[3][1][1]=0.004;  
        r[3][2][0]=0.844; r[3][2][1]=0.004;  */
        r[3][0][0]=0.7828; r[3][0][1]=0.0039;  
        r[3][1][0]=0.8194; r[3][1][1]=0.0032;  
        r[3][2][0]=0.8553; r[3][2][1]=0.0026; 


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
        r[1][0][0]=0.8077; r[1][0][1]=0.0043;  
        r[1][1][0]=0.8408; r[1][1][1]=0.0032;  
        r[1][2][0]=0.8741; r[1][2][1]=0.0026; 
        
    }if(ensembles>0){
        //A53
        /*r[0][0][0]=0.811; r[0][0][1]=0.003;  
        r[0][1][0]=0.843; r[0][1][1]=0.002;  
        r[0][2][0]=0.878; r[0][2][1]=0.002;  */
        r[0][0][0]=0.8093; r[0][0][1]=0.0030;  
        r[0][1][0]=0.8436; r[0][1][1]=0.0023;  
        r[0][2][0]=0.8777; r[0][2][1]=0.0018;

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
        Mw2*=2*Bw*mw;
        Mw2*=(1-0.25 *Delta)*(1-0.25 *Delta);
        
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




double fit_Fpi_and_Mpi_GL_NL0_am_m2(int n, int Nvar, double *x,int Npar,double  *P){
    
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
        
        Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2 +(1./(w0*w0))*mw*P[6] + mw*mw *P[8];
        Mw2*=2*Bw*mw;
        Mw2*=(1-0.25 *Delta)*(1-0.25 *Delta);
    }
    if (n==1){
        
        Mw2=fw*(1-2*xi*log(xi)+P3*xi)+(1./(w0*w0))*P4+(1./(w0*w0))*mw*P[7] + mw*mw *P[9];
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



double fit_Fpi_and_Mpi_GL_NL0_am(int n, int Nvar, double *x,int Npar,double  *P){
    
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
        
        Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2 +(1./(w0*w0))*mw*P[6] ;
        Mw2*=2*Bw*mw;
        Mw2*=(1-0.25 *Delta)*(1-0.25 *Delta);
    }
    if (n==1){
        
        Mw2=fw*(1-2*xi*log(xi)+P3*xi)+(1./(w0*w0))*P4+(1./(w0*w0))*mw*P[7] ;
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



double fit_Fpi_and_Mpi_GL_NL0_am_fonly(int n, int Nvar, double *x,int Npar,double  *P){
    
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
        
        Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2  ;
        Mw2*=2*Bw*mw;
        Mw2*=(1-0.25 *Delta)*(1-0.25 *Delta);
    }
    if (n==1){
        
        Mw2=fw*(1-2*xi*log(xi)+P3*xi)+(1./(w0*w0))*P4+(1./(w0*w0))*mw*P[6] ;
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



double fit_FpiMpi4_and_Mpi2_noGL(int n, int Nvar, double *x,int Npar,double  *P){
    
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
        Mw2*=2*Bw*mw;
        Mw2*=(1-0.25 *Delta)*(1-0.25 *Delta);
        
    }
    if (n==1){// fit fot w0^5fpi Mpi^4
        double ML=sqrt(2*Bw*mw)*L_w;
        Mw2= (4.*pi)* (4.*pi)* (4.*pi)* (4.*pi);
        Mw2*= fw * fw *fw * fw * fw *xi*xi;
        Mw2*= (1.+ 2.*(l4b-l3b) *xi + xi*xi *P[6] +(1./(w0*w0))*P4 );
        //Mw2*=(1.+ P[7]*xi*xi *exp(-ML)/ pow(ML,3./2.));
        
    }
    if (n==2){
        Mw2=(1-0.25 *Delta);  //KM   M(inf)=M(L)/KM
        
    }
    if (n==3){
        Mw2=(1+Delta);  //KfMpi4   M(inf)=M(L)/Kf
        
    }
    
    return Mw2;
    
}


double fit_FpiMpi4_and_Mpi2_noGL_noA2(int n, int Nvar, double *x,int Npar,double  *P){
    
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
        Mw2*=2*Bw*mw;
        Mw2*=(1-0.25 *Delta)*(1-0.25 *Delta);
        
    }
    if (n==1){// fit fot w0^5fpi Mpi^4
        double ML=sqrt(2*Bw*mw)*L_w;
        Mw2= (4.*pi)* (4.*pi)* (4.*pi)* (4.*pi);
        Mw2*= fw * fw *fw * fw * fw *xi*xi;
        Mw2*= (1.+ 2.*(l4b-l3b) *xi  +(1./(w0*w0))*P4 );
        //Mw2*=(1.+ P[7]*xi*xi *exp(-ML)/ pow(ML,3./2.));
        
    }
    if (n==2){
        Mw2=(1-0.25 *Delta);  //KM   M(inf)=M(L)/KM
        
    }
    if (n==3){
        Mw2=(1+Delta);  //KfMpi4   M(inf)=M(L)/Kf
        
    }
    
    return Mw2;
    
}

double FpiMpi4_fromM(int n, int Nvar, double *x,int Npar,double  *P){
    
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
    double A1=l4b+2*log(   v_Mpiw0/(4*pi*fw) );
    double Delta=FVE_GL_fast( L_w, mw, fw, Bw);
    
    xi=dmpi2*w0*w0/(16.*pi*pi*fw*fw);
    
    Mw2= (4.*pi)* (4.*pi)* (4.*pi)* (4.*pi);
    Mw2*= fw * fw *fw * fw * fw *xi*xi;
    Mw2*= (1.-2.*xi*log(xi)+ 2.*A1 *xi   );
    Mw2/=(w0*w0*w0*w0*w0);
    
    return Mw2;
}


double fit_FpiMpi4_fromM_and_Mpi2_noGL_noA2(int n, int Nvar, double *x,int Npar,double  *P){
    
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
    double A1=l4b+2*log(   v_Mpiw0/(4*pi*fw) );
    
    double Delta=FVE_GL_fast( L_w, mw, fw, Bw);
    
    
    if (n==0){
        xi=2*Bw*mw/(16.*pi*pi*fw*fw);
        
        Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2;
        Mw2*=2*Bw*mw;
        Mw2*=(1-0.25 *Delta)*(1-0.25 *Delta);
        
        Mw2=2*Bw*mw*(1+ P1*mw+(1./(w0*w0))*P2);
    }
    if (n==1){// fit fot w0^5fpi Mpi^4
        double ML=sqrt(2*Bw*mw)*L_w;
        xi=dmpi2*w0*w0/(16.*pi*pi*fw*fw);
        //xi/=(1-0.25 *Delta)*(1-0.25 *Delta);
        Mw2= (4.*pi)* (4.*pi)* (4.*pi)* (4.*pi);
        Mw2*= fw * fw *fw * fw * fw *xi*xi;
        Mw2*= (1.-2.*xi*log(xi)+ 2.*A1 *xi  +(1./(w0*w0))*P4 );
        //Mw2*=(1.+ P[7]*xi*xi *exp(-ML)/ pow(ML,3./2.));
        
    }
    if (n==2){
        Mw2=(1-0.25 *Delta);  //KM   M(inf)=M(L)/KM
        
    }
    if (n==3){
        Mw2=(1+Delta);  //KfMpi4   M(inf)=M(L)/Kf
        
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
        Mw2*=2*Bw*mw;
        Mw2*=(1-0.25 *Delta)*(1-0.25 *Delta);
        
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
        Mw2=(1+Delta);  //KfMpi4   M(inf)=M(L)/Kf
        
    }
    
     return Mw2;
    
}


double fit_FpiMpi4_and_Mpi2_GL_am(int n, int Nvar, double *x,int Npar,double  *P){
    
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
        
        Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2 +(1./(w0*w0))*mw*P[8];
        Mw2*=2*Bw*mw;
        Mw2*=(1-0.25 *Delta)*(1-0.25 *Delta);
        
    }
    if (n==1){// fit fot w0^5fpi Mpi^4
        double ML=sqrt(2*Bw*mw)*L_w;
        Mw2= (4.*pi)* (4.*pi)* (4.*pi)* (4.*pi);
        Mw2*= fw * fw *fw * fw * fw *xi*xi;
        Mw2*= (1.+ 2.*(l4b-l3b) *xi + xi*xi *P[6] +(1./(w0*w0))*P4+(1./(w0*w0))*mw*P[9] );
        Mw2*=(1.+ P[7]*xi*xi *exp(-ML)/ pow(ML,3./2.));
        
    }
    if (n==2){
        Mw2=(1-0.25 *Delta);  //KM   M(inf)=M(L)/KM
        
    }
    if (n==3){
        double ML=sqrt(2*Bw*mw)*L_w;
        Mw2=(1.+ P[7]*xi*xi *exp(-ML)/ pow(ML,3./2.));  //KfMpi4   M(inf)=M(L)/Kf
        
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
        Mw2*=2*Bw*mw;
        Mw2/=(1-0.25 *Delta)*(1-0.25 *Delta);
        
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



double fit_FK_and_MK_GL_noP2(int n, int Nvar, double *x,int Npar,double  *P){
    
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
        
        Mw2=1+P[1]*mw + P[2]/(w0*w0);
        Mw2*=P[0]*(mw+msw);//*(1-0.25 *Delta)*(1-0.25 *Delta);
        
    }
    if (n==1){
        
        Mw2=P[3]*(1-(3./4.)*xi*log(xi)+P[4]*xi+(1./(w0*w0))*P[5]);
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


double fit_MD_fD(int n, int Nvar, double *x,int Npar,double  *P){
    
    double MKw=0,xi;
    double pi=3.141592653589793;
    
    double mlw=x[0], w0=x[1];
    
    double     P0=P[0], P1w=P[1], P2ww=P[2],   P3ww=P[3];
    double     P0f, P1fw, P2fww,   P3fww;
    if (n==1){
        P0f=P[4]; P1fw=P[5]; P2fww=P[6];   P3fww=P[7];
    }
    
    if (n==0)
        MKw=P0+P1w*mlw+P2ww*mlw*mlw+P3ww*(1./(w0*w0));
    else if (n==1){
        MKw=P0f+P1fw*mlw+P2fww*mlw*mlw+P3fww*(1./(w0*w0));
    }
    else if (n==2){
        MKw=1;
    }
    else if (n==3){
        MKw=1;
    }
    
    
    
    return MKw;
    
}

double fit_MD_fD_noP2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double MKw=0,xi;
    double pi=3.141592653589793;
    
    double mlw=x[0], w0=x[1];
    
    double     P0=P[0], P1w=P[1],    P3ww=P[2];
    double P0f, P1fw, P2fww,   P3fww;
    if (n==1){
             P0f=P[3]; P1fw=P[4]; P2fww=P[5];   P3fww=P[6];
    }
    
    
    if (n==0)
        MKw=P0+P1w*mlw  +P3ww*(1./(w0*w0));
    else if (n==1){
        MKw=P0f+P1fw*mlw+P2fww*mlw*mlw+P3fww*(1./(w0*w0));
    }
    else if (n==2){
        MKw=1;
    }
    else if (n==3){
        MKw=1;
    }
    
    
    
    return MKw;
    
}

double fit_Mpi_MK(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    
    double     P0=P[0], P1=P[1], P2=P[2];// fw=P[7];
    
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3],  msw=x[5],   MK2w2=x[6] , Bw=x[8], fw=x[9], Mpiw=x[10], fpiw=x[11] ;
    
    int Lsize=(int(x[4]));
    double L_w=(x[4]) /w0;
    double KM,Kf;
    
    //double Delta=FVE_GL_fast( L_w, mw, fw, Bw);
    xi=Mpiw*Mpiw/(16.*pi*pi*fpiw*fpiw);
    double Delta=FVE_GL_Mpi( L_w /* L/w0 */,  xi,fpiw);
    //Delta=0;
    //xi*= ((1+ Delta)*(1+ Delta)) /((1-0.25 *Delta)*(1-0.25 *Delta));
    
    //printf("delta = %f\n",Delta);
    double DeltaK=Delta*(3./8.);
    
    
    if (n==0){
        //Mw2 = 0.5*(1. + msw/mw) *P0* (1. - xi* log(xi) + P1* xi + P2 *xi*xi + P3/(w0*w0) );
        //Mw2*=1./((1-0.25 *Delta)*(1-0.25 *Delta));
        //Mw2=P0*xi*(1+P1*xi-4*xi* log(xi)+ P2/(w0*w0));
        Mw2=P0*xi*(1+P1*xi+ P2/(w0*w0)   +P[3]*xi*xi) ;
        Mw2*=((1-0.25 *Delta)*(1-0.25 *Delta))/1.;
    }
    if (n==1){
        Mw2 = 1 ;
       // Mw2*=(1+DeltaK)/(1+Delta);
    }
    if (n==2){
        Mw2=((1-0.25 *Delta)*(1-0.25 *Delta))/1.;
        //Mw2=1;
    }
    if (n==3){
        //Mw2=(1+DeltaK)/(1+Delta);
        Mw2=1;
    }
    return Mw2;
    
}
double fit_MK_Mpi_FK_Fpi_GL_noP2(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double     P0, P1,  P3 ,P4, P5 ,P6;// fw=P[7];
    
    P0=P[0]; P1=P[1]; P3=P[2];
    
    if (n>=1){
        P0=P[0]; P1=P[1];   P3=P[2]; P4=P[3]; P5=P[4];P6=P[5];// fw=P[7];
    }
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3],  msw=x[5],   MK2w2=x[6] , Bw=x[8], fw=x[9] ;
    
    int Lsize=(int(x[4]));
    double L_w=(x[4]) /w0;
    double KM,Kf;
    
    double Delta=FVE_GL_fast( L_w, mw, fw, Bw);
    double DeltaK=Delta*(3./8.);
    
    xi=2*Bw*mw/(16.*pi*pi*fw*fw);
    if (n==0){
        Mw2 = 0.5*(1. + msw/mw) *P0* (1. - xi* log(xi) + P1* xi +  P3/(w0*w0) );
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

double fit_MK_Mpi_FK_Fpi_GL(int n, int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    
    double     P0=P[0], P1=P[1], P2=P[2],  P3=P[3] ,P4=P[4], P5=P[5] ,P6=P[6];// fw=P[7];
   
    double mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3],  msw=x[5],   MK2w2=x[6] , Bw=x[8], fw=x[9] ;

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
    
    setup_reading_jack( argv,&((*jack_files)[0]),&((*head)[0]),"/media/marco/data/PRACE/beta1.726/cA211ab.53.24/analysis/main/jackknife");  
    (*jack_files)[0].a=0.096;
    if (ensembles>1){
        setup_reading_jack(argv, &((*jack_files)[1]),&((*head)[1]),"/media/marco/data/PRACE/beta1.726/cA211ab.40.24/analysis/main/jackknife");  
        (*jack_files)[1].a=0.096;
    }
    if (ensembles>2){
        setup_reading_jack(argv, &((*jack_files)[2]),&((*head)[2]),"/media/marco/data/PRACE/beta1.726/cA211ab.30.32/analysis/main/jackknife");  
        (*jack_files)[2].a=0.096;
    }
    if (ensembles>3){
     //   setup_reading_jack(argv, &jack_files[3],&head[3],"../../beta1.726/cA211ab.12.48/analysis/main/jackknife");  
        setup_reading_jack(argv, &((*jack_files)[3]),&((*head)[3]),"/media/marco/data/PRACE/beta1.726/cA211ab.12.48_no_rew/analysis/main/jackknife");  
        (*jack_files)[3].a=0.096;
    }
    if (ensembles>4){
        setup_reading_jack(argv, &((*jack_files)[4]),&((*head)[4]),"/media/marco/data/PRACE/beta1.778/cB211ab.25.48/analysis/main/jackknife");  
        (*jack_files)[4].a=0.081;
    }
    if (ensembles>5){
        setup_reading_jack(argv,&((*jack_files)[5]),&((*head)[5]),"/media/marco/data/PRACE/beta1.778/cB211ab.14.64/analysis/main/jackknife");  
        (*jack_files)[5].a=0.070;
    }
    if (ensembles>6){
        setup_reading_jack(argv,&((*jack_files)[6]),&((*head)[6]),"/media/marco/data/PRACE/beta1.778/cB211ab.072.64/analysis/main/jackknife");  
        (*jack_files)[6].a=0.081;
    }
    
    if (ensembles>7){
        setup_reading_jack(argv, &((*jack_files)[7]),&((*head)[7]),"/media/marco/data/PRACE/beta1.836/cC211ab.06.80/analysis/main/jackknife");  
        (*jack_files)[7].a=0.070;
    }
    if (ensembles>8){
        setup_reading_jack(argv, &((*jack_files)[8]),&((*head)[8]),"/media/marco/data/PRACE/beta1.778/cB211ab.25.32/analysis/main/jackknife");  
        (*jack_files)[8].a=0.081;
    }
    if (ensembles>9){
        setup_reading_jack(argv, &((*jack_files)[9]),&((*head)[9]),"/media/marco/data/PRACE/beta1.836/cC211ab.20.48/analysis/main/jackknife");  
        (*jack_files)[9].a=0.081;
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



void  create_fake_distribution(const char *jackboot,double **w0A,double **w0B,double **w0C,double **ZpA,double **ZpB,double **ZpC,int jack_tot, const char *scaletype, const char *Mtype, double ***w0, double **scalefm){
    
      if (strcmp(scaletype,"w0")==0 || strcmp(scaletype,"w0xensemble")==0 ){
          //(*w0A)=fake_sampling(jackboot,1.8381,  0.0037,jack_tot,123);  // petros-from jacob 31/07/20 fit   in M_PS^2/f_PS^2
          //(*w0B)=fake_sampling(jackboot,2.1316 ,0.0024 ,jack_tot,1234);// petros-from jacob 31/07/20 fit  M_PS^2/f_PS^2
          //(*w0C)=fake_sampling(jackboot,2.5039, 0.0017,jack_tot,12345);//  petros-from jacob 31/07/20 fit  M_PS^2/f_PS^2
          (*w0A)=fake_sampling(jackboot,1.83548,  0.00353,jack_tot,123);  // petros-from jacob 31/07/20 fit   in M_PS^2/f_PS^2
          (*w0B)=fake_sampling(jackboot,2.12997 ,0.00157 ,jack_tot,1234);// petros-from jacob 31/07/20 fit  M_PS^2/f_PS^2
          (*w0C)=fake_sampling(jackboot,2.50451, 0.00172,jack_tot,12345);//  petros-from jacob 31/07/20 fit  M_PS^2/f_PS^2
          (*scalefm)=fake_sampling(jackboot,0.17383, 0.00063,jack_tot,1);
          
      }
      else if(strcmp(scaletype,"t0_w0")==0){
          (*w0A)=fake_sampling(jackboot,1.33590,0.00120,jack_tot,123); //t/w0
          (*w0B)=fake_sampling(jackboot,1.52789  ,0.00033,jack_tot,1234);// t/w0
          (*w0C)=fake_sampling(jackboot,1.77671,0.00037,jack_tot,12345);
          (*scalefm)=fake_sampling(jackboot,0.11969, 0.00062,jack_tot,1);
     }
     else if(strcmp(scaletype,"sqrtt0")==0){
         (*w0A)=fake_sampling(jackboot,1.5662,0.00220,jack_tot,123); //t
         (*w0B)=fake_sampling(jackboot,1.80397 ,0.00068,jack_tot,1234);// t
         (*w0C)=fake_sampling(jackboot,2.10945,0.00082,jack_tot,12345);
         (*scalefm)=fake_sampling(jackboot,0.14436, 0.00061,jack_tot,1);
     }
     else {
         printf("error no scale selected\n"); exit(1);
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
         //(*ZpA)=fake_sampling(jackboot,0.477,0.002,jack_tot,321);//M2 2Gev  mattero report Feb/2021  
         //(*ZpB)=fake_sampling(jackboot,0.4781,0.003,jack_tot,3214);//M2 2Gev mattero report Feb/2021 
         //(*ZpC)=fake_sampling(jackboot,0.490,0.002,jack_tot,32145);//M2 2Gev  mattero report Feb/2021 
         
     }
     if(strcmp(Mtype,"M1b")==0){
         (*ZpA)=fake_sampling(jackboot,0.4748,0.0026,jack_tot,321);//M2 2Gev  petros report Oct/2020  
         (*ZpB)=fake_sampling(jackboot,0.4775,0.0027,jack_tot,3214);//M2 2Gev petros report Oct/2020
         (*ZpC)=fake_sampling(jackboot,0.4873,0.0027,jack_tot,32145);//M2 2Gev  petros report Oct/2020
         //(*ZpA)=fake_sampling(jackboot,0.492,0.002,jack_tot,321);//M2 2Gev  mattero report Feb/2021  
         //(*ZpB)=fake_sampling(jackboot,0.494,0.002,jack_tot,3214);//M2 2Gev mattero report Feb/2021
         //(*ZpC)=fake_sampling(jackboot,0.505,0.002,jack_tot,32145);//M2 2Gev   mattero report Feb/2021
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
         //(*ZpA)=fake_sampling(jackboot,0.507,0.002,jack_tot,321);//M2 2Gev   mattero report Feb/2021
         //(*ZpB)=fake_sampling(jackboot,0.504,0.002,jack_tot,3214);//M2 2Gev  mattero report Feb/2021
         //(*ZpC)=fake_sampling(jackboot,0.505,0.002,jack_tot,32145);//M2 2Gev  mattero report Feb/2021
         
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
         //(*ZpA)=fake_sampling(jackboot,0.530,0.002,jack_tot,321);//M2 2Gev  RIMOM  
         //(*ZpB)=fake_sampling(jackboot,0.526,0.002,jack_tot,3214);//M2 2Gev RIMOM
         //(*ZpC)=fake_sampling(jackboot,0.524,0.002,jack_tot,32145);//M2 2Gev  RIMOM
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
      (*w0)=(double**) malloc(sizeof(double*)*ensembles);
      if(ensembles>0) (*w0)[0]=fake_sampling(jackboot,1.75971827958,  0.0043117955,jack_tot,31+0);
      if(ensembles>1) (*w0)[1]=fake_sampling(jackboot,1.7766037849,  0.00326168171497332,jack_tot,31+1);
      if(ensembles>2) (*w0)[2]=fake_sampling(jackboot,1.7928047479464,  0.00166108579455367,jack_tot,31+2);
      if(ensembles>3) (*w0)[3]=fake_sampling(jackboot,1.82491342151081,  0.0032611047018913,jack_tot,31+3);
      if(ensembles>4) (*w0)[4]=fake_sampling(jackboot,2.0981815470160,  0.001880560896,jack_tot,31+4);
      if(ensembles>5) (*w0)[5]=fake_sampling(jackboot,2.1176064583,  0.0013130969466,jack_tot,31+5);
      if(ensembles>6) (*w0)[6]=fake_sampling(jackboot,2.1271529070228,  0.001943835318,jack_tot,31+6);
      if(ensembles>7) (*w0)[7]=fake_sampling(jackboot,2.50451148219,  0.0017199806,jack_tot,31+7);
     
      if(ensembles>8) (*w0)[8]=fake_sampling(jackboot,2.098181547016053,  0.0018805608969,jack_tot,31+8);
      if(ensembles>8) (*w0)[9]=fake_sampling(jackboot, 2.4682126       ,  0.0047903,jack_tot,31+9);
      


    
}



void init_Z( struct database_file_jack *jack_files, struct header *head,int jack_tot, struct data_jack **gJ, const char *scaletype, const char *Mtype){
      int j;
      double *w0A,*w0B,*w0C, *ZpA,*ZpB,*ZpC, **w0,*scalefm;
      create_fake_distribution(jack_files[0].sampling, &w0A, &w0B, &w0C, &ZpA, &ZpB, &ZpC,jack_tot,scaletype,Mtype, &w0, &scalefm);
      printf("ensembles=%d\n",ensembles);
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
      if (ensembles>9){
          (*gJ)[9].w0[j]=w0C[j];
          (*gJ)[9].Zp[j]=ZpC[j];
      } 
      
      
      if (strcmp(scaletype,"w0xensemble")==0){
          for(int e=0;e<ensembles;e++)
             (*gJ)[e].w0[j]=w0[e][j];
      }
      (*gJ)[0].scalefm[j]=scalefm[j];
    }
    free(w0A);free(w0B);free(w0C); free(ZpA);free(ZpB);free(ZpC);
    free_2(ensembles,w0);
    free(scalefm);
    
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
          gJ[e].scalefm=(double*) calloc(*jack_tot,sizeof(double));
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
          gJ[e].scalefm=(double*) malloc((*jack_tot)*sizeof(double));
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

void print_fit_K_info(char **argv,int jack_tot, double **fit, struct fit_type fit_info, double **phys_point, struct result_jack &r1, struct data_jack *grephJ, struct header *head , const char *AV,const char *namefit, const char *save,  std::string M, std::string GF , store_fit_clover mud, std::vector<store_fit_clover> &ms){
    
    int N=fit_info.N;
    error(M!=mud.M,1, "print_fit_K_info"," %s   != %s", M.c_str(),mud.M.c_str());
    char nametex[NAMESIZE];
    mysprintf(nametex,NAMESIZE,"%s/%s.tex",argv[2],namefit);
    
    printf("\n%s\n",namefit);
    FILE *f_fit=open_file(nametex,"w+");
    double *msMeV=(double*) malloc(sizeof(double)*jack_tot);
    double *ms_mud=(double*) malloc(sizeof(double)*jack_tot);
    double *fk_fpi=(double*) malloc(sizeof(double)*jack_tot);
    for (int j=0;j<jack_tot;j++){
        if (strcmp(AV,"K")==0)
            msMeV[j]=fit[0][j]/(mud.w0[j]/197.326963);
        else if (strcmp(AV,"pi/K")==0)
            msMeV[j]=fit[0][j]/(grephJ[0].scalefm[j]/197.326963);
            
        ms_mud[j]=msMeV[j]/(mud.jack_m[j]);
        if (N==0)
            fk_fpi[j]=0;
        if (N>1)
            fk_fpi[j]=fit[1][j];
        
        //fit[0][j]=fit[0][j]/r1.w0MeV[j];
        //fit[1][j]=fit[1][j]/r1.w0MeV[j];
    }
    double *Ci=mean_and_error(argv[1],jack_tot,msMeV);
    
    if (strcmp(AV,"K")==0)
        fprintf(f_fit,"Imposing $M_K =%.2f \\pm %.2g$ and $w_0=%.4f$ fm   $m_{ud}=%.4f \\pm %.2g$\n",v_MKMeV,err_MKMeV,mud.w0[jack_tot-1],mud.jack_m[jack_tot-1],
                error_jackboot( argv[1] ,jack_tot,mud.jack_m)        );
    else if (strcmp(AV,"pi/K")==0)
        fprintf(f_fit,"Imposing $M_K =%.2f \\pm %.2g$ and $w_0=%.4f$ fm \n",v_MKMeV,err_MKMeV,grephJ[0].scalefm[jack_tot-1]);
    
    fprintf(f_fit,"\\begin{gather}\n   m_{s}=(%g\\pm%.2g) MeV   \\\\ \n",Ci[0],Ci[1]);
    //printf("\\end{gather} \n");
    
    free(Ci);
    
    Ci=mean_and_error(argv[1],jack_tot,ms_mud);
    fprintf(f_fit,"  m_{s}/m_{ub}=(%g\\pm%.2g)   \n\\end{gather} \n",Ci[0],Ci[1]);
    free(Ci);
    
    if(N>1){
        Ci=mean_and_error(argv[1],jack_tot,fk_fpi);
        fprintf(f_fit,"\\begin{gather}\n f_{K}/f_{\\pi}=(%g\\pm%.2g)   \n\\end{gather} \n",Ci[0],Ci[1]);
        free(Ci);
    }
    
    fclose(f_fit);
    
    
    if(strcmp(save,"yes")==0){
        store_fit_clover  tmpm;
        tmpm.M=M;
        tmpm.GF=GF;
        tmpm.name=namefit;
        tmpm.jack_m=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.ms=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.ms_mud=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.fk=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.fk_fpi=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.chi2=(double*) malloc(sizeof(double)*jack_tot );
        for (int j=0;j<jack_tot;j++){
            tmpm.jack_m[j]=mud.jack_m[j];
            tmpm.ms[j]=msMeV[j];  
            tmpm.ms_mud[j]=ms_mud[j];
            if (N==1){
                tmpm.fk[j]=0;
                tmpm.fk_fpi[j]=0;
            }
            if (N>1){
                tmpm.fk[j]=fit[1][j]*result.fpiMeV_exp[j];
                tmpm.fk_fpi[j]=fit[1][j];
            }
            tmpm.chi2[j]=fit[N][j];
            //tmpm.chi2[j]=fit_out.chi2[j];
            
        }
        tmpm.Njack=jack_tot;
        
        ms.emplace_back(tmpm);
    }
    
    free(msMeV);free(ms_mud);
        free(fk_fpi);
    free_2(fit_info.N+1,fit);
}


void print_fit_Ds_info(char **argv,int jack_tot, double **fit, struct fit_type fit_info, double **phys_point, struct result_jack &r1, struct data_jack *grephJ, struct header *head , const char *AV,const char *namefit, const char *save,  std::string M, std::string GF , store_fit_clover mud,store_fit_clover ms,store_fit_clover mc, std::vector<store_fit_clover> &mc_Ds){
    
    int N=fit_info.N;
    error(M!=mud.M,1, "print_fit_Ds_info"," %s   != %s", M.c_str(),mud.M.c_str());
    char nametex[NAMESIZE];
    mysprintf(nametex,NAMESIZE,"%s/%s.tex",argv[2],namefit);
    
    printf("\n%s\n",namefit);
    FILE *f_fit=open_file(nametex,"w+");
    double *mcMeV=(double*) malloc(sizeof(double)*jack_tot);
    double *mc_ms=(double*) malloc(sizeof(double)*jack_tot);
    double *fDs=(double*) malloc(sizeof(double)*jack_tot);
    for (int j=0;j<jack_tot;j++){
        
        mcMeV[j]=fit[0][j]/(mud.w0[j]/197.326963);
        
        mc_ms[j]=mcMeV[j]/(ms.ms[j]);
        mcMeV[j]*=m_from2to3GEV;
        
        if (N>1)
        fDs[j]=fit[1][j]/(mud.w0[j]/197.326963);
        
        //fit[0][j]=fit[0][j]/r1.w0MeV[j];
        //fit[1][j]=fit[1][j]/r1.w0MeV[j];
    }
    double *Ci=mean_and_error(argv[1],jack_tot,mcMeV);
    //printf("grep mc: %g     GF=%g      m_ud=%g     m_s=%g    fDs=%g   \n",mcMeV[jack_tot-1], mud.w0[jack_tot-1], mud.jack_m[jack_tot-1], ms.ms[jack_tot-1], fDs[jack_tot-1] );
    
    fprintf(f_fit,"Imposing $M_D =%.2f \\pm %.2g$ and $w_0=%.4f$ fm (%f MeV)\n",v_MDMeV,err_MDMeV,mud.w0[jack_tot-1],mud.w0[jack_tot-1]/197.326963);
    fprintf(f_fit,"\\begin{gather}\n   m_{c}=(%g\\pm%.2g) MeV   \\\\ \n",Ci[0],Ci[1]);
    //printf("\\end{gather} \n");
    
    free(Ci);
    
    Ci=mean_and_error(argv[1],jack_tot,mc_ms);
    
    fprintf(f_fit," m_{c}=(%g\\pm%.2g)  \n\\end{gather} \n",Ci[0],Ci[1]);
    free(Ci);
    
    
    
    if (N>1){
        Ci=mean_and_error(argv[1],jack_tot,fDs);
        
        fprintf(f_fit,"\\begin{gather}\n f_{D_s}=(%g\\pm%.2g)  \n\\end{gather} \n",Ci[0],Ci[1]);
        free(Ci);
    }
    fclose(f_fit);
    
    
    if(strcmp(save,"yes")==0){
        store_fit_clover  tmpm;
        tmpm.M=M;
        tmpm.GF=GF;
        tmpm.name=namefit;
        //tmpm.jack_m=(double*) malloc(sizeof(double)*jack_tot );
        //tmpm.w0=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.ms=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.mc=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.mc_ms=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.fDs=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.fDs_fD=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.chi2=(double*) malloc(sizeof(double)*jack_tot );
        for (int j=0;j<jack_tot;j++){
            //tmpm.jack_m[j]=mud.jack_m[j];
            //tmpm.w0[j]=mud.w0[j];
            tmpm.mc[j]=mcMeV[j];
            tmpm.ms[j]=ms.ms[j];  
            tmpm.mc_ms[j]=mc_ms[j];
            if (N==1){
                tmpm.fDs[j]=0;
                tmpm.fDs_fD[j]=0;
            }
            else if (N>1){
                tmpm.fDs[j]=fDs[j];
                tmpm.fDs_fD[j]=fDs[j]/mc.fD[j];
            }
            tmpm.chi2[j]=fit[N][j];
            
        }
        tmpm.Njack=jack_tot;
        
        mc_Ds.emplace_back(tmpm);
    }
    
    free(mcMeV);free(mc_ms);free(fDs);
    free_2(fit_info.N,fit);
}

void print_fit_D_info(char **argv,int jack_tot, double **fit, struct fit_type fit_info, double **phys_point, struct result_jack &r1, struct data_jack *grephJ, struct header *head , const char *AV,const char *namefit, const char *save,  std::string M, std::string GF , store_fit_clover mud,store_fit_clover ms, std::vector<store_fit_clover> &mc){

    int N=fit_info.N;
    error(M!=mud.M,1, "print_fit_D_info"," %s   != %s", M.c_str(),mud.M.c_str());
    char nametex[NAMESIZE];
    mysprintf(nametex,NAMESIZE,"%s/%s.tex",argv[2],namefit);
    
    printf("\n%s\n",namefit);
    FILE *f_fit=open_file(nametex,"w+");
    double *mcMeV=(double*) malloc(sizeof(double)*jack_tot);
    double *mc_ms=(double*) malloc(sizeof(double)*jack_tot);
    double *fD=(double*) malloc(sizeof(double)*jack_tot);
    for (int j=0;j<jack_tot;j++){
        
        mcMeV[j]=fit[0][j]/(mud.w0[j]/197.326963);
        
        mc_ms[j]=mcMeV[j]/(ms.ms[j]);
        
        mcMeV[j]*=m_from2to3GEV;
        if (N==1)
            fD[j]=0;
        else if (N>1)
            fD[j]=fit[1][j]/(mud.w0[j]/197.326963);
        
        
        //fit[0][j]=fit[0][j]/r1.w0MeV[j];
        //fit[1][j]=fit[1][j]/r1.w0MeV[j];
    }
    double *Ci=mean_and_error(argv[1],jack_tot,mcMeV);
    printf("grep mc: %g     GF=%g      m_ud=%g     m_s=%g\n",mcMeV[jack_tot-1], mud.w0[jack_tot-1], mud.jack_m[jack_tot-1], ms.ms[jack_tot-1] );
    
    fprintf(f_fit,"Imposing $M_D =%.2f \\pm %.2g$ and $w_0=%.4f$ fm (%f MeV)\n",v_MDMeV,err_MDMeV,mud.w0[jack_tot-1],mud.w0[jack_tot-1]/197.326963);
    fprintf(f_fit,"\\begin{gather}\n   m_{c}=(%g\\pm%.2g) MeV   \\\\ \n",Ci[0],Ci[1]);
    //printf("\\end{gather} \n");
    
    free(Ci);
    
    Ci=mean_and_error(argv[1],jack_tot,mc_ms);
    fprintf(f_fit," m_{c}=(%g\\pm%.2g)  \n\\end{gather} \n",Ci[0],Ci[1]);
    
    if (N>1){
        Ci=mean_and_error(argv[1],jack_tot,fD);
        fprintf(f_fit," \\begin{gather}\n \n  f_{D}=(%g\\pm%.2g)   \n\\end{gather} \n",Ci[0],Ci[1]);
    }
    fclose(f_fit);
    free(Ci);
    if(strcmp(save,"yes")==0){
        store_fit_clover  tmpm;
        tmpm.M=M;
        tmpm.GF=GF;
        tmpm.name=namefit;
        //tmpm.jack_m=(double*) malloc(sizeof(double)*jack_tot );
        //tmpm.w0=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.ms=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.mc=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.mc_ms=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.fD=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.fD_fk=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.chi2=(double*) malloc(sizeof(double)*jack_tot );
        for (int j=0;j<jack_tot;j++){
            //tmpm.jack_m[j]=mud.jack_m[j];
            //tmpm.w0[j]=mud.w0[j];
            tmpm.mc[j]=mcMeV[j];
            tmpm.ms[j]=ms.ms[j];  
            tmpm.mc_ms[j]=mc_ms[j];
            if (N==1){
                tmpm.fD[j]=0;
                tmpm.fD_fk[j]=0;
            }
            if (N>1){
                tmpm.fD[j]=fD[j];
                tmpm.fD_fk[j]=fD[j]/ms.fk[j];
            }
            tmpm.chi2[j]=fit[N][j];
            //tmpm.chi2[j]=fit_out.chi2[j];
            
        }
        tmpm.Njack=jack_tot;
        
        mc.emplace_back(tmpm);
    }
    
    free(mcMeV);free(mc_ms);free(fD);
    free_2(fit_info.N+1,fit);
}

void  print_fit_info(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, struct result_jack &r1, struct data_jack *grephJ, struct header *head , const char *AV,const char *namefile, const char *save,  std::string M, std::string GF , std::vector<store_fit_clover> &mud){
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
        /* if (strcmp(namefile,"fit_Mpi_Fpi")==0 ||  strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M1")==0  ||  strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M2a")==0 ||  strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M2b")==0 ||strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M1a")==0 || strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M1b")==0  ||
           strcmp(namefile,"fit_FpiMpi4_GL_w0_M1a")==0        ){*/
             if(i==0)     fprintf(ftex,"& Bw_{0}= %+.5f \\pm \t%.2g  ",Ci[i][0],Ci[i][1]);
             else if(i==1)     fprintf(ftex,"& fw_{0}= %+.5f \\pm \t%.2g  ",Ci[i][0],Ci[i][1]);
             else if(i==2)     fprintf(ftex,"& \\bar{\\ell_3}= %+.5f \\pm \t%.2g   ",Ci[i][0],Ci[i][1]);
             else if(i==3)     fprintf(ftex,"& P_2= %+.5f \\pm \t%.2g   ",Ci[i][0],Ci[i][1]);
             else if(i==4)     fprintf(ftex,"& \\bar{\\ell_4}= %+.5f \\pm \t%.2g   ",Ci[i][0],Ci[i][1]);
             else if(i==5)     fprintf(ftex,"& P_4= %+.5f \\pm \t%.2g   ",Ci[i][0],Ci[i][1]);
             else
                fprintf(ftex,"& P_{%d}= %+.5f \\pm \t%.2g   ",i,Ci[i][0],Ci[i][1]);
             
             if (i== fit_info.Npar-1)
                 fprintf(ftex," \n");
             else
                 fprintf(ftex,"\\\\ \n");
       /* }
        else
                fprintf(ftex,"& P_{%d}= %+.5f \\pm \t%.2g   \\\\ \n",i,Ci[i][0],Ci[i][1]);
        */
       
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
    */
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
        x[j][2]=r1.MpiMeV[j]*r1.MpiMeV[j];//Mpi2
        x[j][3]=r1.fpiMeV_exp[j];//fpi
        x[j][4]=1e+10;//L such that L/w=1e+4
        
    }
    
    for (j=0;j<jack_tot;j++){
        if(strcmp(AV,"pion")==0){
       
       in=r1.MpiMeV[j]*r1.MpiMeV[j]/(r1.fpiMeV_exp[j]*r1.fpiMeV_exp[j]);
       xi[j]=rtbis(Mw2_over_fw2_chiral_FVE_a0_minus,in,fit_info.Npar,tif[j], 0.0001, 0.01, 1e-10);//gives mw
       phys_point[j][0]=xi[j];
       r1.fpiw[j]=fPSw_chiral_FVE(  xi[j],fit_info.Npar,tif[j]);
       w0_estimate[j]=r1.fpiw[j]/(v_fpiMeV_exp/197.326963);
      /* if(strcmp(namefile,"fit_FpiMpi4_GL_w0_M1a")==0 ){
         //  xi[j]=rtbis_func_eq_input(fit_info.function ,
           //                          0,//double n
            //                         fit_info.Nvar, x,fit_info.Npar, tif[j], 
             //                        0,//ivar
              //                       in, 0.0001, 0.01, 1e-10);//gives mw  
           
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
       */
       r1.w0fm[j]= w0_estimate[j] ;
       r1.w0MeV[j]= w0_estimate[j]/197.326963 ;
       r1.mlw[j]=xi[j];
       xi[j]=xi[j]/(r1.w0MeV[j]);
       r1.Bw[j]=fit_out.P[0][j];
       r1.fw[j]=fit_out.P[1][j];
       
       
       tmp3[j]=r1.fpiMeV_exp[j]/(r1.fw[j] /r1.w0MeV[j]);
       
       
        }
       else if(strcmp(AV,"pionM")==0){
           
           in=r1.MpiMeV[j]*r1.MpiMeV[j]*r1.MpiMeV[j]*r1.MpiMeV[j]*r1.fpiMeV_exp[j];
           w0_estimate[j]=rtbis_func_eq_input(FpiMpi4_fromM, 0 , fit_info.Nvar, x[j] ,fit_info.Npar, tif[j], 1 , in, 0, 1, 1e-7);//gives w0MeV
           in=r1.MpiMeV[j]*r1.MpiMeV[j]/(r1.fpiMeV_exp[j]*r1.fpiMeV_exp[j]);
           xi[j]=rtbis(Mw2_over_fw2_chiral_FVE_a0_minus,in,fit_info.Npar,tif[j], 0.0001, 0.01, 1e-10);//gives mw
           xi[j]=xi[j]/(w0_estimate[j]);
           
           w0_estimate[j]*=197.326963; //convert w0 in fm
           
           
       }
       
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
   
    if(strcmp(save,"yes")==0){
        store_fit_clover  tmpm;
        tmpm.M=M;
        tmpm.GF=GF;
        tmpm.name=namefile;
        tmpm.jack_m=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.jack_B=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.jack_f=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.chi2=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.w0=(double*) malloc(sizeof(double)*jack_tot );
        tmpm.Sigma_13=(double*) malloc(sizeof(double)*jack_tot );
        for ( j=0;j<jack_tot;j++){
            tmpm.jack_m[j]=xi[j];
            tmpm.jack_B[j]=fit_out.P[0][j]/(w0_estimate[j]/197.326963);
            tmpm.jack_f[j]=fit_out.P[1][j]/(w0_estimate[j]/197.326963);
            tmpm.chi2[j]=fit_out.chi2[j];
            tmpm.w0[j]=w0_estimate[j];
            
            tmpm.Sigma_13[j]=pow(  tmpm.jack_B[j]*tmpm.jack_f[j]*tmpm.jack_f[j]/2 ,  1./3.);
            
        }
        tmpm.Njack=jack_tot;
        
        mud.emplace_back(tmpm);
    }
   /////// lattice spacing /// 
   
   div_jackboot(jack_tot , tmp,w0_estimate, grephJ[0].w0 );
   fprintf(ftex,"a=(%g\\pm%.2g)  \\\\ \n ",tmp[jack_tot-1],error_jackboot(argv[1], jack_tot,tmp ));
   div_jackboot(jack_tot , tmp,w0_estimate, grephJ[6].w0 );
   fprintf(ftex,"(%g\\pm%.2g) \\\\ \n",tmp[jack_tot-1],error_jackboot(argv[1], jack_tot,tmp ));
   div_jackboot(jack_tot , tmp,w0_estimate, grephJ[7].w0 );
   fprintf(ftex,"(%g\\pm%.2g) fm \\\\\n ",tmp[jack_tot-1],error_jackboot(argv[1], jack_tot,tmp ));
   free(tmp);
    
    free(C1[1]);
    C1[1]=mean_and_error(argv[1],jack_tot,tmp3);
    fprintf(ftex,"f_\\pi/f=(%g\\pm%.2g)    \n\\end{gather}\n",C1[1][0],C1[1][1]);
    printf("f_\\pi/f=(%g\\pm%.2g)    \n\\end{gather}\n",C1[1][0],C1[1][1]);
    
    free(tmp3);

    free_2(2,C1);
    free(w0_estimate);
    
    ////////////////////////print fit end
    //using w0 from silvano
    if (strcmp(AV,"pion")==0){
        for (j=0;j<jack_tot;j++){
            in=r1.MpiMeV[j]*r1.MpiMeV[j]/(r1.fpiMeV_exp[j]*r1.fpiMeV_exp[j]);
            xi[j]=rtbis(Mw2_over_fw2_chiral_FVE_a0_minus,in,fit_info.Npar,tif[j], 0.0001, 0.01, 1e-10);//gives mw
            xi[j]=xi[j]/(gjack[0].scalefm[j]/197.326963);
        }
        fprintf(ftex,"Using scale $ =%g \\pm %g$ fm\n",gjack[0].scalefm[jack_tot-1],error_jackboot(argv[1],jack_tot,gjack[0].scalefm));
        printf("Using $scale  =%g \\pm %g$ fm\n",gjack[0].scalefm[jack_tot-1],error_jackboot(argv[1],jack_tot,gjack[0].scalefm));
        
        fprintf(ftex,"\\begin{gather}\n   m_{ud}=(%g\\pm%.2g) MeV   \n \\end{gather}\n",xi[jack_tot-1],error_jackboot(argv[1], jack_tot, xi));
    }
    
    
    
    free(xi);
    
    
    ////////
    
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
    fclose(ftex);fclose(fgp);free(chi2m); free(Ci);
    
}


double **init_phys_point(int jack_tot){
    int j;
    double **phys_point=(double**) malloc(sizeof(double*)*jack_tot);

    for(j=0;j<jack_tot;j++){
        phys_point[j]=(double*) malloc(sizeof(double)*13);
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


void compute_systematic_D( char **argv, std::vector<store_fit_clover>  mc, std::string suffix=""){
    char name[NAMESIZE];
    mysprintf(name,NAMESIZE,"%s/systematics_mc%s.tex", argv[2],suffix.c_str());
    
    
    FILE *ftex=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/systematics_mc%s.txt", argv[2],suffix.c_str());
    FILE *fdata=open_file(name,"w+");
    
    int Njack=mc[0].Njack;
    double ave_ms=0,ave_ms_mc=0,ave_fk=0, ave_fk_fpi=0;
    double sigma_ms=0, sigma_ms_mc=0,sigma_fk=0,sigma_fk_fpi=0;
    int Nobs=4;
    
    std::vector<double>  sigma(Nobs);
    std::vector<double>  ave(Nobs);
    std::vector<int>     Ngood;
    std::vector<double>  stat(Nobs);
    std::vector<double>  sist(Nobs);
    
    for(int i=0;i<Nobs; i++){
        ave[i]=0;
        sigma[i]=0;
    }
    
    for (int i=0; i<  mc.size(); i++){
        //if (mc[i].chi2[Njack-1]<2){
            ave[0]+=mc[i].mc[Njack-1];
            ave[1]+=mc[i].mc_ms[Njack-1];
            
            ave[2]+=mc[i].fD[Njack-1];
            ave[3]+=mc[i].fD_fk[Njack-1];
            
            //error_i.emplace_back(error_jackboot(argv[1],Njack,mc[i].jack_m));
            
            sigma[0]+= pow(error_jackboot(argv[1],Njack,mc[i].mc),2);
            sigma[1]+=pow(error_jackboot(argv[1],Njack,mc[i].mc_ms),2);
            sigma[2]+=pow(error_jackboot(argv[1],Njack,mc[i].fD),2);
            sigma[3]+=pow(error_jackboot(argv[1],Njack,mc[i].fD_fk),2);
            Ngood.emplace_back(i);
        //}
    }
    for(int i=0;i<Nobs; i++){
        ave[i]/=Ngood.size();
        stat[i]=sigma[i]/Ngood.size();
        
    }
    
    for (int i: Ngood){
        
        sigma[0]+=pow(  ave[0]-mc[i].mc[Njack-1]  ,2);
        sigma[1]+=pow(  ave[1]-mc[i].mc_ms[Njack-1]  ,2);
        sigma[2]+=pow(  ave[2]-mc[i].fD[Njack-1]  ,2);
        sigma[3]+=pow(  ave[3]-mc[i].fD_fk[Njack-1]  ,2);
        
    }
    
    for(int i=0;i<Nobs; i++){
        sigma[i]/=Ngood.size();
        sist[i]=sigma[i]-stat[i];
        sigma[i]=sqrt(sigma[i]);
        stat[i]=sqrt(stat[i]);
        sist[i]=sqrt(sist[i]);
        
    }
    
    
    //fprintf(ftex,"\\begin{tabular}{|c|c|c|c|c|}  \n\\hline\n");
    //fprintf(ftex," $Z_p$ &  scale  &  fit & $m_ud$ [MeV]    & err   \\\\ \n\\hline \n");
    for (int i=0; i<  mc.size(); i++){
        //fprintf(ftex," %s  & %s & %s & %g &  %g    \\\\ \n\\hline \n", mc[i].M.c_str(), mc[i].GF.c_str(), mc[i].name.c_str(), mc[i].jack_m[Njack-1],  error_jackboot(argv[1],Njack,mc[i].jack_m) );
        
        fprintf(fdata,"%s  %s  %s \t", mc[i].M.c_str(), mc[i].GF.c_str(), mc[i].name.c_str());
        
        fprintf(fdata,"%g  %g %g  %g   \t", mc[i].mc[Njack-1], error_jackboot(argv[1],Njack,mc[i].mc),ave[0],sigma[0]);
        fprintf(fdata,"%g  %g %g  %g   \t", mc[i].mc_ms[Njack-1], error_jackboot(argv[1],Njack,mc[i].mc_ms),ave[1],sigma[1]);
        fprintf(fdata,"%g  %g %g  %g   \t", mc[i].fD[Njack-1], error_jackboot(argv[1],Njack,mc[i].fD),ave[2],sigma[2]);
        fprintf(fdata,"%g  %g %g  %g   \t", mc[i].fD_fk[Njack-1], error_jackboot(argv[1],Njack,mc[i].fD_fk),ave[3],sigma[3]);
        fprintf(fdata,"%g     \t", mc[i].chi2[Njack-1]);
        fprintf(fdata,"%g  %g \t %g  %g \t %g  %g \t %g  %g\n", stat[0],sist[0],stat[1],sist[1],stat[2],sist[2],stat[3],sist[3]);
        
        free(mc[i].mc);
        //free(mc[i].jack_m);
        free(mc[i].mc_ms);
        free(mc[i].fD);
        free(mc[i].fD_fk);
        free(mc[i].chi2);
        
    }
    
    //fprintf(ftex,"\\end{tabular}\n");    
    fprintf(ftex, "final value\n");
    fprintf(ftex, "\\begin{gather}\n  m_{c}=%g   \\pm  %g  \\,, \\\\ \n", ave[0],sigma[0]);
    fprintf(ftex, "  m_{c}/m_{s}=%g \\pm  %g  \\,, \\\\ \n", ave[1],sigma[1]);
    fprintf(ftex, "  f_D=%g   \\pm  %g \\,, \\\\ \n", ave[2],sigma[2]);
    fprintf(ftex, "  f_{D}/f_{K}=%g   \\pm  %g  \n \\end{gather}\n", ave[3],sigma[3]);
    
}


void compute_systematic_Ds( char **argv, std::vector<store_fit_clover>  mc, std::string suffix=""){
    char name[NAMESIZE];
    mysprintf(name,NAMESIZE,"%s/systematics_mc_Ds%s.tex", argv[2],suffix.c_str());
    
    FILE *ftex=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/systematics_mc_Ds%s.txt", argv[2],suffix.c_str());
    FILE *fdata=open_file(name,"w+");
    
    int Njack=mc[0].Njack;
    double ave_ms=0,ave_ms_mc=0,ave_fk=0, ave_fk_fpi=0;
    double sigma_ms=0, sigma_ms_mc=0,sigma_fk=0,sigma_fk_fpi=0;
    int Nobs=4;
    
    std::vector<double>  sigma(Nobs);
    std::vector<double>  ave(Nobs);
    std::vector<double>  stat(Nobs);
    std::vector<double>  sist(Nobs);
    
    std::vector<int>     Ngood;
    
    for(int i=0;i<Nobs; i++){
        ave[i]=0;
        sigma[i]=0;
    }
    
    for (int i=0; i<  mc.size(); i++){
        //if (mc[i].chi2[Njack-1]<2){
            ave[0]+=mc[i].mc[Njack-1];
            ave[1]+=mc[i].mc_ms[Njack-1];
            ave[2]+=mc[i].fDs[Njack-1];
            ave[3]+=mc[i].fDs_fD[Njack-1];
            //error_i.emplace_back(error_jackboot(argv[1],Njack,mc[i].jack_m));
            
            sigma[0]+= pow(error_jackboot(argv[1],Njack,mc[i].mc),2);
            sigma[1]+=pow(error_jackboot(argv[1],Njack,mc[i].mc_ms),2);
            sigma[2]+=pow(error_jackboot(argv[1],Njack,mc[i].fDs),2);
            sigma[3]+=pow(error_jackboot(argv[1],Njack,mc[i].fDs_fD),2);
            Ngood.emplace_back(i);
        //}
        
    }
    for(int i=0;i<Nobs; i++){
        ave[i]/=Ngood.size();
        stat[i]=sigma[i]/Ngood.size();
        
    }
    
    for (int i: Ngood){
        
        sigma[0]+=pow(  ave[0]-mc[i].mc[Njack-1]  ,2);
        sigma[1]+=pow(  ave[1]-mc[i].mc_ms[Njack-1]  ,2);
        sigma[2]+=pow(  ave[2]-mc[i].fDs[Njack-1]  ,2);
        sigma[3]+=pow(  ave[3]-mc[i].fDs_fD[Njack-1]  ,2);
        
    }
    
    for(int i=0;i<Nobs; i++){
        sigma[i]/=Ngood.size();
        sist[i]=sigma[i]-stat[i];
        sigma[i]=sqrt(sigma[i]);
        stat[i]=sqrt(stat[i]);
        sist[i]=sqrt(sist[i]);
        
    }
    
    
    //fprintf(ftex,"\\begin{tabular}{|c|c|c|c|c|}  \n\\hline\n");
    //fprintf(ftex," $Z_p$ &  scale  &  fit & $m_ud$ [MeV]    & err   \\\\ \n\\hline \n");
    for (int i=0; i<  mc.size(); i++){
        //fprintf(ftex," %s  & %s & %s & %g &  %g    \\\\ \n\\hline \n", mc[i].M.c_str(), mc[i].GF.c_str(), mc[i].name.c_str(), mc[i].jack_m[Njack-1],  error_jackboot(argv[1],Njack,mc[i].jack_m) );
        
        fprintf(fdata,"%s  %s  %s \t", mc[i].M.c_str(), mc[i].GF.c_str(), mc[i].name.c_str());
        
        fprintf(fdata,"%g  %g %g  %g   \t", mc[i].mc[Njack-1], error_jackboot(argv[1],Njack,mc[i].mc),ave[0],sigma[0]);
        fprintf(fdata,"%g  %g %g  %g   \t", mc[i].mc_ms[Njack-1], error_jackboot(argv[1],Njack,mc[i].mc_ms),ave[1],sigma[1]);
        fprintf(fdata,"%g  %g %g  %g   \t", mc[i].fDs[Njack-1], error_jackboot(argv[1],Njack,mc[i].fDs),ave[2],sigma[2]);
        fprintf(fdata,"%g  %g %g  %g   \t", mc[i].fDs_fD[Njack-1], error_jackboot(argv[1],Njack,mc[i].fDs_fD),ave[3],sigma[3]);
        fprintf(fdata,"%g \t",mc[i].chi2[Njack-1]);
        fprintf(fdata,"%g  %g \t %g  %g \t %g  %g \t %g  %g\n", stat[0],sist[0],stat[1],sist[1],stat[2],sist[2],stat[3],sist[3]);
        
        
        free(mc[i].mc);
        //free(mc[i].jack_m);
        free(mc[i].mc_ms);
        free(mc[i].fDs);
        free(mc[i].fDs_fD);
        free(mc[i].chi2);
        
    }
    
    //fprintf(ftex,"\\end{tabular}\n");    
    fprintf(ftex, "final value\n");
    fprintf(ftex, "\\begin{gather}\n  m_{c}=%g   \\pm  %g  \\,, \\\\ \n", ave[0],sigma[0]);
    fprintf(ftex, "  m_{c}/m_{s}=%g \\pm  %g  \\,, \\\\ \n", ave[1],sigma[1]);
    fprintf(ftex, "  f_{D_s}=%g   \\pm  %g \\,, \\\\ \n", ave[2],sigma[2]);
    fprintf(ftex, "  f_{D_s}/f_{D}=%g   \\pm  %g  \n \\end{gather}\n", ave[3],sigma[3]);
    
}



void compute_systematic_K( char **argv, std::vector<store_fit_clover>  mud, std::string suffix=""){
    char name[NAMESIZE];
    mysprintf(name,NAMESIZE,"%s/systematics_ms%s.tex", argv[2],suffix.c_str());
    
    FILE *ftex=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/systematics_ms%s.txt", argv[2],suffix.c_str());
    FILE *fdata=open_file(name,"w+");
    
    int Njack=mud[0].Njack;
    double ave_ms=0,ave_ms_mud=0,ave_fk=0, ave_fk_fpi=0;
    double sigma_ms=0, sigma_ms_mud=0,sigma_fk=0,sigma_fk_fpi=0;
    int Nobs=4;

    std::vector<double>  sigma(Nobs);
    std::vector<double>  stat(Nobs);
    std::vector<double>  sist(Nobs);
    std::vector<double>  ave(Nobs);
    std::vector<int> Ngood;
    
    for(int i=0;i<Nobs; i++){
        ave[i]=0;
        sigma[i]=0;
    }
    
    for (int i=0; i<  mud.size(); i++){
        
        //if (mud[i].chi2[Njack-1]<2){
            ave[0]+=mud[i].ms[Njack-1];
            ave[1]+=mud[i].ms_mud[Njack-1];
            ave[2]+=mud[i].fk[Njack-1];
            ave[3]+=mud[i].fk_fpi[Njack-1];
            //error_i.emplace_back(error_jackboot(argv[1],Njack,mud[i].jack_m));
            
            sigma[0]+= pow(error_jackboot(argv[1],Njack,mud[i].ms),2);
            sigma[1]+=pow(error_jackboot(argv[1],Njack,mud[i].ms_mud),2);
            sigma[2]+=pow(error_jackboot(argv[1],Njack,mud[i].fk),2);
            sigma[3]+=pow(error_jackboot(argv[1],Njack,mud[i].fk_fpi),2);
            Ngood.emplace_back(i);
        //}
        
    }
    for(int i=0;i<Nobs; i++){
        ave[i]/=Ngood.size();
        stat[i]=sigma[i]/Ngood.size();
    }
     
    for (int i:Ngood){
        
        sigma[0]+=pow(  ave[0]-mud[i].ms[Njack-1]  ,2);
        sigma[1]+=pow(  ave[1]-mud[i].ms_mud[Njack-1]  ,2);
        sigma[2]+=pow(  ave[2]-mud[i].fk[Njack-1]  ,2);
        sigma[3]+=pow(  ave[3]-mud[i].fk_fpi[Njack-1]  ,2);
        
    }
    
    for(int i=0;i<Nobs; i++){
        
        sigma[i]/=Ngood.size();
        sist[i]=sigma[i]-stat[i];
        sigma[i]=sqrt(sigma[i]);
        stat[i]=sqrt(stat[i]);
        sist[i]=sqrt(sist[i]);
    }
    
    
    //fprintf(ftex,"\\begin{tabular}{|c|c|c|c|c|}  \n\\hline\n");
    //fprintf(ftex," $Z_p$ &  scale  &  fit & $m_ud$ [MeV]    & err   \\\\ \n\\hline \n");
    for (int i=0; i<  mud.size(); i++){
        //fprintf(ftex," %s  & %s & %s & %g &  %g    \\\\ \n\\hline \n", mud[i].M.c_str(), mud[i].GF.c_str(), mud[i].name.c_str(), mud[i].jack_m[Njack-1],  error_jackboot(argv[1],Njack,mud[i].jack_m) );
        
        fprintf(fdata,"%s  %s  %s \t", mud[i].M.c_str(), mud[i].GF.c_str(), mud[i].name.c_str());
        
        fprintf(fdata,"%g  %g %g  %g   \t", mud[i].ms[Njack-1], error_jackboot(argv[1],Njack,mud[i].ms),ave[0],sigma[0]);
        fprintf(fdata,"%g  %g %g  %g   \t", mud[i].ms_mud[Njack-1], error_jackboot(argv[1],Njack,mud[i].ms_mud),ave[1],sigma[1]);
        fprintf(fdata,"%g  %g %g  %g   \t", mud[i].fk[Njack-1], error_jackboot(argv[1],Njack,mud[i].fk),ave[2],sigma[2]);
        fprintf(fdata,"%g  %g %g  %g   \t", mud[i].fk_fpi[Njack-1], error_jackboot(argv[1],Njack,mud[i].fk_fpi),ave[3],sigma[3]);
        fprintf(fdata,"%g\t", mud[i].chi2[Njack-1]);
        
        fprintf(fdata,"%g  %g \t %g  %g \t %g  %g \t %g  %g\n", stat[0],sist[0],stat[1],sist[1],stat[2],sist[2],stat[3],sist[3]);
        
        free(mud[i].ms);free(mud[i].jack_m);
        free(mud[i].ms_mud);
        free(mud[i].fk);
        free(mud[i].fk_fpi);
        free(mud[i].chi2);
        
    }
    
    //fprintf(ftex,"\\end{tabular}\n");    
    fprintf(ftex, "final value\n");
    fprintf(ftex, "\\begin{equation}\n  m_{s}=%g   \\pm  %g  \n \\end{equation}\n", ave[0],sigma[0]);
    fprintf(ftex, "\\begin{equation}\n  m_{s}/m_{ud}=%g   \\pm  %g  \n \\end{equation}\n", ave[1],sigma[1]);
    fprintf(ftex, "\\begin{equation}\n  f_K=%g   \\pm  %g  \n \\end{equation}\n", ave[2],sigma[2]);
    fprintf(ftex, "\\begin{equation}\n  f_{K}/f_{\\pi}=%g   \\pm  %g  \n \\end{equation}\n", ave[3],sigma[3]);
    
}



void compute_systematic( char **argv, std::vector<store_fit_clover>  mud, std::string suffix=""){
    char name[NAMESIZE];
    mysprintf(name,NAMESIZE,"%s/systematics_mud%s.tex", argv[2],suffix.c_str());
    
    FILE *ftex=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/systematics_mud%s.txt", argv[2], suffix.c_str());
    FILE *fdata=open_file(name,"w+");
    
    int Njack=mud[0].Njack;
    double ave=0,ave_B=0,ave_f=0, ave_chi2=0;
    double ave_S13=0;
    double sigma=0, sigma_B=0,sigma_f=0;
    double sigma_S13=0;
    double stat=0, stat_B=0,stat_f=0, stat_S13=0;
    double sist=0, sist_B=0,sist_f=0, sist_S13=0;
    std::vector<double>  error_i;
    std::vector<int> Ngood;
    for (int i=0; i<  mud.size(); i++){
        error_i.emplace_back(error_jackboot(argv[1],Njack,mud[i].jack_m));
        //if (mud[i].chi2[Njack-1]<2){
            ave+=mud[i].jack_m[Njack-1];
            ave_B+=mud[i].jack_B[Njack-1];
            ave_f+=mud[i].jack_f[Njack-1];
            ave_chi2+=mud[i].chi2[Njack-1];
            ave_S13+=mud[i].Sigma_13[Njack-1];
            
            sigma+=error_i[i]*error_i[i];
            sigma_B+=error_jackboot(argv[1],Njack,mud[i].jack_B)*error_jackboot(argv[1],Njack,mud[i].jack_B);
            sigma_f+=error_jackboot(argv[1],Njack,mud[i].jack_f)*error_jackboot(argv[1],Njack,mud[i].jack_f);
            sigma_S13+=error_jackboot(argv[1],Njack,mud[i].Sigma_13)*error_jackboot(argv[1],Njack,mud[i].Sigma_13);
            
            stat+=error_i[i]*error_i[i];
            stat_B+=error_jackboot(argv[1],Njack,mud[i].jack_B)*error_jackboot(argv[1],Njack,mud[i].jack_B);
            stat_f+=error_jackboot(argv[1],Njack,mud[i].jack_f)*error_jackboot(argv[1],Njack,mud[i].jack_f);
            stat_S13+=error_jackboot(argv[1],Njack,mud[i].Sigma_13)*error_jackboot(argv[1],Njack,mud[i].Sigma_13);
            
            Ngood.emplace_back(i);
        //}
        
        
    }
    ave=ave/Ngood.size();
    ave_B=ave_B/Ngood.size();
    ave_f=ave_f/Ngood.size();
    ave_chi2=ave_chi2/Ngood.size();
    ave_S13=ave_S13/Ngood.size();
    for (int i:Ngood){
        sist+=  (ave-mud[i].jack_m[Njack-1] )* (ave-mud[i].jack_m[Njack-1] ) ;
        sist_B+=  (ave_B-mud[i].jack_B[Njack-1] )* (ave_B-mud[i].jack_B[Njack-1] ) ;
        sist_f+=  (ave_f-mud[i].jack_f[Njack-1] )* (ave_f-mud[i].jack_f[Njack-1] ) ;
        sist_S13+=  (ave_S13-mud[i].Sigma_13[Njack-1] )* (ave_S13-mud[i].Sigma_13[Njack-1] ) ;
        
        sigma+=  (ave-mud[i].jack_m[Njack-1] )* (ave-mud[i].jack_m[Njack-1] ) ;
        sigma_B+=  (ave_B-mud[i].jack_B[Njack-1] )* (ave_B-mud[i].jack_B[Njack-1] ) ;
        sigma_f+=  (ave_f-mud[i].jack_f[Njack-1] )* (ave_f-mud[i].jack_f[Njack-1] ) ;
        sigma_S13+=  (ave_S13-mud[i].Sigma_13[Njack-1] )* (ave_S13-mud[i].Sigma_13[Njack-1] ) ;
    }
    sigma=sigma/Ngood.size();
    sigma=sqrt(sigma);
    
    sigma_B=sigma_B/Ngood.size();
    sigma_B=sqrt(sigma_B);
    
    sigma_f=sigma_f/Ngood.size();
    sigma_f=sqrt(sigma_f);
    
    sigma_S13=sigma_S13/Ngood.size();
    sigma_S13=sqrt(sigma_S13);
    
    
    stat=stat/Ngood.size();
    stat=sqrt(stat);
    stat_B=stat_B/Ngood.size();
    stat_B=sqrt(stat_B);
    stat_f=sigma_f/Ngood.size();
    stat_f=sqrt(stat_f);
    stat_S13=sigma_S13/Ngood.size();
    stat_S13=sqrt(stat_S13);
    
    
    sist=sist/Ngood.size();
    sist=sqrt(sist);
    sist_B=sist_B/Ngood.size();
    sist_B=sqrt(sist_B);
    sist_f=sigma_f/Ngood.size();
    sist_f=sqrt(sist_f);
    sist_S13=sigma_S13/Ngood.size();
    sist_S13=sqrt(sist_S13);
    
    
    
    fprintf(ftex,"\\begin{tabular}{|c|c|c|c|c|}  \n\\hline\n");
    fprintf(ftex," $Z_p$ &  scale  &  fit & $m_ud$ [MeV]    & err   \\\\ \n\\hline \n");
    for (int i=0; i<  mud.size(); i++){
        fprintf(ftex," %s  & %s & %s & %g &  %g    \\\\ \n\\hline \n", mud[i].M.c_str(), mud[i].GF.c_str(), mud[i].name.c_str(), mud[i].jack_m[Njack-1],  error_i[i]);
        
        fprintf(fdata,"%s \t  %s \t %s  %g   %g  %g   %g\t", mud[i].M.c_str(), mud[i].GF.c_str(), mud[i].name.c_str(), mud[i].jack_m[Njack-1],  error_i[i], ave, sigma );
        
        fprintf(fdata,"%g  %g %g  %g   \t", mud[i].jack_B[Njack-1], error_jackboot(argv[1],Njack,mud[i].jack_B),ave_B,sigma_B);
        fprintf(fdata,"%g  %g %g  %g   \t", mud[i].jack_f[Njack-1], error_jackboot(argv[1],Njack,mud[i].jack_f),ave_f,sigma_f);
        fprintf(fdata,"%g    \t", mud[i].chi2[Njack-1]);
        fprintf(fdata,"%g  %g \t %g  %g \t %g  %g\t", stat,sist,stat_B,sist_B,stat_f,sist_f);
        fprintf(fdata,"%g  %g %g  %g   \t", mud[i].Sigma_13[Njack-1], error_jackboot(argv[1],Njack,mud[i].Sigma_13),ave_S13,sigma_S13);
        fprintf(fdata,"%g  %g \n",stat_S13,sist_S13);
        
        
        free(mud[i].jack_m);
        free(mud[i].jack_B);
        free(mud[i].jack_f);
        free(mud[i].chi2);
        free(mud[i].Sigma_13);
    }
    
    fprintf(ftex,"\\end{tabular}\n");    
    fprintf(ftex, "final value\n");
    fprintf(ftex, "\\begin{equation }\n  m_{ud}=%g   \\pm  %g  \n \\end{equation}\n", ave,sigma);
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
    
    for (e=0;e<ensembles;e++){
     tmp1=(double*) malloc(sizeof(double)* jack_tot); 
     //printf("L=%d\n",head[e].l1);
     for(j=0;j<jack_tot;j++){
         
            double xi=gjack[e].M_PS_jack[0][j]/(4*pi_greco*gjack[e].f_PS_jack[0][j]);
            xi=xi*xi;
            double delta=FVE_GL_Mpi(head[e].l1 , xi,   gjack[e].f_PS_jack[0][j]     );
            tmp1[j]= gjack[e].M_PS_jack[0][j]*(1+delta)/(gjack[e].f_PS_jack[0][j]* (1-0.25*delta));
            tmp1[j]=tmp1[j]*tmp1[j];
        }
        tmp=mean_and_error(argv[1],jack_tot,   tmp1 );
        printf("  M_PS^2/f_PS^2=%g  +-  %g  \n ",  tmp[0],tmp[1] );
        free(tmp);free(tmp1);
    }
          
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
    struct fit_all fit_chi2_good;
    fit_info.Nvar=13;
    
    double    **phys_point=init_phys_point(jack_tot);  
    
    printf("\n\n//////////////////////////////////////////////////////////////\n");
    printf("              OMEGA MK                             \n");
    printf("\n\n//////////////////////////////////////////////////////////////\n");
    //init_Z( jack_files, head, jack_tot, &gjack, "w0xensemble","M1");

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
    
    
    printf("\n\n//////////////////////////////////////////////////////////////\n");
    printf("              OMEGA MK   new                            \n");
    printf("\n\n//////////////////////////////////////////////////////////////\n");
    //init_Z( jack_files, head, jack_tot, &gjack, "w0xensemble","M1");

    Ci=(double**) malloc(sizeof(double*)*2);
    
    fit=fit_Omegaw0_from_M_new(jack_files, head , jack_tot, mass_index, gjack,  &result );
    
    tmp1=(double*) malloc(sizeof(double)*jack_tot);
    for (j=0;j<jack_tot;j++){
        tmp1[j]=fit[0][j]*197.326963;
    }

    Ci[0]=mean_and_error(argv[1],jack_tot,tmp1);
    printf("Imposing $M_\\K =%.2f \\pm %.2g$ , $M_\\pi=%.4f\\pm %2g$ fm  and M_\\Omega=%g \\pm %g\n",v_MKMeV,err_MKMeV,v_MpiMeV,err_MpiMeV,v_MOmegaMeV,err_MOmegaMeV);
    printf("   w_0=(%g\\pm%.2g) fm \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1]);
     
    free(Ci[0]);
    for (j=0;j<jack_tot;j++){
        tmp1[j]=fit[0][j]*197.326963 /gjack[0].w0[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,tmp1);
    printf("   a(A)=(%g\\pm%.2g) fm  \n",Ci[0][0],Ci[0][1]);
    free(Ci[0]);
    for (j=0;j<jack_tot;j++){
        tmp1[j]=fit[0][j]*197.326963 /gjack[4].w0[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,tmp1);
    printf("   a(B)=(%g\\pm%.2g) fm  \n",Ci[0][0],Ci[0][1]);
    free(Ci[0]);
    for (j=0;j<jack_tot;j++){
        tmp1[j]=fit[0][j]*197.326963 /gjack[7].w0[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,tmp1);
    printf("   a(C)=(%g\\pm%.2g) fm  \n",Ci[0][0],Ci[0][1]);
    

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
    std::vector<store_fit_clover>  notused;
   /* printf("\n\n///////////////////////////////////////Pion Mpi^2 and fpiMpi^4  GL  w0 M1a ///////////////////////\n");
    fit_info.Npar=8;
    fit_info.N=2;
    fit_info.function=fit_FpiMpi4_and_Mpi2_GL;
        
    init_Z( jack_files, head, jack_tot, &gjack, "w0","M1a");

        
    fit_out=fit_Mpi_fwMpi4_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info);
 
    //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion","fit_FpiMpi4_GL_w0_M1a", "no","M1a", "w0", notused );
 
    
    printf("\n\n///////////////////////////////////////Pion Mpi^2 and fpiMpi^4  GL  w0 M1a ///////////////////////\n");
    fit_info.Npar=10;
    fit_info.N=2;
    fit_info.function=fit_FpiMpi4_and_Mpi2_GL_am;
    
    init_Z( jack_files, head, jack_tot, &gjack, "w0","M1a");
    
    
    fit_out=fit_Mpi_fwMpi4_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info);
    
    //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion","fit_FpiMpi4_GL_w0_M1a_am", "no","M1a", "w0", notused );
    
    
    printf("\n\n///////////////////////////////////////Pion Mpi^2 and fpiMpi^4  linear  w0 M1a ///////////////////////\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=fit_FpiMpi4_and_Mpi2_linear;
        
    init_Z( jack_files, head, jack_tot, &gjack, "w0","M1a");

    double threshold_Mpiw=0.20 ;//0.164 ~19MeV 
    fit_out=fit_Mpi_fwMpi4_chiral_FVE_clover_threshold(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info, threshold_Mpiw);
 
    //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion","fit_FpiMpi4_GL_w0_M1a_threshold","no","M1a", "w0", notused);
 */
   
   
   /*
    printf("\n\n///////////////////////////////////////Pion of m_l GL   w0 M1a ///////////////////////\n");
    fit_info.Npar=6;
    fit_info.N=2;
    fit_info.function=fit_Fpi_and_Mpi_GL;
        
    init_Z( jack_files, head, jack_tot, &gjack, "w0","M1a");

        
    fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info);
 
    //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion","fit_Mpi_Fpi_GL_w0_M1a");
    
    threshold_Mpiw=0.24;
    fit_out=fit_Mpi_fw_chiral_FVE_clover_treshold(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info, threshold_Mpiw);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion","fit_Mpi_Fpi_GL_w0_M1a_260MeV");
    
    threshold_Mpiw=0.20;
    fit_out=fit_Mpi_fw_chiral_FVE_clover_treshold(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info, threshold_Mpiw);
    print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion","fit_Mpi_Fpi_GL_w0_M1a_190MeV");
    */
    
    std::vector<std::string>   M_Zp;
    M_Zp.emplace_back("M1a");
    M_Zp.emplace_back("M2a");
    M_Zp.emplace_back("M1b");
    M_Zp.emplace_back("M2b");
    std::vector<std::string>   GF_scale;
    GF_scale.emplace_back("w0");
    //GF_scale.emplace_back("t0_w0");
    //GF_scale.emplace_back("sqrtt0");
    fit_chi2_good.Nfits=0;
    
    std::vector<store_fit_clover>  mud;
    std::vector<store_fit_clover>  ms;
    std::vector<store_fit_clover>  ms_Mk;
    std::vector<store_fit_clover>  mc;
    std::vector<store_fit_clover>  mc_Ds;
    
    for (auto M : M_Zp){
        for (auto GF : GF_scale){
            init_Z( jack_files, head, jack_tot, &gjack, GF.c_str(), M.c_str());
            printf("\n\n %s    %s \n\n", M.c_str(),GF.c_str());
            char nameout[NAMESIZE];
       /*     
            printf("\n\n///////////////////////////////////////Pion Mpi^2 and fpiMpi^4  GL  w0 M1a ///////////////////////\n");
            fit_info.Npar=8;
            fit_info.N=2;
            fit_info.function=fit_FpiMpi4_and_Mpi2_GL;
            mysprintf(nameout,NAMESIZE,"fit_FpiMpi4_GL_%s_%s",GF.c_str(),M.c_str() );
            
        
            fit_out=fit_Mpi_fwMpi4_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info, argv,nameout );
            
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion",nameout, "no", M.c_str(), GF.c_str(), notused );
            
            printf("\n\n///////////////////////////////////////Pion Mpi^2 and fpiMpi^4  noGL    ///////////////////////\n");
            fit_info.Npar=7;
            fit_info.N=2;
            fit_info.function=fit_FpiMpi4_and_Mpi2_noGL;
            mysprintf(nameout,NAMESIZE,"fit_FpiMpi4_noGL_%s_%s",GF.c_str(),M.c_str() );
            
            
            fit_out=fit_Mpi_fwMpi4_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info, argv,nameout );
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion",nameout, "no", M.c_str(), GF.c_str(), notused );
            
            printf("\n\n///////////////////////////////////////Pion Mpi^2 and fpiMpi^4  noGL noA2    ///////////////////////\n");
            fit_info.Npar=6;
            fit_info.N=2;
            fit_info.function=fit_FpiMpi4_and_Mpi2_noGL_noA2;
            mysprintf(nameout,NAMESIZE,"fit_FpiMpi4_noGL_noA2_%s_%s",GF.c_str(),M.c_str() );
            
            
            fit_out=fit_Mpi_fwMpi4_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info, argv,nameout );
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion",nameout, "no", M.c_str(), GF.c_str(), notused );
            
            printf("\n\n///////////////////////////////////////Pion Mpi^2 and fpiMpi^4 from M  noGL noA2    ///////////////////////\n");
            fit_info.Npar=6;
            fit_info.N=2;
            fit_info.function=fit_FpiMpi4_fromM_and_Mpi2_noGL_noA2;
            mysprintf(nameout,NAMESIZE,"fit_FpiMpi4_fromM_noGL_noA2_%s_%s",GF.c_str(),M.c_str() );
            
            
            fit_out=fit_Mpi_fwMpi4_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info, argv,nameout );
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pionM",nameout, "no", M.c_str(), GF.c_str(), notused );
            
            
            
            printf("\n\n///////////////////////////////////////Pion Mpi^2 and fpiMpi^4 am  GL    ///////////////////////\n");
            fit_info.Npar=10;
            fit_info.N=2;
            fit_info.function=fit_FpiMpi4_and_Mpi2_GL_am;
            
            mysprintf(nameout,NAMESIZE,"fit_FpiMpi4_GL_am_%s_%s",GF.c_str(),M.c_str() );
            
            
            fit_out=fit_Mpi_fwMpi4_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout);
            
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion",nameout, "no", M.c_str(), GF.c_str(), notused );
            
            
            
            
            fit_info.Npar=6;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL;
            
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout);
            
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion",nameout , "no",M,GF, mud);
            
            
            
            double threshold_Mpiw=0.24;
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_%s_%s_260MeV",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover_treshold(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info, threshold_Mpiw);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion",nameout, "no",M,GF, mud);
            
            threshold_Mpiw=0.20;
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_%s_%s_190MeV",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover_treshold(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info, threshold_Mpiw);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result , gjack, head, "pion",nameout, "no",M,GF, mud);
            
            
            fit_info.Npar=8;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL_NL0_am;
            //double ***scale=double_malloc_3(4,3,jack_tot);
            
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_NLO_am_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "no",M,GF, mud);
            */
            
           
            fit_info.Npar=7;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL_NL0_am_fonly;
            //double ***scale=double_malloc_3(4,3,jack_tot);
            std::vector<int> myen={0,1,2,3,   8,4,5,6,   7,9};
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_NLO_am_fonly_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout,myen);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "yes",M,GF, mud);
            
            
                printf("\n\n///////////////////////////////////////MK     ///////////////////////\n");
                fit_info.Npar=4;
                fit_info.N=1;
                
                fit_info.function=fit_FK_and_MK_GL;
                char namefit[NAMESIZE];
                mysprintf(namefit,NAMESIZE,"MK_fK_GL_%s_%s",GF.c_str(),M.c_str() );
                fit=fit_MK_fK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk);
            
                
                printf("\n\n///////////////////////////////////////Mpi/MK     ///////////////////////\n");
                fit_info.Npar=4;
                fit_info.N=1;
                
                fit_info.function=fit_Mpi_MK;
                mysprintf(namefit,NAMESIZE,"Mpi_MK_%s_%s",GF.c_str(),M.c_str() );
                fit=fit_Mpi_MK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "pi/K",namefit, "no",   M,  GF , mud[mud.size()-1],ms_Mk);
                
                printf("\n\n///////////////////////////////////////MK/Mpi     ///////////////////////\n");
                fit_info.Npar=3;
                fit_info.N=1;
                
                fit_info.function=fit_MK_Mpi_FK_Fpi_GL_noP2;
                mysprintf(namefit,NAMESIZE,"MK_Mpi_%s_%s",GF.c_str(),M.c_str() );
                fit=fit_MK_Mpi_fK_fpi_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "no",   M,  GF , mud[mud.size()-1],ms_Mk);
                
                
                //////////////////////////   BC ////////////////////////// 
                myen={   8,4,5,6,   7};
                printf("\n\n///////////////////////////////////////Mpi/MK  (Mpi)   ///////////////////////\n");
                fit_info.Npar=4;
                fit_info.N=1;
                
                fit_info.function=fit_Mpi_MK;
                mysprintf(namefit,NAMESIZE,"Mpi_MK_BC_mud_ABC_%s_%s",GF.c_str(),M.c_str() );
                fit=fit_Mpi_MK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "pi/K",namefit, "no",   M,  GF , mud[mud.size()-1],ms_Mk);
                
                
                printf("\n\n///////////////////////////////////////MK/Mpi     ///////////////////////\n");
                fit_info.Npar=3;
                fit_info.N=1;
                
                fit_info.function=fit_MK_Mpi_FK_Fpi_GL_noP2;
                mysprintf(namefit,NAMESIZE,"MK_Mpi_BC_mud_ABC_%s_%s",GF.c_str(),M.c_str() );
                fit=fit_MK_Mpi_fK_fpi_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "no",   M,  GF , mud[mud.size()-1],ms_Mk);
                
                //////////////////////////   BC end //////////////////////////   
                myen={0,1,2,3,   8,4,5,6,   7,9};
                
                    printf("\n\n///////////////////////////////////////MD fD   ///////////////////////\n");
                    
                    fit_info.Npar=4;
                    fit_info.N=1;
                    fit_info.function=fit_MD_fD;
                    myen={0,1,2,3,   4,5,6,   7,9};
                    mysprintf(namefit,NAMESIZE,"MD_fD_%s_%s",GF.c_str(),M.c_str() );
                    
                    fit=fit_MD_fD_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                    printf("\n\nHERE HERE\n\n");
                    print_fit_D_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "D",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1], mc);
                    
                    printf("\n\n///////////////////////////////////////MDs fDs   ///////////////////////\n");
                    
                    fit_info.Npar=4;
                    fit_info.N=1;
                    //fit_info.Npar=8;
                    //fit_info.N=2;
                    
                    
                    fit_info.function=fit_MD_fD;
                    
                    mysprintf(namefit,NAMESIZE,"MDs_fDs_%s_%s",GF.c_str(),M.c_str() );
                    
                    fit=fit_MDs_fDs_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1], ms_Mk[ms_Mk.size()-1], myen);
                    print_fit_Ds_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "Ds",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1],mc[mc.size()-1], mc_Ds);
                    
                /*
            printf("\n\n///////////////////////////////////////Mpi fpi   ///////////////////////\n");
            
            fit_info.Npar=7;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL_NL0_am_fonly;
            //double ***scale=double_malloc_3(4,3,jack_tot);
            myen={0,1,2,   8,4,5,6,   7};
            
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_NLO_am_fonly_noA12_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout,myen);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "yes",M,GF, mud);
            
                printf("\n\n///////////////////////////////////////MK     ///////////////////////\n");
                fit_info.Npar=4;
                fit_info.N=1;
                
                fit_info.function=fit_FK_and_MK_GL;
                
                mysprintf(namefit,NAMESIZE,"MK_fK_GL_noA12_%s_%s",GF.c_str(),M.c_str() );
                fit=fit_MK_fK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk);
                
                    printf("\n\n///////////////////////////////////////MD fD   ///////////////////////\n");
                    
                    fit_info.Npar=4;
                    fit_info.N=1;
                    fit_info.function=fit_MD_fD;
                    myen={0,1,2,   4,5,6,   7};
                    mysprintf(namefit,NAMESIZE,"MD_fD_noA12_%s_%s",GF.c_str(),M.c_str() );
                    
                    fit=fit_MD_fD_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                    print_fit_D_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "D",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1], mc);
                    
                    printf("\n\n///////////////////////////////////////MDs fDs   ///////////////////////\n");
                    
                    fit_info.Npar=4;
                    fit_info.N=1;
                    
                    fit_info.function=fit_MD_fD;
                    
                    mysprintf(namefit,NAMESIZE,"MDs_fDs_noA12_%s_%s",GF.c_str(),M.c_str() );
                    
                    fit=fit_MDs_fDs_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1], ms_Mk[ms_Mk.size()-1], myen);
                    print_fit_Ds_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "Ds",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1],mc[mc.size()-1], mc_Ds);
                    
                    
               
            printf("\n\n///////////////////////////////////////Mpi fpi   ///////////////////////\n");
            
            
            fit_info.Npar=7;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL_NL0_am_fonly;
            //double ***scale=double_malloc_3(4,3,jack_tot);
            myen={0,1,   8,4,5,6,   7};
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_NLO_am_fonly_noA12A30_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout,myen);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "yes",M,GF, mud);
            
                printf("\n\n///////////////////////////////////////MK     ///////////////////////\n");
                fit_info.Npar=4;
                fit_info.N=1;
                
                fit_info.function=fit_FK_and_MK_GL;
                
                mysprintf(namefit,NAMESIZE,"MK_fK_GL_noA12A30_%s_%s",GF.c_str(),M.c_str() );
                fit=fit_MK_fK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk);
            
                    printf("\n\n///////////////////////////////////////MD fD   ///////////////////////\n");
                    
                    fit_info.Npar=4;
                    fit_info.N=1;
                    fit_info.function=fit_MD_fD;
                    myen={0,1,   4,5,6,   7};
                    
                    mysprintf(namefit,NAMESIZE,"MD_fD_noA12A30_%s_%s",GF.c_str(),M.c_str() );
                    
                    fit=fit_MD_fD_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                    print_fit_D_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "D",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1], mc);
                    
                    printf("\n\n///////////////////////////////////////MDs fDs   ///////////////////////\n");
                    
                    fit_info.Npar=4;
                    fit_info.N=1;
                    
                    fit_info.function=fit_MD_fD;
                    
                    mysprintf(namefit,NAMESIZE,"MDs_fDs_noA12A30_%s_%s",GF.c_str(),M.c_str() );
                    
                    fit=fit_MDs_fDs_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1], ms_Mk[ms_Mk.size()-1], myen);
                    print_fit_Ds_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "Ds",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1],mc[mc.size()-1], mc_Ds);
                   */
            /*
            
            fit_info.Npar=7;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL_NL0_am_fonly;
            //double ***scale=double_malloc_3(4,3,jack_tot);
            myen={0,1,2,3,   8,4,5,6 };
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_NLO_am_fonly_AB_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout,myen);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "yes",M,GF, mud);
            
                printf("\n\n///////////////////////////////////////MK     ///////////////////////\n");
                fit_info.Npar=4;
                fit_info.N=1;
                
                fit_info.function=fit_FK_and_MK_GL;
                
                mysprintf(namefit,NAMESIZE,"MK_fK_GL_AB_%s_%s",GF.c_str(),M.c_str() );
                fit=fit_MK_fK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk);
            
                    printf("\n\n///////////////////////////////////////MD fD   ///////////////////////\n");
                    
                    fit_info.Npar=4;
                    fit_info.N=1;
                    fit_info.function=fit_MD_fD;
                    
                    mysprintf(namefit,NAMESIZE,"MD_fD_AB_%s_%s",GF.c_str(),M.c_str() );
                    myen={0,1,2,3,   4,5,6 };
                    fit=fit_MD_fD_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                    print_fit_D_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "D",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1], mc);
                    
                   printf("\n\n///////////////////////////////////////MDs fDs   ///////////////////////\n");
                    
                    fit_info.Npar=4;
                    fit_info.N=1;
                    
                    fit_info.function=fit_MD_fD;
                    
                    mysprintf(namefit,NAMESIZE,"MDs_fDs_AB_%s_%s",GF.c_str(),M.c_str() );
                    
                    fit=fit_MDs_fDs_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1], ms_Mk[ms_Mk.size()-1], myen);
                    print_fit_Ds_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "Ds",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1],mc[mc.size()-1], mc_Ds);
                */
            /*
            printf("\n\n///////////////////////////////////////Mpi fpi   ///////////////////////\n");
            
            fit_info.Npar=6;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL;
            //double ***scale=double_malloc_3(4,3,jack_tot);
            myen={   8,4,5,6,   7};
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_BC_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout,myen);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "yes",M,GF, mud);
            
                printf("\n\n///////////////////////////////////////MK     ///////////////////////\n");
                fit_info.Npar=3;
                fit_info.N=1;
                
                fit_info.function=fit_FK_and_MK_GL_noP2;
                
                mysprintf(namefit,NAMESIZE,"MK_fK_GL_noP2_BC_%s_%s",GF.c_str(),M.c_str() );
                fit=fit_MK_fK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk);
                
                    printf("\n\n///////////////////////////////////////MD fD   ///////////////////////\n");
                    
                    fit_info.Npar=3;
                    fit_info.N=1;
                    fit_info.function=fit_MD_fD_noP2;
                    myen={   4,5,6,   7};
                    mysprintf(namefit,NAMESIZE,"MD_fD_noP2_BC_%s_%s",GF.c_str(),M.c_str() );
                    
                    fit=fit_MD_fD_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                    print_fit_D_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "D",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1], mc);
                
                    printf("\n\n///////////////////////////////////////MDs fDs   ///////////////////////\n");
                    
                    fit_info.Npar=3;
                    fit_info.N=1;
                    
                    fit_info.function=fit_MD_fD_noP2;
                    
                    mysprintf(namefit,NAMESIZE,"MDs_fDs_noP2_BC_%s_%s",GF.c_str(),M.c_str() );
                    
                    fit=fit_MDs_fDs_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1], ms_Mk[ms_Mk.size()-1], myen);
                    print_fit_Ds_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "Ds",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1],mc[mc.size()-1], mc_Ds);
                
            
            
            printf("\n\n///////////////////////////////////////Mpi fpi   ///////////////////////\n");
            
            fit_info.Npar=6;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL;
            //double ***scale=double_malloc_3(4,3,jack_tot);
            myen={3,   5,6,   7};
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_190MeV_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout,myen);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "yes",M,GF, mud);
            
                printf("\n\n///////////////////////////////////////MK     ///////////////////////\n");
                fit_info.Npar=3;
                fit_info.N=1;
                
                fit_info.function=fit_FK_and_MK_GL_noP2;
                
                mysprintf(namefit,NAMESIZE,"MK_fK_GL_noP2_190MeV_%s_%s",GF.c_str(),M.c_str() );
                fit=fit_MK_fK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk);
                
                    printf("\n\n///////////////////////////////////////MD fD   ///////////////////////\n");
                    
                    fit_info.Npar=3;
                    fit_info.N=1;
                    fit_info.function=fit_MD_fD_noP2;
                    
                    mysprintf(namefit,NAMESIZE,"MD_fD_noP2_190MeV_%s_%s",GF.c_str(),M.c_str() );
                    myen={3,   5,6,   7};
                    fit=fit_MD_fD_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                    print_fit_D_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "D",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1], mc);
                    
                    printf("\n\n///////////////////////////////////////MDs fDs   ///////////////////////\n");
                    
                    fit_info.Npar=3;
                    fit_info.N=1;
                    
                    fit_info.function=fit_MD_fD_noP2;
                    
                    mysprintf(namefit,NAMESIZE,"MDs_fDs_noP2_190MeV_%s_%s",GF.c_str(),M.c_str() );
                    
                    fit=fit_MDs_fDs_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1], ms_Mk[ms_Mk.size()-1], myen);
                    print_fit_Ds_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "Ds",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1],mc[mc.size()-1], mc_Ds);
         */ /*         
            printf("\n\n///////////////////////////////////////Mpi fpi   ///////////////////////\n");
                        
            fit_info.Npar=6;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL_NL0_am_fonly;
            //double ***scale=double_malloc_3(4,3,jack_tot);
            myen={2,3,  8,4 ,5,6,   7};
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_NL0_am_fonly_260MeV_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout,myen);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "yes",M,GF, mud);
            
                printf("\n\n///////////////////////////////////////MK     ///////////////////////\n");
                fit_info.Npar=3;
                fit_info.N=1;
                
                fit_info.function=fit_FK_and_MK_GL;
                
                mysprintf(namefit,NAMESIZE,"MK_fK_GL_260MeV_%s_%s",GF.c_str(),M.c_str() );
                fit=fit_MK_fK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk);
                        
                    printf("\n\n///////////////////////////////////////MD fD   ///////////////////////\n");
                    
                    fit_info.Npar=3;
                    fit_info.N=1;
                    fit_info.function=fit_MD_fD;
                    
                    mysprintf(namefit,NAMESIZE,"MD_fD_noP2_260MeV_%s_%s",GF.c_str(),M.c_str() );
                    myen={2,3,  4 ,5,6,   7};
                    fit=fit_MD_fD_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
                    print_fit_D_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "D",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1], mc);
                    
                    printf("\n\n///////////////////////////////////////MDs fDs   ///////////////////////\n");
                    
                    fit_info.Npar=3;
                    fit_info.N=1;
                    
                    fit_info.function=fit_MD_fD;
                    
                    mysprintf(namefit,NAMESIZE,"MDs_fDs_noP2_260MeV_%s_%s",GF.c_str(),M.c_str() );
                    
                    fit=fit_MDs_fDs_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1], ms_Mk[ms_Mk.size()-1], myen);
                    print_fit_Ds_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "Ds",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1],mc[mc.size()-1], mc_Ds);
            */    
            
       /*     
            fit_info.Npar=10;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL_NL0_am_m2;
            
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_NLO_am_m2_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout);
            
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "no",M,GF, mud);
            
           */ 
            /*
            ////////////////K
            printf("\n\n///////////////////////////////////////MKoverMpi fKoverfpi   ///////////////////////\n");
            
            fit_info.Npar=7;
            fit_info.N=2;
            fit_info.function=fit_MK_Mpi_FK_Fpi_GL;
            
            char nametex[NAMESIZE];
            //mysprintf(namefit,NAMESIZE,"MK_Mpi_fK_fpi_GL_M1a");
            
            mysprintf(namefit,NAMESIZE,"MK_Mpi_fK_fpi_GL_%s_%s",GF.c_str(),M.c_str() );
              
            Ci=(double**) malloc(sizeof(double*)*2);
            
            fit=fit_MK_Mpi_fK_fpi_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1] );
            //fit_MK_double_chiral_FVE_P40(jack_files, head , jack_tot, mass_index, gjack,  &result );
            //print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point, result, gjack, head, "pion",nameout, "no",M,GF, mud);
            print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "no",   M,  GF , mud[mud.size()-1],ms);
                
            */
            
        }       
    }
    
    printf("\n\n computig systematics pion \n\n");
    compute_systematic(argv, mud);
    printf("\n\n computig systematics K\n\n");
    compute_systematic_K(argv, ms_Mk);
    printf("\n\n computig systematics D\n\n");
    compute_systematic_D(argv, mc);
    printf("\n\n computig systematics Ds\n\n");
    compute_systematic_Ds(argv, mc_Ds);
    
    std::vector<store_fit_clover>().swap(mud);
    std::vector<store_fit_clover>().swap(ms);
    std::vector<store_fit_clover>().swap(ms_Mk);
    std::vector<store_fit_clover>().swap(mc);
    std::vector<store_fit_clover>().swap(mc_Ds);
    mud.resize(0);
    ms.resize(0);
    ms_Mk.resize(0);
    mc.resize(0);
    mc_Ds.resize(0);
    GF_scale[0]="t0_w0";
    for (auto M : M_Zp){
        for (auto GF : GF_scale){
            init_Z( jack_files, head, jack_tot, &gjack, GF.c_str(), M.c_str());
            printf("\n\n %s    %s \n\n", M.c_str(),GF.c_str());
            char nameout[NAMESIZE];
            char namefit[NAMESIZE];
            std::vector<int> myen={0,1,2,3,   8,4,5,6,   7,9};
            
            
            fit_info.Npar=7;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL_NL0_am_fonly;
            //double ***scale=double_malloc_3(4,3,jack_tot);
            myen={0,1,2,3,   8,4,5,6,   7,9};
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_NLO_am_fonly_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout,myen);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "yes",M,GF, mud);
            
            
            printf("\n\n///////////////////////////////////////MK     ///////////////////////\n");
            fit_info.Npar=4;
            fit_info.N=1;
            
            fit_info.function=fit_FK_and_MK_GL;
            mysprintf(namefit,NAMESIZE,"MK_fK_GL_%s_%s",GF.c_str(),M.c_str() );
            fit=fit_MK_fK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
            print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk);
            
            
            
            printf("\n\n///////////////////////////////////////MD fD   ///////////////////////\n");
            
            fit_info.Npar=4;
            fit_info.N=1;
            fit_info.function=fit_MD_fD;
            myen={0,1,2,3,   4,5,6,   7,9};
            mysprintf(namefit,NAMESIZE,"MD_fD_%s_%s",GF.c_str(),M.c_str() );
            
            fit=fit_MD_fD_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
            printf("\n\nHERE HERE\n\n");
            print_fit_D_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "D",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1], mc);
            
            printf("\n\n///////////////////////////////////////MDs fDs   ///////////////////////\n");
            
            fit_info.Npar=4;
            fit_info.N=1;
            //fit_info.Npar=8;
            //fit_info.N=2;
            
            fit_info.function=fit_MD_fD;
            mysprintf(namefit,NAMESIZE,"MDs_fDs_%s_%s",GF.c_str(),M.c_str() );
            
            fit=fit_MDs_fDs_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1], ms_Mk[ms_Mk.size()-1], myen);
            print_fit_Ds_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "Ds",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1],mc[mc.size()-1], mc_Ds);
            
        }
        
    }
    
    printf("\n\n computig systematics pion \n\n");
    compute_systematic(argv, mud, GF_scale[0].c_str());
    printf("\n\n computig systematics K\n\n");
    compute_systematic_K(argv, ms_Mk, GF_scale[0].c_str());
    printf("\n\n computig systematics D\n\n");
    compute_systematic_D(argv, mc, GF_scale[0].c_str());
    printf("\n\n computig systematics Ds\n\n");
    compute_systematic_Ds(argv, mc_Ds, GF_scale[0].c_str());
    
    std::vector<store_fit_clover>().swap(mud);
    std::vector<store_fit_clover>().swap(ms);
    std::vector<store_fit_clover>().swap(ms_Mk);
    std::vector<store_fit_clover>().swap(mc);
    std::vector<store_fit_clover>().swap(mc_Ds);
    mud.resize(0);
    ms.resize(0);
    ms_Mk.resize(0);
    mc.resize(0);
    mc_Ds.resize(0);
    GF_scale[0]="sqrtt0";
    for (auto M : M_Zp){
        for (auto GF : GF_scale){
            init_Z( jack_files, head, jack_tot, &gjack, GF.c_str(), M.c_str());
            printf("\n\n %s    %s \n\n", M.c_str(),GF.c_str());
            char nameout[NAMESIZE];
            char namefit[NAMESIZE];
            std::vector<int> myen={0,1,2,3,   8,4,5,6,   7,9};
            
            
            fit_info.Npar=7;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL_NL0_am_fonly;
            //double ***scale=double_malloc_3(4,3,jack_tot);
            myen={0,1,2,3,   8,4,5,6,   7,9};
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_NLO_am_fonly_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout,myen);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "yes",M,GF, mud);
            
            
            printf("\n\n///////////////////////////////////////MK     ///////////////////////\n");
            fit_info.Npar=4;
            fit_info.N=1;
            
            fit_info.function=fit_FK_and_MK_GL;
            mysprintf(namefit,NAMESIZE,"MK_fK_GL_%s_%s",GF.c_str(),M.c_str() );
            fit=fit_MK_fK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
            print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk);
            
            
            printf("\n\n///////////////////////////////////////MD fD   ///////////////////////\n");
            
            fit_info.Npar=4;
            fit_info.N=1;
            fit_info.function=fit_MD_fD;
            myen={0,1,2,3,   4,5,6,   7,9};
            mysprintf(namefit,NAMESIZE,"MD_fD_%s_%s",GF.c_str(),M.c_str() );
            
            fit=fit_MD_fD_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
            printf("\n\nHERE HERE\n\n");
            print_fit_D_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "D",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1], mc);
            
            printf("\n\n///////////////////////////////////////MDs fDs   ///////////////////////\n");
            
            fit_info.Npar=4;
            fit_info.N=1;
            //fit_info.Npar=8;
            //fit_info.N=2;
            
            fit_info.function=fit_MD_fD;
            mysprintf(namefit,NAMESIZE,"MDs_fDs_%s_%s",GF.c_str(),M.c_str() );
            
            fit=fit_MDs_fDs_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1], ms_Mk[ms_Mk.size()-1], myen);
            print_fit_Ds_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "Ds",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1],mc[mc.size()-1], mc_Ds);
            
        }
        
    }
    
    printf("\n\n computig systematics pion \n\n");
    compute_systematic(argv, mud, GF_scale[0].c_str() );
    printf("\n\n computig systematics K\n\n");
    compute_systematic_K(argv, ms_Mk, GF_scale[0].c_str());
    printf("\n\n computig systematics D\n\n");
    compute_systematic_D(argv, mc, GF_scale[0].c_str());
    printf("\n\n computig systematics Ds\n\n");
    compute_systematic_Ds(argv, mc_Ds, GF_scale[0].c_str());
    
    
    std::vector<store_fit_clover>().swap(mud);
    std::vector<store_fit_clover>().swap(ms);
    std::vector<store_fit_clover>().swap(ms_Mk);
    std::vector<store_fit_clover>().swap(mc);
    std::vector<store_fit_clover>().swap(mc_Ds);
    mud.resize(0);
    ms.resize(0);
    ms_Mk.resize(0);
    mc.resize(0);
    mc_Ds.resize(0);
    GF_scale[0]="w0";
    for (auto M : M_Zp){
        for (auto GF : GF_scale){
            init_Z( jack_files, head, jack_tot, &gjack, GF.c_str(), M.c_str());
            printf("\n\n %s    %s \n\n", M.c_str(),GF.c_str());
            char nameout[NAMESIZE];
            char namefit[NAMESIZE];
            std::vector<int> myen={0,1,2,3,   8,4,5,6,   7,9};
            
            printf("\n\n///////////////////////////////////////Mpi fpi   ///////////////////////\n");
            
            fit_info.Npar=6;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL;
            //double ***scale=double_malloc_3(4,3,jack_tot);
            myen={   8,4,5,6,   7,9};
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_BC_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout,myen);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "yes",M,GF, mud);
            
            printf("\n\n///////////////////////////////////////MK     ///////////////////////\n");
            fit_info.Npar=3;
            fit_info.N=1;
            
            fit_info.function=fit_FK_and_MK_GL_noP2;
            
            mysprintf(namefit,NAMESIZE,"MK_fK_GL_noP2_BC_%s_%s",GF.c_str(),M.c_str() );
            fit=fit_MK_fK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
            print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk);
            
            printf("\n\n///////////////////////////////////////Mpi/MK  (Mpi)   ///////////////////////\n");
            fit_info.Npar=4;
            fit_info.N=1;
            
            fit_info.function=fit_Mpi_MK;
            mysprintf(namefit,NAMESIZE,"Mpi_MK_BC_%s_%s",GF.c_str(),M.c_str() );
            fit=fit_Mpi_MK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
            print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "pi/K",namefit, "no",   M,  GF , mud[mud.size()-1],ms_Mk);
            
            
            printf("\n\n///////////////////////////////////////MK/Mpi     ///////////////////////\n");
            fit_info.Npar=3;
            fit_info.N=1;
            
            fit_info.function=fit_MK_Mpi_FK_Fpi_GL_noP2;
            mysprintf(namefit,NAMESIZE,"MK_Mpi_BC_%s_%s",GF.c_str(),M.c_str() );
            fit=fit_MK_Mpi_fK_fpi_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
            print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "no",   M,  GF , mud[mud.size()-1],ms_Mk);
            
            
            printf("\n\n///////////////////////////////////////MD fD   ///////////////////////\n");
            
            fit_info.Npar=3;
            fit_info.N=1;
            fit_info.function=fit_MD_fD_noP2;
            myen={   4,5,6,   7,9};
            mysprintf(namefit,NAMESIZE,"MD_fD_noP2_BC_%s_%s",GF.c_str(),M.c_str() );
            
            fit=fit_MD_fD_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
            print_fit_D_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "D",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1], mc);
            
            printf("\n\n///////////////////////////////////////MDs fDs   ///////////////////////\n");
            
            fit_info.Npar=3;
            fit_info.N=1;
            
            fit_info.function=fit_MD_fD_noP2;
            
            mysprintf(namefit,NAMESIZE,"MDs_fDs_noP2_BC_%s_%s",GF.c_str(),M.c_str() );
            
            fit=fit_MDs_fDs_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1], ms_Mk[ms_Mk.size()-1], myen);
            print_fit_Ds_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "Ds",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1],mc[mc.size()-1], mc_Ds);
            
            
        }
        
    }
    
    printf("\n\n computig systematics pion \n\n");
    compute_systematic(argv, mud, "BC" );
    printf("\n\n computig systematics K\n\n");
    compute_systematic_K(argv, ms_Mk, "BC");
    printf("\n\n computig systematics D\n\n");
    compute_systematic_D(argv, mc, "BC");
    printf("\n\n computig systematics Ds\n\n");
    compute_systematic_Ds(argv, mc_Ds, "BC");
    
    
    
    std::vector<store_fit_clover>().swap(mud);
    std::vector<store_fit_clover>().swap(ms);
    std::vector<store_fit_clover>().swap(ms_Mk);
    std::vector<store_fit_clover>().swap(mc);
    std::vector<store_fit_clover>().swap(mc_Ds);
    mud.resize(0);
    ms.resize(0);
    ms_Mk.resize(0);
    mc.resize(0);
    mc_Ds.resize(0);
    
    
    for (auto M : M_Zp){
        for (auto GF : GF_scale){
            init_Z( jack_files, head, jack_tot, &gjack, GF.c_str(), M.c_str());
            printf("\n\n %s    %s \n\n", M.c_str(),GF.c_str());
            char nameout[NAMESIZE];
            char namefit[NAMESIZE];
            std::vector<int> myen={0,1,2,3,   8,4,5,6,   7,9};
            
            printf("\n\n///////////////////////////////////////Mpi fpi   ///////////////////////\n");
            
            fit_info.Npar=6;
            fit_info.N=2;
            fit_info.function=fit_Fpi_and_Mpi_GL;
            //double ***scale=double_malloc_3(4,3,jack_tot);
            myen={3,   5,6,   7,9};
            mysprintf(nameout,NAMESIZE,"fit_Mpi_Fpi_GL_190MeV_%s_%s",GF.c_str(),M.c_str() );
            printf("\n%s\n",nameout);
            fit_out=fit_Mpi_fw_chiral_FVE_clover(jack_files,  head ,jack_tot, mass_index,gjack ,fit_info,argv,nameout,myen);
            //fit_chi2_good=save_fit(fit_chi2_good,fit_info,fit_out);
            print_fit_info( argv,jack_tot,  fit_out,  fit_info, phys_point,result, gjack, head, "pion",nameout, "yes",M,GF, mud);
            
            printf("\n\n///////////////////////////////////////MK     ///////////////////////\n");
            fit_info.Npar=3;
            fit_info.N=1;
            
            fit_info.function=fit_FK_and_MK_GL_noP2;
            
            mysprintf(namefit,NAMESIZE,"MK_fK_GL_noP2_190MeV_%s_%s",GF.c_str(),M.c_str() );
            fit=fit_MK_fK_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
            print_fit_K_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "K",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk);
            
            printf("\n\n///////////////////////////////////////MD fD   ///////////////////////\n");
            
            fit_info.Npar=3;
            fit_info.N=1;
            fit_info.function=fit_MD_fD_noP2;
            
            mysprintf(namefit,NAMESIZE,"MD_fD_noP2_190MeV_%s_%s",GF.c_str(),M.c_str() );
            myen={3,   5,6,   7,9};
            fit=fit_MD_fD_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1],myen );
            print_fit_D_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "D",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1], mc);
            
            printf("\n\n///////////////////////////////////////MDs fDs   ///////////////////////\n");
            
            fit_info.Npar=3;
            fit_info.N=1;
            
            fit_info.function=fit_MD_fD_noP2;
            
            mysprintf(namefit,NAMESIZE,"MDs_fDs_noP2_190MeV_%s_%s",GF.c_str(),M.c_str() );
            
            fit=fit_MDs_fDs_chiral_FVE_clover(jack_files,    head , jack_tot, mass_index,  gjack ,  fit_info ,   &result, namefit ,argv,mud[mud.size()-1], ms_Mk[ms_Mk.size()-1], myen);
            print_fit_Ds_info(argv,jack_tot, fit    ,  fit_info, phys_point, result, gjack, head , "Ds",namefit, "yes",   M,  GF , mud[mud.size()-1],ms_Mk[ms_Mk.size()-1],mc[mc.size()-1], mc_Ds);
        }
    }
    printf("\n\n computig systematics pion \n\n");
    compute_systematic(argv, mud,"190MeV");
    printf("\n\n computig systematics K\n\n");
    compute_systematic_K(argv, ms_Mk,"190MeV");
    printf("\n\n computig systematics D\n\n");
    compute_systematic_D(argv, mc,"190MeV");
    printf("\n\n computig systematics Ds\n\n");
    compute_systematic_Ds(argv, mc_Ds,"190MeV");
    
    
    
    
    return 0;
   
    
    
}
