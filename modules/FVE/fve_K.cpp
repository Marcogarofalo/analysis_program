#define FVE_K_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex>
#include "fve.hpp"
#include "bessel.hpp"
#include "global.hpp" 


double   dl1r = 3.80000e-04  ,   dl2r = 1.59000e-03  ,   dl3r =-2.91000e-03  ,   dl4r = 0.00000e+00;
double   dl5r = 1.46000e-03  ,   dl6r = 0.00000e+00  ,   dl7r =-4.90000e-04  ,   dl8r = 1.00000e-03;
double   mu=0.770;//M_pho

void volume(double ampi2,double amK2,double ameta2,double xl,double *g1tf,double *g2tm,double *g2tf,double dn,double dl1m,double dl2m,double dl1f,double dl2f,double dmu0);
void DJKMPQ(double ampi2,double amK2,double ameta2,double **DJPQ,double **DKPQ,double **DMPQ);
void DJKMloop(double amp2,double amq2,double z,double *DJloop,double *DKloop,double *DMloop);


void  FVE_K(double Bw,double fw,double frac_Lw, double mlw, double msw ,double  dmpi2, double dfpi,double  dmK2, double dfK,double *KM, double *Kf){
     double pigreco=3.141592653589790;
     double dn=(4.0*pigreco)*(4.0*pigreco);
     double sizeL=frac_Lw;
     
     double d2v=1.0;
     double dl1m=4.0*dl1r+dl3r-4.0*dl4r-dl5r+4.0*dl6r+2.0*dl8r;
     double dl2m=4.0*dl2r+dl3r;
     double dl1f=4.0*dl1r+dl3r-2.0*dl4r;
     double dl2f=dl5r;
     double dmu0=v_w0GeV*mu;
    
     
     double     ampi2=2.0*Bw*mlw;
     double     amK2=Bw*(msw+mlw);
     double     ameta2=2.0*Bw*(2.0*msw+mlw)/3.0;
     double     csipi=2.0*ampi2/(dn*fw*fw);
     double     csiK=amK2*csipi/ampi2;
     double     xl=sizeL*sqrt(ampi2);
     
     double g1tf, g2tm,g2tf;
     volume(ampi2,amK2,ameta2,xl,&g1tf,&g2tm,&g2tf,dn,dl1m,dl2m,dl1f,dl2f,dmu0);
         
     double     deltam=-d2v*sqrt(ampi2/amK2)*csiK*csiK*g2tm;
     double     deltaf=csipi*dfpi*(g1tf+d2v*csiK*g2tf)/dfK;
     deltaf=1.0/(1.0+deltaf)-1.0;
     //printf("d2v=%f   ampi2=%f   amk2=%f    csiK=%f  g2tm=%f \n",d2v,ampi2,amK2,csiK,g2tm);
     *KM=(1.0+deltam); 
     *Kf=1.0/(1.0+deltaf);
     
}


void volume(double ampi2,double amK2,double ameta2,double xl,double *g1tf,double *g2tm,double *g2tf,double dn,double dl1m,double dl2m,double dl1f,double dl2f,double dmu0){
    int nterm=20,j;
    int mul[20]={6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24};
    int nn,ipq;
    double **DJPQ,**DKPQ, **DMPQ;
    DJPQ=(double**) malloc(sizeof(double*)*3);
    DKPQ=(double**) malloc(sizeof(double*)*3);
    DMPQ=(double**) malloc(sizeof(double*)*3);
    for(ipq=0;ipq<3;ipq++){
        DJPQ[ipq]=(double*) malloc(sizeof(double)*4);
        DKPQ[ipq]=(double*) malloc(sizeof(double)*4);
        DMPQ[ipq]=(double*) malloc(sizeof(double)*4);
    }
    
    double GPQ[3][7],GPQ1[3][7],GPQ2[3][7],SPQ[3][3][7];
    double xpiK=ampi2/amK2;
    double xpieta=ampi2/ameta2;
    double xetaK=ameta2/amK2;
    double dlpi=log(ampi2/(dmu0*dmu0));
    double dlK=log(amK2/(dmu0*dmu0));
    double dleta=log(ameta2/(dmu0*dmu0));
    
    DJKMPQ(ampi2,amK2,ameta2,DJPQ,DKPQ,DMPQ);
    for(ipq=1;ipq<=2;ipq++){
         GPQ[ipq][1]=DJPQ[ipq][0];
         GPQ1[ipq][1]=DJPQ[ipq][1];
         GPQ2[ipq][1]=DJPQ[ipq][2];
         GPQ[ipq][2]=amK2*DJPQ[ipq][1];
         GPQ1[ipq][2]=amK2*DJPQ[ipq][2];
         GPQ2[ipq][2]=amK2*DJPQ[ipq][3];
         GPQ[ipq][3]=DKPQ[ipq][0];
         GPQ1[ipq][3]=DKPQ[ipq][1];
         GPQ2[ipq][3]=DKPQ[ipq][2];
         GPQ[ipq][4]=amK2*DKPQ[ipq][1];
         GPQ1[ipq][4]=amK2*DKPQ[ipq][2];
         GPQ2[ipq][4]=amK2*DKPQ[ipq][3];
         GPQ[ipq][5]=DMPQ[ipq][0];
         GPQ1[ipq][5]=DMPQ[ipq][1];
         GPQ2[ipq][5]=DMPQ[ipq][2];
         GPQ[ipq][6]=amK2*DMPQ[ipq][1];
         GPQ1[ipq][6]=amK2*DMPQ[ipq][2];
         GPQ2[ipq][6]=amK2*DMPQ[ipq][3];
    }

    *g1tf=0.0;
    *g2tm=0.0;
    double s4m=0;
    double s4tm=0.0;
    *g2tf=0.0;
    double s4f=0;
    double s4tf=0.0;
    double y,dmn,b0,b1,b2;
    double cpi1,cpi2,cK1,ceta1,cK2,ceta2;
    
    for(nn=0;nn<nterm;nn++){
         y=xl*sqrt(double(nn+1));
         dmn=double(mul[nn]);
         b0=dbesk0(y);
         b1=dbesk1(y);      
         b2=(b0+2.0*b1/y)/y;
         *g1tf=(*g1tf)+dmn*b1/y;
         cpi1=xpiK*xpiK/(4.0*(1.0-xpiK));
         cK1=((7.0+xpiK)/2.0+(1.0-10.0*xpiK+xpiK*xpiK)/(6.0* (xetaK-1.0))-4.0/(1.0-xpiK))/16.0;
         ceta1=(2.0/3.0+(1.0-xpiK)*(xetaK-1.0)+53.0*xpiK/9.0-xpiK*xpiK/3.0-(1.0-10.0*xpiK+xpiK*xpiK)/(3.0* (xetaK-1.0)))/32.0;
         cpi2=-5.0*xpiK*xpiK/(2.0*(1.0-xpiK));
         cK2=xpiK*(5.0/(1.0-xpiK)-1.0/(xetaK-1.0))/2.0;
         ceta2=xpiK*xetaK/(2.0*(xetaK-1.0));
         *g2tm=(*g2tm)+dmn*(b1*(xpiK/9.0+8.0*dn*xpiK*dl1m+cpi1*dlpi+  cK1*dlK+ceta1*dleta)+b2*(-8.0*dn*xpiK*dl2m+cpi2*dlpi+    cK2*dlK+ceta2*dleta))/y;
         for(ipq=1;ipq<=2;ipq++){
            for(j=1;j<=6;j++){
                SPQ[ipq][0][j]=sqrt(xpiK)*(GPQ[ipq][j]*b1-2.0*ampi2*amK2* GPQ2[ipq][j]*b2);
                SPQ[ipq][1][j]=2.0*xpiK*sqrt(ampi2*amK2)*GPQ1[ipq][j]*b2;
                SPQ[ipq][2][j]=xpiK*sqrt(xpiK)*GPQ[ipq][j]*b2;
            }
         }
         s4m=3.0*(1.0+xpiK)*(1.0+xpiK)  *SPQ[1][0][1]/32.0-5.0*(1.0+xpiK)* SPQ[1][1][1]/8.0-19.0*SPQ[1][2][1]/8.0-3.0*(1.0-xpiK*xpiK)*SPQ[1][0][3]/16.0+13.0*(1.0-xpiK)* SPQ[1][1][3]/8.0-3.0*(xpiK*SPQ[1][0][5]+SPQ[1][2][5])/2.0+      (1.0+xpiK)*(1.0+xpiK)*SPQ[2][0][1]/96.0-(1.0+xpiK)*SPQ[2][1][1]/    8.0-3.0*SPQ[2][2][1]/8.0+(1.0+xpiK)*(3.0*xetaK+       2.0*xpiK-5.0)*SPQ[2][0][3]/16.0+3.0*(1.0-2.0*xpiK+    xetaK)*SPQ[2][1][3]-3.0*(xpiK*SPQ[2][0][5]+SPQ[2][2][5])/2.0;
         
         s4tm=(s4tm)+dmn*s4m/y;
         
         *g2tf=(*g2tf)+dmn*(b1*(12.0*dn*xpiK*dl1f-3.0*dn*dl2f*    (1.0+xpiK)+3.0*xpiK*((dlK-xpiK*dlpi)/(1.0-xpiK)+     (xetaK*dleta-dlK)/(xetaK-1.0)+2.0*dlpi*(xpieta-     9.0/4.0))/16.0+3.0*(2.0*dlK+3.0*xetaK*dleta)/       32.0)+b2*xpiK*(-24.0*dn*dl2m+15.0*(dlK-xpiK*dlpi)/       (2.0*(1.0-xpiK))+3.0*(xetaK*dleta-dlK)/(2.0*     (xetaK-1.0))))/y;
         
         s4f=-15.0*(1.0+xpiK)*SPQ[1][1][1]/16.0-57.0*SPQ[1][2][1]/     8.0-9.0*(1.0+xpiK)*(1.0+xpiK)*SPQ[1][0][2]/64.0+15.0*(1.0+     xpiK)*SPQ[1][1][2]/16.0+57.0*SPQ[1][2][2]/16.0+9.0*(1.0-    5.0*xpiK)*SPQ[1][0][3]/32.0+(6.0-15.0*xpiK/8.0)*      SPQ[1][1][3]+15.0*SPQ[1][2][3]/4.0+9.0*(1.0-xpiK*xpiK)*    SPQ[1][0][4]/32.0-39.0*(1.0-xpiK)*SPQ[1][1][4]/16.0-9.0*     (xpiK*SPQ[1][0][5]+2.0*SPQ[1][2][5]-xpiK*SPQ[1][0][6]-      SPQ[1][2][6])/4.0-3.0*(1.0+xpiK)*SPQ[2][1][1]/16.0-9.0*   SPQ[2][2][1]/8.0-(1.0+xpiK)*(1.0+xpiK)*SPQ[2][0][2]/64.0+3.0*      (1.0+xpiK)*SPQ[2][1][2]/16.0+9.0*SPQ[2][2][2]/16.0+       (27.0*(xetaK-1.0)/32.0-3.0*(1.0+xpiK)/16.0)*    SPQ[2][0][3]+3.0*(4.0-2.0*xpiK+3.0*xetaK)*SPQ[2][1][3]/      8.0+9.0*SPQ[2][2][3]/4.0+3.0*(5.0-3.0*xetaK)*(1.0-    xpiK*xpiK)*SPQ[2][0][4]/32.0-9.0*(1.0-2.0*xpiK+xetaK)*   SPQ[2][1][4]/16.0-9.0*(xpiK*SPQ[2][0][5]+2.0*SPQ[2][2][5]-  xpiK*SPQ[2][0][6]-SPQ[2][2][6])/4.0;
         s4tf=(s4tf)+dmn*s4f/y;
     }
      *g1tf=-3.0*(*g1tf)/2.0;
      *g2tm=3.0*(sqrt(xpiK)*(*g2tm)+(s4tm));
      *g2tf=2.0*((*g2tf)+sqrt(amK2/ampi2)*(s4tf));
      
      
      
      for(ipq=0;ipq<3;ipq++){
        free(DJPQ[ipq]);
        free(DKPQ[ipq]);
        free(DMPQ[ipq]);
      }
      free(DJPQ);free(DKPQ);free(DMPQ);
}





void DJKMPQ(double ampi2,double amK2,double ameta2,double **DJPQ,double **DKPQ,double **DMPQ){
      
      double DJloop[4],DKloop[4],DMloop[4];
      double z=amK2+ampi2;
      int kk;
      DJKMloop(amK2,ampi2,z,DJloop,DKloop,DMloop);
      for( kk=0;kk<=3;kk++){
         DJPQ[1][kk]=DJloop[kk];
         DKPQ[1][kk]=DKloop[kk];
         DMPQ[1][kk]=DMloop[kk];
      }
      DJKMloop(ameta2,amK2,z,DJloop,DKloop,DMloop);
      for( kk=0;kk<=3;kk++){
         DJPQ[2][kk]=DJloop[kk];
         DKPQ[2][kk]=DKloop[kk];
         DMPQ[2][kk]=DMloop[kk];
      }
}
 
void DJKMloop(double amp2,double amq2,double z,double *DJloop,double *DKloop,double *DMloop){

      
      using namespace std::complex_literals;
      std::complex<double> Rho,dJlogz;
      double  pigreco=3.14159265358979;
      double Delta=amp2-amq2;
      double Sigma=amp2+amq2;
      double dlqp=log(amq2/amp2);
      double A=Sigma+2.0*amp2*amq2*dlqp/Delta;
      double Rho0=(z+Delta)*(z+Delta)-4.0*amp2*z;
      if (Rho0<0){    
       Rho=sqrt(fabs(Rho0))*1i;
       dJlogz=pigreco*1i;
      }
      else{
       Rho=sqrt(Rho0);
       dJlogz=log((z+Delta-2.0*amp2+Rho0)/(z+Delta-2.0*amp2-Rho0));
      }
      DJloop[0]=real(0.50*(2.0+(Delta/z-Sigma/Delta)*dlqp-Rho*     dJlogz/z));
      DJloop[1]=real((-2.0*z-Delta*dlqp+(-2.0*amp2*z+Delta*(z+      Delta))*dJlogz/Rho)/(2.0*z*z));
      DJloop[2]=real((z*Rho*(-6.0*amp2*z+(z+Delta)*(z+2.0*Delta))+    Delta*(pow(Rho,3))*dlqp-(6.0*amp2*amp2*z*z+Delta*pow((z+   Delta),3)-2.0*amp2*z*(z*z+3.0*z*Delta+3.0*  pow(Delta,2)))*dJlogz)/(pow((z*Rho),3)));
      DJloop[3]=real((-z*Rho*(60.0*amp2*amp2*z*z+pow((z+Delta),2)*(2.0* z*z+9.0*z*Delta+6.0*pow(Delta,2))-2.0*amp2*z*(13.0*z* z+30.0*z*Delta+21.0*pow(Delta,2)))-3.0*Delta*(pow(Rho,5))* dlqp+3.0*(Delta*pow((z+Delta),5)-20.0*pow((amp2*z),3)+6.0*  amp2*amp2*z*z*(2.0*z*z+5.0*z*Delta+5.0*pow(Delta,2))-  2.0*amp2*z*(z+Delta)*(pow(z,3)+5.0*z*z*Delta+10.0*z*   pow(Delta,2)+5.0*pow(Delta,3)))*dJlogz)/(pow(z,4)*pow(Rho,5)));
      DKloop[0]=Delta*DJloop[0]/(2.0*z);
      DKloop[1]=Delta*(z*DJloop[1]-DJloop[0])/(2.0*z*z);
      DKloop[2]=Delta*(2.0*DJloop[0]+z*(z*DJloop[2]-2.0*DJloop[1]))/  (2.0*pow(z,3));
      DKloop[3]=Delta*(-6.0*DJloop[0]+z*(6.0*DJloop[1]+z*(z*    DJloop[3]-3.0*DJloop[2])))/(2.0*pow(z,4));
      DMloop[0]=(z/4.0-Sigma/2.0+pow(Delta,2)/z)*DJloop[0]/(3.0*z)+  1.0/18.0-A/(6.0*z);
      DMloop[1]=(2.0*A*z+2.0*(-4.0*pow(Delta,2)+z*Sigma)*DJloop[0]+z*  (z*z+4.0*pow(Delta,2)-2.0*z*Sigma)*DJloop[1])/  (12.0*pow(z,3));
      DMloop[2]=(4.0*(6.0*pow(Delta,2)-z*Sigma)*DJloop[0]+z*(-4.0*A+    4.0*(z*Sigma-4.0*pow(Delta,2))*DJloop[1]+z*(z*z+4.0*    pow(Delta,2)-2.0*z*Sigma)*DJloop[2]))/(12.0*pow(z,4));
      DMloop[3]=(12.0*(z*Sigma-8.0*pow(Delta,2))*DJloop[0]+6.0*z*   (2.0*A+2.0*(6.0*pow(Delta,2)-z*Sigma)*DJloop[1]+z*  (z*Sigma-4.0*pow(Delta,2))*DJloop[2])+z*z*z*(z*z+4.0* pow(Delta,2)-2.0*z*Sigma)*DJloop[3])/(12.0*pow(z,5));
}

