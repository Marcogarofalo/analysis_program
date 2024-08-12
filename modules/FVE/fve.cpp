#define FVE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "fve.hpp"
#include "bessel.hpp"
#include "global.hpp"


void volume(double xl,double *g1t,double *g2tm, double *g2tf, double *g2tv ,double g0,double g1,double g2,double g3,double coef1,double coef2,double coef3 ){
      *g1t=0.0;
      *g2tm=0.0;
      *g2tf=0.0;
      *g2tv=0.0;
      int mul[20]={6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24};
      int nn;
      double y,dmn,b0,b1,b2;
      
      for (nn=0;nn<20;nn++){
      
         y=xl*sqrt(double(nn+1));
         dmn=double(mul[nn]);
         b0=dbesk0(y);
         b1=dbesk1(y);
         b2=(b0+2.0*b1/y)/y;
         *g1t=*g1t+4.0*dmn*b1/y;
         *g2tm=*g2tm+dmn*(coef1*b1+coef2*b2+13.0*g0*b1/3.0-(40.0*g0+32.0*g1+26.0*g2)*b2/3.0)/y;
         *g2tf=*g2tf+2.0*dmn*(coef3*b1+coef2*b2+(8.0*g0-13.0*g1)*b1/6.0-(40.0*g0-12.0*g1-8.0*g2-13.0*g3)*b2/3.0)/y;
         *g2tv=*g2tv+dmn*b0;
      }
    
}

//dmpi2=a^2M_pi^2 ,  dfpi=af_pi,  Bw=b_0 w_0,  fw=f_0w_0, 
void  FVE(double w0 ,double w0a,double dl3b, double dl4b,double Bw,double fw,int Lsize, double mw ,double  dmpi2, double dfpi,double *KM, double *Kf){

    
     double r0=w0,r0a=w0a;
     double dl1b=-0.4,dl2b=4.3;

     double B0=Bw/w0;
     double fpi0=fw/w0;
     
      double pigreco=3.14159265358979;
      double dn=(4.0*pigreco)*(4.0*pigreco);
      double dmpiex=0.13498;
      double g0=2-0.5*pigreco;
      double g1=0.25*pigreco-0.5;
      double g2=0.5-0.125*pigreco;
      double g3=3*pigreco/16-0.5;
//ccccccc
//ccccccc                    approximation S^{(4)}=0
//ccccccc
      g0=0;
      g1=0;
      g2=0;
      g3=0;

      int nloopv=2,imPi=0;
      double d2v=0;
      d2v=1.0; 
      double ampiex2=(r0*dmpiex)*(r0*dmpiex);
      double coef10=-55.0/18.0+4.0*dl1b+8.0*dl2b/3.0-2.50*dl3b-2.0*dl4b;   
      double coef20=112.0/9.0-8.0*dl1b/3.0-32.0*dl2b/3.0;
      double coef30=-7.0/9.0+2.0*dl1b+4.0*dl2b/3.0-3.0*dl4b;
 
      double corr0=0.0;
     // if(imPi.lt.0) corr0=1.d0
     // do im=1,nmasses
       //  if(imPi.le.0) then       // in this case imPi=0 then it is this "if" that get activated
          double ampi2=2.0*B0*r0*mw;
          double csi=2.0*ampi2/(dn*fpi0*r0*fpi0*r0);
   /*      else
          ampi2=dmpi2(im)
          csi=2.d0*ampi2/(dn*dfpi(im)*dfpi(im))
         endif*/
         double xl=(double(Lsize))*sqrt(ampi2)/r0a;
         double dlpi=log(ampiex2/ampi2);
         double coef1=coef10+13.0*dlpi/6.0;
         double coef2=coef20-40.0*dlpi/3.0;
         double coef3=coef30+dlpi/3.0;
         
         double g1t,g2tm,g2tf,g2tv;
         volume(xl,&g1t,&g2tm,&g2tf,&g2tv,g0,g1,g2,g3,coef1,coef2,coef3);
         
         double corrg1=corr0*(-2.0*dl4b*g1t+dl3b*g2tv);
         double deltam=0.250*csi*g1t-d2v*csi*csi*(g2tm-0.250*corrg1);
         
         double deltaf=-csi*g1t+d2v*csi*csi*(g2tf-corrg1);
         double deltam2=1.0/(1.0+deltam)*(1.0/(1.0+deltam))  -1.0;
         deltaf=1.0/(1.0+deltaf)-1.0;
         double dmpi2c=dmpi2*(1.0+deltam2);
         double dfpic=dfpi*(1.0+deltaf);
         *KM=(1+deltam);
         *Kf=1.0/(1.0+deltaf);

}


void  FVE_Mpi(double afm,double dl1phys, double dl2phys,double dl3phys, double dl4phys,int Lsize,double  ampi, double afpi,int nloop, int  ischeme, double *KM, double *Kf){

    
     
      double pigreco=3.14159265358979;
      double dn=(4.0*pigreco)*(4.0*pigreco);
      double dmpi0=0.13498;
      double fpi0=0.1225;
      if (ischeme!=1) ischeme=0;
      //ischeme=0 csi = MPi^2/(4pi fPi)^2
      //ischeme=1 csi = MPi^2/(4pi f0)^2
      double am1g=0.19731/afm;//a in GeV
      double ampi0=dmpi0/am1g;
      double afpi0=fpi0/am1g;

      double ampi2=ampi*ampi;
      double csi;
      if(ischeme==0) 
          csi=2.0*ampi2/(dn*afpi*afpi);
      else
          csi=2.0*ampi2/(dn*afpi0*afpi0);
      double ampiL=((double) Lsize)*ampi;
      double dmpi=ampi*am1g;
      double dmpiex=dmpi0;
      double g0=2-0.5*pigreco;
      double g1=0.25*pigreco-0.5;
      double g2=0.5-0.125*pigreco;
      double g3=3*pigreco/16-0.5;

      double d2v=0;
      if (nloop>1)
        d2v=1.0; 
      double dlpi=2.0*log(dmpiex/dmpi);
      double dl1b=dl1phys+dlpi;
      double dl2b=dl2phys+dlpi;
      double dl3b=dl3phys+dlpi;
      double dl4b=dl4phys+dlpi;
      double coef1=-55.0/18.0+4.0*dl1b+8.0*dl2b/3.0-2.50*dl3b-2.0*dl4b;
      double coef2=112.0/9.0-8.0*dl1b/3.0-32.0*dl2b/3.0;
      double coef3=-7.0/9.0+2.0*dl1b+4.0*dl2b/3.0-3.0*dl4b;
        
      double g1t,g2tm,g2tf,g2tv;
      volume(ampiL,&g1t,&g2tm,&g2tf,&g2tv,g0,g1,g2,g3,coef1,coef2,coef3);
      
      double corrgm=0.0;
      double corrgf=0.0;
      if(ischeme==1) 
       corrgf=-2.0*dl4b*g1t;
      else if(ischeme==2) 
       corrgm=dl3b*g2tv;
      else if(ischeme>2) {
       corrgm=dl3b*g2tv;
       corrgf=-2.0*dl4b*g1t;
      }
      double corrg=corrgm+corrgf ;
      double deltam=0.250*csi*g1t-d2v*csi*csi*(g2tm-0.250*corrg);
      double deltaf=-csi*g1t+d2v*csi*csi*(g2tf-corrg);

      double fsemPi=1.0+deltam;
      double fsefPi=1.0+deltaf;
      *KM=fsemPi;
      *Kf=fsefPi;
      /*
      ampic=ampi(im)/fsemPi;
      afpic=afpi(im)/fsefPi;
      errmpic=errampi(im)/fsemPi;
      errfpic=errafpi(im)/fsefPi;
*/
      
}


double FVE_GL(double  Lsize_w,double  aml, double af, double aB){

        int j;
        double pi=3.14159265358979;
        double xi=2*aB*aml/(4*pi*4*pi*af*af);
        double lambda=sqrt(2*aB*aml)*((double) Lsize_w);
        double m[5];
        
        m[0]=0.;
        m[1]=6.;
        m[2]=12.;
        m[3]=8.;
        m[4]=6.;
        
        double Sum=0;
        for (j=1;j<5;j++)
             Sum=Sum+m[j]*dbesk1(lambda*sqrt(j))/sqrt(j);
        
        Sum=Sum*4./(lambda);
        double Delta=-2*Sum* xi ;
        
        return Delta;
        // Mpi^2(L)=Mpi^2(inf)[1-1/4 Delta]^2
        //Mpi^2(inf)=2Bm [1+xi log xi +C1 xi]  //C1 coefficiente
        
        //fpi(L)=fpi(inf)[1+Delta]
        //fpi(inf)=f[1-2xi log xi +2 A xi] // A coefficente
}

double FVE_GL_fast(double Lsize_w /* L/w0 */,double  aml, double af, double aB){

        int j;
        double pi=3.14159265358979;
        double xi=2*aB*aml/(4*pi*4*pi*af*af);
        double lambda=sqrt(2*aB*aml)*((double) Lsize_w);
        double m[5];

        //printf("L=%d   aml=%g    af=%g     aB=%g\n",Lsize,aml,af,aB);
        m[0]=0.;
        m[1]=6.;
        m[2]=12.;
        m[3]=8.;
        m[4]=6.;
    
       
        double Sum=0;
        for (j=1;j<5;j++)
             Sum=Sum+m[j]*exp(-lambda*sqrt(j))/pow(j,3./4.);
        //Sum=Sum*4*sqrt(pi/2.)/pow(lambda,3./2.);
        Sum=Sum*sqrt(pi/2.)*4./(pow(lambda,3./2.));
        //Sum=2*Sum* (M[i]*M[i]/(pow(4.*pi*f0w0/w0a,2)  ))  ;
        double Delta=-2*Sum* xi  ;
        //printf("DElta=%g\n",Delta);
        return Delta;
        // Mpi^2(L)=Mpi^2(inf)[1-1/4 Delta]^2
        //Mpi^2(inf)=2Bm [1+xi log xi +C1 xi]  //C1 coefficiente
        
        //fpi(L)=fpi(inf)[1+Delta]
        //fpi(inf)=f[1-2xi log xi +2 A xi] // A coefficente
}


double FVE_GL_Mpi(double L_a /* L/w0 */,double  xi,double af_PS){

        int j;
        double pi=3.14159265358979;
        
        double lambda=sqrt(xi)*4*pi*af_PS*((double) L_a);
        double m[5];

        //printf("L=%d   aml=%g    af=%g     aB=%g\n",Lsize,aml,af,aB);
        m[0]=0.;
        m[1]=6.;
        m[2]=12.;
        m[3]=8.;
        m[4]=6.;
    
       
        double Sum=0;
        for (j=1;j<5;j++)
             Sum=Sum+m[j]*exp(-lambda*sqrt(j))/pow(j,3./4.);
        //Sum=Sum*4*sqrt(pi/2.)/pow(lambda,3./2.);
        Sum=Sum*sqrt(pi/2.)*4./(pow(lambda,3./2.));
        //Sum=2*Sum* (M[i]*M[i]/(pow(4.*pi*f0w0/w0a,2)  ))  ;
        double Delta=-2*Sum* xi  ;
        //printf("DElta=%g\n",Delta);
        return Delta;
        
}
