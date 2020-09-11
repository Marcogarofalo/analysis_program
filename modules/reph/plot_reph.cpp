#define plot_reph_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"

#include "resampling.hpp"
#include "read.hpp"
//#include "m_eff.h"
#include "gnuplot.hpp"
//#include "eigensystem.h"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "rand.hpp"
#include "non_linear_fit.hpp"
#include "continuum_reph.hpp"
#include "fve.hpp"
#include "indices.hpp"
#include "global_reph.hpp"
#include "mutils.hpp"

#include <unistd.h>
#include <omp.h> 
 
double derivative_xi(int ind,int j ,struct fit_result fit_out, struct fit_type fit_info, double *x){
    
    double df;
    int i;
    int Nvar=fit_info.Nvar;
    int Npar=fit_info.Npar;
    double h=1e-5;
    int n=1;
    double *P=(double *) malloc(sizeof(double)*Npar) ;  
    double *Ph=(double*) malloc(sizeof(double)*Nvar);

    if (x[ind]==0){ 
        free(Ph);
        free(P);
        return 0;
    }
    
    for (i=0;i<Nvar;i++)
        Ph[i]=x[i];
    for (i=0;i<Npar;i++)
        P[i]=fit_out.P[i][j];
    
    Ph[ind]=x[ind]-2.*h;
    df=fit_info.function(n,Nvar,Ph,Npar,P);
    
    Ph[ind]=x[ind]-h;
    df-=8*fit_info.function(n,Nvar,Ph,Npar,P);
     
    Ph[ind]=x[ind]+h;
    df+=8*fit_info.function(n,Nvar,Ph,Npar,P);
    
    Ph[ind]=x[ind]+2.*h;
    df-=fit_info.function(n,Nvar,Ph,Npar,P);
    
    df/=(12.*h);
    
    df/=(2*x[ind]);
    free(Ph);
    free(P);
    
    return df;

} 
 
double **subtract_non_Lorentz_and_int(char **argv, const char  *string,int ikt,int iks,double **kp_tot,double xmin,double xmax,  struct header *head ,int Njack ,int npar_fun,double *fit_function(int,double),struct fit_result fit_out, struct fit_type fit_info  , struct reph_jack *gJ,int e,const char *AV){
   double ***y,*x,**r,*chi2,**tmp,**fit; 
   int i,j,ii;  
   int Ntot=head[e].nmoms*head[e].nmoms*head[e].nmoms;
   int Neff=0;
   int imom0,imoms,imomt;
   double *xp=(double*) malloc(sizeof(double)*fit_info.Nvar); 

   for(imom0=0;imom0<head[e].nmoms;imom0++){       
   for(imomt=0;imomt<head[e].nmoms;imomt++){
   for(imoms=0;imoms<head[e].nmoms;imoms++){
       i=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);
       if (kp_tot[i][Njack-1]>xmin &&  kp_tot[i][Njack-1]<xmax){
           Neff=Neff+1;
       }
   }}}
   
   error(Neff==0,1,"interpolate_FX","no data to interpolate ikt=%d  iks=%d  e=%d",ikt,iks,e);
       
   double **F=(double**) malloc(sizeof(double*)*Neff);
   double **kp=(double**) malloc(sizeof(double*)*Neff);
   ii=0;
   int     i_mDs=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);
   int     i_mD=index_ensemble_twopt_fit(head[e],ikt,0,0,0);
   int     i_mK=index_ensemble_twopt_fit(head[e],0,iks,0,0);
   int     i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
   double cp2,ck2;
   for(imom0=0;imom0<head[e].nmoms;imom0++){       
   for(imomt=0;imomt<head[e].nmoms;imomt++){
   for(imoms=0;imoms<head[e].nmoms;imoms++){
      i=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);
      if (kp_tot[i][Njack-1]>xmin &&  kp_tot[i][Njack-1]<xmax){
           F[ii]=(double*) malloc(sizeof(double)*Njack);
           kp[ii]=(double*) malloc(sizeof(double)*Njack);

           xp[0]=head[e].k[head[e].nk+ikt]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
           xp[1]=gJ[e].w0[Njack-1];//w0
           xp[2]=gJ[e].xG[i][Njack-1];//x_gamma
           xp[3]=gJ[e].M_PS[i_mp][Njack-1]*gJ[e].w0[Njack-1];//x_gamma
                 
           xp[4]=head[e].k[head[e].nk+iks]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ms*w0
           xp[5]=gJ[e].M_PS[i_mK][Njack-1]*gJ[e].w0[Njack-1];//M_K r0
                    
           xp[6]=head[e].k[head[e].nk+ikt]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//mc*w0
           xp[7]=gJ[e].M_PS[i_mD][Njack-1]*gJ[e].w0[Njack-1];//M_D r0
           xp[8]=gJ[e].M_PS[i_mDs][Njack-1]*gJ[e].w0[Njack-1];//M_Ds r0
            
           xp[9]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3]);
           xp[9]=2*sin(xp[9]/2.);
               
           xp[10]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3]);
           xp[10]=2*sin(xp[10]/2.);
    
           for (j=0;j<Njack;j++){ 
               cp2=derivative_xi(9, j , fit_out,  fit_info, xp);
               ck2=derivative_xi(10, j , fit_out,  fit_info, xp);
                if(strcmp(AV,"A")==0)  F[ii][j]=gJ[e].FA[i][j] * gJ[e].ZV[j]; 
                else if(strcmp(AV,"V")==0) F[ii][j]=gJ[e].FV[i][j] * gJ[e].ZA[j];
                else if(strcmp(AV,"{A^H}")==0) F[ii][j]=gJ[e].FA_from_H0[i][j] ;
                else if(strcmp(AV,"{V^H}")==0) F[ii][j]=gJ[e].FV_from_H0[i][j] * gJ[e].ZA[j];
                else if(strcmp(AV,"{V^HA}")==0) F[ii][j]=gJ[e].FV_from_H0_HA[i][j] * gJ[e].ZAV[j];
                else error(0==0,1,"argument of function reweighting_for_plot AV is nor A neither V","");
                F[ii][j]=F[ii][j]-cp2*xp[9]*xp[9]-ck2*xp[10]*xp[10];
                kp[ii][j]=kp_tot[i][j];
           }
           ii=ii+1;
      }
   }}}

   
   x=(double*) malloc(sizeof(double)*(Neff));

   fit=(double**) malloc(sizeof(double*)*(Neff));

   r=(double**) malloc(sizeof(double*)*(npar_fun));
   for(i=0;i<npar_fun;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
  

   chi2=(double*) malloc(sizeof(double)*Njack); 
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(Neff));
        for (i=0;i<Neff;i++){
            y[j][i]=(double*) malloc(sizeof(double)*2);
        }
   }    

   for (i=0;i<Neff;i++){
        fit[i]=mean_and_error_jack_biased(Njack, F[i]);
        //printf("%g   %g   %g\n",head.k[i+head.nk+kmin],fit[i][0],fit[i][1]);
        for (j=0;j<Njack;j++){ 
            y[j][i][0]=F[i][j];
            y[j][i][1]=fit[i][1];
        }
        x[i]=kp[i][Njack-1];
        
   }

   for (j=0;j<Njack;j++){
        tmp=linear_fit( Neff, x, y[j],  npar_fun,fit_function );
        chi2[j]=compute_chisqr( Neff, x, y[j],npar_fun, tmp ,fit_function );
        for(i=0;i<npar_fun;i++){
            r[i][j]=tmp[i][0];
            free(tmp[i]);
        }
        free(tmp);
    }    
    
  
    //////free////
    
    for (i=0;i<Neff;i++){
        free(fit[i]);  free(F[i]);  free(kp[i]);
    }
    free(F), free(kp);
    free(fit);
    free(chi2);
    free(x);
    for (j=0;j<Njack;j++){
            for (i=0;i<Neff;i++){
                free(y[j][i]);
            }
            free(y[j]);
    }
    free(y);
    free(xp);
    return r;

}


double **subtract_non_Lorentz_and_int_K(char **argv, const char  *string,double **phys_point,double xmin,double xmax,  struct header *head ,int Njack ,int npar_fun,double *fit_function(int,double),struct fit_result fit_out, struct fit_type fit_info  , struct reph_jack *gJ,int e,const char *AV){
   double ***y,*x,**r,*chi2,**tmp,**fit; 
   int i,j,ii;  
   int Ntot=head[e].nmoms*head[e].nmoms*head[e].nmoms;
   int Neff=0;
   int imom0,imoms,imomt;
   double *xp=(double*) malloc(sizeof(double)*fit_info.Nvar); 
   int ikt=0 ,iks=1;
   
   ii=0;
   int     i_mDs=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);
   int     i_mD=index_ensemble_twopt_fit(head[e],ikt,0,0,0);
   int     i_mK=index_ensemble_twopt_fit(head[e],0,iks,0,0);
   int     i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
   double cp2,ck2;
   double ***F2=(double***) malloc(sizeof(double**)*Ntot);
   double **kp_tot=(double**) malloc(sizeof(double*)*Ntot);
   
   ii=0;
   for(imom0=0;imom0<head[e].nmoms;imom0++){       
   for(imomt=0;imomt<head[e].nmoms;imomt++){
   for(imoms=0;imoms<head[e].nmoms;imoms++){
        F2[ii]=(double**) malloc(sizeof(double*)*2);
        F2[ii][0]=(double*) malloc(sizeof(double)*Njack);
        F2[ii][1]=(double*) malloc(sizeof(double)*Njack);
        kp_tot[ii]=(double*) malloc(sizeof(double)*Njack);
        ii++;
   }}}
   
   for (int ns=0;ns<2;ns++){
    ii=0   ;
   for(imom0=0;imom0<head[e].nmoms;imom0++){       
   for(imomt=0;imomt<head[e].nmoms;imomt++){
   for(imoms=0;imoms<head[e].nmoms;imoms++){
           i=index_ensemble_twopt_G_fit(head[e],ikt,iks+ns,imom0,imomt,imoms);
           i_mK=index_ensemble_twopt_fit(head[e],0,iks+ns,0,0);
           //F[ii]=(double*) malloc(sizeof(double)*Njack);
           //kp[ii]=(double*) malloc(sizeof(double)*Njack);

           xp[0]=head[e].k[head[e].nk+ikt]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
           xp[1]=gJ[e].w0[Njack-1];//w0
           xp[2]=gJ[e].xG[i][Njack-1];//x_gamma
           xp[3]=gJ[e].M_PS[i_mp][Njack-1]*gJ[e].w0[Njack-1];//x_gamma
                 
           xp[4]=head[e].k[head[e].nk+iks+ns]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ms*w0
           xp[5]=gJ[e].M_PS[i_mK][Njack-1]*gJ[e].w0[Njack-1];//M_K r0
                    
           xp[6]=head[e].k[head[e].nk+ikt]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//mc*w0
           xp[7]=gJ[e].M_PS[i_mD][Njack-1]*gJ[e].w0[Njack-1];//M_D r0
           xp[8]=gJ[e].M_PS[i_mDs][Njack-1]*gJ[e].w0[Njack-1];//M_Ds r0
            
           xp[9]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3]);
           xp[9]=2*sin(xp[9]/2.);
               
           xp[10]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3]);
           xp[10]=2*sin(xp[10]/2.);
    
           for (j=0;j<Njack;j++){ 
               cp2=derivative_xi(9, j , fit_out,  fit_info, xp);
               ck2=derivative_xi(10, j , fit_out,  fit_info, xp);

                if(strcmp(AV,"A")==0)  F2[ii][ns][j]=gJ[e].FA[i][j] * gJ[e].ZV[j]; 
                else if(strcmp(AV,"V")==0) F2[ii][ns][j]=gJ[e].FV[i][j] * gJ[e].ZA[j];
                else if(strcmp(AV,"{A^H}")==0) F2[ii][ns][j]=gJ[e].FA_from_H0[i][j] ;
                else if(strcmp(AV,"{V^H}")==0) F2[ii][ns][j]=gJ[e].FV_from_H0[i][j] * gJ[e].ZA[j];
                else if(strcmp(AV,"{V^HA}")==0) F2[ii][ns][j]=gJ[e].FV_from_H0_HA[i][j] * gJ[e].ZAV[j];
                else error(0==0,1,"argument of function reweighting_for_plot AV is nor A neither V","");
                F2[ii][ns][j]=F2[ii][ns][j]-cp2*xp[9]*xp[9]-ck2*xp[10]*xp[10];
           }
           ii=ii+1;
   }}}
   }
   double xg;
   double E_p,E_k,p,k,kdp,m,b;
   ii=0;
   int is=0,ib=0;
   int i_mK0=index_ensemble_twopt_fit(head[e],0,iks,0,0);
   int i_mK1=index_ensemble_twopt_fit(head[e],0,iks+1,0,0);
   ib=0;
   int i0;
   for(imom0=0;imom0<head[e].nmoms;imom0++){       
   for(imomt=0;imomt<head[e].nmoms;imomt++){
   for(imoms=0;imoms<head[e].nmoms;imoms++){
       i0=index_ensemble_twopt_fit(head[e],0,iks,imom0,imoms);
       for (j=0;j<Njack;j++){ 
            p=2.*sin(  (2.*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3])/2.);
            k=2.*sin(  (2.*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3])/2.);
            E_p=phys_point[j][5]/gJ[e].w0[j];
            E_p*=E_p;
            E_p+=p*p;
            E_p=sqrt(E_p);
            //E_p=gJ[e].M_PS[i0][j];
            E_k=2.*asinh(   sqrt(k*k)/2.   );
            kdp=E_p*E_k- p*k; 
            kp_tot[ib][j]=2*kdp*gJ[e].w0[j]*gJ[e].w0[j]/(phys_point[j][5]*phys_point[j][5]);//  /M_K^2w0^2
            //kp_tot[ib][j]=2.*kdp/(gJ[e].M_PS[i_mK0][j]*gJ[e].M_PS[i_mK0][j]);//  /M_K^2w0^2
       }
       /*i=index_ensemble_twopt_G_fit(head[e],ikt,iks+0,imom0,imomt,imoms);
       printf("1=%g\t",gJ[e].xG[i][Njack-1]);
       i=index_ensemble_twopt_G_fit(head[e],ikt,iks+1,imom0,imomt,imoms);
       printf("2=%g  \t   ave=%g   \n",gJ[e].xG[i][Njack-1],kp_tot[ib][Njack-1]);*/
       if (kp_tot[ib][Njack-1]>xmin &&  kp_tot[ib][Njack-1]<xmax){
           Neff=Neff+1;
       }
       ib++;
   }}}
//   printf("\n\n");
   error(Neff==0,1,"interpolate_FX","no data to interpolate ikt=%d  iks=%d  e=%d",ikt,iks,e);
   double **F=(double**) malloc(sizeof(double*)*Neff);
   double **kp=(double**) malloc(sizeof(double*)*Neff);
   
   ib=0;is=0;
   for(imom0=0;imom0<head[e].nmoms;imom0++){       
   for(imomt=0;imomt<head[e].nmoms;imomt++){
   for(imoms=0;imoms<head[e].nmoms;imoms++){
       if (kp_tot[ib][Njack-1]>xmin &&  kp_tot[ib][Njack-1]<xmax){
          F[is]=(double*) malloc(sizeof(double)*Njack);
          kp[is]=(double*) malloc(sizeof(double)*Njack);
          for (j=0;j<Njack;j++){ 
              m=( F2[ib][0][j]-F2[ib][1][j] )/(   gJ[e].M_PS[i_mK0][j]*gJ[e].M_PS[i_mK0][j]  -  gJ[e].M_PS[i_mK1][j]*gJ[e].M_PS[i_mK1][j]   );
              b=F2[ib][0][j]-m*gJ[e].M_PS[i_mK0][j]*gJ[e].M_PS[i_mK0][j];
              F[is][j]=m*(phys_point[j][5]*phys_point[j][5]/(gJ[e].w0[j]*gJ[e].w0[j]))+b;
              kp[is][j]=kp_tot[ib][j];
          }
          //printf("F1=%g   F2=%g     F_extr=%g\n",F2[ib][0][Njack-1],F2[ib][1][Njack-1],F[is][Njack-1]);
          is=is+1;
          
       }
       free(F2[ib][0]);free(F2[ib][1]);free(F2[ib]);
       free(kp_tot[ib]);
       ib=ib+1;
   }}}
   free(F2);
   free(kp_tot);
   
   x=(double*) malloc(sizeof(double)*(Neff));
   fit=(double**) malloc(sizeof(double*)*(Neff));
   r=(double**) malloc(sizeof(double*)*(npar_fun));
   
   for(i=0;i<npar_fun;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
  
   chi2=(double*) malloc(sizeof(double)*Njack); 
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(Neff));
        for (i=0;i<Neff;i++){
            y[j][i]=(double*) malloc(sizeof(double)*2);
        }
   }    

   for (i=0;i<Neff;i++){
        fit[i]=mean_and_error_jack_biased(Njack, F[i]);
        //printf("%g   %g   %g\n",head.k[i+head.nk+kmin],fit[i][0],fit[i][1]);
        for (j=0;j<Njack;j++){ 
            y[j][i][0]=F[i][j];
            y[j][i][1]=fit[i][1];
        }
        x[i]=kp[i][Njack-1];
        
   }

   for (j=0;j<Njack;j++){
        tmp=linear_fit( Neff, x, y[j],  npar_fun,fit_function );
        chi2[j]=compute_chisqr( Neff, x, y[j],npar_fun, tmp ,fit_function );
        for(i=0;i<npar_fun;i++){
            r[i][j]=tmp[i][0];
            free(tmp[i]);
        }
        free(tmp);
    }    
    
  
    //////free////
    
    for (i=0;i<Neff;i++){
        free(fit[i]);  free(F[i]);  free(kp[i]);
    }
    free(F), free(kp);
    free(fit);
    free(chi2);
    free(x);
    for (j=0;j<Njack;j++){
            for (i=0;i<Neff;i++){
                free(y[j][i]);
            }
            free(y[j]);
    }
    free(y);
    free(xp);
    return r;

}

/////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


double **subtract_non_Lorentz_and_int_Ds(char **argv, const char  *string,double **phys_point,double xmin,double xmax,  struct header *head ,int Njack ,int npar_fun,double *fit_function(int,double),struct fit_result fit_out, struct fit_type fit_info  , struct reph_jack *gJ,int e,const char *AV){
   double ***y,*x,**r,*chi2,**tmp,**fit; 
   int i,j,ii;  
   int Ntot=head[e].nmoms*head[e].nmoms*head[e].nmoms;
   int Neff=0;
   int imom0,imoms,imomt;
   double *xp=(double*) malloc(sizeof(double)*fit_info.Nvar); 
   int ikt=3 ,iks=1;
   
   ii=0;
   int     i_mDs=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);
   int     i_mD=index_ensemble_twopt_fit(head[e],ikt,0,0,0);
   int     i_mK=index_ensemble_twopt_fit(head[e],0,iks,0,0);
   int     i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
   double cp2,ck2;
   double ***F2=(double***) malloc(sizeof(double**)*Ntot);
   double **kp_tot=(double**) malloc(sizeof(double*)*Ntot);
   
   ii=0;
   for(imom0=0;imom0<head[e].nmoms;imom0++){       
   for(imomt=0;imomt<head[e].nmoms;imomt++){
   for(imoms=0;imoms<head[e].nmoms;imoms++){
        F2[ii]=(double**) malloc(sizeof(double*)*4);
        F2[ii][0]=(double*) malloc(sizeof(double)*Njack);
        F2[ii][1]=(double*) malloc(sizeof(double)*Njack);
        F2[ii][2]=(double*) malloc(sizeof(double)*Njack);
        F2[ii][3]=(double*) malloc(sizeof(double)*Njack);
        kp_tot[ii]=(double*) malloc(sizeof(double)*Njack);
        ii++;
   }}}

   for (int nt=0;nt<2;nt++){
   for (int ns=0;ns<2;ns++){
    ii=0   ;
   for(imom0=0;imom0<head[e].nmoms;imom0++){       
   for(imomt=0;imomt<head[e].nmoms;imomt++){
   for(imoms=0;imoms<head[e].nmoms;imoms++){
           i=index_ensemble_twopt_G_fit(head[e],ikt+nt,iks+ns,imom0,imomt,imoms);
           i_mK=index_ensemble_twopt_fit(head[e],0,iks+ns,0,0);
           i_mDs=index_ensemble_twopt_fit(head[e],ikt+nt,iks+ns,0,0);
           i_mD=index_ensemble_twopt_fit(head[e],ikt+nt,0,0,0);
           i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
   
           //F[ii]=(double*) malloc(sizeof(double)*Njack);
           //kp[ii]=(double*) malloc(sizeof(double)*Njack);

           xp[0]=head[e].k[head[e].nk+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
           xp[1]=gJ[e].w0[Njack-1];//w0
           xp[2]=gJ[e].xG[i][Njack-1];//x_gamma
           xp[3]=gJ[e].M_PS[i_mp][Njack-1]*gJ[e].w0[Njack-1];//x_gamma
                 
           xp[4]=head[e].k[head[e].nk+iks+ns]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ms*w0
           xp[5]=gJ[e].M_PS[i_mK][Njack-1]*gJ[e].w0[Njack-1];//M_K r0
                    
           xp[6]=head[e].k[head[e].nk+ikt+nt]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//mc*w0
           xp[7]=gJ[e].M_PS[i_mD][Njack-1]*gJ[e].w0[Njack-1];//M_D r0
           xp[8]=gJ[e].M_PS[i_mDs][Njack-1]*gJ[e].w0[Njack-1];//M_Ds r0
            
           xp[9]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3]);
           xp[9]=2*sin(xp[9]/2.);
               
           xp[10]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3]);
           xp[10]=2*sin(xp[10]/2.);
    
           for (j=0;j<Njack;j++){ 
               cp2=derivative_xi(9, j , fit_out,  fit_info, xp);
               ck2=derivative_xi(10, j , fit_out,  fit_info, xp);

                if(strcmp(AV,"A")==0)  F2[ii][ns+nt*2][j]=gJ[e].FA[i][j] * gJ[e].ZV[j]; 
                else if(strcmp(AV,"V")==0) F2[ii][ns+nt*2][j]=gJ[e].FV[i][j] * gJ[e].ZA[j];
                else if(strcmp(AV,"{A^H}")==0) F2[ii][ns+nt*2][j]=gJ[e].FA_from_H0[i][j] ;
                else if(strcmp(AV,"{V^H}")==0) F2[ii][ns+nt*2][j]=gJ[e].FV_from_H0[i][j] * gJ[e].ZA[j];
                else if(strcmp(AV,"{V^HA}")==0) F2[ii][ns+nt*2][j]=gJ[e].FV_from_H0_HA[i][j] * gJ[e].ZAV[j];
                else error(0==0,1,"argument of function reweighting_for_plot AV is nor A neither V","");
                F2[ii][ns+nt*2][j]=F2[ii][ns+nt*2][j]-cp2*xp[9]*xp[9]-ck2*xp[10]*xp[10];
           }
           ii=ii+1;
   }}}
   }}

   double xg;
   double E_p,E_k,p,k,kdp,m,b;
   ii=0;
   int is=0,ib=0;
   int i_mK0=index_ensemble_twopt_fit(head[e],0,iks,0,0);
   int i_mK1=index_ensemble_twopt_fit(head[e],0,iks+1,0,0);
   int i_mD0=index_ensemble_twopt_fit(head[e],3,0,0,0);
   int i_mD1=index_ensemble_twopt_fit(head[e],4,0,0,0);

   ib=0;
   int i0;
   for(imom0=0;imom0<head[e].nmoms;imom0++){       
   for(imomt=0;imomt<head[e].nmoms;imomt++){
   for(imoms=0;imoms<head[e].nmoms;imoms++){
       i0=index_ensemble_twopt_fit(head[e],0,iks,imom0,imoms);
       for (j=0;j<Njack;j++){ 
            p=2.*sin(  (2.*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3])/2.);
            k=2.*sin(  (2.*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3])/2.);
            E_p=phys_point[j][8]/gJ[e].w0[j];
            E_p*=E_p;
            E_p+=p*p;
            E_p=sqrt(E_p);
            //E_p=gJ[e].M_PS[i0][j];
            E_k=2.*asinh(   sqrt(k*k)/2.   );
            kdp=E_p*E_k- p*k; 
            kp_tot[ib][j]=2*kdp*gJ[e].w0[j]*gJ[e].w0[j]/(phys_point[j][8]*phys_point[j][8]);//  /M_K^2w0^2
            //kp_tot[ib][j]=2.*kdp/(gJ[e].M_PS[i_mK0][j]*gJ[e].M_PS[i_mK0][j]);//  /M_K^2w0^2
       }
       /*i=index_ensemble_twopt_G_fit(head[e],ikt,iks+0,imom0,imomt,imoms);
       printf("1=%g\t",gJ[e].xG[i][Njack-1]);
       i=index_ensemble_twopt_G_fit(head[e],ikt,iks+1,imom0,imomt,imoms);
       printf("2=%g  \t   ave=%g   \n",gJ[e].xG[i][Njack-1],kp_tot[ib][Njack-1]);*/
       if (kp_tot[ib][Njack-1]>xmin &&  kp_tot[ib][Njack-1]<xmax){
           Neff=Neff+1;
       }
       ib++;
   }}}
   error(Neff==0,1,"interpolate_FX","no data to interpolate ikt=%d  iks=%d  e=%d",ikt,iks,e);
   double **F=(double**) malloc(sizeof(double*)*Neff);
   double **kp=(double**) malloc(sizeof(double*)*Neff);
   double m1,b1,FD0,FD1;
   ib=0;is=0;
   for(imom0=0;imom0<head[e].nmoms;imom0++){       
   for(imomt=0;imomt<head[e].nmoms;imomt++){
   for(imoms=0;imoms<head[e].nmoms;imoms++){
       if (kp_tot[ib][Njack-1]>xmin &&  kp_tot[ib][Njack-1]<xmax){
          F[is]=(double*) malloc(sizeof(double)*Njack);
          kp[is]=(double*) malloc(sizeof(double)*Njack);
          for (j=0;j<Njack;j++){ 
              m=( F2[ib][0][j]-F2[ib][1][j] )/(   gJ[e].M_PS[i_mK0][j]*gJ[e].M_PS[i_mK0][j]  -  gJ[e].M_PS[i_mK1][j]*gJ[e].M_PS[i_mK1][j]   );
              b=F2[ib][0][j]-m*gJ[e].M_PS[i_mK0][j]*gJ[e].M_PS[i_mK0][j];
              FD0=m*(phys_point[j][5]*phys_point[j][5]/(gJ[e].w0[j]*gJ[e].w0[j]))+b;
              
              m1=( F2[ib][2][j]-F2[ib][3][j] )/(   gJ[e].M_PS[i_mK0][j]*gJ[e].M_PS[i_mK0][j]  -  gJ[e].M_PS[i_mK1][j]*gJ[e].M_PS[i_mK1][j]   );
              b1=F2[ib][2][j]-m*gJ[e].M_PS[i_mK0][j]*gJ[e].M_PS[i_mK0][j];
              FD1=m*(phys_point[j][5]*phys_point[j][5]/(gJ[e].w0[j]*gJ[e].w0[j]))+b;

              m=( FD0-FD1 )/(   gJ[e].M_PS[i_mD0][j]*gJ[e].M_PS[i_mD0][j]  -  gJ[e].M_PS[i_mD1][j]*gJ[e].M_PS[i_mD1][j]   );
              b=FD0-m*gJ[e].M_PS[i_mD0][j]*gJ[e].M_PS[i_mD0][j];              
              F[is][j]=m*(phys_point[j][7]*phys_point[j][7]/(gJ[e].w0[j]*gJ[e].w0[j]))+b;
              
              
              
              kp[is][j]=kp_tot[ib][j];
          }
          //printf("F1=%g   F2=%g     F_extr=%g\n",F2[ib][0][Njack-1],F2[ib][1][Njack-1],F[is][Njack-1]);
          is=is+1;
          
       }
       free(F2[ib][0]);free(F2[ib][1]);free(F2[ib][2]);free(F2[ib][3]); free(F2[ib]);
       free(kp_tot[ib]);
       ib=ib+1;
   }}}
   free(F2);free(kp_tot);
   
   x=(double*) malloc(sizeof(double)*(Neff));
   fit=(double**) malloc(sizeof(double*)*(Neff));
   r=(double**) malloc(sizeof(double*)*(npar_fun));
   
   for(i=0;i<npar_fun;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
  
   chi2=(double*) malloc(sizeof(double)*Njack); 
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(Neff));
        for (i=0;i<Neff;i++){
            y[j][i]=(double*) malloc(sizeof(double)*2);
        }
   }    

   for (i=0;i<Neff;i++){
        fit[i]=mean_and_error_jack_biased(Njack, F[i]);
        //printf("%g   %g   %g\n",head.k[i+head.nk+kmin],fit[i][0],fit[i][1]);
        for (j=0;j<Njack;j++){ 
            y[j][i][0]=F[i][j];
            y[j][i][1]=fit[i][1];
        }
        x[i]=kp[i][Njack-1];
        
   }

   for (j=0;j<Njack;j++){
        tmp=linear_fit( Neff, x, y[j],  npar_fun,fit_function );
        chi2[j]=compute_chisqr( Neff, x, y[j],npar_fun, tmp ,fit_function );
        for(i=0;i<npar_fun;i++){
            r[i][j]=tmp[i][0];
            free(tmp[i]);
        }
        free(tmp);
    }    
    
  
    //////free////
    
    for (i=0;i<Neff;i++){
        free(fit[i]);  free(F[i]);  free(kp[i]);
    }
    free(F), free(kp);
    free(fit);
    free(chi2);
    free(x);
    for (j=0;j<Njack;j++){
            for (i=0;i<Neff;i++){
                free(y[j][i]);
            }
            free(y[j]);
    }
    free(y);
    free(xp);
    return r;

}
