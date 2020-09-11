#define various_fit_C
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"

#include "resampling.hpp"
#include "gnuplot.hpp"
#include "linear_fit.hpp"
#include "indices.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"

double *polynimial_degree_n(int np1, double in){
    double *r;
    int i,j;
    
    r=(double*) malloc(sizeof(double)*(np1));

    for (i=0;i<np1;i++){
        r[i]=1.;
        for (j=0;j<i;j++){
            r[i]*=in;
            
        }
    }

    return r;
}

double *poles_degree_n(int np1, double in){
    double *r;
    int i,j;
    
    r=(double*) malloc(sizeof(double)*(np1));

    for (i=0;i<np1;i++){
        r[i]=1.;
        for (j=0;j<i;j++){
            r[i]/=in;
            
        }
    }

    return r;
}


double  **fit_polynomial(char **argv,const char* string, int n, int kmin,int ik1, char* name, double ***mass_jack_fit_GEVP,int  Njack , FILE *outfile ){
   double ***y,*x,**r,*chi2,**tmp,**rm,*chi2m,**fit; 
   int i,j;  

   error(kmin<ik1, 1 , "various_fit", "kmin must be smaller than ik1");
   
   x=(double*) malloc(sizeof(double)*(file_head.nk-kmin));

   chi2m=(double*) malloc(sizeof(double)*(n+1));
   rm=(double**) malloc(sizeof(double*)*(n+1));
   fit=(double**) malloc(sizeof(double*)*(file_head.nk-kmin));

   r=(double**) malloc(sizeof(double*)*(n+1));
   for(i=0;i<=n;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(file_head.nk-kmin));
        for (i=0;i<file_head.nk-kmin;i++){
            y[j][i]=(double*) malloc(sizeof(double)*2);
        }
   }
   for (i=0;i<file_head.nk-kmin;i++){
        fit[i]=mean_and_error_jack(Njack, mass_jack_fit_GEVP[i+kmin][ik1]);
        printf("%g   %g   %g\n",file_head.k[i+file_head.nk+kmin],fit[i][0],fit[i][1]);
        for (j=0;j<Njack;j++){ 
            y[j][i][0]=mass_jack_fit_GEVP[i+kmin][ik1][j];
            y[j][i][1]=fit[i][1];
        }
        x[i]=file_head.k[i+file_head.nk+kmin];
        
   }
   
   for (j=0;j<Njack;j++){
        tmp=linear_fit( file_head.nk-kmin, x, y[j],  n+1,polynimial_degree_n );

        chi2[j]=compute_chisqr( file_head.nk-kmin, x, y[j],n, tmp ,polynimial_degree_n );
        for(i=0;i<=n;i++){
            r[i][j]=tmp[i][0];
            free(tmp[i]);
        }
        free(tmp);
    }    
    
    for(i=0;i<=n;i++){
        rm[i]=mean_and_error_jack(Njack, r[i]);
    }
    chi2m=mean_and_error_jack(Njack, chi2);
    
    plotting_fit_deg_n_pdf(argv, string,file_head.nk-kmin, n, x,fit , rm,chi2m );
    //////free////
    for(i=0;i<=n;i++){
        free(r[i]);       
    }
    for (i=0;i<file_head.nk-kmin;i++)
        free(fit[i]);
    free(fit);
    free(r);free(chi2);
    free(x);
    for (j=0;j<Njack;j++){
            for (i=0;i<file_head.nk-kmin;i++){
                free(y[j][i]);
            }
            free(y[j]);
    }
    free(y);
    return rm;
}


double **fit_FX(char **argv, const char  *string,int ikt,int iks,double **kp_tot, double **FX,  struct header file_head ,int Njack ,int npar_fun,double *fit_function(int,double) ){
   double ***y,*x,**r,*chi2,**tmp,**rm,*chi2m,**fit; 
   int i,j,ii;  
   int Ntot=file_head.nmoms*file_head.nmoms*file_head.nmoms;
   int Neff=0;
   int imom0,imoms,imomt;
   char name_RA[500],name_RV[500];
   sprintf(name_RA,"FA-ikt%d-iks%d",ikt,iks); 
   sprintf(name_RV,"FV-ikt%d-iks%d",ikt,iks); 

   for(imom0=0;imom0<file_head.nmoms;imom0++){       
   for(imomt=0;imomt<file_head.nmoms;imomt++){
   for(imoms=0;imoms<file_head.nmoms;imoms++){
       i=index_n_twopt_G_fit(ikt,iks,imom0,imomt,imoms);
       if (kp_tot[i][Njack-1]>1e-6){
           Neff=Neff+1;
       }
   }}}
   
   double **F=(double**) malloc(sizeof(double*)*Neff);
   double **kp=(double**) malloc(sizeof(double*)*Neff);
   
   ii=0;
   for(imom0=0;imom0<file_head.nmoms;imom0++){       
   for(imomt=0;imomt<file_head.nmoms;imomt++){
   for(imoms=0;imoms<file_head.nmoms;imoms++){
      i=index_n_twopt_G_fit(ikt,iks,imom0,imomt,imoms);
      if (kp_tot[i][Njack-1]>1e-6){
           F[ii]=(double*) malloc(sizeof(double)*Njack);
           kp[ii]=(double*) malloc(sizeof(double)*Njack);
           for (j=0;j<Njack;j++){ 
                F[ii][j]=FX[i][j];
                kp[ii][j]=kp_tot[i][j];
           }
           ii=ii+1;
       }
   }}}
   
              

   x=(double*) malloc(sizeof(double)*(Neff));

   chi2m=(double*) malloc(sizeof(double)*(npar_fun));
   rm=(double**) malloc(sizeof(double*)*(npar_fun));
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
        fit[i]=mean_and_error_jack(Njack, F[i]);
        //printf("%g   %g   %g\n",file_head.k[i+file_head.nk+kmin],fit[i][0],fit[i][1]);
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
  
    for(i=0;i<npar_fun;i++){
        rm[i]=mean_and_error_jack(Njack, r[i]);
    }
    chi2m=mean_and_error_jack(Njack, chi2);
    
    if ( strcmp(argv[5],"only_fit_R")==0){
        plotting_fit_deg_n_pdf(argv, string,Neff, npar_fun-1, x,fit , rm,chi2m );///!!!!! put nparameters-1   because  it need the order of the polynomial
        /*if (strcmp(name_RV,string)==0)
            plotting_fit_deg_n_pdf(argv, string,Neff, npar_fun-1, x,fit , rm,chi2m );///!!!!! put nparameters-1   because  it need the order of the polynomial
        else if (strcmp(name_RA,string)==0)
            plotting_fit_poles_deg_n_pdf(argv, string,Neff, npar_fun-1, x,fit , rm,chi2m );///!!!!! put nparameters-1   because  it need the order of the polynomial
            */
    }
    //////free////
    for(i=0;i<npar_fun;i++){
        free(rm[i]);       
    }
    for (i=0;i<Neff;i++){
        free(fit[i]);  free(F[i]);  free(kp[i]);
    }
    free(F), free(kp);
    free(fit);
    free(rm);
    free(chi2);
    free(x);
    for (j=0;j<Njack;j++){
            for (i=0;i<Neff;i++){
                free(y[j][i]);
            }
            free(y[j]);
    }
    free(y);
    return r;

}




double **interpolate_FX(char **argv, const char  *string,int ikt,int iks,double **kp_tot,double xmin,double xmax, double **FX,  struct header head ,int Njack ,int npar_fun,double *fit_function(int,double) ){
   double ***y,*x,**r,*chi2,**tmp,**fit; 
   int i,j,ii;  
   int Ntot=head.nmoms*head.nmoms*head.nmoms;
   int Neff=0;
   int imom0,imoms,imomt;
   char name_RA[500],name_RV[500];
   sprintf(name_RA,"FA-ikt%d-iks%d",ikt,iks); 
   sprintf(name_RV,"FV-ikt%d-iks%d",ikt,iks); 

   for(imom0=0;imom0<head.nmoms;imom0++){       
   for(imomt=0;imomt<head.nmoms;imomt++){
   for(imoms=0;imoms<head.nmoms;imoms++){
       i=index_ensemble_twopt_G_fit(head,ikt,iks,imom0,imomt,imoms);
       if (kp_tot[i][Njack-1]>xmin &&  kp_tot[i][Njack-1]<xmax){
           Neff=Neff+1;
       }
   }}}
   
   error(Neff==0,1,"interpolate_FX","no data to interpolate ikt=%d  iks=%d",ikt,iks);
       
   double **F=(double**) malloc(sizeof(double*)*Neff);
   double **kp=(double**) malloc(sizeof(double*)*Neff);
   ii=0;
   for(imom0=0;imom0<head.nmoms;imom0++){       
   for(imomt=0;imomt<head.nmoms;imomt++){
   for(imoms=0;imoms<head.nmoms;imoms++){
      i=index_ensemble_twopt_G_fit(head,ikt,iks,imom0,imomt,imoms);
      if (kp_tot[i][Njack-1]>xmin &&  kp_tot[i][Njack-1]<xmax){
           F[ii]=(double*) malloc(sizeof(double)*Njack);
           kp[ii]=(double*) malloc(sizeof(double)*Njack);
           for (j=0;j<Njack;j++){ 
                F[ii][j]=FX[i][j];
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
        fit[i]=mean_and_error_jack(Njack, F[i]);
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
    return r;

}


