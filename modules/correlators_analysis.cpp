#define correlators_analysis_C

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
#include "correlators_analysis.hpp"
//#include "eigensystem.hpp"



double  sinh_mass(int n, int Nvar, double *x,int Npar,double  *P){
    
    double t=P[0];
    double T=P[1];
    double cp=sinh(x[0]*(t+1.+(1-T)/2));
    double c=sinh(x[0]*(t+(1-T)/2));
    //printf("inside sinh = %f   %f \n", c,cp);
    return c/cp;
}

//in [varialbe] [t] [re_im]
double two_particle_energy(int t,int T , double **in){
    double c[3];
    double P[2]; 
    double x[1];
    
    P[0]=t;
    P[1]=T;
    c[0]=in[t][0]-in[t+1][0];
    c[1]=in[t+1][0]-in[t+2][0];
    double mass=log(c[0]/c[1]);
    //printf(" %d %g  %g   mass=%g\n", t,c[0],c[1],mass);
    /*double r=rtbis_func_eq_input(sinh_mass ,
                                     0 , 1, x ,//double n [n functions], Nvar, varibles
                                     0, P , 0, // fit_info.Npar, Parameters, ivar 
                                     c[0]/c[1], mass-0.1, mass+0.1, mass*0.0001);
    */
    c[0]=in[t-1][0]-in[t][0];
    c[1]=in[t][0]-in[t+1][0];
    c[2]=in[t+1][0]-in[t+2][0];
    return acosh(  (c[0]+c[2])/ (2.*c[1])  );
    //rp/r = sin(E2(t+1+(1-T/2)) /sin(E2(t+(1-T/2))
    
}


double M_eff_T(  int t, int T, double **in){
    double mass;
 
    double ct[1],ctp[1],res,tmp_mass, u,d ;
    int i,L0;
    
    
    ct[0]=in[t][0];
    ctp[0]=in[t+1][0];


    mass=log(ct[0]/ctp[0]);

    res=1;
    i=t;
    while(res>1e-12){
             u=1.+exp(-mass*(T-2.*i-2.));
             d=1.+exp(-mass*(T-2.*i));
             tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
             res=fabs(tmp_mass - mass);
             mass=tmp_mass;
    }
  
    return mass;

}


double M_eff_sinh_T(  int t, int T, double **in){
    double mass;
 
    double ct[1],ctp[1],res,tmp_mass, u,d ;
    int i,L0;
    
    L0=file_head.l0;
    ct[0]=in[t][0]-in[t+1][0];
    ctp[0]=in[t+1][0]-in[t+2][0];


    mass=log(ct[0]/ctp[0]);

    res=1;
    i=t;
    while(res>1e-12){
             u=-1.+exp(-mass*(L0-1-2.*i-2.));
             d=-1.+exp(-mass*(L0-1-2.*i));
             tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
             res=fabs(tmp_mass - mass);
             mass=tmp_mass;
    }
  
    return mass;

}



double   *plateau_correlator_function(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,FILE **plateaux_masses,FILE *outfile,  int index , const char *description , double (*fun)(int ,int  , double ** )){
   int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(*plateaux_masses,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);


   fprintf(outfile,"#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
              //shift 
              r[i][j]=fun( i,file_head.l0 , conf_jack[j][index]);
             
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);

   }

   fit=fit_plateaux(option, kinematic_2pt ,  name,description/*"M_{PS}^{ll}"*/,mt,r,  Njack,*plateaux_masses,outfile);
   write_jack_bin(Njack,fit,file_jack.M_PS);
     
   for(i=1;i<file_head.l0/2;i++)
       free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
       free(r[i]);
   free(r);

   fflush(outfile);
   
    if ( strcmp(option[1],"read_plateaux")==0 ){
     fclose(*plateaux_masses);
     *plateaux_masses=open_file(kinematic_2pt.plateau_m_ll,"r");

    }
    return fit;    
    
}
