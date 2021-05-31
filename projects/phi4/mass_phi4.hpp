#ifndef mass_phi4_H
#define mass_phi4_H

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
#include <cmath>

extern "C" { 
    #include "../external/rzeta/src/dzeta_function.h"
}

double energy_CM(  double E  , int  *p,int L);

double phase_shift(double E2,double mass, int  *dvec,int L ){
    double delta;
    double z[2];
    
    
    double ECM=energy_CM(E2,dvec,L);
    double k=sqrt(ECM*ECM/4. - mass*mass);
    double q=k*L/(2*pi_greco);
    double gamma=E2/ECM;
    double A=1.;
// printf("k=%g gamma=%g  mass=%g L=%d\n",k,gamma,mass,L);
    //dzeta_function(z,  q*q,0 , 0, dvec, gamma, A, 0.000001, 1.e10,5);
    //printf("q=%g   dvec=(%d,%d,%d)    gamma=%g   E2=%.8g  ECM=%.8g  k=%g L=%d  mass=%.8g\n",q,dvec[0] ,dvec[1] ,dvec[2], gamma,E2,ECM,k,L,mass);
    dzeta_function(z,  q*q,0 , 0, dvec, gamma, A, 1.e-3, 1.e6 ,3);
    std::complex<double>  zc(z[0],z[1]);
    delta=real(  std::atan(  (pow(pi_greco,3./2.) *q*gamma )/zc   ));
    
    
    return delta;
}

double energy_CM(  double E  , int  *p,int L){
    
    double normp=p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
    normp*=2*pi_greco/((double) L);
    normp*=2*pi_greco/((double) L);
    return sqrt(E*E- normp);
}

void phase_shift(double *E2,double *mass, int  *dvec,int L,FILE *outfile, int Njack, char *jackboot ){
    fprintf(outfile,"#E2_CM   err  E2_CM/M   err   re(delta)   err\n " );
    double *E2_CM=(double*) malloc(sizeof(double)*Njack);
    //int dvec[3]= {1,1,1};
    double *k=(double*) malloc(sizeof(double)*Njack);
    for (int j=0;j< Njack;j++){
        E2_CM[j]=energy_CM(E2[j],dvec,L);
        k[j]=sqrt(E2_CM[j]*E2_CM[j]/4. -mass[j]*mass[j]);
    }
    fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(jackboot,Njack,E2_CM ) );
    for (int j=0;j< Njack;j++)
        E2_CM[j]/=mass[j];
    fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(jackboot,Njack,E2_CM ) );
    double *delta=(double*) malloc(sizeof(double)*Njack);
    double *kcotd=(double*) malloc(sizeof(double)*Njack);
    for (int j=0;j< Njack;j++){
        delta[j]=phase_shift( E2[j], mass[j],dvec, L );
        //k[j]=sqrt(E2[j]*E2[j]/4.-mass[j]*mass[j]);
        kcotd[j]=k[j]/std::tan(delta[j]);
    }
    fprintf(outfile,"%.12g  %.12g %.12g  %.12g    %.12g  %.12g\n ", delta[Njack-1],error_jackboot(jackboot,Njack,delta )  , k[Njack-1],  error_jackboot(jackboot,Njack,k ),       kcotd[Njack-1],  error_jackboot(jackboot,Njack,kcotd ));
    free(delta);free(kcotd);
    free(E2_CM);free(k);
    
}

void E3_print_extra(double *E3,double *mass, int  *dvec,int L,FILE *outfile, int Njack, char *jackboot){
    fprintf(outfile,"#E2_CM   err  E2_CM/M   err  k errk \n " );
    double *E3_CM=(double*) malloc(sizeof(double)*Njack);
    //int dvec[3]= {1,1,1};
    double *k=(double*) malloc(sizeof(double)*Njack);
    for (int j=0;j< Njack;j++){
        E3_CM[j]=energy_CM(E3[j],dvec,L);
        k[j]=sqrt(E3_CM[j]*E3_CM[j]/4. -mass[j]*mass[j]);
    }
    fprintf(outfile,"%.12g  %.12g  \t ", E3_CM[Njack-1],error_jackboot(jackboot,Njack,E3_CM ) );
    for (int j=0;j< Njack;j++)
        E3_CM[j]/=mass[j];
    fprintf(outfile,"%.12g  %.12g  \t ", E3_CM[Njack-1],error_jackboot(jackboot,Njack,E3_CM ) );
    
    fprintf(outfile,"\n");
}




double quantization_condition_2particle(int n  , int Nvar , double* x,int Npar,double*P){
    
    double L=P[3];
    double DE= P[2] -P[0]-P[1];
    double mu01=  P[0]*P[1]/(P[0]+P[1]); 
    
    double r=-(2.*pi_greco * x[0])  / (mu01*L*L*L);
    r *=  1 - 2.837297  *(x[0]/L) + 6.375183 *(x[0]/L) *(x[0]/L);
    
    return r;
    
}


double *scattering_len_luscher( int Njack, double *mass0, double *mass1, double *E2 , int L){
    double *a=(double* ) malloc(sizeof(double)*Njack);
    for(int j=0;j<Njack;j++){
        double deltaE2=E2[j]-mass0[j]-mass1[j];
        double *P=(double*) malloc(sizeof(double)*4);
        P[0]=mass0[j];
        P[1]=mass1[j];
        P[2]=E2[j];
        P[3]=(double) L;
        double x[1]={0};
        
        double raw_a=L*L*L *(mass0[j]*mass1[j]/(mass0[j]+mass1[j]) ) *(deltaE2)/(2*pi_greco);
        double xmin=-raw_a-0.2;
        double xmax=-raw_a+0.2;
        a[j]= rtbis_func_eq_input(quantization_condition_2particle, 313 , 1, x ,4, P, 0 , deltaE2, xmin, xmax, 1e-7);
        //a[j]=deltaE2;
        free(P);
    }
   
    return a;
}


double *mass_gamma(int var, int order,int flow ,double *ah){
    double *r=(double*) calloc((1),sizeof(double)); 
    
    //use flow at time of the correlator
    // we need to pass to M_eff a correlator so we create a double **c with has only the correlator at time=flow and flow+1
    double **c=double_malloc_2(flow+2,2);
    c[flow][0]=ah[0];
    c[flow+1][0]=ah[1];
    r[0]=M_eff(flow,c);
    free_2(flow+2,c);
    return r;
}


double   *effective_mass_phi4(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,FILE **plateaux_masses,FILE *outfile,  int index , const char *description ){
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
              r[i][j]=M_eff(i,conf_jack[j][index]);
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);

   }

   char tmp_name[NAMESIZE];
   mysprintf(tmp_name,NAMESIZE,option[1]);
   mysprintf(option[1],NAMESIZE,"see");
   fit=fit_plateaux(option, kinematic_2pt ,  name,description/*"M_{PS}^{ll}"*/,mt,r,  Njack,*plateaux_masses,outfile);
   mysprintf(option[1],NAMESIZE,tmp_name);
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


//double   *effective_mass_phi4_gamma(char **option ,struct kinematic kinematic_2pt , char* name, double ****data, int Confs ,FILE **plateaux_masses,FILE *outfile,  int index , const char *description ){
void effective_mass_phi4_gamma(char **option ,struct kinematic kinematic_2pt , char* name, double ****data, int Confs ,const char *plateaux_masses,FILE *outfile,  int index , const char *description ){

    //int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   //if ( strcmp(option[1],"read_plateaux")==0 )
   //	go_to_line(*plateaux_masses,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Confs);
   //mt=(double**) malloc(sizeof(double*)*file_head.l0);

   
    

   fprintf(outfile,"#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(i=1;i<file_head.l0/2;i++){  
           double *datag;
           // 2 component ot compute the mass
           datag=(double *) malloc(sizeof(double) * Confs *2 ); 
           for (j=0;j<Confs;j++){
              //store in the format for the gamma analysis
              //two variables c(t) and c(t+1), order=0 
              datag[j*2]=data[j][index][i][0];
              datag[1+j*2]=data[j][index][i+1][0];
           }
           
           double *obs=analysis_gamma(  2 , 1, Confs, i//time
                                     , datag  , mass_gamma);
          // mt[i][0]=obs[0];
          // mt[i][0]=obs[1];
           fprintf(outfile,"%d   %.15e    %.15e   %.15e  %.15e  %.15e\n",i,obs[0],obs[1],obs[2],obs[3],obs[4]);
           free(obs);
           free(datag);

   }

   //fit=fit_plateaux(option, kinematic_2pt ,  name,description/*"M_{PS}^{ll}"*/,mt,r,  Confs,*plateaux_masses,outfile);
  // write_jack_bin(Confs,fit,file_jack.M_PS);
     
  /* for(i=1;i<file_head.l0/2;i++)
       free(mt[i]);
   free(mt);*/
   for(i=0;i<file_head.l0;i++)
       free(r[i]);
   free(r);

   fflush(outfile);
   
    /*if ( strcmp(option[1],"read_plateaux")==0 ){
     fclose(*plateaux_masses);
     *plateaux_masses=open_file(kinematic_2pt.plateau_m_ll,"r");

    }*/
    
    //return fit;    
    
}


#endif
