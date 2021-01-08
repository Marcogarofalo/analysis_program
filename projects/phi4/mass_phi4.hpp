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
void effective_mass_phi4_gamma(char **option ,struct kinematic kinematic_2pt , char* name, double ****data, int Confs ,FILE **plateaux_masses,FILE *outfile,  int index , const char *description ){

    //int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   //if ( strcmp(option[1],"read_plateaux")==0 )
   //	go_to_line(*plateaux_masses,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Confs);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);

   
    

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