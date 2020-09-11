#define m_eff1_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "linear_fit.hpp"
#include "mutils.hpp"
#include "resampling.hpp"
#include "m_eff1.hpp"
#include "gnuplot.hpp"
#include "global.hpp"
#include "eigensystem.hpp"


double M_eff1(  int t, double **in){
    double mass;
 
    double ct[1],ctp[1],res,tmp_mass, u,d ;
    int i,L0;
    
    L0=file_head.l0;
    ct[0]=in[t][0];
    ctp[0]=in[t+1][0];


    mass=log(ct[0]/ctp[0]);

    res=1;
    i=t;
    while(res>1e-12){
             u=1.+exp(-mass*(L0-2.*i-2.));
             d=1.+exp(-mass*(L0-2.*i));
             tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
             res=fabs(tmp_mass - mass);
             mass=tmp_mass;
    }
  /*  if(t==0)
        mass=acosh((in[L0-1][0]+in[t+1][0])/(2.*in[t][0]));
    else
        mass=acosh((in[t-1][0]+in[t+1][0])/(2.*in[t][0]));
    */
    return mass;

}


void replace_plateau(struct kinematic kinematic_2pt ,struct observable *obs ,const char* description,int tmin,int tmax,int sep){
    char instruction[NAMESIZE];
    int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1)+1;
    mysprintf(instruction,NAMESIZE,"sed -i \"%ds/.*/%d  %d  %d/\"  %s",line,tmin,tmax,sep,(*obs).name_plateaux);

    system(instruction);
    
}

double *fit_plateaux(char **option,struct kinematic kinematic_2pt , char* name,const char *description,double **mt,double **r, int Njack,struct observable *obs)
{
   int yn;
   int tmin=15, tmax=20, sep=1;
   double *m,*fit;
   char jack_name[NAMESIZE];
   double *chi2;
   
   
   yn=1;
   if ( strcmp(option[1],"see")==0 ){
            while(yn>0){
                printf("#%s 1) mu %g r %d theta %g 2) mu %g r %d theta %g \n",
                       description,
                kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
                kinematic_2pt.k1,kinematic_2pt.r1,kinematic_2pt.mom1);   
                plotting( file_head.l0, mt , &tmin,&tmax, &sep);
                m=try_linear_fit(option ,tmin, tmax, sep , mt, r, Njack,&chi2 );
                yn=plotting_fit(file_head.l0, mt , tmin,tmax,m,chi2);
                free(m);free(chi2);
            }
   }
   if ( strcmp(option[1],"read_plateaux")==0 ){
        fscanf((*obs).f_plateaux,"%d  %d  %d\n",&tmin,&tmax,&sep);
        //m=try_linear_fit(option, tmin,  tmax,sep , mt, r, Njack );    
        m=try_linear_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2);
        yn=1;
        if(chi2[0]>5.5){
            free(m);free(chi2);
             while(yn>0){
                   m=try_linear_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2);
                   printf("#%s from %s 1) mu %g r %d theta %g 2) mu %g r %d theta %g  line_plateaux=%d\n",
                        description,name,
                    kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
                    kinematic_2pt.k1,kinematic_2pt.r1,kinematic_2pt.mom1,kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1)+1 );   
                    printf("The chi2 in the range [%d,%d] is %g  consider changing time interval\n",tmin,tmax,chi2[0]);
                    yn=plotting_fit(file_head.l0, mt , tmin,tmax,m,chi2);
                    if (yn>0){
                        plotting( file_head.l0, mt , &tmin,&tmax, &sep);
                        free(m);free(chi2);
                    }
            }
            replace_plateau(  kinematic_2pt,obs ,  description,tmin,tmax,sep);

         }
          free(m);  free(chi2);


   }
   m=try_linear_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2);              
   fit=give_jack_linear_fit( tmin,  tmax,sep , mt, r, Njack );    
   
    fprintf((*obs).f_out,"\n\n #%s fit in [%d,%d] chi2=%.5f\n  %.15g    %.15g\n\n\n",description,tmin,tmax,chi2[0],m[0],m[1]);
   printf("#%s (mu_h=%.4f, mu_l=%.4f) fit in [%d,%d]:  %.15g    %.15g\n",description,kinematic_2pt.k2,kinematic_2pt.k1 ,tmin,tmax,m[0],m[1]);

   if ( strcmp(option[5],"pdf")==0 ){
           plotting_fit_pdf(option,description,file_head.l1, file_head.l0,file_head.beta,file_head.ksea,file_head.musea, mt , tmin,tmax,m,name, kinematic_2pt );
   }
    
   free(m);free(chi2);
   return fit;
}



double   *compute_effective_mass1(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,  int index , struct observable *obs , const char *description){
    
   int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
        go_to_line((*obs).f_plateaux,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);


   fprintf((*obs).f_out,"#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
              r[i][j]=M_eff1(i,conf_jack[j][index]);
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           fprintf((*obs).f_out,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);

   }

   fit=fit_plateaux(option, kinematic_2pt ,  name,description,mt,r,  Njack,obs);
   write_jack_bin(Njack,fit,(*obs).name_jack);
     
   for(i=1;i<file_head.l0/2;i++)
       free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
       free(r[i]);
   free(r);

   fflush((*obs).f_out);
   
    if ( strcmp(option[1],"read_plateaux")==0 ){
     fclose((*obs).f_plateaux);
     (*obs).f_plateaux=open_file((*obs).name_plateaux,"r");

    }
    return fit;    
    
}
