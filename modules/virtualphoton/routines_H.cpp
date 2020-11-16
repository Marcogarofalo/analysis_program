#define routines_H_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "linear_fit.hpp"
#include "mutils.hpp"
#include "resampling.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "global.hpp"



double   *H_AV(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1,  double* mass_rest, double *oPp, int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
{
   int line=kinematic_2pt_G.i;
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);

   double **r,*m,**mt,*fit;
   int i,j,yn,ii;
   char label[NAMESIZE];
   double kp,xG; 
   
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);
      
    double  E_gT=kinematic_2pt_G.E_gT;
    double  E_gT0=0.5;  
      
   

   fprintf(outfile,"##m_eff_H(t) from %s  propagators:1) mu %.5f r %d theta0 %.5f thetat %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
           kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );
   kp=mass_jack_fit_k2k1[Njack-1]*kinematic_2pt_G.E_g- kinematic_2pt_G.kp;
   xG=2*kp/(mass_rest[Njack-1]*mass_rest[Njack-1]);
   fprintf(outfile,"## E_g=%g      E_gT=%g     E=%g     xG=%g  index=%d\n",kinematic_2pt_G.E_g,kinematic_2pt_G.E_gT,mass_jack_fit_k2k1[Njack-1],xG ,line);
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){            
              ii=sym[index];
              r[i][j]=(conf_jack[j][index][i][ii]);
              r[i][j]/=exp(-i*mass_jack_fit_k2k1[j]-(file_head.l0/2-i)*kinematic_2pt_G.E_g );
              r[i][j]/=oPp[j];
              r[i][j]*=2*mass_jack_fit_k2k1[j];
           }
           
   }
   //double norm=1;//fabs(r[file_head.l0/4][Njack-1]);
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){            
              //r[i][j]/=norm;//;
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0);   mysprintf(label,NAMESIZE,"HA");
   if (sym[index]==1);   mysprintf(label,NAMESIZE,"HV");
   fit=fit_plateaux_G(option, kinematic_2pt_G ,  name,label,mt,r,  Njack,plateaux_masses,outfile);
  // write_jack_bin(Njack,fit,file_jack.f_PS);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
       
   fflush(outfile);
     
    
    return fit;    
    
}

double   *meffH(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
{
   int line=kinematic_2pt_G.i;
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);

   double **r,*m,**mt,*fit;
   int i,j,yn,ii;
   char label[NAMESIZE];
   double kp,xG; 
   
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);
      
    double  E_gT=kinematic_2pt_G.E_gT;
    double  E_gT0=0.5;  
      
   

   fprintf(outfile,"##m_eff_H(t) from %s  propagators:1) mu %.5f r %d theta0 %.5f thetat %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
           kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );
   kp=mass_jack_fit_k2k1[Njack-1]*kinematic_2pt_G.E_g- kinematic_2pt_G.kp;
   xG=2*kp/(mass_rest[Njack-1]*mass_rest[Njack-1]);
   fprintf(outfile,"## E_g=%g      E_gT=%g     E=%g     xG=%g  index=%d\n",kinematic_2pt_G.E_g,kinematic_2pt_G.E_gT,mass_jack_fit_k2k1[Njack-1],xG ,line);
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){            
              ii=sym[index];
              r[i][j]=log(conf_jack[j][index][i][ii]/conf_jack[j][index][i+1][ii]);
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0);   mysprintf(label,NAMESIZE,"meff_HA");
   if (sym[index]==1);   mysprintf(label,NAMESIZE,"m_eff_HV");
   fit=fit_plateaux_G(option, kinematic_2pt_G ,  name,label,mt,r,  Njack,plateaux_masses,outfile);
  // write_jack_bin(Njack,fit,file_jack.f_PS);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
       
   fflush(outfile);
     
    
    return fit;    
    
}



double   *H_over_H0_vir(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
{
   int line=kinematic_2pt_G.i;
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);

   double **r,*m,**mt,*fit;
   int i,j,yn,ii;
   char label[NAMESIZE];
   double kp,xG; 
   
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);
      
    double  E_gT=kinematic_2pt_G.E_gT;
    double  E_gT0=0.5;  
      
   

   fprintf(outfile,"##R_A_or_V(t) from %s  propagators:1) mu %.5f r %d theta0 %.5f thetat %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
           kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );
   kp=mass_jack_fit_k2k1[Njack-1]*kinematic_2pt_G.E_g- kinematic_2pt_G.kp;
   xG=2*kp/(mass_rest[Njack-1]*mass_rest[Njack-1]);
   fprintf(outfile,"## E_g=%g      E_gT=%g     E=%g     xG=%g  index=%d\n",kinematic_2pt_G.E_g,kinematic_2pt_G.E_gT,mass_jack_fit_k2k1[Njack-1],xG ,line);
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){            
              ii=sym[index];
              r[i][j]=conf_jack[j][index][i][ii]/conf_jack[j][index+2][i][ii];    // index+2 is the correlator at p =0
              r[i][j]=r[i][j]/exp(-(file_head.l0/2-i)*kinematic_2pt_G.E_g );
              //r[i][j]=r[i][j]*E_gT/E_gT0;
              //r[i][j]=log(conf_jack[j][index][i][ii]/conf_jack[j][index][i+1][ii]);
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0);   mysprintf(label,NAMESIZE,"H_H0_A_mu");
   if (sym[index]==1);   mysprintf(label,NAMESIZE,"H_H0_V_mu");
   fit=fit_plateaux_G(option, kinematic_2pt_G ,  name,"H_H0_A",mt,r,  Njack,plateaux_masses,outfile);
  // write_jack_bin(Njack,fit,file_jack.f_PS);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
       
   fflush(outfile);
     
    
    return fit;    
    
}

double   *H_minus_H0_HA_vir(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest,  double *Zf_PS_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
{
   int line=kinematic_2pt_G.i;
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn,ii,ia,jr;
   char label[NAMESIZE];
   double kp,xG,tmp; 
   
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);
      
    double  E_gT=kinematic_2pt_G.E_gT;
    double  E_gT0=0.5;  
    double  E_g=kinematic_2pt_G.E_g;
    double  E_g0=0;
    double L[4],tmp1;
    //Assuming L=T/2;
    L[0]=file_head.l0; L[1]=file_head.l0/2.; L[2]=file_head.l0/2.; L[3]=file_head.l0/2.;   

   fprintf(outfile,"##R_A_or_V(t) from %s  propagators:1) mu %.5f r %d theta0 %.5f thetat %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
           kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );
   kp=mass_jack_fit_k2k1[Njack-1]*kinematic_2pt_G.E_g- kinematic_2pt_G.kp;
   xG=2*kp/(mass_rest[Njack-1]*mass_rest[Njack-1]);
   double kin=kinematic_2pt_G.E_g*kinematic_2pt_G.eps2_curl_p[1]-mass_jack_fit_k2k1[Njack-1]*kinematic_2pt_G.eps2_curl_k[1];
   fprintf(outfile,"## E_g=%g      E_gT=%g     E=%g     xG=%g  index=%d   kin=%g\n",kinematic_2pt_G.E_g,kinematic_2pt_G.E_gT,mass_jack_fit_k2k1[Njack-1],xG ,line,kin);
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){  
              r[i][j]=0;
          

              
              ii=sym[index];
              ia=sym[index+1];
              r[i][j]=conf_jack[j][index][i][ii];
              tmp=conf_jack[j][index+2][i][ii]/(   exp(-(L[0]/2.-i)*E_g0 ));    // index+2 is the correlator at p =0
              r[i][j]=r[i][j]/(   exp(-(L[0]/2.-i)*E_g ));
              r[i][j]=r[i][j]-tmp;

              r[i][j]/=conf_jack[j][index+1][i][ia];
              r[i][j]=r[i][j]*mass_rest[j]*Zf_PS_jack_fit[j];
              
              r[i][j]/=(kinematic_2pt_G.E_g*kinematic_2pt_G.eps2_curl_p[1]-mass_jack_fit_k2k1[j]*kinematic_2pt_G.eps2_curl_k[1]);
              r[i][j]*=kinematic_2pt_G.eps1[1];
              //r[i][j]=log(conf_jack[j][index][i][ii]/conf_jack[j][index][i+1][ii]);
        
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0);   mysprintf(label,NAMESIZE,"H_H0_A_mu");
   if (sym[index]==1);   mysprintf(label,NAMESIZE,"HmH0_V_mu");
   fit=fit_plateaux_G(option, kinematic_2pt_G ,  name,"HmH0_V_HA",mt,r,  Njack,plateaux_masses,outfile);
  // write_jack_bin(Njack,fit,file_jack.f_PS);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
       
   fflush(outfile);
     
    
    return fit;    
    
}
