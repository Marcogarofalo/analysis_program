#define automatic_plateau_R_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <memory>

#include "linear_fit.hpp"
#include "mutils.hpp"
#include "resampling.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "global.hpp"
#include "continuum_reph.hpp"
 /*
double Rmur(int t, double ***in,double mass, double oPp,struct kinematic_G kinematic_2pt_G,int index,int *sym){
      double RA=0;
      double Pi=3.141592653589793;
      double vol=0,k[4],E_g,E_gT;
      double L[4];
      int i;
       
      L[0]=file_head.l0; L[1]=file_head.l1; L[2]=file_head.l2; L[3]=file_head.l3;         
          
      E_g=kinematic_2pt_G.E_g;
      E_gT=kinematic_2pt_G.E_gT;

      i=sym[index];
      RA=in[index][t][i]*4*mass*E_gT/(   exp(-t*mass-(L[0]/2-t)*E_g )* oPp      );
        
      
   
      return  RA;
}
*/
 
struct  p_for_plateau
{
    int min,max;
    double chi2;
}; 
 
void auto_plateau(double **m_eff_exp,double *m_eff_th ,int *tmin, int *tmax){
    int i,t,count,yn,count_r,l,l1;
    double m,diff,err,chi2;
    struct p_for_plateau pf;
    int sep=4;
    
    std::vector< p_for_plateau > p(file_head.l0/2);
    std::vector< p_for_plateau > p1(file_head.l0/2);
    //p=(struct p_for_plateau*) malloc(sizeof(struct p_for_plateau)*file_head.l0/2);
    //p1=(struct p_for_plateau*) malloc(sizeof(struct p_for_plateau)*file_head.l0/2);
    
    
    count=0;yn=0;
    double err_th;
    if (m_eff_th[1]!=m_eff_th[1]) err_th=0;
    else err_th=m_eff_th[1];
    
    for(t=1;t<file_head.l0/2;t++){
        diff=fabs(m_eff_exp[t][0]-m_eff_th[0]);
        err=(m_eff_exp[t][1]+err_th);
       // err=(m_eff_exp[t][1]);
        chi2=diff*diff/(err*err);
//printf("%d   %f    %f   - %f   %f\n",t,m_eff_exp[t][0],m_eff_exp[t][1],m_eff_th[0],m_eff_th[1]);
        if (yn==0){
            if(diff<err){
                p[count].min=t;
                p[count].chi2=0;
                p[count].chi2+=chi2;
                yn=1;
            }
        }
        else if (yn==1){
            if(diff<err){
                p[count].chi2+=chi2;
            }
            else if(diff>=err || t==file_head.l0/2-1){
                p[count].max=t;
                p[count].chi2/=((double) (p[count].max-p[count].min));
                yn=0;
                count++;
            }
            
        }
        
    }

     
    //error(count==0,1,"auto_plateau","no plateau found");
    if(count==0){
     // printf("no plateau found. \n plateau set to  [T/4,T/4]\n");
      *tmin=file_head.l0/4;
      *tmax=file_head.l0/4;
      return ;
      
    }
    
    //printf("possible plateaux\n");
    for(t=0;t<count;t++){
    //  printf("%d   %d   %f\n", p[t].min,p[t].max,p[t].chi2);   
     }
  //   printf("\n");
    count_r=1;
    p1[0].min=p[0].min;
    p1[0].max=p[0].max;
    
    
    for(t=1;t<count;t++){
        if(p[t].min-p1[count_r-1].max<sep){
            p1[count_r-1].max=p[t].max;
        }
        else{
            count_r++;
            p1[count_r-1].min=p[t].min;
            p1[count_r-1].max=p[t].max;
        }
    }
    
   // printf("glue plateaux\n");
    for(t=0;t<count_r;t++){
    //  printf("%d   %d   \n", p1[t].min,p1[t].max);   
     }
  //   printf("\n");
    
    
    pf.min=p1[0].min;
    pf.max=p1[0].max;
    l=pf.max-pf.min;
    for(t=1;t<count_r;t++){
        l1=p1[t].max-p1[t].min;
        if(l1>l){
            pf.min=p1[t].min;
            pf.max=p1[t].max;
            l=l1;
        }
            
    }
    *tmin=pf.min;
    *tmax=pf.max-1;
    
    //printf("final plateau: %d   %d   \n", *tmin,*tmax+1);   
     
   //  printf("\n");
    
    //free(p);
    //free(p1);
} 
 
double *fit_plateaux_G_int(char **option,struct kinematic_G kinematic_2pt_G , char* name,const char *description,double **mt,double **r, int Njack,FILE *plateaux_masses,FILE *outfile,int tmin,int tmax)
{
   int yn;
   int  sep=1;
   double *m,*fit, *chi2;
   char jack_name[500];
   yn=1;
   
   m=try_linear_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2);  
   fit=give_jack_linear_fit( tmin,  tmax,sep , mt, r, Njack );   
   sprintf(jack_name,"jackknife/%s_jack.txt",description);


   fprintf(outfile,"\n\n #%s fit in [%d,%d] chi2=%.5f\n  %.15g    %.15g\n\n\n",description,tmin,tmax,chi2[0],m[0],m[1]);
   printf("#%s (mu_h=%.4f, mu_l=%.4f) fit in [%d,%d]:  %.15g    %.15g\n",description,kinematic_2pt_G.kt,kinematic_2pt_G.ks ,tmin,tmax,m[0],m[1]);

   if ( strcmp(option[5],"pdf")==0 ){
           plotting_fit_pdf_G(option,description, mt , tmin,tmax,m,name, kinematic_2pt_G );
   }
    
   free(m);free(chi2);
   return fit;
}


double   *compute_Rmur_auto_plateau(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,double *mass_rest,double *oPp_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
{
   int line=kinematic_2pt_G.i;
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   double **r,*m,**mt,*fit,*m_th;
   int i,j,yn,ii,t;
   char label[500];
   double kp,xG; 
   
   int tmin, tmax;
   
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);


   fprintf(outfile,"##R_A_or_V(t) from %s  propagators:1) mu %.5f r %d theta0 %.5f thetat %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
           kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );
   kp=mass_jack_fit_k2k1[Njack-1]*kinematic_2pt_G.E_g- kinematic_2pt_G.kp;
   xG=2*kp/(mass_rest[Njack-1]*mass_rest[Njack-1]);
   fprintf(outfile,"## E_g=%g      E_gT=%g     E=%g     xG=%g  index=%d\n",kinematic_2pt_G.E_g,kinematic_2pt_G.E_gT,mass_jack_fit_k2k1[Njack-1],xG ,line);
   
   m=(double*) malloc(sizeof(double)*Njack);
   for (j=0;j<Njack;j++){            
        m[j]=mass_jack_fit_k2k1[j]-kinematic_2pt_G.E_g;
   }
   m_th=mean_and_error_jack(Njack, m);
   
   i=sym[index];
   for(t=1;t<file_head.l0/2;t++){    
           for (j=0;j<Njack;j++){ 
              r[t][j]= log(conf_jack[j][index][t][i]/ conf_jack[j][index][t+1][i]);
           }
           mt[t]=mean_and_error_jack(Njack, r[t]);
   }
   auto_plateau(mt ,m_th ,&tmin, &tmax);
   free(m);
   for(t=1;t<file_head.l0/2;t++)
       free(mt[t]);

   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){            

              r[i][j]=Rmur(i,conf_jack[j],mass_jack_fit_k2k1[j],oPp_jack_fit[j],kinematic_2pt_G,index,sym);
//               ii=sym[index];
//              r[i][j]=conf_jack[j][index][i][ii];    // these two lines gives the correlator without any manipulations
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0)   sprintf(label,"R_A_mu");
   if (sym[index]==1)   sprintf(label,"R_V_mu");
   fit=fit_plateaux_G_int(option, kinematic_2pt_G ,  name,"R_{mu}",mt,r,  Njack,plateaux_masses,outfile,tmin,tmax);
  // write_jack_bin(Njack,fit,file_jack.f_PS);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
   free(m_th);
       
   fflush(outfile);
     
    
    return fit;    
    
}



double   *H_over_H0_autoplateaux(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,double *mass_rest,double *oPp_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
{
   int line=kinematic_2pt_G.i;
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   double **r,*m,**mt,*fit,*m_th;
   int i,j,yn,ii,t;
   char label[500];
   double kp,xG; 
   
   int tmin, tmax;
   
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
   
   m=(double*) malloc(sizeof(double)*Njack);
   for (j=0;j<Njack;j++){            
        m[j]=-kinematic_2pt_G.E_g;
   }
   //m_th=mean_and_error(option[4],Njack, m);
   m_th=(double*) malloc(sizeof(double)*2);
   m_th[0]=-kinematic_2pt_G.E_g;
   m_th[1]=fabs(kinematic_2pt_G.E_g/100.0);

   i=sym[index];
   for(t=1;t<file_head.l0/2;t++){    
           for (j=0;j<Njack;j++){ 
              ii=sym[index];
              r[t][j]=conf_jack[j][index][t][ii]/conf_jack[j][index+2][t][ii];    // index+2 is the correlator at p =0
              r[t+1][j]=conf_jack[j][index][t+1][ii]/conf_jack[j][index+2][t+1][ii];
              r[t][j]= log(r[t][j]/ r[t+1][j]);
           }
           mt[t]=mean_and_error_jack(Njack, r[t]);
   }
   auto_plateau(mt ,m_th ,&tmin, &tmax);

   free(m);
   for(t=1;t<file_head.l0/2;t++)
       free(mt[t]);

   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){            

              ii=sym[index];
              r[i][j]=conf_jack[j][index][i][ii]/conf_jack[j][index+2][i][ii];    // index+2 is the correlator at p =0
              r[i][j]=r[i][j]/exp(-(file_head.l0/2-i)*kinematic_2pt_G.E_g );
              r[i][j]=r[i][j]*E_gT/E_gT0;
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0)   sprintf(label,"R_A_mu");
   if (sym[index]==1)   sprintf(label,"R_V_mu");
   fit=fit_plateaux_G_int(option, kinematic_2pt_G ,  name,"FA_from_H0_auto",mt,r,  Njack,plateaux_masses,outfile,tmin,tmax);
  // write_jack_bin(Njack,fit,file_jack.f_PS);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
   free(m_th);
       
   fflush(outfile);
     
    
    return fit;    
    
}


double   *H_minus_H0_autoplateaux(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,double *mass_rest,double *oPp_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
{
   int line=kinematic_2pt_G.i;
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   double **r,*m,**mt,*fit,*m_th;
   int i,j,yn,ii,t;
   char label[500];
   double kp,xG,tmp; 
   
   int tmin, tmax;
   double  E_gT=kinematic_2pt_G.E_gT;
   double  E_gT0=0.5;  
   double  E_g=kinematic_2pt_G.E_g;
   double  E_g0=0;
   double L[4];
   L[0]=file_head.l0; L[1]=file_head.l1; L[2]=file_head.l2; L[3]=file_head.l3;  

   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);


   fprintf(outfile,"##R_A_or_V(t) from %s  propagators:1) mu %.5f r %d theta0 %.5f thetat %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
           kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );
   kp=mass_jack_fit_k2k1[Njack-1]*kinematic_2pt_G.E_g- kinematic_2pt_G.kp;
   xG=2*kp/(mass_rest[Njack-1]*mass_rest[Njack-1]);
   fprintf(outfile,"## E_g=%g      E_gT=%g     E=%g     xG=%g  index=%d\n",kinematic_2pt_G.E_g,kinematic_2pt_G.E_gT,mass_jack_fit_k2k1[Njack-1],xG ,line);
   
   m=(double*) malloc(sizeof(double)*Njack);
   for (j=0;j<Njack;j++){            
        m[j]=mass_jack_fit_k2k1[j]-kinematic_2pt_G.E_g;
   }
   m_th=mean_and_error_jack(Njack, m);
   
   i=sym[index];
   for(t=1;t<file_head.l0/2;t++){    
           for (j=0;j<Njack;j++){ 
              r[t][j]= log(conf_jack[j][index][t][i]/ conf_jack[j][index][t+1][i]);
           }
           mt[t]=mean_and_error_jack(Njack, r[t]);
   }
   auto_plateau(mt ,m_th ,&tmin, &tmax);
   free(m);
   for(t=1;t<file_head.l0/2;t++)
       free(mt[t]);

   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){            
              ii=sym[index];
              r[i][j]=conf_jack[j][index][i][ii]*4*mass_jack_fit_k2k1[j]*E_gT/(   exp(-i*mass_jack_fit_k2k1[j]-(L[0]/2.-i)*E_g )* oPp_jack_fit[j]);
              tmp=conf_jack[j][index+2][i][ii]*4*mass_jack_fit_k2k1[j]*E_gT0/(   exp(-i*mass_jack_fit_k2k1[j]-(L[0]/2.-i)*E_g0 )* oPp_jack_fit[j]);    // index+2 is the correlator at p =0
              r[i][j]=r[i][j]-tmp;
            }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0);   sprintf(label,"R_A_mu");
   if (sym[index]==1);   sprintf(label,"R_V_mu");
   fit=fit_plateaux_G_int(option, kinematic_2pt_G ,  name,"HmH0_{mu}_auto",mt,r,  Njack,plateaux_masses,outfile,tmin,tmax);
  // write_jack_bin(Njack,fit,file_jack.f_PS);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
   free(m_th);
       
   fflush(outfile);
     
    
    return fit;    
    
}


