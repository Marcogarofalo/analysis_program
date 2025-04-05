#define m_eff_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cmath>

#include "linear_fit.hpp"
#include "mutils.hpp"
#include "resampling.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "global.hpp"
#include "eigensystem.hpp"



double M_eff(  int t, double **in){
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
/*  
    if(t==0)
        mass=acosh((in[L0-1][0]+in[t+1][0])/(2.*in[t][0]));
    else
        mass=acosh((in[t-1][0]+in[t+1][0])/(2.*in[t][0]));
  */  
    return mass;
}
double M_eff_in_inp(  int t, double in, double inp){
    double mass;
    
    double ct[1],ctp[1],res,tmp_mass, u,d ;
    int i,L0;
    
    L0=file_head.l0;
    ct[0]=in;
    ctp[0]=inp;
    
    
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
     *       mass=acosh((in[L0-1][0]+in[t+1][0])/(2.*in[t][0]));
     *   else
     *       mass=acosh((in[t-1][0]+in[t+1][0])/(2.*in[t][0]));
     */
    return mass;
}

double matrix_element_sl_ss(int t, double **in_sl,double **in_ss,double mass){
      double me_s,me_l;
      
      me_s=in_ss[t][0]/( exp(-mass*t)+exp(-(file_head.l0-t)*mass) );
      me_s*=2*mass;
      me_s=sqrt(me_s);

      me_l=in_sl[t][0]/( exp(-mass*t)+exp(-(file_head.l0-t)*mass) );
      me_l*=2*mass;
      me_l/=me_s;

      return  me_l;
}
double matrix_element_ll(int t, double **in_ll,double mass){
      double me_l;
  /*   
//extra
    double ct[1], ctp[1], res , u , d, tmp_mass;
    int i,L0=file_head.l0;
    ct[0]=in_ll[t][0];
    ctp[0]=in_ll[t+1][0];


    mass=log(ct[0]/ctp[0]);

    res=1;
    i=t;
    while(res>1e-19){
             u=1.+exp(-mass*(L0-2*i-2));
             d=1.+exp(-mass*(L0-2*i));
             tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
             res=tmp_mass - mass;
             mass=tmp_mass;
    }

//fine
*/


      me_l=in_ll[t][0]/( exp(-mass*t)+exp(-(file_head.l0-t)*mass) );
      me_l*=2*mass;
      me_l=sqrt(me_l);
   
      return  me_l;
}

double f_PS_from_Amu_ll(int t, double ***in,double mass,double oPp ,struct kinematic kinematic_2pt){
      double me_oPp,f_PS=0;
      double Pi=3.141592653589793;
      double vol=0,pi;
      double L[4];
      int i;
      //in[0]= oPPo
      //in[1]=oA0Po
      //in[2]=oA1Po
      L[0]=file_head.l0; L[1]=file_head.l1; L[2]=file_head.l2; L[3]=file_head.l3;
 
 //extra     
/*
    double ct[1], ctp[1], res , u , d, tmp_mass;
    int L0=file_head.l0;
    ct[0]=in[0][t][0];
    ctp[0]=in[0][t+1][0];


    mass=log(ct[0]/ctp[0]);

    res=1;
    i=t;
    while(res>1e-19){
             u=1.+exp(-mass*(L0-2*i-2));
             d=1.+exp(-mass*(L0-2*i));
             tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
             res=tmp_mass - mass;
             mass=tmp_mass;
    }


     //mass=acosh((in[0][t-1][0]+in[0][t+1][0])/(2.*in[0][t][0]));
      
     oPp=matrix_element_ll(t,in[0],mass);
*/ //fine extra    
     
     f_PS=(oPp/(mass))    *  ((in[0+1][t][0])/in[0][t][0]);
     f_PS/=sinh(mass); 
         
      vol++;
     
      return  f_PS;
}






double *constant_fit_m_eff(int M, double in){
    double *r;
    
    r=(double*) malloc(sizeof(double)*M);
    r[0]=1.;
    
    return r;
}
double M_eff_t( int t,int L0, double ***in, int num_corr){
    double mass;
 

    double ct[1],ctp[1],res,tmp_mass, u,d ;
    int i;
    ct[0]=in[num_corr][t][0];
    ctp[0]=in[num_corr][t+1][0];


    mass=log(ct[0]/ctp[0]);

    res=1;
    i=t;
    while(res>1e-12){
	     u=1.+exp(-mass*(L0-2*i-2));
	     d=1.+exp(-mass*(L0-2*i));
	     tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
	     res=fabs(tmp_mass - mass);
	     mass=tmp_mass;
    }
     return mass;

}


double *create_jmeff_from_jcorr(int Njack,int L0, int tmin ,int tmax, double ****in, int num_corr ){
   double **r,***y,**mt;
   double *fit,**tmp,*x;
   int i,j;
    
 
   fit=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*L0);
   r=(double**) malloc(sizeof(double*)*L0);
   x=(double*) malloc(sizeof(double)*(tmax-tmin));
   for(i=0;i<L0;i++)
        r[i]=(double*) malloc(sizeof(double)*Njack);
    
   y=(double***) malloc(sizeof(double**)*Njack);
   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(tmax-tmin));
        for (i=tmin;i<tmax;i++){
	       y[j][i-tmin]=(double*) malloc(sizeof(double)*2);
        }
   }

   for(i=1;i<L0/2;i++){    
        for (j=0;j<Njack;j++){
            r[i][j]=M_eff_t(i,L0,in[j],num_corr);
            if (i>=tmin && i<tmax){
                y[j][i-tmin][0]=r[i][j];
            }

        }
        mt[i]=mean_and_error_jack(Njack, r[i]);
        
        if (i>=tmin && i<tmax){
            for (j=0;j<Njack;j++){
                y[j][i-tmin][1]=mt[i][1];
            }
        }
    free(mt[i]);
    }
    free(mt);

    for (j=0;j<Njack;j++){
        tmp=linear_fit( tmax-tmin, x, y[j],  1,constant_fit_m_eff );
        fit[j]=tmp[0][0];
        free(tmp[0]);free(tmp);
    }
    

   for(i=0;i<L0;i++)
      free(r[i]);
   free(r);

   for (j=0;j<Njack;j++){
        for (i=tmin;i<tmax;i++){
		free(y[j][i-tmin]);
        }
        free(y[j]);
   }
   free(y);
   free(x);
   return fit;
}

double ratio_i( int t, double ***in){
    double mass;
    int i;    
    mass=0;
 
      mass=( (in[4][t][1]+in[8][t][1]+in[12][t][1])/(3.*in[0][t][0]) );// *sqrt(in[0][t][0]*in[1][t][0]/(in[2][t][0]*in[3][t][0]));
    return mass;
}
double ratio_0( int t, double ***in){
    double mass;
     
    
    mass=in[0][t][0]*in[1][t][0]/(in[2][t][0]*in[3][t][0]);
   // mass=log(in[0][t][0]/in[0][t+1][0]);
   //mass=( (in[4][t][1]+in[8][t][1]+in[12][t][1])/(in[0][t][0]) );
    return mass;
}

void replace_plateau(struct kinematic kinematic_2pt ,const char* description,int tmin,int tmax,int sep){
    int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1)+1;
    
    char instruction[NAMESIZE];
    if (strcmp(description,"M_{PS}^{ll}")==0 || strcmp(description,"Zf_{PS}^{ll}")==0 || strcmp(description,"oPp")==0 ){
        mysprintf(instruction,NAMESIZE,"sed -i \"%ds/.*/%d  %d  %d/\"  %s",line,tmin,tmax,sep,kinematic_2pt.plateau_m_ll);
    }
    else if(strcmp(description,"M_{PS}^{GEVP}")==0){
        line=kinematic_2pt.ik2+1;
        mysprintf(instruction,NAMESIZE,"sed -i \"%ds/.*/%d  %d  %d/\"  %s",line,tmin,tmax,sep,kinematic_2pt.plateau_m_GEVP);
    }
    else if(strcmp(description,"f_{PS}^{ls-ss}")==0   || strcmp(description,"f_{PS}")==0   ){
        mysprintf(instruction,NAMESIZE,"sed -i \"%ds/.*/%d  %d  %d/\"  %s",line,tmin,tmax,sep,kinematic_2pt.plateau_f);
    }

    int fi=system(instruction);
    
}

double *fit_plateaux(char **option,struct kinematic kinematic_2pt , char* name,const char *description,double **mt,double **r, int Njack,FILE *plateaux_masses,FILE *outfile)
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
        fscanf(plateaux_masses,"%d  %d  %d\n",&tmin,&tmax,&sep);
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
            replace_plateau(  kinematic_2pt ,  description,tmin,tmax,sep);

         }
          free(m);  free(chi2);

      /*      
        yn=1;
        while(yn>0){
            yn=0;
            m=try_linear_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2);
            if(chi2[0]>2){
                    printf("#%s from %s 1) mu %g r %d theta %g 2) mu %g r %d theta %g \n",
                        description,name,
                    kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
                    kinematic_2pt.k1,kinematic_2pt.r1,kinematic_2pt.mom1);   
                    printf("The chi2 in the range [%d,%d] is %g > 1 try an other interval\n",tmin,tmax,chi2[0]);
                    yn=plotting_fit(file_head.l0, mt , tmin,tmax,m,chi2);
                    if (yn>0) plotting( file_head.l0, mt , &tmin,&tmax, &sep);

            }
            free(m);free(chi2);
       }*/
   }
   m=try_linear_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2);              
   fit=give_jack_linear_fit( tmin,  tmax,sep , mt, r, Njack );    
   
   
   //mysprintf(jack_name,NAMESIZE,"jackknife/%s_jack.txt",description);


   //fprintf(outfile,"\n\n #%s fit in [%d,%d]\n  %.15g    %.15g\n\n\n",description,tmin,tmax,m[0],m[1]);
   fprintf(outfile,"\n\n #%s fit in [%d,%d] chi2=%.5f\n  %.15g    %.15g    %d   %d\n\n\n",description,tmin,tmax,chi2[0],m[0],m[1],tmin,tmax);
   printf("#%s (mu_h=%.4f, mu_l=%.4f) fit in [%d,%d]:  %.15g    %.15g\n",description,kinematic_2pt.k2,kinematic_2pt.k1 ,tmin,tmax,m[0],m[1]);

   if ( strcmp(option[5],"pdf")==0 ){
           plotting_fit_pdf(option,description,file_head.l1, file_head.l0,file_head.beta,file_head.ksea,file_head.musea, mt , tmin,tmax,m,name, kinematic_2pt );
   }
    
   free(m);free(chi2);
   return fit;
}

void replace_plateau_G(struct kinematic_G kinematic_2pt_G , char* name,int tmin,int tmax,int sep){
    int line=kinematic_2pt_G.i+1;
    char instruction[NAMESIZE];
    
    if (strcmp(name,"oVmuGPo")==0 || strcmp(name,"HmH0_V")==0 || strcmp(name,"HmH0_HA_V")==0){
        mysprintf(instruction,NAMESIZE,"sed -i \"%ds/.*/%d  %d  %d/\"  %s",line,tmin,tmax,sep,kinematic_2pt_G.plateau_RV);
        printf("replacing file %s",kinematic_2pt_G.plateau_RV);
    }
    else if(strcmp(name,"oAmuGPo")==0){
        mysprintf(instruction,NAMESIZE,"sed -i \"%ds/.*/%d  %d  %d/\"  %s",line,tmin,tmax,sep,kinematic_2pt_G.plateau_RA);
        printf("replacing file %s",kinematic_2pt_G.plateau_RA);
    }
    else if(strcmp(name,"HmH0_HA_V")==0){
        mysprintf(instruction,NAMESIZE,"sed -i \"%ds/.*/%d  %d  %d/\"  %s",line,tmin,tmax,sep,kinematic_2pt_G.plateau_H_H0_A);
        printf("replacing file %s",kinematic_2pt_G.plateau_H_H0_A);
    }
    else {
        printf("\n!!!!!!!!  NO REPLACEMENT OF PLATEAU !!!!!!!!\n");
    }
    system(instruction);
    
}


double *fit_plateaux_G(char **option,struct kinematic_G kinematic_2pt_G , char* name,const char *description,double **mt,double **r, int Njack,FILE *plateaux_masses,FILE *outfile)
{
   int yn;
   int tmin=15, tmax=20, sep=1;
   double *m,*fit, *chi2;
   char jack_name[NAMESIZE];
   yn=1;
   if ( strcmp(option[1],"see")==0 ){
            while(yn>0){
                printf("#%s from %s %d-%d 1) mu %g r %d theta0 %g thetat %g 2) mu %g r %d theta %g \n",
                        description,name,kinematic_2pt_G.ikt,kinematic_2pt_G.iks,
                    kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
                    kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );  
                plotting( file_head.l0, mt , &tmin,&tmax, &sep);
                m=try_linear_fit(option, tmin, tmax, sep , mt, r, Njack ,&chi2);
                yn=plotting_fit(file_head.l0, mt , tmin,tmax,m,chi2);
                free(chi2);
                free(m);
            }
   }
   if ( strcmp(option[1],"read_plateaux")==0 ){
        fscanf(plateaux_masses,"%d  %d  %d\n",&tmin,&tmax,&sep);
        yn=1;
        m=try_linear_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2);
        yn=1;
        if (kinematic_2pt_G.ikt!=1 && kinematic_2pt_G.iks!=1)
        if (kinematic_2pt_G.ikt!=2 && kinematic_2pt_G.iks!=2)
        if(chi2[0]>3 &&  strcmp(name,"oAmuGPo")!=0  &&  strcmp(name,"oVmuGPo")!=0 &&  strcmp(name,"HmH0_V")!=0  ){
            free(m);free(chi2);
             while(yn>0){
                   m=try_linear_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2);
                   printf("#%s from %s %d-%d 1) mu %g r %d theta0 %g thetat %g 2) mu %g r %d theta %g line_plateaux=%d\n",
                        description,name,kinematic_2pt_G.ikt,kinematic_2pt_G.iks,
                    kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
                    kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] ,kinematic_2pt_G.i+1);  
                    printf("The chi2 in the range [%d,%d] is %g > 1 try an other interval\n",tmin,tmax,chi2[0]);
                    yn=plotting_fit(file_head.l0, mt , tmin,tmax,m,chi2);
                    if (yn>0){ 
                        plotting( file_head.l0, mt , &tmin,&tmax, &sep);
                        free(m);free(chi2);
                    }
            }
            replace_plateau_G(  kinematic_2pt_G ,  name,tmin,tmax,sep);
        }
        
        /*
        while(yn>0){
            yn=0;
                m=try_linear_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2);              
                
                if(chi2[0]>2){
                    printf("#%s from %s 1) mu %g r %d theta0 %g thetat %g 2) mu %g r %d theta %g\n",
                            description,name,
                    kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
                    kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );  
                    printf("The chi2 in the range [%d,%d] is %g > 1 try an other interval\n",tmin,tmax,chi2[0]);
                    yn=plotting_fit(file_head.l0, mt , tmin,tmax,m,chi2);
                    if (yn>0) plotting( file_head.l0, mt , &tmin,&tmax, &sep);

                }
                free(m);free(chi2);
        }*/
        free(m);free(chi2);
   }
   m=try_linear_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2);  
   fit=give_jack_linear_fit( tmin,  tmax,sep , mt, r, Njack );   
   mysprintf(jack_name,NAMESIZE,"jackknife/%s_jack.txt",description);


   fprintf(outfile,"\n\n #%s fit in [%d,%d] chi2=%.5f\n  %.15g    %.15g\n\n\n",description,tmin,tmax,chi2[0],m[0],m[1]);
   printf("#%s (mu_h=%.4f, mu_l=%.4f) fit in [%d,%d]:  %.15g    %.15g\n",description,kinematic_2pt_G.kt,kinematic_2pt_G.ks ,tmin,tmax,m[0],m[1]);

   if ( strcmp(option[5],"pdf")==0 ){
           plotting_fit_pdf_G(option,description, mt , tmin,tmax,m,name, kinematic_2pt_G );
   }
    
   free(m);free(chi2);
   return fit;
}




double   *compute_effective_mass(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,FILE **plateaux_masses,FILE *outfile,  int index , const char *description ){
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

   fit=fit_plateaux(option, kinematic_2pt ,  name,description/*"M_{PS}^{ll}"*/,mt,r,  Njack,*plateaux_masses,outfile);
   if( strcmp(description,"M_{PS}^{ll}")==0)
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


double   *compute_effective_mass_GEVP(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack, int var ,FILE *plateaux_masses,FILE *outfile ){
   int line=kinematic_2pt.ik2;
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);

   double **r,*m,**mt,*fit;
   int i,j,yn,t;
   int tmin=11, tmax=15, sep=1;
   int t0=3;
   double ****M,****vec,****lambda,****lambda0;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);
  
   M=(double****) malloc(sizeof(double***)*file_head.l0+1);
   vec=(double****) malloc(sizeof(double***)*file_head.l0);
   lambda=(double****) malloc(sizeof(double***)*file_head.l0);
   for (i=0;i<file_head.l0;i++){
        M[i]=(double***) malloc(sizeof(double**)*Njack);
        vec[i]=(double***) malloc(sizeof(double**)*Njack);
        lambda[i]=(double***) malloc(sizeof(double**)*Njack);
        for (j=0;j<Njack;j++){
            M[i][j]=(double**) malloc(sizeof(double*)*var*var);
            vec[i][j]=(double**) malloc(sizeof(double*)*var*var);
            lambda[i][j]=(double**) malloc(sizeof(double*)*var*var);
            for(t=0;t<var*var;t++){
                M[i][j][t]=(double*) malloc(2*sizeof(double));
                vec[i][j][t]=(double*) malloc(2*sizeof(double));
                lambda[i][j][t]=(double*) malloc(2*sizeof(double));
            }
        }
   }
    
   lambda0=(double****) malloc(sizeof(double***)*Njack);
   for (i=0;i<Njack;i++){
        lambda0[i]=(double***) malloc(sizeof(double**)*var);
        for (j=0;j<var;j++){
            lambda0[i][j]=(double**) malloc(sizeof(double*)*file_head.l0+1);
            for(t=0;t<file_head.l0;t++){
                lambda0[i][j][t]=(double*) malloc(2*sizeof(double));
    	    }
        }
   }


   for(t=1;t<=file_head.l0/2;t++){    
           for (j=0;j<Njack;j++){
                for (i=0;i<var*var;i++){
                     M[t][j][i][0]=conf_jack[j][i][t][0];
                     M[t][j][i][1]=conf_jack[j][i][t][1];
                }
                M[t][j][1][0]=(M[t][j][1][0]+M[t][j][2][0])/2.;
                M[t][j][2][0]=M[t][j][1][0];
                
           }
   }
      
   for(t=1;t<=file_head.l0/2;t++){
           for (j=0;j<Njack;j++){
  	        generalysed_Eigenproblem(M[t][j],M[t0][j],var,&lambda[t][j],&vec[t][j]);     
                for (i=0;i<var;i++){
                     lambda0[j][i][t][0]=lambda[t][j][i][0];// lambda[t][jack][ index of the eigenstate ][ reim]
                     lambda0[j][i][t][1]=lambda[t][j][i][1];// lambda0[jack] [indec eigenstate] [t ] [reim] 
                }
                if (lambda0[j][0][t][0]<lambda0[j][1][t][0] ) printf("ERROR t=%d   0=%g  1=%g\n",t,lambda0[j][1][t][0],lambda0[j][0][t][0]);
           }
   }
  /* to print th GEVP correlator: change the loop in t in this function from 0 to file_head.l0
     FILE *petros=fopen("M_PS_GEVP_A30.32","w");
for (j=0;j<Njack;j++)
for(t=0;t<file_head.l0;t++)
    fprintf(petros,"%d   %.15f  %.15f\n",t,lambda0[j][0][t][0],lambda0[j][0][t][1]);
fflush(petros);
*/
   fprintf(outfile,"#m_eff(t) of %s  propagators:1) mu %g r %d theta %g 2) mu %g r %d theta %g\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(t=1;t<file_head.l0/2;t++){ 
           for (j=0;j<Njack;j++){
                //r[t][j]=M_eff(t, lambda0[j][0]);
                
                if((t-t0)>=0){
                    //r[t][j]=M_eff_in_inp(t, lambda0[j][0][t][0], lambda0[j][0][t+1][0]);
                    r[t][j]=M_eff(t, lambda0[j][0]);
                }
                else{ 
                    //r[t][j]=M_eff_in_inp(t, lambda0[j][1][t][0], lambda0[j][1][t+1][0]);
                    r[t][j]=M_eff(t, lambda0[j][1]);
                }
           }
           if( strcmp(option[4],"jack")==0)
               mt[t]=mean_and_error_jack(Njack, r[t]);
           if( strcmp(option[4],"boot")==0)
               mt[t]=mean_and_error_boot(Njack, r[t]);
           fprintf(outfile,"%d   %.15e    %.15e\n",t,mt[t][0],mt[t][1]);
   }

   fit=fit_plateaux(option, kinematic_2pt ,  name,"M_{PS}^{GEVP}",mt,r,  Njack,plateaux_masses,outfile);
   write_jack_bin(Njack,fit,file_jack.M_PS_GEVP);

       /////free memory
   for(i=1;i<file_head.l0/2;i++)
       free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
       free(r[i]);
   free(r);
   free_jack( Njack, var , file_head.l0, lambda0);
   free_jack( file_head.l0, Njack, var*var , lambda);
   free_jack( file_head.l0, Njack, var*var , vec);
   free_jack( file_head.l0, Njack, var*var , M);

       fflush(outfile);
    
    return fit;    
    
}
double   *compute_f_PS_ls_ss(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,int Njack ,FILE *plateaux_masses,FILE *outfile )
{
   int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);

   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);


   fprintf(outfile,"##decay constant f_PS(t) from %s  propagators:1) mu %g r %d theta %g 2) mu %g r %d theta %g\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
              r[i][j]=matrix_element_sl_ss(i,conf_jack[j][1], conf_jack[j][3], mass_jack_fit_k2k1[j]);
              r[i][j]/=  (mass_jack_fit_k2k1[j]*sinh(mass_jack_fit_k2k1[j])  );
              r[i][j]*= (kinematic_2pt.k2+kinematic_2pt.k1);
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           else if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           else error(0==0,1,"compute_f_PS_ls_ss","argv[4] is %s while the only options supported are jack or boot",option[4]);
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   fit=fit_plateaux(option, kinematic_2pt ,  name,"f_{PS}^{ls-ss}",mt,r,  Njack,plateaux_masses,outfile);
   write_jack_bin(Njack,fit,file_jack.f_PS_ls_ss);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
       
   fflush(outfile);
     
    
    return fit;    
    
}

double   *compute_f_PS_GEVP(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,int Njack, int var ,FILE *plateaux_masses,FILE *outfile ){

   double **r,*m,**mt,*fit;
   int i,j,yn,t;
   int tmin=11, tmax=15, sep=1;
   int t0=2;
   double ****M,****vec,****lambda,****lambda0;
   double norm;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);
  
   M=(double****) malloc(sizeof(double***)*file_head.l0);
   vec=(double****) malloc(sizeof(double***)*file_head.l0);
   lambda=(double****) malloc(sizeof(double***)*file_head.l0);
   for (i=0;i<file_head.l0;i++){
        M[i]=(double***) malloc(sizeof(double**)*Njack);
        vec[i]=(double***) malloc(sizeof(double**)*Njack);
        lambda[i]=(double***) malloc(sizeof(double**)*Njack);
        for (j=0;j<Njack;j++){
            M[i][j]=(double**) malloc(sizeof(double*)*var*var);
            vec[i][j]=(double**) malloc(sizeof(double*)*var*var);
            lambda[i][j]=(double**) malloc(sizeof(double*)*var*var);
            for(t=0;t<var*var;t++){
                M[i][j][t]=(double*) malloc(2*sizeof(double));
                vec[i][j][t]=(double*) malloc(2*sizeof(double));
                lambda[i][j][t]=(double*) malloc(2*sizeof(double));
            }
        }
   }
    
   lambda0=(double****) malloc(sizeof(double***)*Njack);
   for (i=0;i<Njack;i++){
        lambda0[i]=(double***) malloc(sizeof(double**)*var);
        for (j=0;j<var;j++){
            lambda0[i][j]=(double**) malloc(sizeof(double*)*file_head.l0);
            for(t=0;t<file_head.l0;t++){
                lambda0[i][j][t]=(double*) malloc(2*sizeof(double));
    	    }
        }
   }


   for(t=1;t<file_head.l0/2;t++){    
           for (j=0;j<Njack;j++){
                for (i=0;i<var*var;i++){
                     M[t][j][i][0]=conf_jack[j][i][t][0];
                     M[t][j][i][1]=conf_jack[j][i][t][1];
                }
                M[t][j][1][0]=(M[t][j][1][0]+M[t][j][2][0])/2.;
                M[t][j][2][0]=M[t][j][1][0];
           }
   }
      
   for(t=1;t<file_head.l0/2;t++){
           for (j=0;j<Njack;j++){
  	        generalysed_Eigenproblem(M[t][j],M[t0][j],var,&lambda[t][j],&vec[t][j]);     
            lambda0[j][0][t][0]=0;
            norm=0;
            for (i=0;i<var;i++){
                   lambda0[j][0][t][0]+=vec[t][j][i+0*var][0]*M[t][j][0+i*var][0];
                   lambda0[j][0][t][1]+=vec[t][j][i+0*var][1]*M[t][j][0+i*var][0];
            }
            for (i=0;i<var;i++){
                norm+=lambda0[j][0][t][0]*vec[t][j][i+0*var][0]+lambda0[j][0][t][1]*vec[t][j][i+0*var][1];
            }
            norm=(  exp(-mass_jack_fit_k2k1[j]*t)+exp(  -mass_jack_fit_k2k1[j]*(file_head.l0-t)  )   )/norm;
            norm=sqrt(norm);
            lambda0[j][0][t][0]=0;
            for (i=0;i<var;i++){
                   //lambda0[j][0][t][0]+=vec[t][j][i+0*var][0]*M[t][j][0+i*var][0]/norm;
                   lambda0[j][0][t][0]+=vec[t][j][i+0*var][0]*M[t][j][0+i*var][0]/norm;
            }
           }
   }
   fprintf(outfile,"#f_PS(t) of %s  propagators:1) mu %g r %d theta %g 2) mu %g r %d theta %g\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(t=1;t<file_head.l0/2;t++){
           for (j=0;j<Njack;j++){
                r[t][j]=lambda0[j][0][t][0];
                r[t][j]/=  (mass_jack_fit_k2k1[j]*sinh(mass_jack_fit_k2k1[j])  );
                r[t][j]*= (kinematic_2pt.k2+kinematic_2pt.k1);
           }
           if( strcmp(option[4],"jack")==0)
               mt[t]=mean_and_error_jack(Njack, r[t]);
           if( strcmp(option[4],"boot")==0)
               mt[t]=mean_and_error_boot(Njack, r[t]);
           fprintf(outfile,"%d   %.15e    %.15e\n",t,mt[t][0],mt[t][1]);
   }

   fit=fit_plateaux(option, kinematic_2pt ,  name,"f_{PS}^{GEVP}",mt,r,  Njack,plateaux_masses,outfile);
   write_jack_bin(Njack,fit,file_jack.M_PS_GEVP);

       /////free memory
   for(i=1;i<file_head.l0/2;i++)
       free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
       free(r[i]);
   free(r);
   free_jack( Njack, var , file_head.l0, lambda0);
   free_jack( file_head.l0, Njack, var*var , lambda);
   free_jack( file_head.l0, Njack, var*var , vec);
   free_jack( file_head.l0, Njack, var*var , M);

       fflush(outfile);
    
    return fit;    
    
}  
  
  
double   *compute_Zf_PS_ll(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,int Njack ,FILE *plateaux_masses,FILE *outfile )
{
    
   int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);

   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);


   fprintf(outfile,"##decay constant Zf_PS(t) from %s  propagators:1) mu %g r %d theta %g 2) mu %g r %d theta %g\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
              r[i][j]=matrix_element_ll(i,conf_jack[j][0],mass_jack_fit_k2k1[j]);
              r[i][j]/=  (mass_jack_fit_k2k1[j]*sinh(mass_jack_fit_k2k1[j])  );
              r[i][j]*= (kinematic_2pt.k2+kinematic_2pt.k1);
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   fit=fit_plateaux(option, kinematic_2pt ,  name,"Zf_{PS}^{ll}",mt,r,  Njack,plateaux_masses,outfile);
   write_jack_bin(Njack,fit,file_jack.Zf_PS);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
       
   fflush(outfile);
     
    
    return fit;    
    
}
double   *compute_oPp_ll(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,int Njack ,FILE *plateaux_masses,FILE *outfile , int id)
{
   int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);

   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);


   fprintf(outfile,"##matrix elemnt oPp(t) from %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   fprintf(outfile,"#mass=%g\n",mass_jack_fit_k2k1[Njack-1]);
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
              r[i][j]=matrix_element_ll(i,conf_jack[j][id],mass_jack_fit_k2k1[j]);

           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   fit=fit_plateaux(option, kinematic_2pt ,  name,"oPp",mt,r,  Njack,plateaux_masses,outfile);
   write_jack_bin(Njack,fit,file_jack.f_PS);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
       
   fflush(outfile);
     
    
    return fit;    
    
}
   
   
double   *compute_oPp_s(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *oPp_l,int Njack ,FILE *plateaux_masses,FILE *outfile , int i_ll, int i_ls)
{
   int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);

   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);


   fprintf(outfile,"##matrix elemnt oPp(t) from %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   //fprintf(outfile,"#mass=%g\n",mass_jack_fit_k2k1[Njack-1]);
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
              //r[i][j]=matrix_element_ll(i,conf_jack[j][0],mass_jack_fit_k2k1[j]);
              r[i][j]=conf_jack[j][i_ls][i][0]/conf_jack[j][i_ll][i][0];
              r[i][j]*=oPp_l[j];

           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   fit=fit_plateaux(option, kinematic_2pt ,  name,"oPp",mt,r,  Njack,plateaux_masses,outfile);
   write_jack_bin(Njack,fit,file_jack.f_PS);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
       
   fflush(outfile);
     
    
    return fit;    
    
}   

double   *compute_f_PS_ll(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1, double *oPp_jack_fit ,int Njack ,FILE *plateaux_masses,FILE *outfile )
{
    int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);

   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);


   fprintf(outfile,"##decay constant f_PS(t) from %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
              r[i][j]=f_PS_from_Amu_ll(i,conf_jack[j],mass_jack_fit_k2k1[j],oPp_jack_fit[j],kinematic_2pt);
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   fit=fit_plateaux(option, kinematic_2pt ,  name,"f_{PS}",mt,r,  Njack,plateaux_masses,outfile);
   write_jack_bin(Njack,fit,file_jack.f_PS);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
       
   fflush(outfile);
     
    
    return fit;    
    
}



double Rmur(int t, double ***in,double mass, double oPp,struct kinematic_G kinematic_2pt_G,int index,int *sym){
      double RA=0;
     
 double Pi=3.141592653589793;
      double vol=0,k[4],E_g,E_gT;
      double L[4];
      int i;
       
      L[0]=file_head.l0; L[1]=file_head.l1; L[2]=file_head.l2; L[3]=file_head.l3;         
          
      E_g=kinematic_2pt_G.E_g;
      E_gT=kinematic_2pt_G.E_gT;
//m_eff
      /*
    double ct[1], ctp[1], res , u , d, tmp_mass;
    int L0=file_head.l0;
    ct[0]=in[0][t][0];
    ctp[0]=in[0][t+1][0];

    mass=log(ct[0]/ctp[0]);

    res=1;
    while(res>1e-12){
             u=1.+exp(-mass*(L0-2*t-2));
             d=1.+exp(-mass*(L0-2*t));
             tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
             res=fabs(tmp_mass - mass);
             mass=tmp_mass;
    }*/
///
      if (E_gT==0) E_gT=0.5;
      if (fabs(kinematic_2pt_G.Momt-kinematic_2pt_G.Mom0)<1e-6) E_gT=0.5;
      //E_gT=0.5;
      i=sym[index];
      RA=in[index][t][i]*4*mass*E_gT/(   exp(-t*mass-(L[0]/2-t)*E_g )* oPp      );
     // RA=in[index][t][i]*4*mass*E_gT/(   exp(-t*mass-(L[0]/2-t)*E_g )      );
      return  RA;
}








double Rmur_from_meff(int t, double ***in,struct kinematic_G kinematic_2pt_G,int index,int *sym){
      double RA=0;
      double Pi=3.141592653589793;
      double vol=0,k[4],E_g,E_gT;
      double L[4];
      int i;
      if (E_gT==0) E_gT=0.5;
       
      double me_l;
     
//m_eff
    double ct[1], ctp[1], res , u , d, tmp_mass,mass;
    int L0=file_head.l0;
    ct[0]=in[0][t][0];
    ctp[0]=in[0][t+1][0];

    mass=log(ct[0]/ctp[0]);

    res=1;
    while(res>1e-12){
             u=1.+exp(-mass*(L0-2*t-2));
             d=1.+exp(-mass*(L0-2*t));
             tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
             res=fabs(tmp_mass - mass);
             mass=tmp_mass;
    }

    //mass=acosh((in[0][t-1][0]+in[0][t+1][0])/(2.*in[0][t][0]));
// oPp

      me_l=in[0][t][0]/( exp(-mass*t)+exp(-(file_head.l0-t)*mass) );
      me_l*=2*mass;
      me_l=sqrt(me_l);
//RA
      
      L[0]=file_head.l0; L[1]=file_head.l1; L[2]=file_head.l2; L[3]=file_head.l3;         
          
      E_g=kinematic_2pt_G.E_g;
      E_gT=kinematic_2pt_G.E_gT;

      i=sym[index];
      RA=in[index][t][i]*4*mass*E_gT/(   exp(-t*mass-(L[0]/2-t)*E_g )* me_l     );
        
      return  RA;
}

double   *compute_Rmur(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,  double *mass_rest,double *oPp_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
{
    gdb_hook();
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


   fprintf(outfile,"##R_A_or_V(t) from %s  propagators:1) mu %.5f r %d theta0 %.5f thetat %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
           kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );
   kp=mass_jack_fit_k2k1[Njack-1]*kinematic_2pt_G.E_g- kinematic_2pt_G.kp;
   xG=2*kp/(mass_rest[Njack-1]*mass_rest[Njack-1]);
   fprintf(outfile,"## E_g=%g      E_gT=%g     E=%g     xG=%g  index=%d\n",kinematic_2pt_G.E_g,kinematic_2pt_G.E_gT,mass_jack_fit_k2k1[Njack-1],xG ,line);
  // printf("##R_A_or_V(t) from %s  propagators:1) mu %.5f r %d theta0 %.5f thetat %.5f 2) mu %.5f r %d theta %.5f\n",name,
   //        kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
    //       kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );

   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){            

              r[i][j]=Rmur(i,conf_jack[j],mass_jack_fit_k2k1[j],oPp_jack_fit[j],kinematic_2pt_G,index,sym);
  //            ii=sym[index];
  //            r[i][j]=conf_jack[j][index][i][ii];    // these two lines gives the correlator without any manipulations
           }
           
    //       printf("%d  %g\n",4*mass_jack_fit_k2k1[Njack-1]*kinematic_2pt_G.E_gT/(   exp(-i*mass_jack_fit_k2k1[Njack-1]-(file_head.l0[0]/2-i)*kinematic_2pt_G.E_g )* oPp_jack_fit[Njack-1]);

           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0)   mysprintf(label,NAMESIZE,"R_A_mu");
   if (sym[index]==1)   mysprintf(label,NAMESIZE,"R_V_mu");
   fit=fit_plateaux_G(option, kinematic_2pt_G ,  name,"R_{mu}",mt,r,  Njack,plateaux_masses,outfile);
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

double   *H_over_H0(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest, double *oPp_PS_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
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
              r[i][j]=r[i][j]*E_gT/E_gT0;
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0)   mysprintf(label,NAMESIZE,"H_H0_A_mu");
   if (sym[index]==1)   mysprintf(label,NAMESIZE,"H_H0_V_mu");
   fit=fit_plateaux_G(option, kinematic_2pt_G ,  name,"H_H0_{mu}",mt,r,  Njack,plateaux_masses,outfile);
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
/*
double   *H_over_H0(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest, double *oPp_PS_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
{
   int line=kinematic_2pt_G.i;
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn,ii,jr;
   char label[NAMESIZE];
   double kp,xG,tmp; 
   
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
               r[i][j]=0;
               for(jr=0;jr<4;jr++){
                    ii=sym[6+jr];
                    tmp=(conf_jack[j][6+jr][i][ii]+conf_jack[j][6+jr][file_head.l0-i][ii])/(conf_jack[j][10+jr][i][ii]  +conf_jack[j][10+jr][file_head.l0-i][ii]);    // index+2 is the correlator at p =0
                    tmp=tmp/exp(-(file_head.l0/2-i)*kinematic_2pt_G.E_g );
                    tmp=tmp*E_gT/E_gT0;
                    r[i][j]+=tmp;
               }
               r[i][j]/=4;
               //r[i][j]-=1;
           }
              
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0);   mysprintf(label,NAMESIZE,"H_H0_A_mu");
   if (sym[index]==1);   mysprintf(label,NAMESIZE,"H_H0_V_mu");
   fit=fit_plateaux_G(option, kinematic_2pt_G ,  name,"H_H0_{mu}",mt,r,  Njack,plateaux_masses,outfile);
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

*/
double   *H_minus_H0(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest, double *oPp_PS_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
{
   int line=kinematic_2pt_G.i;
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn,ii,t;
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
    double L[4];
    double mass;
    L[0]=file_head.l0; L[1]=file_head.l1; L[2]=file_head.l2; L[3]=file_head.l3;  

   fprintf(outfile,"##R_A_or_V(t) from %s  propagators:1) mu %.5f r %d theta0 %.5f thetat %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
           kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );
   kp=mass_jack_fit_k2k1[Njack-1]*kinematic_2pt_G.E_g- kinematic_2pt_G.kp;
   xG=2*kp/(mass_rest[Njack-1]*mass_rest[Njack-1]);
   fprintf(outfile,"## E_g=%g      E_gT=%g     E=%g     xG=%g  index=%d\n",kinematic_2pt_G.E_g,kinematic_2pt_G.E_gT,mass_jack_fit_k2k1[Njack-1],xG ,line);
   for(i=1;i<file_head.l0/2;i++){    
           t=i;
           for (j=0;j<Njack;j++){            
              ii=sym[index];
              
              mass=mass_jack_fit_k2k1[j];
              //m_eff
          /*          double ct[1], ctp[1], res , u , d, tmp_mass;
                    int L0=file_head.l0;
                    ct[0]=conf_jack[j][0][t][0];
                    ctp[0]=conf_jack[j][0][t+1][0];

                    mass=log(ct[0]/ctp[0]);

                    res=1;
                    while(res>1e-12){
                            u=1.+exp(-mass*(L0-2*t-2));
                            d=1.+exp(-mass*(L0-2*t));
                            tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
                            res=fabs(tmp_mass - mass);
                            mass=tmp_mass;
                    }*/
              ///mass
              r[i][j]=conf_jack[j][index][i][ii]*4*mass*E_gT/(   exp(-i*mass-(L[0]/2.-i)*E_g )* oPp_PS_jack_fit[j]);
//              r[i][j]/=(kinematic_2pt_G.E_g*kinematic_2pt_G.eps1_curl_p[1]-mass_jack_fit_k2k1[j]*kinematic_2pt_G.eps1_curl_k[1]);
              tmp=conf_jack[j][index+2][i][ii]*4*mass*E_gT0/(   exp(-i*mass-(L[0]/2.-i)*E_g0 )* oPp_PS_jack_fit[j]);    // index+2 is the correlator at p =0
              r[i][j]=r[i][j]-tmp;
              //r[i][j]/=(   exp(-i*mass_jack_fit_k2k1[j]-(L[0]/2.-i)*E_g )* oPp_PS_jack_fit[j]);
              r[i][j]=-r[i][j]*mass_rest[j];
              r[i][j]/=(kinematic_2pt_G.E_g*kinematic_2pt_G.eps1_curl_p[1]-mass*kinematic_2pt_G.eps1_curl_k[1]);
              
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0)   mysprintf(label,NAMESIZE,"H_H0_A_mu");
   if (sym[index]==1)   mysprintf(label,NAMESIZE,"HmH0_V_mu");
   fit=fit_plateaux_G(option, kinematic_2pt_G ,  name,"HmH0_{mu}",mt,r,  Njack,plateaux_masses,outfile);
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


double   *H_minus_H0_HA(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1, double* mass_rest, double *oPp_PS_jack_fit, double *Zf_PS_jack_fit,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
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
    L[0]=file_head.l0; L[1]=file_head.l1; L[2]=file_head.l2; L[3]=file_head.l3;  

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
              /*
              for(jr=0;jr<4;jr++){
                    ii=sym[index];
                    ia=sym[index+1];
                    tmp1=(conf_jack[j][14+jr][i][ii]-conf_jack[j][14+jr][file_head.l0-i][ii])*E_gT;
                    tmp=(conf_jack[j][18+jr][i][ii]-conf_jack[j][18+jr][file_head.l0-i][ii])*E_gT0;    // index+2 is the correlator at p =0
                    tmp1=tmp1/(   exp(-(L[0]/2.-i)*E_g ));
                    //tmp1=tmp1-tmp;
                    tmp1/=(conf_jack[j][10+jr][i][ia]+conf_jack[j][10+jr][file_head.l0-i][ia])*E_gT0;
                    
                    tmp1=tmp1*mass_rest[j]*Zf_PS_jack_fit[j];
                    if(jr<2){
                        tmp1/=(kinematic_2pt_G.E_g*kinematic_2pt_G.eps1_curl_p[1+jr]-mass_jack_fit_k2k1[j]*kinematic_2pt_G.eps1_curl_k[1+jr]);
                        tmp1*=kinematic_2pt_G.eps1[1+jr];
                        if (j==0 && i ==1)printf("fattore (%d )= %g eps=%g\n",jr,(kinematic_2pt_G.E_g*kinematic_2pt_G.eps1_curl_p[1+jr]-mass_jack_fit_k2k1[j]*kinematic_2pt_G.eps1_curl_k[1+jr]),kinematic_2pt_G.eps1[1+jr]);
                    }
                    else{
                        tmp1/=(kinematic_2pt_G.E_g*kinematic_2pt_G.eps2_curl_p[1+jr-2]-mass_jack_fit_k2k1[j]*kinematic_2pt_G.eps2_curl_k[1+jr-2]);
                        tmp1*=kinematic_2pt_G.eps2[1+jr-2];
                        if (j==0 && i ==1)printf("fattore (%d )= %g eps=%g\n",jr,(kinematic_2pt_G.E_g*kinematic_2pt_G.eps2_curl_p[1+jr-2]-mass_jack_fit_k2k1[j]*kinematic_2pt_G.eps2_curl_k[1+jr-2]),kinematic_2pt_G.eps2[1+jr-2]);

                    }
                    r[i][j]+=tmp1;
                    
              }
              r[i][j]/=4.;
              */
              

              
              ii=sym[index];
              ia=sym[index+1];
              r[i][j]=conf_jack[j][index][i][ii];
              tmp=conf_jack[j][index+2][i][ii]*E_gT0/(   exp(-(L[0]/2.-i)*E_g0 ));    // index+2 is the correlator at p =0
              r[i][j]=r[i][j]*E_gT/(   exp(-(L[0]/2.-i)*E_g ));
              r[i][j]=r[i][j]-tmp;

              r[i][j]/=conf_jack[j][index+1][i][ia]*E_gT0;
              r[i][j]=r[i][j]*mass_rest[j]*Zf_PS_jack_fit[j];
              
              r[i][j]/=(kinematic_2pt_G.E_g*kinematic_2pt_G.eps2_curl_p[1]-mass_jack_fit_k2k1[j]*kinematic_2pt_G.eps2_curl_k[1]);
              r[i][j]*=kinematic_2pt_G.eps1[1];
        
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0)   mysprintf(label,NAMESIZE,"H_H0_A_mu");
   if (sym[index]==1)   mysprintf(label,NAMESIZE,"HmH0_V_mu");
   fit=fit_plateaux_G(option, kinematic_2pt_G ,  name,"HmH0_HA_{mu}",mt,r,  Njack,plateaux_masses,outfile);
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

double   *compute_Rmur_from_meff(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,double *mass_jack_fit_k2k1,double *mass_rest,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
{
     int line=kinematic_2pt_G.i;
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
 
   double **r,*m,**mt,*fit;
   int i,j,yn;
   char label[NAMESIZE];
    double kp,xG; 
  
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);


   fprintf(outfile,"##R_A_or_V(t) from %s  propagators:1) mu %.5f r %d theta0 %.5f thetat %.4f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
           kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );
   kp=mass_jack_fit_k2k1[Njack-1]*kinematic_2pt_G.E_g- kinematic_2pt_G.kp;
   xG=2*kp/(mass_rest[Njack-1]*mass_rest[Njack-1]);
   fprintf(outfile,"## E_g=%g      E_gT=%g     E=%g     xG=%g  index=%d\n",kinematic_2pt_G.E_g,kinematic_2pt_G.E_gT,mass_jack_fit_k2k1[Njack-1],xG ,line);
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){            

              r[i][j]=Rmur_from_meff(i,conf_jack[j],kinematic_2pt_G,index,sym);
          //    r[i][j]=conf_jack[j][0][i][0];
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (sym[index]==0)   mysprintf(label,NAMESIZE,"R_A_mu");
   if (sym[index]==1)   mysprintf(label,NAMESIZE,"R_V_mu");
   fit=fit_plateaux_G(option, kinematic_2pt_G ,  name,"R_{mu}",mt,r,  Njack,plateaux_masses,outfile);
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
double   *compute_CAmur(char **option ,struct kinematic_G kinematic_2pt_G , char* name, double ****conf_jack,int Njack ,FILE *plateaux_masses,FILE *outfile,int index,int *sym )
{
   double **r,*m,**mt,*fit;
   int i,j,yn;
   char label[NAMESIZE];
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);


   fprintf(outfile,"##R_A(t) from %s  propagators:1) mu %g r %d theta0 %g thetat %g 2) mu %g r %d theta %g\n",name,
           kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],
           kinematic_2pt_G.ks,kinematic_2pt_G.rs, kinematic_2pt_G.Moms[3] );
   fprintf(outfile,"## E_g=%g      E_gT=%g     \n",kinematic_2pt_G.E_g,kinematic_2pt_G.E_gT );
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){            

              r[i][j]=conf_jack[j][index][i][0];
          //    r[i][j]=conf_jack[j][0][i][0];
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   if (index==0)   mysprintf(label,NAMESIZE,"R_A_mu");
   if (index==1)   mysprintf(label,NAMESIZE,"R_V_mu");
   fit=fit_plateaux_G(option, kinematic_2pt_G ,  name,"R_{mu}",mt,r,  Njack,plateaux_masses,outfile);
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

double *get_ZA(char **option ,int Njack){
   double *ZA;
   if(file_head.beta==1.726){
  //	 ZA=fake_jack(0.751,0.005,Njack);//M1 a^2 g_*^2  16/2/2019
  //	 ZA=fake_jack(0.735,0.004,Njack);//M2 a^2 g_*^2  16/2/2019
  	 
   //	  ZA=fake_jack(0.747,0.003,Njack);//M1 a^inf g_0^2  16/2/2019
  //	 ZA=fake_jack(0.755,0.004,Njack);//M2 a^inf g_0^2  16/2/2019
     
   //  ZA=fake_jack(0.748,0.005,Njack);//M1 a^2 g_0^2   tilde p  16/2/2019
       ZA=fake_sampling(option[4],0.735,0.004,Njack,123);//M2 a^2 g_+^2    p  29/03/2019
   }
   if(file_head.beta==1.778){
     
       ZA=fake_sampling(option[4],0.755,0.003,Njack,1234);//M2 a^2 g_+^2    p  29/03/2019
   }
   
   return ZA;
}

/*
double   *compute_effective_mass_out_max_twist(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,FILE *plateaux_masses,FILE *outfile,  int index  ){
   int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);

   double mu=kinematic_2pt.k1;
   double mpcac,*ZA;

   ZA=get_ZA(option,Njack);
   
   
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
		         mpcac=conf_jack[j][4][i][1]/(2.*conf_jack[j][5][i][0]);
	             r[i][j]=mpcac;
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           //fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   fit=fit_plateaux(option, kinematic_2pt ,  name,"m_{PCAC}",mt,r,  Njack,plateaux_masses,outfile);
   //free(fit);
   //fit=fake_sampling(option[4],0.0024,0.0005,Njack,666);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   fprintf(outfile,"#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
		         mpcac=fit[j];
	             r[i][j]=M_eff(i,conf_jack[j][index])*pow(1+mpcac*mpcac*ZA[j]*ZA[j]/(mu*mu),-1./4.);
                 //if (j==Njack-1) printf("t=%d   correction=%g     mpcac=%g\n",i,pow(1+mpcac*mpcac*ZA[j]*ZA[j]/(mu*mu),-1./4.),mpcac);
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }
   free(fit);
   fit=fit_plateaux(option, kinematic_2pt ,  name,"M_{PS}^{ll}",mt,r,  Njack,plateaux_masses,outfile);
   write_jack_bin(Njack,fit,file_jack.M_PS);
     
   for(i=1;i<file_head.l0/2;i++)
       free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
       free(r[i]);
   free(r);
   
   free(ZA);
   fflush(outfile);
     
    return fit;    
    
}
*/

double   *compute_mpcac(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,FILE *plateaux_masses,FILE *outfile,  int index  ){
   int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);

   double mu1=kinematic_2pt.k1;
   double mu2=kinematic_2pt.k2;
   double cl,cs,tmp;

  
   
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
		         r[i][j]=-conf_jack[j][4][i][1]/(2.*conf_jack[j][5][i][0]);
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }
   fit=fit_plateaux(option, kinematic_2pt ,  name,"m_{PCAC}",mt,r,  Njack,plateaux_masses,outfile);
   
   write_jack_bin(Njack,fit,file_jack.mpcac);
     
   for(i=1;i<file_head.l0/2;i++)
       free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
       free(r[i]);
   free(r);
   
   fflush(outfile);
     
    return fit;    
    
}
double   *compute_effective_mass_out_max_twist(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,FILE *plateaux_masses,FILE *outfile,  int index  ){
   int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);

   double mu1=kinematic_2pt.k1;
   double mu2=kinematic_2pt.k2;
   double mpcac,*ZA;
   double cl,cs,tmp;

   ZA=get_ZA(option,Njack);
   
   
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
		         mpcac=conf_jack[j][4][i][1]/(2.*conf_jack[j][5][i][0]);
	             r[i][j]=mpcac;
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           //fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }
   fit=fit_plateaux(option, kinematic_2pt ,  name,"m_{PCAC}",mt,r,  Njack,plateaux_masses,outfile);
 free(fit);
   fit=fake_sampling(option[4],0.00024,0.00005,Njack,666);
  if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   fprintf(outfile,"#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
		         mpcac=fit[j];
                 cl=pow(1+mpcac*mpcac*ZA[j]*ZA[j]/(mu1*mu1),-1./2.);
                 cs=pow(1+mpcac*mpcac*ZA[j]*ZA[j]/(mu2*mu2),-1./2.);
                 tmp=mu1/cl + mu2/cs;
	             r[i][j]=M_eff(i,conf_jack[j][index])*sqrt( (mu1+mu2) /tmp );
                 //if (j==Njack-1) printf("t=%d   correction=%g     mpcac=%g\n",i,pow(1+mpcac*mpcac*ZA[j]*ZA[j]/(mu*mu),-1./4.),mpcac);
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }
   free(fit);
   fit=fit_plateaux(option, kinematic_2pt ,  name,"M_{PS}^{ll}",mt,r,  Njack,plateaux_masses,outfile);
   write_jack_bin(Njack,fit,file_jack.M_PS);
     
   for(i=1;i<file_head.l0/2;i++)
       free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
       free(r[i]);
   free(r);
   
   free(ZA);
   fflush(outfile);
     
    return fit;    
    
}


double   *compute_Zf_PS_ll_out_max_twist(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,int Njack ,FILE *plateaux_masses,FILE *outfile )
{
   int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
	   go_to_line(plateaux_masses,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);

   double mu=kinematic_2pt.k1;
   double mu1=kinematic_2pt.k1,  mu2=kinematic_2pt.k2;
   double mpcac,*ZA,cl, cs,mass ,tmp;
   ZA=get_ZA(option,Njack);

   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
		         mpcac=conf_jack[j][4][i][1]/(2.*conf_jack[j][5][i][0]);
	             r[i][j]=mpcac;
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           //fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }
   fit=fit_plateaux(option, kinematic_2pt ,  name,"m_{PCAC}",mt,r,  Njack,plateaux_masses,outfile);
   free(fit);
   fit=fake_sampling(option[4],0.00024,0.00005,Njack,666);
    if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   fprintf(outfile,"##decay constant Zf_PS(t) from %s  propagators:1) mu %g r %d theta %g 2) mu %g r %d theta %g\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
	          mpcac=fit[j];
              //cl=pow(1+mpcac*mpcac*ZA[j]*ZA[j]/(mu*mu),-1./2.);
              cl=pow(1+mpcac*mpcac*ZA[j]*ZA[j]/(mu1*mu1),-1./2.);
              cs=pow(1+mpcac*mpcac*ZA[j]*ZA[j]/(mu2*mu2),-1./2.);
              double Kl=1./cl;
              double Ks=1.;// /cs; 
              tmp=mu1/cl + mu2/cs;
              mass=mass_jack_fit_k2k1[j]/sqrt(  (mu1+mu2)/tmp);
              r[i][j]=matrix_element_ll(i,conf_jack[j][0],mass);
              r[i][j]/=  (mass*sinh(mass)  );
              r[i][j]*= (kinematic_2pt.k2+kinematic_2pt.k1);
              if (kinematic_2pt.ik2==0 && kinematic_2pt.ik1==0)
                    r[i][j]/=cl;
              if (kinematic_2pt.ik2>0 || kinematic_2pt.ik1 >0)
                    r[i][j]*=sqrt(2*Ks*Kl/(Ks+Kl)  );
              if (j==Njack-1) printf("cl=%g   mpcac=%g    mu=%g   ZA=%g\n",cl,mpcac,mu,ZA[j]);
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   fit=fit_plateaux(option, kinematic_2pt ,  name,"Zf_{PS}^{ll}",mt,r,  Njack,plateaux_masses,outfile);
   write_jack_bin(Njack,fit,file_jack.Zf_PS);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
       
   free(ZA);
   fflush(outfile);
     
    
    return fit;    
    
}

double   *compute_f_PS_ls_ss_out_max_twist(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, double *mass_jack_fit_k2k1,int Njack ,FILE *plateaux_masses,FILE *plateaux_mpcac,FILE *outfile )
{
   int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_mpcac,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);

   double mu=kinematic_2pt.k1;
   double mpcac,*ZA;
   ZA=get_ZA(option,Njack);

   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
		         mpcac=conf_jack[j][4][i][1]/(2.*conf_jack[j][5][i][0]);
	             r[i][j]=mpcac;
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           //fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }
   fit=fit_plateaux(option, kinematic_2pt ,  name,"m_{PCAC}",mt,r,  Njack,plateaux_mpcac,outfile);
   free(fit);
   fit=fake_sampling(option[4],0.00024,0.00005,Njack,666);
    if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(plateaux_masses,line);
   
   fprintf(outfile,"##decay constant Zf_PS(t) from %s  propagators:1) mu %g r %d theta %g 2) mu %g r %d theta %g\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
	          mpcac=fit[j];
              r[i][j]=matrix_element_sl_ss(i,conf_jack[j][1], conf_jack[j][3], mass_jack_fit_k2k1[j]);
              r[i][j]/=  (mass_jack_fit_k2k1[j]*sinh(mass_jack_fit_k2k1[j])  );
              r[i][j]*= (kinematic_2pt.k2+kinematic_2pt.k1);
	          //r[i][j]*=(1+mpcac*mpcac*ZA[j]*ZA[j]/(mu*mu));
           }
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }

   fit=fit_plateaux(option, kinematic_2pt ,  name,"f_{PS}^{ls-ss}",mt,r,  Njack,plateaux_masses,outfile);
   write_jack_bin(Njack,fit,file_jack.f_PS_ls_ss);

     
   for(i=1;i<file_head.l0/2;i++)
      free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r);
       
   free(ZA);
   fflush(outfile);
     
    
    return fit;    
    
}
