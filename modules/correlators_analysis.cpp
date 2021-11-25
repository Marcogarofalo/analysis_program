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
#include "tower.hpp"
#include "read.hpp"
//#include "eigensystem.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <iterator>
#include <random>

void check_correlatro_counter(int i){
    if( corr_counter!=i){
        printf("correlator counter wrong\n");
        printf("corr_counter=%d     expected=%d\n",i,corr_counter);
        exit(-1);
    }
}

double constant_fit(int n, int Nvar, double *x,int Npar,double  *P){
    
    return P[0];
}


double  sinh_mass(int n, int Nvar, double *x,int Npar,double  *P){
    
    double t=P[0];
    double T=P[1];
    double cp=sinh(x[0]*(t+1.+(1-T)/2));
    double c=sinh(x[0]*(t+(1-T)/2));
    //printf("inside sinh = %f   %f \n", c,cp);
    return c/cp;
}

//in [varialbe] [t] [re_im]
double shift_and_M_eff_acosh(int t,int T , double **in){
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

inline double second_der_time(  int t, int T, double **in){
    return (6*in[t][0]-4*in[(t+1)%T][0]+in[(t+2)%T][0]-4*in[(t-1+T)%T][0]+in[(t-2+T)%T][0]);
}


double laplacian_M_eff_T(  int t, int T, double **in){
    double mass;
    double ct[1],ctp[1],res,tmp_mass, u,d ;
    int i,L0;
    ct[0]=second_der_time( t,  T, in);
    ctp[0]=second_der_time( (t+1)%T,  T, in);
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


inline double der2sym(  int t, int T, double **in){
    return (in[(t+1)%T][0]-2*in[(t)][0]+in[(t-1+T)%T][0]);
}
double der2corr_M_eff_T(  int t, int T, double **in){
    double mass;
    double ct[1],ctp[1],res,tmp_mass, u,d ;
    int i,L0;
    ct[0]=der2sym( t,  T, in);
    ctp[0]=der2sym( (t+1)%T,  T, in);
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
    
    L0=T;
    ct[0]=in[t][0];
    ctp[0]=in[t+1][0];


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

double M_eff_sinh_T_ct_ctp(  int t, int T, double ct,double ctp){
    double mass;
 
    double res,tmp_mass, u,d ;
    

    mass=log(ct/ctp);
    res=1;
    while(res>1e-12){
             u=-1.+exp(-mass*(T-1-2.*t-2.));
             d=-1.+exp(-mass*(T-1-2.*t));
             tmp_mass=log( (ct/ctp) * (u/d)) ;
             res=fabs(tmp_mass - mass);
             mass=tmp_mass;

    }
    return mass;

}


double M_eff_t0_sinh_T(  int t,int t0, int T, double **in){
    double mass;
 
    double ct[1],ctp[1],res,tmp_mass, u,d ;
    int i,L0;
    
    L0=T;
    ct[0]=in[t][0];
    ctp[0]=in[t+1][0];


    mass=log(ct[0]/ctp[0]);

    res=1;
    i=t-t0;
    while(res>1e-12){
             u=-1.+exp(-mass*(L0-1-2.*i-2.));
             d=-1.+exp(-mass*(L0-1-2.*i));
             tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
             res=fabs(tmp_mass - mass);
             mass=tmp_mass;
    }
  
    return mass;

}


double shift_and_M_eff_sinh_T(  int t, int T, double **in){
    double mass;
 
    double ct[1],ctp[1],res,tmp_mass, u,d ;
    int i,L0;
    
    L0=T;
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


double M_eff_log(  int t, int T, double **in){
    double mass;
 
    double ct[1],ctp[1],res,tmp_mass, u,d ;
    int i,L0;
        
    ct[0]=in[t][0];
    ctp[0]=in[t+1][0];

    mass=log(ct[0]/ctp[0]);
   
    return mass;

}


double identity(int t, int T, double ** in){ 
    return in[t][0];

}


double shift_corr(int t, int T, double ** in){ 
    return in[t][0]-in[t+1][0];

}


double M_eff_log_shift(  int t, int T, double **in){
    double mass;
 
    double ct[1],ctp[1],res,tmp_mass, u,d ;
    int i,L0;
        
    ct[0]=in[t][0]-in[t+1][0];
    ctp[0]=in[t+1][0]-in[t+2][0];

    mass=log(ct[0]/ctp[0]);
   
    return mass;

}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//return the line where to read the plateau
int   line_read_plateaux(char **option, const char *corr , int &tmin, int &tmax, int &sep , const  char *namefile_plateaux ){
    
    int line=0;
   //option[3] path   
   //option[6] file
   char namefile[NAMESIZE]; 
   mysprintf(namefile,NAMESIZE,"%s/%s",option[3],namefile_plateaux);
   std::fstream newfile;
   
   newfile.open(namefile,std::ios::in); //open a file to perform read operation using file object
   int match=0;
   if (newfile.is_open()){ //checking whether the file is open
      std::string tp;
      while(getline(newfile, tp)){ //read data from file object and put it into string.
         line++;
         std::vector<std::string> x = split(tp,' ');
         
         std::string  name=option[6]   ;
         std::string  correlator=corr   ;
         if (x.empty()==0){
         if (x[0].compare(name)==0  &&  x[1].compare(correlator)==0){
             tmin=stoi(x[2]);
             tmax=stoi(x[3]);
             sep=stoi(x[4]);
             printf("correlator %s  plateaux %d  %d %d\n", correlator.c_str(),tmin,tmax,sep);
             match++;
             break;
         } 
         }
         
      }
      newfile.close(); //close the file object.
   }
   else{
      error(0==0,1,"correlators_analysis.cpp line_read_plateaux",
            "unable to open %s",namefile);
   }
   //error(match==0,1,"correlators_analysis.cpp line_read_plateaux",
   //      "no match for plateau %s   %s \n in the file %s ",option[6], corr,namefile);
   if (match==0){
      printf("no plateau found for %s in plateau file %s\n", corr,namefile);
      printf("looking for a line:\n %s  %s\n", option[6] ,corr);
      fprintf(stderr,"Please enter the plateau interval and separation:\n"); 
      myscanf(3, (char*) "%d %d %d",&tmin,&tmax,&sep);
      while (tmin>tmax && tmax>= file_head.l0/2 &&  tmin <0){
            fprintf(stderr,"please enter a valid number for tmin and tmax\n");
            myscanf(3,(char*)"%d %d %d",&tmin,&tmax,&sep);
    }
   }
   if (tmin>tmax){
      printf("\n\nerror tmin can not be more than tmax set tmin=tmax-1\n\n");
      tmin=tmax-1;
   }
   if (tmax>file_head.l0){
       printf("\n\n error tmax over T, set  tmax=T/2-1\n\n");
       tmin=file_head.l0/2-1;
   }
   
   
   return line;
    
}




struct fit_result try_fit(char **option,int tmin, int tmax, int sep ,double **corr_ave, double **corr_J,int Njack ,double **chi2, struct fit_type fit_info){
        
    double ***y,***x,*r,**tmp; 
    int i,j;  
    double *chi2j;
    int Npar=fit_info.Npar;  //parameters
    int Nvar=fit_info.Nvar;
    int N=fit_info.N; //functions to fit
    int *en=(int*) malloc(sizeof(int)*fit_info.N); // data to fit for each function
    for (i=0;i< fit_info.N; i++)
        en[0]=(tmax-tmin)/sep +1;
    int en_tot=0;  // total data to fit
    for ( int n=0;n<N;n++)   {  en_tot+=en[n];   }
    
    
    double *guess=(double*) malloc(sizeof(double)*fit_info.Npar);// initial guess for the parameter
    for (i=0;i< fit_info.Npar; i++)
        guess[i]=1.;
    
    if(fit_info.guess.size()>0){
        error(fit_info.guess.size()>Npar,1,"try_fit", " more guess than parametes");
        for (int i=0;i<fit_info.guess.size();i++)
            guess[i]=fit_info.guess[i];
        for (int i=fit_info.guess.size();i<Npar;i++)
            guess[i]=1;
    }
    
    
    chi2j=(double*) malloc(sizeof(double)*Njack);
    
    y=double_malloc_3(Njack, tmax-tmin+1, 2);
    x=double_malloc_3(Njack,en_tot,Nvar);
    
    for (j=0;j<Njack;j++){
        int count=0;
        for (int n=0;n<N;n++){
                for (int e=0;e<en[n];e++){
                    x[j][count][0]=e+tmin;
                    for(int i=0 ; i< fit_info.n_ext_P; i++){
                        x[j][count][i+1]=fit_info.ext_P[i][j];
                    }

                    count++;
                }
        }
    }

   
    double **fit=(double**) malloc(sizeof(double*)*Njack);//result of the fit, the other dimension is allocated by the function non_linear_fit_Nf()

    int count=0;
    for (int n=0;n<N;n++){
        for (i=tmin;i<=tmax;i+=sep){
            for (j=0;j<Njack;j++){
                y[j][count][0]=corr_J[i][j];
                y[j][count][1]= error_jackboot(option[4], Njack, corr_J[i]);   ;//corr_ave[i][1];
            }
            count++;
        }
    }
//     if(fit_info.guess.size()==0)
    guess=guess_for_non_linear_fit_Nf(N, en,x[Njack-1], y[Njack-1] , Nvar,  Npar, fit_info.function,guess ,fit_info);
    //for (j=0;j<Njack;j++){
    for (j=Njack-1;j>=0;j--){    
        //tmp=linear_fit( (tmax-tmin)/sep +1, x, y[j],  1,constant_fit_to_try );

        fit[j]=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess, fit_info );
        if (fit_info.precision_sum >0 )
            chi2j[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],fit[j] , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        else 
            chi2j[j]=compute_chi_non_linear_Nf_kahan(N, en,x[j], y[j],fit[j] , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        int max=0;
        
        while (fabs(chi2j[j]-chi2j[Njack-1])/chi2j[Njack-1]>fit_info.chi2_gap_jackboot && max<fit_info.guess_per_jack){
            std::mt19937 mt_rand(max);
            printf("jack %d has a chi2= %g   while the mean has chi2=%g \n retry\n",j,chi2j[j],chi2j[Njack-1]);
            
            double *guess1=(double*) malloc(sizeof(double)*Npar);
            for (int i=0;i< Npar;i++) guess1[i]=guess[i]+guess[i]* mt_rand()/((double) 10*mt_rand.max() ) ;
            guess1=guess_for_non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess1 ,fit_info);
            
            double *tmp_fit=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess1, fit_info );
            double tmp_chi2;
            if (fit_info.precision_sum>0 )
                tmp_chi2=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp_fit , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
            else 
                tmp_chi2=compute_chi_non_linear_Nf_kahan(N, en,x[j], y[j],tmp_fit , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
//             double tmp_chi2=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp_fit , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
            
            if (tmp_chi2<chi2j[j]){
                chi2j[j]=tmp_chi2;
                for (int i=0;i< Npar;i++) fit[j][i]=tmp_fit[i];
            }
             
            printf("%d  chi= %g  P=\t",j,   chi2j[j] );
            for (int i =0; i< fit_info.Npar;i++){
                printf("%g\t",fit[j][i]);
            }
            printf("\n");
            
                
            free(tmp_fit);
            free(guess1);
            
            max++;
        }
        //chi2j[j]=compute_chisqr((tmax-tmin)/sep +1, x, y[j],  1, tmp, constant_fit_to_try )/((tmax-tmin)/sep +1);
       
    } 
    
        
    struct fit_result fit_out=malloc_fit(fit_info);
    for(i=0;i<Npar;i++){
        for (j=0;j<Njack;j++){
            fit_out.P[i][j]=fit[j][i];
        }
       
       /*fit_out.C[i]=(double**) malloc(sizeof(double*)*Npar);
       for(n=0;n<Npar;n++){     
           fit_out.C[i][n]=(double*) malloc(sizeof(double)*Njack);
           for (j=0;j<Njack;j++){
                fit_out.C[i][n][j]=(*C)[j][i][n];
           }
       }*/
    }
    for (j=0;j<Njack;j++){
        fit_out.chi2[j]=chi2j[j];
    }
    

    free_2(Njack,fit);
    
    free_3(Njack,en_tot,x);
    free_3(Njack, tmax-tmin+1, y);
    free(chi2j);
    
    free(en);free(guess);
    return fit_out;    
}


struct fit_result fit_fun_to_corr(char **option,struct kinematic kinematic_2pt , char* name,const char *description,double **mt,double **r, int Njack,const char *plateaux_masses, FILE *outfile,struct fit_type fit_info)
{
   int yn;
   int tmin=1, tmax=3, sep=1;
   double *m,*fit;
   struct fit_result fit_out;
   char jack_name[NAMESIZE];
   double *chi2;
   
   if (fit_info.plateaux_scan){
       fprintf(fit_info.f_plateaux_scan,"\n\n#tmin tmax chi2 P1 P1err ...\n");
       for(tmin=1;tmin<file_head.l0/2-fit_info.Npar;tmin++){
            for(tmax=tmin+fit_info.Npar;tmax<file_head.l0/2;tmax++){
               fit_result tmp=try_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2,fit_info);
               fprintf(fit_info.f_plateaux_scan,"%d   %d ",tmin,tmax);
               fprintf(fit_info.f_plateaux_scan," %.5g \t",tmp.chi2[Njack-1]);
               for (int i=0; i< fit_info.Npar; i++){
                   fprintf(fit_info.f_plateaux_scan,"%.15g    %.15g    \t",tmp.P[i][Njack-1],error_jackboot(option[4],Njack,tmp.P[i])  );
               }
               free_fit_result( fit_info, tmp);
               fprintf(fit_info.f_plateaux_scan,"\n");               
           }
       }
       
   }
   
   yn=1;
   if ( strcmp(option[1],"see")==0 ){
            while(yn>0){
                printf("#%s 1) mu %g r %d theta %g 2) mu %g r %d theta %g \n",
                       description,
                kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
                kinematic_2pt.k1,kinematic_2pt.r1,kinematic_2pt.mom1);   
                plotting( file_head.l0, mt , &tmin,&tmax, &sep);
                fit_out=try_fit(option ,tmin, tmax, sep , mt, r, Njack,&chi2, fit_info );
                //yn=plotting_fit(file_head.l0, mt , tmin,tmax,m,chi2);
                // to_do:    add test of the fit
                //yn=0;
                free_fit_result( fit_info, fit_out);
            }
   }
   if ( strcmp(option[1],"read_plateaux")==0 ){
        int l=line_read_plateaux(option, description ,  tmin,  tmax,  sep ,plateaux_masses);
        //m=try_linear_fit(option, tmin,  tmax,sep , mt, r, Njack );    
//         fit_out=try_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2,fit_info);
//         yn=1;
//         if(fit_out.chi2[Njack-1]>5.5){
//              free_fit_result( fit_info, fit_out);
//              while(yn>0){
//                    fit_out=try_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2 ,fit_info);
//                    printf("#%s from %s 1) mu %g r %d theta %g 2) mu %g r %d theta %g  line_plateaux=%d\n",
//                         description,name,
//                     kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
//                     kinematic_2pt.k1,kinematic_2pt.r1,kinematic_2pt.mom1,kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1)+1 );   
//                     printf("The chi2 in the range [%d,%d] is %g  consider changing time interval\n",tmin,tmax,fit_out.chi2[Njack-1]);
//                     //plotting( file_head.l0, mt , &tmin,&tmax, &sep);
//                     //yn=plotting_fit(file_head.l0, mt , tmin,tmax,m,chi2);
//                      // to_do:    add test of the fit
//                     yn=0;
//                     if (yn>0){
//                         plotting( file_head.l0, mt , &tmin,&tmax, &sep);
//                         free_fit_result( fit_info, fit_out);
//                     }
//             }
//             //to_do
//             //replace_plateau(  kinematic_2pt ,  description,tmin,tmax,sep);
// 
//          }
//          free_fit_result( fit_info, fit_out);


      
   }
   
   
   fit_out=try_fit(option, tmin,  tmax,sep , mt, r, Njack ,&chi2,fit_info);              
   //fit=give_jack_linear_fit( tmin,  tmax,sep , mt, r, Njack );    
   //printing the correlator and the function   
   double **tif=swap_indices(fit_info.Npar,Njack, fit_out.P);
   double *tmp=(double*) malloc(sizeof(double)*Njack);
   double *x=(double*) malloc(sizeof(double)*fit_info.Nvar);
   int xs=fit_info.Nvar-fit_info.n_ext_P;
   for(int t=1;t<file_head.l0/2;t++){  
           fprintf(outfile,"%d   %.15e    %.15e\t",t,mt[t][0],mt[t][1]);
           // variables and external parameters 
           x[0]=t;
           for(int i=0 ; i< fit_info.n_ext_P; i++)
               x[i+xs]=fit_info.ext_P[i][Njack-1];
           
           for (int j=0;j<Njack;j++)
               tmp[j]= fit_info.function(0 ,fit_info.Nvar,x,fit_info.Nvar, tif[j] );
           double *f_val=mean_and_error(option[4],Njack,tmp );
           
           fprintf(outfile,"   %.15e    %.15e\t",f_val[0],f_val[1]);
           free(f_val);
           
           x[0]=t+0.33;
           for (int j=0;j<Njack;j++)
               tmp[j]= fit_info.function(0 ,fit_info.Nvar,x,fit_info.Nvar, tif[j] );
           f_val=mean_and_error(option[4],Njack,tmp );
           fprintf(outfile,"%.15e   %.15e    %.15e\t",x[0],f_val[0],f_val[1]);
           free(f_val);
           x[0]=t+0.66;
           for (int j=0;j<Njack;j++)
               tmp[j]= fit_info.function(0 ,fit_info.Nvar,x,fit_info.Nvar, tif[j] );
           f_val=mean_and_error(option[4],Njack,tmp );
           fprintf(outfile,"%.15e   %.15e    %.15e\n",x[0],f_val[0],f_val[1]);
           free(f_val);

   }
   free(tmp);free(x);
   free_2(Njack,tif);
   
   m=mean_and_error(option[4],Njack, fit_out.chi2);
   fprintf(outfile,"\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n",description,tmin,tmax,m[0],m[1]);
   if (fit_info.plateaux_scan)fprintf(fit_info.f_plateaux_scan,"\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n",description,tmin,tmax,m[0],m[1]);
   free(m);
   for (int i=0; i< fit_info.Npar; i++){
       m=mean_and_error(option[4],Njack, fit_out.P[i]);
       fprintf(outfile,"%.15g    %.15g    \t",m[0],m[1]);
       if (fit_info.plateaux_scan) fprintf(fit_info.f_plateaux_scan,"%.15g    %.15g    \t",m[0],m[1]);
       if (i==0){
           if (fit_info.verbosity>=0)
               printf("#%s  fit in [%d,%d]:  %.15g    %.15g\n",description,tmin,tmax,m[0],m[1]);
       }
       free(m);

   }
   fprintf(outfile,"\n");
   if (fit_info.plateaux_scan) fprintf(fit_info.f_plateaux_scan,"\n");
   
   
   corr_counter++;
    
   return fit_out;
}



 double   *plateau_correlator_function(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,const char  *plateaux_masses,FILE *outfile,  int index , const char *description , double (*fun)(int ,int  , double ** ),  FILE * file_jack, struct fit_type fit_info){
// jackknife plateau_correlator_function(char **option ,struct kinematic kinematic_2pt , char* name, double ****conf_jack, int Njack ,const char  *plateaux_masses,FILE *outfile,  int index , const char *description , double (*fun)(int ,int  , double ** ),  FILE * file_jack){    
   /*int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(*plateaux_masses,line);
   */
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=double_malloc_2(file_head.l0, Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);

   fprintf(outfile," \n\n");
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
           //fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);

   }

   //struct fit_type fit_info;
   fit_info.Nvar=1;
   fit_info.Npar=1;
   fit_info.N=1;
   fit_info.Njack=Njack;
   fit_info.function=constant_fit;
   fit_info.n_ext_P=0;
  
   //fit=fit_plateaux(option, kinematic_2pt ,  name,description/*"M_{PS}^{ll}"*/,mt,r,  Njack,*plateaux_masses,outfile);
   struct fit_result fit_out=fit_fun_to_corr(option, kinematic_2pt ,  name, description, mt, r,  Njack, plateaux_masses, outfile, fit_info);
   
   fit=(double*) malloc(sizeof(double)*Njack);
   for (j=0;j<Njack;j++)
       fit[j]=fit_out.P[0][j];
   

   fwrite(fit_out.P[0],sizeof(double),Njack, file_jack );
   
   free_fit_result(fit_info,fit_out);  

   for(i=1;i<file_head.l0/2;i++)
       free(mt[i]);
   free(mt);
   free_2(file_head.l0,r);
  
   fflush(outfile);
   /*
    if ( strcmp(option[1],"read_plateaux")==0 ){
     fclose(*plateaux_masses);
     *plateaux_masses=open_file(kinematic_2pt.plateau_m_ll,"r");

    }*/
//    jackknife tmp(option[4], Njack, fit);
//    free(fit);
//    return tmp;
   return fit;    
    
}

struct fit_result fit_function_to_corr(char **option ,struct kinematic kinematic_2pt ,  char* name, double ****conf_jack ,const char  *plateaux_masses,FILE *outfile,  int index, int re_im , const char *description , struct fit_type fit_info,  FILE * file_jack ){

   int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   /*if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(*plateaux_masses,line);
   */
   int Njack=fit_info.Njack;
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);
   
   fprintf(outfile," \n\n");
   fprintf(outfile,"#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
              
              r[i][j]= conf_jack[j][index][i][re_im];
             
           }
           
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           //fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);

   }
   
   struct fit_result fit_out=fit_fun_to_corr(option, kinematic_2pt ,  name, description, mt, r,  Njack, plateaux_masses, outfile, fit_info);
   //fit=fit_plateaux(option, kinematic_2pt ,  name,description/*"M_{PS}^{ll}"*/,mt,r,  Njack,*plateaux_masses,outfile);
   
   fwrite(fit_out.P[0],sizeof(double),Njack, file_jack );
   
   
   
   for(i=1;i<file_head.l0/2;i++)
       free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
       free(r[i]);
   free(r);
   
   fflush(outfile);
   
    /*if ( strcmp(option[1],"read_plateaux")==0 ){
     fclose(*plateaux_masses);
     *plateaux_masses=open_file(kinematic_2pt.plateau_m_ll,"r");

    }*/
    return fit_out;   
    
} 




//fit the function of fit info.function() to a combination of correlator  y= fun_of_corr( jack_index, data[jack][corr][time][reim], time, fit_info)
//  
struct fit_result fit_fun_to_fun_of_corr(char **option ,struct kinematic kinematic_2pt ,  char* name, double ****conf_jack ,const char  *plateaux_masses,FILE *outfile,  double fun_of_corr(int, double****,int, struct fit_type ) , const char *description , struct fit_type fit_info,  FILE *file_jack ){
/*
int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   if ( strcmp(option[1],"read_plateaux")==0 )
   	go_to_line(*plateaux_masses,line);
  */ 
   int Njack=fit_info.Njack;
   double **r,*m,**mt,*fit;
   int i,j,yn;
   
   error(fit_info.N!=1,1,"fit_fun_to_fun_of_corr", "multiple correlator fit not implemented");
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);
   fprintf(outfile," \n\n");
   fprintf(outfile,"#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
              
              //r[i][j]= conf_jack[j][index][i][re_im];
              r[i][j]= fun_of_corr(j,conf_jack,i,fit_info);
             
           }
           
           if( strcmp(option[4],"jack")==0)
               mt[i]=mean_and_error_jack(Njack, r[i]);
           if( strcmp(option[4],"boot")==0)
               mt[i]=mean_and_error_boot(Njack, r[i]);
           //fprintf(outfile,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);

   }
   
   struct fit_result fit_out=fit_fun_to_corr(option, kinematic_2pt ,  name, description, mt, r,  Njack, plateaux_masses, outfile, fit_info);
   //fit=fit_plateaux(option, kinematic_2pt ,  name,description/*"M_{PS}^{ll}"*/,mt,r,  Njack,*plateaux_masses,outfile);
   
   
   fwrite(fit_out.P[0],sizeof(double),Njack, file_jack );
   
   
   for(i=1;i<file_head.l0/2;i++)
       free(mt[i]);
   free(mt);
   for(i=0;i<file_head.l0;i++)
       free(r[i]);
   free(r);
   
   fflush(outfile);
   
   /* if ( strcmp(option[1],"read_plateaux")==0 ){
     fclose(*plateaux_masses);
     *plateaux_masses=open_file(kinematic_2pt.plateau_m_ll,"r");

    }*/
    return fit_out;   
    
} 


void add_correlators(char **option , int& ncorr_conf_jack, double ****&conf_jack , double **fun_of_corr(int, double****,int, struct fit_type ) ,  struct fit_type fit_info ){
    
    int correlators_out=ncorr_conf_jack+fit_info.N;
    int Njack=fit_info.Njack;
    
    double ****corr_out=calloc_corr(Njack, correlators_out, file_head.l0 );
    //copy the first part
    for(int j=0; j<Njack; j++){
        for (int v=0; v<ncorr_conf_jack; v++){
            for(int t=0; t<file_head.l0; t++){
                corr_out[j][v][t][0]=conf_jack[j][v][t][0];
                corr_out[j][v][t][1]=conf_jack[j][v][t][1];
            }
        }
    }
    
    // for(int n=0; n<fit_info.N; n++){
    //     fit_info.n=n; // this goes in to fun_of_corr
    //     for(int j=0; j<Njack; j++){
    //         for(int t=0; t<file_head.l0; t++){
    //             corr_out[j][ncorr_conf_jack+ n][t][0]=fun_of_corr(j,conf_jack,t,fit_info);
    //             corr_out[j][ncorr_conf_jack+ n][t][1]=0;
    //         }
    //     }
    // }
    
    
    for(int j=0; j<Njack; j++){
        for(int t=0; t<file_head.l0; t++){
            for(int n=0; n<fit_info.N; n++){
                double **r=fun_of_corr(j,conf_jack,t,fit_info);
                corr_out[j][ncorr_conf_jack+ n][t][0]=r[n][0];
                corr_out[j][ncorr_conf_jack+ n][t][1]=r[n][1];
                free(r);
            }
        }
    }
    
    free_corr(Njack, ncorr_conf_jack, file_head.l0 ,conf_jack);
    conf_jack=corr_out;
    ncorr_conf_jack=correlators_out;
}


void zero_corr(double *zeros, int Njack, FILE *jack_file){
    fwrite(zeros,sizeof(double),Njack, jack_file );
    corr_counter++;
}
