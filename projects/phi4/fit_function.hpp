#ifndef fit_function_H
#define fit_function_H

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

#include "header_phi4.hpp"
using namespace std;



double M_finite_volume(int n, int Nvar, double *x,int Npar,double  *P){
    double M=P[0]  ;// mass at inif volume
    double L=x[0];
    //int K1=dbesk1(M*L);
    double K1=exp(-M*L)/sqrt(M*L);
    double r=M+P[1]*K1/(M*L) ;
    return r;
}

double a_luscher(int n, int Nvar, double *x,int Npar,double  *P){
    double a=P[0]  ;// mass at inif volume
    double L=x[0];
    double mu=x[1]*x[2]/(x[1]+x[2]);
    //int K1=dbesk1(M*L);
    double r=-(2.*pi_greco * a)  / (mu*L*L*L);
    r *=  1 - 2.837297  *(a/L) + 6.375183 *(a/L) *(a/L);
    return r;
    
}


void print_fit_output(char **argv,vector<data_phi> gjack ,struct fit_type fit_info, const char* label, struct fit_result fit_out, int *en, double ***x, double ***y){
    int Npar=fit_info.Npar;
    int Nvar=fit_info.Nvar;
    int Njack=gjack[0].Njack;
    int N=fit_info.N;
   
    char namefile[NAMESIZE];
    FILE *f;
    
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_data.txt",argv[3], label);
    f=open_file(namefile,"w+");
    int count=0;
    for (int n=0;n<N;n++){
        for (int e=0;e<en[n];e++){
            fprintf(f," %g   %g   %g\n ",x[Njack-1][e+count][0], y[Njack-1][e+count][0], y[Njack-1][e+count][1] );
        }
        count+=en[n];
        fprintf(f,"\n\n");
    }
            
    fclose(f);
    ////////////// parameter print and correlation matrix
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_P.tex",argv[3], label);
    f=open_file(namefile,"w+");
    fprintf(f,"\\begin{gather}\n");
    fprintf(f,"\\chi^2/d.o.f.=%g\\pm %.2g \\\\ \n", fit_out.chi2[Njack-1], error_jackboot(argv[1], Njack, fit_out.chi2));
    for (int i=0;i<Npar;i++){
        fprintf(f,"P[%d]=%g\\pm %.2g \\\\ \n",i, fit_out.P[i][Njack-1], error_jackboot(argv[1], Njack, fit_out.P[i]));
    }
    fprintf(f,"\\end{gather}\n");
    double **cov=covariance(argv[1], Npar, Njack, fit_out.P);
    fprintf(f,"{\\tiny\\begin{gather}\n C=\\begin{pmatrix}\n");
    for (int i=0;i<fit_info.Npar;i++){
        for (int k=0;k<i;k++)
                    cov[i][k]/=sqrt(cov[i][i]*cov[k][k]);
       for (int k=i+1;k<fit_info.Npar;k++)
                    cov[i][k]/=sqrt(cov[i][i]*cov[k][k]);
    }
    

     for (int i=0;i<fit_info.Npar;i++){
        for (int j=0;j<fit_info.Npar;j++){
            if (j==0)  fprintf(f,"%.3g",  cov[i][j] );
            else       fprintf(f,"& %.3g", cov[i][j] );
        }
        if (i!=fit_info.Npar) fprintf(f,"\\\\ \n");
        else fprintf(f,"\n");
    }
    fprintf(f,"\\end{pmatrix}\n\\end{gather}}\n");
    free_2(Npar,cov);
    fclose(f);
    /////////fit band
    
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_out.txt",argv[3], label);
    f=open_file(namefile,"w+");
    
    double **tif=swap_indices(fit_info.Npar,Njack,fit_out.P);
    for (int i=0 ; i<100; i++){
        double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
        double *tmpy=(double*) malloc(sizeof(double*)* Njack);
        tmpx[0]=10+i*0.5;
        for(int i=0 ; i< fit_info.n_ext_P; i++)
            tmpx[i+1]=fit_info.ext_P[i][Njack-1];
                
        for (int j=0;j<Njack;j++){
            
            tmpy[j]=fit_info.function(N,Nvar,tmpx,Npar,tif[j]);//N, Nvar, x ,Npar,P
        }
        fprintf(f,"%g  \t %g  %g\n",tmpx[0],tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
        free(tmpy);free(tmpx);
    }
    free_2(Njack,tif);
    fclose(f);  
       
    
}

struct fit_result fit_m0_phi4(char **argv, vector<cluster::IO_params> params ,vector<data_phi> gjack ,struct fit_type fit_info, const char* label){
    int Npar=fit_info.Npar;
    int Nvar=fit_info.Nvar;
    int Njack=gjack[0].Njack;
    int N=fit_info.N;
    ////// allocation
    int *en=(int*) malloc(sizeof(int)*fit_info.N);// we need to init en and en_tot to allocate the other 
        for (int e=0;e< fit_info.N; e++){     en[e]=gjack.size();}
        int en_tot=0;      for ( int n=0;n<N;n++)   {  en_tot+=en[n];   }// total data to fit
    double ***y=double_malloc_3(Njack, en_tot, 2);// 2 is mean/error
    double ***x=double_malloc_3(Njack,en_tot,Nvar);
    struct fit_result fit_out=malloc_fit(fit_info);
    double *guess=(double*) malloc(sizeof(double)*Npar);
    double **fit=(double**) malloc(sizeof(double*)*Njack);//result of the fit, the other dimension is allocated by the function non_linear_fit_Nf()

    ////// allocation end
    /////////////// init
    
    for (int i=0;i<Npar;i++)
        guess[i]=1;//rand()/((double)RAND_MAX);
    
    
    //init x
    for (int j=0;j<Njack;j++){
       int count=0;
       for (int n=0;n<N;n++){
            for (int e=0;e<en[n];e++){
                x[j][count][0]=params[e].data.L[1];
                for(int i=0 ; i< fit_info.n_ext_P; i++){
                    x[j][count][i+1]=fit_info.ext_P[i][j];
                }
                count++;
            }
       }
    }
    
    
    int count=0;
    for (int n=0;n<N;n++){
        for (int e=0;e<en[n];e++){
            for (int j=0;j<Njack;j++){
                y[j][e+count][0]=gjack[e].jack[1][j];
                y[j][e+count][1]= error_jackboot(argv[1], gjack[e].Njack, gjack[e].jack[1]);   ;//corr_ave[i][1];
                if (j==Njack-1) {printf(" %g   %g   %g\n",x[j][e+count][0], y[j][e+count][0], y[j][e+count][1] );}
                
                
            }
        }
        count+=en[n];
    }
    //////  init end
    
    ///////////////// the fit 
    // scan the parameter of the fit with the last jack
    guess=guess_for_non_linear_fit_Nf(N, en,x[Njack-1], y[Njack-1] , Nvar,  Npar, fit_info.function,guess );
    
    
    for (int j=0;j<Njack;j++){
        
        fit[j]=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess );
        //tmp=non_linear_fit_Nf_sigmax( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_iterative( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_covariance( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess ,cov_yx1);
        //tmp=non_linear_fit_Nf_covariance(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess ,cov1);
       
        
        fit_out.chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],fit[j] , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        
        // we do not need the covariance of the fit, it will be computed with jackboot
        //double **C=covariance_non_linear_fit_Nf(N, en,x[j], y[j],fit[j] , Nvar,  Npar, fit_info.function );            
        //for(int i=0;i<Npar;i++)
        //    for(int k=0;k<Npar;k++)
        //        fit_out.C[i][k][j]=C[i][k];
        //free_2(Npar, C);
        
        
    }
    for(int i=0;i<Npar;i++)
       for (int j=0;j<Njack;j++)
           fit_out.P[i][j]=fit[j][i];
    
     
    /////////////////////////////////////////////////////////////////////writing the result
    print_fit_output(argv,   gjack , fit_info,  label,  fit_out , en,x,y);

    ////// free
    free(en);
    //free(chi2j);
    free_3(Njack, en_tot, y);
    free_3(Njack, en_tot, x);
    free(guess);
    free_2(Njack, fit);

    ////// free end
     
       
    return fit_out;   
   
  
}

double M0_finite_volume_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack ){
    return gjack[e].jack[1][j];
}

double M1_finite_volume_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack ){
    return gjack[e].jack[2][j];
}

double DE_00_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack ){
    return gjack[e].jack[4][j]-2*gjack[e].jack[1][j];
}


struct fit_result fit_data(char **argv, vector<cluster::IO_params> params ,vector<data_phi> gjack, double lhs_fun(int, int, int ,vector<cluster::IO_params> ,vector<data_phi> ) , struct fit_type fit_info, const char* label){
    int Npar=fit_info.Npar;
    int Nvar=fit_info.Nvar;
    int Njack=gjack[0].Njack;
    int N=fit_info.N;
    ////// allocation
    int *en=(int*) malloc(sizeof(int)*fit_info.N);// we need to init en and en_tot to allocate the other 
        for (int e=0;e< fit_info.N; e++){     en[e]=gjack.size();}
        int en_tot=0;      for ( int n=0;n<N;n++)   {  en_tot+=en[n];   }// total data to fit
    double ***y=double_malloc_3(Njack, en_tot, 2);// 2 is mean/error
    double ***x=double_malloc_3(Njack,en_tot,Nvar);
    struct fit_result fit_out=malloc_fit(fit_info);
    double *guess=(double*) malloc(sizeof(double)*Npar);
    double **fit=(double**) malloc(sizeof(double*)*Njack);//result of the fit, the other dimension is allocated by the function non_linear_fit_Nf()

    ////// allocation end
    /////////////// init
    
    for (int i=0;i<Npar;i++)
        guess[i]=1;//rand()/((double)RAND_MAX);
    
    
    //init x
    for (int j=0;j<Njack;j++){
       int count=0;
       for (int n=0;n<N;n++){
            for (int e=0;e<en[n];e++){
                x[j][count][0]=(double) params[e].data.L[1];
                //other var
                //
                for(int i=0 ; i< fit_info.n_ext_P; i++){
                    x[j][count][i+1]=fit_info.ext_P[i][j];
                }
                count++;
            }
       }
    }
    
    
    int count=0;
    for (int n=0;n<N;n++){
        for (int e=0;e<en[n];e++){
            double *tmpj=(double*) malloc(sizeof(double)*Njack);
            for (int j=0;j<Njack;j++){
                //y[j][e+count][0]=gjack[e].jack[1][j];
                tmpj[j]=lhs_fun(n,e,j,params,gjack);
            }
            double err=error_jackboot(argv[1], Njack, tmpj); 
            for (int j=0;j<Njack;j++){
                y[j][e+count][0]=tmpj[j];
                y[j][e+count][1]= err;
            }
            printf(" %g   %g   %g\n",x[Njack-1][e+count][0], y[Njack-1][e+count][0], y[Njack-1][e+count][1] );   
            free(tmpj);
        }
        count+=en[n];
    }
    //////  init end
    
    ///////////////// the fit 
    // scan the parameter of the fit with the last jack
    guess=guess_for_non_linear_fit_Nf(N, en,x[Njack-1], y[Njack-1] , Nvar,  Npar, fit_info.function,guess );
    
    
    for (int j=0;j<Njack;j++){
        
        fit[j]=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess );
        //tmp=non_linear_fit_Nf_sigmax( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_iterative( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_covariance( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess ,cov_yx1);
        //tmp=non_linear_fit_Nf_covariance(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess ,cov1);
       
        
        fit_out.chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],fit[j] , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        
        // we do not need the covariance of the fit, it will be computed with jackboot
        //double **C=covariance_non_linear_fit_Nf(N, en,x[j], y[j],fit[j] , Nvar,  Npar, fit_info.function );            
        //for(int i=0;i<Npar;i++)
        //    for(int k=0;k<Npar;k++)
        //        fit_out.C[i][k][j]=C[i][k];
        //free_2(Npar, C);
        
        
    }
    for(int i=0;i<Npar;i++)
       for (int j=0;j<Njack;j++)
           fit_out.P[i][j]=fit[j][i];
    
     
    /////////////////////////////////////////////////////////////////////writing the result
    print_fit_output(argv,   gjack , fit_info,  label,  fit_out , en,x,y);

    ////// free
    free(en);
    //free(chi2j);
    free_3(Njack, en_tot, y);
    free_3(Njack, en_tot, x);
    free(guess);
    free_2(Njack, fit);

    ////// free end
     
       
    return fit_out;   
   
  
}


#endif
