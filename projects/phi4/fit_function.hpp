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




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// functions lhs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double M0_finite_volume_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack,  struct fit_type fit_info){
    return gjack[e].jack[1][j];
}

double M1_finite_volume_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack, struct fit_type fit_info ){
    return gjack[e].jack[2][j];
}

double DE_00_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack , struct fit_type fit_info){
    return gjack[e].jack[4][j]-2*gjack[e].jack[1][j];
}

double muDE_00_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack, struct fit_type fit_info ){
    double DE=gjack[e].jack[4][j]-2*gjack[e].jack[1][j];
    return DE*gjack[e].jack[1][j]/2.;
}

double muDE_00_infm_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack, struct fit_type fit_info ){
    double DE=gjack[e].jack[4][j]-2*gjack[e].jack[1][j];
    double mu=fit_info.ext_P[0][j]*fit_info.ext_P[1][j]/(fit_info.ext_P[0][j]+fit_info.ext_P[1][j]);
    return DE*mu;
}


double muDE_01_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack, struct fit_type fit_info ){
    double DE=gjack[e].jack[19][j]-gjack[e].jack[1][j]-gjack[e].jack[2][j];
    double mu=gjack[e].jack[1][j]*gjack[e].jack[2][j]/(gjack[e].jack[1][j]+gjack[e].jack[2][j]);
    return DE*mu;
}

double muDE_01_div_shift_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack, struct fit_type fit_info ){
    double DE=gjack[e].jack[85][j]-gjack[e].jack[1][j]-gjack[e].jack[2][j];
    double mu=gjack[e].jack[1][j]*gjack[e].jack[2][j]/(gjack[e].jack[1][j]+gjack[e].jack[2][j]);
    return DE*mu;
    
}



template<int id>
double lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack, struct fit_type fit_info ){
    return gjack[e].jack[id][j];
}
    
double a_01_BH_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack, struct fit_type fit_info ){
    //return gjack[e].jack[12][j];// 0t_8tT_2
    //return gjack[e].jack[23][j];//   03t16
    return gjack[e].jack[53][j];//   03t16_shifted1
}
double a_00_BH_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack, struct fit_type fit_info ){
    //return gjack[e].jack[10][j];
    return gjack[e].jack[21][j];//   03t16
}
double a_11_BH_lhs(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack, struct fit_type fit_info ){
    return gjack[e].jack[11][j];
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// functions rhs
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double M_finite_volume(int n, int Nvar, double *x,int Npar,double  *P){
    double M=P[0]  ;// mass at inif volume
    double L=x[0];
    //int K1=dbesk1(M*L);
    double K1=exp(-M*L)/sqrt(M*L);
    double r=M+P[1]*K1/(M*L) ;
    return r;
}


double muDE_rhs(int n, int Nvar, double *x,int Npar,double  *P){
    double a=P[0]  ;// mass at inif volume
    double L=x[0];
    //double mu=x[Nvar]*x[Nvar+1]/(x[Nvar]+x[Nvar+1]);
    //int K1=dbesk1(M*L);
    double r=-(2.*pi_greco * a)  / (L*L*L);
    r *=  1 - 2.837297  *(a/L) + 6.375183 *(a/L) *(a/L);
    return r;
    
}

double a_luscher_infm(int n, int Nvar, double *x,int Npar,double  *P){
    double a=P[0]  ;// mass at inif volume
    double L=x[0];
    double mu=x[Nvar]*x[Nvar+1]/(x[Nvar]+x[Nvar+1]);
    //int K1=dbesk1(M*L);
    double r=-(2.*pi_greco * a)  / (mu*L*L*L);
    r *=  1 - 2.837297  *(a/L) + 6.375183 *(a/L) *(a/L);
    return r;
    
}
double a_luscher(int n, int Nvar, double *x,int Npar,double  *P){
    double a=P[0]  ;// mass at inif volume
    double L=x[0];
    double mu=x[1]*x[1]/(x[1]+x[1]);
    //int K1=dbesk1(M*L);
    double r=-(2.*pi_greco * a)  / (mu*L*L*L);
    r *=  1 - 2.837297  *(a/L) + 6.375183 *(a/L) *(a/L);
    return r;
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// print output
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void print_fit_output(char **argv,vector<data_phi> gjack ,struct fit_type fit_info, const char* label, struct fit_result fit_out, int *en, double ***x, double ***y,  vector<cluster::IO_params> params, std::vector<int> myen){
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
            fprintf(f," %g   %g   %g   %d\n ",x[Njack-1][e+count][0], y[Njack-1][e+count][0], y[Njack-1][e+count][1] , params[myen[e]].data.L[0]);
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
    /////////fit band L
    
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_L.txt",argv[3], label);
    f=open_file(namefile,"w+");
    
    double **tif=swap_indices(fit_info.Npar,Njack,fit_out.P);
    for (int i=0 ; i<100; i++){
        double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
        double *tmpy=(double*) malloc(sizeof(double*)* Njack);
                
        for (int j=0;j<Njack;j++){
            tmpx[0]=10+i*0.5;
            tmpx[1]=gjack[0].jack[1][j];//m0
            tmpx[2]=gjack[0].jack[2][j];//m1
            tmpx[3]=gjack[0].jack[4][j];//E20
            tmpx[4]=gjack[0].jack[5][j];//E21
            tmpx[5]=(double) params[0].data.L[0];//T
            
            for(int i=fit_info.Nvar ; i< fit_info.n_ext_P; i++)
                tmpx[i]=fit_info.ext_P[fit_info.Nvar][j];
            
            tmpy[j]=fit_info.function(N,Nvar,tmpx,Npar,tif[j]);//N, Nvar, x ,Npar,P
        }
        fprintf(f,"%g  \t %g  %g\n",tmpx[0],tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
        free(tmpy);free(tmpx);
    }
    free_2(Njack,tif);
    fclose(f); 
    /////////fit band L
    
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_T.txt",argv[3], label);
    f=open_file(namefile,"w+");
    
    tif=swap_indices(fit_info.Npar,Njack,fit_out.P);
    for (int i=0 ; i<100; i++){
        double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
        double *tmpy=(double*) malloc(sizeof(double*)* Njack);
        
        for (int j=0;j<Njack;j++){
            tmpx[0]=(double) params[0].data.L[1];//L
            tmpx[1]=gjack[0].jack[1][j];//m0
            tmpx[2]=gjack[0].jack[2][j];//m1
            tmpx[3]=gjack[0].jack[4][j];//E20
            tmpx[4]=gjack[0].jack[5][j];//E21
            tmpx[5]=10+i*1;//T
            
            for(int i=fit_info.Nvar ; i< fit_info.n_ext_P; i++)
                tmpx[i]=fit_info.ext_P[fit_info.Nvar][j];
            
            tmpy[j]=fit_info.function(N,Nvar,tmpx,Npar,tif[j]);//N, Nvar, x ,Npar,P
        }
        fprintf(f,"%g  \t %g  %g\n",tmpx[5],tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
        free(tmpy);free(tmpx);
    }
    free_2(Njack,tif);
    fclose(f);  
    ////////// end fit band T
       
    
}


struct fit_result fit_data(char **argv, vector<cluster::IO_params> params ,vector<data_phi> gjack, double lhs_fun(int, int, int ,vector<cluster::IO_params> ,vector<data_phi>,struct fit_type ) , struct fit_type fit_info, const char* label, std::vector<int> myen){
    int Npar=fit_info.Npar;
    int Nvar=fit_info.Nvar+fit_info.n_ext_P;
    int Njack=gjack[0].Njack;
    int N=fit_info.N;
    ////// allocation
    int *en=(int*) malloc(sizeof(int)*fit_info.N);// we need to init en and en_tot to allocate the other 
    for (int e=0;e< fit_info.N; e++){     en[e]=myen.size();}
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
                x[j][count][0]=(double) params[myen[e]].data.L[1];//L
                x[j][count][1]=gjack[myen[e]].jack[1][j];//m0
                x[j][count][2]=gjack[myen[e]].jack[2][j];//m1
                x[j][count][3]=gjack[myen[e]].jack[4][j];//E20
                x[j][count][4]=gjack[myen[e]].jack[5][j];//E21
                x[j][count][5]=(double) params[myen[e]].data.L[0];//T
                for(int i=fit_info.Nvar ; i< fit_info.n_ext_P; i++){
                    x[j][count][i]=fit_info.ext_P[i-fit_info.Nvar][j];
                }
                //copy the other variables after 
                 //other var
                //x[j][count][fit_info.n_ext_P]=(double) params[e].data.L[1];
                
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
                tmpj[j]=lhs_fun(n,myen[e],j,params,gjack,fit_info);
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
    guess=guess_for_non_linear_fit_Nf(N, en,x[Njack-1], y[Njack-1] , fit_info.Nvar,  Npar, fit_info.function,guess );
    
    
    for (int j=0;j<Njack;j++){
        
        fit[j]=non_linear_fit_Nf(N, en,x[j], y[j] , fit_info.Nvar,  Npar, fit_info.function,guess );
        //tmp=non_linear_fit_Nf_sigmax( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_iterative( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_covariance( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess ,cov_yx1);
        //tmp=non_linear_fit_Nf_covariance(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess ,cov1);
       
        
        fit_out.chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],fit[j] , fit_info.Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        
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
    print_fit_output(argv,   gjack , fit_info,  label,  fit_out , en,x,y, params,myen);

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
