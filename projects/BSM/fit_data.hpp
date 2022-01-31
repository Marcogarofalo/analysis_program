#ifndef fit_data_H
#define fit_data_H

#include "header_BSM.hpp"
#include "resampling.hpp"
#include <random>



double rhs_critical_eta_mu_m0(int n, int Nvar, double *x,int Npar,double  *P){
    double a0=P[0], r0=P[1], P2=P[2];
    double eta=x[3], mu=x[5], m0=x[6];
    double r;
    if (n==0) //r_AWI
        r= P[0]+P[1]*eta+P[2]*mu+P[3]*m0;
    else if (n==1) //m_pcac
        r= P[4]+P[5]*eta+P[6]*mu+P[7]*m0;
    else{ r=0; exit(1);}
        
    return r;
    
}


double rhs_critical_eta_mu_m0_shifted(int n, int Nvar, double *x,int Npar,double  *P){
    double a0=P[0], r0=P[1], P2=P[2];
    double eta=x[3], mu=x[5], m0=x[6];
    double r;
    double eta_cr=P[0];
    double m0_cr=P[1];
    double eta_sub=eta-eta_cr;
    double m0_sub=m0-m0_cr;
    if (n==0) //r_AWI
        r= P[2]*eta_sub+P[3]*m0_sub+P[6]*mu;
    else if (n==1) //m_pcac
        r= P[4]*eta_sub+P[5]*m0_sub+P[7]*mu;
    else{ r=0; exit(1);}
        
    return r;
    
}


double rhs_NG_mpcac_MPS2(int n, int Nvar, double *x,int Npar,double  *P){
    double eta=x[3], mu=x[5], m0=x[6];
    double r;
    double eta_cr=x[Nvar-2];
    double m0_cr=x[Nvar-1];
    double eta_sub=eta-eta_cr;
    double m0_sub=m0-m0_cr;
    if (n==0) //m_pcac
        r= P[0]+P[2]*eta_sub+P[4]*mu;
    else if (n==1) //MPS^2
        r= P[1]+P[3]*eta_sub+P[5]*mu+P[6]*eta_sub*eta_sub;
    else{ r=0; exit(1);}
        
    return r;
    
}

double rhs_critical_eta_mu_m0_simple(int n, int Nvar, double *x,int Npar,double  *P){
    double eta=x[3], mu=x[5], m0=x[6];
    double r;
    double eta_cr=P[0];
    double m0_cr=P[1];
    double eta_sub=eta-eta_cr;
    double m0_sub=m0-m0_cr;
    if (n==0) //r_AWI
        r= P[2]*eta_sub+P[3]*m0_sub;
    else if (n==1) //m_pcac
        r= P[4]*m0_sub+P[5]*mu;
    else{ r=0; exit(1);}
        
    return r;
    
}


double lhs_fit_two_func(int n, int e , int j , vector<header_BSM> params,vector<data_BSM> gjack, struct fit_type fit_info ){
    double r;
    int n1=fit_info.corr_id[0];
    int n2=fit_info.corr_id[1];
    if(n==0)
        r= gjack[e].jack[n1][j]; //r_AWI
    else if( n==1)
        r= gjack[e].jack[n2][j]; //m_pcac
    else{ r=0; exit(1);}    
    return r;
}

double lhs_critical_eta_mu_m0(int n, int e , int j , vector<header_BSM> params,vector<data_BSM> gjack, struct fit_type fit_info ){
    double r;
    if(n==0)
        r= gjack[e].jack[3][j]; //r_AWI
    else if( n==1)
        r= gjack[e].jack[4][j]; //m_pcac
    else{ r=0; exit(1);}    
    return r;
}

double lhs_critical_eta_mu_m0_loc(int n, int e , int j , vector<header_BSM> params,vector<data_BSM> gjack, struct fit_type fit_info ){
    double r;
    if(n==0)
        r= gjack[e].jack[7][j]; //r_AWI
    else if( n==1)
        r= gjack[e].jack[8][j]; //m_pcac
    else{ r=0; exit(1);}    
    return r;
}
 
double lhs_mpcac_MPS2(int n, int e , int j , vector<header_BSM> params,vector<data_BSM> gjack, struct fit_type fit_info ){
    double r;
    if(n==0)
        r= gjack[e].jack[4][j]; //m_pcac
    else if( n==1){
        r= gjack[e].jack[0][j]; //MPS
        r=r*r;  //MPS2
    }
    else{ r=0; exit(1);}    
    return r;
}   

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// print output
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_fit_output(char **argv,vector<data_BSM> gjack ,struct fit_type fit_info, const char* label, struct fit_result fit_out, int *en, double ***x, double ***y,  vector<header_BSM> params, std::vector<int> myen){
    int Npar=fit_info.Npar;
    int Nvar=fit_info.Nvar+fit_info.n_ext_P;
    int Njack=gjack[0].Njack;
    int N=fit_info.N;
    
    char namefile[NAMESIZE];
    FILE *f;
    
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_data.txt",argv[3], label);
    f=open_file(namefile,"w+");
    int count=0;
    for (int n=0;n<N;n++){
        for (int e=0;e<en[n];e++){
            fprintf(f," %d   \t ",n);
            fprintf(f," %g   %g   \t ", y[Njack-1][e+count][0], y[Njack-1][e+count][1] );
            for (int v=0;v<Nvar;v++){
                fprintf(f," %g   \t ",x[Njack-1][e+count][v]);
            }
            fprintf(f,"    \n ");
        }
        count+=en[n];
        fprintf(f,"\n\n");
    }
    
    fclose(f);
    ////////////// parameter print and correlation matrix
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_P.tex",argv[3], label);
    f=open_file(namefile,"w+");
    fprintf(f,"\\begin{gather}\n");
    fprintf(f,"\\chi^2/d.o.f.=%g \\\\ \n", fit_out.chi2[Njack-1]);//error_jackboot(argv[1], Njack, fit_out.chi2)
    for (int i=0;i<Npar;i++){
        fprintf(f,"P[%d]=%g\\pm (%.2g) \\\\ \n",i, fit_out.P[i][Njack-1], error_jackboot(argv[1], Njack, fit_out.P[i]));
    }
    fprintf(f,"\\end{gather}\n");
    double **cov=covariance(argv[1], Npar, Njack, fit_out.P);
    fprintf(f,"{\\tiny\\begin{gather}\n C=\\begin{pmatrix}\n");
    for (int i=0;i<fit_info.Npar;i++){
        for (int k=0;k<i;k++)
            cov[i][k]/=sqrt(cov[i][i]*cov[k][k]);
        for (int k=i;k<fit_info.Npar;k++)
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
    
    double **tif=swap_indices(fit_info.Npar,Njack,fit_out.P);
    
    /////////fit band L
    //    print_fit_band_L( argv, gjack , fit_info,  label,   fit_out, en, x,  y,   params,  myen);
     
    
    free_2(Njack,tif);
    
}
    
    
struct fit_result fit_data(char **argv, vector<header_BSM> params ,vector<data_BSM> gjack, double lhs_fun(int, int, int ,vector<header_BSM> ,vector<data_BSM>,struct fit_type ) , struct fit_type fit_info, const char* label, std::vector<int> myen ){
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
    
    printf("///// fit name:  %s \n",label);
    ////// allocation end
    /////////////// init
    std::mt19937 mt_rand(123);
    if (fit_info.guess.size()==0){
        for (int i=0;i<Npar;i++)
            guess[i]=mt_rand()/((double)mt_rand.max() );//rand()/((double)RAND_MAX);
    }
    else{
        for (int i=0;i<fit_info.guess.size();i++)
            guess[i]=fit_info.guess[i];
        for (int i=fit_info.guess.size();i<Npar;i++)
            guess[i]=1;
    }
    
    
    
    //init x
    for (int j=0;j<Njack;j++){
        int count=0;
        for (int n=0;n<N;n++){
            for (int e=0;e<en[n];e++){
                x[j][count][0]=(double) params[myen[e]].L;//L
                x[j][count][1]=(double) params[myen[e]].T;//T
                x[j][count][2]=(double) params[myen[e]].rho;
                x[j][count][3]=(double) params[myen[e]].eta;
                x[j][count][4]=(double) params[myen[e]].csw;
                x[j][count][5]=(double) params[myen[e]].mu03;
                x[j][count][6]=(double) params[myen[e]].m0;
                
                
                
                for(int i=0 ; i< fit_info.n_ext_P; i++){
                    x[j][count][i+fit_info.Nvar]=fit_info.ext_P[i][j];
                }
                
                count++;
            }
        }
    }
    
    ////////////////////////////////////////// y
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
    if (fit_info.guess.size()==0){
        guess=guess_for_non_linear_fit_Nf(N, en,x[Njack-1], y[Njack-1] , Nvar,  Npar, fit_info.function,guess );
    }
    
    for (int j=Njack-1;j>=0;j--){
        
        fit[j]=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function, guess, fit_info);
        
        
        fit_out.chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],fit[j] , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        
                    
        
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
