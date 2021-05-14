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
#include "mass_phi4.hpp"
using namespace std;


///////////////////////////////////////////////////////////////////////////////////////////
///////kcotd
//////////////////////////////////////////////////////////////////////////////////////////
inline double kcotd(double E2,double mass, int  *dvec,int L  ){
    double E2_CM=energy_CM(E2,dvec,L);
    double k=sqrt(E2_CM*E2_CM/4. -mass*mass);
    double delta;
    double kcotd;
    
    delta=phase_shift( E2, mass,dvec, L );
    kcotd=k/std::tan(delta);
    return kcotd;
}

double lhs_kcotd(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack, struct fit_type fit_info ){
    double r;
    if(n==0){//E2_0
        int dvec[3]= {0,0,0};
        r=kcotd( gjack[e].jack[4][j] ,gjack[e].jack[1][j], dvec, params[e].data.L[1]  );
    }
    if(n==1){//E2_0_p1
        int dvec[3]= {1,0,0};
        r=kcotd( gjack[e].jack[100][j] ,gjack[e].jack[1][j], dvec, params[e].data.L[1]  );
    }
    if(n==2){//E2_0_p1
        int dvec[3]= {1,1,0};
        r=kcotd( gjack[e].jack[102][j] ,gjack[e].jack[1][j], dvec, params[e].data.L[1]  );
    }
    if(n==3){//E2_0_p1
        int dvec[3]= {1,1,1};
        r=kcotd( gjack[e].jack[104][j] ,gjack[e].jack[1][j], dvec, params[e].data.L[1]  );
    }
    if(n==4){//E2_0_A1
        int dvec[3]= {0,0,0};
        r=kcotd( gjack[e].jack[80][j] ,gjack[e].jack[1][j], dvec, params[e].data.L[1]  );
    }
    
    return r;
    
}
double compute_k(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack, struct fit_type fit_info ){
    double L=params[e].data.L[1];
    double mass=gjack[e].jack[1][j];
    double E2=gjack[e].jack[4][j];
    double E2_0_p1=gjack[e].jack[100][j];
    double E2_0_p11=gjack[e].jack[102][j];
    double E2_0_p111=gjack[e].jack[104][j];
    double E2_0_A1=gjack[e].jack[80][j];
    double E2_CM;
    if(n==0){//E2_0
        int dvec[3]= {0,0,0};
        E2_CM=energy_CM(E2,dvec,L);
    }
    else if(n==1){//E2_0_p1
        int dvec[3]= {1,0,0};
        E2_CM=energy_CM(E2_0_p1,dvec,L);
    }
    else if(n==2){//E2_0_p11
        int dvec[3]= {1,1,0};
        E2_CM=energy_CM(E2_0_p11,dvec,L);
    }
    else if(n==3){//E2_0_p111
        int dvec[3]= {1,1,1};
        E2_CM=energy_CM(E2_0_p111,dvec,L);
    }
    else if(n==4){//E2_0_A1
        int dvec[3]= {0,0,0};
        E2_CM=energy_CM(E2_0_A1,dvec,L);
    }
    else {
        exit(1);
    }
    return sqrt(E2_CM*E2_CM/4. -mass*mass);
    
}
double rhs_kcotd(int n, int Nvar, double *x,int Npar,double  *P){
    double a0=P[0], r0=P[1], P2=P[2];
    double k=x[6];
    

    return 1.0/a0  + r0*k*k/2. - P2*r0*r0*r0*k*k*k*k;

}
void init_dvec(int n,int *dvec,int *dvec1, int *dvec2){
    
    if(n==0){//E2_0
        dvec[0]=0; dvec[1]=0; dvec[2]=0;
        dvec1[0]=0; dvec1[1]=0; dvec1[2]=0;
        dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
    }
    else if(n==1){//E2_0_p1
        dvec[0]=1; dvec[1]=0; dvec[2]=0;
        dvec1[0]=1; dvec1[1]=0; dvec1[2]=0;
        dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
    }
    else if(n==2){//E2_0_p11
        dvec[0]=1; dvec[1]=1; dvec[2]=0;
        dvec1[0]=1; dvec1[1]=1; dvec1[2]=0;
        dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
    }
    else if(n==3){//E2_0_p111
        dvec[0]=1; dvec[1]=1; dvec[2]=1;
        dvec1[0]=1; dvec1[1]=1; dvec1[2]=1;
        dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
    }
    else if(n==4){//E2_0_A1
        dvec[0]=0; dvec[1]=0; dvec[2]=0;
        dvec1[0]=1; dvec1[1]=0; dvec1[2]=0;
        dvec2[0]=-1; dvec2[1]=0; dvec2[2]=0;
    }
    else {
        exit(1);
    }
}


double to_invert_k_from_phase_shift(int n, int Nvar, double *x,int Npar,double  *P){
    double a0=P[0], r0=P[1], P2=P[2];
    int L=x[0];
    double mass=x[1];
    double k=x[2];
    int dvec[3],dvec1[3],dvec2[3];
    
    init_dvec(n,dvec,dvec1,dvec2);
    
    double E2_CM=sqrt(k*k+mass*mass)*2.;
    double E2=sqrt( (k*k+mass*mass)*4.+ (2*pi_greco/L)*(2*pi_greco/L)*( dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2]   )  );
    
    double qsq=k*k * (L/(2*pi_greco))*  (L/(2*pi_greco));
    double gamma=E2/E2_CM;
    double A=1.;
    double z[2];
//     printf("a0=%.12g   r0=%.12g    P2=%.12g  dvec=(%d,%d,%d)\n",a0,r0,P2,dvec[0] ,dvec[1] ,dvec[2]);
//     printf("zeta input:  qsq=%g   dvec=(%d,%d,%d) gamma=%g\n ",qsq,dvec[0] ,dvec[1] ,dvec[2], gamma);
    dzeta_function(z,  qsq,0 , 0, dvec, gamma, A, 1.e-3, 1.e6 ,2);
    std::complex<double>  zc(z[0],z[1]);
//     std::complex<double>   q=std::sqrt(std::complex<double> (qsq,0));
    
    double r=real(zc*2.*pi_greco/(pow(pi_greco,3./2.) *L*gamma )    );
    double kcotdelta=1.0/a0  + r0*k*k/2. - P2*r0*r0*r0*k*k*k*k;
//     printf("kcotdelta =%g  k=%g   r=%g  qsq=%g  gamma=%g  z=%g, %g\n",kcotdelta,k,r,qsq,gamma,z[0],z[1]);
//    std::cout<< "ZC " << zc << "       q "<< q<<endl;
//     printf(" k=%g   fit=%g        mass=%g    dvec=(%d,%d,%d) E2CM=%g\n ",k, r/kcotdelta  ,mass, dvec[0] ,dvec[1] ,dvec[2], E2_CM);
    return         r/kcotdelta;
    
    
}

double rhs_k_from_phase_shift(int n, int Nvar, double *x,int Npar,double  *P){
    double E2;
    double xx[3]={x[0],x[1],1} ; //{L,mass,   k to be find by the bisection}
    //double p[5]={0 ,2.*pi_greco/x[0], sqrt(2)*2.*pi_greco/x[0] , sqrt(3)*2.*pi_greco/x[0],0 };
    
//     printf("before   n=%d   L=%g   m=%g  a=%g  r=%g   P=%g\n",n,x[0],x[1],P[0],P[1],P[2]);
    int dvec[3],dvec1[3],dvec2[3];
    init_dvec(n,dvec,dvec1,dvec2);
    double mass=x[1];
    int L=x[0];
    // xmin have to be k in free theory
    double E1f=sqrt(mass*mass+(2*pi_greco/L)*(2*pi_greco/L)*(dvec1[0]*dvec1[0]+dvec1[1]*dvec1[1]+dvec1[2]*dvec1[2])   );
    double E2f=sqrt(mass*mass+(2*pi_greco/L)*(2*pi_greco/L)*(dvec2[0]*dvec2[0]+dvec2[1]*dvec2[1]+dvec2[2]*dvec2[2])   );
    double Ef=E1f+E2f;
    double ECMfsq=Ef*Ef-(2*pi_greco/L)*(2*pi_greco/L)*(dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
    double kf=sqrt(ECMfsq/4.-mass*mass);
    
    
    E1f=sqrt(mass*mass+(2*pi_greco/L)*(2*pi_greco/L)*((dvec1[0]+1)*(dvec1[0]+1)+dvec1[1]*dvec1[1]+dvec1[2]*dvec1[2])   );
    Ef=E1f+E2f;
    ECMfsq=Ef*Ef-(2*pi_greco/L)*(2*pi_greco/L)*((dvec[0]+1)*(dvec[0]+1)+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
    double kf1=sqrt(ECMfsq/4.-mass*mass);
    

    double xmin=kf+ (kf1-kf)/1e+3;
    double xmax=kf1- (kf1-kf)/1e+3;
//     printf("kf=%g   kf1=%g      x=%g   %g\n",kf,kf1,xmin,xmax);
     if (n==0)  { xmax=0.07;}//0.015
     if (n==1)  { xmax=0.1;}//0.7
     if (n==2)  { xmax=0.1;}//0.8
     if (n==3)  { xmax=0.1;}//0.8
     if (n==4)  { xmax=0.1;}//0.8
     
   /* for(int ik =0;ik<100;ik++){
        double k=0+0.0001*ik;
        xx[2]=k;
        double f=to_invert_k_from_phase_shift( n, 3, xx, Npar,P);
        printf("k=%g   fun(to be 1)=%g\n ",k,f);
        
    }
  */  
    
    E2=rtbis_func_eq_input( to_invert_k_from_phase_shift , n,  3, xx, Npar, P, 2/*ivar*/,1. /*input*/,   xmin,   xmax , 1e-3,  1);
    
//     printf("after  found root k=%.8g   n=%d\n",E2,n);
    return E2;
}

double lhs_k(int n, int e , int j , vector<cluster::IO_params> params,vector<data_phi> gjack, struct fit_type fit_info ){
    double E2;
    int dvec[3];
    double mass=gjack[e].jack[1][j];
    if(n==0){//E2_0
        E2=gjack[e].jack[4][j];
        dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else if(n==1){//E2_0_p1
        E2=gjack[e].jack[100][j];
        dvec[0]=1; dvec[1]=0; dvec[2]=0;
    }
    else if(n==2){//E2_0_p11
        E2=gjack[e].jack[102][j];
        dvec[0]=1; dvec[1]=1; dvec[2]=0;
    }
    else if(n==3){//E2_0_p111
        E2=gjack[e].jack[104][j];
        dvec[0]=1; dvec[1]=1; dvec[2]=1;
    }
    else if(n==4){//E2_0_A1
        E2=gjack[e].jack[80][j];
        dvec[0]=0; dvec[1]=0; dvec[2]=0;
    }
    else{exit(1); E2=0 ; dvec[0]=0; dvec[1]=0; dvec[2]=0;}
    double E2_CM=sqrt(E2*E2-(2.*pi_greco/params[e].data.L[1])*(2.*pi_greco/params[e].data.L[1])*( dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2]   ));
    double k=sqrt((E2_CM*E2_CM/4. - mass*mass));
    return k;
    
}




///////////////////////////////////////////////////////////////////////////////////////////
///////kcotd
//////////////////////////////////////////////////////////////////////////////////////////



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
            fprintf(f," %g   %g   %g   %d\t ",x[Njack-1][e+count][0], y[Njack-1][e+count][0], y[Njack-1][e+count][1] , params[myen[e]].data.L[0]);
            fprintf(f," %g   \n ",x[Njack-1][e+count][6]);
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
    
    double **tif=swap_indices(fit_info.Npar,Njack,fit_out.P);
    /////////fit band L
    
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_L.txt",argv[3], label);
    f=open_file(namefile,"w+");
    double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
    double *tmpy=(double*) malloc(sizeof(double*)* Njack);
    
    for (int i=0 ; i<100; i++){
                 
        for (int j=0;j<Njack;j++){
            tmpx[0]=10+i*0.5;
            tmpx[1]=gjack[0].jack[1][j];//m0
            tmpx[2]=gjack[0].jack[2][j];//m1
            tmpx[3]=gjack[0].jack[4][j];//E20
            tmpx[4]=gjack[0].jack[5][j];//E21
            tmpx[5]=(double) params[0].data.L[0];//T
            tmpx[6]=x[j][0][6];
            for(int i=fit_info.Nvar ; i<fit_info.Nvar+ fit_info.n_ext_P; i++)
                tmpx[i]=fit_info.ext_P[fit_info.Nvar][j];
            
            tmpy[j]=fit_info.function(0,Nvar,tmpx,Npar,tif[j]);//N, Nvar, x ,Npar,P
        }
        fprintf(f,"%g  \t %g  %g\n",tmpx[0],tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
        
    }
    free(tmpy);free(tmpx);
    //free_2(Njack,tif);
    fclose(f); 
    
    /////////fit band L for each n
    for (int n=0;n< N; n++){

        mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_n%d_L.txt",argv[3], label,n);
        f=open_file(namefile,"w+");
        double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
        double *tmpy=(double*) malloc(sizeof(double*)* Njack);
        
        for (int i=0 ; i<100; i++){
            
            for (int j=0;j<Njack;j++){
                tmpx[0]=10+i*0.5;
                tmpx[1]=gjack[myen.back()].jack[1][j];//m0   put for each n the mass of the last ensemble
                tmpx[2]=gjack[0].jack[2][j];//m1  //
                tmpx[3]=gjack[0].jack[4][j];//E20
                tmpx[4]=gjack[0].jack[5][j];//E21
                tmpx[5]=(double) params[0].data.L[0];//T
                tmpx[6]=x[j][0][6];
                for(int i=fit_info.Nvar ; i<fit_info.Nvar+ fit_info.n_ext_P; i++)
                    tmpx[i]=fit_info.ext_P[fit_info.Nvar][j];
                
                tmpy[j]=fit_info.function(n,Nvar,tmpx,Npar,tif[j]);//N, Nvar, x ,Npar,P
            }
            fprintf(f,"%g  \t %g  %g\n",tmpx[0],tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
            
        }
        free(tmpy);free(tmpx);
        fclose(f); 
        
    }
    /////////fit band T
    
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_T.txt",argv[3], label);
    f=open_file(namefile,"w+");
    
    tmpx=(double*) malloc(sizeof(double*)* Nvar);
    tmpy=(double*) malloc(sizeof(double*)* Njack);
    
    for (int i=0 ; i<100; i++){
        
        for (int j=0;j<Njack;j++){
            tmpx[0]=(double) params[0].data.L[1];//L
            tmpx[1]=gjack[0].jack[1][j];//m0
            tmpx[2]=gjack[0].jack[2][j];//m1
            tmpx[3]=gjack[0].jack[4][j];//E20
            tmpx[4]=gjack[0].jack[5][j];//E21
            tmpx[5]=10+i*1;//T
            tmpx[6]=x[j][0][6];
            for(int i=fit_info.Nvar ; i< fit_info.n_ext_P; i++)
                tmpx[i]=fit_info.ext_P[fit_info.Nvar][j];
            
            tmpy[j]=fit_info.function(0,Nvar,tmpx,Npar,tif[j]);//N, Nvar, x ,Npar,P
        }
        fprintf(f,"%g  \t %g  %g\n",tmpx[5],tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
       
    }
    free(tmpy);free(tmpx);
    fclose(f);  
    ////////// end fit band T
    /////////fit band T
    
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_k.txt",argv[3], label);
    f=open_file(namefile,"w+");
    tmpx=(double*) malloc(sizeof(double*)* Nvar);
    tmpy=(double*) malloc(sizeof(double*)* Njack);
    
    for (int i=0 ; i<100; i++){
         
        for (int j=0;j<Njack;j++){
            tmpx[0]=(double) params[0].data.L[1];//L
            tmpx[1]=gjack[0].jack[1][j];//m0
            tmpx[2]=gjack[0].jack[2][j];//m1
            tmpx[3]=gjack[0].jack[4][j];//E20
            tmpx[4]=gjack[0].jack[5][j];//E21
            tmpx[5]=params[0].data.L[1];//T
            tmpx[6]=0+i*0.004;//k
            for(int i=fit_info.Nvar ; i< fit_info.n_ext_P; i++)
                tmpx[i]=fit_info.ext_P[fit_info.Nvar][j];
            
            tmpy[j]=fit_info.function(0,Nvar,tmpx,Npar,tif[j]);//N, Nvar, x ,Npar,P
        }
        fprintf(f,"%g  \t %g  %g\n",tmpx[6],tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
        
    }
    free(tmpy);free(tmpx);
    fclose(f);  
    ////////// end fit band k
    free_2(Njack,tif);
    
    
}


struct fit_result fit_data(char **argv, vector<cluster::IO_params> params ,vector<data_phi> gjack, double lhs_fun(int, int, int ,vector<cluster::IO_params> ,vector<data_phi>,struct fit_type ) , struct fit_type fit_info, const char* label, std::vector<int> myen, std::vector<double> start_point={}){
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
    if (start_point.size()==0){
        for (int i=0;i<Npar;i++)
            guess[i]=-1;//rand()/((double)RAND_MAX);
    }
    else{
        for (int i=0;i<start_point.size();i++)
            guess[i]=start_point[i];
        for (int i=start_point.size();i<Npar;i++)
            guess[i]=1;
    }
    
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
                
                x[j][count][6]=compute_k(n,myen[e],j,params,gjack,fit_info);//k
    
                
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
    if (start_point.size()==0){
        guess=guess_for_non_linear_fit_Nf(N, en,x[Njack-1], y[Njack-1] , fit_info.Nvar,  Npar, fit_info.function,guess );
    }
    
    for (int j=0;j<Njack;j++){
        if (strcmp(label,"k_from_phase_shift")==0 ){
            printf("jack =%d\n",j);
            fit[j]=non_linear_fit_Nf(N, en,x[j], y[j] , fit_info.Nvar,  Npar, fit_info.function,guess,0.0001, 0.01 );
        }
        else 
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
