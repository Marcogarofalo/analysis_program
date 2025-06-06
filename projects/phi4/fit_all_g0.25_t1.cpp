#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <memory>
#include <cmath>
#include <vector>

#include "global.hpp"
#include "resampling.hpp"
#include "linear_fit.hpp"
#include "mutils.hpp"
#include "header_phi4.hpp"
#include "zeta_interpolation.hpp"
#include "fit_function.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"



// using namespace std;
int Ne=0;


void print_fit_band_L_M(char **argv,  std::vector<data_phi> gjack ,struct fit_type fit_info,struct fit_type fit_info_m0 , const char* label, struct fit_result fit_out, struct fit_result fit_out_m0,    std::vector<cluster::IO_params> params, std::vector<int> myen, std::vector<int> Lrange={16,50}){
    int Npar=fit_info.Npar;
    int Nvar=fit_info.Nvar+fit_info.n_ext_P;
    int Njack=gjack[0].Njack;
    int N=fit_info.N;
    char namefile[NAMESIZE];
    FILE *f;
    
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_k.txt",argv[3], label);
    f=open_file(namefile,"w+");
    double **tif=swap_indices(fit_info.Npar,Njack,fit_out.P);
    double **tif_m0=swap_indices(fit_info_m0.Npar,Njack,fit_out_m0.P);
    double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
    double *tmpy=(double*) malloc(sizeof(double*)* Njack);
    printf("writing: %s\n",namefile);
    
    for (int n=0;n< N; n++){
        
        mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_n%d_L.txt",argv[3], label,n);
        f=open_file(namefile,"w+");
        double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
        double *tmpy=(double*) malloc(sizeof(double*)* Njack);
        printf("writing: %s\n",namefile);
        
        for (int i=Lrange[0] ; i<Lrange[1]; i++){
            double finalL=i;
            tmpx[0]=finalL;
            double *E3_m=(double*) malloc(sizeof(double) *Njack);
            for (int j=0;j<Njack;j++){
                E3_m[j]=lhs_E3_m(n,0,j,params,gjack,fit_info);
            }
            double E3_m_err=error_jackboot(argv[1], Njack, E3_m );
            for (int j=0;j<Njack;j++){
                tmpx[1]=fit_info_m0.function(0,fit_info_m0.Nvar,tmpx,fit_info_m0.Npar,tif_m0[j]); //m0   put for each n the mass of the last ensemble
                
                tmpx[2]=gjack[0].jack[2][j];//m1  //
                tmpx[3]=gjack[0].jack[4][j];//E20
                tmpx[4]=gjack[0].jack[5][j];//E21
                tmpx[5]=(double) params[0].data.L[0];//T
                tmpx[6]=1;//k
                tmpx[7]=1;//MpiL
                
                double *x=(double*) malloc(sizeof(double)*myen.size());
                double *y=(double*) malloc(sizeof(double)*myen.size());
                for (int e=0; e<myen.size();e++){
                    y[e]= lhs_E3_m(n,myen[e],j,params,gjack,fit_info);// E3( \vec{n} )/mass
                    x[e]=params[myen[e]].data.L[1];
                }
                tmpx[8]=inter_spline( tmpx[0], myen.size(), x, y   ); // tmpx[0]=L 
                                     
                free(x);free(y);
                //tmpx[8]=E3_m[j];// E3( \vec{n} )/mass
                
                // as error on E3/m  we take the error on the ensemble=0
                tmpx[9]=E3_m_err;
                
                
                for(int i=fit_info.Nvar ; i<fit_info.Nvar+ fit_info.n_ext_P; i++)
                    tmpx[i]=fit_info.ext_P[i-fit_info.Nvar][j];
                
                
                
                tmpy[j]=fit_info.function(n,Nvar,tmpx,Npar,tif[j]);
//                 tmpy[j]=tmpx[8];
//                 printf("L=%g  j=%d  E=%g\n",finalL,j,tmpx[8]);
//                 if (fabs(finalL-36)<1e-5  && n==2)
//                     printf("%g   j=%d \t k=%g  m=%g P0=%g  P1=%g\n",finalL,j,tmpy[j], tmpx[1], tif[j][0] ,tif[j][1]  );
//                 
            }
            /*if (fabs(finalL-36)<1e-5  && n==2){
                printf("%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
            }*/
            fprintf(f,"%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
            
            free(E3_m);
        }
        
        free(tmpy);free(tmpx);
        fclose(f); 
    }
    free_2(Njack,tif);
    free_2(Njack,tif_m0);
    
}





void print_fit_band_E3_vs_L(char **argv, std::vector<data_phi> gjack , struct fit_type fit_info, struct fit_type fit_info_m0 , const char* label, struct fit_result fit_out, struct fit_result fit_out_m0,    std::vector<cluster::IO_params> params, std::vector<int> myen,  struct fit_type fit_info_E3_poly ,fit_result fit_E3_poly, std::vector<int> Lrange={16,50}){
    int Npar=fit_info.Npar;
    int Nvar=fit_info.Nvar+fit_info.n_ext_P;
    int Njack=gjack[0].Njack;
    int N=fit_info.N;
    char namefile[NAMESIZE];
    FILE *f;
    
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_k.txt",argv[3], label);
    f=open_file(namefile,"w+");
    double **tif=swap_indices(fit_info.Npar,Njack,fit_out.P);
    double **tif_m0=swap_indices(fit_info_m0.Npar,Njack,fit_out_m0.P);
    double **tif_E3_poly=swap_indices(fit_info_E3_poly.Npar, Njack, fit_E3_poly.P);
    double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
    double *tmpy=(double*) malloc(sizeof(double*)* Njack);
    printf("writing: %s\n",namefile);
    
    for (int n=0;n< N; n++){
        
        mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_n%d_L.txt",argv[3], label,n);
        f=open_file(namefile,"w+");
        double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
        double *tmpy=(double*) malloc(sizeof(double*)* Njack);
        printf("writing: %s\n",namefile);
        
        for (int i=Lrange[0] ; i<Lrange[1]; i++){
            double finalL=i;
            tmpx[0]=finalL;
            double *E3_m=(double*) malloc(sizeof(double) *Njack);
            for (int j=0;j<Njack;j++){
                E3_m[j]=fit_info_E3_poly.function(n,fit_info_E3_poly.Nvar , tmpx, fit_info_E3_poly.Npar, tif_E3_poly[j]);
            }
            double E3_m_err=error_jackboot(argv[1], Njack, E3_m );
            
            
            
            for (int j=0;j<Njack;j++){
                tmpx[1]=fit_info_m0.function(0,fit_info_m0.Nvar,tmpx,fit_info_m0.Npar,tif_m0[j]); //m0   put for each n the mass of the last ensemble
                
                tmpx[2]=gjack[0].jack[2][j];//m1  //
                tmpx[3]=gjack[0].jack[4][j];//E20
                tmpx[4]=gjack[0].jack[5][j];//E21
                tmpx[5]=(double) params[0].data.L[0];//T
                tmpx[6]=1;//k
                tmpx[7]=1;//MpiL
                
//                 double *x=(double*) malloc(sizeof(double)*myen.size());
//                 double *y=(double*) malloc(sizeof(double)*myen.size());
//                 for (int e=0; e<myen.size();e++){
//                     y[e]= lhs_E3_m(n,myen[e],j,params,gjack,fit_info);// E3( \vec{n} )/mass
//                     x[e]=params[myen[e]].data.L[1];
//                 }
//                 tmpx[8]=inter_spline( tmpx[0], myen.size(), x, y   ); // tmpx[0]=L 
//                 free(x);free(y);
                
                tmpx[8]=E3_m[Njack-1];// E3( \vec{n} )/mass
                // as error on E3/m  we take the error on the ensemble=0
                tmpx[9]=E3_m_err;
                
                tmpx[11]=E3_m[Njack-1];// E3( \vec{n} )/mass
                // as error on E3/m  we take the error on the ensemble=0
                tmpx[12]=E3_m_err;
                
                
                
                for(int i=fit_info.Nvar ; i<fit_info.Nvar+ fit_info.n_ext_P; i++)
                    tmpx[i]=fit_info.ext_P[i-fit_info.Nvar][j];
                
                
                
                tmpy[j]=fit_info.function(n,Nvar,tmpx,Npar,tif[j]);
                //                 tmpy[j]=tmpx[8];
                //                 printf("L=%g  j=%d  E=%g\n",finalL,j,tmpx[8]);
                //                 if (fabs(finalL-36)<1e-5  && n==2)
                //                     printf("%g   j=%d \t k=%g  m=%g P0=%g  P1=%g\n",finalL,j,tmpy[j], tmpx[1], tif[j][0] ,tif[j][1]  );
                //                 
            }
            /*if (fabs(finalL-36)<1e-5  && n==2){
             *                printf("%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
        }*/
            fprintf(f,"%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
            
            free(E3_m);
        }
        
        free(tmpy);free(tmpx);
        fclose(f); 
    }
    free_2(Njack,tif);
    free_2(Njack,tif_m0);
    free_2(Njack,tif_E3_poly);
}





void print_kiso_P0_inf_L_M(char **argv, std::vector<data_phi> gjack , struct fit_type fit_info, struct fit_type fit_info_m0 , const char* label, struct fit_result fit_out, struct fit_result fit_out_m0,    std::vector<cluster::IO_params> params, std::vector<int> myen,  struct fit_type fit_info_E3_poly ,fit_result fit_E3_poly, std::vector<int> Lrange={16,50}){
    int Npar=fit_info.Npar;
    int Nvar=fit_info.Nvar+fit_info.n_ext_P;
    int Njack=gjack[0].Njack;
    int N=fit_info.N;
    char namefile[NAMESIZE];
    FILE *f;
    
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_k.txt",argv[3], label);
    f=open_file(namefile,"w+");
    double **tif=double_malloc_2( Njack,Npar);
    for (int i=0;i<Npar; i++)
        for (int j=0;j<Njack;j++)
            tif[j][i]=0;
    double **tif_m0=swap_indices(fit_info_m0.Npar,Njack,fit_out_m0.P);
    double **tif_E3_poly=swap_indices(fit_info_E3_poly.Npar, Njack, fit_E3_poly.P);
    double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
    double *tmpy=(double*) malloc(sizeof(double*)* Njack);
    printf("writing: %s\n",namefile);
    
    for (int n=0;n< N; n++){
        
        mysprintf(namefile,NAMESIZE,"%s/kiso_P0_n%d_L.txt",argv[3],n);
        f=open_file(namefile,"w+");
        double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
        double *tmpy=(double*) malloc(sizeof(double*)* Njack);
        printf("writing: %s\n",namefile);
        
        for (int i=Lrange[0] ; i<Lrange[1]; i++){
            double finalL=i;
            tmpx[0]=finalL;
            double *E3_m=(double*) malloc(sizeof(double) *Njack);
            for (int j=0;j<Njack;j++){
                E3_m[j]=fit_info_E3_poly.function(n,fit_info_E3_poly.Nvar , tmpx, fit_info_E3_poly.Npar, tif_E3_poly[j]);
            }
            double E3_m_err=error_jackboot(argv[1], Njack, E3_m );
            for (int j=0;j<Njack;j++){
                tmpx[1]=fit_info_m0.function(0,fit_info_m0.Nvar,tmpx,fit_info_m0.Npar,tif_m0[j]); //m0   put for each n the mass of the last ensemble
                
                tmpx[2]=gjack[0].jack[2][j];//m1  //
                tmpx[3]=gjack[0].jack[4][j];//E20
                tmpx[4]=gjack[0].jack[5][j];//E21
                tmpx[5]=(double) params[0].data.L[0];//T
                tmpx[6]=1;//k
                tmpx[7]=1;//MpiL
                
                    
                tmpx[8]=E3_m[Njack-1];// E3( \vec{n} )/mass
                
                tmpx[9]=E3_m_err;
                
                
                for(int i=fit_info.Nvar ; i<fit_info.Nvar+ fit_info.n_ext_P; i++)
                    tmpx[i]=fit_info.ext_P[i-fit_info.Nvar][j];
                
                
                
                tmpy[j]=fit_info.function(n,Nvar,tmpx,Npar,tif[j]);
                              
            }
            
            fprintf(f,"%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
            
            free(E3_m);
        }
        
        free(tmpy);free(tmpx);
        fclose(f); 
    }
    for (int i=0;i<Npar; i++)
        for (int j=0;j<Njack;j++)
            tif[j][i]=-1e+3;
    for (int n=0;n< N; n++){
        
        mysprintf(namefile,NAMESIZE,"%s/kiso_P-1e+3_n%d_L.txt",argv[3],n);
        f=open_file(namefile,"w+");
        double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
        double *tmpy=(double*) malloc(sizeof(double*)* Njack);
        printf("writing: %s\n",namefile);
        
        for (int i=Lrange[0] ; i<Lrange[1]; i++){
            double finalL=i;
            tmpx[0]=finalL;
            double *E3_m=(double*) malloc(sizeof(double) *Njack);
            for (int j=0;j<Njack;j++){
                E3_m[j]=fit_info_E3_poly.function(n,fit_info_E3_poly.Nvar , tmpx, fit_info_E3_poly.Npar, tif_E3_poly[j]);
            }
            double E3_m_err=error_jackboot(argv[1], Njack, E3_m );
            for (int j=0;j<Njack;j++){
                tmpx[1]=fit_info_m0.function(0,fit_info_m0.Nvar,tmpx,fit_info_m0.Npar,tif_m0[j]); //m0   put for each n the mass of the last ensemble
                
                tmpx[2]=gjack[0].jack[2][j];//m1  //
                tmpx[3]=gjack[0].jack[4][j];//E20
                tmpx[4]=gjack[0].jack[5][j];//E21
                tmpx[5]=(double) params[0].data.L[0];//T
                tmpx[6]=1;//k
                tmpx[7]=1;//MpiL
                
                
                tmpx[8]=E3_m[Njack-1];// E3( \vec{n} )/mass
                
                tmpx[9]=E3_m_err;
                
                
                for(int i=fit_info.Nvar ; i<fit_info.Nvar+ fit_info.n_ext_P; i++)
                    tmpx[i]=fit_info.ext_P[i-fit_info.Nvar][j];
                
                
                
                tmpy[j]=fit_info.function(n,Nvar,tmpx,Npar,tif[j]);
                
            }
            
            fprintf(f,"%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
            
            free(E3_m);
        }
        
        free(tmpy);free(tmpx);
        fclose(f); 
    }
    
    free_2(Njack,tif);
    free_2(Njack,tif_m0);
    free_2(Njack,tif_E3_poly);
}

void print_phase_shift(char **argv, std::vector<data_phi> gjack , struct fit_type fit_info , const char* label, struct fit_result fit_out){
    
    char namefile[NAMESIZE];
    mysprintf(namefile,NAMESIZE,"%s/%s_fit_phase_shift.txt",argv[3], label);
    FILE *f=open_file(namefile,"w+");
    int Npar=fit_info.Npar;
    int Njack=gjack[0].Njack;
    
    int ikmax=100;
    double dk=0.03;
    double *delta=(double* ) malloc(sizeof(double)*Njack);
    for (int ik=0;ik<ikmax;ik++){
    
        double k_m=ik*dk;
        for (int j=0; j<Njack;j++){
        
            double a0m0=fit_out.P[0][j];
            double r0m0=fit_out.P[1][j];
            
            double kcotdelta_m=1.0/a0m0+ + r0m0*k_m*k_m/2.;  //   (k cot(d) )/ mass
            if (Npar>=3){
                kcotdelta_m+=fit_out.P[2][j]*r0m0*r0m0*r0m0*k_m*k_m*k_m*k_m;
            }
            delta[j]=std::atan( k_m/kcotdelta_m);
            
        }
        fprintf(f,"%g    %g   %g\n",k_m,delta[Njack-1], error_jackboot( argv[1] , Njack, delta) );
    }
    free(delta);
    fclose(f);
}


int main(int argc, char **argv){
     error(argc!=4,1,"main ",
         "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");   
    
     //int Ne=0;   
     /*cluster::IO_params *params=(cluster::IO_params*) malloc(sizeof(cluster::IO_params)*Ne);
     
     
     data_phi data;
     data_phi *dataj=(data_phi*) malloc(sizeof(data_phi*)*Ne);
     
     
     char namefile[NAMESIZE];
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L32_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     FILE *f=open_file(namefile,"r");
     read_header_phi4( f  , params[0]);
     read_dataj(f,params[0],dataj[0] );
     fclose(f);
     */
    
     std::vector<cluster::IO_params> paramsj;
     std::vector<data_phi> dataj;
     
     int Ne=0; 
     cluster::IO_params params;
     
     char namefile[NAMESIZE];
     char jackboot[NAMESIZE];
     mysprintf(jackboot,NAMESIZE,"%s" ,argv[1] );
     
     //mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L16_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     //emplace_back_par_data(namefile,paramsj,dataj);
    
     //0
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L24_msq0-4.864000_msq1-4.705000_l02.500000_l12.500000_mu5.000000_g0.250000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //1
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L26_msq0-4.864000_msq1-4.705000_l02.500000_l12.500000_mu5.000000_g0.250000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //2
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L28_msq0-4.864000_msq1-4.705000_l02.500000_l12.500000_mu5.000000_g0.250000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //3
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L30_msq0-4.864000_msq1-4.705000_l02.500000_l12.500000_mu5.000000_g0.250000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     
     //4
    //  mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L32_msq0-4.864000_msq1-4.710000_l02.500000_l12.500000_mu5.000000_g0.250000_rep0",argv[2],argv[1]);
    //  emplace_back_par_data(namefile,paramsj,dataj);

    //  mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L34_msq0-4.864000_msq1-4.710000_l02.500000_l12.500000_mu5.000000_g0.250000_rep0",argv[2],argv[1]);
    //  emplace_back_par_data(namefile,paramsj,dataj);
     
    //  mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L34_msq0-4.864000_msq1-4.710000_l02.500000_l12.500000_mu5.000000_g0.250000_rep0",argv[2],argv[1]);
    //  emplace_back_par_data(namefile,paramsj,dataj);
     
    //  mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L36_msq0-4.864000_msq1-4.710000_l02.500000_l12.500000_mu5.000000_g0.250000_rep0",argv[2],argv[1]);
    //  emplace_back_par_data(namefile,paramsj,dataj);
     
    //  mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L40_msq0-4.864000_msq1-4.710000_l02.500000_l12.500000_mu5.000000_g0.250000_rep0",argv[2],argv[1]);
    //  emplace_back_par_data(namefile,paramsj,dataj);
     
     
     
     printf("E1_0 =%f   %f\n", dataj[0].jack[1][dataj[0].Njack-1], error_jackboot(argv[1],dataj[0].Njack, dataj[0].jack[1]  ) );
    
     
     std::vector<data_phi> gjack= create_generalised_resampling_phi( dataj, paramsj );
     printf("GEVP_E2_01 =%f   %f\n", gjack[1].jack[19][gjack[1].Njack-1], error_jackboot(argv[1],gjack[1].Njack, gjack[1].jack[19]  ) );
     Ne=gjack.size();
     printf("number of ensembles = %d\n",Ne);
     
     std::vector<int> myen(Ne);
     for(int i=0;i< Ne; i++)  myen[i]=i;
    
     int Njack=gjack[0].Njack;
     /*std::vector<double> masses(Ne);
     printf("masses:");
     for(int e=0;e< Ne; e++){ masses[e]=gjack[e].jack[1][Njack-1];
        printf(" %g\t",masses[e]);
     }
     printf(" \n");
     std::sort(masses.begin(),masses.end());
     for(int e=0;e< Ne; e++){ 
         printf(" %g\t",masses[e]);
     }
     printf(" \n");
     
     std::vector<int> Ls(Ne);
     for(int e=0;e< Ne; e++) Ls[e]= paramsj[e].data.L[1];
     std::vector< std::vector<int> > momenta(5,std::vector<int>(3));
     
     momenta[0][0]=0; momenta[0][1]=0; momenta[0][2]=0; 
     momenta[1][0]=1; momenta[1][1]=0; momenta[1][2]=0; 
     momenta[2][0]=1; momenta[2][1]=1; momenta[2][2]=0; 
     momenta[3][0]=0; momenta[3][1]=0; momenta[3][2]=0; 
     momenta[4][0]=1; momenta[4][1]=1; momenta[4][2]=1; 
     */
//      std::vector<int> Ls={16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56};
     std::vector<int> Ls;//={16,17,18,19,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56};
     for (int  i=20;i<42;i++){Ls.emplace_back(i);}
//      std::vector<int> Ls={28,30,32,36};
     std::vector<double> masses;
     std::vector<double> err_mass;
     int e1;
     for (int L: Ls){
         e1=Ne-1;
         for(int e=0; e< Ne;e++){
             if (L==paramsj[e].data.L[1]){
                 e1=e;
                 printf("L=%d found   mass=%g\n",L,gjack[e1].jack[443][Njack-1]);
             }
         }
         masses.emplace_back(gjack[e1].jack[443][Njack-1]);// 443 ->GEVP_0_3_1_001_011_meffl0
         err_mass.emplace_back(  error_jackboot(argv[1],Njack, gjack[e1].jack[443] )  );
     }
     
    // zeta.Init_Lmq(Ls, masses, err_mass  );
    // zeta.write("zeta_g0.25.dat");
    zeta.read("zeta_g0.25.dat");
    // //   zeta_interpolation   zz;
    // //  zz.Init(argv[1] , myen , paramsj ,gjack  );
    // //  for(int e=0; e< Ne;e++){
    // //      printf("%g  %g\t",zeta.compute(28,0 ,12 *2* pi_greco/28., 0.5), zz.compute(28,0 ,12 *2* pi_greco/28., 0.5));
    // //      printf("%g  %g\t",zeta.compute(28,1 ,12 *2* pi_greco/28., 0.5), zz.compute(28,1 ,12 *2* pi_greco/28., 0.5));
    // //      printf("%g  %g\t",zeta.compute(28,2 ,12 *2* pi_greco/28., 0.5), zz.compute(28,2 ,12 *2* pi_greco/28., 0.5));
    // //      printf("%g  %g  %d %d\n",zeta.mass(e,0) , zz.mass(e,0),zeta.Ls(e),zz.Ls(e) );
    // //      printf(" %d %d %d %d %d %d\n",zeta.moms(0,0) , zz.moms(0,0),zeta.moms(0,1) , zz.moms(0,1),zeta.moms(0,1) , zz.moms(0,1) );
    // //      printf(" %d %d %d %d %d %d\n",zeta.moms(1,0) , zz.moms(1,0),zeta.moms(1,1) , zz.moms(1,1),zeta.moms(1,1) , zz.moms(1,1) );
    // // }
   
    //  for(int e=0+12; e< Ne+20;e++){
    //      for(int im=0;im<10;im++){
    //          printf("e=%d  L=%d im=%d  krange=[%g,%g] \n",e,zeta.Ls(e),im,
    //                 sqrt(zeta.kmax(e,2,im))*(2.*pi_greco)/zeta.Ls(e),
    //                 sqrt(zeta.kmin(e,2,im))*(2.*pi_greco)/zeta.Ls(e)
    //                 );
    //     }
    //  }
//      zeta_qsqg.Init(  );
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     //////////////////////////////////////////////////////////////////////////////////////////////////
     struct fit_type fit_info, fit_info_m0;
     
     struct fit_result  fit_m1, fit_m0;
     fit_info.Nvar=14;                   fit_info_m0.Nvar=14;            
     fit_info.Npar=2;                   fit_info_m0.Npar=2;           
     fit_info.N=1;                      fit_info_m0.N=1;
     fit_info.Njack=gjack[0].Njack;     fit_info_m0.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;                fit_info_m0.n_ext_P=0;
     fit_info.function=M_finite_volume; fit_info_m0.function=M_finite_volume;

     
     fit_m0=fit_data(argv,  paramsj ,gjack, M0_finite_volume_lhs ,fit_info, "M0_finite_vol_m1-4.705" ,myen);
    
     
     fit_m1=fit_data(argv,  paramsj ,gjack, M1_finite_volume_lhs ,fit_info, "M1_finite_vol_m1-4.705" ,myen);
     
     // corr 75 =E1_1_px
     fit_result fit_m1_p1=fit_data(argv,  paramsj ,gjack, M1_p_finite_volume_lhs<75> ,fit_info, "M1_p1_finite_vol" ,myen);
     free_fit_result(fit_info,fit_m1_p1);
     printf("\n/////////////////////////////////     E2_0//////////////////\n");
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // E20
     //////////////////////////////////////////////////////////////////////////////////////////////////
    
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     
     fit_info.Npar=1;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=muDE_rhs;
     
     struct fit_result fit_a_00=fit_data(argv,  paramsj ,gjack, muDE_00_lhs ,fit_info, "a_00_luscher_m1-4.705",myen );
     
     printf("\n/////////////////////////////////     k cot delta    //////////////////\n");
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // kcot
     //////////////////////////////////////////////////////////////////////////////////////////////////
     
    
     fit_info.Npar=3;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=rhs_kcotd_g;
     
     struct fit_result fit_kcotd=fit_data(argv,  paramsj ,gjack, lhs_kcotd_g ,fit_info, "kcotd_m1-4.705",myen );
     
      
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     printf("\n/////////////////////////////////   fit  k  form from_phase_shift    //////////////////\n");
     //////////////////////////////////////////////////////////////////////////////////////////////////
     
     
    //  fit_info.Npar=2;
    //  fit_info.N=1;
    //  fit_info.Njack=gjack[0].Njack;
    //  fit_info.n_ext_P=0;
    //  //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
    //  fit_info.function=rhs_k_from_phase_shift_g;
    //  fit_info.lambda=0.001;
    //  fit_info.acc=0.01;
    //  fit_info.h=1e-3;
    //  fit_info.Prange={100,100000};
    //  fit_info.devorder=2;
    //  fit_info.guess= {-0.121902,-100};
     
    //  struct fit_result k_from_phase_shift=fit_data(argv,  paramsj ,gjack, lhs_k_g ,fit_info, "k_from_phase_shift_m1-4.705",myen  );// {-0.948817,-114.788,0.0003987}
    //  print_fit_band_L_M( argv, gjack , fit_info,fit_info_m0 ,  "k_from_phase_shift_m1-4.705",   k_from_phase_shift ,fit_m0,    paramsj,  myen, {22,40});
    //  fit_info.guess=std::vector<double>();
     
     
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     printf("\n/////////////////////////////////   fit  deltaE2_m_quant_cond  //////////////////\n");
     //////////////////////////////////////////////////////////////////////////////////////////////////
    //  fit_info.restore_default();
     
    //  fit_info.Npar=2;
    //  fit_info.N=1;
    //  fit_info.Njack=gjack[0].Njack;
    //  fit_info.n_ext_P=0;
    //  //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
    //  fit_info.function=rhs_deltaE2_m_quant_cond_g;
     
    //  fit_info.lambda=0.001;
    //  fit_info.acc=0.01;
    //  fit_info.h=1e-3;
    //  fit_info.Prange={100,100000};
    //  fit_info.devorder=2;
    //  fit_info.guess={-4069.,-90} ;
     
    //  //      the zeta is computed analytically, use the interpolated one for faster result!!!!!!!!!
    //  //      struct fit_result k_from_phase_shift_3par=fit_data(argv,  paramsj ,gjack, lhs_k ,fit_info, "k_from_phase_shift_n5_3par",myen ,  {-0.11,-950, 6.4e-6} );// {-0.948817,-114.788,0.0003987}
    //  struct fit_result deltaE2_m_quant_cond=fit_data(argv,  paramsj ,gjack, lhs_deltaE2_m_latt_g ,fit_info, "deltaE2_m_quant_cond_m1-4.705",myen );
    //  print_fit_band_L_M( argv, gjack , fit_info,fit_info_m0 ,  "deltaE2_m_quant_cond_m1-4.705",   deltaE2_m_quant_cond ,fit_m0,    paramsj,  myen);
     
    //  print_phase_shift(argv, gjack ,  fit_info , "deltaE2_m_quant_cond_m1-4.705", deltaE2_m_quant_cond);
    //  fit_info.restore_default();
     

    printf("\n/////////////////////////////////     delta 2par   //////////////////\n");
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // kcot
     //////////////////////////////////////////////////////////////////////////////////////////////////
     
     
     fit_info.Npar=2;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=rhs_delta;
     
     struct fit_result fit_delta_2par=fit_data(argv,  paramsj ,gjack, lhs_delta_g ,fit_info, "delta_2par_m1-4.705",myen );
     free_fit_result(fit_info,fit_delta_2par);

        #ifdef not_compile


     ///////////////////////////////////////////////////////////////////////////////////////////////////
     printf("\n/////////////////////////////////   fit  E3 quant cond  //////////////////\n");
     //////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef PYTHON
     //// we need python
     wchar_t *program = Py_DecodeLocale(argv[0], NULL);
     if (program == NULL) {
         fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
         exit(1);
     }  
     Py_SetProgramName(program);  /* optional but recommended */
     
     Py_Initialize();
     ///////////// end python init
     
     
     printf("//////////////////// poly fit E3   ////////////////////////////////////\n");
     fit_type fit_info_E3_poly;
     fit_info_E3_poly.N=6;
     
     fit_info_E3_poly.Njack=gjack[0].Njack;
     fit_info_E3_poly.n_ext_P=2;
     fit_info_E3_poly.ext_P=(double**) malloc(sizeof(double*)*fit_info_E3_poly.n_ext_P);
     
     fit_info_E3_poly.ext_P[0]=deltaE2_m_quant_cond.P[0];
     fit_info_E3_poly.ext_P[1]=deltaE2_m_quant_cond.P[1];
     fit_info_E3_poly.Npar=fit_info_E3_poly.N*4;
     fit_info_E3_poly.function=rhs_poly_order_E3_m<4>;
     
     fit_info_E3_poly.lambda=0.001;
     fit_info_E3_poly.acc=0.01;
     fit_info_E3_poly.h=1e-3;
     fit_info_E3_poly.Prange={1000,10000};
     fit_info_E3_poly.devorder=2;
     
     mysprintf(namefile,NAMESIZE,"poly_E3andE1_N%d",fit_info_E3_poly.N );
    // struct fit_result fit_QC3_poly=fit_data(argv,  paramsj ,gjack, lhs_E3_m ,fit_info_E3_poly, namefile,myen   );
     
     struct fit_result fit_QC3_poly=fit_data(argv,  paramsj ,gjack, lhs_E3orE1_m ,fit_info_E3_poly, namefile,myen   );
     print_fit_band_L_M( argv, gjack , fit_info_E3_poly,fit_info_m0 ,  namefile,   fit_QC3_poly ,fit_m0,    paramsj,  myen, {23,41});
//      free_fit_result(fit_info,fit_QC3_poly);
//      fit_info_E3_poly.restore_default();
     
     
     printf("////////////////////  kiso pole fit   ////////////////////////////////////\n");
     init_python_detQC();
      init_python_detQC_kcot_kiso("kcot_2par", "kiso_pole");
//     init_python_detQC_kcot_kiso("kcot_2par", "kiso_2par");
//      init_python_detQC_kcot_kiso("kcot_2par", "kiso_1par");
     fit_info.Npar=2;
     fit_info.N=2;
     fit_info.Njack=gjack[0].Njack;
     
     
     
     fit_info.function=rhs_E3_m_QC3_pole;
     
     fit_info.n_ext_P=2;
     fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.ext_P[0]=deltaE2_m_quant_cond.P[0];
     fit_info.ext_P[1]=deltaE2_m_quant_cond.P[1];
     
     
     
     fit_info.lambda=0.001;
     fit_info.acc=0.01;
     fit_info.h=1e-3;
     //fit_info.Prange={1000,10000};
     fit_info.devorder=2;
     fit_info.verbosity=1;
     //fit_info.repeat_start=2;
     fit_info.guess={-22533.5, 80.5154};
     mysprintf(namefile,NAMESIZE,"QC3_N%d_%dpar_pole",fit_info.N, fit_info.Npar);
     struct fit_result fit_QC3_1par=fit_data(argv,  paramsj ,gjack, lhs_E3orE1_m ,fit_info, namefile,myen   );
     print_fit_band_E3_vs_L( argv, gjack , fit_info,fit_info_m0 ,  namefile,   fit_QC3_1par ,fit_m0,    paramsj,  myen,  fit_info_E3_poly, fit_QC3_poly, {23,41});
     
     fit_info.restore_default();
    
    
     ///// close python
     python_detQC_write_database();
     python_detQC_free();
     if (Py_FinalizeEx() < 0) {
         exit(120);
     }
     PyMem_RawFree(program);
     
#endif
     #endif
     return 0;
}
