#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"

#include "resampling.hpp"
#include "read.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"
#include <string>
#include "correlators_analysis.hpp"
#include "eigensystem.hpp"

#include <cstring> 
#include <string>
#include <fstream>
#include <memory>


//local folder
#include "header_BSM.hpp"
#include "lhs_functions.hpp"
#include "fit_data.hpp"



using namespace std;
int Ne=0;


void print_fit_band_eta(char **argv,vector<data_BSM> gjack ,struct fit_type fit_info , const char* label, struct fit_result fit_out,     vector<header_BSM> params, std::vector<int> myen){
    int Npar=fit_info.Npar;
    int Nvar=fit_info.Nvar+fit_info.n_ext_P;
    int Njack=gjack[0].Njack;
    int N=fit_info.N;
    char namefile[NAMESIZE];
    FILE *f;
    
    double **tif=swap_indices(fit_info.Npar,Njack,fit_out.P);
    double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
    double *tmpy=(double*) malloc(sizeof(double*)* Njack);
    
    for (int n=0;n< N; n++){
        
        mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_n%d_eta.txt",argv[3], label,n);
        f=open_file(namefile,"w+");
        double *tmpx=(double*) malloc(sizeof(double*)* Nvar);
        double *tmpy=(double*) malloc(sizeof(double*)* Njack);
        printf("writing: %s\n",namefile);
        
        for (int i=0 ; i<100; i++){
            for (int j=0;j<Njack;j++){
                double finalL=i;
                
                tmpx[0]=(double) params[myen[0]].L;
                tmpx[1]=(double) params[myen[0]].T;
                tmpx[2]=(double) params[myen[0]].rho;
                tmpx[3]=-1.5+ i/100.0;
                tmpx[4]=(double) params[myen[0]].csw;
                tmpx[5]=(double) params[myen[0]].mu03;
                tmpx[6]=(double) params[myen[0]].m0;
                //m0   put for each n the mass of the last ensemble
                
                    
                for(int i=fit_info.Nvar ; i<fit_info.Nvar+ fit_info.n_ext_P; i++)
                    tmpx[i]=fit_info.ext_P[i-fit_info.Nvar][j];
                
                            
                tmpy[j]=fit_info.function(n,Nvar,tmpx,Npar,tif[j]);
            }  
            fprintf(f,"%g  \t %g  %g\n",tmpx[3],tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
        }
        free(tmpy);free(tmpx);
        fclose(f); 
            
    }
    
    free_2(Njack,tif);
    
}






int main(int argc, char **argv){
    error(argc!=4,1,"main ",
          "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");   
    
    vector<header_BSM> paramsj;
    vector<data_BSM> dataj;
    
    int Ne=0; 
    header_BSM params;
    
    char namefile[NAMESIZE];
    
    //mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L16_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
    //emplace_back_par_data(namefile,paramsj,dataj);
    
     //0
    mysprintf(namefile,NAMESIZE,"%s/b5.85/L20T40/eta_m1.0983_M02_-0.010396_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.098300_csw1.000000_mu030.022400_m0-0.010396",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //1
     mysprintf(namefile,NAMESIZE,"%s/b5.85/L20T40/eta_m1.0983_M02_-0.024604_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.098300_csw1.000000_mu030.022400_m0-0.024604",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //1
     mysprintf(namefile,NAMESIZE,"%s/b5.85/L20T40/eta_m1.0983_M02_-0.040000_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.098300_csw1.000000_mu030.022400_m0-0.040000",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     
     
     mysprintf(namefile,NAMESIZE,"%s/b5.85/L20T40/eta_m1.2944_M02_-0.024604_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.294400_csw1.000000_mu030.022400_m0-0.024604",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
//      mysprintf(namefile,NAMESIZE,"%s/b5.85/L20T40/eta_m1.2944_M02_-0.040000_mu03_0.0224_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.294400_csw1.000000_mu030.022400_m0-0.040000",argv[2],argv[1]);
//      emplace_back_par_data(namefile,paramsj,dataj);
     mysprintf(namefile,NAMESIZE,"%s/b5.85/L20T40/eta_m1.2944_M02_-0.024604_mu03_0.0120_csw_1.0_rho1.96/jackknife/%s_T40_L20_rho1.960000_eta-1.294400_csw1.000000_mu030.012000_m0-0.024604",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     
     
     
     printf("%g   %g\n",dataj[0].jack[3][dataj[0].Njack-3] , error_jackboot(argv[1], dataj[0].Njack, dataj[0].jack[3]));
    /* for (int j=0; j<dataj[0].Njack; j++)
         printf("%d   %g\n",j,dataj[0].jack[3][j]);
    */ 
     vector<data_BSM> gjack= create_generalised_resampling( dataj );
     
     Ne=gjack.size();
     int Njack=gjack[0].Njack;
     printf("number of ensembles = %d\n",Ne);
     
     vector<int> myen(Ne);
     
     
     for(int i=0;i< Ne; i++)  myen[i]=i;
    
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     //////////////////////////////////////////////////////////////////////////////////////////////////
     struct fit_type fit_info;
     struct fit_result  fit_critical;
     fit_info.Nvar=7;
     fit_info.Npar=8;
     fit_info.N=2;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     fit_info.function=rhs_critical_eta_mu_m0;
     
     
     fit_critical=fit_data(argv,  paramsj ,gjack, lhs_critical_eta_mu_m0 ,fit_info, "eta_m0_critical_b585" ,myen);
     print_fit_band_eta( argv, gjack , fit_info ,  "eta_m0_critical_b585",   fit_critical ,    paramsj,  myen);
     
     
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // E20
     //////////////////////////////////////////////////////////////////////////////////////////////////
     
     return 0;
}
