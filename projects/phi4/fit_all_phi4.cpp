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
#include "mass_phi4.hpp"
#include "fit_function.hpp"
#include "header_phi4.hpp"


using namespace std;
int Ne=0;





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
    
     vector<cluster::IO_params> paramsj;
     vector<data_phi> dataj;
     
     int Ne=0; 
     cluster::IO_params params;
     
     char namefile[NAMESIZE];
     
     //mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L16_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     //emplace_back_par_data(namefile,paramsj,dataj);
    
     //0
    /* mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //1
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L24_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //2
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L26_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //3
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L32_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
   */  //4
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T48_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //5
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //6
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T96_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //7
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T128_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     
     
     //8
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T96_L22_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //9
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T96_L24_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //10
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T96_L26_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //11
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T96_L32_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     
     
     
     
     printf("E1_0 =%f   %f\n", dataj[0].jack[1][dataj[0].Njack-1], error_jackboot(argv[1],dataj[0].Njack, dataj[0].jack[1]  ) );
    
     
     vector<data_phi> gjack= create_generalised_resampling( dataj );
     printf("GEVP_E2_01 =%f   %f\n", gjack[1].jack[19][gjack[1].Njack-1], error_jackboot(argv[1],gjack[1].Njack, gjack[1].jack[19]  ) );
     Ne=gjack.size();
     printf("number of ensembles = %d\n",Ne);
     
     vector<int> myen(Ne);
     for(int i=0;i< Ne; i++)  myen[i]=i;
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     //////////////////////////////////////////////////////////////////////////////////////////////////
     int Njack=gjack[0].Njack;
     struct fit_type fit_info;
     struct fit_result  fit_m1, fit_m0;
     fit_info.Nvar=11;
     fit_info.Npar=2;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     fit_info.function=M_finite_volume;

     
     fit_m0=fit_data(argv,  paramsj ,gjack, M0_finite_volume_lhs ,fit_info, "M0_finite_vol" ,myen);
    
     
     fit_m1=fit_data(argv,  paramsj ,gjack, M1_finite_volume_lhs ,fit_info, "M1_finite_vol" ,myen);
     
     
     printf("\n/////////////////////////////////     E2_0//////////////////\n");
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // E20
     //////////////////////////////////////////////////////////////////////////////////////////////////
     
     fit_info.Npar=1;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=2;
     fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=muDE_rhs;
     
     fit_info.ext_P[0]=fit_m0.P[0];
     fit_info.ext_P[1]=fit_m0.P[0];
     struct fit_result fit_a_00_infm=fit_data(argv,  paramsj ,gjack, muDE_00_infm_lhs ,fit_info, "a_00_luscher_infm" ,myen);
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     
     fit_info.Npar=1;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=muDE_rhs;
     
     struct fit_result fit_a_00=fit_data(argv,  paramsj ,gjack, muDE_00_lhs ,fit_info, "a_00_luscher",myen );
     
     
    
     printf("\n/////////////////////////////////     E2_01//////////////////\n");
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     
     fit_info.Npar=1;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=muDE_rhs;
     
     struct fit_result fit_a_01=fit_data(argv,  paramsj ,gjack, muDE_01_lhs ,fit_info, "a_01_luscher",myen );
   
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     
     fit_info.Npar=1;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=constant_fit;
     
     struct fit_result fit_a_01_const=fit_data(argv,  paramsj ,gjack, a_01_luescher_lhs ,fit_info, "a_01_luscher_const",myen );
     
     
     fit_info.Npar=1;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=muDE_rhs;
     
     struct fit_result fit_a_01_div_shift=fit_data(argv,  paramsj ,gjack, muDE_01_div_shift_lhs ,fit_info, "a_01_luscher_div_shift",myen );
     
     
     /////////////reload the data
    /* cout << endl <<"reload the data" << endl;
     vector<data_phi>().swap(dataj);
     vector<data_phi>().swap(gjack);
     vector<cluster::IO_params>().swap(paramsj);
     dataj.resize(0);
     gjack.resize(0);
     paramsj.resize(0);
     
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T48_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T96_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     gjack= create_generalised_resampling( dataj );
     myen.resize(gjack.size());
     for(int i=0;i< gjack.size(); i++)  myen[i]=i;
     */
     
     printf("size=%ld\n",myen.size());
//      myen={4,5,6,7,8,9,10};
     printf("size=%ld\n",myen.size());
     
     printf("\n/////////////////////////////////     a_00_BH//////////////////\n");
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     
     fit_info.Npar=1;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=constant_fit;
     
     struct fit_result fit_a_00_BH=fit_data(argv,  paramsj ,gjack, a_00_BH_lhs ,fit_info, "a_00_BH" ,myen);
   
     
     printf("\n/////////////////////////////////     a_01_BH_03t16//////////////////\n");
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     
     fit_info.Npar=1;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=constant_fit;
     
     //lhs<53>  E4_03t16_shifted
     struct fit_result fit_a_01_BH_03t16=fit_data(argv,  paramsj ,gjack, lhs<53> ,fit_info, "a_01_BH_03t16_shifted" ,myen);
     
     
     printf("\n/////////////////////////////////     a_01_BH_02t10//////////////////\n");
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     
     fit_info.Npar=1;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=constant_fit;
     
     //lhs<53>  E4_03t16_shifted
     struct fit_result fit_a_01_BH_02t10=fit_data(argv,  paramsj ,gjack, lhs<73> ,fit_info, "a_01_BH_02t10_shifted" ,myen);
     
     
     printf("\n/////////////////////////////////     a_01_BH_04t16//////////////////\n");
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     
     fit_info.Npar=1;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=constant_fit;
     //lhs<53>  E4_04t16_shifted
     struct fit_result fit_a_01_BH_04t16=fit_data(argv,  paramsj ,gjack, lhs<56> ,fit_info, "a_01_BH_04t16_shifted" ,myen);
     
     mysprintf(namefile,NAMESIZE,"%s/%s.tex",argv[3], "PT_test_BH");
     FILE *f=open_file(namefile,"w+");
    
     fprintf(f,"\\begin{gather}\n");
     double *tmp=(double*) malloc(sizeof(double)*gjack[0].Njack);
     for (int j=0;j<gjack[0].Njack; j++)
         tmp[j]=fit_a_00_BH.P[0][j]*fit_m0.P[0][j];
     fprintf(f,"a_{00}  m_0=%f  \\pm %f \\\\ \n", tmp[Njack-1],error_jackboot(argv[1],Njack,tmp )  );
     for (int j=0;j<gjack[0].Njack; j++)
         tmp[j]=fit_a_01_BH_03t16.P[0][j]*(fit_m0.P[0][j]+fit_m1.P[0][j]);
     fprintf(f,"a_{01}  (m_0+m_1)=%f  \\pm %f\\\\ \n", tmp[Njack-1],error_jackboot(argv[1],Njack,tmp )  );
     for (int j=0;j<gjack[0].Njack; j++)
         tmp[j]=fit_a_01_BH_03t16.P[0][j]*(fit_m0.P[0][j]+fit_m1.P[0][j])   /(  fit_a_00_BH.P[0][j]*fit_m0.P[0][j]   );
     fprintf(f,"\\frac{a_{01}  (m_0+m_1)}{a_{00}  m_0}=%f  \\pm %f \n", tmp[Njack-1],error_jackboot(argv[1],Njack,tmp )  );
     fprintf(f,"\\end{gather}\n");

     free(tmp);
     
     return 0;
}
