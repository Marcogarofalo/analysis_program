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
#include "header_phi4.hpp"
#include "zeta_interpolation.hpp"
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>


//local folder
#include "mass_phi4.hpp"
#include "fit_function.hpp"



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
     char jackboot[NAMESIZE];
     mysprintf(jackboot,NAMESIZE,"%s" ,argv[1] );
     
     //mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L16_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     //emplace_back_par_data(namefile,paramsj,dataj);
    
     //0
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L28_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //1
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L30_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //2
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L32_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     //3
     mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L36_msq0-4.900000_msq1-4.650000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
     emplace_back_par_data(namefile,paramsj,dataj);
     
     
     printf("E1_0 =%f   %f\n", dataj[0].jack[1][dataj[0].Njack-1], error_jackboot(argv[1],dataj[0].Njack, dataj[0].jack[1]  ) );
    
     
     vector<data_phi> gjack= create_generalised_resampling( dataj );
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
     std::vector<int> Ls={16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56};
//      std::vector<int> Ls={28,30,32,36};
     std::vector<double> masses;
     std::vector<double> err_mass;
     int e1;
     for (int L: Ls){
         e1=Ne-1;
         for(int e=0; e< Ne;e++){
             if (L==paramsj[e].data.L[1]){
                 e1=e;
                 printf("L=%d found\n",L);
             }
         }
         masses.emplace_back(gjack[e1].jack[1][Njack-1]);
         err_mass.emplace_back(  error_jackboot(argv[1],Njack, gjack[e1].jack[1] )  );
     }
     
//       zeta.Init_Lmq(Ls, masses, err_mass  );
//      zeta_interpolation   zz;
//      zz.Init(argv[1] , myen , paramsj ,gjack  );
//      for(int e=0; e< Ne;e++){
//          printf("%g  %g\t",zeta.compute(28,0 ,12 *2* pi_greco/28., 0.5), zz.compute(28,0 ,12 *2* pi_greco/28., 0.5));
//          printf("%g  %g\t",zeta.compute(28,1 ,12 *2* pi_greco/28., 0.5), zz.compute(28,1 ,12 *2* pi_greco/28., 0.5));
//          printf("%g  %g\t",zeta.compute(28,2 ,12 *2* pi_greco/28., 0.5), zz.compute(28,2 ,12 *2* pi_greco/28., 0.5));
//          printf("%g  %g  %d %d\n",zeta.mass(e,0) , zz.mass(e,0),zeta.Ls(e),zz.Ls(e) );
//          printf(" %d %d %d %d %d %d\n",zeta.moms(0,0) , zz.moms(0,0),zeta.moms(0,1) , zz.moms(0,1),zeta.moms(0,1) , zz.moms(0,1) );
//          printf(" %d %d %d %d %d %d\n",zeta.moms(1,0) , zz.moms(1,0),zeta.moms(1,1) , zz.moms(1,1),zeta.moms(1,1) , zz.moms(1,1) );
//     }
//       zeta.write();
      zeta.read();
//      zeta_qsqg.Init(  );
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     //////////////////////////////////////////////////////////////////////////////////////////////////
     struct fit_type fit_info;
     struct fit_result  fit_m1, fit_m0;
     fit_info.Nvar=8;
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
    
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     
     fit_info.Npar=1;
     fit_info.N=1;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=muDE_rhs;
     
     struct fit_result fit_a_00=fit_data(argv,  paramsj ,gjack, muDE_00_lhs ,fit_info, "a_00_lusher",myen );
     
     printf("\n/////////////////////////////////     k cot delta    //////////////////\n");
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // kcot
     //////////////////////////////////////////////////////////////////////////////////////////////////
     
    
     fit_info.Npar=3;
     fit_info.N=5;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=rhs_kcotd;
     
     struct fit_result fit_kcotd=fit_data(argv,  paramsj ,gjack, lhs_kcotd ,fit_info, "kcotd",myen );
     
     ///////////////////////////////////////////////
     /// compute covariance energies
     int Nmom=5;
     double **E2corr=double_malloc_2(Nmom*gjack.size(),gjack[0].Njack);
     int **momlist=int_malloc_2(Nmom*gjack.size(),3);
     for(int e=0; e<gjack.size();e++ ){
         for (int j=0;j<gjack[0].Njack;j++ ){
             E2corr[e*Nmom+0][j]=gjack[e].jack[4][j];
             E2corr[e*Nmom+1][j]=gjack[e].jack[100][j];
             E2corr[e*Nmom+2][j]=gjack[e].jack[102][j];
             E2corr[e*Nmom+3][j]=gjack[e].jack[104][j];
             E2corr[e*Nmom+4][j]=gjack[e].jack[80][j];
         }
         momlist[e*Nmom+0][0]=0; momlist[e*Nmom+0][1]=0; momlist[e*Nmom+0][2]=0;
         momlist[e*Nmom+1][0]=1; momlist[e*Nmom+1][1]=0; momlist[e*Nmom+1][2]=0;
         momlist[e*Nmom+2][0]=1; momlist[e*Nmom+2][1]=1; momlist[e*Nmom+2][2]=0;
         momlist[e*Nmom+3][0]=1; momlist[e*Nmom+3][1]=1; momlist[e*Nmom+3][2]=1;
         momlist[e*Nmom+4][0]=0; momlist[e*Nmom+4][1]=0; momlist[e*Nmom+4][2]=0;
     }
     mysprintf(namefile,NAMESIZE,"%s/two_particle_energies.txt",argv[3] );
     FILE *f_two_particle=open_file(namefile,"w+");
     fprintf(f_two_particle,"#E2 dE2  P1  P2   P3   L  T\n");
     for(int e=0; e<gjack.size();e++ ){
         for(int n=0; n<Nmom;n++ ){
             fprintf(f_two_particle,"%.12g  %.12g  %d  %d   %d   %d  %d\n",
                     E2corr[e*Nmom+n][gjack[0].Njack-1],
                     error_jackboot(jackboot, gjack[0].Njack, E2corr[e*Nmom+n]),
                    momlist[e*Nmom+n][0],momlist[e*Nmom+n][1],momlist[e*Nmom+n][2],
                    paramsj[e].data.L[1],paramsj[e].data.L[0]     );
         }
     }
     fclose(f_two_particle);
     free(E2corr);
     ///////CM
     
     E2corr=double_malloc_2(Nmom*gjack.size(),gjack[0].Njack);
     for(int e=0; e<gjack.size();e++ ){
         for (int j=0;j<gjack[0].Njack;j++ ){
             E2corr[e*Nmom+0][j]=energy_CM(gjack[e].jack[4][j], momlist[e*Nmom+0],paramsj[e].data.L[1]  )/gjack[e].jack[1][j];
             E2corr[e*Nmom+1][j]=energy_CM(gjack[e].jack[100][j], momlist[e*Nmom+1],paramsj[e].data.L[1]  )/gjack[e].jack[1][j];
             E2corr[e*Nmom+2][j]=energy_CM(gjack[e].jack[102][j], momlist[e*Nmom+2],paramsj[e].data.L[1]  )/gjack[e].jack[1][j];
             E2corr[e*Nmom+3][j]=energy_CM(gjack[e].jack[104][j], momlist[e*Nmom+3],paramsj[e].data.L[1]  )/gjack[e].jack[1][j];
             E2corr[e*Nmom+4][j]=energy_CM(gjack[e].jack[80][j], momlist[e*Nmom+4],paramsj[e].data.L[1] )/gjack[e].jack[1][j];
         }
         
     }
     mysprintf(namefile,NAMESIZE,"%s/two_particle_energies_CM.txt",argv[3] );
     f_two_particle=open_file(namefile,"w+");
     fprintf(f_two_particle,"#E2CM/m0 dE2CM/m0  P1  P2   P3   L  T\n");
     for(int e=0; e<gjack.size();e++ ){
         for(int n=0; n<Nmom;n++ ){
             fprintf(f_two_particle,"%.12g  %.12g  %d  %d   %d   %d  %d\n",
                     E2corr[e*Nmom+n][gjack[0].Njack-1],
                     error_jackboot(jackboot, gjack[0].Njack, E2corr[e*Nmom+n]),
                     momlist[e*Nmom+n][0],momlist[e*Nmom+n][1],momlist[e*Nmom+n][2],
                     paramsj[e].data.L[1],paramsj[e].data.L[0]     );
         }
     }
     fclose(f_two_particle);
     double **cov=covariance(jackboot,Nmom*gjack.size() , gjack[0].Njack, E2corr);
     /*double **err_cov=error_covariance(jackboot,Nmom*gjack.size() , gjack[0].Njack, E2corr);
     mysprintf(namefile,NAMESIZE,"%s/two_particle_energies_covariance.txt",argv[3] );
     FILE *f_two_particle_covariance=open_file(namefile,"w+");
     for(int ne=0; ne<gjack.size()*Nmom;ne++ ){
         for(int ne1=0; ne1<gjack.size()*Nmom;ne1++ ){
             //fprintf( f_two_particle_covariance,"%.12g\t",cov[ne][ne1]/sqrt(cov[ne][ne]*cov[ne1][ne1]) );
             fprintf( f_two_particle_covariance,"%.12g  %.12g\t",cov[ne][ne1],err_cov[ne][ne1] );
         }
         fprintf(f_two_particle_covariance,"\n");
     }
     free_2(Nmom,cov);
     free_2(Nmom,err_cov);
     fclose(f_two_particle_covariance);
     */
//      ///////////////////////////////////////////////////////////////////////////////////////////////////
//      printf("\n/////////////////////////////////   fit  q  form from_phase_shift    //////////////////\n");
//      //////////////////////////////////////////////////////////////////////////////////////////////////
//      
//      
//      fit_info.Npar=2;
//      fit_info.N=2;
//      fit_info.Njack=gjack[0].Njack;
//      fit_info.n_ext_P=0;
//      //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
//      fit_info.function=rhs_q_from_phase_shift;
//      
//      struct fit_result q_from_phase_shift=fit_data(argv,  paramsj ,gjack, lhs_q ,fit_info, "k_from_phase_shift",myen ,  {-0.121902,8.20332} );// {-0.948817,-114.788,0.0003987}
//      
//      
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     printf("\n/////////////////////////////////   fit  k  form from_phase_shift    //////////////////\n");
     //////////////////////////////////////////////////////////////////////////////////////////////////
     
     
     fit_info.Npar=2;
     fit_info.N=3;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=rhs_k_from_phase_shift;
     
     struct fit_result k_from_phase_shift=fit_data(argv,  paramsj ,gjack, lhs_k ,fit_info, "k_from_phase_shift",myen ,  {-0.121902,8.20332} );// {-0.948817,-114.788,0.0003987}
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     printf("\n/////////////////////////////////   fit  k  form from_phase_shift   n4 //////////////////\n");
     //////////////////////////////////////////////////////////////////////////////////////////////////
     
     
     fit_info.Npar=2;
     fit_info.N=4;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=rhs_k_from_phase_shift;
     
     struct fit_result k_from_phase_shift_n4=fit_data(argv,  paramsj ,gjack, lhs_k ,fit_info, "k_from_phase_shift_n4",myen ,  {-0.121902,8.20332} );// {-0.948817,-114.788,0.0003987}
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     printf("\n/////////////////////////////////   fit  k  form from_phase_shift   n5  //////////////////\n");
     //////////////////////////////////////////////////////////////////////////////////////////////////
     
     
     fit_info.Npar=2;
     fit_info.N=5;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=rhs_k_from_phase_shift;
     
     struct fit_result k_from_phase_shift_n5=fit_data(argv,  paramsj ,gjack, lhs_k ,fit_info, "k_from_phase_shift_n5",myen ,  {-0.121902,8.20332} );// {-0.948817,-114.788,0.0003987}
     
     
     
     return 0;
}
