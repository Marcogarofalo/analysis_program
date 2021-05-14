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


void read_Njack_Nobs( FILE *stream, cluster::IO_params params, int &Njack, int &Nobs ){

   long int tmp;
   int s=params.data.header_size;
    
   fread(&Njack, sizeof(int), 1, stream );
   
   
   fseek(stream, 0, SEEK_END);
   tmp = ftell(stream);
   tmp-= params.data.header_size+sizeof(int) ;
   
   s=Njack;

   Nobs= (tmp)/ ((s)*sizeof(double) );

   fseek(stream, params.data.header_size+sizeof(int), SEEK_SET);
   
  

}

void read_dataj(FILE *stream,cluster::IO_params params, data_phi &dj){
    read_Njack_Nobs(stream, params, dj.Njack, dj.Nobs );
    //printf("Nobs=%d   Njack=%d\n",dj.Nobs,dj.Njack)
    dj.jack=double_malloc_2( dj.Nobs, dj.Njack);
    
    for (int obs=0; obs<dj.Nobs; obs++ ){
        fread(dj.jack[obs], sizeof(double ), dj.Njack, stream );
    }
    
}

void emplace_back_par_data( char *namefile , vector<cluster::IO_params> &paramsj, vector<data_phi> &dataj){
    cluster::IO_params params;
    data_phi  data;
    FILE *f=open_file(namefile,"r");
    read_header( f  , params);
    read_dataj(f,params,data );
    fclose(f);

    paramsj.emplace_back(params);
    dataj.emplace_back(data);
    //printf("E1=%g    %g\n",dataj[0].jack[1][   dataj[0].Njack-1 ],    data.jack[1][data.Njack-1]);
    printf("ending the function: emplace_back_par_data\n");
}

vector<data_phi> create_generalised_resampling(  vector<data_phi> &dataj ){
    // if the length is the same return dataj
    int same=0;
    for( auto &d :dataj){
        printf("jacks=%d\n",d.Njack);
        if (d.Njack==dataj[0].Njack)
            same++;
    }
    if (same==dataj.size()){
            cout << "all the files have the same number of jack/boot , do nothing"<<endl;
            return dataj;
    }
    else{
        cout << "creating generalised jack"<<endl;
        vector<data_phi> gjack;
        //gjack.resize(dataj.size());
        //jac_tot is the summ of all jackknife +1 
        //remember alle the dataj have one extra entry for the mean
        int jack_tot=0;
        for( int e=0 ; e<dataj.size();e++)
            jack_tot+=dataj[e].Njack;
        jack_tot=jack_tot-dataj.size()+1;
        cout << "ensembles "<< dataj.size() << endl;
        cout<< "jack tot= "<< jack_tot<< endl;
        
        //get Nobs the minimum number of observable between the diles
        int Nobs=1000;
        for( auto &d :dataj)
            if(Nobs> d.Nobs )
                Nobs=d.Nobs;
       
        
        
        
        for(int e=0;e<dataj.size();e++){
            data_phi tmp;
            tmp.Njack=jack_tot;
            tmp.Nobs=Nobs;
            cout<< Nobs<<" " <<jack_tot<< endl;
            tmp.jack=double_malloc_2( Nobs, jack_tot);
            int counter=0;
            for(int e1=0;e1<dataj.size();e1++){
                for(int j=0;j<(dataj[e1].Njack-1);j++){
                    for(int o=0;o<Nobs;o++){ 
                        if (e==e1){
                            tmp.jack[o][j+counter]=dataj[e].jack[o][j];
                        }
                        else{
                            tmp.jack[o][j+counter]=dataj[e].jack[o][ dataj[e].Njack-1  ];
                        }
                    }
                }
                counter+=dataj[e1].Njack-1 ;
            }
            for(int o=0;o<Nobs;o++)
                tmp.jack[o][ jack_tot-1 ]=dataj[e].jack[o][ dataj[e].Njack-1  ];
            
            free_2(dataj[e].Nobs, dataj[e].jack);
            gjack.emplace_back(tmp);
        }
        vector<data_phi>().swap(dataj);
        
        return gjack;
        
    }

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
     read_header( f  , params[0]);
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
     
     vector<int> myen(Ne);
     for(int i=0;i< Ne; i++)  myen[i]=i;
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     // start fitting
     //////////////////////////////////////////////////////////////////////////////////////////////////
     int Njack=gjack[0].Njack;
     struct fit_type fit_info;
     struct fit_result  fit_m1, fit_m0;
     fit_info.Nvar=7;
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
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     printf("\n/////////////////////////////////   fit  k  form from_phase_shift    //////////////////\n");
     //////////////////////////////////////////////////////////////////////////////////////////////////
     
     
     fit_info.Nvar=7;
     fit_info.Npar=3;
     fit_info.N=5;
     fit_info.Njack=gjack[0].Njack;
     fit_info.n_ext_P=0;
     //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
     fit_info.function=rhs_k_from_phase_shift;
     
     struct fit_result k_from_phase_shift=fit_data(argv,  paramsj ,gjack, lhs_k ,fit_info, "k_from_phase_shift",myen , {-0.948817,-114.788,0.0003987});
     
     
     
     return 0;
}
