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
   ;
  

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
    
}


int main(int argc, char **argv){
    vector<cluster::IO_params> paramsj;
     vector<data_phi> dataj;
     
     int Ne=0; 
     cluster::IO_params params;
     
     char namefile[NAMESIZE];
     
     mysprintf(namefile,NAMESIZE,"../../tests/phi4/jackknife/jack_G2t_T32_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_reference");
     emplace_back_par_data(namefile,paramsj,dataj); 
     mysprintf(namefile,NAMESIZE,"../../tests/phi4/jackknife/jack_G2t_T32_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0");
     emplace_back_par_data(namefile,paramsj,dataj); 
     printf("observables reference: %d\n",dataj[0].Nobs);
     printf("observables test:      %d\n",dataj[1].Nobs);
     int nerrors=0;
     for(int i=0;i< dataj[0].Nobs;i++){
        for(int j=0;j< dataj[0].Njack;j++){
            double c=fabs( (dataj[0].jack[i][j]- dataj[1].jack[i][j] ) /dataj[0].jack[i][j]  );
            if(c>1e-6) { printf("error at   obs=%d  jack=%d     ref=%g  test=%g \n",i,j,dataj[0].jack[i][j],dataj[1].jack[i][j]);  nerrors++; }
        }
     }
     if (nerrors>0){
         exit(1);
    }
     
}

