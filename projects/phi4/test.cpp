#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <memory>
#include <vector>

#include "global.hpp"
#include "resampling.hpp"
#include "mutils.hpp"
#include "header_phi4.hpp"
#include "tower.hpp"

using namespace std;

int main(int argc, char **argv){
    vector<cluster::IO_params> paramsj;
     vector<data_phi> dataj;
     
     int Ne=0; 
     cluster::IO_params params;
     
     char namefile[NAMESIZE];
     error(argc!=3,1,"usage:","./test file1 file2");
     
//     mysprintf(namefile,NAMESIZE,"../../tests/phi4/jackknife/jack_G2t_T48_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0_reference");
     emplace_back_par_data(argv[1],paramsj,dataj); 
//     mysprintf(namefile,NAMESIZE,"../../tests/phi4/jackknife/jack_G2t_T48_L20_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0");
     emplace_back_par_data(argv[2],paramsj,dataj); 
     printf("observables reference: %d\n",dataj[0].Nobs);
     printf("observables test:      %d\n",dataj[1].Nobs);
     int nerrors=0;
     for(int i=0;i< dataj[0].Nobs;i++){
        for(int j=0;j< dataj[0].Njack;j++){
            double err=error_jackboot("jack", dataj[0].Njack, dataj[0].jack[i]  );
            double c=fabs( (dataj[0].jack[i][j]- dataj[1].jack[i][j] ) / dataj[0].jack[i][j] );
            if(c>1e-6) { printf("error at   obs=%d  jack=%d     ref=%.10g  test=%.10g err=%.10g \n",i,j,dataj[0].jack[i][j],dataj[1].jack[i][j], err);  nerrors++; }
        }
     }
     if (nerrors>0){
         exit(1);
    }
    free_2(dataj[0].Nobs,dataj[0].jack);
    free_2(dataj[1].Nobs,dataj[1].jack);
     
}

