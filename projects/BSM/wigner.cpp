#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>
#include <random>



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
#include "lhs_functions.hpp"
#include "header_BSM.hpp"

using namespace std;

struct  kinematic kinematic_2pt;


int r_value(int r)
{
    int vr;
    if (r==0) vr= 1;
    else if (r==1) vr= -1;
    else{ vr=-10; error(0==0,0,"r_value","r value is not 0 neither 1\n");}
    return vr;
}

void get_kinematic( int ik2,int r2, int ik1,int r1,int imom2, int imom1 ){
    kinematic_2pt.ik2=ik2;
    kinematic_2pt.ik1=ik1;
    kinematic_2pt.k2=file_head.k[ik2+file_head.nk];
    kinematic_2pt.k1=file_head.k[ik1+file_head.nk];
    kinematic_2pt.r2=-r_value(r2);
    kinematic_2pt.r1=r_value(r1);
    kinematic_2pt.mom2=-file_head.mom[imom2][1];
    kinematic_2pt.mom1=file_head.mom[imom1][1];
    
    kinematic_2pt.mom02=file_head.mom[imom2][0];
    kinematic_2pt.mom01=file_head.mom[imom1][0];
    
    kinematic_2pt.line=1;
    printf("%g %g  %d  %d\n",kinematic_2pt.k1,kinematic_2pt.k2,ik1+file_head.nk,ik2+file_head.nk);
    
}


int main(int argc, char **argv){
    int size;
    
    
    double ****data, **out,**tmp; 
    char c;
    clock_t t1,t2;
    double *in;
    
    char namefile[NAMESIZE];
    
    
    //double ****conf_jack,**r,**mt,**met;
    int Ncorr=1;
    int t0=2;
    
    
    
    FILE *plateaux_masses=NULL, *plateaux_masses_GEVP=NULL; 
    char namefile_plateaux[NAMESIZE];
    
    mysprintf(namefile_plateaux,NAMESIZE,"plateaux.txt");
    
    error(argc!=7,1,"main ",
          "usage:./phi4  blind/see/read_plateaux -p path  -bin $bin  jack/boot  \n separate path and file please");
    error(strcmp(argv[1],"blind")!=0 &&  strcmp(argv[1],"see")!=0 && strcmp(argv[1],"read_plateaux")!=0   ,1,"main ",
          "argv[1] only options:  blind/see/read_plateaux ");
    error(strcmp(argv[4],"-bin")!=0 ,1,"main", "argv[4] must be: -bin");
    error(strcmp(argv[6],"jack")!=0 &&  strcmp(argv[6],"boot")!=0,1,"main",
          "argv[6] only options: jack/boot");
    char **option;
    option=(char **) malloc(sizeof(char*)*7);
    option[0]=(char *) malloc(sizeof(char)*NAMESIZE);
    option[1]=(char *) malloc(sizeof(char)*NAMESIZE);
    option[2]=(char *) malloc(sizeof(char)*NAMESIZE);
    option[3]=(char *) malloc(sizeof(char)*NAMESIZE);
    option[4]=(char *) malloc(sizeof(char)*NAMESIZE);
    option[5]=(char *) malloc(sizeof(char)*NAMESIZE);
    option[6]=(char *) malloc(sizeof(char)*NAMESIZE);
    
    mysprintf(option[1],NAMESIZE,argv[1]); // blind/see/read_plateaux
    mysprintf(option[2],NAMESIZE,"-p"); // -p
    mysprintf(option[3],NAMESIZE,argv[3]); // path
    mysprintf(option[4],NAMESIZE,argv[6]); //resampling
    mysprintf(option[5],NAMESIZE,"no"); // pdf
    mysprintf(option[6],NAMESIZE,"plateau"); // infile
    
    std::vector<std::string>  correlators;
    correlators.emplace_back("JTILDEA1P1TRIVIAL");
    correlators.emplace_back("P1P1TRIVIAL");
    correlators.emplace_back("P2P2TRIVIAL");
    correlators.emplace_back("P3P3TRIVIAL");
    correlators.emplace_back("S0S0TRIVIAL");
    correlators.emplace_back("P1DP1NONSMEAREDNONTRIVIAL");
    correlators.emplace_back("VECTORDENSITY3DENSITY3NONTRIVIAL");
    correlators.emplace_back("phit");
    std::vector<FILE*>  f_correlators(correlators.size());
    header_BSM header;
    for (int i =0;i<correlators.size();i++){
        mysprintf(namefile,NAMESIZE,"%s/data/%s",  argv[3], correlators[i].c_str() );
        f_correlators[i]=open_file(namefile,"r+");
        printf("reading file: %s\n",namefile);
        if (i==0)
            read_header_BSM(header,f_correlators[i]);
        else 
            check_header_BSM(header,f_correlators[i],correlators[i]);
            
    }
    
    
    
    int T=header.T;
    file_head.l0=T;
    file_head.l1=header.L;file_head.l2=header.L;file_head.l3=header.L;
    file_head.nk=2;
    file_head.k= (double*) malloc(sizeof(double )*file_head.nk*2);
    file_head.k[0]=0;file_head.k[1]=0;
    file_head.k[2]=0;
    file_head.k[3]=0;
    
    file_head.nmoms=1;
    file_head.mom=(double**) malloc(sizeof(double*)*file_head.nmoms);
    for(int i=0;i<file_head.nmoms;i++) {
        file_head.mom[i]=(double*) malloc(sizeof(double)*4);
        file_head.mom[i][0]=0;
        file_head.mom[i][1]=0;
        file_head.mom[i][2]=0;
        file_head.mom[i][3]=0;
    }
    
    mysprintf(namefile,NAMESIZE,"%s/out/T%d_L%d_output",  argv[3], T,file_head.l1 );
    printf("writing output in :\n %s \n",namefile);
    FILE *outfile=open_file(namefile,"w+");   
    mysprintf(namefile,NAMESIZE,"%s/out/T%d_L%d_meff_correlators",  argv[3], T,file_head.l1 );
    FILE *outfile_meff_corr=open_file(namefile,"w+");    
    mysprintf(namefile,NAMESIZE,"%s/out/T%d_L%d_raw_correlators",  argv[3], T,file_head.l1 );
    FILE *outfile_raw_corr=open_file(namefile,"w+"); 
    mysprintf(namefile,NAMESIZE,"%s/out/T%d_L%d_shifted_correlators",  argv[3], T,file_head.l1 );
    FILE *outfile_shifted_corr=open_file(namefile,"w+"); 
    mysprintf(namefile,NAMESIZE,"%s/out/T%d_L%d_log_meff_shifted",  argv[3], T,file_head.l1 );
    FILE *outfile_log_meff_shifted=open_file(namefile,"w+"); 
    mysprintf(namefile,NAMESIZE,"%s/out/T%d_L%d_gamma",  argv[3], T,file_head.l1 );
    FILE *out_gamma=open_file(namefile,"w+");     
    
    
    
    int bin=atoi(argv[5]);
    int Neff=header.confs/bin;
    int confs=header.confs;
    error(confs<=0,1, "main","nconf<0");
    if( strcmp(argv[6],"jack")==0)
        header.Njack=Neff+1;
    else if( strcmp(argv[6],"boot")==0)
        header.Njack=Nbootstrap+1;
    else 
        error(1==1,1,"main","argv[6]= %s is not jack or boot",argv[6]);
    
    printf("configurations: %d , binning size: %d , resulting confs: %d ,  Njackboot: %d \n",header.confs, bin, Neff, header.Njack );
    // prepare jack file
    mysprintf(namefile,NAMESIZE,"%s/jackknife/%s_T%d_L%d_rho%.6f_eta%.6f_csw%.6f_mu03%.6f_m0%.6f",
              argv[3],option[4],
              header.T, header.L, header.rho,header.eta,header.csw,header.mu03,header.m0);
    FILE *jack_file=open_file(namefile,"w+");
    write_header_BSM_bin(header,jack_file);
    
    int var_to_read=correlators.size();
    correlators.emplace_back("JTILDEA1P1TRIVIALphi");
    correlators.emplace_back("P1DP1NONSMEAREDNONTRIVIALphi");
    int var=correlators.size();
    data=calloc_corr(confs, var,  header.T );
    int tau=5;
    
    for (int iconf=0; iconf< confs ;iconf++){
        
        for(int i =0 ; i< 8; i++){//var_to_read
            for (int t =0; t< header.T;t++)
                fscanf(f_correlators[i],"%lf  %lf\n",&data[iconf][i][t][0],&data[iconf][i][t][1]);
        }
                //read_twopt(f_correlators[i], iconf, &data[iconf][i], params,i);
        for (int t =0; t< header.T;t++){
            int tptau=(t+tau)%T;
            data[iconf][8][t][0]= (data[iconf][0][(t+1)%T][0]-data[iconf][0][t][0])*data[iconf][7][tptau][0];
            data[iconf][8][t][1]=0;
            data[iconf][9][t][0]= data[iconf][5][t][0]*data[iconf][7][tptau][0];
            data[iconf][9][t][1]=0;
        }
        for(int i =8 ; i<var_to_read ; i++){//var_to_read
            for (int t =0; t< header.T;t++)
                fscanf(f_correlators[i],"%lf  %lf\n",&data[iconf][i+2][t][0],&data[iconf][i+2][t][1]);
        }
        
    }
    
    
    
//forward_derivative_corr(confs,0,header.T,data);
    symmetrise_corr(confs, 1, header.T,data);    
    symmetrise_corr(confs, 2, header.T,data);
    symmetrise_corr(confs, 3, header.T,data);
    symmetrise_corr(confs, 4, header.T,data);
    
    
    ///resampling
    double ****data_bin=binning(confs, var, header.T ,data, bin);
    //if you want to do the gamma analysis you need to do before freeing the raw data
    //effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data_bin,  Neff ,namefile_plateaux,out_gamma,0,"M_{PS}^{ll}");
    //effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data,  confs ,namefile_plateaux,out_gamma,3,"M_{PS}^{ll}");
    
    free_corr(confs, var, header.T ,data);
    
    printf("before =%g\n", data_bin[0][0][0][0]);
    if(Neff==1){
        int Nfake=9;
        printf("only one conf found: generating a second one \n");
        Neff+=Nfake;
        double ****data_x2=calloc_corr(Neff, var,  header.T );
        std::mt19937 mt_rand(123);
            for(int k=0;k<var;k++)
                for (int i=0;i<header.T;i++)
                    for (int l=0;l<2;l++){
                        double r=mt_rand()/((double)mt_rand.max() );
                        r=r/10.+0.9;
                        for(int j=0;j<Neff;j++){
                            data_x2[j][k][i][l]=data_bin[0][k][i][l]*(1+ (j- Neff*(Neff-1)/(2.0*Neff))*r*0.1);
                        }
                        
                    }
        if( strcmp(argv[6],"jack")==0)
            header.Njack=Neff+1;
        else if( strcmp(argv[6],"boot")==0)
            header.Njack=Nbootstrap+1;
        else 
            error(1==1,1,"main","argv[6]= %s is not jack or boot",argv[6]);
        free_corr(Neff-10, var, header.T ,data_bin);// Neff-1 because Neff changed
        data_bin=data_x2;
    }
    
    double ****conf_jack=create_resampling(option[4],Neff, var, header.T, data_bin);
    free_corr(Neff, var, header.T ,data_bin);
    
    
    double *zeros=(double*) calloc(header.Njack,sizeof(double));
    printf("Neff=%d\n",Neff);
    
    double *tmp1=(double*) malloc(header.Njack);
    
//     error_jackboot(argv[6],header.Njack,conf_jack[0][1][1]) );
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //print all the effective masses correlators
    //set the option to not read for a plateaux
    char  save_option[NAMESIZE];
    sprintf(save_option,"%s",option[1]);
    sprintf(option[1],"blind");
    FILE *dev_null=open_file("/dev/null","w");
    get_kinematic( 0,0,  1, 0,0,  0 );
    
    for(int icorr=0; icorr<correlators.size(); icorr++ ){
        //log effective mass
        double *tmp_meff_corr  =plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  header.Njack
        ,namefile_plateaux,outfile_meff_corr,icorr,(char*) correlators[icorr].c_str(), M_eff_log,dev_null);
        free(tmp_meff_corr);
        //raw correlator
        tmp_meff_corr  =plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,   header.Njack ,
                                                      namefile_plateaux,outfile_raw_corr,icorr,(char*) correlators[icorr].c_str(), identity,dev_null);
        free(tmp_meff_corr);
        // shifted correlator
        tmp_meff_corr  =plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,   header.Njack ,
                                                      namefile_plateaux,outfile_shifted_corr,icorr,(char*) correlators[icorr].c_str(), shift_corr,dev_null);
        free(tmp_meff_corr);
        // log_meff shifted correlator
        tmp_meff_corr  =plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,   header.Njack
        ,namefile_plateaux,outfile_log_meff_shifted,icorr,(char*) correlators[icorr].c_str(), M_eff_log_shift,dev_null);
        free(tmp_meff_corr);
    }
    sprintf(option[1],"%s",save_option);// restore option
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // 0
    double *mass=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,   header.Njack ,namefile_plateaux,outfile,1,"m_PS", M_eff_T,jack_file);
    // 1
    double *massDs3=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,   header.Njack ,namefile_plateaux,outfile,6,"m_DS3", M_eff_T,jack_file);
    // 2
    double *masss0=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,   header.Njack ,namefile_plateaux,outfile,4,"m_S0", M_eff_T,jack_file);
    
    
    fit_type fit_info;

    fit_info.Nvar=1;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=header.Njack;
    fit_info.n_ext_P=0;
    fit_info.function=constant_fit;
    
    //file_head.k[2]=mu1;    file_head.k[3]=mu2;
    
    //c++ 3 || r 2
    fit_result fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "A1P1phi", conf_jack ,namefile_plateaux, outfile, 
                                              r_AWI, "r_AWI",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    // 4
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "A1P1", conf_jack ,namefile_plateaux, outfile, 
                                              m_PCAC, "m_PCAC",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
      
    
    
    
}
