#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"
#include "read.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"
#include "correlators_analysis.hpp"
#include "gamma_analysis.hpp"
#include "tower.hpp"
#include "gamma_analysis.hpp"


#include "correlators_analysis.hpp"
#include "eigensystem.hpp"


#include <string>
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>
// using namespace std;
struct  kinematic kinematic_2pt;

double afm=0.0908026;
double aMeV=0.0908026/ 197.326963;

generic_header read_head(FILE* stream) {
    generic_header header;
    size_t fi = 0;
    fi += fscanf(stream, "%d", &header.L);
    fi += fscanf(stream, "%d", &header.T);
    header.mus = std::vector<double>(2);
    fi += fscanf(stream, "%lf", &header.mus[0]);
    fi += fscanf(stream, "%lf", &header.kappa);

    int a, b;
    fi += fscanf(stream, "%d", &a);
    fi += fscanf(stream, "%d", &b);
    header.ncorr = a * b;
    fi += fscanf(stream, "%d", &header.Njack);
    header.size = header.ncorr * 2 * header.T;
    return header;
}
void  write_header_g2(FILE* jack_file, generic_header head) {
    int fi = 0;
    fi += fwrite(&head.T, sizeof(int), 1, jack_file);
    fi += fwrite(&head.L, sizeof(int), 1, jack_file);
    int nmu = head.mus.size();
    fi += fwrite(&nmu, sizeof(int), 1, jack_file);
    for (double mu : head.mus) {
        fi += fwrite(&mu, sizeof(double), 1, jack_file);
    }

}




void read_twopt(FILE* stream, double*** to_write, generic_header head) {
    int fi = 0;
    for (int k = 0; k < head.ncorr; k++) {
        for (int t = 0; t < head.T;t++) {
            fi += fscanf(stream, "%lf  %lf", &to_write[k][t][0], &to_write[k][t][1]);
            // printf(" corr=%d t=%d %.12g   %.12g\n", k, t, to_write[k][t][0], to_write[k][t][1]);
        }
    }

}



double DeltaM(int j, double**** in, int t, struct fit_type fit_info) {

    double m = fit_info.ext_P[0][j];
    double A = fit_info.ext_P[1][j];
    double T = fit_info.T;
    // double f = A * (t * exp(-m * t) + (T - t) * exp(-m * (T - t)));
    // double fp = A * ((t + 1.) * exp(-m * (t + 1)) + (T - (t + 1.)) * exp(-m * (T - (t + 1))));
    double f = (t * exp(-m * t) + (T - t) * exp(-m * (T - t))) / (exp(-m * t) + exp(-m * (T - t)));
    double fp = ((t + 1.) * exp(-m * (t + 1)) + (T - (t + 1.)) * exp(-m * (T - (t + 1)))) / (exp(-m * (t + 1)) + exp(-m * (T - (t + 1))));


    double r = -(in[j][3][t + 1][1] / in[j][0][t + 1][0] - in[j][3][t][1] / in[j][0][t][0]) / (fp - f);

    return r;
}


double me_P5P5(int j, double**** in, int t, struct fit_type fit_info) {

    double m = fit_info.ext_P[0][j];
    double T = fit_info.T;
    double r = in[j][0][t][0] / (exp(-m * t) + exp(-m * (T - t)));

    return r;
}

double DeltaMt(int j, double**** in, int t, struct fit_type fit_info) {

    double T = fit_info.T;

    double m = M_eff_T(t, T, in[j][0]);
    double mp = M_eff_T(t + 1, T, in[j][0]);
    double A = me_P5P5(j, in, t, fit_info);
    double Ap = me_P5P5(j, in, t + 1, fit_info);
    double f = A * (t * exp(-m * t) + (T - t) * exp(-m * (T - t)));
    double fp = Ap * ((t + 1) * exp(-mp * (t + 1)) + (T - (t + 1)) * exp(-mp * (T - (t + 1))));

    double r = -(in[j][3][t + 1][1] / in[j][0][t + 1][0] - in[j][3][t][1] / in[j][0][t][0]) / (fp - f);

    return r;
}
double mpcac(int j, double**** in, int t, struct fit_type fit_info) {


    double r = -(in[j][1][t + 1][1] - in[j][1][t][1]) / (2. * in[j][0][t][0]);

    return r;
}

double deltam(int j, double**** in, int t, struct fit_type fit_info) {


    double r = in[j][1][t][1] / (2. * in[j][4][t][0]);

    return r;
}

double deltam_sub(int j, double**** in, int t, struct fit_type fit_info) {
    double T = fit_info.T;
    double m = fit_info.ext_P[0][j];
    double DM = fit_info.ext_P[1][j];
    double g = DM * (t * exp(-m * t) - (T - t) * exp(-m * (T - t))) / (exp(-m * t) - exp(-m * (T - t)));
    double r =  (2. * in[j][4][t][0])/in[j][1][t][1]  -2* g;

    return 1/r;
}


// double *mass_gamma(int var, int order,int flow ,double *ah){
//     double *r=(double*) calloc((1),sizeof(double)); 
    
//     //use flow at time of the correlator
//     // we need to pass to M_eff a correlator so we create a double **c with has only the correlator at time=flow and flow+1
//     double **c=double_malloc_2(flow+2,2);
//     c[flow][0]=ah[0];
//     c[flow+1][0]=ah[1];
//     r[0]=M_eff(flow,c);
//     free_2(flow+2,c);
//     return r;
// }

void effective_mass_phi4_gamma(char **option ,struct kinematic kinematic_2pt , char* name, double ****data, int Confs ,const char *plateaux_masses,FILE *outfile,  int index , const char *description ){

   
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Confs);
    

   fprintf(outfile,"#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n",name,
           kinematic_2pt.k2,kinematic_2pt.r2,kinematic_2pt.mom2,
           kinematic_2pt.k1,kinematic_2pt.r1, kinematic_2pt.mom1 );
   for(i=1;i<file_head.l0/2;i++){  
           double *datag;
           // 2 component ot compute the mass
           datag=(double *) malloc(sizeof(double) * Confs *2 ); 
           for (j=0;j<Confs;j++){
              //store in the format for the gamma analysis
              //two variables c(t) and c(t+1), order=0 
              datag[j*2]=data[j][index][i][0];
              datag[1+j*2]=data[j][index][i+1][0];
           }
           
           double *obs=analysis_gamma(  2 , 1, Confs, i//time
                                     , datag  , mass_gamma);
          
           fprintf(outfile,"%d   %.15e    %.15e   %.15e  %.15e  %.15e\n",i,obs[0],obs[1],obs[2],obs[3],obs[4]);
           free(obs);
           free(datag);

   }

  
   for(i=0;i<file_head.l0;i++)
       free(r[i]);
   free(r);

   fflush(outfile);
    
}

int main(int argc, char** argv) {
    error(argc != 7, 1, "nissa_mpcac ",
        "usage:./nissa_mpcac -p path file -bin $bin  jack/boot \n separate path and file please");

    char** option;
    option = (char**)malloc(sizeof(char*) * 7);
    option[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, "read_plateaux"); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p"); // -p
    mysprintf(option[3], NAMESIZE, argv[2]); // path
    mysprintf(option[4], NAMESIZE, argv[6]); //resampling
    mysprintf(option[5], NAMESIZE, "no"); // pdf
    mysprintf(option[6], NAMESIZE, argv[3]); // infile

    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], option[6]);

    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");


    FILE* infile = open_file(namefile, "r");
    generic_header head = read_head(infile);


    file_head.l1 = head.L;
    file_head.l0 = head.T;
    file_head.l2 = head.L;
    file_head.l3 = head.L;
    file_head.nk = 2;
    file_head.musea = head.mus[1];
    file_head.k = (double*)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;file_head.k[1] = 0;
    file_head.k[2] = head.mus[1];
    file_head.k[3] = head.mus[1];

    file_head.nmoms = 1;
    file_head.mom = (double**)malloc(sizeof(double*) * file_head.nmoms);
    for (int i = 0;i < file_head.nmoms;i++) {
        file_head.mom[i] = (double*)malloc(sizeof(double) * 4);
        file_head.mom[i][0] = 0;
        file_head.mom[i][1] = 0;
        file_head.mom[i][2] = 0;
        file_head.mom[i][3] = 0;
    }


    int confs = head.Njack;
    int bin = atoi(argv[5]);
    int Neff = confs / bin;
    int Njack;
    if (strcmp(argv[6], "jack") == 0) {
        Njack = Neff + 1;
        myres = new resampling_jack(Neff);
    }
    else if (strcmp(argv[6], "boot") == 0) {
        Njack = (Neff * 2 + 1);
        myres = new resampling_boot(Neff * 2);
    }
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[7]= %s is not jack or boot", argv[7]);
    }
    mysprintf(namefile, NAMESIZE, "%s/out/%s_output", option[3], option[6]);
    printf("writing output in :\n %s \n", namefile);
    FILE* outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s", option[3], option[4], option[6]);
    FILE* jack_file = open_file(namefile, "w+");
    write_header_g2(jack_file, head);


    double**** data = calloc_corr(confs, head.ncorr, head.T);


    printf("confs=%d\n", confs);
    printf("ncorr=%d\n", head.ncorr);
    printf("kappa=%g\n", head.kappa);
    for (int iconf = 0; iconf < confs;iconf++) {
        read_twopt(infile, data[iconf], head);
    }

    // symmetrise_corr(confs, 0, head.T, data);
    // antisymmetrise_corr(confs, 1, head.T, data);
    // antisymmetrise_corr(confs, 2, head.T, data);



    double**** data_bin = binning(confs, head.ncorr, head.T, data, bin);
    double**** conf_jack = myres->create(Neff, head.ncorr, head.T, data_bin);

    free_corr(confs, head.ncorr, head.T, data);



    /////////////////////////////////////////////////////////////////////////////////////////////////////////
           //print all the effective masses correlators
           //set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/%s_meff_correlators", option[3], option[6]);
    FILE* outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_raw_correlators", option[3], option[6]);
    FILE* outfile_raw_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_shifted_correlators", option[3], option[6]);
    FILE* outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_log_meff_shifted", option[3], option[6]);
    FILE* outfile_log_meff_shifted = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_gamma", option[3], option[6]);
    FILE* out_gamma = open_file(namefile, "w+");

    char  save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");
    FILE* dev_null = open_file("/dev/null", "w");
    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    fit_info_silent.chi2_gap_jackboot = 1e+6;
    fit_info_silent.guess_per_jack = 0;

    for (int icorr = 0; icorr < head.ncorr; icorr++) {
        //log effective mass
        double* tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack
            , namefile_plateaux, outfile_meff_corr, icorr, "log", M_eff_log, dev_null, fit_info_silent);
        free(tmp_meff_corr);
        //raw correlator
        file_head.l0 = head.T * 2;
        tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity, dev_null, fit_info_silent);
        free(tmp_meff_corr);
        file_head.l0 = head.T;
        // shifted correlator
        tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_shifted_corr, icorr, "shift_cor", shift_corr, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        // log_meff shifted correlator
        tmp_meff_corr = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack
            , namefile_plateaux, outfile_log_meff_shifted, icorr, "log_shift", M_eff_log_shift, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);

        effective_mass_phi4_gamma(option, kinematic_2pt, (char*)"P5P5", data_bin, Neff, namefile_plateaux, out_gamma, 0, "M_{PS}");

    }
    free_corr(Neff, head.ncorr, head.T, data_bin);

    symmetrise_jackboot(Njack, 0, head.T, conf_jack);
    symmetrise_jackboot(Njack, 1, head.T, conf_jack, -1);
    symmetrise_jackboot(Njack, 2, head.T, conf_jack, -1);
    //insP
    // symmetrise_jackboot(Njack, 3, head.T, conf_jack);
    // symmetrise_jackboot(Njack, 4, head.T, conf_jack, -1);
    // symmetrise_jackboot(Njack, 5, head.T, conf_jack, -1);

    fit_info_silent.restore_default();
    sprintf(option[1], "%s", save_option);// restore option
    corr_counter = -1;

    double* M_PS = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 0, "M_{PS}", M_eff_T, jack_file);
    check_correlatro_counter(0);
    // double *MPSMEV=(double*) malloc(sizeof(double)*Njack);
    // for(int j=0; j<Njack;j++){
    //     MPSMEV[j]=M_PS[j]/aMeV;
    // }
    // fprintf(outfile,"%g  %g %g %g\n",aMeV,0.0,MPSMEV[Njack-1], myres->comp_error(MPSMEV));
 
    struct fit_type fit_info;
    struct fit_result  fit_out;

    // fit_info.Nvar = 1;
    // fit_info.Npar = 1;
    // fit_info.N = 1;
    // fit_info.Njack = Njack;
    // fit_info.n_ext_P = 1;
    // fit_info.ext_P = (double**)malloc(sizeof(double*) * 1);
    // fit_info.ext_P[0] = M_PS;
    // fit_info.function = constant_fit;
    // fit_info.linear_fit = true;
    // fit_info.T = head.T;

    // //c++ 1 || r 2
    // struct fit_result fit_me_P5P5 = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, me_P5P5, "me_P5P5", fit_info, jack_file);
    // // free_fit_result(fit_info, fit_out);
    // fit_info.restore_default();


    // fit_info.Nvar = 1;
    // fit_info.Npar = 1;
    // fit_info.N = 1;
    // fit_info.Njack = Njack;
    // fit_info.n_ext_P = 2;
    // fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);
    // fit_info.ext_P[0] = M_PS;
    // fit_info.ext_P[1] = fit_me_P5P5.P[0];

    // fit_info.function = constant_fit;
    // fit_info.linear_fit = true;
    // fit_info.T = head.T;

    // //c++ 1 || r 2
    // struct fit_result fit_DeltaM = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, DeltaM, "DeltaM", fit_info, jack_file);
    // // free_fit_result(fit_info, fit_out);
    // fit_info.restore_default();


    // fit_info.Nvar = 1;
    // fit_info.Npar = 1;
    // fit_info.N = 1;
    // fit_info.Njack = Njack;
    // fit_info.n_ext_P = 2;
    // fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);
    // fit_info.ext_P[0] = M_PS;
    // fit_info.ext_P[1] = fit_me_P5P5.P[0];

    // fit_info.function = constant_fit;
    // fit_info.linear_fit = true;
    // fit_info.T = head.T;

    // //c++ 1 || r 2
    // struct fit_result fit_DeltaMt = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, DeltaMt, "DeltaMt", fit_info, jack_file);
    // // free_fit_result(fit_info, fit_out);
    // fit_info.restore_default();




    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 0;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * 0);
    fit_info.function = constant_fit;
    fit_info.linear_fit = true;

    //c++ 1 || r 2
    fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, mpcac, "mpcac", fit_info, jack_file);
    free_fit_result(fit_info, fit_out);
    fit_info.restore_default();





    // fit_info.Nvar = 1;
    // fit_info.Npar = 1;
    // fit_info.N = 1;
    // fit_info.Njack = Njack;
    // fit_info.n_ext_P = 0;
    // fit_info.ext_P = (double**)malloc(sizeof(double*) * 0);
    // fit_info.function = constant_fit;
    // fit_info.linear_fit = true;

    // //c++ 1 || r 2
    // fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, deltam, "deltam", fit_info, jack_file);
    // double* true_kappa = (double*)malloc(sizeof(double) * Njack);
    // for (int j = 0; j < Njack;j++) {
    //     true_kappa[j] = head.kappa / (1. + 2. * fit_out.P[0][j] * head.kappa);
    // }
    // fprintf(outfile, "%.12g   %.12g\n", true_kappa[Njack - 1], myres->comp_error(true_kappa));
    // printf("true_kappa=%.12g   %.12g\n", true_kappa[Njack - 1], myres->comp_error(true_kappa));
    // free_fit_result(fit_info, fit_out);
    // fit_info.restore_default();

    // fit_info.Nvar = 1;
    // fit_info.Npar = 1;
    // fit_info.N = 1;
    // fit_info.Njack = Njack;
    // fit_info.n_ext_P = 2;
    // fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
    // fit_info.ext_P[0] = M_PS;
    // fit_info.ext_P[1] = fit_DeltaM.P[0];
    // fit_info.T = head.T;
    // fit_info.function = constant_fit;
    // fit_info.linear_fit = true;

    // //c++ 1 || r 2
    // fit_out = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux, outfile, deltam_sub, "deltam_sub", fit_info, jack_file);
    // double* true_kappa_sub = (double*)malloc(sizeof(double) * Njack);
    // for (int j = 0; j < Njack;j++) {
    //     true_kappa_sub[j] = head.kappa / (1. + 2. * fit_out.P[0][j] * head.kappa);
    // }
    // fprintf(outfile, "%.12g   %.12g\n", true_kappa_sub[Njack - 1], myres->comp_error(true_kappa_sub));
    // printf("true_kappa_sub=%.12g   %.12g\n", true_kappa_sub[Njack - 1], myres->comp_error(true_kappa_sub));
    // free_fit_result(fit_info, fit_out);
    // fit_info.restore_default();


}