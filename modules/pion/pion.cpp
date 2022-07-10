#define pion_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"

#include "resampling.hpp"
#include "read.hpp"
//#include "m_eff.hpp"
#include "gnuplot.hpp"
//#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "rand.hpp"
#include "non_linear_fit.hpp"
#include "pion.hpp"



double Mw2_fw(int n, int Nvar, double* x, int Npar, double* P) {

    double Mw2 = 0, xi;
    double pi = 3.141592653589793;
    double Bw = P[0], fw = P[1], P1 = P[2], P2 = P[3], P3 = P[4], P4 = P[5];
    double mw = x[0], w0 = x[1];// KM=x[2],Kf=x[3];

    xi = 2 * Bw * mw / (16. * pi * pi * fw * fw);

    if (n == 0) {

        Mw2 = 1 + xi * log(xi) + P1 * xi + (1. / (w0 * w0)) * P2;
        Mw2 *= 2 * Bw * mw;

    }
    if (n == 1) {

        Mw2 = 1 - 2 * xi * log(xi) + P3 * xi + (1. / (w0 * w0)) * P4;
        Mw2 *= fw;


    }
    return Mw2;
}


double Mw2_fw_a0_minus(double x, int Npar, double* P) {

    int n = 0;//Mw2


    double Mw2 = 0, xi;
    double pi = 3.141592653589793;
    double Bw = P[0], fw = P[1], P1 = P[2], P2 = P[3], P3 = P[4], P4 = P[5];
    double mw = x;// KM=x[2],Kf=x[3];

    xi = 2 * Bw * mw / (16. * pi * pi * fw * fw);

    if (n == 0) {

        Mw2 = 1 + xi * log(xi) + P1 * xi;
        Mw2 *= 2 * Bw * mw;
        Mw2 -= pow(134.98 * 0.1714 / 197.326963, 2);

    }
    if (n == 1) {

        Mw2 = 1 - 2 * xi * log(xi) + P3 * xi;
        Mw2 *= fw;
        Mw2 -= 130.41 * 0.1714 / 197.326963;

    }
    return Mw2;
}



double Mw2_fw_chiral_c2(int n, int Nvar, double* x, int Npar, double* P) {

    double Mw2 = 0, xi;
    double pi = 3.141592653589793;
    double Bw = P[0], fw = P[1], P1 = P[2], P2 = P[3], P3 = P[4], P4 = P[5], c2ww = P[6];
    double mw = x[0], w0 = x[1];// KM=x[2],Kf=x[3];

    xi = 2 * Bw * mw / (16. * pi * pi * fw * fw);

    if (n == 0) {

        Mw2 = 1 + xi * log(xi) + P1 * xi + (1. / (w0 * w0)) * (P2 + 4 * c2ww * log(xi) / (4 * pi * 4 * pi * fw * fw));
        Mw2 *= 2 * Bw * mw;
    }
    if (n == 1) {

        Mw2 = 1 - 2 * xi * log(xi) + P3 * xi + (1. / (w0 * w0)) * (P4 - 4 * c2ww * log(xi) / (4 * pi * 4 * pi * fw * fw));
        Mw2 *= fw;
    }
    return Mw2;
}

double Mw2_fw_chiral_c2_a0_minus(double x, int Npar, double* P) {

    int n = 0;//Mw2

    double Mw2 = 0, xi;
    double pi = 3.141592653589793;
    double Bw = P[0], fw = P[1], P1 = P[2], P2 = P[3], P3 = P[4], P4 = P[5], c2ww = P[6];
    double mw = x;// KM=x[2],Kf=x[3];

    xi = 2 * Bw * mw / (16. * pi * pi * fw * fw);

    if (n == 0) {

        Mw2 = 1 + xi * log(xi) + P1 * xi;
        Mw2 *= 2 * Bw * mw;
        Mw2 -= pow(134.98 * 0.1714 / 197.326963, 2);

    }
    if (n == 1) {

        Mw2 = 1 - 2 * xi * log(xi) + P3 * xi;
        Mw2 *= fw;
        Mw2 -= 130.41 * 0.1714 / 197.326963;

    }
    return Mw2;
}



double m_over_f_xi(double xi, int Npar, double* P) {
    double Mw2 = 0;
    double pi = 3.141592653589793;
    double Bw = P[0], fw = P[1], P1 = P[2], P2 = P[3], P3 = P[4], P4 = P[5];
    double tmp;

    Mw2 = 1 + xi * log(xi) + P1 * xi;
    tmp = 1 - 2 * xi * log(xi) + P3 * xi;
    tmp = 16 * pi * pi * xi * Mw2 / (tmp * tmp);

    return tmp - 1.07131468442;//value of M_pi^2/f^2

}

double m_over_f(double x, int Npar, double* P) {
    double Bw = P[0], fw = P[1], P1 = P[2], P2 = P[3], P3 = P[4], P4 = P[5];

    double xi = 2. * Bw * x / (16 * pi_greco * pi_greco * fw * fw);


    return m_over_f_xi(xi, Npar, P);

}






double Mw2_fw_polynomial(int n, int Nvar, double* x, int Npar, double* P) {

    double Mw2 = 0, xi;
    double pi = 3.141592653589793;
    double Bw = P[0], fw = P[1], P1 = P[2], P2 = P[3], P3 = P[4], P4 = P[5], P5 = P[6], P6 = P[7];
    double mw = x[0], w0 = x[1];// KM=x[2],Kf=x[3];

    //xi=2*Bw*ml*w0/(16.*pi*pi*fw*fw);

    if (n == 0) {

        Mw2 = 1 + P1 * mw + (1. / (w0 * w0)) * P2 + P3 * mw * mw;
        Mw2 *= 2 * Bw * mw;

    }
    if (n == 1) {

        Mw2 = 1 + P4 * mw + (1. / (w0 * w0)) * P5 + P6 * mw * mw;
        Mw2 *= fw;

    }
    return Mw2;
}



double m_over_f_pol(double mw, int Npar, double* P) {
    double Bw = P[0], fw = P[1], P1 = P[2], P2 = P[3], P3 = P[4], P4 = P[5], P5 = P[6], P6 = P[7];
    double r = 0;
    double pi = 3.141592653589793;

    r = 2 * Bw * mw * (1 + P1 * mw + P3 * mw * mw);
    r /= (fw * (1 + P4 * mw + P6 * mw * mw));
    r /= (fw * (1 + P4 * mw + P6 * mw * mw));


    return r - 1.07131468442;


}



double Mw2_fw_polynomial_6(int n, int Nvar, double* x, int Npar, double* P) {

    double Mw2 = 0, xi;
    double pi = 3.141592653589793;
    double Bw = P[0], fw = P[1], P1 = P[2], P2 = P[3], P4 = P[4], P5 = P[5];
    double mw = x[0], w0 = x[1];// KM=x[2],Kf=x[3];

    //xi=2*Bw*ml*w0/(16.*pi*pi*fw*fw);

    if (n == 0) {

        Mw2 = 1 + P1 * mw + (1. / (w0 * w0)) * P2;
        Mw2 *= 2 * Bw * mw;

    }
    if (n == 1) {

        Mw2 = 1 + P4 * mw + (1. / (w0 * w0)) * P5;
        Mw2 *= fw;

    }
    return Mw2;
}
double m_over_f_pol_6(double mw, int Npar, double* P) {
    double Bw = P[0], fw = P[1], P1 = P[2], P2 = P[3], P4 = P[4], P5 = P[5]
        ;
    double r = 0;
    double pi = 3.141592653589793;

    r = 2 * Bw * mw * (1 + P1 * mw);
    r /= (fw * (1 + P4 * mw));
    r /= (fw * (1 + P4 * mw));

    return r - 1.07131468442;


}



double** fit_Mpi_fw(struct database_file_jack* jack_files, struct header* head, int Njack, int*** mass_index, struct data_jack* gJ) {
    double*** y, ** x, ** r, * chi2, * tmp, * rm, * chi2m, ** fit;
    int i, j, e, im;
    int Npar = 6;
    int Nvar = 2;//m_l, w0,KM
    int ik1 = 0, ik2 = 0;
    int n, count, N = 2;
    int* en = (int*)malloc(sizeof(int) * N);
    en[0] = ensembles;
    en[1] = ensembles;
    int en_tot = 0;

    for (n = 0;n < N;n++)
        en_tot += en[n];
    double* guess = (double*)malloc(sizeof(double) * Npar);
    guess[0] = 2.05478;
    guess[1] = 0.113021;
    guess[2] = 2.54451;
    guess[3] = 0.410103;
    guess[4] = 4.61247;
    guess[5] = -0.272884;

    x = (double**)malloc(sizeof(double*) * (en_tot));

    chi2m = (double*)malloc(sizeof(double) * (Npar));
    rm = (double*)malloc(sizeof(double*) * (Njack));
    fit = (double**)malloc(sizeof(double*) * (en_tot));

    r = (double**)malloc(sizeof(double*) * (Npar));
    for (i = 0;i < Npar;i++) {
        r[i] = (double*)malloc(sizeof(double) * Njack);
    }

    chi2 = (double*)malloc(sizeof(double) * Njack);
    y = (double***)malloc(sizeof(double**) * Njack);

    for (j = 0;j < Njack;j++) {
        y[j] = (double**)malloc(sizeof(double*) * (en_tot));
        count = 0;
        for (n = 0;n < N;n++) {
            for (i = 0;i < en[n];i++) {
                y[j][i + count] = (double*)malloc(sizeof(double) * 2);
            }
            count += en[n];
        }
    }

    printf("w0/a[fm]     m*w0/aZp[]      (M_Pi w0/KM2)^2 or fw/Kf [MeV] err\n");
    count = 0;
    for (n = 0;n < N;n++) {
        printf("#function %d\n", n);
        for (e = 0;e < en[n];e++) {
            im = mass_index[e][ik2][ik1];
            if (n == 0) {
                for (j = 0;j < Njack;j++) {
                    rm[j] = gJ[e].M_PS_GEVP_jack[im][j] * gJ[e].w0[j] / (gJ[e].KM[im]);
                    rm[j] *= rm[j];
                }
                fit[e + count] = mean_and_error(jack_files[0].sampling, Njack, rm);
            }
            if (n == 1) {
                for (j = 0;j < Njack;j++) {
                    rm[j] = gJ[e].f_PS_jack[im][j] * gJ[e].w0[j] / gJ[e].Kf[im];
                }
                fit[e + count] = mean_and_error(jack_files[0].sampling, Njack, rm);
            }

            for (j = 0;j < jack_tot;j++) {
                y[j][e + count][0] = rm[j];
                y[j][e + count][1] = fit[e + count][1];
            }


            x[e + count] = (double*)malloc(sizeof(double) * Nvar);
            x[e + count][0] = head[e].k[head[e].nk + ik2] * gJ[e].w0[Njack - 1] / gJ[e].Zp[Njack - 1];//ml*w0
            x[e + count][1] = gJ[e].w0[Njack - 1];//w0
            //x[e+count][2]=gJ[e].KM[im];//KM
            //x[e+count][3]=gJ[e].Kf[im];//Kf
            printf("%g     %g      %g   %g\n", x[e + count][1], x[e + count][0], fit[e + count][0], fit[e + count][1]);
        }
        count += en[n];
    }



    for (j = 0;j < Njack;j++) {
        tmp = non_linear_fit_Nf(N, en, x, y[j], Nvar, Npar, Mw2_fw, guess).P;

        chi2[j] = compute_chi_non_linear_Nf(N, en, x, y[j], tmp, Nvar, Npar, Mw2_fw);

        for (i = 0;i < Npar;i++) {
            r[i][j] = tmp[i];
        }
        free(tmp);

    }

    chi2m = mean_and_error(jack_files[0].sampling, Njack, chi2);
    printf("$\\chi^2=%f+-%f$\n", chi2m[0], chi2m[1]);
    free(rm);free(chi2m);

    // make_plots_MPi_fPi(en,);

    for (e = 0;e < en_tot;e++) {
        free(x[e]);   free(fit[e]);

    }

    free(fit);

    free(x);
    for (j = 0;j < Njack;j++) {
        for (e = 0;e < en_tot;e++) {
            free(y[j][e]);
        }
        free(y[j]);
    }
    free(y); free(guess);
    return r;

}

double** fit_Mpi_fw_chiral_c2(struct database_file_jack* jack_files, struct header* head, int Njack, int*** mass_index, struct data_jack* gJ) {
    double*** y, ** x, ** r, * chi2, * tmp, * rm, * chi2m, ** fit;
    int i, j, e, im;
    int Npar = 7;
    int Nvar = 2;//m_l, w0,KM
    int ik1 = 0, ik2 = 0;
    int n, count, N = 2;
    int* en = (int*)malloc(sizeof(int) * N);
    en[0] = ensembles;
    en[1] = ensembles;
    int en_tot = 0;

    for (n = 0;n < N;n++)
        en_tot += en[n];

    double* guess = (double*)malloc(sizeof(double) * Npar);
    guess[0] = 2.05478;
    guess[1] = 0.113021;
    guess[2] = 2.54451;
    guess[3] = 0.410103;
    guess[4] = 4.61247;
    guess[5] = -0.272884;
    guess[6] = -0.0412976;
    x = (double**)malloc(sizeof(double*) * (en_tot));

    chi2m = (double*)malloc(sizeof(double) * (Npar));
    rm = (double*)malloc(sizeof(double*) * (Njack));
    fit = (double**)malloc(sizeof(double*) * (en_tot));

    r = (double**)malloc(sizeof(double*) * (Npar));
    for (i = 0;i < Npar;i++) {
        r[i] = (double*)malloc(sizeof(double) * Njack);
    }

    chi2 = (double*)malloc(sizeof(double) * Njack);
    y = (double***)malloc(sizeof(double**) * Njack);

    for (j = 0;j < Njack;j++) {
        y[j] = (double**)malloc(sizeof(double*) * (en_tot));
        count = 0;
        for (n = 0;n < N;n++) {
            for (i = 0;i < en[n];i++) {
                y[j][i + count] = (double*)malloc(sizeof(double) * 2);
            }
            count += en[n];
        }
    }

    printf("w0/a[fm]     m*w0/aZp[]      (M_Pi w0/KM2)^2 or fw/Kf [MeV] err\n");
    count = 0;
    for (n = 0;n < N;n++) {
        printf("#function %d\n", n);
        for (e = 0;e < en[n];e++) {
            im = mass_index[e][ik2][ik1];
            if (n == 0) {
                for (j = 0;j < Njack;j++) {
                    rm[j] = gJ[e].M_PS_GEVP_jack[im][j] * gJ[e].w0[j] / (gJ[e].KM[im]);
                    rm[j] *= rm[j];
                }
                fit[e + count] = mean_and_error(jack_files[0].sampling, Njack, rm);
            }
            if (n == 1) {
                for (j = 0;j < Njack;j++) {
                    rm[j] = gJ[e].f_PS_jack[im][j] * gJ[e].w0[j] / gJ[e].Kf[im];
                }
                fit[e + count] = mean_and_error(jack_files[0].sampling, Njack, rm);
            }

            for (j = 0;j < jack_tot;j++) {
                y[j][e + count][0] = rm[j];
                y[j][e + count][1] = fit[e + count][1];
            }


            x[e + count] = (double*)malloc(sizeof(double) * Nvar);
            x[e + count][0] = head[e].k[head[e].nk + ik2] * gJ[e].w0[Njack - 1] / gJ[e].Zp[Njack - 1];//ml*w0
            x[e + count][1] = gJ[e].w0[Njack - 1];//w0
            //x[e+count][2]=gJ[e].KM[im];//KM
            //x[e+count][3]=gJ[e].Kf[im];//Kf
            printf("%g     %g      %g   %g\n", x[e + count][1], x[e + count][0], fit[e + count][0], fit[e + count][1]);
        }
        count += en[n];
    }



    for (j = 0;j < Njack;j++) {
        tmp = non_linear_fit_Nf(N, en, x, y[j], Nvar, Npar, Mw2_fw_chiral_c2, guess).P;

        chi2[j] = compute_chi_non_linear_Nf(N, en, x, y[j], tmp, Nvar, Npar, Mw2_fw_chiral_c2);

        for (i = 0;i < Npar;i++) {
            r[i][j] = tmp[i];
        }
        free(tmp);

    }

    chi2m = mean_and_error(jack_files[0].sampling, Njack, chi2);
    printf("$\\chi^2=%f+-%f$\n", chi2m[0], chi2m[1]);
    free(rm);free(chi2m);

    // make_plots_MPi_fPi(en,);

    for (e = 0;e < en_tot;e++) {
        free(x[e]);   free(fit[e]);

    }

    free(fit);

    free(x);
    for (j = 0;j < Njack;j++) {
        for (e = 0;e < en_tot;e++) {
            free(y[j][e]);
        }
        free(y[j]);
    }
    free(y); free(guess);
    return r;

}




double** fit_Mpi_fw_polynomial(struct database_file_jack* jack_files, struct header* head, int Njack, int*** mass_index, struct data_jack* gJ) {
    double*** y, ** x, ** r, * chi2, * tmp, * rm, * chi2m, ** fit;
    int i, j, e, im;
    int Npar = 8;
    int Nvar = 2;//m_l, w0,KM
    int ik1 = 0, ik2 = 0;
    int n, count, N = 2;
    int* en = (int*)malloc(sizeof(int) * N);
    en[0] = ensembles;
    en[1] = ensembles;
    int en_tot = 0;

    for (n = 0;n < N;n++)
        en_tot += en[n];
    double* guess = (double*)malloc(sizeof(double) * Npar);
    guess[0] = 0.1;
    guess[1] = 0.1;
    guess[2] = 0.1;
    guess[3] = 0.1;
    guess[4] = 0.1;
    guess[5] = 0.1;
    guess[6] = 0.1;
    guess[7] = 0.1;


    x = (double**)malloc(sizeof(double*) * (en_tot));

    chi2m = (double*)malloc(sizeof(double) * (Npar));
    rm = (double*)malloc(sizeof(double*) * (Njack));
    fit = (double**)malloc(sizeof(double*) * (en_tot));

    r = (double**)malloc(sizeof(double*) * (Npar));
    for (i = 0;i < Npar;i++) {
        r[i] = (double*)malloc(sizeof(double) * Njack);
    }

    chi2 = (double*)malloc(sizeof(double) * Njack);
    y = (double***)malloc(sizeof(double**) * Njack);

    for (j = 0;j < Njack;j++) {
        y[j] = (double**)malloc(sizeof(double*) * (en_tot));
        count = 0;
        for (n = 0;n < N;n++) {
            for (i = 0;i < en[n];i++) {
                y[j][i + count] = (double*)malloc(sizeof(double) * 2);
            }
            count += en[n];
        }
    }

    // printf("w0/a[fm]     am/aZp[]      (M_Pi w0/KM2)^2 or fw/Kf [MeV] err\n");
    count = 0;
    for (n = 0;n < N;n++) {
        //    printf("#function %d\n",n);
        for (e = 0;e < en[n];e++) {
            im = mass_index[e][ik2][ik1];
            if (n == 0) {
                for (j = 0;j < Njack;j++) {
                    rm[j] = gJ[e].M_PS_GEVP_jack[im][j] * gJ[e].w0[j] / (gJ[e].KM[im] * gJ[e].KM[im]);
                    rm[j] *= rm[j];
                }
                fit[e + count] = mean_and_error(jack_files[0].sampling, Njack, rm);
            }
            if (n == 1) {
                for (j = 0;j < Njack;j++) {
                    rm[j] = gJ[e].f_PS_jack[im][j] * gJ[e].w0[j] / gJ[e].Kf[im];
                }
                fit[e + count] = mean_and_error(jack_files[0].sampling, Njack, rm);
            }

            for (j = 0;j < jack_tot;j++) {
                y[j][e + count][0] = rm[j];
                y[j][e + count][1] = fit[e + count][1];
            }


            x[e + count] = (double*)malloc(sizeof(double) * Nvar);
            x[e + count][0] = head[e].k[head[e].nk + ik2] * gJ[e].w0[Njack - 1] / gJ[e].Zp[Njack - 1];//ml
            x[e + count][1] = gJ[e].w0[Njack - 1];//w0
            //x[e+count][2]=gJ[e].KM[im];//KM
            //x[e+count][3]=gJ[e].Kf[im];//Kf
  //          printf("%g     %g      %g   %g\n",x[e+count][1],x[e+count][0],fit[e+count][0],fit[e+count][1]);
        }
        count += en[n];
    }



    for (j = 0;j < Njack;j++) {
        tmp = non_linear_fit_Nf(N, en, x, y[j], Nvar, Npar, Mw2_fw_polynomial, guess).P;

        chi2[j] = compute_chi_non_linear_Nf(N, en, x, y[j], tmp, Nvar, Npar, Mw2_fw_polynomial);

        for (i = 0;i < Npar;i++) {
            r[i][j] = tmp[i];
        }
        free(tmp);

    }

    chi2m = mean_and_error(jack_files[0].sampling, Njack, chi2);
    printf("$\\chi^2=%f+-%f$\n", chi2m[0], chi2m[1]);
    free(rm);free(chi2m);

    // make_plots_MPi_fPi(en,);

    for (e = 0;e < en_tot;e++) {
        free(x[e]);   free(fit[e]);

    }

    free(fit);

    free(x);
    for (j = 0;j < Njack;j++) {
        for (e = 0;e < en_tot;e++) {
            free(y[j][e]);
        }
        free(y[j]);
    }
    free(y); free(guess);
    return r;

}




double** fit_Mpi_fw_polynomial_6(struct database_file_jack* jack_files, struct header* head, int Njack, int*** mass_index, struct data_jack* gJ) {
    double*** y, ** x, ** r, * chi2, * tmp, * rm, * chi2m, ** fit;
    int i, j, e, im;
    int Npar = 6;
    int Nvar = 2;//m_l, w0,KM
    int ik1 = 0, ik2 = 0;
    int n, count, N = 2;
    int* en = (int*)malloc(sizeof(int) * N);
    en[0] = ensembles;
    en[1] = ensembles;
    int en_tot = 0;
    FILE* temp = fopen("data.temp", "w");
    FILE* gnuplotPipe = popen("gnuplot ", "w");
    char name[500], instruction[1000];
    double* m_pi;

    for (n = 0;n < N;n++)
        en_tot += en[n];

    double* guess = (double*)malloc(sizeof(double) * Npar);
    guess[0] = 0.1;
    guess[1] = 0.1;
    guess[2] = 0.1;
    guess[3] = 0.1;
    guess[4] = 0.1;
    guess[5] = 0.1;


    x = (double**)malloc(sizeof(double*) * (en_tot));

    chi2m = (double*)malloc(sizeof(double) * (Npar));
    rm = (double*)malloc(sizeof(double*) * (Njack));
    fit = (double**)malloc(sizeof(double*) * (en_tot));

    r = (double**)malloc(sizeof(double*) * (Npar));
    for (i = 0;i < Npar;i++) {
        r[i] = (double*)malloc(sizeof(double) * Njack);
    }

    chi2 = (double*)malloc(sizeof(double) * Njack);
    y = (double***)malloc(sizeof(double**) * Njack);

    for (j = 0;j < Njack;j++) {
        y[j] = (double**)malloc(sizeof(double*) * (en_tot));
        count = 0;
        for (n = 0;n < N;n++) {
            for (i = 0;i < en[n];i++) {
                y[j][i + count] = (double*)malloc(sizeof(double) * 2);
            }
            count += en[n];
        }
    }

    // printf("w0/a[fm]     am/aZp[]      (M_Pi w0/KM2)^2 or fw/Kf [MeV] err\n");
    count = 0;
    for (n = 0;n < N;n++) {
        //    printf("#function %d\n",n);
        for (e = 0;e < en[n];e++) {
            im = mass_index[e][ik2][ik1];
            if (n == 0) {
                for (j = 0;j < Njack;j++) {
                    rm[j] = gJ[e].M_PS_GEVP_jack[im][j] * gJ[e].w0[j] / (gJ[e].KM[im] * gJ[e].KM[im]);
                    rm[j] *= rm[j];
                }
                fit[e + count] = mean_and_error(jack_files[0].sampling, Njack, rm);
            }
            if (n == 1) {
                for (j = 0;j < Njack;j++) {
                    rm[j] = gJ[e].f_PS_jack[im][j] * gJ[e].w0[j] / gJ[e].Kf[im];
                }
                fit[e + count] = mean_and_error(jack_files[0].sampling, Njack, rm);
            }

            for (j = 0;j < jack_tot;j++) {
                y[j][e + count][0] = rm[j];
                y[j][e + count][1] = fit[e + count][1];
            }


            x[e + count] = (double*)malloc(sizeof(double) * Nvar);
            x[e + count][0] = head[e].k[head[e].nk + ik2] * gJ[e].w0[Njack - 1] / gJ[e].Zp[Njack - 1];//ml
            x[e + count][1] = gJ[e].w0[Njack - 1];//w0

            //  printf("%g     %g      %g   %g\n",x[e+count][1],x[e+count][0],fit[e+count][0],fit[e+count][1]);
    /*          if(n==0)
                  fprintf(temp,"%g \t %.15g \t %.15g\n",x[e+count][1]*x[e+count][0],fit[e+count][0]/(x[e+count][0]*x[e+count][1]),fit[e+count][1]/(x[e+count][0]*x[e+count][1]));
              if(n==1)
                  fprintf(temp,"%g \t %.15g \t %.15g\n",x[e+count][1]*x[e+count][0],fit[e+count][0],fit[e+count][1]);*/
        }
        fprintf(temp, "\n\n");
        count += en[n];
    }
    /*fclose(temp);
     fprintf(gnuplotPipe,"set sample 200\n");
     fprintf(gnuplotPipe,"set key right \n");
     fprintf(gnuplotPipe,"set term epslatex  color colortext standalone size 16.5cm,11.5cm linewidth 2 \n");
     fprintf(gnuplotPipe,"set linetype 1 lw 1.5 pt 5 ps 2 lc 'blue'\n");
     fprintf(gnuplotPipe,"set linetype 2 lw 1.5 pt 7 ps 2 lc 'red'\n");
     fprintf(gnuplotPipe,"set mytics 2\n");
     fprintf(gnuplotPipe,"set key spacing 1.2\n");
     */
    for (j = 0;j < Njack;j++) {
        tmp = non_linear_fit_Nf(N, en, x, y[j], Nvar, Npar, Mw2_fw_polynomial_6, guess).P;

        chi2[j] = compute_chi_non_linear_Nf(N, en, x, y[j], tmp, Nvar, Npar, Mw2_fw_polynomial_6);

        for (i = 0;i < Npar;i++) {
            r[i][j] = tmp[i];
        }
        free(tmp);

    }

    chi2m = mean_and_error(jack_files[0].sampling, Njack, chi2);
    printf("$\\chi^2=%f+-%f$\n", chi2m[0], chi2m[1]);
    free(rm);free(chi2m);

    /*
    sprintf(name,"M_pi_polynomial6");
    fprintf(gnuplotPipe,"set xlabel '$w_0 m_l$'\n");
    fprintf(gnuplotPipe,"set ylabel '$w_0 M_{\\pi}^2/m_l$'\n");
    fprintf(gnuplotPipe,"set output '%s.tex'\n",name);
    fprintf(gnuplotPipe,"f(x,y)=2*Bw*x*(1+P1*x+ (1./(y*y))*P2)\n");
    fprintf(gnuplotPipe,"g(x,y)=f(x,y)/(x)\n");
    fprintf(gnuplotPipe,"Bw=%g\nP1=%g\nP2=%g\n",r[0][Njack-1],r[2][Njack-1],r[3][Njack-1]);
       fprintf(gnuplotPipe,"plot [0.002:0.02] 'data.temp' i 0 u 1:2:3 w e t '$\\chi^2=%.2f\\pm %.f$', g(x,1.760)  t 'A', g(x,2.1330729)  t 'B', g(x,1e+10)  t 'continuum'\n",chi2m[0],chi2m[1]);
   // make_plots_MPi_fPi(en,);
     fprintf(gnuplotPipe,"exit \n");

     pclose(gnuplotPipe);
 //    sprintf(instruction,"pdflatex -output-directory=./ %s", name);
 //    sprintf(instruction,"pdflatex  %s", name);
     system(instruction);
     */
     //    sprintf(instruction,"rm %s.tex %s.log %s.aux  %s-inc.eps %s-inc-eps-converted-to.pdf",name,name,name,name,name);
    system(instruction);

    for (e = 0;e < en_tot;e++) {
        free(x[e]);   free(fit[e]);

    }

    free(fit);

    free(x);
    for (j = 0;j < Njack;j++) {
        for (e = 0;e < en_tot;e++) {
            free(y[j][e]);
        }
        free(y[j]);
    }
    free(y); free(guess);
    return r;

}









//w0=1/Mss
double** fit_Mpi_fw_Mss(struct database_file_jack* jack_files, struct header* head, int Njack, int*** mass_index, struct data_jack* gJ) {
    double*** y, ** x, ** r, * chi2, * tmp, * rm, * chi2m, ** fit;
    int i, j, e, im;
    int Npar = 6;
    int Nvar = 2;//m_l, w0,KM
    int ik1 = 0, ik2 = 0;
    int n, count, N = 2;
    int* en = (int*)malloc(sizeof(int) * N);
    en[0] = ensembles;
    en[1] = ensembles;
    int en_tot = 0;

    for (n = 0;n < N;n++)
        en_tot += en[n];
    double* guess = (double*)malloc(sizeof(double) * Npar);
    for (n = 0;n < Npar;n++)
        guess[n] = 0.1;

    x = (double**)malloc(sizeof(double*) * (en_tot));

    chi2m = (double*)malloc(sizeof(double) * (Npar));
    rm = (double*)malloc(sizeof(double*) * (Njack));
    fit = (double**)malloc(sizeof(double*) * (en_tot));

    r = (double**)malloc(sizeof(double*) * (Npar));
    for (i = 0;i < Npar;i++) {
        r[i] = (double*)malloc(sizeof(double) * Njack);
    }

    chi2 = (double*)malloc(sizeof(double) * Njack);
    y = (double***)malloc(sizeof(double**) * Njack);

    for (j = 0;j < Njack;j++) {
        y[j] = (double**)malloc(sizeof(double*) * (en_tot));
        count = 0;
        for (n = 0;n < N;n++) {
            for (i = 0;i < en[n];i++) {
                y[j][i + count] = (double*)malloc(sizeof(double) * 2);
            }
            count += en[n];
        }
    }

    printf("w0/a[fm]     m*w0/aZp[]      (M_Pi w0/KM2)^2 or fw/Kf [MeV] err\n");
    count = 0;
    for (n = 0;n < N;n++) {
        printf("#function %d\n", n);
        for (e = 0;e < en[n];e++) {
            im = mass_index[e][ik2][ik1];
            if (n == 0) {
                for (j = 0;j < Njack;j++) {
                    rm[j] = gJ[e].M_PS_GEVP_jack[im][j] * gJ[e].w0[j] / (gJ[e].KM[im]);
                    rm[j] *= rm[j];
                }
                fit[e + count] = mean_and_error(jack_files[0].sampling, Njack, rm);
            }
            if (n == 1) {
                for (j = 0;j < Njack;j++) {
                    rm[j] = gJ[e].f_PS_jack[im][j] * gJ[e].w0[j] / gJ[e].Kf[im];
                }
                fit[e + count] = mean_and_error(jack_files[0].sampling, Njack, rm);
            }

            for (j = 0;j < jack_tot;j++) {
                y[j][e + count][0] = rm[j];
                y[j][e + count][1] = fit[e + count][1];
            }


            x[e + count] = (double*)malloc(sizeof(double) * Nvar);
            x[e + count][0] = head[e].k[head[e].nk + ik2] * gJ[e].w0[Njack - 1] / gJ[e].Zp[Njack - 1];//ml*w0
            x[e + count][1] = gJ[e].w0[Njack - 1];//w0
            //x[e+count][2]=gJ[e].KM[im];//KM
            //x[e+count][3]=gJ[e].Kf[im];//Kf
            printf("%g     %g      %g   %g\n", x[e + count][1], x[e + count][0], fit[e + count][0], fit[e + count][1]);
        }
        count += en[n];
    }



    for (j = 0;j < Njack;j++) {
        tmp = non_linear_fit_Nf(N, en, x, y[j], Nvar, Npar, Mw2_fw, guess).P;

        chi2[j] = compute_chi_non_linear_Nf(N, en, x, y[j], tmp, Nvar, Npar, Mw2_fw);

        for (i = 0;i < Npar;i++) {
            r[i][j] = tmp[i];
        }
        free(tmp);

    }

    chi2m = mean_and_error(jack_files[0].sampling, Njack, chi2);
    printf("$\\chi^2=%f+-%f$\n", chi2m[0], chi2m[1]);
    free(rm);free(chi2m);

    // make_plots_MPi_fPi(en,);

    for (e = 0;e < en_tot;e++) {
        free(x[e]);   free(fit[e]);

    }

    free(fit);

    free(x);
    for (j = 0;j < Njack;j++) {
        for (e = 0;e < en_tot;e++) {
            free(y[j][e]);
        }
        free(y[j]);
    }
    free(y); free(guess);
    return r;

}


