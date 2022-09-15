#define common_integral_eq_C
#include <complex.h>
#include <vector>
#include "mutils.hpp"
#include "global.hpp"
#include "tower.hpp"
#include "resampling.hpp"
#include "fit_all.hpp"
#include "non_linear_fit.hpp"
#include "integral_eq_QC3.hpp"

using namespace std::complex_literals;

data_all setup_data_for_fits(int NE, int Njack, double*** M, double*** F) {
    data_all jackall;
    jackall.ens = NE;
    jackall.resampling = "jack";
    jackall.en = new data_single[jackall.ens];
    for (int i = 0;i < NE;i++) {
        jackall.en[i].Njack = Njack;
        jackall.en[i].Nobs = 4;
        jackall.en[i].resampling = jackall.resampling;
        jackall.en[i].jack = double_malloc_2(jackall.en[i].Nobs, Njack);
        for (int j = 0;j < Njack;j++) {
            jackall.en[i].jack[0][j] = M[i][0][j];
            jackall.en[i].jack[1][j] = M[i][1][j];
            jackall.en[i].jack[2][j] = F[i][0][j];
            jackall.en[i].jack[3][j] = F[i][1][j];

        }
    }
    return jackall;
}

void write_M3(int NE, int Njack, std::vector<double>& E3, double*** M3, double*** F, std::string namef, const char* prefix) {
    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", prefix, namef.c_str());
    FILE* f = open_file(namefile, "w+");
    fprintf(f, "%d\n", NE);
    fprintf(f, "%d\n", Njack);
    for (int i = 0;i < NE;i++) {
        for (int j = 0;j < Njack;j++) {
            // printf("%.12g  %.12g  %.12g\n", E3[i], M3[i][0][j], M3[i][1][j]);
            fprintf(f, "%.12g  %.12g  %.12g  %.12g  %.12g\n", E3[i], M3[i][0][j], M3[i][1][j], F[i][0][j], F[i][1][j]);
        }
    }
    fclose(f);

}

void read_M3(int& NE, int& Njack, std::vector<double>& E3, double***& M3, double***& F, std::string namef, const char* prefix) {
    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", prefix, namef.c_str());
    FILE* f = open_file(namefile, "r+");
    int is = 0;
    is += fscanf(f, "%d\n", &NE);
    is += fscanf(f, "%d\n", &Njack);
    printf("reading NE=%d  Njack=%d\n", NE, Njack);

    M3 = double_malloc_3(NE, 2, Njack);
    F = double_malloc_3(NE, 2, Njack);
    for (int i = 0;i < NE;i++) {
        for (int j = 0;j < Njack;j++) {
            double a, b, c, d;
            is += fscanf(f, "%lf %lf  %lf %lf  %lf\n", &E3[i], &a, &b, &c, &d);
            M3[i][0][j] = a;
            M3[i][1][j] = b;
            F[i][0][j] = c;
            F[i][1][j] = d;

            // printf("%d %d  %g %g\n",i,j,M3[i][0][j], M3[i][1][j]);
        }
    }
    error(is != 5 * Njack * NE + 2, 1, "read_M3", "numer of elements read = %d not correct, expected=%d", is, 2 * Njack * NE + 3);
    fclose(f);

}

void compute_M3(int NE, double Emin, double dE, int Njack, std::vector<double>& E3, double***& M3, double***& F, int N, int Npar, double** P,
    double compute_kcot(int, double*, int, double*), double** PKiso, double compute_kiso(double, double*), double eps) {

    M3 = double_malloc_3(NE, 2, Njack);
    F = double_malloc_3(NE, 2, Njack);
    // #pragma omp parallel for   shared(NE,Njack,Emin,dE,Npar,P,eps) 
    for (int i = 0;i < NE;i++) {
        #pragma omp parallel for   shared(NE,Njack,Emin,dE,Npar,P,eps) 
        for (int j = 0;j < Njack;j++) {
            E3[i] = Emin + i * dE;
            // E3[i] = 3.02;printf("\n\n MODIFY HERE \n\n");


            // std::complex<double> m3 = compute_M3_sym(E3[i], N, Npar, P[j], compute_kcot, PKiso[j], compute_kiso, eps);
            

            // std::complex<double> Kdf = compute_kiso(E3[i], PKiso[j]);
            Eigen::MatrixXcd D = compute_D(E3[i], N, Npar, P[j], compute_kcot, eps);
            std::complex<double> Finf = comput_Finf(E3[i], D, N, Npar, P[j], compute_kcot, eps);
            F[i][0][j] = Finf.real(); F[i][1][j] = Finf.imag();
            Eigen::VectorXcd  L = cal_L1(E3[i], D, N, Npar, P[j], compute_kcot, eps);


            std::complex<double> kiso = compute_kiso(E3[i], PKiso[j]);
            std::complex<double> m3 = std::complex<double>(0, 0);
            for (int i = 0;i < 3; i++) {
                for (int j = 0;j < 3; j++) {
                    // std::cout << j <<"   " << L(j)<<"   " << (1. / (1. / kiso + Finf))<<   std::endl;
                    m3 += L(i) * (1. / (1. / kiso + Finf)) * L(j);
                    // M3 += (1. / (1. / kiso + Finf)) ;
                }
            }
            M3[i][0][j] = m3.real(); M3[i][1][j] = m3.imag();
            // printf("jack =%-4d%-18.8g%-14g%-18g%-14g%-18g||%-25g%-25g%-25g%-25g\n",
            //     j, E3[i], real(m3), imag(m3), real(Kdf), imag(Kdf), PKiso[Njack-1][0], PKiso[Njack-1][1], PKiso[Njack-1][2], P[Njack-1][3]);
            // printf("%-18.8g%-14g%-18g%-14g%-18g%-18.12g%-20.12g\n", E3[i], real(m3), imag(m3), real(Kdf), imag(Kdf), real(Finf), imag(Finf));
            // printf("%-18.8g%-14g%-18g%\n", E3[i], real(m3), imag(m3));
        }

        // printf("EMFP:%-20.8g%-20.12g%-18g%-20.12g%-18g%-20.12g%-18g%-22.12g%-20g\n", E3[i], M3[i][0][Njack - 1], error_jackboot("jack", Njack, M3[i][0]),
        //     M3[i][1][Njack - 1], error_jackboot("jack", Njack, M3[i][1]),
        //     F[i][0][Njack - 1], error_jackboot("jack", Njack, F[i][0]), F[i][1][Njack - 1], error_jackboot("jack", Njack, F[i][1])
        // );
        // printf("%-10g%-10d%-10g%-20.12g%-18g%-22.12g%-20g\n",E3[i],N,eps, F[i][0][Njack - 1], error_jackboot("jack", Njack, F[i][0]), F[i][1][Njack - 1], error_jackboot("jack", Njack, F[i][1]));
        // exit(0);
    }
}

double rhs_laurent_pole(int n, int Nvar, double* x, int Npar, double* P) {
    error(Npar % 2 != 0, 1, "rhs_laurent_pole:", "Npar=%d but it must be multiple of two since the parameters are complex", Npar);
    std::complex<double> z(x[0], x[1]);

    std::complex<double> p(P[0], P[1]);
    std::complex<double> am1(P[2], P[3]);

    std::complex<double> r = am1 / (z * z - p * p);

    if (Npar >= 6) {
        std::complex<double> a0(P[4], P[5]);
        r += a0;
    }
    r = 1. / r;

    if (n == 0)      return real(r);
    else if (n == 1) return imag(r);
    else { printf("%s\n", __func__);  exit(1); }
}


double rhs_BW(int n, int Nvar, double* x, int Npar, double* P) {
    error(Npar % 2 != 0, 1, "rhs_laurent_pole:", "Npar=%d but it must be multiple of two since the parameters are complex", Npar);
    std::complex<double> E(x[0], x[1]);

    std::complex<double> r = (P[2] + 1i * P[3]) / (E - P[0] + 1i * P[1] / 2.0);

    if (Npar >= 5) {
        r += P[4] + 1i * P[5];
    }
    r = 1. / r;

    if (n == 0)      return real(r);
    else if (n == 1) return imag(r);
    else { printf("%s\n", __func__);  exit(1); }
}


double rhs_absBW(int n, int Nvar, double* x, int Npar, double* P) {
    error(Npar % 2 != 0, 1, "rhs_laurent_pole:", "Npar=%d but it must be multiple of two since the parameters are complex", Npar);
    std::complex<double> E(x[0], x[1]);

    std::complex<double> r = P[2] / (E - P[0] + 1i * P[1] / 2.0);

    if (Npar >= 4) {
        r += P[3];
    }
    r = 1. / r;

    if (n == 0)      return std::abs(r);
    else { printf("%s\n", __func__);  exit(1); }
}

double rhs_F(int n, int Nvar, double* x, int Npar, double* P) {

    std::complex<double> F = P[0] + 1i * P[1] + (P[2] + 1i * P[3]) * (x[0] + 1i * x[1]) * (x[0] + 1i * x[1]);
    if (n == 0)      return F.real();
    else if (n == 1) return F.imag();
    else { printf("%s\n", __func__);  exit(1); }
}


double compute_kiso(double E3_m, double* P) {
    return -P[0] / (E3_m * E3_m - P[1]) + P[2];
}
std::complex<double> compute_kiso_complex(std::complex<double> E3_m, double* P) {
    return -P[0] / (E3_m * E3_m - P[1]) + P[2];
}

double compute_kcot(int Nvar, double* x, int Npar, double* P) {
    double r;
    r = 1. / P[0];

    return r;
}
double denom_M(int n, int Nvar, double* x, int Npar, double* P) {


    //to minimize we need change parameters with variables
    double* Pf = (double*)malloc(sizeof(double) * 4);
    double* xf = (double*)malloc(sizeof(double) * 2);
    Pf[0] = x[0];
    Pf[1] = x[1];
    Pf[2] = x[2];
    Pf[3] = x[3];

    xf[0] = P[0];
    xf[1] = P[1];
    double Fre = rhs_F(0, Nvar, xf, Npar, Pf);
    double Fim = rhs_F(1, Nvar, xf, Npar, Pf);
    std::complex<double> F(Fre, Fim);

    double* PK = (double*)malloc(sizeof(double) * 4);
    double* xK = (double*)malloc(sizeof(double) * 1);
    PK[0] = x[4];
    PK[1] = x[5];
    PK[2] = x[6];
    xK[0] = P[0];

    std::complex<double> K = compute_kiso_complex(std::complex<double>(P[0], P[1]), PK);
    // printf("F=%g  %g   E=%g   pF=%g   %g  %g  %g\n", real(F ), imag(F), P[0], Pf[0],Pf[1],Pf[2],Pf[3]);
    F = F + 1. / K;
    // printf("D=%g   E=%g   F=%g   %g  K=%g\n", real(F * conj(F)), P[0], F.real(), F.imag(),K);
    free(xK);free(PK);free(xf);free(Pf);

    return sqrt(real(F * conj(F)));
}

double lhs_M3(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;

    std::complex<double> M3(gjack.en[e].jack[0][j], gjack.en[e].jack[1][j]);
    M3 = 1.0 / M3;

    if (n == 0) {
        r = M3.real();
    }
    else if (n == 1) {
        r = M3.imag();
    }
    else {
        r = 0;  printf("lhs_M3 n=%d not implemented\n", n); exit(1);
    }
    return  r;
}


double lhs_absM3(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;

    std::complex<double> M3(gjack.en[e].jack[0][j], gjack.en[e].jack[1][j]);
    M3 = 1.0 / M3;

    if (n == 0) {
        r = std::abs(M3);
    }

    else {
        r = 0;  printf("lhs_M3 n=%d not implemented\n", n); exit(1);
    }
    return  r;
}



double lhs_F(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;

    std::complex<double> F(gjack.en[e].jack[2][j], gjack.en[e].jack[3][j]);


    if (n == 0) {
        r = F.real();
    }
    else if (n == 1) {
        r = F.imag();
    }
    else {
        r = 0;  printf("lhs_F3 n=%d not implemented\n", n); exit(1);
    }
    return  r;
}

