#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "HLT.hpp"
#include "mutils.hpp"


#include "arb.h"
#include "acb_modular.h"
#include "acb_calc.h"


double no_smearing(double x, double* p) {
    return 1.0;
};

int main() {

    int T = 24;
    int tmax = 10;
    double E0 = 0.5;
    int Njack = 10;
    {
        HLT_type_d HLT_space(tmax, T, E0, Njack, HLT_EXP_b);

        printf("testing integration without smearing\n");
        wrapper_smearing_d Delta0(no_smearing, std::vector<double> {}, & HLT_space);
        HLT_space.compute_f_EXP_b(Delta0, 1e-10);
        for (int t = 0;t < HLT_space.Tmax;t++) {
            double analytic = exp(-0.5 * (t + 1)) / (t + 1) + exp(-0.5 * (T - (t + 1))) / (T - (t + 1));
            double diff = fabs(HLT_space.f[t] - analytic);
            printf("f[%d]=%.12g  analytic=%.12g  diff=%.12g\n", t, HLT_space.f[t], analytic, diff);
            error(diff > 1e-7, 1, "integration f not working", "t=%d diff=%.12g", t, diff);
        }
        printf("testing integration with smearing\n");
        std::vector<double>  gauss_p = { 0,1 };
        printf("gaussian mu=%.12g  sigma=%.12g\n", gauss_p[0], gauss_p[1]);
        wrapper_smearing_d Delta(gaussian_for_HLT_d, gauss_p, &HLT_space);
        HLT_space.compute_f_EXP_b(Delta, 1e-10);
        for (int t = 0;t < HLT_space.Tmax;t++) {
            double a = t + 1;
            double b = (T - (t + 1));
            double analytic = -0.5 * exp(0.5 * a * a) * (-1. + erf(M_SQRT1_2f64 / 2. + M_SQRT1_2f64 * (a)))
                - 0.5 * exp(0.5 * b * b) * (-1. + erf(M_SQRT1_2f64 / 2. + M_SQRT1_2 * b));
            double diff = fabs(HLT_space.f[t] - analytic);
            printf("f[%d]=%.12g  analytic=%.12g  diff=%.12g\n", t, HLT_space.f[t], analytic, diff);
            error(diff > 1e-3, 1, "integration f not working", "t=%d diff=%.12g", t, diff);
        }

        {
            double prec = 50 * 3.33;
            HLT_type HLT_space(tmax, T, E0, HLT_EXP_b, prec);
            printf("testing integration with arb\n");
            std::vector<double>  gauss_p = { 0,1 };
            printf("gaussian mu=%.12g  sigma=%.12g\n", gauss_p[0], gauss_p[1]);
            wrapper_smearing Delta(gaussian_for_HLT, gauss_p, &HLT_space);
            HLT_space.compute_f_EXP_b(Delta, 1e-10);
            for (int t = 0;t < HLT_space.Tmax;t++) {
                double a = t + 1;
                double b = (T - (t + 1));
                printf("f[%d]=", t); arb_printn(arb_mat_entry(HLT_space.f,t,0), prec / 3.33, 0); flint_printf("\n");
                // printf("f[%d]=%.12g  analytic=%.12g  diff=%.12g\n", t, HLT_space.f[t], analytic, diff);
                // error(diff > 1e-3, 1, "integration f not working", "t=%d diff=%.12g", t, diff);
            }
        }

    }
    {

        double prec = 50 * 3.33;
        int n = 2;
        arb_mat_t A;
        arb_mat_init(A, n, n);
        arb_mat_t R;
        arb_mat_init(R, n, 1);
        arb_t det;
        arb_init(det);
        arb_set_d(arb_mat_entry(A, 0, 0), 1);
        arb_set_d(arb_mat_entry(A, 1, 1), 2);
        arb_set_d(arb_mat_entry(A, 0, 1), 2);
        arb_set_d(arb_mat_entry(A, 1, 0), 2);
        arb_mat_det(det, A, prec);
        arb_printn(det, prec / 3.33, 0); flint_printf("\n");

        arb_set_d(arb_mat_entry(R, 0, 0), 1);
        arb_set_d(arb_mat_entry(R, 1, 0), 2);
        arb_mat_mul(R,A,R,prec);
        arb_printn(arb_mat_entry(R, 0, 0), prec / 3.33, 0); flint_printf("\n");

    }
}

