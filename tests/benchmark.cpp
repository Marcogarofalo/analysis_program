#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include <unistd.h>
#include <sys/time.h>
#include <fcntl.h>
#include "global.hpp"
#include "mutils.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"


double new_error_jack_biased(int Np1, double* in) {
    double r[2]={0,0};
    int i, N;

    N = Np1 - 1;

    for (i = 0;i < N;i++)
        r[0] += in[i];

    r[0] /= ((double)N);


    for (i = 0;i < N;i++) {
        r[1] += (in[i] - r[0]) * (in[i] - r[0]);
        //error(in[i]!=in[i],1,"mean_and_error_jack_biased","errore jack=%d is nan",i);
    }
    r[1] *= (N - 1.) / ((double)N);
    r[1] = sqrt(r[1]);


    return r[1];
}

int main() {
    int Njack = 20;
    int Nr = 1e+4;
    double* jack = new double[Njack];
    for (int j = 0; j < Njack; j++) {
        jack[j] = j;
    }
    {
        double a = timestamp();
        for (int r = 0; r < Nr; r++) {
            double da = error_jackboot("jack", Njack, jack);
        }
        printf("time general = %g\n", timestamp() - a);
    }
    {
        double a = timestamp();
        for (int r = 0; r < Nr; r++) {
            if (r == 99);
            else {
                double* da = mean_and_error_jack_biased(Njack, jack);
                free(da);
            }
        }
        printf("time jack = %g\n", timestamp() - a);
    }

    {
        double a = timestamp();
        myres = new resampling_jack(Njack);
        for (int r = 0; r < Nr; r++) {

            double da = myres->comp_error(jack);
        }
        printf("time class virtual func = %g\n", timestamp() - a);
    }
    {
        double a = timestamp();
        double (*p)(int, double*);
        p = new_error_jack_biased;
        for (int r = 0; r < Nr; r++) {

            double da = p(Njack, jack);

        }
        printf("time func pointer = %g\n", timestamp() - a);
    }

    delete[] jack;

}