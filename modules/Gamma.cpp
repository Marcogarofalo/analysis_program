#define GAMMA_C

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "global.hpp"
#include "tower.hpp"
#include "m_eff.hpp"

static double order = 0;
static double obs[1], dobs[1], ddobs[1], taubb_intF[1], dtau[1];
static double* abb;
static double CbbF[1], w[1];
static double** gammaFbb;

double** Gamma(int t, int var, int order, int rep, int nconf, double* a, double* bba) {
    double** r;
    int i0, i1, i2, i3, alpha, N;

    alpha = (order + 1) * var;
    N = alpha * rep;
    r = (double**)malloc(sizeof(double*) * (alpha));
    for (i0 = 0;i0 < alpha;i0++)
        r[i0] = (double*)calloc(alpha, sizeof(double));

    for (i0 = 0;i0 < alpha;i0++)
        for (i1 = 0;i1 < alpha;i1++)
            for (i2 = 0;i2 < rep;i2++)
                for (i3 = 0;i3 < (nconf - t);i3++)
                    r[i0][i1] += (a[i0 + i2 * alpha + i3 * N] - bba[i0]) * (a[i1 + i2 * alpha + (i3 + t) * N] - bba[i1]);
    // if(t==0) printf("r[0][0]=%f\n",r[0][0]);  


          /*printf("a-bba=%f\n",a[0]-bba[0]);*/
    for (i0 = 0;i0 < alpha;i0++)
        for (i1 = 0;i1 < alpha;i1++)
            r[i0][i1] /= ((double)(rep * nconf - rep * t));



    return r;
}




double** barf(int var, int order, int rep, int nconf, int flow, double* bba, double** ga, double* function(int var, int order, int flow, double* ah)) {
    double** r, * tmp, * tmp1;
    double* h, * ah;
    int i, j, N, alpha;

    N = rep * nconf;
    alpha = (order + 1) * var;
    ah = (double*)malloc(alpha * sizeof(double));
    h = (double*)malloc(alpha * sizeof(double));
    r = (double**)malloc(alpha * sizeof(double*));
    for (i = 0;i < alpha;i++) {
        r[i] = (double*)calloc(order + 1, sizeof(double));
        h[i] = sqrt(ga[i][i] / ((double)N * 4.));
        ah[i] = bba[i];
    }

    for (i = 0;i < alpha;i++) {
        ah[i] += h[i];
        tmp1 = function(var, order, flow, ah);
        ah[i] -= h[i];ah[i] -= h[i];
        tmp = function(var, order, flow, ah);
        //sub_pseries(order,tmp1,tmp,r[i]);
        for (j = 0;j <= order;j++) {
            r[i][j] = tmp1[j] - tmp[j];
            r[i][j] /= 2. * h[i];
        }
        free(tmp);free(tmp1);ah[i] += h[i];
        //scale_pseries(order,0,1./(2.*h[i]),r[i]  );
    }



    free(h);free(ah);
    return r;
}

// *a is an array of data a[ a+ rep *alpha + conf* Nv ] 
// a =0 ... (alpha-1)
//alpha =(order+1)*var
void mean_value(int var, int order, int rep, int nconf, int flow, double* a, double* function(int var, int order, int flow, double* ah)) {
    int i0, i1, i2, j, alpha, imax;
    double** ab, N;
    double* Fbb, * Fb, * tmp;

    alpha = (order + 1) * var;
    N = nconf * rep;
    imax = alpha * rep;

    for (i0 = 0;i0 < alpha;i0++)
        abb[i0] = 0;
    int i, Nfit = 1;
    // printf("a=%f\n",a[1+38*2]); 
    ab = (double**)malloc(rep * sizeof(double*));
    for (i1 = 0;i1 < rep;i1++)
        ab[i1] = (double*)calloc(alpha, sizeof(double));

    Fb = (double*)calloc(order + 1, sizeof(double));

    for (i0 = 0;i0 < alpha;i0++)
        for (i1 = 0;i1 < rep;i1++)
            for (i2 = 0;i2 < nconf;i2++)
                ab[i1][i0] += a[i0 + i1 * alpha + i2 * imax];
    //printf("nconf=%d\t replicas=%d\t alpha=%d\n",nconf,rep,alpha);

    for (i0 = 0;i0 < alpha;i0++)
        for (i1 = 0;i1 < rep;i1++)
            abb[i0] += ab[i1][i0];

    for (i0 = 0;i0 < alpha;i0++)
        for (i1 = 0;i1 < rep;i1++)
            ab[i1][i0] /= (double)nconf;

    for (i0 = 0;i0 < alpha;i0++)
        abb[i0] /= (double)(((double)nconf) * rep);


    Fbb = function(var, order, flow, abb);
    /* printf("Fbb=%f\n",Fbb[0]);*/
    if (rep == 1)   for (i0 = 0;i0 <= order;i0++) obs[i0] = Fbb[i0];
    else {
        for (i1 = 0;i1 < rep;i1++) {
            printf(" ok untill here  replica %d\n", i1);
            tmp = function(var, order, flow, ab[i1]);
            for (j = 0;j <= order;j++) {
                tmp[j] /= (double)nconf;
                Fb[j] += tmp[j];
            }
            //scale_pseries(order,0,nconf,tmp);
            //add_pseries(order,tmp,Fb,Fb);
        }
        for (j = 0;j <= order;j++)
            Fb[j] /= (double)N;
        //scale_pseries(order,0,1./((double)N),Fb);

        for (i0 = 0;i0 <= order;i0++)
            obs[i0] = (((double)rep) * Fbb[i0] - Fb[i0]) / (((double)rep) - 1.);
    }

    //free_dpseries(rep-1,ab);
    for (j = 0;j < rep;j++) {
        free(ab[j]);
    }
    free(ab);
    free(Fb);free(Fbb);
}

double* Gammaf(int var, int order, double** ga, double** fa) {
    int i, j, k, alpha;
    double* r;
    r = (double*)calloc(order + 1, sizeof(double));
    alpha = var * (order + 1);

    for (i = 0;i <= order;i++)
        for (j = 0;j < alpha;j++)
            for (k = 0;k < alpha;k++)
                r[i] += fa[j][i] * fa[k][i] * ga[j][k];

    return r;
}

void windowing(int var, int order, int rep, int nconf, int flow, double* a, double* bba, double* function(int var, int order, int flow, double* ah)) {
    double** fbba, ** tmp, * g, * tau, Caa = 0;
    int count = 0, i, j, i1, N, alpha;
    double S = 1.5;

    FILE* file_tau = NULL;
    file_tau = fopen("tau_int", "w+");
    if (file_tau == NULL) {
        printf("unable to open analysis file"); exit(0);
    }


    alpha = (order + 1) * var;

    g = (double*)calloc(order + 1, sizeof(double));
    tau = (double*)calloc(order + 1, sizeof(double));

    N = rep * nconf;

    for (i = 0;i <= order;i++) {
        CbbF[i] = 0;
        w[i] = -1;
    }

    tmp = Gamma(0, var, order, rep, nconf, a, bba);
    fbba = barf(var, order, rep, nconf, flow, bba, tmp, function);
    gammaFbb[0] = Gammaf(var, order, tmp, fbba);

    for (i1 = 0;i1 <= order;i1++)
        CbbF[i1] += gammaFbb[0][i1];
    Caa += tmp[0][0];

    for (i1 = 0;i1 < alpha;i1++)
        free(tmp[i1]);
    free(tmp);
    fprintf(file_tau, "%d   \t %d %g\n", flow, 0, 0.5);

    for (i = 1;i < nconf;i++) {

        tmp = Gamma(i, var, order, rep, nconf, a, bba);
        gammaFbb[i] = Gammaf(var, order, tmp, fbba);
        if (i == 1) for (j = 0;j <= order;j++) {
            if (gammaFbb[1][j] < 0) {
                /*printf("there is no autocorrelation\n");*/
                w[j] = 0;count++;
            }
        }
        if (w[0] == -1)  Caa += 2 * tmp[0][0];
        //free_dpseries(alpha-1,tmp);
        for (i1 = 0;i1 < alpha;i1++)
            free(tmp[i1]);
        free(tmp);
        for (j = 0;j <= order;j++) {

            if (w[j] == -1) {
                CbbF[j] += 2. * gammaFbb[i][j];
                taubb_intF[j] = CbbF[j] / (2. * gammaFbb[0][j]);
                fprintf(file_tau, "%d   \t %d %g\n", flow, i, taubb_intF[j]);
                tau[j] = 0.6;
                if (taubb_intF[j] > 0.5)
                    tau[j] = S / (log((2. * taubb_intF[j] + 1.) / (2. * taubb_intF[j] - 1.)));
                g[j] = exp(-((double)i) / tau[j]) - (tau[j] / (sqrt((double)(i * N))));
                /*printf("g=%f",g[j]);*/
                if (g[j] < 0) {
                    count++;  w[j] = i;
                }

            }
            /*if(j==0) printf("gammaFbb[%d]=%0.10f\n",i,gammaFbb[i][0]);*/
            if (count == order + 1) break;
        }
        free(gammaFbb[i]);
        if (count == order + 1) break;
    }



    //free_dpseries(alpha-1,fbba);
    for (i1 = 0;i1 < alpha;i1++)
        free(fbba[i1]);
    free(fbba);
    free(g);free(tau);

    for (j = 0;j <= order;j++) {
        if (w[j] == -1) {
            printf("Windowing condition order %d failed up to W = %d\n", j, nconf - 1);
            w[j] = nconf - 1;
        }
    }


    for (j = 0;j <= order;j++) {

        gammaFbb[0][j] += CbbF[j] / ((double)N);
        CbbF[j] += CbbF[j] * (2. * w[j] + 1) / ((double)N);
        taubb_intF[j] = CbbF[j] / (2. * gammaFbb[0][j]);
    }
    free(gammaFbb[0]);
    fprintf(file_tau, "\n");
    fclose(file_tau);

}


void return_answer(int var, int order, int rep, int nconf) {
    int i, N;


    N = rep * nconf;
    for (i = 0;i <= order;i++) {
        dobs[i] = CbbF[i] / ((double)N);
        dobs[i] = sqrt(dobs[i]);

        ddobs[i] = dobs[i] * sqrt((w[i] + 0.5) / N);
        dtau[i] = sqrt((w[i] + 0.5 - taubb_intF[i]) / ((double)N)) * 2. * taubb_intF[i];
    }

}

//replicas =rep 
// flow is an extra parameter to pass to the function
double* analysis_gamma(int var, int rep, int nconf, int flow, double* a, double* function(int var, int order, int flow, double* ah)) {
    double* r = (double*)malloc(sizeof(double) * 5);

    abb = (double*)calloc(var, sizeof(double));
    gammaFbb = (double**)malloc(sizeof(double*) * nconf);

    // inside they use gammaFbb
    mean_value(var, order, rep, nconf, flow, a, function);
    windowing(var, order, rep, nconf, flow, a, abb, function);
    return_answer(var, order, rep, nconf);//printf("HERE  %d %d %f\n",k,i+L0,abb[2*39]);
    //fprintf(file_out,"%d \t %.15f  %.15f %.15f %.15f %.15f \n",i,obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);


    r[0] = obs[0];
    r[1] = dobs[0];
    r[2] = ddobs[0];
    r[3] = taubb_intF[0];
    r[4] = dtau[0];

    free(abb);
    free(gammaFbb);
    return r;
}


double* identity_func_gamma(int var, int order, int flow, double* ah) {
    double* r = (double*)calloc((1), sizeof(double));

    r[0] = ah[0];
    return r;
}
//double   *effective_mass_phi4_gamma(char **option ,struct kinematic kinematic_2pt , char* name, double ****data, int Confs ,FILE **plateaux_masses,FILE *outfile,  int index , const char *description ){
void gamma_correlator(char** option, struct kinematic kinematic_2pt,
    char* name, double**** data, int Confs, const char* plateaux_masses, FILE* outfile,
    int index, const char* description) {

    //int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   //if ( strcmp(option[1],"read_plateaux")==0 )
   //	go_to_line(*plateaux_masses,line);

    double** r, * m, ** mt, * fit;
    int i, j, yn;

    r = (double**)malloc(sizeof(double*) * file_head.l0);
    for (i = 0;i < file_head.l0;i++)
        r[i] = (double*)malloc(sizeof(double) * Confs);
    //mt=(double**) malloc(sizeof(double*)*file_head.l0);


    fprintf(outfile, "#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n", name,
        kinematic_2pt.k2, kinematic_2pt.r2, kinematic_2pt.mom2,
        kinematic_2pt.k1, kinematic_2pt.r1, kinematic_2pt.mom1);
    double* datag;
    int Nvar=1;
    datag = (double*)malloc(sizeof(double) * Confs * Nvar);
    for (i = 1;i < file_head.l0 / 2;i++) {

        for (j = 0;j < Confs;j++) {
            //store in the format for the gamma analysis
            //two variables c(t) and c(t+1), order=0 
            datag[j *Nvar+0] = data[j][index][i][0];
            // datag[1 + j * 2] = data[j][index][i + 1][0];
        }

        double* obs = analysis_gamma(1, 1, Confs, i//time
            , datag, identity_func_gamma);

        fprintf(outfile, "%d   %.15e    %.15e   %.15e  %.15e  %.15e\n", i, obs[0], obs[1], obs[2], obs[3], obs[4]);
        free(obs);

    }

    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, 1, 1, 0.0, 0.0);
    // fprintf(outfile, "%.15g    %.15g    \n", 0.0, 0.0);

    free(datag);
    for (i = 0;i < file_head.l0;i++)
        free(r[i]);
    free(r);

    fflush(outfile);



}
double *mass_gamma(int var, int order,int flow ,double *ah){
    double *r=(double*) calloc((1),sizeof(double)); 
    
    //use flow at time of the correlator
    // we need to pass to M_eff a correlator so we create a double **c with has only the correlator at time=flow and flow+1
    double **c=double_malloc_2(flow+2,2);
    c[flow][0]=ah[0];
    c[flow+1][0]=ah[1];
    r[0]=M_eff(flow,c);
    free_2(flow+2,c);
    return r;
}
//double   *effective_mass_phi4_gamma(char **option ,struct kinematic kinematic_2pt , char* name, double ****data, int Confs ,FILE **plateaux_masses,FILE *outfile,  int index , const char *description ){
void gamma_correlator_func(char** option, struct kinematic kinematic_2pt,
    char* name, double**** data, int Confs, const char* plateaux_masses, FILE* outfile,
    int index, const char* description, double* function(int var, int order, int flow, double* ah)) {
    //int line=kinematic_2pt.ik2+kinematic_2pt.ik1*(file_head.nk+1);
   //if ( strcmp(option[1],"read_plateaux")==0 )
   //	go_to_line(*plateaux_masses,line);
   
   double **r,*m,**mt,*fit;
   int i,j,yn;
    
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Confs);
   //mt=(double**) malloc(sizeof(double*)*file_head.l0);

   
    

   fprintf(outfile,"\n\n#m_eff(t) of %s  propagators:1) mu %.5f r %d theta %.5f 2) mu %.5f r %d theta %.5f\n",name,
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
                                     , datag  , function);
          // mt[i][0]=obs[0];
          // mt[i][0]=obs[1];
           fprintf(outfile,"%d   %.15e    %.15e   %.15e  %.15e  %.15e\n",i,obs[0],obs[1],obs[2],obs[3],obs[4]);
           free(obs);
           free(datag);

   }

   //fit=fit_plateaux(option, kinematic_2pt ,  name,description/*"M_{PS}^{ll}"*/,mt,r,  Confs,*plateaux_masses,outfile);
  // write_jack_bin(Confs,fit,file_jack.M_PS);
     

  fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, 1, 1, 0.0, 0.0);
  fprintf(outfile, "%.15g    %.15g    \n", 0.0, 0.0);
  /* for(i=1;i<file_head.l0/2;i++)
       free(mt[i]);
   free(mt);*/
   for(i=0;i<file_head.l0;i++)
       free(r[i]);
   free(r);

   fflush(outfile);
   
    /*if ( strcmp(option[1],"read_plateaux")==0 ){
     fclose(*plateaux_masses);
     *plateaux_masses=open_file(kinematic_2pt.plateau_m_ll,"r");

    }*/
    
    //return fit;    
    
}