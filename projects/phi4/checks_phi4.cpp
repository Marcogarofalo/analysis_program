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
#include "mass_phi4.hpp"
#include "header_phi4.hpp"
#include "correlators_analysis.hpp"
#include "eigensystem.hpp"
#include "lhs_functions.hpp"

#include <cstring> 
#include <string>
#include <fstream>
#include <memory>
using namespace std;

struct  kinematic kinematic_2pt;

int r_value(int r)
{
    int vr;
    if (r==0) vr= 1;
    else if (r==1) vr= -1;
    else error(0==0,0,"r_value","r value is not 0 neither 1\n");
    return vr;
}



double C3(int n, int Nvar, double *x,int Npar,double  *P){
    
    double C3;
    
    double E3=P[0];
    double A3=P[1];
    double A12=P[2];
    double t= x[0];
    double M=x[1];
    double E2=x[2];
    //check if file_head.l0 arrives here
    double T=(double)file_head.l0;
    C3=A3*A3 *  exp(-E3*T/2.) * cosh( E3 *(t -T/2.));
    C3+=A12*A12 *  exp(-(E3+M)*T/2.) * cosh( (E2-M) *(t -T/2.));
    
    return C3;
    
}

double C2_diff_masses(int n, int Nvar, double *x,int Npar,double  *P){
    
    double C2;
    
    double E2=P[0];
    double A2=P[1];
    double A12=P[2];
    
    double t= x[0];
    double M0=x[1];
    double M1=x[2];
    //check if file_head.l0 arrives here
    double T=(double)file_head.l0;
    C2=A2*A2 * ( exp(-t*E2)+ exp(-(T-t)*E2) );
    C2+=A12*A12 *( exp(-T*M0-t*M1 +t*M0)+ exp(-T*M1-t*M0 +t*M1) ) ;
    
    return C2;
    
}



double four_pt_BH_t_T_8(int n, int Nvar, double *x,int Npar,double  *P){
    
    double C4;
    
    double aN=P[0];
    //double norm=P[1];
    
    
    double T=(double)file_head.l0;
    double t= x[0]- T/8.;
    double M0=x[1];
    double M1=x[2];
    
    C4 = 8. *pi_greco*(M0+M1)*aN*t;
    C4 -= 16*aN*aN * sqrt( 2.*pi_greco *(M0+M1)* M0*M1*t );
    //C4 *= norm*exp(M1*t) ;
    C4 /= 8*M0*M1*t ;
    
    return C4;
    
}
double four_pt_BH_t_t3(int n, int Nvar, double *x,int Npar,double  *P){
    double C4;
    double aN=P[0];
    double T=(double)file_head.l0;
    double t= x[0]- 3.;
    double M0=x[1];
    double M1=x[2];
    
    C4 = 8. *pi_greco*(M0+M1)*aN*t;
    C4 -= 16*aN*aN * sqrt( 2.*pi_greco *(M0+M1)* M0*M1*t );
    C4 /= 8*M0*M1*t ;
    
    return C4;
}

double four_pt_BH_t_t4(int n, int Nvar, double *x,int Npar,double  *P){
    double C4;
    double aN=P[0];
    double T=(double)file_head.l0;
    double t= x[0]- 4.;
    double M0=x[1];
    double M1=x[2];
    
    C4 = 8. *pi_greco*(M0+M1)*aN*t;
    C4 -= 16*aN*aN * sqrt( 2.*pi_greco *(M0+M1)* M0*M1*t );
    C4 /= 8*M0*M1 *t;
    
    return C4;
}
double four_pt_BH_t_t5(int n, int Nvar, double *x,int Npar,double  *P){
    double C4;
    double aN=P[0];
    double T=(double)file_head.l0;
    double t= x[0]- 5.;
    double M0=x[1];
    double M1=x[2];
    
    C4 = 8. *pi_greco*(M0+M1)*aN*t;
    C4 -= 16*aN*aN * sqrt( 2.*pi_greco *(M0+M1)* M0*M1*t );
    C4 /= 8*M0*M1*t ;
    
    return C4;
}
double four_pt_BH_00_t_t3(int n, int Nvar, double *x,int Npar,double  *P){
    double C4;
    double aN=P[0];
    double T=(double)file_head.l0;
    double t= x[0]- 3.;
    double M0=x[1];
    double M1=x[2];
    
    C4 = 16. *pi_greco*(M0+M1)*aN*t;
    C4 -= 16*aN*aN * sqrt( 2.*pi_greco *(M0+M1)* M0*M1*t );
    C4 /= 8*M0*M1*t ;
    
    return C4;
}

double four_pt_BH_00_t_t4(int n, int Nvar, double *x,int Npar,double  *P){
    double C4;
    double aN=P[0];
    double T=(double)file_head.l0;
    double t= x[0]- 4.;
    double M0=x[1];
    double M1=x[2];
    
    C4 = 16. *pi_greco*(M0+M1)*aN*t;
    C4 -= 16*aN*aN * sqrt( 2.*pi_greco *(M0+M1)* M0*M1*t );
    C4 /= 8*M0*M1 *t;
    
    return C4;
}
double four_pt_BH_00_t_t5(int n, int Nvar, double *x,int Npar,double  *P){
    double C4;
    double aN=P[0];
    double T=(double)file_head.l0;
    double t= x[0]- 5.;
    double M0=x[1];
    double M1=x[2];
    
    C4 = 16. *pi_greco*(M0+M1)*aN*t;
    C4 -= 16*aN*aN * sqrt( 2.*pi_greco *(M0+M1)* M0*M1*t );
    C4 /= 8*M0*M1*t ;
    
    return C4;
}
double four_pt_BH_t_T_8_const(int n, int Nvar, double *x,int Npar,double  *P){
    double r=four_pt_BH_t_T_8( n,  Nvar, x, Npar, P);
    return r+P[1]/(x[0]- file_head.l0/8.); 
}

double four_pt_BH_t_t3_const(int n, int Nvar, double *x,int Npar,double  *P){
    double r=four_pt_BH_t_t3( n,  Nvar, x, Npar, P);
    return r+P[1]/(x[0]- 3.); 
}

double four_pt_BH_t_t4_const(int n, int Nvar, double *x,int Npar,double  *P){
    double r=four_pt_BH_t_t4( n,  Nvar, x, Npar, P);
    return r+P[1]/(x[0]- 4.); 
}
double four_pt_BH_t_t5_const(int n, int Nvar, double *x,int Npar,double  *P){
    double r=four_pt_BH_t_t5( n,  Nvar, x, Npar, P);
    return r+P[1]/(x[0]- 5.); 
}

double four_pt_BH_00_t_t3_const(int n, int Nvar, double *x,int Npar,double  *P){
    double r=four_pt_BH_00_t_t3( n,  Nvar, x, Npar, P);
    return r+P[1]/(x[0]- 3.); 
}
double four_pt_BH_00_t_t4_const(int n, int Nvar, double *x,int Npar,double  *P){
    double r=four_pt_BH_00_t_t4( n,  Nvar, x, Npar, P);
    return r+P[1]/(x[0]- 4.); 
}
double four_pt_BH_00_t_t5_const(int n, int Nvar, double *x,int Npar,double  *P){
    double r=four_pt_BH_00_t_t5( n,  Nvar, x, Npar, P);
    return r+P[1]/(x[0]- 5.); 
}


template<int tx,int delta  >
double four_pt_BH_t_tx_shifted(int n, int Nvar, double *x,int Npar,double  *P){
    double C4;
    double aN=P[0];
    double T=(double)file_head.l0;
    double t= x[0]- tx;
    double M0=x[1];
    double M1=x[2];
    
    C4 = 8. *pi_greco*(M0+M1)*aN;
    C4 -= 16*aN*aN * sqrt( 2.*pi_greco *(M0+M1)* M0*M1 )*(sqrt(t+delta)-sqrt(t))/delta;
    C4 /= 8*M0*M1;
    
    return C4;
}

template<int tx,int delta  >
double four_pt_BH_00_t_tx_shifted(int n, int Nvar, double *x,int Npar,double  *P){
    double C4;
    double aN=P[0];
    double T=(double)file_head.l0;
    double t= x[0]- tx;
    double M0=x[1];
    double M1=x[2];
    
    C4 = 16. *pi_greco*(M0+M1)*aN;
    C4 -= 16*aN*aN * sqrt( 2.*pi_greco *(M0+M1)* M0*M1 )*(sqrt(t+delta)-sqrt(t))/delta;
    C4 /= 8*M0*M1;
    
    return C4;
}


double four_pt_BH_line(int n, int Nvar, double *x,int Npar,double  *P){
    
    double C4;
    
    double aN=P[0];
    //double norm=P[1];
    
    
    double T=(double)file_head.l0;
    double t= x[0]- T/8.;
    double M0=x[1];
    double M1=x[2];
    
    C4 = 8. *pi_greco*(M0+M1)*aN*t;
    C4 /= 8*M0*M1*t ;
    return C4;
    
}

double four_pt_BH_2par(int n, int Nvar, double *x,int Npar,double  *P){
    
    double C4;
    
    double A=P[0];
    double B=P[1];
    
    
    double T=(double)file_head.l0;
    double t= x[0]- T/8.;
    double M0=x[1];
    double M1=x[2];
    
    C4 = A*8. *pi_greco*(M0+M1)*t;
    C4 += B * sqrt( 2.*pi_greco *(M0+M1)* M0*M1*t );
    //C4 *= norm*exp(M1*t) ;
    //C4 *= norm ;
    C4 /= 8*M0*M1*t ;
    return C4;
    
}


double four_pt_BH_3par(int n, int Nvar, double *x,int Npar,double  *P){
    
    double C4;
    
    double A=P[0];
    double B=P[1];
    double C=P[2];
    
    
    double T=(double)file_head.l0;
    double t= x[0]- T/8.;
    double M0=x[1];
    double M1=x[2];
    
    C4 = A*t;
    C4 += B * sqrt( 2.*pi_greco *(M0+M1)* M0*M1*t );
    C4+=C;
    //C4 *= norm*exp(M1*t) ;
    //C4 *= norm ;
    C4 /= 8*M0*M1*t ;
    return C4;
    
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


static int ****mass_index;

void init_mass_index()
{
    int k1, k2,r1,r2,i;
    int nk=file_head.nk;
    
    mass_index=(int****) malloc(sizeof(int***)*nk);
    for (k1=0;k1<nk;k1++){
        mass_index[k1]=(int***) malloc(sizeof(int**)*2);
        for (r1=0;r1<2;r1++){
            mass_index[k1][r1]=(int**) malloc(sizeof(int*)*(k1+1));
            for (k2=0;k2<=k1;k2++){
                mass_index[k1][r1][k2]=(int*) malloc(sizeof(int)*2);		
            }
        }
    }
    
    i=0;
    for (k1=0;k1<nk;k1++)
        for (r1=0;r1<2;r1++)
            for (k2=0;k2<=k1;k2++)
                for (r2=0;r2<2;r2++){
                    mass_index[k1][r1][k2][r2]=i;
                    i++;
                }
                
                
}

static void  print_file_head(FILE *stream)
{
    int i;
    
    fprintf(stream,"twist= %d\n",file_head.twist);
    fprintf(stream,"nf=%d\n",file_head.nf);
    fprintf(stream,"nsrc=%d\n",file_head.nsrc);
    fprintf(stream,"L0=%d\n",file_head.l0);
    fprintf(stream,"L1=%d\n",file_head.l1);
    fprintf(stream,"L2=%d\n",file_head.l2);
    fprintf(stream,"L3=%d\n",file_head.l3);
    fprintf(stream,"mus=%d\n",file_head.nk);
    fprintf(stream,"moms=%d\n",file_head.nmoms);
    
    fprintf(stream,"beta=%f\n",file_head.beta);
    fprintf(stream,"ksea=%f\n",file_head.ksea);
    fprintf(stream,"musea=%f\n",file_head.musea);
    fprintf(stream,"csw=%f\n",file_head.csw);
    
    fprintf(stream,"masses=");
    for(i=0;i<2*file_head.nk;++i)
        fprintf(stream,"%f\t",file_head.k[i]);
    fprintf(stream,"\n");
    
    fprintf(stream,"momenta=");
    for(i=0;i<file_head.nmoms;++i)
        fprintf(stream,"%f  %f   %f   %f\n",file_head.mom[i][0],file_head.mom[i][1],file_head.mom[i][2],file_head.mom[i][3]);
}



int read_nconfs( FILE *stream, cluster::IO_params params){
    
    long int tmp;
    int s=params.data.header_size;
    
    
    
    fseek(stream, 0, SEEK_END);
    tmp = ftell(stream);
    tmp-= params.data.header_size ;
    
    s=params.data.size;
    std::cout<< "size="<<s<<std::endl;
    
    int c= (tmp)/ ( sizeof(int)+(s)*sizeof(double) );
    
    
    std::cout<< "confs="<<c<<std::endl;
    fseek(stream, params.data.header_size, SEEK_SET);
    
    return c;
    
}
static void  write_file_head(FILE *stream)
{
    int i,dsize;
    double *dstd;
    
    fwrite(&file_head.twist,sizeof(int),1,stream);
    fwrite(&file_head.nf,sizeof(int),1,stream);
    fwrite(&file_head.nsrc,sizeof(int),1,stream);
    fwrite(&file_head.l0,sizeof(int),1,stream);
    fwrite(&file_head.l1,sizeof(int),1,stream);
    fwrite(&file_head.l2,sizeof(int),1,stream);
    fwrite(&file_head.l3,sizeof(int),1,stream);
    fwrite(&file_head.nk,sizeof(int),1,stream);
    fwrite(&file_head.nmoms,sizeof(int),1,stream);
    
    fwrite(&file_head.beta,sizeof(double),1,stream);
    fwrite(&file_head.ksea,sizeof(double),1,stream);
    fwrite(&file_head.musea,sizeof(double),1,stream);
    fwrite(&file_head.csw,sizeof(double),1,stream);
    
    fwrite(file_head.k,sizeof(double),2*file_head.nk,stream);
    
    for(i=0;i<file_head.nmoms;i++)  
        fwrite(file_head.mom[i],sizeof(double),4,stream);
}



double matrix_element_GEVP(int t, double **cor,double mass){
    double me;
    
    me=cor[t][0]/sqrt( exp(-mass*t)+exp(-(file_head.l0-t)*mass) );
    me*=2*mass;
    
    return  me;
}




void read_twopt(FILE *stream, int iconf , double ***to_write ,cluster::IO_params params, int index ){
    
    int tmp=params.data.header_size;// 
    tmp+=sizeof(double)*iconf*params.data.size+sizeof(int)*(iconf+1);
    
    
    double *obs=(double*) malloc(params.data.size*sizeof(double)); 
    
    fseek(stream, tmp, SEEK_SET);
    fread(obs,sizeof(double),params.data.size,stream); 
    
    for(int t=0 ;t<params.data.L[0];t++){
        size_t  id=index+ t*params.data.ncorr;
        (*to_write)[t][0]=obs[id];
        
    }
    free(obs);
    
    
}



void setup_single_file_jack(char  *save_name,char **argv, const char  *name,int Njack){
    FILE *f;
    mysprintf(save_name,NAMESIZE,"/dev/null");
    f=fopen(save_name,"w+");
    error(f==NULL,1,"setup_file_jack ",
          "Unable to open output file /dev/null");
    write_file_head(f);
    fwrite(&Njack,sizeof(int),1,f);
    fclose(f);
}

void setup_single_file_jack_ASCI(char  *save_name, char **argv,const char  *name,int Njack){
    FILE *f;
    mysprintf(save_name,NAMESIZE,"/dev/null");
    f=fopen(save_name,"w+");
    error(f==NULL,1,"setup_file_jack ",
          "Unable to open output file /dev/null");
    fclose(f);
}

void setup_file_jack(char **argv,int Njack){
    if( strcmp(argv[4],"jack")==0){
        setup_single_file_jack(file_jack.M_PS,argv,"/dev/null",Njack);
        setup_single_file_jack(file_jack.f_PS,argv,"/dev/null",Njack);
        setup_single_file_jack(file_jack.Zf_PS,argv,"/dev/null",Njack);
        
        setup_single_file_jack(file_jack.M_PS_GEVP,argv,"/dev/null",Njack);
        setup_single_file_jack(file_jack.f_PS_ls_ss,argv,"/dev/null",Njack);
        
    }
    
    if( strcmp(argv[4],"boot")==0){
        
        setup_single_file_jack(file_jack.M_PS,argv,"/dev/null",Njack);
        setup_single_file_jack(file_jack.f_PS,argv,"/dev/null",Njack);
        setup_single_file_jack(file_jack.Zf_PS,argv,"/dev/null",Njack);
        
        setup_single_file_jack(file_jack.M_PS_GEVP,argv,"/dev/null",Njack);
        setup_single_file_jack(file_jack.f_PS_ls_ss,argv,"/dev/null",Njack);
    }
}

int main(int argc, char **argv){
    int size;
    int i,j,t;
    
    int *iconf,confs;
    double ****data,****data_bin, **out,**tmp; 
    char c;
    clock_t t1,t2;
    double *in;
    
    double ****M,****vec,****projected_O;
    double  ****lambda,****lambda0;
    
    double *fit,***y,*x,*m,*me;
    
    
    double ****conf_jack,**r,**mt,**met;
    int Ncorr=1;
    int t0=2;
    
    FILE  *f_ll=NULL, *f_sl=NULL,*f_ls=NULL,*f_ss=NULL;
    
    FILE *plateaux_masses=NULL, *plateaux_masses_GEVP=NULL; 
    FILE *plateaux_f=NULL;   
    char namefile[NAMESIZE];
    srand(1);
    
    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux,NAMESIZE,"checks_plateaux.txt");
    
    
    error(argc!=8,1,"main ",
          "usage:./phi4  blind/see/read_plateaux -p path file -bin $bin  jack/boot \n separate path and file please");
    error(strcmp(argv[1],"blind")!=0 &&  strcmp(argv[1],"see")!=0 && strcmp(argv[1],"read_plateaux")!=0   ,1,"main ",
          "argv[1] only options:  blind/see/read_plateaux ");
    
    cluster::IO_params params;
    mysprintf(namefile,NAMESIZE,"%s/%s",argv[3],argv[4]);
    FILE *infile=open_file(namefile,"r+");
    read_header(infile,params);
    
    
    
    
    error(strcmp(argv[5],"-bin")!=0 ,1,"main", "argv[4] must be: -bin");
    error(strcmp(argv[7],"jack")!=0 &&  strcmp(argv[7],"boot")!=0,1,"main",
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
    mysprintf(option[4],NAMESIZE,argv[7]); //resampling
    mysprintf(option[5],NAMESIZE,"no"); // pdf
    mysprintf(option[6],NAMESIZE,argv[4]); // infile
    
    printf("resampling %s\n",option[4] );
    int T=params.data.L[0];
    
    double mu1=params.data.msq0;
    double mu2=params.data.msq1;
    printf("mu=%g  %g\n",mu1,mu2);
    //mysprintf(argv[4],NAMESIZE,"jack");
    
    file_head.l0=T;
    file_head.l1=params.data.L[1];file_head.l2=params.data.L[2];file_head.l3=params.data.L[3];
    file_head.nk=2;
    file_head.k= (double*) malloc(sizeof(double )*file_head.nk*2);
    file_head.k[0]=0;file_head.k[1]=0;
    file_head.k[2]=mu1;
    file_head.k[3]=mu2;
    
    file_head.nmoms=1;
    file_head.mom=(double**) malloc(sizeof(double*)*file_head.nmoms);
    for(i=0;i<file_head.nmoms;i++) {
        file_head.mom[i]=(double*) malloc(sizeof(double)*4);
        file_head.mom[i][0]=0;
        file_head.mom[i][1]=0;
        file_head.mom[i][2]=0;
        file_head.mom[i][3]=0;
    }
    
    
    
    
    
    mysprintf(namefile,NAMESIZE,"%s/out/checks_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_output",
              argv[3], T, params.data.L[1],params.data.msq0, params.data.msq1,
              params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
    printf("writing output in :\n %s \n",namefile);
    FILE *outfile=open_file(namefile,"w+");   
    
    mysprintf(namefile,NAMESIZE,"%s/out/checks_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_meff_correlators",
              argv[3], T, params.data.L[1],params.data.msq0, params.data.msq1,
              params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
    FILE *outfile_meff_corr=open_file(namefile,"w+");    
    
    mysprintf(namefile,NAMESIZE,"%s/out/checks_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_raw_correlators",
              argv[3], T, params.data.L[1],params.data.msq0, params.data.msq1,
              params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
    FILE *outfile_raw_corr=open_file(namefile,"w+"); 
    
    mysprintf(namefile,NAMESIZE,"%s/out/checks_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_shifted_correlators",
              argv[3], T, params.data.L[1],params.data.msq0, params.data.msq1,
              params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
    FILE *outfile_shifted_corr=open_file(namefile,"w+"); 
    
    mysprintf(namefile,NAMESIZE,"%s/out/checks_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_log_meff_shifted",
              argv[3], T, params.data.L[1],params.data.msq0, params.data.msq1,
              params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
    FILE *outfile_log_meff_shifted=open_file(namefile,"w+"); 
    
    mysprintf(namefile,NAMESIZE,"%s/out/checks_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_gamma",
              argv[3], T, params.data.L[1],params.data.msq0, params.data.msq1,
              params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
    
    FILE *out_gamma=open_file(namefile,"w+");      
    
    mysprintf(namefile,NAMESIZE,"%s/jackknife/%s_checks_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d",
              argv[3],option[4],
              T, params.data.L[1],params.data.msq0, params.data.msq1,
              params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
    
    FILE *jack_file=open_file(namefile,"w+");
    write_header(jack_file,params);
    
    // open infile and count the lines
    //
    
    
    int count=0;
    confs=read_nconfs( infile,  params);
    cout << "correlators ="<< params.data.ncorr<< endl;
    // compute what will be the neff after the binning 
    int bin=atoi(argv[6]);
    int Neff=confs/bin;
    cout << "effective configurations after binning (" << bin  <<"):  "<<Neff << endl;
    
    int Njack;
    if( strcmp(argv[7],"jack")==0)
        Njack=Neff+1;
    else if( strcmp(argv[7],"boot")==0)
        Njack=Nbootstrap+1;
    else
        error(1==1,1,"main","argv[7]= %s is not jack or boot",argv[7]);
    fwrite(&Njack,sizeof(int),1,jack_file );
    
    int var=params.data.ncorr;
    data=calloc_corr(confs, var,  file_head.l0 );
    
    setup_file_jack(option,Njack);
    
    get_kinematic( 0,0,  1, 0,0,  0 );
    printf("option[4]=%s\n",option[4]);
    
    for (int iconf=0; iconf< confs ;iconf++){
        
        for(int i =0 ; i< params.data.ncorr; i++)
            read_twopt(infile, iconf, &data[iconf][i], params,i);
        
        
    }
   
    symmetrise_corr(confs, 0, file_head.l0,data);
    symmetrise_corr(confs, 1, file_head.l0,data);
    
    symmetrise_corr(confs, 2, file_head.l0,data);
     
    if (params.data.ncorr>11){
        symmetrise_corr(confs, 11, file_head.l0,data);
    }
    
    data_bin=binning(confs, var, file_head.l0 ,data, bin);
    //if you want to do the gamma analysis you need to do before freeing the raw data
    effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data_bin,  Neff ,namefile_plateaux,out_gamma,0,"M_{PS}^{ll}");
    //effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data,  confs ,namefile_plateaux,out_gamma,3,"M_{PS}^{ll}");
    
    free_corr(confs, var, file_head.l0 ,data);
    
    conf_jack=create_resampling(option[4],Neff, var, file_head.l0, data_bin);
    double *zeros=(double*) calloc(Njack,sizeof(double));
    
    double **mass=(double**) malloc(sizeof(double*)*3);
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //print all the effective masses correlators
    //set the option to not read for a plateaux
    char  save_option[NAMESIZE];
    sprintf(save_option,"%s",option[1]);
    sprintf(option[1],"blind");
    
    FILE *dev_null=open_file("/dev/null","w");
    for(int icorr=0; icorr<params.data.ncorr; icorr++ ){
        //log effective mass
        double *tmp_meff_corr  =plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack
        ,namefile_plateaux,outfile_meff_corr,icorr,"meff_corr", M_eff_log,dev_null);
        free(tmp_meff_corr);
        //raw correlator
        tmp_meff_corr  =plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,
                                                      namefile_plateaux,outfile_raw_corr,icorr,"raw_corr", identity,dev_null);
        free(tmp_meff_corr);
        // shifted correlator
        tmp_meff_corr  =plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,
                                                      namefile_plateaux,outfile_shifted_corr,icorr,"shifted_corr", shift_corr,dev_null);
        free(tmp_meff_corr);
        // log_meff shifted correlator
        tmp_meff_corr  =plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack
        ,namefile_plateaux,outfile_log_meff_shifted,icorr,"log_meff_shifted", M_eff_log_shift,dev_null);
        free(tmp_meff_corr);
        
        
        
    }
    sprintf(option[1],"%s",save_option);// restore option
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // c++ 0  || r 1
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    double *E2_0_pxmpx=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,0,"E2_0_pxmpx", shift_and_M_eff_sinh_T,jack_file);
    free(E2_0_pxmpx);
    
    // c++ 1  || r 2
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    double *E2_0_pympy=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,1,"E2_0_pympy", shift_and_M_eff_sinh_T,jack_file);
    free(E2_0_pympy);
    
    // c++ 2  || r 3
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    double *E2_0_pzmpz=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,2,"E2_0_pzmpz", shift_and_M_eff_sinh_T,jack_file);
    free(E2_0_pzmpz);
    
    
    
    //free_2(3,mass);
    //free_2(3,E2);
    free_corr(Neff, var, file_head.l0 ,data_bin);
    free_jack(Njack,var , file_head.l0, conf_jack);
    
    fclose(out_gamma);
   // free(fit_info.ext_P);
    free_tower(7,(void**)option);
    
    for(i=0;i<file_head.nmoms;i++) 
        free(file_head.mom[i]);
    free(file_head.mom);
    fclose(infile); fclose(outfile);
    
    free(file_head.k);
    fclose(jack_file);
    return 0;   
}


