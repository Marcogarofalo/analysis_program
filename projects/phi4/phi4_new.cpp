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

extern "C" { 
    #include "dzeta_function.h"
}

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
    else{ vr=-10; error(0==0,0,"r_value","r value is not 0 neither 1\n");}
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
   char namefile_plateaux[NAMESIZE];
   mysprintf(namefile_plateaux,NAMESIZE,"plateaux.txt");
   FILE *plateaux_f=NULL;   
   char namefile[NAMESIZE];
   srand(1);
   
   
   
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
   
  
   
   
   
   mysprintf(namefile,NAMESIZE,"%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_output",
             argv[3], T, params.data.L[1],params.data.msq0, params.data.msq1,
             params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
   printf("writing output in :\n %s \n",namefile);
   FILE *outfile=open_file(namefile,"w+");   

   mysprintf(namefile,NAMESIZE,"%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_meff_correlators",
             argv[3], T, params.data.L[1],params.data.msq0, params.data.msq1,
             params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
   FILE *outfile_meff_corr=open_file(namefile,"w+");    
   
   mysprintf(namefile,NAMESIZE,"%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_raw_correlators",
             argv[3], T, params.data.L[1],params.data.msq0, params.data.msq1,
             params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
   FILE *outfile_raw_corr=open_file(namefile,"w+"); 
   
   mysprintf(namefile,NAMESIZE,"%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_shifted_correlators",
             argv[3], T, params.data.L[1],params.data.msq0, params.data.msq1,
             params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
   FILE *outfile_shifted_corr=open_file(namefile,"w+"); 
   
   mysprintf(namefile,NAMESIZE,"%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_log_meff_shifted",
             argv[3], T, params.data.L[1],params.data.msq0, params.data.msq1,
             params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
   FILE *outfile_log_meff_shifted=open_file(namefile,"w+"); 
  
   mysprintf(namefile,NAMESIZE,"%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_gamma",
             argv[3], T, params.data.L[1],params.data.msq0, params.data.msq1,
             params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
   
   FILE *out_gamma=open_file(namefile,"w+");      
   
   mysprintf(namefile,NAMESIZE,"%s/jackknife/%s_G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d",
             argv[3],option[4],
             T, params.data.L[1],params.data.msq0, params.data.msq1,
             params.data.lambdaC0, params.data.lambdaC1, params.data.muC, params.data.gC,params.data.replica);
   
   FILE *jack_file=open_file(namefile,"w+");
   write_header(jack_file,params);
   
   // open infile and count the lines
   //
   
   
   int count=0;
   confs=read_nconfs( infile,  params);
   //confs=confs/2;
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
        
        
            
        //read_twopt(infile, iconf, &data[iconf][0], params,0);//2pt 0
        //read_twopt(infile, iconf, &data[iconf][1], params,1);//2pt 1
        
        //read_twopt(infile, iconf, &data[iconf][2], params,2);//C2t0
        //read_twopt(infile, iconf, &data[iconf][3], params,3);//C2t1
        //read_twopt(infile, iconf, &data[iconf][4], params,4);//C2t
        
        //read_twopt(infile, iconf, &data[iconf][5], params,5);//C3t0
        //read_twopt(infile, iconf, &data[iconf][6], params,6);//C3t1
        //read_twopt(infile, iconf, &data[iconf][7], params,7);//C3t
        
        //read_twopt(infile, iconf, &data[iconf][8], params,8);//C4t0
        //read_twopt(infile, iconf, &data[iconf][9], params,9);//C4t1
        //read_twopt(infile, iconf, &data[iconf][10], params,10);//C4t
        
        
        //read_twopt(infile, iconf, &data[iconf][11], params,11);//C201
        
    }
/*
    for(int i =0 ; i< 2; i++){
        std::string myname="one_to_one_"+ std::to_string(i) + ".txt";
        FILE *fi=open_file(myname.c_str() ,"w+" );
        for (int iconf=0; iconf< confs ;iconf++){
            for(int t=0;t<params.data.L[0];t++)
                fprintf(fi,"%d   %d  %.12g\n",iconf,t,data[iconf][i][t][0]);
        }
        fclose(fi);
    }
  */  
    symmetrise_corr(confs, 0, file_head.l0,data);
    symmetrise_corr(confs, 1, file_head.l0,data);
    
    symmetrise_corr(confs, 2, file_head.l0,data);
    symmetrise_corr(confs, 3, file_head.l0,data);
    symmetrise_corr(confs, 4, file_head.l0,data);
    
    symmetrise_corr(confs, 5, file_head.l0,data);
    symmetrise_corr(confs, 6, file_head.l0,data);
    symmetrise_corr(confs, 7, file_head.l0,data);
    
    if (params.data.ncorr>11){
        symmetrise_corr(confs, 11, file_head.l0,data);
    }
    if(params.data.ncorr>33){
        for(int i =33 ;i< 48; i++)
            symmetrise_corr(confs, i, file_head.l0,data);
    }
   
   
   FILE *f3t16=fopen("E1_1.txt","w+");
   for (int iconf=0; iconf< confs ;iconf++){
       for (int t =0; t< T ;t++){
           fprintf(f3t16,"%.12g  ",data[iconf][1][t][0]);
       }
       fprintf(f3t16,"\n");
   }
   fclose(f3t16);
   
    data_bin=binning(confs, var, file_head.l0 ,data, bin);
    //if you want to do the gamma analysis you need to do before freeing the raw data
    effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data_bin,  Neff ,namefile_plateaux,out_gamma,0,"M_{PS}^{ll}");
    //effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data,  confs ,namefile_plateaux,out_gamma,3,"M_{PS}^{ll}");
    
    free_corr(confs, var, file_head.l0 ,data);
    
    conf_jack=create_resampling(option[4],Neff, var, file_head.l0, data_bin);
    double *zeros=(double*) calloc(Njack,sizeof(double));
    
    double **mass=(double**) malloc(sizeof(double*)*3);
    fprintf(outfile,"#correlator\n");
    for (int t =1; t< T/2 ;t++){
        for (int v=0;v< 2;v++){
            double *mj0=(double*) malloc(sizeof(double)*Njack); 
            for (int j=0 ;j<Njack;j++){
                  mj0[j]=conf_jack[j][v][t][0]-conf_jack[j][v][t+1][0];
            }
            double *m=mean_and_error(option[4],Njack,mj0);
            if (v==0){        fprintf(outfile,"%d  %g  %g\t",t,m[0],m[1]);   free(m);}
            else {            fprintf(outfile,"%g  %g\t",m[0],m[1]);   free(m);}
            free(mj0);
        }
        fprintf(outfile,"\n");
    }
    
    /////////////////// some declarations///////////
    double *E1_0_px,*E1_1_px,*E1_0_py,*E1_1_py,*E1_0_pz,*E1_1_pz;
    
    
    
    /////!!!!!!!!!!!!!!!!! write 0 in the first jackknife!!!!!!!!!!!!!
    double *mj0=(double*) calloc(Njack,sizeof(double)); 
    fwrite(mj0,sizeof(double),Njack,jack_file);
    free(mj0);
    fprintf(outfile,"\n\n #%s fit in [%d,%d] chi2=%.5f\n  %.15g    %.15g    %d   %d\n\n\n","#need_for_gnuplot",0,0,0.0,0.0,0.0,0,0);
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
    
    file_head.k[2]=mu1;
    file_head.k[3]=mu1;

    //c++ 1 || r 2
    mass[0]=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,0,"E1_0", M_eff_T,jack_file);
    //mass=compute_effective_mass(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,0,"M_{PS}^{ll}");
    
    file_head.k[2]=mu2;
    file_head.k[3]=mu2;
    
    //c++ 2 || r 3
    mass[1]=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,1,"E1_1", M_eff_T,jack_file);
    
    //!!!!!
    //there is not this correlation function  <phi0 phi1>
    //!!!!//c++ 3 || r 4
    mass[2]=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,1,"E1", M_eff_T,jack_file);
    //compute_effective_mass(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,1,"M_{PS}^{ll}");

    
    //c++ 4 || r 5
    double **E2=(double**) malloc(sizeof(double*)*3);
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    E2[0]=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,2,"E2_0", shift_and_M_eff_sinh_T,jack_file);
    double *a_0=scattering_len_luscher(  Njack,  mass[0], mass[0], E2[0] ,params.data.L[3]);
    double *tmpj=(double*) malloc(sizeof(double)*Njack);
    sub_jackboot(Njack,  tmpj, E2[0], mass[0] );
    sub_jackboot(Njack,  tmpj, tmpj, mass[0] );
    fprintf(outfile,"#scattering length  a  err deltaE2 err    mu00 err    deltaE2*mu00  err   a_00*(m0)=-3lambda/4pi  err\n %g  %g     %g  %g\t",
           a_0[Njack-1], error_jackboot(option[4],Njack,a_0),   tmpj[Njack-1],  error_jackboot(option[4],Njack,tmpj));
    
    
    double *tmp_muj=(double*) malloc(sizeof(double)*Njack);
      
    //reduced mass
    for (j=0; j<Njack;j++)
        tmp_muj[j]=mass[0][j]*mass[0][j]/(mass[0][j]+mass[0][j]);
    fprintf(outfile,"%g   %g\t", tmp_muj[Njack-1], error_jackboot(option[4],Njack,tmp_muj) );
    
    //reduced mass time DeltaE2
    for (j=0; j<Njack;j++)
        tmp_muj[j]=tmp_muj[j]*tmpj[j];
     fprintf(outfile,"%g   %g\t", tmp_muj[Njack-1], error_jackboot(option[4],Njack,tmp_muj)  );
    
    double *a0m0=(double*) malloc(sizeof(double)*Njack);
    for (j=0; j<Njack;j++)
        a0m0[j]=a_0[j]*(mass[0][j]);
     fprintf(outfile,"%g   %g\n", a0m0[Njack-1], error_jackboot(option[4],Njack,a0m0)  );
    
    
    free(tmpj); free(tmp_muj);
    ///////////////////////////////////////////////////////////////////////////////////////
    
    
    //c++ 5 || r 6
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    E2[1]=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,3,"E2_1", shift_and_M_eff_sinh_T, jack_file);
    double *a_1=scattering_len_luscher(  Njack,  mass[1], mass[1], E2[1] ,params.data.L[1]);
    tmpj=(double*) malloc(sizeof(double)*Njack);
    sub_jackboot(Njack,  tmpj, E2[1], mass[1] );
    sub_jackboot(Njack,  tmpj, tmpj, mass[1] );
    fprintf(outfile,"#scattering length  a  err deltaE2 err\n %g  %g     %g  %g\n",
           a_1[Njack-1], error_jackboot(option[4],Njack,a_1),   tmpj[Njack-1],  error_jackboot(option[4],Njack,tmpj));
    free(tmpj);
    
    //c++ 6 || r 7
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    E2[2]=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,4,"E2", shift_and_M_eff_sinh_T,jack_file);

    
    struct fit_type fit_info;
    struct fit_result  fit_out;
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.function=C3;
    fit_info.n_ext_P=2;
    fit_info.ext_P=(double**) malloc(sizeof(double*)*2);
    
    //c++ 7 || r 8
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=E2[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  5,0/*reim*/ , "E3_0",  fit_info, jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 8 || r 9
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=E2[1];
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  6,0/*reim*/ , "E3_1",  fit_info ,jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 9 || r 10
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=E2[2];
    fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  7,0/*reim*/ , "E3",  fit_info 
        ,jack_file    );
    free_fit_result(fit_info,fit_out);
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_t_T_8;
    //c++ 10 || r 11
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_0 , "E4_0",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 11 || r 12
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_1 , "E4_1",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    //c++ 12 || r 13
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH , "E4",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    //c++ 13 || r 14
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH , "E4_plat1",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    //c++ 14 || r 15
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH , "E4_plat2",  fit_info, jack_file);
    free_fit_result(fit_info,fit_out);
    //c++ 15 || r 16
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH , "E4_line",  fit_info,jack_file);
    free_fit_result(fit_info,fit_out);
    

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    fit_info.Nvar=3;
    fit_info.Npar=2;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_2par;

    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    //c++ 16 || r 17
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH , "E4_2p",  fit_info, jack_file);
    free_fit_result(fit_info,fit_out);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    fit_info.Nvar=3;
    fit_info.Npar=2;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_t_T_8_const;

    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    //c++ 17 || r 18
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH , "E4_const",  fit_info, jack_file);
    free_fit_result(fit_info,fit_out);

    
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_3par;

    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
        //c++ 18 || r 19
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH , "E4_3p",  fit_info, jack_file);
    free_fit_result(fit_info,fit_out);
    
   
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=C2_diff_masses;

    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
        //c++ 19 || r 20
    fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  11,0/*reim*/ , "E2_01",  fit_info ,jack_file);
    double *a=scattering_len_luscher(  Njack,  mass[0], mass[1], fit_out.P[0] ,params.data.L[1]);
    tmpj=(double*) malloc(sizeof(double)*Njack);
    tmp_muj=(double*) malloc(sizeof(double)*Njack);
    sub_jackboot(Njack,  tmpj, fit_out.P[0], mass[0] );
    sub_jackboot(Njack,  tmpj, tmpj, mass[1] );
    fprintf(outfile,"#scattering length  a  err deltaE2 err    mu01 err    deltaE2*mu01  err   a_01*(m0+m1)=-mu/2pi  err      a_01*(m0+m1)/(a0m0)=4/3 a_01*pi/mu   err    \n %g  %g     %g  %g\t",
           a[Njack-1], error_jackboot(option[4],Njack,a),   tmpj[Njack-1],  error_jackboot(option[4],Njack,tmpj));
    
    /*
    FILE *fm0=open_file("mass_0.txt","w+");
    FILE *fm1=open_file("mass_1.txt","w+");
    FILE *fm2=open_file("mass_2.txt","w+");
    FILE *fE2_01=open_file("E2_01.txt","w+");
    for (j=0; j<Njack;j++){
        fprintf(fm0,"%.12g\n",mass[0][j]);
        fprintf(fm1,"%.12g\n",mass[1][j]);
        fprintf(fm2,"%.12g\n",mass[2][j]);
        fprintf(fE2_01,"%.12g\n",fit_out.P[0][j]);
    }
    */
    //reduced mass
    for (j=0; j<Njack;j++)
        tmp_muj[j]=mass[0][j]* mass[1][j]/(mass[0][j]+ mass[1][j]);
    fprintf(outfile,"%g   %g\t", tmp_muj[Njack-1], error_jackboot(option[4],Njack,tmp_muj) );
    
    //reduced mass time DeltaE2
    for (j=0; j<Njack;j++)
        tmp_muj[j]=tmp_muj[j]*tmpj[j];
     fprintf(outfile,"%g   %g\t", tmp_muj[Njack-1], error_jackboot(option[4],Njack,tmp_muj)  );
    
    
    for (j=0; j<Njack;j++)
        tmp_muj[j]=a[j]*(mass[0][j]+mass[1][j]);
    fprintf(outfile,"%g   %g\t", tmp_muj[Njack-1], error_jackboot(option[4],Njack,tmp_muj)  );
    
    for (j=0; j<Njack;j++)
        tmp_muj[j]=(a[j]*(mass[0][j]+mass[1][j]))/a0m0[j];
    fprintf(outfile,"%g   %g\t", tmp_muj[Njack-1], error_jackboot(option[4],Njack,tmp_muj)  );

    // a_01*pi/mu
    for (j=0; j<Njack;j++)
        tmp_muj[j]=a[j]*pi_greco*(mass[0][j]+mass[1][j])/(mass[0][j]*mass[1][j]) ;
    fprintf(outfile,"%g   %g\n", tmp_muj[Njack-1], error_jackboot(option[4],Njack,tmp_muj)  );
    
    
    free(tmpj); free(tmp_muj);
    
    free_fit_result(fit_info,fit_out);
    fflush(outfile);
    
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //GEVP two particle
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    fit_info.Nvar=1;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=0;
    fit_info.function=constant_fit;

    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
        //c++ 20 || r 21
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, GEVP_shift_matrix , "GEVP_E2_01",  fit_info, jack_file);
    free_fit_result(fit_info,fit_out);
       
    
    
    
if(params.data.ncorr>15){    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// C4_BH_03t16
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_t3;
    //c++ 21 || r 22
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_0_03t16 , "E4_0_03t16",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 22 || r 23
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_1_03t16, "E4_1_03t16",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 23 || r 24
    fit_info.function=four_pt_BH_t_t3;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_03t16 , "E4_03t16",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    ///// const fit
    
    fit_info.Nvar=3;
    fit_info.Npar=2;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_t3_const;
    
    //c++ 24 || r 25
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_0_03t16 , "E4_0_03t16_const",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 25 || r 26
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_1_03t16, "E4_1_03t16_const",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 26 || r 27
    fit_info.function=four_pt_BH_t_t3_const;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_03t16 , "E4_03t16_const",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);


    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// C4_BH_04t16
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_t4;
    //c++ 27 || r 22
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_0_04t16 , "E4_0_04t16",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 28 || r 23
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_1_04t16, "E4_1_04t16",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 29 || r 24
    fit_info.function=four_pt_BH_t_t4;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_04t16 , "E4_04t16",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    ///// const fit
    fit_info.Nvar=3;
    fit_info.Npar=2;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_t4_const;
    
    //c++ 30 || r 25
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_0_04t16 , "E4_0_04t16_const",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 31 || r 26
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_1_04t16, "E4_1_04t16_const",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 32 || r 27
    fit_info.function=four_pt_BH_t_t4_const;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_04t16 , "E4_04t16_const",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// C4_BH_03t20
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_t3;
    //c++ 33 || r 22
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_0_03t20 , "E4_0_03t20",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 34 || r 23
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_1_03t20, "E4_1_03t20",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 35 || r 24
    fit_info.function=four_pt_BH_t_t3;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_03t20, "E4_03t20",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    ///// const fit
    fit_info.Nvar=3;
    fit_info.Npar=2;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_t3_const;
    
    //c++ 36 || r 25
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_0_03t20 , "E4_0_03t20_const",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 37 || r 26
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_1_03t20, "E4_1_03t20_const",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 38 || r 27
    fit_info.function=four_pt_BH_t_t3_const;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_03t20 , "E4_03t20_const",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// C4_BH_04t20
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_t4;
    //c++ 39 || r 22
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_0_04t20 , "E4_0_04t20",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 40 || r 23
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_1_04t20, "E4_1_04t20",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 41 || r 24
    fit_info.function=four_pt_BH_t_t4;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_04t20, "E4_04t20",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    ///// const fit
    fit_info.Nvar=3;
    fit_info.Npar=2;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_t4_const;
    
    //c++ 42 || r 25
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_0_04t20 , "E4_0_04t20_const",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 43 || r 26
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_1_04t20, "E4_1_04t20_const",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 44 || r 27
    fit_info.function=four_pt_BH_t_t4_const;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_04t20 , "E4_04t20_const",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    
     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// C4_BH_05t20
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_t5;
    //c++ 45 || r 22
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_0_05t20 , "E4_0_05t20",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 46 || r 23
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_1_05t20, "E4_1_05t20",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 47 || r 24
    fit_info.function=four_pt_BH_t_t5;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_05t20, "E4_05t20",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    ///// const fit
    fit_info.Nvar=3;
    fit_info.Npar=2;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_t5_const;
    
    //c++ 48 || r 25
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_0_05t20 , "E4_0_05t20_const",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 49 || r 26
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_1_05t20, "E4_1_05t20_const",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 50 || r 27
    fit_info.function=four_pt_BH_t_t5_const;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_05t20 , "E4_05t20_const",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    
    
    ///////////////////////togrep:shifted fit of BH
    //03t16
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_tx_shifted<3,1>;
    
    //c++ 51 || r 52
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   lhs_four_BH_00_tx_tf_shifetd<1, 3,16,0,0, 15> , "E4_0_03t16_shifted",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 52 || r 53
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   lhs_four_BH_01_tx_tf_shifetd<1, 3,16,1,1, 16>, "E4_1_03t16_shifted",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    fit_info.function=four_pt_BH_t_tx_shifted<3,1>;
    //c++ 53 || r 54
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   lhs_four_BH_01_tx_tf_shifetd<1, 3,16,0,1, 17> , "E4_03t16_shifted",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    //////////////////////
    //04t16
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_tx_shifted<4,1>;
    
    //c++ 54 || r 55
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   lhs_four_BH_00_tx_tf_shifetd<1, 4,16,0,0, 18> , "E4_0_04t16_shifted",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 55 || r 56
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,
                                   lhs_four_BH_01_tx_tf_shifetd<1, 4,16,1,1, 19>, "E4_1_04t16_shifted",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 56 || r 57
    fit_info.function=four_pt_BH_t_tx_shifted<4,1>;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   lhs_four_BH_01_tx_tf_shifetd<1, 4,16,0,1, 20> , "E4_04t16_shifted",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //////////////////////
    //03t20
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_tx_shifted<3,1>;
    
    //c++ 57 || r 58
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   lhs_four_BH_00_tx_tf_shifetd<1, 3,20,0,0, 21> , "E4_0_03t20_shifted",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 58 || r 59
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,
                                   lhs_four_BH_01_tx_tf_shifetd<1, 3,20,1,1, 22>, "E4_1_03t20_shifted",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 59 || r 60
    fit_info.function=four_pt_BH_t_tx_shifted<3,1>;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,
                                   lhs_four_BH_01_tx_tf_shifetd<1, 3,20,0,1, 23> , "E4_03t20_shifted",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    
    //////////////////////
    //04t20
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_tx_shifted<4,1>;
    
    //c++ 60 || r 61
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,
                                   lhs_four_BH_00_tx_tf_shifetd<1, 4,20,0,0, 24> , "E4_0_04t20_shifted",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 61 || r 62
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,
                                   lhs_four_BH_01_tx_tf_shifetd<1, 4,20,1,1, 25>, "E4_1_04t20_shifted",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 62 || r 63
    fit_info.function=four_pt_BH_t_tx_shifted<4,1>;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,
                                   lhs_four_BH_01_tx_tf_shifetd<1, 4,20,0,1, 26> , "E4_04t20_shifted",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //////////////////////
    //05t20
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_00_t_tx_shifted<5,1>;
    
    //c++ 63 || r 64
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,
                                   lhs_four_BH_00_tx_tf_shifetd<1, 5,20,0,0, 27> , "E4_0_05t20_shifted",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    //c++ 64 || r 65
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,
                                   lhs_four_BH_01_tx_tf_shifetd<1, 5,20,1,1, 28>, "E4_1_05t20_shifted",  fit_info,  jack_file);
    free_fit_result(fit_info,fit_out);
    
    //c++ 65 || r 66
    fit_info.function=four_pt_BH_t_tx_shifted<5,1>;
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   lhs_four_BH_01_tx_tf_shifetd<1, 5,20,0,1, 29> , "E4_05t20_shifted",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    
    
    /////HERE we try different shifts
    //03t16
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_t_tx_shifted<3,2>;
    
    //c++ 66 || r 67
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   lhs_four_BH_01_tx_tf_shifetd<2, 3,16,0,1, 17> , "E4_03t16_shifted_2",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    fit_info.function=four_pt_BH_t_tx_shifted<3,3>;
    
    //c++ 67 || r 68
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,
                                   lhs_four_BH_01_tx_tf_shifetd<3, 3,16,0,1, 17> , "E4_03t16_shifted_3",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    fit_info.function=four_pt_BH_t_tx_shifted<3,4>;
    
    //c++ 68 || r 69
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,
                                   lhs_four_BH_01_tx_tf_shifetd<4, 3,16,0,1, 17> , "E4_03t16_shifted_4",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    fit_info.function=four_pt_BH_t_tx_shifted<3,5>;
    
    //c++ 69 || r 70
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,
                                   lhs_four_BH_01_tx_tf_shifetd<4, 3,16,0,1, 17> , "E4_03t16_shifted_5",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    
if(params.data.ncorr>30){    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// C4_BH_10_03t16
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_t_t3;
     
    //c++ 70 || r 71
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[0];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_10_03t16 , "E4_10_03t16",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    ///// const fit
    
    fit_info.Nvar=3;
    fit_info.Npar=2;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_t_t3_const;
    
     
    //c++ 71 || r 72
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[1];
    fit_info.ext_P[1]=mass[0];
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, lhs_four_BH_10_03t16 , "E4_10_03t16_const",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
 
if(params.data.ncorr>32){
        
    //02t10
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_t_tx_shifted<2,1>;
    
    //c++ 72 || r 73
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   lhs_four_BH_01_tx_tf_shifetd<1, 2,10,0,1, 31> , "E4_02t10_shifted",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
    
    
    //02t12
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=four_pt_BH_t_tx_shifted<2,1>;
    
    //c++ 73 || r 74
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[0];
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   lhs_four_BH_01_tx_tf_shifetd<1, 2,12,0,1, 32> , "E4_02t12_shifted",  fit_info, jack_file );
    free_fit_result(fit_info,fit_out);
  
if(params.data.ncorr>33){  
    file_head.k[2]=mu1;
    file_head.k[3]=mu1;
    //c++ 74 || r 75
    E1_0_px=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,33,"E1_0_px", M_eff_T,jack_file);
    //free(E1_0_px);
    
    file_head.k[2]=mu2;
    file_head.k[3]=mu2;
    //c++ 75 || r 76
    E1_1_px=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,34,"E1_1_px", M_eff_T,jack_file);
    //free(E1_1_px);
    
    file_head.k[2]=mu1;
    file_head.k[3]=mu1;
    //c++ 76 || r 77
    E1_0_py=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,35,"E1_0_py", M_eff_T,jack_file);
    //free(E1_0_py);
    
    file_head.k[2]=mu2;
    file_head.k[3]=mu2;
    //c++ 77 || r 78
    E1_1_py=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,36,"E1_1_py", M_eff_T,jack_file);
    //free(E1_1_py);
    
    
    
    file_head.k[2]=mu1;
    file_head.k[3]=mu1;
    //c++ 78 || r 79
    E1_0_pz=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,37,"E1_0_pz", M_eff_T,jack_file);
    //free(E1_0_pz);
    
    file_head.k[2]=mu2;
    file_head.k[3]=mu2;
    //c++ 79 || r 80
    E1_1_pz=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,38,"E1_1_pz", M_eff_T,jack_file);
    //free(E1_1_pz);
    
    
    //c++ 80 || r 81
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    double *E2_0_A1=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,39,"E2_0_A1", shift_and_M_eff_sinh_T,jack_file);
    free(E2_0_A1);
    
    //c++ 81 || r 82
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    double *E2_1_A1=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,40,"E2_1_A1", shift_and_M_eff_sinh_T,jack_file);
    free(E2_1_A1);
    
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=C2_diff_masses;
    
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    //c++ 82 || r 83
    fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  41,0/*reim*/ , "E2_01_A1",  fit_info ,jack_file);
    free_fit_result(fit_info,fit_out);
       
    //c++ 83 || r 84
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    double *E2_0_E1=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,42,"E2_0_E1", shift_and_M_eff_sinh_T,jack_file);
    free(E2_0_E1);
    
    //c++ 84 || r 85
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    double *E2_1_E1=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,43,"E2_1_E1", shift_and_M_eff_sinh_T,jack_file);
    free(E2_1_E1);
    
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=C2_diff_masses;
    
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    //c++ 85 || r 86
    fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  44,0/*reim*/ , "E2_01_E1",  fit_info ,jack_file);
    
    //c++ 86 || r 87
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    double *E2_0_E2=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,45,"E2_0_E2", shift_and_M_eff_sinh_T,jack_file);
    free(E2_0_E2);
    
    //c++ 87 || r 88
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    double *E2_1_E2=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,46,"E2_1_E2", shift_and_M_eff_sinh_T,jack_file);
    free(E2_1_E2);
    
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=C2_diff_masses;
    
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    //c++ 88 || r 89
    fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  47,0/*reim*/ , "E2_01_E2",  fit_info ,jack_file);
    
    
} //if ncor>33   
else {    for(int i=74;i < 89;i++ )  fwrite(zeros,sizeof(double),Njack, jack_file );}

} //if ncorr>32
else {    for(int i=72;i < 89;i++ )  fwrite(zeros,sizeof(double),Njack, jack_file );}
    
} //if ncorr>30  
else {    for(int i=70;i < 89;i++ )  fwrite(zeros,sizeof(double),Njack, jack_file );}
}//if ncorr>15  
else {    for(int i=21;i < 89;i++ )  fwrite(zeros,sizeof(double),Njack, jack_file );}

    //E2_01_div_shift
    fit_info.Nvar=3;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=constant_fit;

    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    //c++ 89 || r 90
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   lhs_E2_div_shift<11>, "E2_01_div_shift",  fit_info, jack_file );
   
    free(a);
    //a=scattering_len_luscher(  Njack,  mass[0], mass[1], fit_out.P[0] ,params.data.L[1]);
    a=(double*) calloc(Njack,sizeof(double));
    tmpj=(double*) malloc(sizeof(double)*Njack);
    tmp_muj=(double*) malloc(sizeof(double)*Njack);
    sub_jackboot(Njack,  tmpj, fit_out.P[0], mass[0] );
    sub_jackboot(Njack,  tmpj, tmpj, mass[1] );
    fprintf(outfile,"#scattering length  a  err deltaE2 err    mu01 err    deltaE2*mu01  err   a_01*(m0+m1)=-mu/2pi  err      a_01*(m0+m1)/(a0m0)=4/3 a_01*pi/mu   err    \n %g  %g     %g  %g\t",
            a[Njack-1], error_jackboot(option[4],Njack,a),   tmpj[Njack-1],  error_jackboot(option[4],Njack,tmpj));

    //reduced mass
    for (j=0; j<Njack;j++)
        tmp_muj[j]=mass[0][j]* mass[1][j]/(mass[0][j]+ mass[1][j]);
    fprintf(outfile,"%g   %g\t", tmp_muj[Njack-1], error_jackboot(option[4],Njack,tmp_muj) );

    //reduced mass time DeltaE2
    for (j=0; j<Njack;j++)
        tmp_muj[j]=tmp_muj[j]*tmpj[j];
    fprintf(outfile,"%g   %g\t", tmp_muj[Njack-1], error_jackboot(option[4],Njack,tmp_muj)  );


    for (j=0; j<Njack;j++)
        tmp_muj[j]=a[j]*(mass[0][j]+mass[1][j]);
    fprintf(outfile,"%g   %g\t", tmp_muj[Njack-1], error_jackboot(option[4],Njack,tmp_muj)  );

    for (j=0; j<Njack;j++)
        tmp_muj[j]=(a[j]*(mass[0][j]+mass[1][j]))/a0m0[j];
    fprintf(outfile,"%g   %g\t", tmp_muj[Njack-1], error_jackboot(option[4],Njack,tmp_muj)  );

    // a_01*pi/mu
    for (j=0; j<Njack;j++)
        tmp_muj[j]=a[j]*pi_greco*(mass[0][j]+mass[1][j])/(mass[0][j]*mass[1][j]) ;
    fprintf(outfile,"%g   %g\n", tmp_muj[Njack-1], error_jackboot(option[4],Njack,tmp_muj)  );


    free(tmpj); free(tmp_muj);

    free_fit_result(fit_info,fit_out);
    fflush(outfile);
    //////////////////////////////////////////////////////77
if (params.data.ncorr>48){
    //c++ 90 || r 91
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    double *E2_0_A1=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,48,"E2_0_A1E1", shift_and_M_eff_sinh_T,jack_file);
    free(E2_0_A1);
    
    //c++ 91 || r 92
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    double *E2_1_A1=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,49,"E2_1_A1E1", shift_and_M_eff_sinh_T,jack_file);
    free(E2_1_A1);
    
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=C2_diff_masses;
    
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    //c++ 92 || r 93
    fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  50,0/*reim*/ , "E2_01_A1E1",  fit_info ,jack_file);
    free_fit_result(fit_info,fit_out);
    
}else {    for(int i=90;i < 93;i++ )  fwrite(zeros,sizeof(double),Njack, jack_file );}

if (params.data.ncorr>65){
    //c++ 93 || r 94
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    double *E2_0_A1E2=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,54,"E2_0_A1E2", shift_and_M_eff_sinh_T,jack_file);
    free(E2_0_A1E2);
    
    //c++ 94 || r 95
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    double *E2_1_A1E2=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,55,"E2_1_A1E2", shift_and_M_eff_sinh_T,jack_file);
    free(E2_1_A1E2);
    
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=C2_diff_masses;
    
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    //c++ 95 || r 96
    fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  56,0/*reim*/ , "E2_01_A1E2",  fit_info ,jack_file);
    free_fit_result(fit_info,fit_out);
    
    
    
    ///////////////////////////////////////////////   A1o20   /////////////////////////////////////////////
    //c++ 96 || r 97
    file_head.k[2]=mu1;    file_head.k[3]=mu1;
    double *E2_0_A1o20=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,60,"E2_0_A1o20", shift_and_M_eff_sinh_T,jack_file);
    free(E2_0_A1o20);
    
    //c++ 97 || r 98
    file_head.k[2]=mu2;    file_head.k[3]=mu2;
    double *E2_1_A1o20=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,61,"E2_1_A1o20", shift_and_M_eff_sinh_T,jack_file);
    free(E2_1_A1o20);
    
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=C2_diff_masses;
    
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=mass[1];
    //c++ 98 || r 99
    fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  62,0/*reim*/ , "E2_01_A1o20",  fit_info ,jack_file);
    free_fit_result(fit_info,fit_out);
    
    
    
}else {    for(int i=93;i < 99;i++ )  fwrite(zeros,sizeof(double),Njack, jack_file );}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

if (params.data.ncorr>74){
    
    
    fit_info.Nvar=1;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=0;
    fit_info.function=constant_fit;
    
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    //c++ 99 || r 100
    //fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  66,0/*reim*/ , "E2_0_px",  fit_info ,jack_file);
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   m_eff_of_sum<33,35,37>, "E1_0_p1",  fit_info, jack_file );
    
    //free_fit_result(fit_info,fit_out);
    double *E1_0_p1=fit_out.P[0];
    
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=C2_diff_masses;
    
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=E1_0_p1;
    //c++ 100 || r 101
    //fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  66,0/*reim*/ , "E2_0_px",  fit_info ,jack_file);
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   sum_corr_directions_shift<66,67,68>, "E2_0_p1",  fit_info, jack_file );
    fprintf(outfile,"#E2_CM   err  E2_CM/M   err \n " );
    double *E2_CM=(double*) malloc(sizeof(double)*Njack);
    for (int j=0;j< Njack;j++)
        E2_CM[j]=energy_CM(fit_out.P[0][j],{1,0,0},params.data.L[1]);
    fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
    for (int j=0;j< Njack;j++)
        E2_CM[j]/=mass[0][j];
    fprintf(outfile,"%.12g  %.12g  \n ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
    free(E2_CM);
    
    
    
    free_fit_result(fit_info,fit_out);
    
    free(E1_0_p1);
    
}else {    for(int i=99;i < 101;i++ )  fwrite(zeros,sizeof(double),Njack, jack_file );}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////p=(1,1,0)
if (params.data.ncorr>90){
    
    fit_info.Nvar=1;
    fit_info.Npar=1;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=0;
    fit_info.function=constant_fit;
    
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    //c++ 101 || r 102
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   m_eff_of_sum<75,76,77>, "E1_0_p11",  fit_info, jack_file );
    
    double *E1_0_p11=fit_out.P[0];
    
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=C2_diff_masses;
    
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=E1_0_p11;
    //c++ 102 || r 103
    fit_out=fit_fun_to_fun_of_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile, 
                                   sum_corr_directions_shift<81,82,83>, "E2_0_p11",  fit_info, jack_file );
    fprintf(outfile,"#E2_CM   err  E2_CM/M   err \n " );
    double *E2_CM=(double*) malloc(sizeof(double)*Njack);
    for (int j=0;j< Njack;j++)
        E2_CM[j]=energy_CM(fit_out.P[0][j],{1,1,0},params.data.L[1]);
    fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
    for (int j=0;j< Njack;j++)
        E2_CM[j]/=mass[0][j];
    fprintf(outfile,"%.12g  %.12g  \n ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
    free(E2_CM);
    free_fit_result(fit_info,fit_out);
    
    free(E1_0_p11);
    
}else {    for(int i=101;i < 103;i++ )  fwrite(zeros,sizeof(double),Njack, jack_file );}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////p=(1,1,1)
if (params.data.ncorr>90){
    
    file_head.k[2]=mu1;
    file_head.k[3]=mu1;
    //c++ 103 || r 74
    double *E1_0_p111=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,namefile_plateaux,outfile,90,"E1_0_p111", M_eff_T,jack_file);
    
    
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.n_ext_P=2;
    fit_info.function=C2_diff_masses;
    
    file_head.k[2]=mu1;    file_head.k[3]=mu2;
    fit_info.ext_P[0]=mass[0];
    fit_info.ext_P[1]=E1_0_p111;
    //c++ 104 || r 105
    fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,namefile_plateaux, outfile,  92,0/*reim*/ , "E2_0_p111",  fit_info ,jack_file);
    fprintf(outfile,"#E2_CM   err  E2_CM/M   err \n " );
    double *E2_CM=(double*) malloc(sizeof(double)*Njack);
    for (int j=0;j< Njack;j++)
        E2_CM[j]=energy_CM(fit_out.P[0][j],{1,1,1},params.data.L[1]);
    fprintf(outfile,"%.12g  %.12g  \t ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
    for (int j=0;j< Njack;j++)
        E2_CM[j]/=mass[0][j];
    fprintf(outfile,"%.12g  %.12g  \n ", E2_CM[Njack-1],error_jackboot(option[4],Njack,E2_CM ) );
    free(E2_CM);
    free_fit_result(fit_info,fit_out);
    
    free(E1_0_p111);
    
}else {    for(int i=103;i < 105;i++ )  fwrite(zeros,sizeof(double),Njack, jack_file );}


    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    free_2(3,mass);free_2(3,E2);
    free_corr(Neff, var, file_head.l0 ,data_bin);
    free_jack(Njack,var , file_head.l0, conf_jack);
 
    fclose(out_gamma);
    free(fit_info.ext_P);
    free_tower(7,(void**)option);
    
    for(i=0;i<file_head.nmoms;i++) 
    	free(file_head.mom[i]);
    free(file_head.mom);
    fclose(infile); fclose(outfile);
   
    free(file_head.k);
    fclose(jack_file);
    double z[2]={1,2};
    double q2=0;
    int dvec[3]={0,0,1};
    //double abc=dintegrand_12( q2, z);
    int aaa= dzeta_function(z,  q2,0 , 0, dvec, 0.1, 0.1, 0.1, 0.1,100);
    return 0;   
}

