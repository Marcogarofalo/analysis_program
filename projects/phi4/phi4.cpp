 
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

#include <cstring> 
#include <string>
#include <fstream>
// using namespace std;

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


/*
void read_nconfs(int *s, int *c, FILE *stream){

   FILE *f1;
   long int tmp;

   fread(s,sizeof(int),1,stream);
   f1=stream;
   
   fseek(stream, 0, SEEK_END);
   tmp = ftell(stream);
   tmp-= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*9 ;
   tmp-= sizeof(int)*2;
   (*c)= tmp/ (sizeof(int)+ (*s)*sizeof(double) );
 
  rewind(stream);
  read_file_head_bin(stream);
  fread(s,sizeof(int),1,stream);

}
*/

double *constant_fit(int M, double in){
    double *r;
    
    r=(double*) malloc(sizeof(double)*M);
    r[0]=1.;
    
    return r;
}

double matrix_element_GEVP(int t, double **cor,double mass){
      double me;
      
      me=cor[t][0]/sqrt( exp(-mass*t)+exp(-(file_head.l0-t)*mass) );
      me*=2*mass;

      return  me;
}



int index_minus_r( int r)
{
    int mr;
    if (r==0) mr= 1;
    else if (r==1) mr= 0;
    else error(0==0,0,"index_minus_r","r is not 0 or 1\n");
    return mr;
}

static int index_twopt(int si,int ii,int ix0,int imom2,int imom1,int ik2,int r2,int ik1,int r1)
{
    int nk,nmoms;

    nk=file_head.nk;
    nmoms=file_head.nmoms;

    return ii+si*(ix0+file_head.l0*(imom2+nmoms*(imom1+nmoms*(mass_index[ik2][r2][ik1][r1]))));
}
static int index_threept(int si,int ii,int ix0,int imom1,int imom2,int ik1,int ik2,int ik3)
{
    int nk,nmoms;

    nk=file_head.nk;
    nmoms=file_head.nmoms;

    return ii+si*(ix0+file_head.l0*(imom1+nmoms*(imom2+nmoms*(ik1+nk*(ik2+nk*ik3)))));
}
int index_minus_kappa(int ik)
{
    int imk,i;
    double mu;
    
    mu=-file_head.k[ file_head.nk+ik ];
    imk=-1;
    for (i=0;i<file_head.nk;i++){
	if ( file_head.k[ file_head.nk+i ]==mu )
            imk=i;
    }

   error(imk==-1,1,"inde_minus_kappa ",  "Unable to find mass=%g",mu);
   return imk; 


}

static int index_minus_theta(int imom)
{
   int i,imth;
   double m0,m1,m2,m3;

   imth=-1;
   m0= file_head.mom[imom][0];
   m1= -file_head.mom[imom][1];
   m2= -file_head.mom[imom][2];
   m3= -file_head.mom[imom][3];
   for(i=0;i<file_head.nmoms;++i)
      if(m0==file_head.mom[i][0] && m1==file_head.mom[i][1] && 
	 m2==file_head.mom[i][2] && m3==file_head.mom[i][3])
	 imth=i;

   error(imth==-1,1,"inde_minus_theta ",  "Unable to find theta=%g",m1);
   return imth;
}


void read_twopt(FILE *stream,int size, int iconf , double **to_write,int si, int ii, int imom2, int imom1, int ik2, int r2, int ik1,int r1 ){
   
   long int tmp;
   int iiconf,N,s;
   double *obs;
   int mik1, mik2;
   int mimom1,mimom2,mr1,mr2;
   int t,vol,index;
   double re,im;
 
   tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(12) ;
   tmp+=sizeof(double)*iconf*size+sizeof(int)*iconf;
   fseek(stream, tmp, SEEK_SET);

   obs=(double*) malloc(size*sizeof(double)); 
   fread(obs,sizeof(double),size,stream);   
   
   mimom1=index_minus_theta(imom1);
   mimom2=index_minus_theta(imom2);
   mr1=index_minus_r(r1);
   mr2=index_minus_r(r2);

   for(t=0;t<file_head.l0;t++){
	   re=0;vol=0;im=0;
	   index=2*index_twopt(si,ii,t,imom2,imom1,ik2,r2,ik1,r1);
	   re+= obs[index];
       im+= obs[index+1];
       vol++;
       index=2*index_twopt(si,ii,t,imom2,imom1,ik2,mr2,ik1,mr1);
	   re+= obs[index];
       im+= obs[index+1];
       vol++;
	   to_write[t][0]=re/( (double) vol );
	   to_write[t][1]=im/( (double) vol );
   }
   free(obs);
}


 
void extract_threept(double *to_read , double **to_write,int si, int ii, int imom1, int imom2, int ik1, int ik2 ,int ik3,int sym ){

   int mik1, mik2,mik3;
   int mimom1,mimom2;
   int t,vol,index;
   double re,im;

   double symm=(double) sym;
   mik1=index_minus_kappa(ik1);
   mik2=index_minus_kappa(ik2);
   mik3=index_minus_kappa(ik3);
   mimom1=index_minus_theta(imom1);
   mimom2=index_minus_theta(imom2);

	for(t=0;t<file_head.l0;t++){
	   re=0;vol=0;im=0;
	   index=2*index_threept(si,ii,t,imom1,imom2,ik1,ik2,ik3);
	   re+= to_read[index];
	   im+= to_read[index+1];
           vol++;
	   index=2*index_threept(si,ii,t,imom1,imom2,mik1,mik2,mik3);
	   re+= to_read[index];
	   im+= to_read[index+1];
           vol++;
        /*   if (  ik1 ==ik2){
		   index=2*index_threept(si,ii,t,imom1,imom2,mik1,ik1,mik3);
		   re+= to_read[index];
                   im+= to_read[index+1];
                   vol++;
		   index=2*index_threept(si,ii,t,imom1,imom2,ik3,mik3,ik1);
		   re+= to_read[index];
	           im+= to_read[index+1];
		   vol++;
           }*/
  

	   if (  mimom1>=0 &&  mimom2 >=0 )
	   { 
		   index=2*index_threept(si,ii,t,mimom1,mimom2,ik1,ik2,ik3);
		   re+=symm* to_read[index];
	           im+=symm* to_read[index+1];
                   vol++;
		   index=2*index_threept(si,ii,t,mimom1,mimom2,mik1,mik2,mik3);
		   re+=symm* to_read[index];
	           im+=symm* to_read[index+1];
		   vol++;
		/*   if (  ik1 ==ik2){
			   index=2*index_threept(si,ii,t,mimom1,mimom2,mik3,ik3,mik1);
			   re+=symm* to_read[index];
	           	   im+= symm*to_read[index+1];
	                   vol++;
			   index=2*index_threept(si,ii,t,mimom1,mimom2,ik3,mik3,ik1);
			   re+=symm* to_read[index];
            		   im+=symm* to_read[index+1];
		   	   vol++;
		   }*/
	   }
	   to_write[t][0]=re/( (double) vol );
	   to_write[t][1]=im/((double) vol);
	}
}

/*
double *jack_mass(int tmin, int tmax, int sep ,double **corr_ave, double **corr_J,int Njack ){
    
   double ***y,*x,**tmp,*fit; 
   int i,j;  

   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(tmax-tmin+1));
        for (i=tmin;i<=tmax;i++){
            y[j][i-tmin]=(double*) malloc(sizeof(double)*2);
        }
   }
   x=(double*) malloc(sizeof(double)*(tmax-tmin+1));
   fit=(double*) malloc(sizeof(double)*Njack);

   for (i=tmin;i<=tmax;i+=sep){
        for (j=0;j<Njack;j++){
            y[j][(i-tmin)/sep][0]=corr_J[i][j];
            y[j][(i-tmin)/sep][1]=corr_ave[i][1];
        }
    }
    for (j=0;j<Njack;j++){
        tmp=linear_fit( (tmax-tmin)/sep +1, x, y[j],  1,constant_fit_to_try );
        fit[j]=tmp[0][0];
        free(tmp[0]);free(tmp);
    }    
    return fit;    
}
*/
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

   
   
   
   error(argc!=15,1,"main ",
         "usage:./phi4  blind/see/read_plateaux path T L  msq0 msq1 l0 l1 mu g -bin $bin  jack/boot ");
   error(strcmp(argv[1],"blind")!=0 &&  strcmp(argv[1],"see")!=0 && strcmp(argv[1],"read_plateaux")!=0   ,1,"main ",
         "argv[1] only options:  blind/see/read_plateaux ");
   error(strcmp(argv[12],"-bin")!=0 ,1,"main", "argv[12] must be: -bin");
   error(strcmp(argv[14],"jack")!=0 &&  strcmp(argv[14],"boot")!=0,1,"main",
         "argv[13] only options: jack/boot");
   
   
    char **option;
    option=(char **) malloc(sizeof(char*)*6);
    option[0]=(char *) malloc(sizeof(char)*NAMESIZE);
    option[1]=(char *) malloc(sizeof(char)*NAMESIZE);
    option[2]=(char *) malloc(sizeof(char)*NAMESIZE);
    option[3]=(char *) malloc(sizeof(char)*NAMESIZE);
    option[4]=(char *) malloc(sizeof(char)*NAMESIZE);
    option[5]=(char *) malloc(sizeof(char)*NAMESIZE);

    mysprintf(option[1],NAMESIZE,argv[1]); // blind/see/read_plateaux
    mysprintf(option[2],NAMESIZE,"-p"); // -p
    mysprintf(option[3],NAMESIZE,argv[2]); // path
    mysprintf(option[4],NAMESIZE,argv[14]); //resampling
    mysprintf(option[5],NAMESIZE,"no"); // pdf
    printf("resampling %s\n",option[4] );
    int T=atoi(argv[3]);
   
   double mu1=atof(argv[5]);
   double mu2=atof(argv[6]);
   printf("mu=%g  %g\n",mu1,mu2);
   char namefile[NAMESIZE];
   //mysprintf(argv[4],NAMESIZE,"jack");
   
   file_head.l0=T;
   file_head.l1=T/2;file_head.l2=T/2;file_head.l3=T/2;
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
             argv[2], atoi(argv[3]), atoi(argv[4]), atof(argv[5]), atof(argv[6]),
             atof(argv[7]), atof(argv[8]), atof(argv[9]), atof(argv[10]),atoi(argv[11]) );
   FILE *outfile=open_file(namefile,"w+");      
   
   mysprintf(namefile,NAMESIZE,"%s/out/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d_gamma",
             argv[2], atoi(argv[3]), atoi(argv[4]), atof(argv[5]), atof(argv[6]),
             atof(argv[7]), atof(argv[8]), atof(argv[9]), atof(argv[10]),atoi(argv[11]) );
   FILE *out_gamma=open_file(namefile,"w+");      
   
   // open infile and count the lines
   //
   mysprintf(namefile,NAMESIZE,"%s/G2t_T%d_L%d_msq0%.6f_msq1%.6f_l0%.6f_l1%.6f_mu%.6f_g%.6f_rep%d",
             argv[2], atoi(argv[3]), atoi(argv[4]), atof(argv[5]), atof(argv[6]),
             atof(argv[7]), atof(argv[8]), atof(argv[9]), atof(argv[10]),atoi(argv[11]));
   printf("opening file: %s \n", namefile);
   FILE *infile=open_file(namefile,"r");
   
   
   int count=0;
   string line;
   ifstream file( namefile);
   while (getline(file, line))
        count++;
 
   std::cout << "Numbers of lines in the file : " << count << std::endl;
   confs=count/T;
   std::cout << "Numbers of configurations in the file : " << confs << std::endl;
   
   // compute what will be the neff after the binning 
   int bin=atoi(argv[13]);
   int Neff=confs/bin;
   std::cout << "effective configurations after binning (" << bin  <<"):  "<<Neff << std::endl;

   int Njack;
   if( strcmp(argv[14],"jack")==0)
                Njack=Neff+1;
   else if( strcmp(argv[14],"boot")==0)
                Njack=Nbootstrap+1;
   else
       error(1==1,1,"main","argv[14]= %s is not jack or boot",argv[14]);
   
   int var=4;
   data=calloc_corr(confs, var,  file_head.l0 );
   
   setup_file_jack(option,Njack);
    
   get_kinematic( 0,0,  1, 0,0,  0 );
   printf("option[4]=%s\n",option[4]);

   for (int iconf=0; iconf< confs ;iconf++){
       for (int t =0; t< T ;t++){
           int tt;
           //fscanf(infile,"%d  %lf",&tt,&data[iconf][0][t][0]);
           //error(t!=tt, 1, "main: reading","time do not match  conf=%d   t=%d  read %d",iconf ,t,tt);
           double a ,b,c;
           fscanf(infile,"%d  %lf %lf",&tt,&data[iconf][0][t][0],&data[iconf][1][t][0]);
           fscanf(infile,"%lf %lf\n",&data[iconf][2][t][0],&data[iconf][3][t][0]);
           error(t!=tt, 1, "main: reading","time do not match  conf=%d   t=%d  read %d",iconf ,t,tt);
           
           //fscanf(infile,"%lf",&data[iconf][0][t][0]);
           //data[iconf][0][t][0]*=-1;
       }
    }

    symmetrise_corr(confs, 0, file_head.l0,data);
    symmetrise_corr(confs, 1, file_head.l0,data);
    symmetrise_corr(confs, 2, file_head.l0,data);
    symmetrise_corr(confs, 3, file_head.l0,data);
    
    //if you want to do the gamma analysis you need to do before freeing the raw data
    effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data,  confs ,&plateaux_masses,out_gamma,0,"M_{PS}^{ll}");
    //effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data,  confs ,&plateaux_masses,out_gamma,3,"M_{PS}^{ll}");
    data_bin=binning(confs, var, file_head.l0 ,data, bin);
    free_corr(confs, var, file_head.l0 ,data);
    
    conf_jack=create_resampling(option[4],Neff, var, file_head.l0, data_bin);
   
    double *mass;
    fprintf(outfile,"#correlator\n");
    for (int t =1; t< T/2 ;t++){
        for (int v=0;v< var;v++){
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
    fprintf(outfile,"\n\n #%s fit in [%d,%d] chi2=%.5f\n  %.15g    %.15g    %d   %d\n\n\n","#need_for_gnuplot",0,0,0.0,0.0,0.0,0,0);
    
   
    
  
    file_head.k[2]=mu1;
    file_head.k[3]=mu1;

    
    mass=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,&plateaux_masses,outfile,0,"E1", M_eff_T);
    //mass=compute_effective_mass(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,&plateaux_masses,outfile,0,"M_{PS}^{ll}");
    free(mass);
    
    file_head.k[2]=mu2;
    file_head.k[3]=mu2;
    
    mass=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,&plateaux_masses,outfile,1,"E1", M_eff_T);
    //compute_effective_mass(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,&plateaux_masses,outfile,1,"M_{PS}^{ll}");

    //Zfpi=compute_Zf_PS_ll(  option, kinematic_2pt, (char*) "P5P5", conf_jack, mass,  Njack ,plateaux_masses,outfile );
    
     double *E2;
    
    E2=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,&plateaux_masses,outfile,2,"E2", two_particle_energy);
    
    //E2=plateau_correlator_function(  option, kinematic_2pt,   (char*) "P5P5", conf_jack,  Njack ,&plateaux_masses,outfile,3,"E2", two_particle_energy);
    //free(E2);
    
    struct fit_type fit_info;
    struct fit_result  fit_out;
    fit_info.Nvar=3;
    fit_info.Npar=3;
    fit_info.N=1;
    fit_info.Njack=Njack;
    fit_info.function=C3;
    fit_info.n_ext_P=2;
    fit_info.ext_P=(double**) malloc(sizeof(double*)*2);
    fit_info.ext_P[0]=mass;
    fit_info.ext_P[1]=E2;
    
    fit_out=fit_function_to_corr(option , kinematic_2pt ,  (char*) "P5P5", conf_jack ,&plateaux_masses, outfile,  3,0/*reim*/ , "E3",  fit_info, file_jack.M_PS );
    
    
    free(mass);free(E2);
    free_corr(Neff, var, file_head.l0 ,data_bin);
    free_jack(Njack,var , file_head.l0, conf_jack);
 
    fclose(out_gamma);
    return 0;   
}

