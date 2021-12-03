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



struct  kinematic kinematic_2pt;

int r_value(int r)
{
    int vr;
    if (r==0) vr= 1;
    else if (r==1) vr= -1;
    else error(0==0,0,"r_value","r value is not 0 neither 1\n");
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

static void  read_file_head_bin(FILE *stream)
{
    int i,dsize;
    double *dstd;
    
    fread(&file_head.twist,sizeof(int),1,stream);
    fread(&file_head.nf,sizeof(int),1,stream);
    fread(&file_head.nsrc,sizeof(int),1,stream);
    fread(&file_head.l0,sizeof(int),1,stream);
    fread(&file_head.l1,sizeof(int),1,stream);
    fread(&file_head.l2,sizeof(int),1,stream);
    fread(&file_head.l3,sizeof(int),1,stream);
    fread(&file_head.nk,sizeof(int),1,stream);
    fread(&file_head.nmoms,sizeof(int),1,stream);
    
    fread(&file_head.beta,sizeof(double),1,stream);
    fread(&file_head.ksea,sizeof(double),1,stream);
    fread(&file_head.musea,sizeof(double),1,stream);
    fread(&file_head.csw,sizeof(double),1,stream);
   
    file_head.k=(double*) malloc(sizeof(double)*2*file_head.nk);
    for(i=0;i<2*file_head.nk;++i)
    	fread(&file_head.k[i],sizeof(double),1,stream);
    
    file_head.mom=(double**) malloc(sizeof(double*)*file_head.nmoms);
    for(i=0;i<file_head.nmoms;i++) {
    	file_head.mom[i]=(double*) malloc(sizeof(double)*4);
        fread(&file_head.mom[i][0],sizeof(double),1,stream);
        fread(&file_head.mom[i][1],sizeof(double),1,stream);
        fread(&file_head.mom[i][2],sizeof(double),1,stream);
        fread(&file_head.mom[i][3],sizeof(double),1,stream);

    }
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


void read_twopt(FILE *stream,int size, int iconf , double **to_write,int si, int ii, int imom2, int imom1, int ik2, int r2, int ik1,int r1,
int sym ){
   
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
	   re+=((double) sym) * obs[index];
       im+=((double) sym) * obs[index+1];
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
     mysprintf(save_name,NAMESIZE,"%s/%s",argv[3],name);
     f=fopen(save_name,"w+");
     error(f==NULL,1,"setup_file_jack ",
         "Unable to open output file %s/%s",argv[3],name);
     write_file_head(f);
     fwrite(&Njack,sizeof(int),1,f);
     fclose(f);
}

void setup_single_file_jack_ASCI(char  *save_name, char **argv,const char  *name,int Njack){
     FILE *f;
     mysprintf(save_name,NAMESIZE,"%s/%s",argv[3],name);
     f=fopen(save_name,"w+");
     error(f==NULL,1,"setup_file_jack ",
         "Unable to open output file %s/%s",argv[3],name);
     fclose(f);
}

void setup_file_jack(char **argv,int Njack){
    if( strcmp(argv[4],"jack")==0){
     setup_single_file_jack(file_jack.mpcac,argv,"jackknife/mpcac_jack",Njack);
     setup_single_file_jack(file_jack.M_PS,argv,"jackknife/M_{PS}_jack",Njack);
     setup_single_file_jack(file_jack.f_PS,argv,"jackknife/f_{PS}_jack",Njack);
     setup_single_file_jack(file_jack.Zf_PS,argv,"jackknife/Zf_{PS}_jack",Njack);
    
     setup_single_file_jack(file_jack.M_PS_GEVP,argv,"jackknife/M_{PS}^{GEVP}_jack",Njack);
     setup_single_file_jack(file_jack.f_PS_ls_ss,argv,"jackknife/f_{PS}_ls_ss_jack",Njack);

    }
               
    if( strcmp(argv[4],"boot")==0){
     setup_single_file_jack(file_jack.mpcac,argv,"jackknife/mpcac_boot",Njack);
     
     setup_single_file_jack(file_jack.M_PS,argv,"jackknife/M_{PS}_boot",Njack);
     setup_single_file_jack(file_jack.f_PS,argv,"jackknife/f_{PS}_boot",Njack);
     setup_single_file_jack(file_jack.Zf_PS,argv,"jackknife/Zf_{PS}_boot",Njack);
     
     setup_single_file_jack(file_jack.M_PS_GEVP,argv,"jackknife/M_{PS}^{GEVP}_boot",Njack);
     setup_single_file_jack(file_jack.f_PS_ls_ss,argv,"jackknife/f_{PS}_ls_ss_boot",Njack);
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

   error(argc!=6,1,"main ",
         "usage: ./form_factors blind/see/read_plateaux -p inpath   jack/boot pdf");

   error(strcmp(argv[2],"-p")!=0,1,"main ",
         "missing -p \n usage: ./form_factors blind/see/read_plateaux -p inpa.hpp");

   error(strcmp(argv[4],"jack")!=0 && strcmp(argv[4],"boot")!=0 ,2,"main ",
         "choose jack or boot \n usage: ./form_factors blind/see/read_plateaux -p inpath   jack/boot");
 
   char namefile[NAMESIZE];
   
   
   FILE *outfile        =NULL; mysprintf(namefile,NAMESIZE,"%s/out/out_E0.txt",argv[3]);        outfile=fopen(namefile,"w+");       error(outfile==NULL,1,"main ", "Unable to open %s file",namefile);
   FILE *outfile_ls        =NULL; mysprintf(namefile,NAMESIZE,"%s/out/out_E0_ls.txt",argv[3]);        outfile_ls=fopen(namefile,"w+");       error(outfile_ls==NULL,1,"main ", "Unable to open %s file",namefile);
   FILE *outfile_ss        =NULL; mysprintf(namefile,NAMESIZE,"%s/out/out_E0_ss.txt",argv[3]);        outfile_ss=fopen(namefile,"w+");       error(outfile_ss==NULL,1,"main ", "Unable to open %s file",namefile);
   
   FILE *m_pcac         =NULL; mysprintf(namefile,NAMESIZE,"%s/out/m_pcac.txt",argv[3]);        m_pcac=fopen(namefile,"w+");        error(m_pcac==NULL,1,"main ", "Unable to open %s file",namefile);
   FILE *outfile_RM     =NULL; mysprintf(namefile,NAMESIZE,"%s/out/Ratio_masses.txt",argv[3]);  outfile_RM=fopen(namefile,"w+");    error(outfile_RM==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_GEVP   =NULL; mysprintf(namefile,NAMESIZE,"%s/out/out_E0_GEVP.txt",argv[3]);       outfile_GEVP=fopen(namefile,"w+");  error(outfile_GEVP==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_RM_GEVP=NULL; mysprintf(namefile,NAMESIZE,"%s/out/Ratio_masses_GEVP.txt",argv[3]);  outfile_RM_GEVP=fopen(namefile,"w+");    error(outfile_RM_GEVP==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_f      =NULL; mysprintf(namefile,NAMESIZE,"%s/out/f_PS.txt",argv[3]);          outfile_f=fopen(namefile,"w+");     error(outfile_f==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_f_ls_ss      =NULL; mysprintf(namefile,NAMESIZE,"%s/out/f_PS_ls_ss.txt",argv[3]);          outfile_f_ls_ss=fopen(namefile,"w+");     error(outfile_f_ls_ss==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_f_GEVP =NULL; mysprintf(namefile,NAMESIZE,"%s/out/f_PS_GEVP.txt",argv[3]);     outfile_f_GEVP=fopen(namefile,"w+");error(outfile_f_GEVP==NULL,1,"main ", "Unable to open %s file",namefile);
   FILE *outfile_Rf     =NULL; mysprintf(namefile,NAMESIZE,"%s/out/Ratio_f_PS.txt",argv[3]);    outfile_Rf=fopen(namefile,"w+");    error(outfile_Rf==NULL,1,"main ",  "Unable to open %s file",namefile);
   
 
   
   double E_B,E_Pi, x_SCHET,q2,vec_pB,vec_pPi;
   int Neff,Njack,bin=1;
   
   t1=clock();
   if ( strcmp(argv[1],"read_plateaux")==0 ){
      //mysprintf(namefile,NAMESIZE,"%s/plateaux_masses.txt",argv[3]);
      mysprintf(kinematic_2pt.plateau_m_ll,NAMESIZE,"%s/plateaux_masses.txt",argv[3]);
      plateaux_masses=     fopen(kinematic_2pt.plateau_m_ll,"r");          error(plateaux_masses==NULL,1,"main ", "Unable to open %s file",kinematic_2pt.plateau_m_ll);
      
      //mysprintf(namefile,NAMESIZE,"%s/plateaux_masses_GEVP.txt",argv[3]);      
      mysprintf(kinematic_2pt.plateau_m_GEVP,NAMESIZE,"%s/plateaux_masses_GEVP.txt",argv[3]); 
      plateaux_masses_GEVP=fopen(kinematic_2pt.plateau_m_GEVP,"r");          error(plateaux_masses_GEVP==NULL,1,"main ", "Unable to open %s file",kinematic_2pt.plateau_m_GEVP);
      
      //mysprintf(namefile,NAMESIZE,"%s/plateaux_f.txt",argv[3]);
      mysprintf(kinematic_2pt.plateau_f,NAMESIZE,"%s/plateaux_f.txt",argv[3]);
      plateaux_f=        fopen(kinematic_2pt.plateau_f,"r");    error(plateaux_f==NULL,1,"main ", "Unable to open %s file",kinematic_2pt.plateau_f); 
        
  }
   // f=fopen("./meas_2pts_bin10.dat","r");
   mysprintf(namefile,NAMESIZE,"%s/data/to_read_ll_bin10.dat",argv[3]);
   f_ll=fopen(namefile,"r"); if (f_ll==NULL) {printf("2pt file not found\n"); exit(0);}
   mysprintf(namefile,NAMESIZE,"%s/data/to_read_sl_bin10.dat",argv[3]);
   f_sl=fopen(namefile,"r"); if (f_sl==NULL) {printf("2pt file not found\n"); exit(0);}
   mysprintf(namefile,NAMESIZE,"%s/data/to_read_ls_bin10.dat",argv[3]);
   f_ls=fopen(namefile,"r"); if (f_ls==NULL) {printf("2pt file not found\n"); exit(0);}
   mysprintf(namefile,NAMESIZE,"%s/data/to_read_ss_bin10.dat",argv[3]);
   f_ss=fopen(namefile,"r"); if (f_ss==NULL) {printf("2pt file not found\n"); exit(0);}
   
   read_file_head_bin(f_ll);
   print_file_head(outfile);
   print_file_head(outfile);
   print_file_head(outfile_ls);
   print_file_head(outfile_ss);
   print_file_head(outfile_GEVP);
  

   fflush(outfile);
   fflush(outfile_ls);
   fflush(outfile_ss);
   
   
   init_mass_index();

   printf("index mass=%d\n",mass_index[1][0][0][0]);
   read_nconfs(&size,&confs,f_ll);   
   int Nconfs;
   fread(&Nconfs,sizeof(int),1,f_ll );
////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//confs=200;
//printf("number of confs set manually to %d\n",confs);
////////////////////////////////////////////////////////////////////
   printf("size=%d  confs=%d  confs=%d\n",size,confs,Nconfs);

   Neff=confs/bin;
   printf("bi_size=%d   confs after binning=%d\n",bin,Neff);
   if( strcmp(argv[4],"jack")==0)
                Njack=Neff+1;
   if( strcmp(argv[4],"boot")==0)
                Njack=Nbootstrap+1;

   setup_file_jack(argv,Njack);
   iconf=(int*) malloc(sizeof(int)*confs);
   out=(double**) malloc(sizeof(double*)*confs);
   
   
   
  /* printf("reading confs:\n");
   for(i=0;i<confs;i++){
       fread(&iconf[i],sizeof(int),1,f_ll);
       out[i]=(double*) malloc(sizeof(double)*size);
       fread(out[i],sizeof(double),size,f_ll);
       printf("%d\t",iconf[i]);
   }
   printf("\n");
   */
 //  fclose(f);
   //copyng the input
   int var=2,si_2pt,si;
   
   
   data=calloc_corr(confs, var*var+2,  file_head.l0 );
   lambda0=calloc_corr(Njack, var*var,  file_head.l0 );
   projected_O=calloc_corr(Njack, var*var,  file_head.l0);
     
   int ii, imom1, imom2, ik1,ik2,ik3;
   ii=0;
   imom1=0;
   imom2=0;
   ik1=0;
   ik2=0;
   ik3=0;
   int tmin=18, tmax=23, sep=1;
   int yn;
  
   char name[NAMESIZE];
   int index;
   double **E0,***Mass;
   E0=(double**) malloc(sizeof(double*)*size);
   for (i=0;i<size ;i++)
       E0[i]=(double*) malloc(sizeof(double)*2); 

    Mass=(double***) malloc(sizeof(double**)*file_head.nk);
   for(ik2=0;ik2<file_head.nk;ik2++){
   	Mass[ik2]=(double**) malloc(sizeof(double*)*file_head.nk);
   	for(ik1=0;ik1<file_head.nk;ik1++){
   		Mass[ik2][ik1]=(double*) malloc(sizeof(double)*Njack);
        }
   }

   double ave=0;
   int vol=0;
   si_2pt=size/(4*file_head.nk*((file_head.nk+1)/2)*file_head.nmoms*file_head.nmoms*file_head.l0*2);
   
   printf("%d\n\n",si_2pt); 
   int mik3, mik2,mik1, mimom1,mimom2,counter,r1,r2,mr1,mr2;
  
   
   
   print_file_head(m_pcac); 
   
   
   double ***mass_jack_fit,***f_jack_fit,***mass_jack_fit_GEVP,***Zf_jack_fit;   
   mass_jack_fit=      (double***) malloc(sizeof(double**)*file_head.nk);
   mass_jack_fit_GEVP= (double***) malloc(sizeof(double**)*file_head.nk);
   f_jack_fit=         (double***) malloc(sizeof(double**)*file_head.nk);
   Zf_jack_fit=         (double***) malloc(sizeof(double**)*file_head.nk);
   for(ik2=0;ik2<file_head.nk;ik2++){
        mass_jack_fit[ik2]=     (double**) malloc(sizeof(double*)*file_head.nk);
        mass_jack_fit_GEVP[ik2]=(double**) malloc(sizeof(double*)*file_head.nk);
        f_jack_fit[ik2]=        (double**) malloc(sizeof(double*)*file_head.nk);
        Zf_jack_fit[ik2]=        (double**) malloc(sizeof(double*)*file_head.nk);
   }
   
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);
   mt=(double**) malloc(sizeof(double*)*file_head.l0);

   
   
   M=calloc_corr(file_head.l0, Njack,  var*var );
   vec=calloc_corr(file_head.l0, Njack,  var*var );
   lambda=calloc_corr(file_head.l0, Njack,  var*var );
   /*
   M=(double****) malloc(sizeof(double***)*file_head.l0);
   vec=(double****) malloc(sizeof(double***)*file_head.l0);
   lambda=(double****) malloc(sizeof(double***)*file_head.l0);
   for (i=0;i<file_head.l0;i++){
        M[i]=(double***) malloc(sizeof(double**)*Njack);
        vec[i]=(double***) malloc(sizeof(double**)*Njack);
        lambda[i]=(double***) malloc(sizeof(double**)*Njack);
        for (j=0;j<Njack;j++){
            M[i][j]=(double**) malloc(sizeof(double*)*var*var);
            vec[i][j]=(double**) malloc(sizeof(double*)*var*var);
            lambda[i][j]=(double**) malloc(sizeof(double*)*var*var);
            for(t=0;t<var*var;t++){
                M[i][j][t]=(double*) malloc(2*sizeof(double));
                vec[i][j][t]=(double*) malloc(2*sizeof(double));
                lambda[i][j][t]=(double*) malloc(2*sizeof(double));
	   }
        }
   }
   
     */
   
   int ik1_max=4;
   if (file_head.nk<4) ik1_max=file_head.nk;
//   for (ii=0;ii<si_2pt;ii++){
   for(ik1=0;ik1<ik1_max;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
   for(r1=0;r1<1;r1++){  // it is not a loop, r1=0
        if ( strcmp(argv[1],"read_plateaux")==0 ){
           // go_to_line(plateaux_masses,ik1);
            go_to_line(plateaux_masses_GEVP,ik1);  
           // if (ik1!=0) move_line(plateaux_f,1);
        }
   for(ik2=ik1;ik2<file_head.nk;ik2++){
   for(r2=0;r2<1;r2++){    // it is not a loop, r2=0
   for(imom2=0;imom2<file_head.nmoms;imom2++){
   for(imom1=0;imom1<file_head.nmoms;imom1++){
       mr2=index_minus_r(r2);
       mr1=index_minus_r(r1);
       mimom1=index_minus_theta(imom1);
       mimom2=index_minus_theta(imom2);
      

       /*   P5P5 mass*/ 
       contraction_index(&ii,"P5P5");
       contraction_name(ii,name);


       for (i=0;i<confs;i++){
      // extract_twopt(out[i] ,data[i][0],si_2pt,ii,imom2,imom1,ik2,r2,ik1,r1 );
          read_twopt(f_ll,size,i ,data[i][0],si_2pt,ii,imom2,imom1,ik2,r2,ik1,r1, 1 );
          read_twopt(f_ls,size,i ,data[i][1],si_2pt,ii,imom2,imom1,ik2,r2,ik1,r1, 1 );
          read_twopt(f_sl,size,i ,data[i][2],si_2pt,ii,imom2,imom1,ik2,r2,ik1,r1, 1 );
          read_twopt(f_ss,size,i ,data[i][3],si_2pt,ii,imom2,imom1,ik2,r2,ik1,r1, 1 );

          //ii=0 V0P5 , ii=1 P5P5
	      read_twopt(f_ll,size,i ,data[i][4],si_2pt,0,imom2,imom1,0,r2,0,r1,-1 );
          read_twopt(f_ll,size,i ,data[i][5],si_2pt,ii,imom2,imom1,0,r2,0,r1,1 );
     // if(  strcmp(name,"P5P5")==0 ) if(ik2==0) if (r2==0) if (ik1==0) if (r1==0) if (imom1==0) if (imom2==0) for(t=0;t<file_head.l0;t++){ 
       //   fprintf(petros,"%d   %.15f   %.15f\n",t,data[i][0][t][0],data[i][0][t][1]);}
       }
       symmetrise_corr(confs, 0, file_head.l0,data);
       symmetrise_corr(confs, 1, file_head.l0,data);
       symmetrise_corr(confs, 2, file_head.l0,data);
       symmetrise_corr(confs, 3, file_head.l0,data);

       antisymmetrise_corr(confs, 4, file_head.l0,data);
       forward_derivative_corr(confs, 4, file_head.l0,data);
       symmetrise_corr(confs, 5, file_head.l0,data);

       data_bin=binning(confs, var*var+2, file_head.l0 ,data, bin);

       conf_jack=create_resampling(argv[4],Neff, var*var+2, file_head.l0, data_bin);
      
       get_kinematic( ik2,r2,  ik1, r1,imom2,  imom1 );
       
       /*if you want to print the GEVP correlator
  get_kinematic( ik2,r2,  ik1, r1,imom2,  imom1 );
mass_jack_fit_GEVP[ik2][ik1]=compute_effective_mass_GEVP(  argv, kinematic_2pt,  name, data,  confs, var ,plateaux_masses_GEVP,outfile_GEVP );
  */  
       
       mass_jack_fit[ik2][ik1]=compute_effective_mass(  argv, kinematic_2pt,  name, conf_jack,  Njack ,&plateaux_masses,outfile,0,"M_{PS}^{ll}");
       
       double *tmp1=compute_mpcac(  argv, kinematic_2pt,  name, conf_jack,  Njack ,plateaux_masses,m_pcac,0);
       free(tmp1);
       
       double *tmp=compute_effective_mass(  argv, kinematic_2pt,  name, conf_jack,  Njack ,&plateaux_masses,outfile_ls,1,"M_{PS}^{ls}");
       free(tmp);
       tmp=compute_effective_mass(  argv, kinematic_2pt,  name, conf_jack,  Njack ,&plateaux_masses,outfile_ss,3,"M_{PS}^{ss}");
       free(tmp);
       //plateaux_masses=open_file(kinematic_2pt.plateau_m_ll,"r");
       
       mass_jack_fit_GEVP[ik2][ik1]=compute_effective_mass_GEVP(  argv, kinematic_2pt,  name, conf_jack,  Njack, var ,plateaux_masses_GEVP,outfile_GEVP );
       f_jack_fit[ik2][ik1]=compute_f_PS_ls_ss(  argv, kinematic_2pt,  name, conf_jack, mass_jack_fit_GEVP[ik2][ik1],  Njack ,plateaux_f,outfile_f_ls_ss );
       Zf_jack_fit[ik2][ik1]=compute_Zf_PS_ll(  argv, kinematic_2pt,  name, conf_jack, mass_jack_fit[ik2][ik1],  Njack ,plateaux_masses,outfile_f );
       //f_jack_fit[ik2][ik1]=compute_f_PS_GEVP(  argv, kinematic_2pt,  name, conf_jack,  mass_jack_fit_GEVP[ik2][ik1], Njack, var ,plateaux_f,outfile_f );
       
       free_corr(Neff, var*var, file_head.l0 ,data_bin);
       free_jack(Njack,var*var , file_head.l0, conf_jack);
      ////////////////////GEVP f_PS
  /*   
       for(t=1;t<file_head.l0/2;t++){    
           for (j=0;j<Njack;j++){
                for (i=0;i<var;i++){
                     projected_O[j][0][t][0]+=vec[t][j][i+0*var][0]*M[t][j][0+i*var][0];  //vec_0[i] M_0i
                }
               // printf("norm=%g\n",sqrt(vec[t][j][0][0]*vec[t][j][0][0]+vec[t][j][2][0]*vec[t][j][2][0] ));
           }
       }
       fprintf(outfile_f_GEVP,"#decay constant f_PS(t) of %s  propagators:1) mu %g r %d theta %g 2) mu %g r %d theta %g\n",name,
           file_head.k[ik2+file_head.nk],r_value(mr2),file_head.mom[mimom2][1],
           file_head.k[ik1+file_head.nk],r_value(r1), file_head.mom[imom1][1] );
       for(i=1;i<file_head.l0/2;i++){    
           for (j=0;j<Njack;j++){
              r[i][j]=matrix_element_GEVP(i,projected_O[j][0], mass_jack_fit_GEVP[ik2][ik1][j]);
              r[i][j]/=  (mass_jack_fit_GEVP[ik2][ik1][j]*mass_jack_fit_GEVP[ik2][ik1][j]  );
              r[i][j]*= (file_head.k[ik2+file_head.nk]+file_head.k[ik1+file_head.nk]);

            }
           mt[i]=mean_and_error_jack(Njack, r[i]);
           fprintf(outfile_f_GEVP,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);                
       }

       yn=1;
       if ( strcmp(argv[1],"see")==0 ){
            while(yn>0){
                printf("#f_PS^GEVP 1) mu %g r %d theta %g 2) mu %g r %d theta %g \n",
                file_head.k[ik2+file_head.nk],r_value(mr2),file_head.mom[mimom1][1],
                file_head.k[ik1+file_head.nk],r_value(r1),file_head.mom[imom2][1]);   
                plotting( file_head.l0, mt , &tmin,&tmax, &sep);
                m=try_linear_fit( tmin, tmax, sep , mt, r, Njack );
                yn=plotting_fit( file_head.l0, mt , tmin,tmax,m);
                free(m);
            }
       }
       if ( strcmp(argv[1],"read_plateaux")==0 ) fscanf(plateaux_f_GEVP,"%d  %d  %d\n",&tmin,&tmax,&sep);


       m=try_linear_fit( tmin,  tmax,sep , mt, r, Njack );
       //f_jack_fit[ik2]=give_jack_linear_fit( tmin,  tmax,sep , mt, r, Njack );

       fprintf(outfile_f_GEVP,"\n\n #f_PS fit in [%d,%d]\n  %.15g    %.15g\n\n\n",tmin,tmax,m[0],m[1]);
       printf("#f_PS GEVP(mu_h=%.4f, mu_l=%.4f) fit in [%d,%d]:  %.15g    %.15g\n",file_head.k[ik2+file_head.nk],file_head.k[ik1+file_head.nk] ,tmin,tmax,m[0],m[1]);
   
       if ( strcmp(argv[1],"see")==0 || strcmp(argv[1],"read_plateaux")==0 ){
           plotting_fit_pdf("f_{PS}^{GEVP}",file_head.l1, file_head.l0,file_head.beta,file_head.ksea,file_head.musea, mt , tmin,tmax,m,name, file_head.k[ik2+file_head.nk],r_value(mr2),
                file_head.k[ik1+file_head.nk],r_value(r1),file_head.mom[mimom1][1],file_head.mom[imom2][1] );
       }
    
       /////free memory
       for(i=1;i<file_head.l0/2;i++)
           free(mt[i]);
       free(m);
       fflush(outfile);
     ////////////////////end  f_PS
    */ 
     
     
     
/////end GEVP

 /*
   if(  strcmp(name,"V0P5")==0 ){
   
   for (i=0;i<confs;i++){
       extract_twopt(out[i] ,data[i][0],si_2pt,ii,imom2,imom1,ik2,r2,ik1,r1 );
       extract_twopt(out[i] ,data[i][1],si_2pt,1,imom2,imom1,ik2,r2,ik1,r1 );
   }
   
   antisymmetrise_corr(Neff, 0, file_head.l0,data);
   symmetric_derivative_corr(Neff, 0, file_head.l0,data);

   symmetrise_corr(Neff, 1, file_head.l0,data);

   conf_jack=create_jack(Neff, 2, file_head.l0, data);
  
      printf("HERE4\n"); 

////////////////////allocation
    
   fprintf(m_pcac,"#m_pcaca of %s  propagators:1) mu %g r %d theta %g 2) mu %g r %d theta %g\n",name,
           file_head.k[ik2+file_head.nk],r_value(mr2),file_head.mom[mimom2][1],
           file_head.k[ik1+file_head.nk],r_value(r1), file_head.mom[imom1][1] );
   for(i=1;i<file_head.l0/2;i++){    
        for (j=0;j<Njack;j++){
            r[i][j]=conf_jack[j][0][i][1]/conf_jack[j][1][i][0];
        }
        mt[i]=mean_and_error_jack(Njack, r[i]);
        fprintf(m_pcac,"%d   %.15e    %.15e\n",i,mt[i][0],mt[i][1]);
   }
    
   yn=1;
   if ( strcmp(argv[1],"see")==0 ){
      while(yn>0){
            printf("# m_pcac 1) mu %g r %d theta %g 2) mu %g r %d theta %g \n",
                file_head.k[ik2+file_head.nk],r_value(mr2),file_head.mom[mimom1][1],
                file_head.k[ik1+file_head.nk],r_value(r1),file_head.mom[imom2][1]);   
            plotting( file_head.l0, mt , &tmin,&tmax, &sep);
            m=try_linear_fit( tmin, tmax, sep , mt, r, Njack );
            yn=plotting_fit( file_head.l0, mt , tmin,tmax,m);
            free(m);
      }
   }

   m=try_linear_fit( tmin,  tmax,sep , mt, r, Njack );

   fprintf(m_pcac,"\n\n #m_pcac fit in [%d,%d]\n  %.15g    %.15g\n\n\n",tmin,tmax,m[0],m[1]);
   
   if ( strcmp(argv[1],"see")==0 ){
   plotting_fit_pdf("m_{PCAC}",file_head.l1, file_head.l0,file_head.beta,file_head.ksea,file_head.musea, mt , tmin,tmax,m,name, file_head.k[ik2+file_head.nk],r_value(mr2),
                file_head.k[ik1+file_head.nk],r_value(r1),file_head.mom[mimom1][1],file_head.mom[imom2][1] );
   }

   
    
   for(i=1;i<file_head.l0-1;i++)
        free(mt[i]);
   free(m);
    fflush(m_pcac);
   ///////free memory
//////////////////////////end m_pcac
   } 
    
  free_jack(Njack,4,file_head.l0,conf_jack); 
    
      printf("NEW mass\n"); 
 */   
}}   }}}} // end loop imom1 imom2 r1 ik1 r2 ik2  
   int kmin=4;
   contraction_index(&ii,"P5P5");
   contraction_name(ii,name);

  /* 
   double **fit3;
   if (file_head.nk>kmin){
        fit3=fit_polynomial( argv,"M_{PS}^{GEVP}" ,3,  kmin,0, name, mass_jack_fit_GEVP,  Njack , outfile );
        
        printf("fit=");
        for (i=0;i<=3;i++){
                printf("+(%g +- %g)x^ %d", fit3[i][0],fit3[i][1], i);
        }
        printf("\n");
   }*/
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////                        RATIOS                                                    //////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   for(ik1=0;ik1<ik1_max;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
   for(r1=0;r1<1;r1++){  // it is not a loop, r1=0
   for(ik2=ik1;ik2<file_head.nk;ik2++){
   for(r2=0;r2<1;r2++){    // it is not a loop, r2=0
   for(imom2=0;imom2<file_head.nmoms;imom2++){
   for(imom1=0;imom1<file_head.nmoms;imom1++){

       mr2=index_minus_r(r2);
       mr1=index_minus_r(r1);
       mimom1=index_minus_theta(imom1);
       mimom2=index_minus_theta(imom2);
//////////////////////////ratio of masses
      if (ik2==ik1)  fprintf(outfile_RM,"\n\n#Ratio mu_h  M(mu_h+1)/M(mu_h)   error  #  mu_l=%g\n",file_head.k[file_head.nk+ik1] );
      if (ik2-ik1-1>=0){ 
	      /* fprintf(outfile_RM,"#Ratio of m_eff(t) from %s :1) mu %g r %d theta %g 2) mu %g r %d theta %g /",name,
		   file_head.k[ik2+file_head.nk],r_value(mr2),file_head.mom[mimom2][1],
		   file_head.k[ik1+file_head.nk],r_value(r1), file_head.mom[imom1][1] );
	       fprintf(outfile_RM,"2) mu %g r %d theta %g 3) mu %g r %d theta %g \n",
		   file_head.k[ik2+file_head.nk-1],r_value(mr2),file_head.mom[mimom2][1],
		   file_head.k[ik1+file_head.nk],r_value(r1), file_head.mom[imom1][1] );*/
	       for (j=0;j<Njack;j++){
		   r[0][j]=mass_jack_fit[ik2][ik1][j]/mass_jack_fit[ik2-1][ik1][j];
	       }
	       m=mean_and_error(argv[4],Njack, r[0]);
	       fprintf(outfile_RM,"%.4f   %.15e    %.15e\n",file_head.k[ik2+file_head.nk-1]  ,m[0],m[1]);
	    
	       free(m);   
	       fflush(outfile_RM);
      }

////////////////////end  
//////////////////////////ratio of masses  GEVP
      if (ik2==ik1)  fprintf(outfile_RM_GEVP,"\n\n#Ratio mu_h  M(mu_h+1)/M(mu_h)    error  #  mu_l=%g\n",file_head.k[file_head.nk+ik1] );
      if (ik2-ik1-1>=0){ 
	       for (j=0;j<Njack;j++){
		   r[0][j]=mass_jack_fit_GEVP[ik2][ik1][j]/mass_jack_fit_GEVP[ik2-1][ik1][j];
	       }
	       m=mean_and_error(argv[4],Njack, r[0]);
	       fprintf(outfile_RM_GEVP,"%.4f   %.15e    %.15e\n",file_head.k[ik2+file_head.nk-1]  ,m[0],m[1]);
	    
	       free(m);   
	       fflush(outfile_RM_GEVP);
      }
////////////////////end  
//////////////////////////ratio of f_PS
      if (ik2==ik1)  fprintf(outfile_Rf,"\n\n#Ratio mu_h  f_{PS}(mu_h+1)/f_{PS}M(mu_h)    error  #  mu_l=%g\n",file_head.k[file_head.nk+ik1] );
      if (ik2-ik1-1>=0){ 
    	      /* fprintf(outfile_RM,"#Ratio of m_eff(t) from %s :1) mu %g r %d theta %g 2) mu %g r %d theta %g /",name,
		   file_head.k[ik2+file_head.nk],r_value(mr2),file_head.mom[mimom2][1],
		   file_head.k[ik1+file_head.nk],r_value(r1), file_head.mom[imom1][1] );
	       fprintf(outfile_RM,"2) mu %g r %d theta %g 3) mu %g r %d theta %g \n",
		   file_head.k[ik2+file_head.nk-1],r_value(mr2),file_head.mom[mimom2][1],
		   file_head.k[ik1+file_head.nk],r_value(r1), file_head.mom[imom1][1] );*/
	       for (j=0;j<Njack;j++){
               r[0][j]=f_jack_fit[ik2][ik1][j]/f_jack_fit[ik2-1][ik1][j];
	       }
	       m=mean_and_error(argv[4],Njack, r[0]);
	       fprintf(outfile_Rf,"%.4f   %.15e    %.15e\n",file_head.k[ik2+file_head.nk-1]  ,m[0],m[1]);
	    
	       free(m);   
	       fflush(outfile_Rf);
      }
////////////////////end  P5P5 mass

}}   }}}} // end loop imom1 imom2 r1 ik1 r2 ik2  


printf("here\n");




   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r); 
printf("here\n");
   
   fclose(outfile);
   fclose(m_pcac);
   fclose(f_ll);
  // free_corr(confs, var*var, file_head.l0 ,data);
printf("here\n");

  // for (i=0;i<confs;i++)
  //     free(out[i]);
free(out);
 printf("here\n");

    return 0;   
}
