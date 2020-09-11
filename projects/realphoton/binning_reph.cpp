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
#include "resampling.hpp"
#include "read.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "indices.hpp"
#include "mutils.hpp"


#include <unistd.h>
#include <omp.h>

int Nboot=100;
int fdA,fdV;
struct  kinematic kinematic_2pt;
struct  kinematic_G kinematic_2pt_G;
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
int r_value(int r)
{   int a;
    if (r==0) a=1;
    if (r==1) a=-1;
    return a;
}
void get_kinematic( int ik2, int ik1,int imom2, int imom1 ){
    kinematic_2pt.k2=file_head.k[ik2+file_head.nk];
    kinematic_2pt.k1=file_head.k[ik1+file_head.nk];
    
    kinematic_2pt.mom2=-file_head.mom[imom2][3];
    if (kinematic_2pt.mom2==0) kinematic_2pt.mom2=0;
    kinematic_2pt.mom1=file_head.mom[imom1][3];

    kinematic_2pt.mom02=file_head.mom[imom2][0];
    kinematic_2pt.mom01=file_head.mom[imom1][0];
    
    kinematic_2pt.r2=-1;
    kinematic_2pt.r1=1;
    
    
    int i ;
    kinematic_2pt.Mom2[0]=file_head.mom[imom2][0];
    kinematic_2pt.Mom1[0]=file_head.mom[imom1][0];
 
    for (i=1;i<4;i++){
        kinematic_2pt.Mom2[i]=-file_head.mom[imom2][i];
        if (kinematic_2pt.Mom2[i]==0) kinematic_2pt.Mom2[i]=0;
        kinematic_2pt.Mom1[i]=file_head.mom[imom1][i];
 
    }
        
 
}

void get_kinematic_G( int ikt, int iks,int imom0, int imomt, int imoms ){
    
    double Pi=3.141592653589793;
    double k[4],p[4],k3;
    double L[4];
    int i;
  
    L[0]=file_head.l0; L[1]=file_head.l1; L[2]=file_head.l2; L[3]=file_head.l3;         

    kinematic_2pt_G.i=index_n_twopt_G_fit(ikt,iks,imom0,imomt,imoms);
    kinematic_2pt_G.kt=file_head.k[ikt+file_head.nk];
    kinematic_2pt_G.ks=file_head.k[iks+file_head.nk];
    
    
    kinematic_2pt_G.rt=1;
    kinematic_2pt_G.rs=-1;
    
    kinematic_2pt_G.Mom0[0]=file_head.mom[imom0][0];
    kinematic_2pt_G.Momt[0]=file_head.mom[imomt][0];
    kinematic_2pt_G.Moms[0]=file_head.mom[imoms][0];
    
    for (i=1;i<4;i++){
        kinematic_2pt_G.Mom0[i]=file_head.mom[imom0][i];
        if (kinematic_2pt_G.Mom0[i]==0) kinematic_2pt_G.Mom0[i]=0;
        kinematic_2pt_G.Momt[i]=file_head.mom[imomt][i];
        if (kinematic_2pt_G.Momt[i]==0) kinematic_2pt_G.Momt[i]=0;
        kinematic_2pt_G.Moms[i]=-file_head.mom[imoms][i];
        if (kinematic_2pt_G.Moms[i]==0) kinematic_2pt_G.Moms[i]=0;
    }
    
    
    for(i=0;i<4;i++){

          k[i]=(2*Pi/L[i])*(kinematic_2pt_G.Mom0[i]-kinematic_2pt_G.Momt[i]);
          k[i]=2*sin(k[i]/2.);
          p[i]=(2*Pi/L[i])*(kinematic_2pt_G.Mom0[i]+kinematic_2pt_G.Moms[i]);
          p[i]=2*sin(p[i]/2.);
	 // p[i]=fabs(p[i]);
    }
    error(fabs(k[0])>0.00001,1,"m_eff.c","error: k0 of the photon is not zero");    
    error(fabs(p[0])>0.00001,1,"m_eff.c","error: p0 of the photon is not zero");    
      
    kinematic_2pt_G.E_g=2.*asinh(   sqrt(k[1]*k[1]+k[2]*k[2]+k[3]*k[3])/2.   );
    kinematic_2pt_G.E_gT=sinh(kinematic_2pt_G.E_g)*(1.-exp(-L[0]*kinematic_2pt_G.E_g));
    
    k3=kinematic_2pt_G.E_g*k[3]/fabs(k[3]);
    kinematic_2pt_G.kp=0;
     for(i=1;i<4;i++){
          kinematic_2pt_G.kp+= k[i]*p[i];
    }
    //kinematic_2pt_G.kp= kinematic_2pt_G.E_g*p[3];
    //kinematic_2pt_G.kp= k3*p[3];

    
    kinematic_2pt_G.eps1_curl_p[0]=0;
    kinematic_2pt_G.eps1_curl_p[1]=kinematic_2pt_G.eps1[2]*p[3]-kinematic_2pt_G.eps1[3]*p[2];
    kinematic_2pt_G.eps1_curl_p[2]=-kinematic_2pt_G.eps1[1]*p[3]+kinematic_2pt_G.eps1[3]*p[1];
    kinematic_2pt_G.eps1_curl_p[3]=kinematic_2pt_G.eps1[1]*p[2]-kinematic_2pt_G.eps1[2]*p[1];
    
    kinematic_2pt_G.eps2_curl_p[0]=0;
    kinematic_2pt_G.eps2_curl_p[1]=kinematic_2pt_G.eps2[2]*p[3]-kinematic_2pt_G.eps2[3]*p[2];
    kinematic_2pt_G.eps2_curl_p[2]=-kinematic_2pt_G.eps2[1]*p[3]+kinematic_2pt_G.eps2[3]*p[1];
    kinematic_2pt_G.eps2_curl_p[3]=kinematic_2pt_G.eps2[1]*p[2]-kinematic_2pt_G.eps2[2]*p[1];
    
    
    kinematic_2pt_G.eps1_curl_k[0]=0;
    kinematic_2pt_G.eps1_curl_k[1]=kinematic_2pt_G.eps1[2]*k[3]-kinematic_2pt_G.eps1[3]*k[2];
    kinematic_2pt_G.eps1_curl_k[2]=-kinematic_2pt_G.eps1[1]*k[3]+kinematic_2pt_G.eps1[3]*k[1];
    kinematic_2pt_G.eps1_curl_k[3]=kinematic_2pt_G.eps1[1]*k[2]-kinematic_2pt_G.eps1[2]*k[1];
    
    kinematic_2pt_G.eps2_curl_k[0]=0;
    kinematic_2pt_G.eps2_curl_k[1]=kinematic_2pt_G.eps2[2]*k[3]-kinematic_2pt_G.eps2[3]*k[2];
    kinematic_2pt_G.eps2_curl_k[2]=-kinematic_2pt_G.eps2[1]*k[3]+kinematic_2pt_G.eps2[3]*k[1];
    kinematic_2pt_G.eps2_curl_k[3]=kinematic_2pt_G.eps2[1]*k[2]-kinematic_2pt_G.eps2[2]*k[1];
    
 /*   kinematic_2pt_G.eps1_curl_k[0]=0;
    kinematic_2pt_G.eps1_curl_k[1]=kinematic_2pt_G.eps1[2]*k3;
    kinematic_2pt_G.eps1_curl_k[2]=-kinematic_2pt_G.eps1[1]*k3;
    kinematic_2pt_G.eps1_curl_k[3]=0;
    
    kinematic_2pt_G.eps2_curl_k[0]=0;
    kinematic_2pt_G.eps2_curl_k[1]=kinematic_2pt_G.eps2[2]*k3;
    kinematic_2pt_G.eps2_curl_k[2]=-kinematic_2pt_G.eps2[1]*k3;
    kinematic_2pt_G.eps2_curl_k[3]=0;
    */
}


double timestamp()
{
	struct timeval tm;
	gettimeofday(&tm, NULL);
	return tm.tv_sec + 1.0e-6 * tm.tv_usec;
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
    for(i=file_head.nk;i<2*file_head.nk;++i)
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
   tmp-= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*10 ;
   //tmp-= sizeof(int)*2;
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



int index_n_minus_r( int r)
{
    int a;
    if (r==0) a=1;
    if (r==1) a=0;
    return a;
}



/*
void read_chunk(FILE *stream,int index, double *obs, double **to_write,int vol, int si){
       int t;  
       fseek(stream, index, SEEK_CUR);
       fread(obs,sizeof(double),2*si*file_head.l0,stream); 
   
       for(t=0;t<file_head.l0;t++){
           to_write[t][0]+=obs[si*t*2];
           to_write[t][1]+=obs[si*t*2+1];
       }
       vol++;

}*/
void read_chunk(FILE *stream,int tmp,int index, double *obs, double **to_write,int *vol, int si){
       int t;  
       fseek(stream,tmp+ index, SEEK_SET);
       fread(obs,sizeof(double),2*si*file_head.l0,stream); 
   
       for(t=0;t<file_head.l0;t++){
           to_write[t][0]+=obs[si*t*2];
           to_write[t][1]+=obs[si*t*2+1];
       }
       *vol+=1;

}
/*
void read_twopt_Nconfs(FILE *stream,int size, int Nconfs,int var , double ****to_write,int si, int ii,int ik1, int ik2, int imom1, int imom2 ){
  long int tmp;
   int iiconf,N,s;
   double *obs;
   int cik1, cik2;
   int cimom1,cimom2,cr1;
   int t,vol=0,index,index1,i;
   double re,im;
   int    vol3=file_head.l1*file_head.l2*file_head.l3; 
   int diff;

   tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
   fseek(stream, tmp, SEEK_SET);
  // tmp+=sizeof(double)*iconf*size+sizeof(int)*iconf;

   obs=(double*) malloc(2*si*file_head.l0*sizeof(double)); 
   
   cimom1=index_n_minus_theta(imom1);
   cimom2=index_n_minus_theta(imom2);
   cik1=index_n_minus_kappa(ik1);
   cik2=index_n_minus_kappa(ik2);
   
   tmp=sizeof(double)*size+sizeof(int);
   for (i=0;i<Nconfs;i++){
   
   index=sizeof(double)*2*index_n_twopt(si,ii,0,ik1,ik2,imom1,imom2);
   fseek(stream, index, SEEK_CUR);

   fread(obs,sizeof(double),2*si*file_head.l0,stream); 
   
   for(t=0;t<file_head.l0;t++){
       to_write[i][var][t][0]=obs[si*t*2];
       to_write[i][var][t][1]=obs[si*t*2+1];
   }
   vol++;
   
   index1=sizeof(double)*2*index_n_twopt(si,ii,0,ik1,ik2,imom2,imom1);
   diff=index1-index;
   index=index1;
   read_chunk(stream,  diff,  obs,  to_write[i][var], vol,  si);
  
   index1=sizeof(double)*2*index_n_twopt(si,ii,0,ik2,ik1,imom1,imom2);
   diff=index1-index;
   index=index1;
   read_chunk(stream,  diff,  obs,  to_write[i][var], vol,  si);

   if (  cimom1>=0 &&  cimom2 >=0 ){ 
       index1=sizeof(double)*2*index_n_twopt(si,ii,0,ik1,ik2,cimom1,cimom2);
       diff=index1-index;
       index=index1;
       read_chunk(stream,  diff,  obs,  to_write[i][var], vol,  si);

       index1=sizeof(double)*2*index_n_twopt(si,ii,0,ik1,ik2,cimom2,cimom1);
       diff=index1-index;
       index=index1;
       read_chunk(stream,  diff,  obs,  to_write[i][var], vol,  si);

       index1=sizeof(double)*2*index_n_twopt(si,ii,0,ik2,ik1,cimom1,cimom2);
       diff=index1-index;
       index=index1;
       read_chunk(stream,  diff,  obs,  to_write[i][var], vol,  si);


   }
   diff=tmp-index;
   fseek(stream, diff, SEEK_CUR);
   
   for(t=0;t<file_head.l0;t++){
  	   to_write[i][var][t][0]/=( (double) vol*vol3 );
	   to_write[i][var][t][1]/=( (double) vol*vol3 );
   }

   }   
   free(obs);

    
}*/
void read_twopt(FILE *stream,int size, int iconf , double **to_write,int si, int ii,int ik1, int ik2, int imom1, int imom2 ){
   
   long int tmp;
   int iiconf,N,s;
   double *obs;
   int cik1, cik2;
   int cimom1,cimom2,cr1;
   int t,vol=0,index;
   double re,im;
   int    vol3=file_head.l1*file_head.l2*file_head.l3; 

   tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
   tmp+=sizeof(double)*iconf*size+sizeof(int)*iconf;

   obs=(double*) malloc(2*si*file_head.l0*sizeof(double)); 
   
   cimom1=index_n_minus_theta(imom1);
   cimom2=index_n_minus_theta(imom2);
   cik1=index_n_minus_kappa(ik1);
   cik2=index_n_minus_kappa(ik2);
   
   index=sizeof(double)*2*index_n_twopt(si,ii,0,ik1,ik2,imom1,imom2);
   fseek(stream, tmp+index, SEEK_SET);

   fread(obs,sizeof(double),2*si*file_head.l0,stream); 
   
   for(t=0;t<file_head.l0;t++){
       to_write[t][0]=obs[si*t*2];
       to_write[t][1]=obs[si*t*2+1];
   }
   vol++;
   
   index=sizeof(double)*2*index_n_twopt(si,ii,0,ik1,ik2,imom2,imom1);
   read_chunk(stream, tmp, index,  obs,  to_write, &vol,  si);
  
   index=sizeof(double)*2*index_n_twopt(si,ii,0,ik2,ik1,imom1,imom2);
   read_chunk(stream, tmp, index,  obs,  to_write, &vol,  si);

   if (  cimom1>=0 &&  cimom2 >=0 ){ 
       index=sizeof(double)*2*index_n_twopt(si,ii,0,ik1,ik2,cimom1,cimom2);
       read_chunk(stream, tmp, index,  obs,  to_write, &vol,  si);

       index=sizeof(double)*2*index_n_twopt(si,ii,0,ik1,ik2,cimom2,cimom1);
       read_chunk(stream, tmp, index,  obs,  to_write, &vol,  si);

       index=sizeof(double)*2*index_n_twopt(si,ii,0,ik2,ik1,cimom1,cimom2);
       read_chunk(stream, tmp, index,  obs,  to_write, &vol,  si);


   }
   
   for(t=0;t<file_head.l0;t++){
  	   to_write[t][0]/=( (double) vol*vol3 );
	   to_write[t][1]/=( (double) vol*vol3 );
   }
   
   


   free(obs);

}

/*

void read_twopt(FILE *stream,int size, int iconf , double **to_write,int si, int ii,int ik1, int ik2, int imom1, int imom2 ){
   
   long int tmp;
   int iiconf,N,s;
   double *obs;
   int cik1, cik2;
   int cimom1,cimom2,cr1;
   int t,vol,index;
   double re,im;
   int    vol3=file_head.l1*file_head.l2*file_head.l3; 

   tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
   tmp+=sizeof(double)*iconf*size+sizeof(int)*iconf;

   obs=(double*) malloc(2*sizeof(double)); 
   
   cimom1=index_n_minus_theta(imom1);
   cimom2=index_n_minus_theta(imom2);
   cik1=index_n_minus_kappa(ik1);
   cik2=index_n_minus_kappa(ik2);

   
	for(t=0;t<file_head.l0;t++){
	   re=0;im=0;vol=0;
	   index=sizeof(double)*2*index_n_twopt(si,ii,t,ik1,ik2,imom1,imom2);
       fseek(stream, tmp+index, SEEK_SET);
       fread(obs,sizeof(double),2,stream); 
       re+= obs[0];
       im+= obs[1];
       vol++;
       index=sizeof(double)*2*index_n_twopt(si,ii,t,ik1,ik2,imom2,imom1);
       fseek(stream, tmp+index, SEEK_SET);
       fread(obs,sizeof(double),2,stream); 
       re+= obs[0];
       im+= obs[1];
       vol++;
       index=sizeof(double)*2*index_n_twopt(si,ii,t,ik2,ik1,imom1,imom2);
       fseek(stream, tmp+index, SEEK_SET);
       fread(obs,sizeof(double),2,stream); 
       re+= obs[0];
       im+= obs[1];
       vol++;
	   if (  cimom1>=0 &&  cimom2 >=0 )
  	   { 
  	       index=sizeof(double)*2*index_n_twopt(si,ii,t,ik1,ik2,cimom1,cimom2);
           fseek(stream, tmp+index, SEEK_SET);
           fread(obs,sizeof(double),2,stream); 
           re+= obs[0];
           im+= obs[1];
           vol++;
  	       index=sizeof(double)*2*index_n_twopt(si,ii,t,ik1,ik2,cimom2,cimom1);
           fseek(stream, tmp+index, SEEK_SET);
           fread(obs,sizeof(double),2,stream); 
           re+= obs[0];
           im+= obs[1];
           vol++;
  	       index=sizeof(double)*2*index_n_twopt(si,ii,t,ik2,ik1,cimom1,cimom2);
           fseek(stream, tmp+index, SEEK_SET);
           fread(obs,sizeof(double),2,stream); 
           re+= obs[0];
           im+= obs[1];
           vol++;
             	
  	   }
	   to_write[t][0]=re/( (double) vol*vol3 );
	   to_write[t][1]=im/( (double) vol*vol3 );
	}


   free(obs);

}*/
void read_twopt_gamma(FILE *stream,int size, int iconf , double **to_write, const char*name,int si,int ikt,int iks,int imom0,int imomt,int imoms ){
   
   long int tmp ;
   int iiconf,N,s;
   double *obs;
   int cikt, ciks;
   int cimom0,cimomt,cimoms;
   int ix0,vol,index;
   int    vol3=file_head.l1*file_head.l2*file_head.l3; 
   double re,im;
   double y0=2./3,ys=-1./3.;
  
   
   vol3=file_head.l1*file_head.l2*file_head.l3; 
   tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
   tmp+=sizeof(double)*iconf*size+sizeof(int)*iconf;
   obs=(double*) malloc(si*file_head.l0*2*sizeof(double)); 


   cimom0=index_n_minus_theta(imom0);
   cimomt=index_n_minus_theta(imomt);
   cimoms=index_n_minus_theta(imoms);
   cikt=index_n_minus_kappa(ikt);
   ciks=index_n_minus_kappa(iks);
  
   
   index=sizeof(double)*2*index_n_twoptgamma( si, 0, 0, ikt, iks, imom0, imomt, imoms);
   fseek(stream, tmp+index, SEEK_SET);
   fread(obs,sizeof(double),si*file_head.l0*2,stream); 
   
    
   for(ix0=0;ix0<file_head.l0;ix0++){
	   re=0;vol=0;im=0;
       if (strcmp(name,"oAmuGPo")==0){
            re+= obs[(1+si*ix0)*2];
            im+= obs[(1+si*ix0)*2+1];
            vol++;
            
            re+= obs[(2+si*ix0)*2];
            im+= obs[(2+si*ix0)*2+1];
            vol++;
           
            re-= obs[(5+si*ix0)*2];
            im-= obs[(5+si*ix0)*2+1];
            vol++;
           
            re+= obs[(6+si*ix0)*2];
            im+= obs[(6+si*ix0)*2+1];
            vol++;
       }
       if (strcmp(name,"oVmuGPo")==0){
           re+= obs[(1+si*ix0)*2];
            im+= obs[(1+si*ix0)*2+1];
            vol++;
            
            re-= obs[(2+si*ix0)*2];
            im-= obs[(2+si*ix0)*2+1];
            vol++;
            
            re+= obs[(5+si*ix0)*2];
            im+= obs[(5+si*ix0)*2+1];
            vol++;
            
            re+= obs[(6+si*ix0)*2];
            im+= obs[(6+si*ix0)*2+1];
            vol++;
       }
 
	   to_write[ix0][0]=re/( (double) vol*vol3 );
	   to_write[ix0][1]=im/( (double) vol*vol3 );
   }
   index=sizeof(double)*2*index_n_twoptgamma( si, 0, 0, iks, ikt, imom0, imomt, imoms);
   fseek(stream, tmp+index, SEEK_SET);
   fread(obs,sizeof(double),si*file_head.l0*2,stream); 
    
   for(ix0=0;ix0<file_head.l0;ix0++){
	   re=0;vol=0;im=0;
       if (strcmp(name,"oAmuGPo")==0){
            re+= obs[(1+si*ix0)*2];
            im+= obs[(1+si*ix0)*2+1];
            vol++;
            
            re+= obs[(2+si*ix0)*2];
            im+= obs[(2+si*ix0)*2+1];
            vol++;
           
            re-= obs[(5+si*ix0)*2];
            im-= obs[(5+si*ix0)*2+1];
            vol++;
           
            re+= obs[(6+si*ix0)*2];
            im+= obs[(6+si*ix0)*2+1];
            vol++;

	    re=-re;im=-im;
       }
       if (strcmp(name,"oVmuGPo")==0){
            re+= obs[(1+si*ix0)*2];
            im+= obs[(1+si*ix0)*2+1];
            vol++;
            
            re-= obs[(2+si*ix0)*2];
            im-= obs[(2+si*ix0)*2+1];
            vol++;
            
            re+= obs[(5+si*ix0)*2];
            im+= obs[(5+si*ix0)*2+1];
            vol++;
            
            re+= obs[(6+si*ix0)*2];
            im+= obs[(6+si*ix0)*2+1];
            vol++;
       }
 
	   to_write[ix0][0]= y0*to_write[ix0][0] +  ys*re/( (double) vol*vol3 );
	   to_write[ix0][1]= y0*to_write[ix0][1] +  ys*im/( (double) vol*vol3 );
   }
   
   free(obs);
}


void read_twoptgamma_Nconfs(FILE *stream,int size, int Nconfs,int var , double ****to_write,int si, int ii,int ikt, int iks, int imom0,int imomt, int imoms ){
  long int tmp;
   int iiconf,N,s;
   double *obs;
   int cik1, cik2;
   int cimom1,cimom2,cr1;
   int t,vol=0,index,index1,i;
   double re,im;
   int    vol3=file_head.l1*file_head.l2*file_head.l3; 
   int diff;

   tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
   fseek(stream, tmp, SEEK_SET);
  // tmp+=sizeof(double)*iconf*size+sizeof(int)*iconf;

   obs=(double*) malloc(2*si*file_head.l0*sizeof(double)); 
   
   cimom1=index_n_minus_theta(imom0);
   cimom2=index_n_minus_theta(imomt);
   cik1=index_n_minus_kappa(ikt);
   cik2=index_n_minus_kappa(iks);
   
   tmp=sizeof(double)*size+sizeof(int);
   for (i=0;i<Nconfs;i++){
   
   index=sizeof(double)*2*index_n_twoptgamma(si,ii,0,ikt,iks,imom0,imomt,imoms);
   fseek(stream, index, SEEK_CUR);

   fread(obs,sizeof(double),2*si*file_head.l0,stream); 
   
   for(t=0;t<file_head.l0;t++){
       to_write[i][var][t][0]=obs[si*t*2];
       to_write[i][var][t][1]=obs[si*t*2+1];
   }
   vol++;
   
  
   diff=tmp-index;
   fseek(stream, diff, SEEK_CUR);
   
   for(t=0;t<file_head.l0;t++){
  	   to_write[i][var][t][0]/=( (double) vol*vol3 );
	   to_write[i][var][t][1]/=( (double) vol*vol3 );
   }

   }   
   free(obs);

    
}
 
void setup_single_file_jack(char  *save_name,char **argv, const char  *name,int Njack){
     FILE *f;
     mysprintf(save_name,NAMESIZE,"%s/%s",argv[3],name);
     f=fopen(save_name,"w+");
     error(f==NULL,1,"setup_file_jack ",
         "Unable to open output jackknife file %s/%s",argv[3],name);
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
    
     setup_single_file_jack(file_jack.M_PS_GEVP,argv,"jackknife/M_{PS}^{GEVP}_jack",Njack);
     setup_single_file_jack(file_jack.M_PS,argv,"jackknife/M_{PS}_jack",Njack);
     setup_single_file_jack(file_jack.f_PS,argv,"jackknife/f_{PS}_jack",Njack);
     setup_single_file_jack(file_jack.FV,argv,"jackknife/FV_jack",Njack);
     setup_single_file_jack(file_jack.FAp1,argv,"jackknife/FAp1_jack",Njack);
     setup_single_file_jack(file_jack.FV_exclude,argv,"jackknife/FV_jack",Njack);
     setup_single_file_jack(file_jack.FAp1_exclude,argv,"jackknife/FAp1_jack",Njack);
     setup_single_file_jack(file_jack.kp,argv,"jackknife/FAp1_jack",Njack);

}

void allocate_jack_fit(struct jack_fit jack_fit, int sizePP){
   jack_fit.m=      (double**) malloc(sizeof(double*)*sizePP/file_head.l0);
   jack_fit.oPp=      (double**) malloc(sizeof(double*)*sizePP/file_head.l0);
   jack_fit.f_PS= (double**) malloc(sizeof(double*)*sizePP/file_head.l0);
   jack_fit.Zf_PS= (double**) malloc(sizeof(double*)*sizePP/file_head.l0);

   jack_fit.RA= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   jack_fit.HA= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   
   jack_fit.RV= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   jack_fit.HV= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);

}

int main(int argc, char **argv){
   int sizePP, sizeAmuP,sizeAmuGP,sizeVV,sizeVmuGP;
   int i,j,t,iG;
   clock_t t1,t2;
   double dt;
   int *sym;

   int *iconf,confs;
   double ****data,****data_bin, **out,**tmp; 
   char c;
   double *in;

   double ****M,****vec,****projected_O;
   double  ****lambda,****lambda0;
   
   double *fit,***y,*x,*m,*me;
   

   double ****conf_jack,**r,**mt,**met;
   int Ncorr=1;
   int t0=2;
   
   FILE  *oPPo=NULL, *oAmuPo=NULL,*oAmuGPo=NULL,*oVVo=NULL, *oVmuGPo=NULL;
   FILE  *oPPo_b=NULL, *oAmuPo_b=NULL,*oAmuGPo_b=NULL,*oVVo_b=NULL, *oVmuGPo_b=NULL;
   
   FILE *plateaux_masses=NULL, *plateaux_masses_GEVP=NULL; 
   FILE *plateaux_f=NULL;   
FILE  *plateaux_RA=NULL,  *plateaux_RV=NULL;
   struct jack_fit jack_fit;
   
   error(argc!=6,1,"main ",
         "usage argc!=5: ./form_factors blind/see/read_plateaux -p inpath   jack/boot  pdf");

   error(strcmp(argv[2],"-p")!=0,1,"main ",
         "missing -p \n usage: ./form_factors blind/see/read_plateaux -p inpath   jack/boot");

    error(strcmp(argv[4],"jack")!=0 && strcmp(argv[4],"boot")!=0 ,2,"main ",
         "choose jack or boot \n usage: ./form_factors blind/see/read_plateaux -p inpath   jack/boot");
    
    
   char namefile[NAMESIZE];
   
   
   
  
   double E_B,E_Pi, x_SCHET,q2,vec_pB,vec_pPi;
   int Neff,Njack,bin=10;
   
   t1=clock();
   
   // f=fopen("./meas_2pts_bin10.dat","r");
   mysprintf(namefile,NAMESIZE,"%s/data/oPPo-ss_conf.realph.dat",argv[3]);
   oPPo=fopen(namefile,"r"); error(oPPo==NULL,1, "main","2pt file %s not found",namefile );
   mysprintf(namefile,NAMESIZE,"%s/data/oAmuPo-ss_conf.realph.dat",argv[3]);
   oAmuPo=fopen(namefile,"r"); error(oAmuPo==NULL,1, "main","2pt file %s not found",namefile ); 
   mysprintf(namefile,NAMESIZE,"%s/data/oAmuGPo-gs_conf.realph.dat",argv[3]);
   oAmuGPo=fopen(namefile,"r"); error(oAmuGPo==NULL,1, "main","2pt file %s not found",namefile );   
   std::ifstream oAmuGPopp (namefile, std::ifstream::binary);  
   //mysprintf(namefile,"%s/data/oVVo-ss_conf.realph.dat",argv[3]);
  // oVVo=fopen(namefile,"r"); error(oVVo==NULL,1, "main","2pt file %s not found",namefile );
   mysprintf(namefile,NAMESIZE,"%s/data/oVmuGPo-gs_conf.realph.dat",argv[3]);
   oVmuGPo=fopen(namefile,"r"); error(oVmuGPo==NULL,1, "main","2pt file %s not found",namefile );

/////////////////////////////////////////////////////////////binned data   
   mysprintf(namefile,NAMESIZE,"%s/data/oPPo-ss_bin%d_conf.realph.dat",argv[3],bin);
   oPPo_b=fopen(namefile,"w+"); error(oPPo==NULL,1, "main","2pt file %s not found",namefile );
   mysprintf(namefile,NAMESIZE,"%s/data/oAmuPo-ss_bin%d_conf.realph.dat",argv[3],bin);
   oAmuPo_b=fopen(namefile,"w+"); error(oAmuPo==NULL,1, "main","2pt file %s not found",namefile ); 
   mysprintf(namefile,NAMESIZE,"%s/data/oAmuGPo-gs_bin%d_conf.realph.dat",argv[3],bin);
   oAmuGPo_b=fopen(namefile,"w+"); error(oAmuGPo==NULL,1, "main","2pt file %s not found",namefile );   
   //mysprintf(namefile,"%s/data/oVVo-ss_conf.realph.dat",argv[3]);
  // oVVo=fopen(namefile,"r"); error(oVVo==NULL,1, "main","2pt file %s not found",namefile );
   mysprintf(namefile,NAMESIZE,"%s/data/oVmuGPo-gs_bin%d_conf.realph.dat",argv[3],bin);
   oVmuGPo_b=fopen(namefile,"w+"); error(oVmuGPo==NULL,1, "main","2pt file %s not found",namefile );

   
   read_file_head_bin(oPPo);
   read_nconfs(&sizePP,&confs,oPPo);   
   printf("oPPo\t size=%d  confs=%d  \n",sizePP,confs);
 
   read_file_head_bin(oAmuPo);
   read_nconfs(&sizeAmuP,&confs,oAmuPo);   
   printf("oAmuPo\t size=%d  confs=%d \n",sizeAmuP,confs);
 
   read_file_head_bin(oAmuGPo);
   read_nconfs(&sizeAmuGP,&confs,oAmuGPo);   
   printf("oAmuGPo\t size=%d  confs=%d \n",sizeAmuGP,confs);
   
   /*  read_file_head_bin(oVVo);
   read_nconfs(&sizeVV,&confs,oVVo);   
   printf("oVVo\t size=%d  confs=%d \n",sizeVV,confs);
*/
   read_file_head_bin(oVmuGPo);
   read_nconfs(&sizeVmuGP,&confs,oVmuGPo);   
   printf("oVmuGPo\t size=%d  confs=%d \n",sizeVmuGP,confs);

 ////////////////////////////////////////////writing file
   int conf_out=confs/bin;
   printf("HERE\n");
   write_file_head(oPPo_b);
   fwrite(&sizePP,sizeof(int),1,oPPo_b);

   write_file_head(oAmuPo_b);
   fwrite(&sizeAmuP,sizeof(int),1,oAmuPo_b);
 
    write_file_head(oAmuGPo_b);
   fwrite(&sizeAmuGP,sizeof(int),1,oAmuGPo_b);
   
 /*   write_file_head(oVVo_b);
   fwrite(&sizeVV,sizeof(int),1,oVVo_b);*/
 
   write_file_head(oVmuGPo_b);
   fwrite(&sizeVmuGP,sizeof(int),1,oVmuGPo_b);
  
   
   
   double *obs_PP=(double*) malloc(sizeof(double)*(sizePP));
   double *sum_PP=(double*) calloc(sizePP,sizeof(double));
   int  iii;
   int ib=0;
   
   ib=0;
   for (i=1;i<=confs;i++){   
         /*tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
         tmp+=sizeof(double)*iconf*sizePP+sizeof(int)*iconf;
         fseek(oPPo, tmp+index, SEEK_SET);*/
         fread(&iii,sizeof(int),1,oPPo);
         fread(obs_PP,sizeof(double),sizePP,oPPo);
         for (j=0;j<sizePP;j++){
             sum_PP[j]+=obs_PP[j];
         }
         if ((i%(bin))==0){
             for (j=0;j<sizePP;j++)
                sum_PP[j]/=(double) bin;
             fwrite(&ib,sizeof(int),1,oPPo_b);
             fwrite(sum_PP,sizeof(double),sizePP,oPPo_b);
              for (j=0;j<sizePP;j++)
                sum_PP[j]=0;
             ib++;
        }

   }
   fclose(oPPo_b);
   printf("PP binned\n");
   
   double *obs_AmuP=(double*) malloc(sizeof(double)*(sizeAmuP));
   double *sum_AmuP=(double*) calloc(sizeAmuP,sizeof(double));
   
   
   ib=0;
   for (i=1;i<=confs;i++){   
         /*tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
         tmp+=sizeof(double)*iconf*sizeAmuP+sizeof(int)*iconf;
         fseek(oAmuPo, tmp+index, SEEK_SET);*/
         fread(&iii,sizeof(int),1,oAmuPo);
         fread(obs_AmuP,sizeof(double),sizeAmuP,oAmuPo); 
         for (j=0;j<sizeAmuP;j++){
             sum_AmuP[j]+=obs_AmuP[j];
         }
         if ((i%(bin))==0){
             for (j=0;j<sizeAmuP;j++)
                sum_AmuP[j]/=(double) bin;
             fwrite(&ib,sizeof(int),1,oAmuPo_b);
             fwrite(sum_AmuP,sizeof(double),sizeAmuP,oAmuPo_b);
              for (j=0;j<sizeAmuP;j++)
                sum_AmuP[j]=0;
             ib++;
        }

   }
   fclose(oAmuPo_b);
    printf("AmuP binned\n");
   double *obs_AmuGP=(double*) malloc(sizeof(double)*(sizeAmuGP));
   double *sum_AmuGP=(double*) calloc((sizeAmuGP),sizeof(double));
   
   
   ib=0;
   for (i=1;i<=confs;i++){   
         /*tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
         tmp+=sizeof(double)*iconf*sizeAmuGP+sizeof(int)*iconf;
         fseek(oAmuGPo, tmp+index, SEEK_SET);*/
         fread(&iii,sizeof(int),1,oAmuGPo);
         fread(obs_AmuGP,sizeof(double),sizeAmuGP,oAmuGPo); 
         for (j=0;j<sizeAmuGP;j++){
             sum_AmuGP[j]+=obs_AmuGP[j];
         }
         if ((i%(bin))==0){
             for (j=0;j<sizeAmuGP;j++)
                sum_AmuGP[j]/=(double) bin;
             fwrite(&ib,sizeof(int),1,oAmuGPo_b);
             fwrite(sum_AmuGP,sizeof(double),sizeAmuGP,oAmuGPo_b);
             for (j=0;j<sizeAmuGP;j++)
                sum_AmuGP[j]=0;
             ib++;
        }

   }
   
   fclose(oAmuGPo_b);
    printf("AmuGP binned\n");

   double *obs_VmuGP=(double*) malloc(sizeof(double)*(sizeVmuGP));
   double *sum_VmuGP=(double*) calloc(sizeVmuGP,sizeof(double));
   
   
   ib=0;
   for (i=1;i<=confs;i++){   
         /*tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
         tmp+=sizeof(double)*iconf*sizeVmuGP+sizeof(int)*iconf;
         fseek(oVmuGPo, tmp+index, SEEK_SET);*/
         fread(&iii,sizeof(int),1,oVmuGPo);
         fread(obs_VmuGP,sizeof(double),sizeVmuGP,oVmuGPo); 
         for (j=0;j<sizeVmuGP;j++){
             sum_VmuGP[j]+=obs_VmuGP[j];
         }
         if ((i%(bin))==0){
             for (j=0;j<sizeVmuGP;j++)
                sum_VmuGP[j]/=(double) bin;
             fwrite(&ib,sizeof(int),1,oVmuGPo_b);
             fwrite(sum_VmuGP,sizeof(double),sizeVmuGP,oVmuGPo_b);
             for (j=0;j<sizeVmuGP;j++)
                sum_VmuGP[j]=0;
             ib++;
        }

   }
   fclose(oVmuGPo_b);
    printf("VmuGP binned\n");

   
 
    return 0;   
}
