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
#include "continuum_reph.hpp"
#include "tower.hpp"
#include "mutils.hpp"

#include <unistd.h>
#include <omp.h>

int head_allocated=0;

int Nboot=100;
int fdA,fdV;
struct  kinematic kinematic_2pt;
struct  kinematic_G kinematic_2pt_G;
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
int r_value(int r)
{
    int vr;
    if (r==0) vr=1;
    else if (r==1) vr=-1;
    else error(0==0,1,"r_value nor 0 neither 1","");
    return vr;
}
void get_kinematic( int ik2, int ik1,int imom2, int imom1 ){
    kinematic_2pt.ik2=ik2;
    kinematic_2pt.ik1=ik1;

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
    double k[4],p[4];
    double L[4];
    int i;
  
    L[0]=file_head.l0; L[1]=file_head.l1; L[2]=file_head.l2; L[3]=file_head.l3;         

    kinematic_2pt_G.i=index_n_twopt_G_fit(ikt,iks,imom0,imomt,imoms);
    kinematic_2pt_G.kt=file_head.k[ikt+file_head.nk];
    kinematic_2pt_G.ks=file_head.k[iks+file_head.nk];
    
    kinematic_2pt_G.ikt=ikt;
    kinematic_2pt_G.iks=iks;
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
          kinematic_2pt_G.k[i]=k[i];
          kinematic_2pt_G.p[i]=p[i];
	 // p[i]=fabs(p[i]);
    }
    error(fabs(k[0])>0.00001,1,"m_eff.c","error: k0 of the photon is not zero");    
    error(fabs(p[0])>0.00001,1,"m_eff.c","error: p0 of the photon is not zero");    
      
    kinematic_2pt_G.E_g=2.*asinh(   sqrt(k[1]*k[1]+k[2]*k[2]+k[3]*k[3])/2.   );
    kinematic_2pt_G.E_gT=sinh(kinematic_2pt_G.E_g)*(1.-exp(-L[0]*kinematic_2pt_G.E_g));
    
    //k3=kinematic_2pt_G.E_g*k[3]/fabs(k[3]);
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
    int i;
    int nk_old,nmoms_old;
    if(file_head.allocated==1){
        nk_old=file_head.nk;
        nmoms_old=file_head.nmoms;
    }
    
    fread(&file_head.twist,sizeof(int),1,stream);
    fread(&file_head.nf,sizeof(int),1,stream);
    fread(&file_head.nsrc,sizeof(int),1,stream);
    fread(&file_head.l0,sizeof(int),1,stream);
    fread(&file_head.l1,sizeof(int),1,stream);
    fread(&file_head.l2,sizeof(int),1,stream);
    fread(&file_head.l3,sizeof(int),1,stream);
    fread(&file_head.nk,sizeof(int),1,stream);
    fread(&file_head.nmoms,sizeof(int),1,stream);
    if(file_head.allocated==1)
        error(nk_old!=file_head.nk || nmoms_old!=file_head.nmoms,1,"read_file_head_jack", " file head nk or nmoms has changed %d %d -> %d %d",nk_old,nmoms_old,file_head.nk,file_head.nmoms );

    fread(&file_head.beta,sizeof(double),1,stream);
    fread(&file_head.ksea,sizeof(double),1,stream);
    fread(&file_head.musea,sizeof(double),1,stream);
    fread(&file_head.csw,sizeof(double),1,stream);
   
    if(file_head.allocated==0){
        file_head.k=(double*)  malloc(sizeof(double)*2*file_head.nk);
        file_head.mom=(double**) malloc(sizeof(double*)*file_head.nmoms);
        for(i=0;i<file_head.nmoms;i++) 
            file_head.mom[i]=(double*) malloc(sizeof(double)*4);
        file_head.allocated=1;
    }
    for(i=0;i<2*file_head.nk;++i)
    	fread(&file_head.k[i],sizeof(double),1,stream);
    
    
    for(i=0;i<file_head.nmoms;i++) {
        fread(&file_head.mom[i][0],sizeof(double),1,stream);
        fread(&file_head.mom[i][1],sizeof(double),1,stream);
        fread(&file_head.mom[i][2],sizeof(double),1,stream);
        fread(&file_head.mom[i][3],sizeof(double),1,stream);

    }
}
static void  write_file_head(FILE *stream)
{
    int i;    
    
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

   long int tmp;

   fread(s,sizeof(int),1,stream);
   
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
    int ir;
    if (r==0) ir=1;
    if (r==1) ir=0;
    return ir;
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
   double *obs;
   //int cik1, cik2;
   int cimom1,cimom2,cr1;
   int t,vol=0,index;
   double re,im;
   int    vol3=file_head.l1*file_head.l2*file_head.l3; 

   tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
   tmp+=sizeof(double)*iconf*size+sizeof(int)*iconf;

   obs=(double*) malloc(2*si*file_head.l0*sizeof(double)); 
   
   cimom1=index_n_minus_theta(imom1);
   cimom2=index_n_minus_theta(imom2);
   //cik1=index_n_minus_kappa(ik1);
   //cik2=index_n_minus_kappa(ik2);
   
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

void ave_polarization_gamma_A(  double **to_write,int si, double y, double *obs)
{
    double re,im;
    int vol,ix0;
    int vol3=file_head.l1*file_head.l2*file_head.l3; 
     for(ix0=0;ix0<file_head.l0;ix0++){
            re=0;vol=0;im=0;
  
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
      
            to_write[ix0][0]=to_write[ix0][0]+y*re/( (double) vol*vol3 );
            to_write[ix0][1]=to_write[ix0][1]+y*im/( (double) vol*vol3 );
   }
}

void ave_polarization_gamma_V(  double **to_write,int si, double y, double *obs)
{
    double re,im;
    int vol,ix0;
    int vol3=file_head.l1*file_head.l2*file_head.l3; 
     for(ix0=0;ix0<file_head.l0;ix0++){
	   re=0;vol=0;im=0;
       
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
       
       to_write[ix0][0]=to_write[ix0][0]+y*re/( (double) vol*vol3 );
	   to_write[ix0][1]=to_write[ix0][1]+y*im/( (double) vol*vol3 );
   }
}
void hypercharge(int ikt, int iks, double *yt, double *ys){
    if(ikt==0 && iks==0){
        *yt=2./3.;*ys=-1./3.;
    }
    
    //K
    else if (ikt==1 && iks==0 ){
        *yt=-1./3.;*ys=2./3.;
    }
    else if (ikt==2 && iks==0 ){
        *yt=-1./3.;*ys=2./3.;
    }
    else if (ikt==0 && iks==1 ){
        *yt=2./3.;*ys=-1./3.;
    }
    else if (ikt==0 && iks==2 ){
        *yt=2./3.;*ys=-1./3.;
    }
    
    //D+ and Ds+
    else if (ikt>=3  ){
        *yt=2./3.;*ys=-1./3.;
    }
    else{
        *yt=2./3.;*ys=-1./3.;
    }
   
}
void charge2(int ikt, int iks, double *yt, double *ys){
    if(ikt==0 && iks==0){
        *yt=-1./3;*ys=-1./3.;
    }
    
    //Ds
    else if (ikt==1 && iks==0 ){
        *yt=2./3.;*ys=-1./3.;
    }
    
    else if (ikt==0 && iks==1 ){
          *yt=-1./3.;*ys=2./3.;
    }
   
   
}


void read_twopt_gamma(FILE *stream,int size, int iconf , double **to_write, const char*name,int si,int ikt,int iks,int imom0,int imomt,int imoms ){
   
   long int tmp ;
   int iiconf,N,s;
   double *obs;
   int cikt, ciks;
   int cimom0,cimomt,cimoms;
   int ix0,vol,index;
   int    vol3=file_head.l1*file_head.l2*file_head.l3; 
   double re,im;
   double yt,ys;
   
   hypercharge( ikt,iks, &yt,&ys);
   if (file_head.nk==2){
        charge2( ikt,iks, &yt,&ys);
   }
   vol3=file_head.l1*file_head.l2*file_head.l3; 
   tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
   tmp+=sizeof(double)*iconf*size+sizeof(int)*iconf;
   obs=(double*) malloc(si*file_head.l0*2*sizeof(double)); 


   cimom0=index_n_minus_theta(imom0);
   cimomt=index_n_minus_theta(imomt);
   cimoms=index_n_minus_theta(imoms);
   cikt=index_n_minus_kappa(ikt);
   ciks=index_n_minus_kappa(iks);
   
   for(ix0=0;ix0<file_head.l0;ix0++){
       to_write[ix0][0]=0;
	   to_write[ix0][1]=0;
   }
   vol=0;
   
   index=sizeof(double)*2*index_n_twoptgamma( si, 0, 0, ikt, iks, imom0, imomt, imoms);
   fseek(stream, tmp+index, SEEK_SET);
   fread(obs,sizeof(double),si*file_head.l0*2,stream); 
   
   if (strcmp(name,"oAmuGPo")==0) 
        ave_polarization_gamma_A(  to_write, si,yt,  obs);
   if (strcmp(name,"oVmuGPo")==0)
        ave_polarization_gamma_V(  to_write, si,yt,  obs);
   
   index=sizeof(double)*2*index_n_twoptgamma( si, 0, 0, iks, ikt, imom0, imomt, imoms);
   fseek(stream, tmp+index, SEEK_SET);
   fread(obs,sizeof(double),si*file_head.l0*2,stream); 
  
   if (strcmp(name,"oAmuGPo")==0) 
        ave_polarization_gamma_A(  to_write, si,-ys,  obs);
   if (strcmp(name,"oVmuGPo")==0)
        ave_polarization_gamma_V(  to_write, si,ys,  obs);
   
   vol++;
   
  /* 
   if (imoms==imomt){
        index=sizeof(double)*2*index_n_twoptgamma( si, 0, 0, ikt, iks, imoms, imom0, imom0);
        fseek(stream, tmp+index, SEEK_SET);
        fread(obs,sizeof(double),si*file_head.l0*2,stream); 
        
            
       if (strcmp(name,"oAmuGPo")==0) 
            ave_polarization_gamma_A(  to_write, si,yt,  obs);
       if (strcmp(name,"oVmuGPo")==0)
            ave_polarization_gamma_V(  to_write, si,-yt,  obs);
   
        
        index=sizeof(double)*2*index_n_twoptgamma( si, 0, 0, iks, ikt, imoms, imom0, imom0);
        fseek(stream, tmp+index, SEEK_SET);
        fread(obs,sizeof(double),si*file_head.l0*2,stream); 
            
        if (strcmp(name,"oAmuGPo")==0) 
            ave_polarization_gamma_A(  to_write, si,-ys,  obs);
       if (strcmp(name,"oVmuGPo")==0)
            ave_polarization_gamma_V(  to_write, si,-ys,  obs);
   
       vol++;
      
   }
   if (imom0==imoms){
       index=sizeof(double)*2*index_n_twoptgamma( si, 0, 0, ikt, iks, imomt, imom0, imomt);
        fseek(stream, tmp+index, SEEK_SET);
        fread(obs,sizeof(double),si*file_head.l0*2,stream); 
        
            
       if (strcmp(name,"oAmuGPo")==0) 
            ave_polarization_gamma_A(  to_write, si,yt,  obs);
       if (strcmp(name,"oVmuGPo")==0)
            ave_polarization_gamma_V(  to_write, si,-yt,  obs);
   
        
        index=sizeof(double)*2*index_n_twoptgamma( si, 0, 0, iks, ikt, imomt, imom0, imomt);
        fseek(stream, tmp+index, SEEK_SET);
        fread(obs,sizeof(double),si*file_head.l0*2,stream); 
            
        if (strcmp(name,"oAmuGPo")==0) 
            ave_polarization_gamma_A(  to_write, si,-ys,  obs);
       if (strcmp(name,"oVmuGPo")==0)
            ave_polarization_gamma_V(  to_write, si,-ys,  obs);
   
       vol++;
   }*/
   for(ix0=0;ix0<file_head.l0;ix0++){
       to_write[ix0][0]/=(double) vol;
	   to_write[ix0][1]/=(double) vol;
   }
     
   
   free(obs);
}

void read_twopt_gamma_jr(FILE *stream,int size, int iconf , double **to_write, const char*name,int si,int ikt,int iks,int imom0,int imomt,int imoms,int jr ){
   
   long int tmp ;
   int iiconf,N,s;
   double *obs;
   int cikt, ciks;
   int cimom0,cimomt,cimoms;
   int ix0,vol,index;
   int    vol3=file_head.l1*file_head.l2*file_head.l3; 
   double re,im;
   double yt,ys;
   
   hypercharge( ikt,iks, &yt,&ys);
   if (file_head.nk==2){
        charge2( ikt,iks, &yt,&ys);
   }
   vol3=file_head.l1*file_head.l2*file_head.l3; 
   tmp= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
   tmp+=sizeof(double)*iconf*size+sizeof(int)*iconf;
   obs=(double*) malloc(si*file_head.l0*2*sizeof(double)); 


   cimom0=index_n_minus_theta(imom0);
   cimomt=index_n_minus_theta(imomt);
   cimoms=index_n_minus_theta(imoms);
   cikt=index_n_minus_kappa(ikt);
   ciks=index_n_minus_kappa(iks);
   
   for(ix0=0;ix0<file_head.l0;ix0++){
       to_write[ix0][0]=0;
	   to_write[ix0][1]=0;
   }
   vol=0;
   
   index=sizeof(double)*2*index_n_twoptgamma( si, 0, 0, ikt, iks, imom0, imomt, imoms);
   fseek(stream, tmp+index, SEEK_SET);
   fread(obs,sizeof(double),si*file_head.l0*2,stream); 
   
   for(ix0=0;ix0<file_head.l0;ix0++){
            re=0;vol=0;im=0;
  
            re+= obs[(jr+si*ix0)*2];
            im+= obs[(jr+si*ix0)*2+1];
      
            to_write[ix0][0]=to_write[ix0][0]+yt*re/( (double) vol3 );
            to_write[ix0][1]=to_write[ix0][1]+yt*im/( (double) vol3 );
   }
   
   
   index=sizeof(double)*2*index_n_twoptgamma( si, 0, 0, iks, ikt, imom0, imomt, imoms);
   fseek(stream, tmp+index, SEEK_SET);
   fread(obs,sizeof(double),si*file_head.l0*2,stream); 
  
   if (strcmp(name,"oAmuGPo")==0) 
       for(ix0=0;ix0<file_head.l0;ix0++){
            re=0;vol=0;im=0;
  
            re+= obs[(jr+si*ix0)*2];
            im+= obs[(jr+si*ix0)*2+1];
      
            to_write[ix0][0]=to_write[ix0][0]-ys*re/( (double) vol3 );
            to_write[ix0][1]=to_write[ix0][1]-ys*im/( (double) vol3 );
         }
   if (strcmp(name,"oVmuGPo")==0)
       for(ix0=0;ix0<file_head.l0;ix0++){
            re=0;vol=0;im=0;
  
            re+= obs[(jr+si*ix0)*2];
            im+= obs[(jr+si*ix0)*2+1];
      
            to_write[ix0][0]=to_write[ix0][0]+ys*re/( (double) vol3 );
            to_write[ix0][1]=to_write[ix0][1]+ys*im/( (double) vol3 );
       }
   
   vol++;
   
  
   for(ix0=0;ix0<file_head.l0;ix0++){
       to_write[ix0][0]/=(double) vol;
	   to_write[ix0][1]/=(double) vol;
   }
     
   
   free(obs);
}
/*
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
*/

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
     setup_single_file_jack(file_jack.M_PS,argv,"jackknife/M_{PS}_jack",Njack);
     setup_single_file_jack(file_jack.f_PS,argv,"jackknife/f_{PS}_jack",Njack);
     setup_single_file_jack(file_jack.Zf_PS,argv,"jackknife/Zf_{PS}_jack",Njack);
    
     setup_single_file_jack(file_jack.FV,argv,"jackknife/FV_jack",Njack);
     setup_single_file_jack(file_jack.FA,argv,"jackknife/FA_jack",Njack);
     
     setup_single_file_jack(file_jack.FV_autoplateaux,argv,"jackknife/FV_autoplateaux_jack",Njack);
     setup_single_file_jack(file_jack.FA_autoplateaux,argv,"jackknife/FA_autoplateaux_jack",Njack);
     
     setup_single_file_jack(file_jack.FAp,argv,"jackknife/FAp_jack",Njack);
     setup_single_file_jack(file_jack.FV_exclude,argv,"jackknife/FV_exclude_jack",Njack);
     setup_single_file_jack(file_jack.FAp1_exclude,argv,"jackknife/FAp1_jack",Njack);
     setup_single_file_jack(file_jack.xG,argv,"jackknife/xG_jack",Njack);

     setup_single_file_jack(file_jack.FA_from_H0,argv,"jackknife/FA_from_H0_jack",Njack);
     setup_single_file_jack(file_jack.FA_from_H0_autoplateaux,argv,"jackknife/FA_from_H0_autoplateaux_jack",Njack);
     setup_single_file_jack(file_jack.FV_from_H0,argv,"jackknife/FV_from_H0_jack",Njack);
     setup_single_file_jack(file_jack.FV_from_H0_autoplateaux,argv,"jackknife/FV_from_H0_autoplateaux_jack",Njack);
     setup_single_file_jack(file_jack.FV_from_H0_HA,argv,"jackknife/FV_from_H0_HA_jack",Njack);


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
   
   jack_fit.RA_autoplateaux= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   jack_fit.RV_autoplateaux= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);


}

int main(int argc, char **argv){
   int sizePP, sizeAmuP,sizeAmuGP,sizeVV,sizeVmuGP;
   int i,j,t,iG,iG0;
   clock_t t1,t2;
   double dt;
   int *sym;

   int *iconf,confs;
   double ****data,****data_bin,**tmp; 
   char c;
   double *in;

   //double ****M,****vec;
   //****projected_O;
   //double  ****lambda,****lambda0;
   
   double *fit,***y,*x,*m,*me;
   

   double ****conf_jack,**r,**met;
   int Ncorr=1;
   int t0=2;
   
   FILE  *oPPo=NULL, *oAmuPo=NULL,*oAmuGPo=NULL,*oVVo=NULL, *oVmuGPo=NULL;
   
   FILE *plateaux_masses=NULL, *plateaux_masses_GEVP=NULL; 
   FILE *plateaux_f=NULL;   
   FILE  *plateaux_RA=NULL,  *plateaux_RV=NULL;
   FILE  *plateaux_H_H0_A=NULL;
   
   char bin_name[NAMESIZE];
   
    struct jack_fit jack_fit;
   
   error(argc!=6,1,"main ",
         "usage argc!=5: ./form_factors blind/see/read_plateaux -p inpath   jack/boot  pdf/only_fit_R/no");

   error(strcmp(argv[2],"-p")!=0,1,"main ",
         "missing -p \n usage: ./form_factors blind/see/read_plateaux -p inpath   jack/boot  pdf/only_fit_R/no");

    error(strcmp(argv[4],"jack")!=0 && strcmp(argv[4],"boot")!=0 ,2,"main ",
         "choose jack or boot \n usage: ./form_factors blind/see/read_plateaux -p inpath   jack/boot   pdf/only_fit_R/no");
    
    
   char namefile[NAMESIZE];
   
   
   
   FILE *outfile        =NULL; mysprintf(namefile,NAMESIZE,"%s/out/out_E0.txt",argv[3]);        outfile=fopen(namefile,"w+");       error(outfile==NULL,1,"main ", "Unable to open %s file",namefile);
   FILE *m_pcac         =NULL; mysprintf(namefile,NAMESIZE,"%s/out/m_pcac.txt",argv[3]);        m_pcac=fopen(namefile,"w+");        error(m_pcac==NULL,1,"main ", "Unable to open %s file",namefile);
   FILE *outfile_oPp     =NULL; mysprintf(namefile,NAMESIZE,"%s/out/oPp.txt",argv[3]);  outfile_oPp=fopen(namefile,"w+");    error(outfile_oPp==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_RA     =NULL; mysprintf(namefile,NAMESIZE,"%s/out/RA.txt",argv[3]);                outfile_RA=fopen(namefile,"w+");    error(outfile_RA==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_CA     =NULL; mysprintf(namefile,NAMESIZE,"%s/out/CA.txt",argv[3]);                outfile_CA=fopen(namefile,"w+");    error(outfile_CA==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_RV     =NULL; mysprintf(namefile,NAMESIZE,"%s/out/RV.txt",argv[3]);                outfile_RV=fopen(namefile,"w+");    error(outfile_RV==NULL,1,"main ",  "Unable to open %s file",namefile);
   
   
   FILE *outfile_RA_autoplateaux     =NULL; mysprintf(namefile,NAMESIZE,"%s/out/RA_autoplateaux.txt",argv[3]);  outfile_RA_autoplateaux=open_file(namefile,"w+");
   FILE *outfile_RV_autoplateaux     =NULL; mysprintf(namefile,NAMESIZE,"%s/out/RV_autoplateaux.txt",argv[3]);  outfile_RV_autoplateaux=open_file(namefile,"w+");

   FILE *outfile_RA_vs_kp =NULL; mysprintf(namefile,NAMESIZE,"%s/out/RA_vs_kp.txt",argv[3]);  outfile_RA_vs_kp=fopen(namefile,"w+");    error(outfile_RA_vs_kp==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_FA_vs_kp =NULL; mysprintf(namefile,NAMESIZE,"%s/out/FA_vs_kp.txt",argv[3]);  outfile_FA_vs_kp=fopen(namefile,"w+");    error(outfile_FA_vs_kp==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_FAp1_vs_kp =NULL; mysprintf(namefile,NAMESIZE,"%s/out/FAp1_vs_kp.txt",argv[3]);  outfile_FAp1_vs_kp=fopen(namefile,"w+");    error(outfile_FAp1_vs_kp==NULL,1,"main ",  "Unable to open %s file",namefile);

   FILE *outfile_FAp1_vs_kp_exclude =NULL; mysprintf(namefile,NAMESIZE,"%s/out/FAp1_vs_kp_exclude.txt",argv[3]);  outfile_FAp1_vs_kp_exclude=fopen(namefile,"w+");    error(outfile_FAp1_vs_kp_exclude==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_FV_vs_kp =NULL; mysprintf(namefile,NAMESIZE,"%s/out/FV_vs_kp.txt",argv[3]);  outfile_FV_vs_kp=fopen(namefile,"w+");    error(outfile_FV_vs_kp==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_FV_vs_kp_exclude =NULL; mysprintf(namefile,NAMESIZE,"%s/out/FV_vs_kp_exclude.txt",argv[3]);  outfile_FV_vs_kp_exclude=fopen(namefile,"w+");    error(outfile_FV_vs_kp_exclude==NULL,1,"main ",  "Unable to open %s file",namefile);

   FILE *outfile_FA_autoplateaux_vs_kp =NULL; mysprintf(namefile,NAMESIZE,"%s/out/FA_autoplateaux_vs_kp.txt",argv[3]);  outfile_FA_autoplateaux_vs_kp=open_file(namefile,"w+");
   FILE *outfile_FV_autoplateaux_vs_kp =NULL; mysprintf(namefile,NAMESIZE,"%s/out/FV_autoplateaux_vs_kp.txt",argv[3]);  outfile_FV_autoplateaux_vs_kp=open_file(namefile,"w+");

   
   FILE *outfile_f       =NULL; mysprintf(namefile,NAMESIZE,"%s/out/f_PS.txt",argv[3]);          outfile_f=fopen(namefile,"w+");     error(outfile_f==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_Zf      =NULL; mysprintf(namefile,NAMESIZE,"%s/out/Zf_PS.txt",argv[3]);         outfile_Zf=fopen(namefile,"w+");    error(outfile_f==NULL,1,"main ", "Unable to open %s file",namefile);
   FILE *outfile_Rf      =NULL; mysprintf(namefile,NAMESIZE,"%s/out/Ratio_f_PS.txt",argv[3]);    outfile_Rf=fopen(namefile,"w+");    error(outfile_Rf==NULL,1,"main ",  "Unable to open %s file",namefile);
   FILE *outfile_FAoverFP=NULL; mysprintf(namefile,NAMESIZE,"%s/out/FAoverFP.txt",argv[3]);     outfile_FAoverFP=fopen(namefile,"w+");  error(outfile_Rf==NULL,1,"main ",  "Unable to open %s file",namefile);

   mysprintf(namefile,NAMESIZE,"%s/out/FAxg.txt",argv[3]);
   FILE *outfile_FAxg=open_file(namefile,"w+");

   mysprintf(namefile,NAMESIZE,"%s/out/FA_from_H0.txt",argv[3]);
   FILE *outfile_FA_from_H0=open_file(namefile,"w+");
   
   mysprintf(namefile,NAMESIZE,"%s/out/FA_from_H0_autoplateaux.txt",argv[3]);
   FILE *outfile_FA_from_H0_autoplateaux=open_file(namefile,"w+");

   mysprintf(namefile,NAMESIZE,"%s/out/FV_from_H0.txt",argv[3]);
   FILE *outfile_FV_from_H0=open_file(namefile,"w+");
   
   mysprintf(namefile,NAMESIZE,"%s/out/FV_from_H0_autoplateaux.txt",argv[3]);
   FILE *outfile_FV_from_H0_autoplateaux=open_file(namefile,"w+");
   
   
   mysprintf(namefile,NAMESIZE,"%s/out/FV_from_H0_HA.txt",argv[3]);
   FILE *outfile_FV_from_H0_HA=open_file(namefile,"w+");
   
   mysprintf(namefile,NAMESIZE,"%s/out/H_H0_A.txt",argv[3]);
   FILE *outfile_H_H0_A=open_file(namefile,"w+");
   
   mysprintf(namefile,NAMESIZE,"%s/out/H_H0_A_autoplateaux.txt",argv[3]);
   FILE *outfile_H_H0_A_autoplateaux=open_file(namefile,"w+");
   
   mysprintf(namefile,NAMESIZE,"%s/out/HmH0_V.txt",argv[3]);
   FILE *outfile_HmH0_V=open_file(namefile,"w+");

   mysprintf(namefile,NAMESIZE,"%s/out/HmH0_V_autoplateaux.txt",argv[3]);
   FILE *outfile_HmH0_V_autoplateaux=open_file(namefile,"w+");

   mysprintf(namefile,NAMESIZE,"%s/out/HmH0_V_HA.txt",argv[3]);
   FILE *outfile_HmH0_V_HA=open_file(namefile,"w+");

   mysprintf(namefile,NAMESIZE,"%s/out/FVxg.txt",argv[3]);
   FILE *outfile_FVxg=open_file(namefile,"w+");
   
   mysprintf(namefile,NAMESIZE,"%s/out/RA_POS_wuhan.txt",argv[3]);
   FILE *outfile_RA_POS_wuhan=open_file(namefile,"w+");
   
   mysprintf(namefile,NAMESIZE,"%s/out/interpolation_FA.txt",argv[3]);
   FILE *interpolation_FA=open_file(namefile,"w+");
   
   mysprintf(namefile,NAMESIZE,"%s/out/interpolation_FV.txt",argv[3]);
   FILE *interpolation_FV=open_file(namefile,"w+");

   mysprintf(namefile,NAMESIZE,"%s/out/interpolation_FA_from_H0.txt",argv[3]);
   FILE *interpolation_FA_from_H0=open_file(namefile,"w+");

   mysprintf(namefile,NAMESIZE,"%s/out/interpolation_FV_from_H0.txt",argv[3]);
   FILE *interpolation_FV_from_H0=open_file(namefile,"w+");

   
   
   double E_B,E_Pi, x_SCHET,q2,vec_pB,vec_pPi;
   int Neff,Njack,bin=1;
   
   t1=clock();
   if ( strcmp(argv[1],"read_plateaux")==0 ){
      mysprintf(namefile,NAMESIZE,"%s/plateaux_masses.txt",argv[3]);
      mysprintf(kinematic_2pt.plateau_m_ll,NAMESIZE,"%s/plateaux_masses.txt",argv[3]);
      plateaux_masses=     fopen(kinematic_2pt.plateau_m_ll,"r");          error(plateaux_masses==NULL,1,"main ", "Unable to open %s file",kinematic_2pt.plateau_m_ll);

      mysprintf(namefile,NAMESIZE,"%s/plateaux_f.txt",argv[3]);
      mysprintf(kinematic_2pt.plateau_f,NAMESIZE,"%s/plateaux_f.txt",argv[3]);
      plateaux_f=        fopen(kinematic_2pt.plateau_f,"r");    error(plateaux_f==NULL,1,"main ", "Unable to open %s file",kinematic_2pt.plateau_f); 

      mysprintf(namefile,NAMESIZE,"%s/plateaux_RA.txt",argv[3]);
      mysprintf(kinematic_2pt_G.plateau_RA,NAMESIZE,"%s/plateaux_RA.txt",argv[3]);      
      plateaux_RA=        fopen(kinematic_2pt_G.plateau_RA,"r");    error(plateaux_f==NULL,1,"main ", "Unable to open %s file",kinematic_2pt_G.plateau_RA); 
      
      mysprintf(namefile,NAMESIZE,"%s/plateaux_RV.txt",argv[3]);
      mysprintf(kinematic_2pt_G.plateau_RV,NAMESIZE,"%s/plateaux_RV.txt",argv[3]);      
      plateaux_RV=        fopen(kinematic_2pt_G.plateau_RV,"r");    error(plateaux_f==NULL,1,"main ", "Unable to open %s file",kinematic_2pt_G.plateau_RV); 
      
      mysprintf(kinematic_2pt_G.plateau_H_H0_A,NAMESIZE,"%s/plateaux_H_H0_A.txt",argv[3]);      
      plateaux_H_H0_A=open_file(kinematic_2pt_G.plateau_H_H0_A,"r");
      
  }
   // f=fopen("./meas_2pts_bin10.dat","r");
   mysprintf(bin_name,NAMESIZE,"_bin10_conf.realph.dat");
   //mysprintf(bin_name,NAMESIZE,"_conf.realph.dat");

   mysprintf(namefile,NAMESIZE,"%s/data/oPPo-ss%s",argv[3],bin_name);
   oPPo=fopen(namefile,"r"); error(oPPo==NULL,1, "main","2pt file %s not found",namefile );
   mysprintf(namefile,NAMESIZE,"%s/data/oAmuPo-ss%s",argv[3],bin_name);
   oAmuPo=fopen(namefile,"r"); error(oAmuPo==NULL,1, "main","2pt file %s not found",namefile ); 
   mysprintf(namefile,NAMESIZE,"%s/data/oAmuGPo-gs%s",argv[3],bin_name);
   oAmuGPo=fopen(namefile,"r"); error(oAmuGPo==NULL,1, "main","2pt file %s not found",namefile );   
   std::ifstream oAmuGPopp (namefile, std::ifstream::binary);  

   //mysprintf(namefile,"%s/data/oVVo-ss_conf.realph.dat",argv[3]);
  // oVVo=fopen(namefile,"r"); error(oVVo==NULL,1, "main","2pt file %s not found",namefile );
   mysprintf(namefile,NAMESIZE,"%s/data/oVmuGPo-gs%s",argv[3],bin_name);
   oVmuGPo=fopen(namefile,"r"); error(oVmuGPo==NULL,1, "main","2pt file %s not found",namefile );
   std::ifstream oVmuGPopp (namefile, std::ifstream::binary);   

   
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
 
   
   
   print_file_head(outfile);
   print_file_head(outfile_oPp);
   print_file_head(outfile_f);
   print_file_head(outfile_RA);
   print_file_head(outfile_H_H0_A);
   print_file_head(outfile_H_H0_A_autoplateaux);
   print_file_head(outfile_HmH0_V);
   print_file_head(outfile_HmH0_V_autoplateaux);
   print_file_head(outfile_HmH0_V_HA);
   print_file_head(outfile_RV);
   print_file_head(outfile_RA_autoplateaux);
   print_file_head(outfile_RV_autoplateaux);

   fflush(outfile);

   bin=1;
   Neff=confs/bin;
   Njack=Neff+1;
 
   setup_file_jack(argv,Njack);
   int ks_min=0,ks_max=3;
   if (ks_max>file_head.nk) ks_max=file_head.nk;
   int kt_min=0,kt_max=file_head.nk;
   
////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//confs=25;
//printf("number of confs set manually to %d\n",confs);
////////////////////////////////////////////////////////////////////
   iconf=(int*) malloc(sizeof(int)*confs);
   
   
   int var=22/*6*/,si_2pt,si;
   
   
   //data=       calloc_corr(confs, var*var,  file_head.l0 );
   //lambda0=    calloc_corr(Njack, var,  file_head.l0 );
   //projected_O=calloc_corr(Njack, var,  file_head.l0);
     
   int ii, imom1, imom2, ik1,ik2,ik3;
   ii=0;
   imom1=0;
   imom2=0;
   ik1=0;
   ik2=0;
   ik3=0;
   int tmin=15, tmax=20, sep=1;
   int yn;
  
   char name[NAMESIZE];
   int index;
   
   int vol=0;
   
   int cik3, cik2,cik1, cimom1,cimom2,counter,r1,r2,cr1,cr2;
  
   
   
   print_file_head(m_pcac); 
   
   
   double **mass_jack_fit,**f_PS_jack_fit,**Zf_PS_jack_fit, **oPp_jack_fit,**RA_jack_fit,**RV_jack_fit,**RA_autoplateaux_jack_fit,**RV_autoplateaux_jack_fit;
   double **FA_jack_fit,**FV_jack_fit,**kp,**xG,**FAp_jack_fit,**FA_autoplateaux_jack_fit,**FV_autoplateaux_jack_fit;
   double **H_H0,**H_H0_autoplateaux,**HmH0,**HmH0_autoplateaux,**HmH0_HA,**FA_from_H0_jack_fit,**FA_from_H0_autoplateaux_jack_fit,**FV_from_H0_jack_fit,**FV_from_H0_HA_jack_fit;
   double **FV_from_H0_autoplateaux_jack_fit;
   double **FAxg_jack_fit,**FVxg_jack_fit;
   double **fp_mx;
   mass_jack_fit=      (double**) malloc(sizeof(double*)*sizePP/file_head.l0);
   oPp_jack_fit=      (double**) malloc(sizeof(double*)*sizePP/file_head.l0);
   f_PS_jack_fit= (double**) malloc(sizeof(double*)*sizePP/file_head.l0);
   Zf_PS_jack_fit= (double**) malloc(sizeof(double*)*sizePP/file_head.l0);

   RA_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   H_H0= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   H_H0_autoplateaux= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   HmH0= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   HmH0_autoplateaux= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   HmH0_HA= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   FAp_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   RA_autoplateaux_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   FA_autoplateaux_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   FA_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   kp= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   xG= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   fp_mx= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);

   
   RV_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   FV_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   RV_autoplateaux_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   FV_autoplateaux_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);

   r=(double**) malloc(sizeof(double*)*file_head.l0);
   
   FAxg_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   FA_from_H0_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   FA_from_H0_autoplateaux_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   FV_from_H0_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   FV_from_H0_autoplateaux_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   FV_from_H0_HA_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);

   

   FVxg_jack_fit= (double**) malloc(sizeof(double*)*sizePP*file_head.nmoms/file_head.l0);
   
   for(i=0;i<file_head.l0;i++)
       r[i]=(double*) malloc(sizeof(double)*Njack);

 double a ,b;
  
   
   //M=calloc_corr(file_head.l0, Njack,  var );
   //vec=calloc_corr(file_head.l0, Njack,  var );
  // lambda=calloc_corr(file_head.l0, Njack,  var );
 t1=clock();       a=timestamp();
 
 
//   for (ii=0;ii<si_2pt;ii++){
  
   for(ik1=0;ik1<1;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
        if ( strcmp(argv[1],"read_plateaux")==0 ){
            go_to_line(plateaux_masses,ik1);
            go_to_line(plateaux_f,ik1); 
            go_to_line(plateaux_RA,ik1);
            go_to_line(plateaux_RV,ik1);
        }
   for(ik2=0;ik2<file_head.nk;ik2++){
   for(imom2=0;imom2<file_head.nmoms;imom2++){
   for(imom1=0;imom1<file_head.nmoms;imom1++){
       
       /*
       cr2=index_n_minus_r(r2);
       cr1=index_n_minus_r(r1);
       cimom1=index_n_minus_theta(imom1);
       cimom2=index_n_minus_theta(imom2);
      
 
       for (i=0;i<confs;i++){
          read_twopt(oPPo,sizePP,i ,data[i][0],1,0,ik1,ik2,imom1,imom2 );
          read_twopt(oAmuPo,sizeAmuP,i ,data[i][1],4,0,ik1,ik2 ,imom1,imom2);
          read_twopt(oAmuPo,sizeAmuP,i ,data[i][2],4,1,ik1,ik2 ,imom1,imom2 );
          read_twopt(oAmuPo,sizeAmuP,i ,data[i][3],4,2,ik1,ik2 ,imom1,imom2 );
          read_twopt(oAmuPo,sizeAmuP,i ,data[i][4],4,3,ik1,ik2 ,imom1,imom2 );

          read_twopt(oVVo,sizeVV,i ,data[i][5],4,0,ik1,ik2,imom1,imom2 );
    
       }
       symmetrise_corr(confs, 0, file_head.l0,data);
       antisymmetrise_corr(confs, 1, file_head.l0,data);
 
       data_bin=binning(confs, var*var, file_head.l0 ,data, bin);

       if( strcmp(argv[4],"jack")==0)
       conf_jack=create_jack(Neff, var*var, file_head.l0, data_bin);
       if( strcmp(argv[4],"boot")==0)
       conf_jack=create_boot(Neff,Nboot, var*var, file_head.l0, data_bin);
       
   
printf("here\n");
       i=index_n_twopt_fit(ik1,ik2,imom1,imom2);
       get_kinematic( ik1,  ik2,imom1,  imom2 );
       
    
       if( strcmp(argv[4],"jack")==0){
            mass_jack_fit[i]=compute_effective_mass(  argv, kinematic_2pt, (char*) "oPPo", conf_jack,  Njack ,plateaux_masses,outfile ,0);
            oPp_jack_fit[i]=compute_oPp_ll(  argv, kinematic_2pt,  (char*) "oPp", conf_jack, mass_jack_fit[i],  Njack ,plateaux_masses,outfile_oPp );
            Zf_PS_jack_fit[i]=compute_Zf_PS_ll(  argv, kinematic_2pt,  (char*) "oPPo", conf_jack, mass_jack_fit[i],  Njack ,plateaux_f,outfile_Zf );
            f_PS_jack_fit[i]=compute_f_PS_ll(  argv, kinematic_2pt,  (char*) "oAmuPo", conf_jack, mass_jack_fit[i],  Njack ,plateaux_f,outfile_f );
       } 
       if( strcmp(argv[4],"boot")==0){
            mass_jack_fit[i]=compute_effective_mass(  argv, kinematic_2pt, (char*) "oPPo", conf_jack,  Nboot+1 ,plateaux_masses,outfile ,0);
            oPp_jack_fit[i]=compute_oPp_ll(  argv, kinematic_2pt,  (char*) "oPp", conf_jack, mass_jack_fit[i],  Nboot+1 ,plateaux_masses,outfile_oPp );
            Zf_PS_jack_fit[i]=compute_Zf_PS_ll(  argv, kinematic_2pt,  (char*) "oPPo", conf_jack, mass_jack_fit[i],  Nboot+1 ,plateaux_f,outfile_Zf );
            f_PS_jack_fit[i]=compute_f_PS_ll(  argv, kinematic_2pt,  (char*) "oAmuPo", conf_jack, mass_jack_fit[i],  Nboot+1 ,plateaux_f,outfile_f );
       }
      
       
       free_corr(Neff, var*var, file_head.l0 ,data_bin);
       free_jack(Njack,var*var , file_head.l0, conf_jack);
        */
    }}  }} //end loop  ik1 ik2 imom2 imom1

     

     t2=clock();       b=timestamp();

       dt=(double)(t2-t1)/(double)(CLOCKS_PER_SEC);
  printf("time to read clock %f, real  %f\n",dt,b-a);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
///////////////////real photon        
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
     t1=clock();
      data=       calloc_corr(confs, var,  file_head.l0 );
  sym=(int*) malloc(sizeof(int)*var);
  int ikt,iks,imom0,imomt,imoms;   
  int i_m;
  double go;
  FILE *f_ODE=NULL;
  for(iks=ks_min;iks<ks_max;iks++){
  for(ikt=kt_min;ikt<kt_max;ikt++){     
   for(imoms=0;imoms<file_head.nmoms;imoms++){
   for(imomt=0;imomt<file_head.nmoms;imomt++){
   for(imom0=0;imom0<file_head.nmoms;imom0++){       
       
     a=timestamp();  
        if(imomt==0){  
            for (i=0;i<confs;i++){
               read_twopt(oPPo,sizePP,i ,data[i][0],1,0,ikt,iks,imom0,imoms );
               read_twopt(oAmuPo,sizeAmuP,i ,data[i][1],4,0,ikt,iks ,imom0,imoms);
           
                    
            }
        }
       
      
       for (i=0;i<confs;i++){
       //   read_twopt(oPPo,sizePP,i ,data[i][0],1,0,ikt,iks,imom0,imoms );
              read_twopt_gamma(oAmuGPo,sizeAmuGP,i ,data[i][2],"oAmuGPo",8,ikt,iks,imom0,imomt,imoms ); sym[2]=0;
              read_twopt_gamma(oAmuGPo,sizeAmuGP,i ,data[i][4],"oAmuGPo",8,ikt,iks,imom0,imom0,imoms ); sym[4]=0;
   //if(ikt==0) if (iks==0)  if (imomt==0) if (imoms==0) printf("%g   %g\n",data[i][0][0][0],data[i][0][0][1]);
              read_twopt_gamma_jr(oAmuGPo,sizeAmuGP,i ,data[i][6],"oAmuGPo",8,ikt,iks,imom0,imomt,imoms,1 ); sym[6]=0;
              read_twopt_gamma_jr(oAmuGPo,sizeAmuGP,i ,data[i][7],"oAmuGPo",8,ikt,iks,imom0,imomt,imoms,2 ); sym[7]=0;
              read_twopt_gamma_jr(oAmuGPo,sizeAmuGP,i ,data[i][8],"oAmuGPo",8,ikt,iks,imom0,imomt,imoms,5 ); sym[8]=0;
              read_twopt_gamma_jr(oAmuGPo,sizeAmuGP,i ,data[i][9],"oAmuGPo",8,ikt,iks,imom0,imomt,imoms,6 ); sym[9]=0;
              
              read_twopt_gamma_jr(oAmuGPo,sizeAmuGP,i ,data[i][10],"oAmuGPo",8,ikt,iks,imom0,imom0,imoms,1 ); sym[10]=0;
              read_twopt_gamma_jr(oAmuGPo,sizeAmuGP,i ,data[i][11],"oAmuGPo",8,ikt,iks,imom0,imom0,imoms,2 ); sym[11]=0;
              read_twopt_gamma_jr(oAmuGPo,sizeAmuGP,i ,data[i][12],"oAmuGPo",8,ikt,iks,imom0,imom0,imoms,5 ); sym[12]=0;
              read_twopt_gamma_jr(oAmuGPo,sizeAmuGP,i ,data[i][13],"oAmuGPo",8,ikt,iks,imom0,imom0,imoms,6 ); sym[13]=0;

       }
      

      
        for (i=0;i<confs;i++){
            read_twopt_gamma(oVmuGPo,sizeVmuGP,i ,data[i][3],"oVmuGPo",8,ikt,iks,imom0,imomt,imoms ); sym[3]=1;
            read_twopt_gamma(oVmuGPo,sizeVmuGP,i ,data[i][5],"oVmuGPo",8,ikt,iks,imom0,imom0,imoms ); sym[5]=1;
            
            read_twopt_gamma_jr(oVmuGPo,sizeVmuGP,i ,data[i][14],"oVmuGPo",8,ikt,iks,imom0,imomt,imoms,1 ); sym[14]=1;
            read_twopt_gamma_jr(oVmuGPo,sizeVmuGP,i ,data[i][15],"oVmuGPo",8,ikt,iks,imom0,imomt,imoms,2 ); sym[15]=1;
            read_twopt_gamma_jr(oVmuGPo,sizeVmuGP,i ,data[i][16],"oVmuGPo",8,ikt,iks,imom0,imomt,imoms,5 ); sym[16]=1;
            read_twopt_gamma_jr(oVmuGPo,sizeVmuGP,i ,data[i][17],"oVmuGPo",8,ikt,iks,imom0,imomt,imoms,6 ); sym[17]=1;
              
            read_twopt_gamma_jr(oVmuGPo,sizeVmuGP,i ,data[i][18],"oVmuGPo",8,ikt,iks,imom0,imom0,imoms,1 ); sym[18]=1;
            read_twopt_gamma_jr(oVmuGPo,sizeVmuGP,i ,data[i][19],"oVmuGPo",8,ikt,iks,imom0,imom0,imoms,2 ); sym[19]=1;
            read_twopt_gamma_jr(oVmuGPo,sizeVmuGP,i ,data[i][20],"oVmuGPo",8,ikt,iks,imom0,imom0,imoms,5 ); sym[20]=1;
            read_twopt_gamma_jr(oVmuGPo,sizeVmuGP,i ,data[i][21],"oVmuGPo",8,ikt,iks,imom0,imom0,imoms,6 ); sym[21]=1;
            
        }
        b=timestamp();
        printf("time=%f\n",b-a);
        
        
     //   read_twoptgamma_Nconfs(oAmuGPo,sizeAmuGP,confs,0 ,data,8,1,ikt,iks,imom0,imomt,imoms );
     //   read_twoptgamma_Nconfs(oVmuGPo,sizeVmuGP,confs,1 ,data,8,1,ikt,iks,imom0,imomt,imoms );

       symmetrise_corr(confs, 0, file_head.l0,data);

       if(imomt==0){  
           antisymmetrise_corr(confs, 1, file_head.l0,data);
           symmetric_derivative_corr(confs, 1, file_head.l0,data);
           //forward_derivative_corr(confs, 1, file_head.l0,data);
       }
       symmetrise_corr(confs, 2, file_head.l0,data);
       symmetrise_corr(confs, 4, file_head.l0,data);
       antisymmetrise_corr(confs, 3, file_head.l0,data);
       antisymmetrise_corr(confs, 5, file_head.l0,data);
       
       
       data_bin=binning(confs, var, file_head.l0 ,data, bin);
       conf_jack=create_jack(Neff, var, file_head.l0, data_bin);

       
       i=index_n_twopt_fit(ikt,iks,imom0,imoms);       
       i_m=index_n_twopt_fit(ikt,iks,0,0);      
       get_kinematic( ikt,  iks,imom0,  imoms );
       
       get_kinematic_G( ikt,  iks,imom0,imomt,  imoms );
       iG=index_n_twopt_G_fit(ikt,iks,imom0,imomt,imoms);
       /*
       if (ikt==0 && iks==0){
           mysprintf(namefile,NAMESIZE,"%s/data_ODE/oAmuGPo-pion-mom0%d-momt%d-moms%d",argv[3],imom0,imomt,imoms);
           f_ODE=open_file(namefile,"w+");
                for (j=0;j<Njack;j++){
                    for (t=0;t<file_head.l0;t++){
                        fprintf(f_ODE,"%d   %.15g\n",t,conf_jack[j][2][t][0]);
                    }
                }
           fclose(f_ODE);    
           mysprintf(namefile,NAMESIZE,"%s/data_ODE/oVmuGPo-pion-mom0%d-momt%d-moms%d",argv[3],imom0,imomt,imoms);
           f_ODE=open_file(namefile,"w+");
                for (j=0;j<Njack;j++){
                    for (t=0;t<file_head.l0;t++){
                        fprintf(f_ODE,"%d   %.15g\n",t,conf_jack[j][3][t][1]);
                    }
                }
           fclose(f_ODE);   
           if(imomt==0){ 
                mysprintf(namefile,NAMESIZE,"%s/data_ODE/oP5P5o-pion-mom0%d-moms%d",argv[3],imom0,imoms);
                f_ODE=open_file(namefile,"w+");
                        for (j=0;j<Njack;j++){
                            for (t=0;t<file_head.l0;t++){
                                fprintf(f_ODE,"%d   %.15g\n",t,conf_jack[j][0][t][0]);
                            }
                        }
                fclose(f_ODE);
           }
       }
       if (ikt==3 && iks==1){
           mysprintf(namefile,NAMESIZE,"%s/data_ODE/oAmuGPo-Ds%d-%d-mom0%d-momt%d-moms%d",argv[3],ikt,iks,imom0,imomt,imoms);
           f_ODE=open_file(namefile,"w+");
                for (j=0;j<Njack;j++){
                    for (t=0;t<file_head.l0;t++){
                        fprintf(f_ODE,"%d   %.15g\n",t,conf_jack[j][2][t][0]);
                    }
                }
           fclose(f_ODE);     
           mysprintf(namefile,NAMESIZE,"%s/data_ODE/oVmuGPo-Ds%d-%d-mom0%d-momt%d-moms%d",argv[3],ikt,iks,imom0,imomt,imoms);
           f_ODE=open_file(namefile,"w+");
                for (j=0;j<Njack;j++){
                    for (t=0;t<file_head.l0;t++){
                        fprintf(f_ODE,"%d   %.15g\n",t,conf_jack[j][3][t][1]);
                    }
                }
           fclose(f_ODE);   
           if(imomt==0){ 
                mysprintf(namefile,NAMESIZE,"%s/data_ODE/oP5P5o-Ds%d-%d-mom0%d-moms%d",argv[3],ikt,iks,imom0,imoms);
                f_ODE=open_file(namefile,"w+");
                        for (j=0;j<Njack;j++){
                            for (t=0;t<file_head.l0;t++){
                                fprintf(f_ODE,"%d   %.15g\n",t,conf_jack[j][0][t][0]);
                            }
                        }
                fclose(f_ODE);
           }
       }
       
       if (ikt==3 && iks==0){
           mysprintf(namefile,NAMESIZE,"%s/data_ODE/oAmuGPo-D%d-%d-mom0%d-momt%d-moms%d",argv[3],ikt,iks,imom0,imomt,imoms);
           f_ODE=open_file(namefile,"w+");
                for (j=0;j<Njack;j++){
                    for (t=0;t<file_head.l0;t++){
                        fprintf(f_ODE,"%d   %.15g\n",t,conf_jack[j][2][t][0]);
                    }
                }
           fclose(f_ODE);     
           mysprintf(namefile,NAMESIZE,"%s/data_ODE/oVmuGPo-D%d-%d-mom0%d-momt%d-moms%d",argv[3],ikt,iks,imom0,imomt,imoms);
           f_ODE=open_file(namefile,"w+");
                for (j=0;j<Njack;j++){
                    for (t=0;t<file_head.l0;t++){
                        fprintf(f_ODE,"%d   %.15g\n",t,conf_jack[j][3][t][1]);
                    }
                }
           fclose(f_ODE);   
           if(imomt==0){ 
                mysprintf(namefile,NAMESIZE,"%s/data_ODE/oP5P5o-D%d-%d-mom0%d-moms%d",argv[3],ikt,iks,imom0,imoms);
                f_ODE=open_file(namefile,"w+");
                        for (j=0;j<Njack;j++){
                            for (t=0;t<file_head.l0;t++){
                                fprintf(f_ODE,"%d   %.15g\n",t,conf_jack[j][0][t][0]);
                            }
                        }
                fclose(f_ODE);
           }
       }*/
       if(imomt==0){  
            mass_jack_fit[i]=compute_effective_mass(  argv, kinematic_2pt, (char*) "oPPo", conf_jack,  Njack ,&plateaux_masses,outfile ,0, "M_{PS}^{ll}");
            oPp_jack_fit[i]=compute_oPp_ll(  argv, kinematic_2pt,  (char*) "oPp", conf_jack, mass_jack_fit[i],  Njack ,plateaux_masses,outfile_oPp, 0 );
            Zf_PS_jack_fit[i]=compute_Zf_PS_ll(  argv, kinematic_2pt,  (char*) "oPPo", conf_jack, mass_jack_fit[i],  Njack ,plateaux_masses,outfile_Zf );// !!!!!Zf_PS=(mu1+mu2) oPp/M^2 anche ad impulso diverso da zero  !!! da correggere
            f_PS_jack_fit[i]=compute_f_PS_ll(  argv, kinematic_2pt,  (char*) "oAmuPo", conf_jack, mass_jack_fit[i], oPp_jack_fit[i] ,Njack ,plateaux_f,outfile_f );
       }
       
       if (ikt==3 && iks==1  && imom0==0 && imomt==2 && imoms==0){
            FILE *aaa;
           aaa=open_file("prova_A.txt","w+");
           for (int t=0;t<file_head.l0;t++)
               fprintf(aaa,"%d   %g    %g  %g    %g\n",t,conf_jack[Njack-1][2][t][0]*2*kinematic_2pt_G.E_gT,conf_jack[Njack-1][2][t][1]*2*kinematic_2pt_G.E_gT, conf_jack[Njack-1][4][t][0],conf_jack[Njack-1][4][t][1]);
           fclose(aaa);
           aaa=open_file("prova_V.txt","w+");
           for (int t=0;t<file_head.l0;t++)
               fprintf(aaa,"%d   %g    %g   %g    %g\n",t,conf_jack[Njack-1][14][t][1]*2*kinematic_2pt_G.E_gT,conf_jack[Njack-1][15][t][1]*2*kinematic_2pt_G.E_gT, conf_jack[Njack-1][16][t][1]*2*kinematic_2pt_G.E_gT,conf_jack[Njack-1][17][t][1]*2*kinematic_2pt_G.E_gT);
           fclose(aaa);
           
           
       }
      // RA_jack_fit[iG]=compute_CAmur(  argv, kinematic_2pt_G,  (char*) "oAmuGPo", conf_jack ,  Njack ,plateaux_f,outfile_CA ,1,sym);
       RA_jack_fit[iG]=compute_Rmur(  argv, kinematic_2pt_G,  (char*) "oAmuGPo", conf_jack, mass_jack_fit[i],  mass_jack_fit[i_m],   oPp_jack_fit[i] ,  Njack ,plateaux_RA,outfile_RA ,2,sym);
       H_H0[iG]=H_over_H0(  argv, kinematic_2pt_G,  (char*) "H_H0_A", conf_jack,  mass_jack_fit[i],  mass_jack_fit[i_m],   oPp_jack_fit[i],  Njack ,plateaux_H_H0_A,outfile_H_H0_A ,2,sym);

       RA_autoplateaux_jack_fit[iG]=compute_Rmur_auto_plateau(  argv, kinematic_2pt_G,  (char*) "oAmuGPo", conf_jack, mass_jack_fit[i],  mass_jack_fit[i_m],   oPp_jack_fit[i] ,  Njack ,plateaux_RA,outfile_RA_autoplateaux ,2,sym);
       H_H0_autoplateaux[iG]=H_over_H0_autoplateaux(  argv, kinematic_2pt_G,  (char*) "H_H0_A", conf_jack, mass_jack_fit[i],  mass_jack_fit[i_m],   oPp_jack_fit[i] ,  Njack ,plateaux_RA,outfile_H_H0_A_autoplateaux ,2,sym);

      //RA_jack_fit[iG]=compute_Rmur_from_meff(  argv, kinematic_2pt_G,  (char*) "oAmuGPo", conf_jack, mass_jack_fit[i],  mass_jack_fit[i_m], Njack ,plateaux_RA,outfile_RA ,2,sym);
           
       //RA_jack_fit[iG]=compute_Rmur_from_meff(  argv, kinematic_2pt_G,  (char*) "oAmuGPo", conf_jack,   Njack ,plateaux_f,outfile_RA ,1,sym);
      // HA_jack_fit[iG]=compute_Rmur_from_meff(  argv, kinematic_2pt_G,  (char*) "oAmuGPo", conf_jack,   Njack ,plateaux_f,outfile_RA ,1,sym);

       
       RV_jack_fit[iG]=compute_Rmur(  argv, kinematic_2pt_G,  (char*) "oVmuGPo", conf_jack, mass_jack_fit[i],  mass_jack_fit[i_m],   oPp_jack_fit[i] ,  Njack ,plateaux_RV,outfile_RV, 3,sym);
       RV_autoplateaux_jack_fit[iG]=compute_Rmur_auto_plateau(  argv, kinematic_2pt_G,  (char*) "oVmuGPo", conf_jack, mass_jack_fit[i],  mass_jack_fit[i_m],   oPp_jack_fit[i] ,  Njack ,plateaux_RV,outfile_RV_autoplateaux, 3,sym);
       // RV_jack_fit[iG]=compute_Rmur_from_meff(  argv, kinematic_2pt_G,  (char*) "oVmuGPo", conf_jack, mass_jack_fit[i],  mass_jack_fit[i_m], Njack ,plateaux_RV,outfile_RV ,3,sym);
       HmH0[iG]=H_minus_H0(  argv, kinematic_2pt_G,  (char*) "HmH0_V", conf_jack,  mass_jack_fit[i],  mass_jack_fit[i_m],   oPp_jack_fit[i],  Njack ,plateaux_RV,outfile_HmH0_V ,3,sym);
       HmH0_autoplateaux[iG]=H_minus_H0_autoplateaux(  argv, kinematic_2pt_G,  (char*) "oVmuGPo", conf_jack, mass_jack_fit[i],  mass_jack_fit[i_m],   oPp_jack_fit[i] ,  Njack ,plateaux_RV,outfile_HmH0_V_autoplateaux ,3,sym);

       HmH0_HA[iG]=H_minus_H0_HA(  argv, kinematic_2pt_G,  (char*) "HmH0_HA_V", conf_jack,  mass_jack_fit[i],  mass_jack_fit[i_m],   oPp_jack_fit[i],Zf_PS_jack_fit[i_m],  Njack ,plateaux_RV,outfile_HmH0_V_HA ,3,sym);

       
       free_corr(Neff, var, file_head.l0 ,data_bin);
       free_jack(Njack,var , file_head.l0, conf_jack);
       
   }}}}} //end loop  ikt iks imom0 imomt imoms
   
   free_corr(confs, var,  file_head.l0 ,data);
   fclose(outfile);     fclose(outfile_oPp);      fclose(outfile_Zf);     fclose(outfile_f);
printf("here\n");

       t2=clock();
       dt=(double)(t2-t1)/(double)(CLOCKS_PER_SEC);
       printf("time to read clock %f, real  %f\n",dt,b-a);

       
   double HA,ecurlp,ecurlk;

       
       
        fprintf(outfile_RA_vs_kp,"#k\\cdot p      RA    err\n");
   for(iks=ks_min;iks<ks_max;iks++){
   for(ikt=kt_min;ikt<kt_max;ikt++){     
      if ( strcmp(argv[1],"read_plateaux")==0 ){
            go_to_line(plateaux_masses,ik1);
            go_to_line(plateaux_f,ik1); 
            go_to_line(plateaux_RA,ik1);
            go_to_line(plateaux_RV,ik1);
        }
   fprintf(outfile_FV_vs_kp,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(outfile_FV_vs_kp_exclude,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(outfile_FA_vs_kp,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(outfile_RA_vs_kp,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(outfile_FAp1_vs_kp,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(outfile_FAp1_vs_kp_exclude,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(outfile_FAoverFP,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );

   fprintf(outfile_RA_POS_wuhan,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );

   fprintf(outfile_FAxg,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(outfile_FVxg,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   
   fprintf(outfile_FA_from_H0,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(outfile_FA_from_H0_autoplateaux,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(outfile_FV_from_H0,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(outfile_FV_from_H0_autoplateaux,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(outfile_FV_from_H0_HA,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );

   
   fprintf(outfile_FA_autoplateaux_vs_kp,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(outfile_FV_autoplateaux_vs_kp,"#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );

   for(imoms=0;imoms<file_head.nmoms;imoms++){
   for(imomt=0;imomt<file_head.nmoms;imomt++){
   for(imom0=0;imom0<file_head.nmoms;imom0++){       
  
        i=index_n_twopt_fit(ikt,iks,imom0,imoms);       
        i_m=index_n_twopt_fit(ikt,iks,0,0);       
        iG=index_n_twopt_G_fit(ikt,iks,imom0,imomt,imoms);
        
        get_kinematic_G( ikt,  iks,imom0,imomt,  imoms );
         
       
       
       /*FA*/
       FA_jack_fit[iG]=(double*) malloc(sizeof(double)*Njack);
       FA_autoplateaux_jack_fit[iG]=(double*) malloc(sizeof(double)*Njack);
       kp[iG]=         (double*) malloc(sizeof(double)*Njack);
       xG[iG]=         (double*) malloc(sizeof(double)*Njack);
       fp_mx[iG]=         (double*) malloc(sizeof(double)*Njack);
       double *ave,sign,*ave1;
       for(j=0;j<Njack;j++){
            sign=RA_jack_fit[iG][j]/fabs(RA_jack_fit[iG][j]);
            FA_jack_fit[iG][j]=RA_jack_fit[iG][j]+sign*f_PS_jack_fit[i][j]*kinematic_2pt_G.eps1[1];
            kp[iG][j]=(mass_jack_fit[i][j]*kinematic_2pt_G.E_g)- kinematic_2pt_G.kp;
            //kp[iG][j]=(0.5.*sinh(mass_jack_fit[i][j]*2)*kinematic_2pt_G.E_g)- kinematic_2pt_G.kp;
            //kp[iG][j]=(sqrt(mass_jack_fit[i_m][j]*mass_jack_fit[i_m][j]+kinematic_2pt_G.p[3]*kinematic_2pt_G.p[3]  )*kinematic_2pt_G.E_g)- kinematic_2pt_G.kp;
            xG[iG][j]=2*kp[iG][j]/(mass_jack_fit[i_m][j]*mass_jack_fit[i_m][j]);
            FA_jack_fit[iG][j]=FA_jack_fit[iG][j]*mass_jack_fit[i_m][j]/(-kp[iG][j]*kinematic_2pt_G.eps1[1]);
            fp_mx[iG][j]=2*Zf_PS_jack_fit[i_m][j]/( mass_jack_fit[i_m][j]  *xG[iG][j]);
            
            sign=RA_autoplateaux_jack_fit[iG][j]/fabs(RA_autoplateaux_jack_fit[iG][j]);
            FA_autoplateaux_jack_fit[iG][j]=RA_autoplateaux_jack_fit[iG][j]+sign*f_PS_jack_fit[i][j]*kinematic_2pt_G.eps1[1];
            FA_autoplateaux_jack_fit[iG][j]=FA_autoplateaux_jack_fit[iG][j]*mass_jack_fit[i_m][j]/(-kp[iG][j]*kinematic_2pt_G.eps1[1]);

       }
      
        ave1=mean_and_error_jack(Njack, fp_mx[iG]);
        ave=mean_and_error_jack(Njack, xG[iG]);
        if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, RA_jack_fit[iG]);
         if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, RA_jack_fit[iG]); 
       
       fprintf(outfile_RA_vs_kp,"%g     %g     %g \n",ave[0],m[0],m[1] );
       free(m);
       
       
       
       if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, FA_jack_fit[iG]);
       if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, FA_jack_fit[iG]); 
       fprintf(outfile_FA_vs_kp,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );
       free(m);

       write_jack_bin(Njack,FA_jack_fit[iG],file_jack.FA);

       m=mean_and_error(argv[4],Njack, FA_autoplateaux_jack_fit[iG]); 

       fprintf(outfile_FA_autoplateaux_vs_kp,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );
       free(m);
       write_jack_bin(Njack,FA_autoplateaux_jack_fit[iG],file_jack.FA_autoplateaux);

       
       
/////////////////////////////////////////////////////////////////////////////////////////////
       ////////////////////////////////////////////
       FAp_jack_fit[iG]=(double*) malloc(sizeof(double)*Njack);

       for(j=0;j<Njack;j++){
            //FAp_jack_fit[iG][j]=-RA_jack_fit[iG][j]/(kinematic_2pt_G.eps1[1] );
            //FAp_jack_fit[iG][j]=FAp_jack_fit[iG][j]*mass_jack_fit[i_m][j]/(kp[iG][j]);
            FAp_jack_fit[iG][j]=(H_H0[iG][j]-1.)*mass_jack_fit[i_m][j]/(kp[iG][j]);
            FAp_jack_fit[iG][j]=FAp_jack_fit[iG][j]*Zf_PS_jack_fit[i_m][j];
            FAp_jack_fit[iG][j]+=fp_mx[iG][j];
       }
      
       if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, FAp_jack_fit[iG]);
       if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, FAp_jack_fit[iG]); 
      
       fprintf(outfile_RA_POS_wuhan,"%g   %g    %g     %g    %g    %g \n",ave[0],ave[1],m[0],m[1],ave1[0],ave1[1] );
       free(m);

       
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
       /*FA x_g/2 +f_p*/

       for(j=0;j<Njack;j++){
            FAp_jack_fit[iG][j]=-RA_jack_fit[iG][j]/(kinematic_2pt_G.eps1[1] );// * f_PS_jack_fit[i][j]);
       
       }
      
      
       if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, FAp_jack_fit[iG]);
       if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, FAp_jack_fit[iG]); 
       fprintf(outfile_FAoverFP,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );
       write_jack_bin(Njack,FAp_jack_fit[iG],file_jack.FAp);
       fprintf(outfile_FAp1_vs_kp,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );
       if ( abs(imom0- imoms)<pdiffmax ){
             if ( abs(imomt- imoms)<pdiffmax ){

       		 write_jack_bin(Njack,FA_jack_fit[iG],file_jack.FAp1_exclude);
		//fprintf(outfile_FAp1_vs_kp_exclude,"#mu1=%g   #mu2=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
       		 fprintf(outfile_FAp1_vs_kp_exclude,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );
             }
       }
       
       free(m);
       
       
////////////////FA*xG 
       
       FAxg_jack_fit[iG]=(double*) malloc(sizeof(double)*Njack);

       for(j=0;j<Njack;j++){

            sign=RA_jack_fit[iG][j]/fabs(RA_jack_fit[iG][j]);
            FAxg_jack_fit[iG][j]=RA_jack_fit[iG][j]+sign*f_PS_jack_fit[i][j]*kinematic_2pt_G.eps1[1];
            FAxg_jack_fit[iG][j]=FAxg_jack_fit[iG][j]*2/(-kinematic_2pt_G.eps1[1]*mass_jack_fit[i_m][j]);
            //   FAxg_jack_fit[iG][j]=FA_jack_fit[iG][j]*xG[iG][j];
        }
        
           
       if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, FAxg_jack_fit[iG]);
       if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, FAxg_jack_fit[iG]); 
       fprintf(outfile_FAxg,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );
       
       free(m);

////////////////// H/H0-1
       FA_from_H0_jack_fit[iG]=(double*) malloc(sizeof(double)*Njack);
       iG0=index_n_twopt_G_fit(ikt,iks,imom0,imom0,imoms);

       for(j=0;j<Njack;j++){

            FA_from_H0_jack_fit[iG][j]=(H_H0[iG][j]-1.)*mass_jack_fit[i_m][j]/(kp[iG][j]);
            FA_from_H0_jack_fit[iG][j]=FA_from_H0_jack_fit[iG][j]*Zf_PS_jack_fit[i_m][j];
        }
       if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, FA_from_H0_jack_fit[iG]);
       if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, FA_from_H0_jack_fit[iG]); 
       fprintf(outfile_FA_from_H0,"%g   %g    %g     %g   \t\t %g    %g \n",ave[0],ave[1],m[0],m[1],kinematic_2pt_G.p[3],kinematic_2pt_G.k[3] );
       
       write_jack_bin(Njack,FA_from_H0_jack_fit[iG],file_jack.FA_from_H0);

       free(m);
       
       ////////////////// H/H0-1
       FA_from_H0_autoplateaux_jack_fit[iG]=(double*) malloc(sizeof(double)*Njack);
       iG0=index_n_twopt_G_fit(ikt,iks,imom0,imom0,imoms);

       for(j=0;j<Njack;j++){

            FA_from_H0_autoplateaux_jack_fit[iG][j]=(H_H0_autoplateaux[iG][j]-1.)*mass_jack_fit[i_m][j]/(kp[iG][j]);
            FA_from_H0_autoplateaux_jack_fit[iG][j]=FA_from_H0_autoplateaux_jack_fit[iG][j]*Zf_PS_jack_fit[i_m][j];
        }
       if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, FA_from_H0_autoplateaux_jack_fit[iG]);
       if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, FA_from_H0_autoplateaux_jack_fit[iG]); 
       fprintf(outfile_FA_from_H0_autoplateaux,"%g   %g    %g     %g   \t\t %g    %g \n",ave[0],ave[1],m[0],m[1],kinematic_2pt_G.p[3],kinematic_2pt_G.k[3] );
       
       write_jack_bin(Njack,FA_from_H0_autoplateaux_jack_fit[iG],file_jack.FA_from_H0_autoplateaux);

       free(m);

///////////////////       /*FV*//////////////////// /////////////////// /////////////////// 
       FV_jack_fit[iG]=(double*) malloc(sizeof(double)*Njack);
       FV_from_H0_jack_fit[iG]=(double*) malloc(sizeof(double)*Njack);
       FV_from_H0_autoplateaux_jack_fit[iG]=(double*) malloc(sizeof(double)*Njack);
       FV_from_H0_HA_jack_fit[iG]=(double*) malloc(sizeof(double)*Njack);
       FV_autoplateaux_jack_fit[iG]=(double*) malloc(sizeof(double)*Njack);
       
       for(j=0;j<Njack;j++){
            
            FV_jack_fit[iG][j]=-RV_jack_fit[iG][j]*mass_jack_fit[i_m][j];
            FV_jack_fit[iG][j]/=(kinematic_2pt_G.E_g*kinematic_2pt_G.eps1_curl_p[1]-mass_jack_fit[i][j]*kinematic_2pt_G.eps1_curl_k[1]);
           // FV_jack_fit[iG][j]*=xG[iG][j];
            FV_autoplateaux_jack_fit[iG][j]=-RV_autoplateaux_jack_fit[iG][j]*mass_jack_fit[i_m][j];
            FV_autoplateaux_jack_fit[iG][j]/=(kinematic_2pt_G.E_g*kinematic_2pt_G.eps1_curl_p[1]-mass_jack_fit[i][j]*kinematic_2pt_G.eps1_curl_k[1]);
            
            //FV_from_H0_jack_fit[iG][j]=-HmH0[iG][j]*mass_jack_fit[i_m][j];
            //FV_from_H0_jack_fit[iG][j]/=(kinematic_2pt_G.E_g*kinematic_2pt_G.eps1_curl_p[1]-mass_jack_fit[i][j]*kinematic_2pt_G.eps1_curl_k[1]);
            FV_from_H0_jack_fit[iG][j]=HmH0[iG][j];
            
            FV_from_H0_autoplateaux_jack_fit[iG][j]=-HmH0_autoplateaux[iG][j]*mass_jack_fit[i_m][j];
            FV_from_H0_autoplateaux_jack_fit[iG][j]/=(kinematic_2pt_G.E_g*kinematic_2pt_G.eps1_curl_p[1]-mass_jack_fit[i][j]*kinematic_2pt_G.eps1_curl_k[1]);
           
            
            //FV_from_H0_HA_jack_fit[iG][j]=-HmH0_HA[iG][j]*mass_jack_fit[i_m][j]*Zf_PS_jack_fit[i][j];
            //FV_from_H0_HA_jack_fit[iG][j]/=(kinematic_2pt_G.E_g*kinematic_2pt_G.eps1_curl_p[1]-mass_jack_fit[i][j]*kinematic_2pt_G.eps1_curl_k[1]);
            FV_from_H0_HA_jack_fit[iG][j]=HmH0_HA[iG][j];
            //FV_from_H0_HA_jack_fit[iG][j]=HmH0_HA[iG][j]*mass_jack_fit[i_m][j]*Zf_PS_jack_fit[i_m][j];

       }
       
       if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, FV_jack_fit[iG]);
       if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, FV_jack_fit[iG]); 
       fprintf(outfile_FV_vs_kp,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );

       write_jack_bin(Njack,FV_jack_fit[iG],file_jack.FV);
       write_jack_bin(Njack,xG[iG],file_jack.xG);
       
       
       if ( abs(imom0- imoms)<pdiffmax ){
                  if ( abs(imomt- imoms)<pdiffmax ){
                    fprintf(outfile_FV_vs_kp_exclude,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );
       	
                   // write_jack_bin(Njack,FV_jack_fit[iG],file_jack.FV_exclude);
                  }
       }

       free(m);
       
       m=mean_and_error(argv[4],Njack, FV_autoplateaux_jack_fit[iG]); 
       fprintf(outfile_FV_autoplateaux_vs_kp,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );
       write_jack_bin(Njack,FV_autoplateaux_jack_fit[iG],file_jack.FV_autoplateaux);
       free(m);
       
       m=mean_and_error(argv[4],Njack, FV_from_H0_jack_fit[iG]); 
       fprintf(outfile_FV_from_H0,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );
       write_jack_bin(Njack,FV_from_H0_jack_fit[iG],file_jack.FV_from_H0);
       free(m);
       
       m=mean_and_error(argv[4],Njack, FV_from_H0_autoplateaux_jack_fit[iG]); 
       fprintf(outfile_FV_from_H0_autoplateaux,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );
       write_jack_bin(Njack,FV_from_H0_autoplateaux_jack_fit[iG],file_jack.FV_from_H0_autoplateaux);
       free(m);
       
       
       m=mean_and_error(argv[4],Njack, FV_from_H0_HA_jack_fit[iG]); 
       fprintf(outfile_FV_from_H0_HA,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );
       write_jack_bin(Njack,FV_from_H0_HA_jack_fit[iG],file_jack.FV_from_H0_HA);
       free(m);
       
////////////////FV*xG       
       FVxg_jack_fit[iG]=(double*) malloc(sizeof(double)*Njack);

       for(j=0;j<Njack;j++)
           FVxg_jack_fit[iG][j]=FV_jack_fit[iG][j]*xG[iG][j];
           
       if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, FVxg_jack_fit[iG]);
       if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, FVxg_jack_fit[iG]); 
       fprintf(outfile_FVxg,"%g   %g    %g     %g \n",ave[0],ave[1],m[0],m[1] );

       free(m);
       free(ave);
       free(ave1);
       
       
   }}} 
   fprintf(outfile_FAxg,"\n\n");
   fprintf(outfile_FA_from_H0,"\n\n");
   fprintf(outfile_FA_from_H0_autoplateaux,"\n\n");
   fprintf(outfile_FVxg,"\n\n");
  
   fprintf(outfile_FV_vs_kp,"\n\n");
   fprintf(outfile_FV_vs_kp_exclude,"\n\n");
   fprintf(outfile_FV_from_H0,"\n\n");
   fprintf(outfile_FV_from_H0_autoplateaux,"\n\n");
   fprintf(outfile_FV_from_H0_HA,"\n\n");

   
   fprintf(outfile_FA_vs_kp,"\n\n");
   fprintf(outfile_RA_vs_kp,"\n\n");
   fprintf(outfile_RA_POS_wuhan,"\n\n");
   fprintf(outfile_FAp1_vs_kp_exclude,"\n\n");
   fprintf(outfile_FAp1_vs_kp,"\n\n");
   fprintf(outfile_FAoverFP,"\n\n");
   
   fprintf(outfile_FV_autoplateaux_vs_kp,"\n\n");
   fprintf(outfile_FA_autoplateaux_vs_kp,"\n\n");


   }} //end loop  ikt iks imom0 imomt imoms
   
   
   double ***fit_FAp1,**int_FA_from_H0,**int_FV_from_H0;
   fit_FAp1=(double***) malloc(sizeof(double**)*file_head.nk*file_head.nk);
    double ***fit_FV;
   fit_FV=(double***) malloc(sizeof(double**)*file_head.nk*file_head.nk);
   double xmin, xmax,dx,xMAX;
   int Neff1;
   
   
   
   for(iks=ks_min;iks<ks_max;iks++){
   for(ikt=kt_min;ikt<kt_max;ikt++){     
   fprintf(interpolation_FA,"\n\n#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(interpolation_FV,"\n\n#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(interpolation_FA_from_H0,"\n\n#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   fprintf(interpolation_FV_from_H0,"\n\n#mut=%g   #mus=%g \n",file_head.k[file_head.nk+ikt],  file_head.k[file_head.nk+iks] );
   i=ikt+iks*file_head.nk;
   i_m=index_n_twopt_fit(ikt,iks,0,0);    

//   mysprintf(namefile,"FA-ikt%d-iks%d",ikt,iks); 
//   fit_FAp1[i]=fit_FX(argv, namefile, ikt, iks,xG,FA_jack_fit,  file_head , Njack ,2/* npar_fun,double*/,polynimial_degree_n);
//   free(fit_FAp1[i]);

   mysprintf(namefile,NAMESIZE,"FAp-ikt%d-iks%d",ikt,iks); 
//   fit_FAp1[i]=fit_FX(argv, namefile, ikt, iks,xG,FA_jack_fit,  file_head , Njack ,2/* npar_fun,double*/,polynimial_degree_n )
  // fit_FAp1[i]=fit_FX(argv, namefile, ikt, iks,xG,FAp_jack_fit,  file_head , Njack ,2/* npar_fun,double*/,polynimial_degree_n/*poles_degree_n*/ );

   mysprintf(namefile,NAMESIZE,"FV-ikt%d-iks%d",ikt,iks); 
   //fit_FV[i]=fit_FX(argv, namefile, ikt, iks,xG,FV_jack_fit,  file_head , Njack ,2/* npar_fun,double*/,polynimial_degree_n );
/*
      if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, fit_FAp1[i][0]);
       if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, fit_FAp1[i][0]); 
       printf("FAoverFP=(%g +-  %g)",m[0],m[1] );
     
      if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, fit_FAp1[i][1]);
       if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, fit_FAp1[1][1]); 
       printf("+(1/xG)*(    %g     %g) \n",m[0],m[1] );
      
    free(fit_FAp1[i][0]);free(fit_FAp1[i][1]);
    free(fit_FV[i][0]);free(fit_FV[i][1]);
  */    
    if (ikt==0 && iks==0){dx=0.2; xMAX=4;}
    else if (ikt==0 && iks==1){dx=0.2; xMAX=1.8;} else if (ikt==1 && iks==0){dx=0.2; xMAX=1.8;}
    else if (ikt==0 && iks==2){dx=0.2; xMAX=1.8;} else if (ikt==2 && iks==0){dx=0.2; xMAX=1.8;}
    else  {dx=0.02; xMAX=0.4;} 
   
    xmin=0;
    xmax=xmin+dx;
    for (j=0;j<((int) (xMAX/dx));j++){
        
        
        Neff1=0;
        for(imom0=0;imom0<file_head.nmoms;imom0++){       
        for(imomt=0;imomt<file_head.nmoms;imomt++){
        for(imoms=0;imoms<file_head.nmoms;imoms++){
            i=index_n_twopt_G_fit(ikt,iks,imom0,imomt,imoms);
            if (xG[i][Njack-1]>xmin &&  xG[i][Njack-1]<xmax){
                Neff1=Neff1+1;
                //printf("%g\t",xG[i][Njack-1]);
            }
        }}}
        //printf("\nikt=%d   iks=%d   xg[%g,%g]  Neff1=%d\n ",ikt,iks,xmin,xmax,Neff1);
        if (Neff1!=0){
            fit_FAp1[0]=interpolate_FX(argv, namefile, ikt, iks,xG,xmin,xmax,FA_jack_fit,  file_head , Njack ,1/* npar_fun,double*/,polynimial_degree_n/*poles_degree_n*/ );
            fit_FV[0]=interpolate_FX(argv, namefile, ikt, iks,xG,xmin,xmax,FV_jack_fit,  file_head , Njack ,1/* npar_fun,double*/,polynimial_degree_n );
            
            int_FA_from_H0=interpolate_FX(argv, namefile, ikt, iks,xG,xmin,xmax,FA_from_H0_jack_fit,  file_head , Njack ,1/* npar_fun,double*/,polynimial_degree_n );
            int_FV_from_H0=interpolate_FX(argv, namefile, ikt, iks,xG,xmin,xmax,FV_from_H0_HA_jack_fit,  file_head , Njack ,1/* npar_fun,double*/,polynimial_degree_n );
            
            m=mean_and_error(argv[4],Njack, fit_FAp1[0][0]);
            //fprintf(stdout,"%g  %g    %g   %g\n",(xmin+xmax)/2.,m[0],m[1],mass_jack_fit[i_m][Njack-1] );
            fprintf(interpolation_FA,"%g  %g    %g   %g\n",(xmin+xmax)/2.,m[0],m[1],mass_jack_fit[0][Njack-1] );
            free(m);
            
            m=mean_and_error(argv[4],Njack, fit_FV[0][0]);
            fprintf(interpolation_FV,"%g  %g    %g   %g  \n",(xmin+xmax)/2.,m[0],m[1],mass_jack_fit[0][Njack-1]);
            free(m);
     
            m=mean_and_error(argv[4],Njack, int_FA_from_H0[0]);
            fprintf(interpolation_FA_from_H0,"%g  %g    %g   %g   %g \n",(xmin+xmax)/2.,m[0],m[1],mass_jack_fit[0][Njack-1],mass_jack_fit[0][Njack-1]/Zf_PS_jack_fit[0][Njack-1]);
            free(m);

            m=mean_and_error(argv[4],Njack, int_FV_from_H0[0]);
            fprintf(interpolation_FV_from_H0,"%g  %g    %g   %g  %g \n",(xmin+xmax)/2.,m[0],m[1],mass_jack_fit[0][Njack-1],mass_jack_fit[0][Njack-1]/Zf_PS_jack_fit[0][Njack-1]);
           
            free(m);
            free(fit_FAp1[0][0]);free(fit_FAp1[0]);
            free(fit_FV[0][0]);free(fit_FV[0]);
            free(int_FA_from_H0[0]);free(int_FA_from_H0);
            free(int_FV_from_H0[0]);free(int_FV_from_H0);
        }
        xmin=xmin+dx;
        xmax=xmin+dx;
    }
   /*   for(j=0;j<Njack;j++){
          fit_FAp1[i][1][j]=fit_FAp1[i][1][j]*f_PS_jack_fit[i_m][j]*mass_jack_fit[i_m][j];
      }
      if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, fit_FAp1[i][1]);
       if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, fit_FAp1[1][1]); 
       printf("FA=(%g   +-  %g) \n",m[0],m[1] );*/
   }}
   
   free(sym);free(iconf);free(fit_FAp1);free(fit_FV);
      
   for(i=0;i<file_head.l0;i++)
      free(r[i]);
   free(r); 
   
  // fclose(outfile);
   fclose(m_pcac);
   fclose(oPPo);
  fclose(interpolation_FA); fclose(interpolation_FV); fclose(interpolation_FA_from_H0); fclose(interpolation_FV_from_H0);
   //free_corr(confs, var*var, file_head.l0 ,data);

   //for (i=0;i<confs;i++)
       //free(out[i]);
  
  for(iks=ks_min;iks<ks_max;iks++){
  for(ikt=kt_min;ikt<kt_max;ikt++){     
  for(imoms=0;imoms<file_head.nmoms;imoms++){
  for(imomt=0;imomt<file_head.nmoms;imomt++){
  for(imom0=0;imom0<file_head.nmoms;imom0++){  
      i=index_n_twopt_fit(ikt,iks,imom0,imoms);       
      iG=index_n_twopt_G_fit(ikt,iks,imom0,imomt,imoms);
      if(imomt==0){  
        
        free(mass_jack_fit[i]);
        free(oPp_jack_fit[i]);
        free(f_PS_jack_fit[i]);
        free(Zf_PS_jack_fit[i]);
      }
      free(RA_jack_fit[iG]);       
      free(H_H0[iG]);   
      free(H_H0_autoplateaux[iG]);   
      free(HmH0[iG]);       
      free(HmH0_autoplateaux[iG]);       
      free(HmH0_HA[iG]);       
      free(FAp_jack_fit[iG]);       
      free(RA_autoplateaux_jack_fit[iG]);       
      free(FA_autoplateaux_jack_fit[iG]);       
      free(FA_jack_fit[iG]);       
      free(kp[iG]);       
      free(xG[iG]);       
      free(fp_mx[iG]);       
      free(RV_jack_fit[iG]);       
      free(FV_jack_fit[iG]);       
      free(RV_autoplateaux_jack_fit[iG]);       
      free(FV_autoplateaux_jack_fit[iG]);    
      
      free(FAxg_jack_fit[iG]);    
      free(FA_from_H0_jack_fit[iG]);    
      free(FVxg_jack_fit[iG]);    
      
 
      
  }}}}}
  
  free(mass_jack_fit);
  free(oPp_jack_fit);
  free(f_PS_jack_fit);
  free(Zf_PS_jack_fit);
  free(RA_jack_fit);       
  free(H_H0);       
  free(H_H0_autoplateaux);       
  free(HmH0);       
  free(HmH0_autoplateaux);       
  free(HmH0_HA);       
  free(FAp_jack_fit);       
  free(RA_autoplateaux_jack_fit);       
  free(FA_autoplateaux_jack_fit);       
  free(FA_jack_fit);       
  free(kp);       
  free(xG);       
  free(fp_mx);       
  free(RV_jack_fit);       
  free(FV_jack_fit);       
  free(RV_autoplateaux_jack_fit);       
  free(FV_autoplateaux_jack_fit);    
      
  free(FAxg_jack_fit);    
  free(FA_from_H0_jack_fit);    
  free(FVxg_jack_fit);    
 
   

    return 0;   
}
