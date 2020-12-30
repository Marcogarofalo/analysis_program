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
#include "header_file_virph.hpp"
#include "virtualphoton_routines.hpp"

#include <unistd.h>
#include <omp.h>
#include <vector> 


struct  kinematic kinematic_2pt;
struct  kinematic_G kinematic_2pt_G;

int head_allocated=0;

int Nboot=100;
int fdA,fdV;

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

  
void get_kinematic( struct combination_t head ){
    kinematic_2pt.ik1=head.i0;
    kinematic_2pt.ik2=head.is;

    kinematic_2pt.k1=head.mu1;
    kinematic_2pt.k2=head.mu2;
    
    kinematic_2pt.mom1=head.th0[2];
    kinematic_2pt.mom2=-head.ths[2];
    if (kinematic_2pt.mom2==0) kinematic_2pt.mom2=0;
    

    kinematic_2pt.mom01=head.th0[0];
    kinematic_2pt.mom02=head.ths[0];
    
    
    kinematic_2pt.r1=1;
    kinematic_2pt.r2=-1;

}


void get_kinematic_G( struct header_virph head ,int icomb ){
    double Pi=3.141592653589793;
    double k[4],p[4];
    double L[4];
    int i;
  
    L[0]=head.tmax; L[1]=head.tmax/2; L[2]=head.tmax/2; L[3]=head.tmax/2;         

    auto c=head.comb[icomb];
    kinematic_2pt_G.i=icomb;//index_n_twopt_G_fit(head.ikt,head.iks,head.imom0,head.imomt,head.imoms);
    kinematic_2pt_G.k0=c.mu1;//file_head.k[ikt+file_head.nk];
    kinematic_2pt_G.kt=c.mu1;//file_head.k[ikt+file_head.nk];
    kinematic_2pt_G.ks=c.mu2;//file_head.k[iks+file_head.nk];
    
    kinematic_2pt_G.ik0=c.i0;//ikt;
    kinematic_2pt_G.ikt=c.it;//ikt;
    kinematic_2pt_G.iks=c.is;//iks;
    kinematic_2pt_G.rt=1;
    kinematic_2pt_G.rs=-1;
    
    if (head.phptype==0)
        kinematic_2pt_G.Twall=head.z0;
    else if (head.phptype==1)
        kinematic_2pt_G.Twall=head.tmax/2;
    
    error(head.phptype !=1 && head.phptype !=0,1,"get_kinematic", "phptype=%d non supported not 0 nor 1",head.phptype );
    
    kinematic_2pt_G.Mom0[0]=0;//file_head.mom[imom0][0];
    kinematic_2pt_G.Momt[0]=0;//file_head.mom[imomt][0];
    kinematic_2pt_G.Moms[0]=0;//file_head.mom[imoms][0];
    
    for (i=1;i<4;i++){
        kinematic_2pt_G.Mom0[i]=c.th0[i-1];//file_head.mom[imom0][i];
        if (kinematic_2pt_G.Mom0[i]==0) kinematic_2pt_G.Mom0[i]=0;
        kinematic_2pt_G.Momt[i]=c.tht[i-1];//file_head.mom[imomt][i];
        if (kinematic_2pt_G.Momt[i]==0) kinematic_2pt_G.Momt[i]=0;
        kinematic_2pt_G.Moms[i]=-c.ths[i-1];//-file_head.mom[imoms][i];
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

static void  print_file_head(FILE *stream,struct header_virph header)
{
    
    fprintf(stream,"tmax =%d\n",header.tmax);
    fprintf(stream,"x0   =%d\n",header.x0);
    fprintf(stream,"stype=%d\n",header.stype);
    fprintf(stream,"ninv =%d\n",header.ninv);
    fprintf(stream,"neinv=%d\n",header.neinv);
    fprintf(stream,"nsolv=%d\n",header.nsolv);
    fprintf(stream,"nhits=%d\n",header.nhits);
    fprintf(stream,"phptype=%d\n",header.phptype);
    fprintf(stream,"z0=%d\n",header.z0);
    fprintf(stream,"ncomb=%d\n",header.ncomb);
    fprintf(stream,"ngsm =%d\n",header.ngsm);
    fprintf(stream,"nqsml=%d\n",header.nqsml);
    fprintf(stream,"nqsm0=%d\n",header.nqsm0);
    fprintf(stream,"nqsm =%d\n",header.nqsm);
    
    fprintf(stream,"epsgsm=%f\n",header.epsgsm);
    fprintf(stream,"epsqsm=%f\n",header.epsqsm);
   
    fprintf(stream,"gflv.qhat    =%d\n",header.gflv.qhat);

    fprintf(stream,"gflv.kappa   =%f\n",header.gflv.kappa);
    fprintf(stream,"gflv.mu      =%f\n",header.gflv.mu);
    fprintf(stream,"gflv.su3csw  =%f\n",header.gflv.su3csw);
    fprintf(stream,"gflv.u1csw   =%f\n",header.gflv.u1csw);
    fprintf(stream,"gflv.cF      =%f\n",header.gflv.cF);
    fprintf(stream,"gflv.cF_prime=%f\n",header.gflv.cF_prime);
    fprintf(stream,"gflv.th1     =%f\n",header.gflv.th1);
    fprintf(stream,"gflv.th2     =%f\n",header.gflv.th2);
    fprintf(stream,"gflv.th3     =%f\n",header.gflv.th3);
    

      
    for(int icomb=0;icomb<header.ncomb;++icomb)
        {
        auto c=header.comb[icomb];
        fprintf(stream,"icomb=%d\n",icomb);
        fprintf(stream,"i0 it is=%d  %d  %d\n",c.i0,c.it,c.is);
        fprintf(stream,"mu0  mus  off=%f  %f  %f\n",c.mu1,c.mu2,c.off);
        fprintf(stream,"th0=%f  %f  %f\n",c.th0[0],c.th0[1],c.th0[2]);
        fprintf(stream,"tht=%f  %f  %f\n",c.tht[0],c.tht[1],c.tht[2]);
        fprintf(stream,"ths=%f  %f  %f\n",c.ths[0],c.ths[1],c.ths[2]);
    
    }
    
    for(int inv=0;inv<header.ninv;++inv)
    {
      auto v=header.inv[inv];
      fprintf(stream,"inversion=%d\n",inv);
      fprintf(stream,"mu=%f\n",v.mu);
      fprintf(stream,"th0=%f  %f  %f\n",v.th[0],v.th[1],v.th[2]);

      
    }
     
}

static void  read_file_head_bin(FILE *stream, struct header_virph &header)
{
    
   
    for(auto& i : {&header.tmax,&header.x0,&header.stype,&header.phptype,&header.z0,&header.ninv,&header.neinv,&header.nsolv,&header.nhits,&header.ncomb,&header.ngsm,&header.nqsml, &header.nqsm0,&header.nqsm})
        bin_read(*i,stream);
 
    header.inv.resize(header.ninv);
    header.comb.resize(header.ncomb);
    header.einv.resize(header.neinv);

    fread(&header.epsgsm,sizeof(double),1,stream);
    fread(&header.epsqsm,sizeof(double),1,stream);
    
    fread(&header.gflv.qhat,sizeof(int),1,stream);
    
    for(auto& d : {&header.gflv.kappa,&header.gflv.mu,&header.gflv.su3csw,&header.gflv.u1csw,&header.gflv.cF,&header.gflv.cF_prime,&header.gflv.th1,&header.gflv.th2,&header.gflv.th3})
        bin_read(*d,stream);
    
    for(int icomb=0;icomb<header.ncomb;++icomb)
        {
        auto& c=header.comb[icomb];
        for(auto& i : {&c.i0,&c.it,&c.is})
            bin_read(*i,stream);
        for(auto& d : {&c.mu1,&c.mu2,&c.off})
            bin_read(*d,stream);
        for(auto& d : {&c.th0[0],&c.th0[1],&c.th0[2]})
            bin_read(*d,stream);
        for(auto& d : {&c.tht[0],&c.tht[1],&c.tht[2]})
            bin_read(*d,stream);
        for(auto& d : {&c.ths[0],&c.ths[1],&c.ths[2]})
            bin_read(*d,stream);
    }
    
    for(int inv=0;inv<header.ninv;++inv)
    {
      auto& v=header.inv[inv];
      for(auto& d : {&v.mu,&v.th[0],&v.th[1],&v.th[2]})
            bin_read(*d,stream);
    }
    header.header_size=ftell(stream);      
        
  
}
static void  write_file_head(FILE *stream, struct header_virph header)
{
    
    for(auto& i : {&header.tmax,&header.x0,&header.stype,&header.phptype,&header.z0,&header.ninv,&header.neinv,&header.nsolv,&header.nhits,&header.ncomb,&header.ngsm,&header.nqsml, &header.nqsm0,&header.nqsm})
        bin_write(*i,stream);
 
   
    fwrite(&header.epsgsm,sizeof(double),1,stream);
    fwrite(&header.epsqsm,sizeof(double),1,stream);
    
    fwrite(&header.gflv.qhat,sizeof(int),1,stream);
    
    for(auto& d : {&header.gflv.kappa,&header.gflv.mu,&header.gflv.su3csw,&header.gflv.u1csw,&header.gflv.cF,&header.gflv.cF_prime,&header.gflv.th1,&header.gflv.th2,&header.gflv.th3})
        bin_write(*d,stream);
    
    for(int icomb=0;icomb<header.ncomb;++icomb)
        {
        auto& c=header.comb[icomb];
        for(auto& i : {&c.i0,&c.it,&c.is})
            bin_write(*i,stream);
        for(auto& d : {&c.mu1,&c.mu2,&c.off})
            bin_write(*d,stream);
        for(auto& d : {&c.th0[0],&c.th0[1],&c.th0[2]})
            bin_write(*d,stream);
        for(auto& d : {&c.tht[0],&c.tht[1],&c.tht[2]})
            bin_write(*d,stream);
        for(auto& d : {&c.ths[0],&c.ths[1],&c.ths[2]})
            bin_write(*d,stream);
    }
    
    for(int inv=0;inv<header.ninv;++inv)
    {
      auto& v=header.inv[inv];
      for(auto& d : {&v.mu,&v.th[0],&v.th[1],&v.th[2]})
            bin_write(*d,stream);
    }
          
}

void read_nconfs( FILE *stream,struct header_virph &header){

   long int tmp;
   int& s=header.file_size;
   int& c=header.file_nconf;

  
   fseek(stream, 0, SEEK_END);
   tmp = ftell(stream);
   tmp-= header.header_size ;
   
   //ncorr=2 HA and HV,  in general there is H_{A/V}^{\mu\alpha}
   int mus=4,alphas=4, ncorr=2;
   //  2 stands for re im;
   s=2*ncorr*mus*alphas*header.tmax*header.ncomb*header.nqsml;
   std::cout<< "size="<<s<<std::endl;

   c= (tmp)/ ( sizeof(int)+(s)*sizeof(double) );

   
   std::cout<< "confs="<<c<<std::endl;
   fseek(stream, header.header_size, SEEK_SET);
 

}



void read_nconfs_2pt( FILE *stream,struct header_virph &header){

   long int tmp;
   int& s=header.file_size;
   int& c=header.file_nconf;

  
   fseek(stream, 0, SEEK_END);
   tmp = ftell(stream);
   tmp-= header.header_size ;
   
   int ncorr=5;
   // 2 stand for re im
   s=2*header.tmax*header.ninv*header.ninv*header.nqsml*ncorr;
   std::cout<< "size="<<s<<std::endl;

   c= (tmp)/ ( sizeof(int)+(s)*sizeof(double));

   
   std::cout<< "confs="<<c<<std::endl;
   fseek(stream, header.header_size, SEEK_SET);
 

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



void read_twopt(FILE *stream,int size, int iconf , double **to_write, int icorr ,struct header_virph header, int icomb, int smearing ){
   
   long int tmp;
   double *obs;
   int vol=0,index;
   //int    vol3=1;//header.tmax*header.tmax*header.tmax/8.0; 
   auto c=header.comb[icomb];
   
   tmp=header.header_size;// sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
   tmp+=sizeof(double)*iconf*size+sizeof(int)*(iconf+1);

   obs=(double*) malloc(2*header.tmax*sizeof(double)); 
    
   int ire=0,ix0=0;
   index=sizeof(double)*(  ire+2*(ix0+header.tmax*(smearing +header.nqsml*(c.i0+header.ninv*(c.is+header.ninv*icorr))))  );
   fseek(stream, tmp+index, SEEK_SET);

   fread(obs,sizeof(double),2*header.tmax,stream); 
   
   for(int t=0;t<header.tmax;t++){
       to_write[t][0]=obs[t*2];
       to_write[t][1]=obs[t*2+1];
   }
   vol++;
  /* 
   index=sizeof(double)*(ire+2*(ix0+header.tmax*(smearing +header.nqsml*(c.is+header.ninv*(c.i0+header.ninv*icorr)))));
   fseek(stream,tmp+ index, SEEK_SET);
   fread(obs,sizeof(double),2*header.tmax,stream); 
   
   for(int t=0;t<header.tmax;t++){
       to_write[t][0]+=obs[t*2];
       to_write[t][1]+=obs[t*2+1];
   }
   vol++;
  */
   
   for(int t=0;t<header.tmax;t++){
  	   to_write[t][0]/=( (double) vol );
	   to_write[t][1]/=( (double) vol );
   }
   
   


   free(obs);

}

void electric_charge(int ikt, int iks, double *yt, double *ys){
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
void electric_charge2(int ikt, int iks, double *yt, double *ys){
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

void ave_polarization_gamma_A(  double **to_write,int si, double y, double *obs,struct header_virph header)
{
    double re,im;
    int vol,ix0;
    int vol3=1;//header.tmax*header.tmax*header.tmax/8.;
    for(ix0=0;ix0<header.tmax;ix0++){
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

void ave_polarization_gamma_V(  double **to_write,int si, double y, double *obs,struct header_virph header)
{
    double re,im;
    int vol,ix0;
    int vol3=1;//header.tmax*header.tmax*header.tmax/8.; 
     for(ix0=0;ix0<header.tmax;ix0++){
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


double  *contract_eps(double *obs ,struct header_virph header){
    double *H=(double*) malloc(sizeof(double)*4*2*2*header.tmax);  //[mu r  reim   ix0]  
    int ndim=4;
    
    for(int ix0=0;ix0<header.tmax;ix0++){
    for(int reim=0;reim<2;reim++){    
        for (int alpha=0;alpha<ndim;alpha++){
            int idH=reim+2*( alpha +ndim*( 0+ 2*ix0));
            H[idH]=0;  
            for (int mu=0;mu<ndim;mu++){
                int index=reim+2*(ix0+header.tmax*(alpha +ndim*(mu)));
                H[idH]+=obs[index]*kinematic_2pt_G.eps1[mu]; 
            }
            idH=reim+2*( alpha +ndim*( 1+ 2*ix0));
            H[idH]=0;  
            for (int mu=0;mu<ndim;mu++){
                int index=reim+2*(ix0+header.tmax*(alpha +ndim*(mu)));
                H[idH]+=obs[index]*kinematic_2pt_G.eps2[mu];         
            }
        
      }
   }}   
   
   return H;
    
}


void read_twopt_gamma(FILE *stream,int size, int iconf , double **to_write ,struct header_virph header, const char* name, int icomb, int smearing ){
   
   long int tmp;
   double *obs,*H;
   int index;
   
   tmp=header.header_size;// sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
   tmp+=sizeof(double)*iconf*size+sizeof(int)*(iconf+1);

   double yt,ys;
   electric_charge( /*ikt*/ header.comb[icomb].it,/*iks*/header.comb[icomb].is , &yt,&ys);
   if (file_head.nk==2){
        electric_charge2(/*ikt*/ header.comb[icomb].it,/*iks*/header.comb[icomb].is, &yt,&ys);
   }
   yt=-1./3.;
   ys=2./3.;
   
   yt=2./3.;
   ys=-1./3.;
   
   
   
   for(int ix0=0;ix0<header.tmax;ix0++){
       to_write[ix0][0]=0;
	   to_write[ix0][1]=0;
   }
   
   int ire=0,ix0=0;
   int alpha=0, mu=0;
   int ndim=4;
   int si=ndim*2; //mu * r
   int ncorr=2;
   int corr;
   if (strcmp(name,"oAmuGPo")==0) corr=0;
   else if (strcmp(name,"oVmuGPo")==0) corr=1;
   else error(1==1,3,"read_twopt_gamma", "name is not oAmuGPo or oVmuGPo \n name=%s",name);
   
   obs=(double*) malloc(2*header.tmax*ndim*ndim*sizeof(double)); 
                        
   index=sizeof(double)*(  ire+2*(ix0 +header.tmax*(alpha+ndim*(mu+ndim*( smearing + header.nqsml*(icomb +header.ncomb*corr  ))))) );
   fseek(stream, tmp+index, SEEK_SET);
   fread(obs,sizeof(double),2*header.tmax*ndim*ndim,stream); 
   
   H=contract_eps(obs , header);

   if (strcmp(name,"oAmuGPo")==0) 
        ave_polarization_gamma_A(  to_write, si,yt,  H,header);
   else if (strcmp(name,"oVmuGPo")==0)
        ave_polarization_gamma_V(  to_write, si,yt,  H,header);
   free(H);
   
   int ci=find_icomb_with_opposite_mu(header,icomb);
   index=sizeof(double)*(  ire+2*(ix0 +header.tmax*(alpha+ndim*(mu+ndim*( smearing + header.nqsml*(ci +header.ncomb*corr  ))))) );
   
   fseek(stream, tmp+index, SEEK_SET);
   fread(obs,sizeof(double),2*header.tmax*ndim*ndim,stream); 
   
   H=contract_eps(obs , header);

   if (strcmp(name,"oAmuGPo")==0) 
        ave_polarization_gamma_A(  to_write, si,-ys,  H,header);
   else if (strcmp(name,"oVmuGPo")==0)
        ave_polarization_gamma_V(  to_write, si,ys,  H,header);
   free(H);

   free(obs);
}


void read_twopt_gamma_jr(FILE *stream,int size, int iconf , double **to_write ,struct header_virph header, const char* name, int icomb, int smearing ,int jr){
   
   long int tmp;
   double *obs,*H;
   int index;
   
   tmp=header.header_size;// sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*(11) ;
   tmp+=sizeof(double)*iconf*size+sizeof(int)*(iconf+1);

   double yt,ys;
   electric_charge( /*ikt*/ header.comb[icomb].it,/*iks*/header.comb[icomb].is , &yt,&ys);
   if (file_head.nk==2){
        electric_charge2(/*ikt*/ header.comb[icomb].it,/*iks*/header.comb[icomb].is, &yt,&ys);
   }
   yt=-1./3.;
   ys=2./3.;
   
   yt=2./3.;
   ys=-1./3.;
   
   
   
   for(int ix0=0;ix0<header.tmax;ix0++){
       to_write[ix0][0]=0;
	   to_write[ix0][1]=0;
   }
   
   int ire=0,ix0=0;
   int alpha=0, mu=0;
   int ndim=4;
   int si=ndim*2; //mu * r
   int ncorr=2;
   int corr;
   if (strcmp(name,"oAmuGPo")==0) corr=0;
   else if (strcmp(name,"oVmuGPo")==0) corr=1;
   else error(1==1,3,"read_twopt_gamma", "name is not oAmuGPo or oVmuGPo \n name=%s",name);
   
   obs=(double*) malloc(2*header.tmax*ndim*ndim*sizeof(double)); 
                        
   index=sizeof(double)*(  ire+2*(ix0 +header.tmax*(alpha+ndim*(mu+ndim*( smearing + header.nqsml*(icomb +header.ncomb*corr  ))))) );
   fseek(stream, tmp+index, SEEK_SET);
   fread(obs,sizeof(double),2*header.tmax*ndim*ndim,stream); 
   
   H=contract_eps(obs , header);
   int vol3=1;
   for(ix0=0;ix0<file_head.l0;ix0++){
            double re=0,im=0;
  
            re+= H[(jr+si*ix0)*2];
            im+= H[(jr+si*ix0)*2+1];
      
            to_write[ix0][0]=to_write[ix0][0]+yt*re/( (double) vol3 );
            to_write[ix0][1]=to_write[ix0][1]+yt*im/( (double) vol3 );
   }
   
   int ci=find_icomb_with_opposite_mu(header,icomb);
   index=sizeof(double)*(  ire+2*(ix0 +header.tmax*(alpha+ndim*(mu+ndim*( smearing + header.nqsml*(ci +header.ncomb*corr  ))))) );
   
   fseek(stream, tmp+index, SEEK_SET);
   fread(obs,sizeof(double),2*header.tmax*ndim*ndim,stream); 
   
   H=contract_eps(obs , header);

   if (strcmp(name,"oAmuGPo")==0) 
        ave_polarization_gamma_A(  to_write, si,-ys,  H,header);
   else if (strcmp(name,"oVmuGPo")==0)
        ave_polarization_gamma_V(  to_write, si,ys,  H,header);
   free(H);
  if (strcmp(name,"oAmuGPo")==0) 
       for(ix0=0;ix0<file_head.l0;ix0++){
            double re=0,im=0;
  
            re+= obs[(jr+si*ix0)*2];
            im+= obs[(jr+si*ix0)*2+1];
      
            to_write[ix0][0]=to_write[ix0][0]-ys*re/( (double) vol3 );
            to_write[ix0][1]=to_write[ix0][1]-ys*im/( (double) vol3 );
         }
   if (strcmp(name,"oVmuGPo")==0)
       for(ix0=0;ix0<file_head.l0;ix0++){
            double re=0,im=0;
  
            re+= obs[(jr+si*ix0)*2];
            im+= obs[(jr+si*ix0)*2+1];
      
            to_write[ix0][0]=to_write[ix0][0]+ys*re/( (double) vol3 );
            to_write[ix0][1]=to_write[ix0][1]+ys*im/( (double) vol3 );
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
   
if (2+2==5){
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
   }
}
   for(ix0=0;ix0<file_head.l0;ix0++){
       to_write[ix0][0]/=(double) vol;
	   to_write[ix0][1]/=(double) vol;
   }
     
   
   free(obs);
}*/
/*
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
*/
 
void setup_single_file_jack(char  *save_name,char **argv, const char  *name,int Njack, struct header_virph header){
     FILE *f;
     mysprintf(save_name,NAMESIZE,"%s/%s",argv[3],name);
     f=fopen(save_name,"w+");
     error(f==NULL,1,"setup_file_jack ",
         "Unable to open output jackknife file %s/%s",argv[3],name);
     write_file_head(f,header);
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

void setup_file_jack(char **argv,int Njack,struct header_virph header){
     setup_single_file_jack(file_jack.M_PS,argv,"jackknife/M_{PS}_jack",Njack,header);
     setup_single_file_jack(file_jack.f_PS,argv,"jackknife/f_{PS}_jack",Njack,header);
     setup_single_file_jack(file_jack.Zf_PS,argv,"jackknife/Zf_{PS}_jack",Njack,header);
    
     setup_single_file_jack(file_jack.FV,argv,"jackknife/FV_jack",Njack,header);
     setup_single_file_jack(file_jack.FA,argv,"jackknife/FA_jack",Njack,header);
     
     setup_single_file_jack(file_jack.FV_autoplateaux,argv,"jackknife/FV_autoplateaux_jack",Njack,header);
     setup_single_file_jack(file_jack.FA_autoplateaux,argv,"jackknife/FA_autoplateaux_jack",Njack,header);
     
     setup_single_file_jack(file_jack.FAp,argv,"jackknife/FAp_jack",Njack,header);
     setup_single_file_jack(file_jack.FV_exclude,argv,"jackknife/FV_exclude_jack",Njack,header);
     setup_single_file_jack(file_jack.FAp1_exclude,argv,"jackknife/FAp1_jack",Njack,header);
     setup_single_file_jack(file_jack.xG,argv,"jackknife/xG_jack",Njack,header);

     setup_single_file_jack(file_jack.FA_from_H0,argv,"jackknife/FA_from_H0_jack",Njack,header);
     setup_single_file_jack(file_jack.FA_from_H0_autoplateaux,argv,"jackknife/FA_from_H0_autoplateaux_jack",Njack,header);
     setup_single_file_jack(file_jack.FV_from_H0,argv,"jackknife/FV_from_H0_jack",Njack,header);
     setup_single_file_jack(file_jack.FV_from_H0_autoplateaux,argv,"jackknife/FV_from_H0_autoplateaux_jack",Njack,header);
     setup_single_file_jack(file_jack.FV_from_H0_HA,argv,"jackknife/FV_from_H0_HA_jack",Njack,header);


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
   clock_t t1,t2;
   double dt;

   int *iconf,confs;
   double **tmp; 
   char c;
   double *in;

   //double ****M,****vec;
   //****projected_O;
   //double  ****lambda,****lambda0;
   
   double *fit,***y,*x,*m,*me;
   

   double **r,**met;
   int Ncorr=1;
   int t0=2;
   
   FILE  *oPPo=NULL, *oAmuPo=NULL,*oAmuGPo=NULL,*oVVo=NULL, *oVmuGPo=NULL;
   
   FILE *plateaux_masses=NULL, *plateaux_masses_GEVP=NULL; 
   FILE *plateaux_f=NULL;   
   FILE  *plateaux_RA=NULL,  *plateaux_RV=NULL;
   FILE  *plateaux_H_H0_A=NULL;
   
   char bin_name[NAMESIZE];
   struct header_virph header,header_2pt;
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
   mysprintf(namefile,NAMESIZE,"%s/out/oPp_s.txt",argv[3]);
   FILE *outfile_oPp_s     =open_file(namefile,"w+");  
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

   mysprintf(namefile,NAMESIZE,"%s/out/meffHA.txt",argv[3]);
   FILE *outfile_meffHA=open_file(namefile,"w+");
   
   mysprintf(namefile,NAMESIZE,"%s/out/meffHV.txt",argv[3]);
   FILE *outfile_meffHV=open_file(namefile,"w+");
   
   mysprintf(namefile,NAMESIZE,"%s/out/HA.txt",argv[3]);
   FILE *outfile_HA=open_file(namefile,"w+");
   mysprintf(namefile,NAMESIZE,"%s/out/HV.txt",argv[3]);
   FILE *outfile_HV=open_file(namefile,"w+");
   
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
/*   // f=fopen("./meas_2pts_bin10.dat","r");
   mysprintf(bin_name,NAMESIZE,"_bin10_conf.realph.dat");
   //mysprintf(bin_name,NAMESIZE,"_conf.realph.dat");

   mysprintf(namefile,NAMESIZE,"%s/data/oPPo-ss%s",argv[3],bin_name);
   oPPo=fopen(namefile,"r"); error(oPPo==NULL,1, "main","2pt file %s not found",namefile );
   mysprintf(namefile,NAMESIZE,"%s/data/oAmuPo-ss%s",argv[3],bin_name);
   oAmuPo=fopen(namefile,"r"); error(oAmuPo==NULL,1, "main","2pt file %s not found",namefile ); 
   mysprintf(namefile,NAMESIZE,"%s/data/oAmuGPo-gs%s",argv[3],bin_name);
   oAmuGPo=fopen(namefile,"r"); error(oAmuGPo==NULL,1, "main","2pt file %s not found",namefile );   
   //std::ifstream oAmuGPopp (namefile, std::ifstream::binary);  

   mysprintf(namefile,NAMESIZE,"%s/data/oVmuGPo-gs%s",argv[3],bin_name);
   oVmuGPo=fopen(namefile,"r"); error(oVmuGPo==NULL,1, "main","2pt file %s not found",namefile );
   //std::ifstream oVmuGPopp (namefile, std::ifstream::binary);   

  
   read_file_head_bin(oPPo);
   read_nconfs(&sizePP,&confs,oPPo);   
   printf("oPPo\t size=%d  confs=%d  \n",sizePP,confs);
 
   read_file_head_bin(oAmuPo);
   read_nconfs(&sizeAmuP,&confs,oAmuPo);   
   printf("oAmuPo\t size=%d  confs=%d \n",sizeAmuP,confs);
 
   read_file_head_bin(oAmuGPo);
   read_nconfs(&sizeAmuGP,&confs,oAmuGPo);   
   printf("oAmuGPo\t size=%d  confs=%d \n",sizeAmuGP,confs);
 

   read_file_head_bin(oVmuGPo);
   read_nconfs(&sizeVmuGP,&confs,oVmuGPo);   
   printf("oVmuGPo\t size=%d  confs=%d \n",sizeVmuGP,confs);
 */
   mysprintf(namefile,NAMESIZE,"%s/data/conf.virtualph.dat2",argv[3]);
   oPPo=open_file(namefile,"r");
   read_file_head_bin(oPPo,header_2pt);
   print_file_head(stdout,header_2pt);
   std::cout << "\n-----end header_2pt-----\n"<<  std::endl;
   
   std::cout << "the header size is " << header_2pt.header_size << " bytes?"<<  std::endl;
   read_nconfs_2pt(oPPo,header_2pt);
   std::cout << "Nconf= " << header_2pt.file_nconf <<  std::endl;
      
   for (int i =0 ; i< header_2pt.file_nconf;i++){
       int iconf;
       bin_read(iconf,oPPo);
       fseek(oPPo,header_2pt.file_size*sizeof(double),SEEK_CUR);
       std::cout << iconf <<"  ";
   }  
   std::cout <<  std::endl <<  std::endl;
   
   mysprintf(namefile,NAMESIZE,"%s/data/conf.virtualph.dat",argv[3]);
   oAmuGPo=open_file(namefile,"r");
   read_file_head_bin(oAmuGPo,header);
   //print_file_head(stdout,header);
   std::cout << "\n-----end header-----\n"<<  std::endl;
   
   std::cout << "the header size is " << header.header_size << " bytes?"<<  std::endl;
   read_nconfs(oAmuGPo,header);
   std::cout << "Nconf= " << header.file_nconf <<  std::endl;
   for (int i =0 ; i< header.file_nconf;i++){
       int iconf;
       bin_read(iconf,oAmuGPo);
       fseek(oAmuGPo,header.file_size*sizeof(double),SEEK_CUR);
       std::cout << iconf <<"  ";
   }  
   std::cout <<  std::endl;
   
   error(header_2pt.file_nconf!=header.file_nconf, 1, "main","Configuration number in file 2pt differs form 2pt_gamma file");
   confs=header_2pt.file_nconf;
   
   print_file_head(outfile,header);
   print_file_head(outfile_Zf,header);
   print_file_head(outfile_oPp,header);
   print_file_head(outfile_f,header);
   print_file_head(outfile_RA,header);
   print_file_head(outfile_H_H0_A,header);
   print_file_head(outfile_H_H0_A_autoplateaux,header);
   print_file_head(outfile_HmH0_V,header);
   print_file_head(outfile_HmH0_V_autoplateaux,header);
   print_file_head(outfile_HmH0_V_HA,header);
   print_file_head(outfile_RV,header);
   print_file_head(outfile_RA_autoplateaux,header);
   print_file_head(outfile_RV_autoplateaux,header);
   print_file_head(outfile_meffHA,header);
   print_file_head(outfile_meffHV,header);
   print_file_head(outfile_HA,header);
   print_file_head(outfile_HV,header);
   
   fflush(outfile);

   bin=1;
   Neff=confs/bin;
   if( strcmp(argv[4],"jack")==0)
                Njack=Neff+1;
   if( strcmp(argv[4],"boot")==0)
                Njack=Nbootstrap+1;

 
   setup_file_jack(argv,Njack,header);
   
   
   int var=11,*sym;
   double ****data,****data_bin,****conf_jack;
   sym=(int*) malloc(sizeof(int)*var);
   data=       calloc_corr(confs, var,  header.tmax );
   
   double **mass_jack_fit=      (double**) malloc(sizeof(double*)*header.ncomb*header.nqsml);
   double **oPp_jack_fit=      (double**) malloc(sizeof(double*)*header.ncomb*header.nqsml);
   double **oPp_s_jack_fit=      (double**) malloc(sizeof(double*)*header.ncomb*header.nqsml);
   double **H_H0=               (double**) malloc(sizeof(double*)*header.ncomb*header.nqsml);
   double **HmH0_HA=            (double**) malloc(sizeof(double*)*header.ncomb*header.nqsml);
   double **Zf_PS_jack_fit=     (double**) malloc(sizeof(double*)*header.ncomb*header.nqsml);
   
   for(int smearing=0; smearing<header.nqsml; smearing++){
      if(smearing==0) mysprintf(namefile,NAMESIZE,"M_{PS}^{ll}");
      if(smearing==1) mysprintf(namefile,NAMESIZE,"M_{PS}^{s_1l}");
   for(int icomb=0; icomb<header.ncomb; icomb++){
      // initialise the old header so that compute_effective_mass works
      file_head.nk=header.ninv;
      file_head.l0=header.tmax;
      get_kinematic(header.comb[icomb]);
      get_kinematic_G(header, icomb);

      int icombk0=find_icomb_k0(header, icomb);
      
      for (int iconf=0; iconf<confs; iconf++ ){
            read_twopt(oPPo, header_2pt.file_size, iconf ,data[iconf][0], 0/*icorr*/, header_2pt, icomb,smearing );//PP
            read_twopt(oPPo, header_2pt.file_size, iconf ,data[iconf][1], 1/*icorr*/, header_2pt, icomb,smearing );//oAmuPo
            read_twopt(oPPo, header_2pt.file_size, iconf ,data[iconf][6], 0/*icorr*/, header_2pt, icomb,0/*smearing*/ ); //PP not smeared
            
            read_twopt_gamma(oAmuGPo, header.file_size, iconf ,data[iconf][2],  header, "oAmuGPo", icomb,smearing); sym[2]=0;
            read_twopt_gamma(oAmuGPo, header.file_size, iconf ,data[iconf][4],  header, "oAmuGPo", icombk0,smearing); sym[4]=0;
            
            read_twopt_gamma(oAmuGPo, header.file_size, iconf ,data[iconf][3],  header, "oVmuGPo", icomb,smearing); sym[3]=1;
            read_twopt_gamma(oAmuGPo, header.file_size, iconf ,data[iconf][5],  header, "oVmuGPo", icombk0,smearing); sym[5]=1;

            
            read_twopt_gamma_jr(oAmuGPo, header.file_size, iconf ,data[iconf][7],  header, "oVmuGPo", icomb,smearing,1); sym[6]=1;
            read_twopt_gamma_jr(oAmuGPo, header.file_size, iconf ,data[iconf][8],  header, "oVmuGPo", icomb,smearing,2); sym[6]=1;
            read_twopt_gamma_jr(oAmuGPo, header.file_size, iconf ,data[iconf][9],  header, "oVmuGPo", icomb,smearing,5); sym[6]=1;
            read_twopt_gamma_jr(oAmuGPo, header.file_size, iconf ,data[iconf][10],  header, "oVmuGPo", icomb,smearing,6); sym[6]=1;
      }
      
      symmetrise_corr(confs, 0, header.tmax,data);
      symmetrise_corr(confs, 6, header.tmax,data);

      //symmetrise_corr(confs, 2, header.tmax,data);
      symmetrise_corr(confs, 4, header.tmax,data);
      //antisymmetrise_corr(confs, 3, header.tmax,data);
      antisymmetrise_corr(confs, 5, header.tmax,data);

      data_bin=binning(confs, var, header.tmax ,data, bin);
      conf_jack=create_resampling(argv[4],Neff, var, header.tmax, data_bin);
     
      if (smearing==0 && icomb==1){
             FILE *aaa;
           aaa=open_file("prova_A.txt","w+");
           for (int t=0;t<file_head.l0;t++)
               fprintf(aaa,"%d   %g    %g  %g    %g\n",t,conf_jack[Njack-1][2][t][0],conf_jack[Njack-1][2][t][1], conf_jack[Njack-1][4][t][0],conf_jack[Njack-1][4][t][1]);
           fclose(aaa);
           aaa=open_file("prova_V.txt","w+");
           for (int t=0;t<file_head.l0;t++)
               fprintf(aaa,"%d   %g    %g  %g    %g\n",t,conf_jack[Njack-1][6][t][1],conf_jack[Njack-1][7][t][1], conf_jack[Njack-1][8][t][1],conf_jack[Njack-1][9][t][1]);
           fclose(aaa);
            
       }
      
      int i=icomb+smearing *header.ncomb;
      int iG=icomb+smearing *header.ncomb;
      int i_m=icomb_2pt_p0k0( header,  icomb)+smearing *header.ncomb;
      int i_ml=icomb_2pt_p0k0( header,  icomb)+0 *header.ncomb;
      mass_jack_fit[i]=compute_effective_mass(  argv, kinematic_2pt, (char*) "oPPo", conf_jack,  Njack ,&plateaux_masses,outfile ,0/*index in data*/, namefile);
      Zf_PS_jack_fit[i]=compute_Zf_PS_ll(  argv, kinematic_2pt,  (char*) "oPPo", conf_jack, mass_jack_fit[i],  Njack ,plateaux_masses,outfile_Zf );
      
      oPp_jack_fit[i]=compute_oPp_ll(  argv, kinematic_2pt,  (char*) "oPp", conf_jack, mass_jack_fit[i],  Njack ,plateaux_masses,outfile_oPp,6 );
      if (smearing >0 ){
          oPp_s_jack_fit[i]=compute_oPp_s(  argv, kinematic_2pt,  (char*) "oPp_s", conf_jack, oPp_jack_fit[i],  Njack ,plateaux_masses,outfile_oPp_s ,6,0);
          for(int j=0 ; j<Njack; j++ ){
                oPp_jack_fit[i][j]=oPp_s_jack_fit[i][j];
                Zf_PS_jack_fit[i][j]= oPp_jack_fit[i][j]/ (mass_jack_fit[i][j]*sinh(mass_jack_fit[i][j])  );
                Zf_PS_jack_fit[i][j]*= (kinematic_2pt.k2+kinematic_2pt.k1);
          }
      }
      double *tmp_mass;
      tmp_mass=H_AV(  argv, kinematic_2pt_G,  (char*) "m_eff_HA", conf_jack,   mass_jack_fit[i],  mass_jack_fit[i_m], oPp_jack_fit[i], Njack, plateaux_H_H0_A, outfile_HA,2,sym);
      free(tmp_mass);
      tmp_mass=H_AV(  argv, kinematic_2pt_G,  (char*) "m_eff_HV", conf_jack,   mass_jack_fit[i],  mass_jack_fit[i_m], oPp_jack_fit[i], Njack, plateaux_RV    , outfile_HV,3,sym);
      free(tmp_mass);
      tmp_mass=meffH(  argv, kinematic_2pt_G,  (char*) "m_eff_HA", conf_jack,   mass_jack_fit[i],  mass_jack_fit[i_m], Njack, plateaux_H_H0_A, outfile_meffHA,2,sym);
      free(tmp_mass);
      tmp_mass=meffH(  argv, kinematic_2pt_G,  (char*) "m_eff_HV", conf_jack,   mass_jack_fit[i],  mass_jack_fit[i_m], Njack, plateaux_RV    , outfile_meffHV,3,sym);
      free(tmp_mass);
      
      
      H_H0[iG]=H_over_H0_vir(  argv, kinematic_2pt_G,  (char*) "H_H0_A", conf_jack,  mass_jack_fit[i],  mass_jack_fit[i_m], Njack ,plateaux_H_H0_A,outfile_H_H0_A ,2,sym);
      HmH0_HA[iG]=H_minus_H0_HA_vir(  argv, kinematic_2pt_G,  (char*) "HmH0_V_HA", conf_jack,  mass_jack_fit[i],  mass_jack_fit[i_m],   Zf_PS_jack_fit[i_ml],  Njack ,plateaux_RV,outfile_HmH0_V_HA ,3,sym);
      //free
      free_corr(Neff, var, header.tmax ,data_bin);
      free_jack(Njack,var , header.tmax, conf_jack);

   }}
   fclose(outfile_meffHA);  fclose(outfile_meffHV);
   fclose(outfile);     fclose(outfile_oPp);  fclose(outfile_oPp_s);     fclose(outfile_Zf);     fclose(outfile_f);
   printf("We have done with plateaux let's move on\n");

   int imu1=-1,imu2=-1;
   for(int smearing=0; smearing<header.nqsml; smearing++){
   for(int icomb=0; icomb<header.ncomb; icomb++){
      
      file_head.nk=header.ninv;
      file_head.l0=header.tmax;
      get_kinematic(header.comb[icomb]);
      get_kinematic_G(header, icomb);
      if(imu1!=kinematic_2pt_G.ik0  || imu2!=kinematic_2pt_G.iks ){
        fprintf(outfile_FA_from_H0,"\n\n#mut=%g   #mus=%g \n",kinematic_2pt_G.k0,  kinematic_2pt_G.ks );
        fprintf(outfile_FV_from_H0_HA,"\n\n#mut=%g   #mus=%g \n",kinematic_2pt_G.k0, kinematic_2pt_G.ks);

        imu1=kinematic_2pt_G.ik0;
        imu2=kinematic_2pt_G.iks;
      }
      
      int iG=icomb+smearing *header.ncomb;
      int i=iG;
      int i_m=icomb_2pt_p0k0( header,  icomb)+smearing *header.ncomb;
      int i_ml=icomb_2pt_p0k0( header,  icomb)+0 *header.ncomb;
      ////////////////// H/H0-1
      double *ave;
      double *xG=(double*) malloc(sizeof(double)*Njack/*header.ncomb*header.nqsml*/);
      double *kp=(double*) malloc(sizeof(double)*Njack/*header.ncomb*header.nqsml*/);
      for(int j=0;j<Njack;j++){
            
            kp[j]=(mass_jack_fit[i][j]*kinematic_2pt_G.E_g)- kinematic_2pt_G.kp;
            xG[j]=2*kp[j]/(mass_jack_fit[i_m][j]*mass_jack_fit[i_m][j]);
      }
      write_jack_bin(Njack,xG,file_jack.xG);
      ave=mean_and_error_jack(Njack, xG);

      double *FA_from_H0_jack_fit=(double*) malloc(sizeof(double)*Njack);

      for(int j=0;j<Njack;j++){
            FA_from_H0_jack_fit[j]=(H_H0[iG][j]-1.)*mass_jack_fit[i_m][j]/(kp[j]);
            FA_from_H0_jack_fit[j]=FA_from_H0_jack_fit[j]*Zf_PS_jack_fit[i_ml][j];
      }
       /*if( strcmp(argv[4],"jack")==0)
               m=mean_and_error_jack(Njack, FA_from_H0_jack_fit);
       if( strcmp(argv[4],"boot")==0)
               m=mean_and_error_boot(Njack, FA_from_H0_jack_fit); 
       */
       m=mean_and_error(argv[4],Njack, FA_from_H0_jack_fit);
       fprintf(outfile_FA_from_H0,"%g   %g    %g     %g   \t\t %g    %g \n",ave[0],ave[1],m[0],m[1],kinematic_2pt_G.p[3],kinematic_2pt_G.k[3] );
       write_jack_bin(Njack,FA_from_H0_jack_fit,file_jack.FA_from_H0);

       free(FA_from_H0_jack_fit);
       free(m);
       //////////////////////////FV
       
       double *FV_from_H0_HA_jack_fit=(double*) malloc(sizeof(double)*Njack);
       for(int j=0;j<Njack;j++){
           FV_from_H0_HA_jack_fit[j]=HmH0_HA[iG][j];
       }
       m=mean_and_error(argv[4],Njack, FV_from_H0_HA_jack_fit);
       fprintf(outfile_FV_from_H0_HA,"%g   %g    %g     %g   \t\t %g    %g \n",ave[0],ave[1],m[0],m[1],kinematic_2pt_G.p[3],kinematic_2pt_G.k[3] );
       write_jack_bin(Njack,FV_from_H0_HA_jack_fit,file_jack.FV_from_H0_HA);

       free(FV_from_H0_HA_jack_fit);
       free(m);
       
       free(ave);
    ///////////////////
   }}
   
    free_corr(confs, var, header.tmax ,data);
   free_2(header.ncomb*header.nqsml,mass_jack_fit);
   return 0;   
}
