#ifndef header_file_virph_H
#define header_file_virph_H


 

#include <complex>
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <typeinfo>
#include <vector>
 
#include "mutils.hpp"
struct flavour_t
{
  int qhat;
  double kappa;
  double mu;
  double su3csw;
  double u1csw;
  double cF;
  double cF_prime;
  double th1;
  double th2;
  double th3;
};

struct inv_t
{
  int isolv;
  double mu,th[3];
};

struct einv_t
{
  int isolv;
  double mu,th0[3],tht[3],off;
};

struct combination_t
{
  int i0,it,is;
  double mu1,mu2,off;
  double th0[3],tht[3],ths[3];
};

struct  header_virph
{
  int twist;
  int nf;
  int nsrc,nsrcd;
  int l0,l1,l2,l3;
  int nk,nmoms;
  double beta,ksea,musea,csw;
  double *k,*mu,**mom;
  int allocated=0; 
  int tmax;
  int x0;
  int stype;
  int ninv;
  int neinv;
  int nsolv;
  int nhits;
  int phptype;
  int z0;
  int ncomb;
  int ngsm;
  double epsgsm;
  int nqsml,nqsm0,nqsm;
  double epsqsm;
  flavour_t gflv;
  std::vector<combination_t> comb;
  std::vector<inv_t> inv;
  std::vector<einv_t> einv;
  
  
  int header_size;
  int file_size;
  int file_nconf;
} ;


  //! binary read, non-vector case
  template <class T>
  auto bin_read(T &out,FILE *file) //const -> enable_if_t<is_pod<T>::value>
  {
    int rc=fread(&out,sizeof(T),1,file);
    error(rc!=1,1,"bin_read","Reading from file  rc: %d",rc);
  }
/*  
  //! binary read, complex case
  template <class T>
  auto bin_read(complex<T> &out,FILE *file) const -> enable_if_t<is_pod<T>::value>
  {
    for(size_t ri=0;ri<2;ri++)
      bin_read(((T*)&out)[ri],file);
  }
  
  //! specialization for vector
  template <class T>
  auto bin_read(T &out,FILE *file) const -> enable_if_t<is_vector<T>::value and not is_pod<T>::value>
  {
    for(auto &it : out)
      bin_read(it,file);
  }
*/  
  //! return what read
  template <class T>
  T bin_read()
  {
    T out;
    FILE *file;
    bin_read(out,file);
    return out;
  }

  
  //! binary read, non-vector case
  template <class T>
  auto bin_write(T &out,FILE *file) //const -> enable_if_t<is_pod<T>::value>
  {
    int rc=fwrite(&out,sizeof(T),1,file);
    error(rc!=1,1,"bin_read","writing into file  rc: %d",rc);
  }
/*  
  //! binary read, complex case
  template <class T>
  auto bin_read(complex<T> &out,FILE *file) const -> enable_if_t<is_pod<T>::value>
  {
    for(size_t ri=0;ri<2;ri++)
      bin_read(((T*)&out)[ri],file);
  }
  
  //! specialization for vector
  template <class T>
  auto bin_read(T &out,FILE *file) const -> enable_if_t<is_vector<T>::value and not is_pod<T>::value>
  {
    for(auto &it : out)
      bin_read(it,file);
  }
*/  
  //! return what read
  template <class T>
  T bin_write()
  {
    T out;
    FILE *file;
    bin_write(out,file);
    return out;
  }
  


int find_icomb_with_opposite_mu(struct  header_virph header, int icomb){
    int ci=-1;
    int i0=header.comb[icomb].i0;
    int is=header.comb[icomb].is;
    int found=0;
    
     
    for (int i =0; i<header.ncomb;i++ ){
        if (header.comb[i].i0==is   &&  header.comb[i].is==i0  )
        if (fabs(header.comb[i].off-header.comb[icomb].off )<1e-10)
        if (fabs(header.comb[i].th0[0]-header.comb[icomb].th0[0] )<1e-10)    
        if (fabs(header.comb[i].tht[0]-header.comb[icomb].tht[0] )<1e-10)        
        if (fabs(header.comb[i].ths[0]-header.comb[icomb].ths[0] )<1e-10)
        if (fabs(header.comb[i].th0[1]-header.comb[icomb].th0[1] )<1e-10)    
        if (fabs(header.comb[i].tht[1]-header.comb[icomb].tht[1] )<1e-10)        
        if (fabs(header.comb[i].ths[1]-header.comb[icomb].ths[1] )<1e-10)
        if (fabs(header.comb[i].th0[2]-header.comb[icomb].th0[2] )<1e-10)    
        if (fabs(header.comb[i].tht[2]-header.comb[icomb].tht[2] )<1e-10)        
        if (fabs(header.comb[i].ths[2]-header.comb[icomb].ths[2] )<1e-10){
                    ci=i;
                    found++;
        }
    }
    
    if (found!=1){
    auto c=header.comb[icomb];
        printf("icomb=%d\n",icomb);
        printf("i0 it is=%d  %d  %d\n",c.i0,c.it,c.is);
        printf("mu0  mut  off=%f  %f  %f\n",c.mu1,c.mu2,c.off);
        printf("th0=%f  %f  %f\n",c.th0[0],c.th0[1],c.th0[2]);
        printf("tht=%f  %f  %f\n",c.tht[0],c.tht[1],c.tht[2]);
        printf("ths=%f  %f  %f\n",c.ths[0],c.ths[1],c.ths[2]);
        printf("find_icomb_with_opposite_mu\n");
        printf("Either there is no combination with opposite mu either there are many\n");
        exit(3);
    }

    return ci;
}  

int find_icomb_k0(struct  header_virph header, int icomb){
    int ci=-1;
    int found=0;
    int foundk;
          
    for (int i =0; i<header.ncomb;i++ ){
        foundk=0;
        for (int k=0; k<3;k++){
            if (fabs(header.comb[i].th0[k]-header.comb[i].tht[k] )<1e-10){
                if (fabs(header.comb[i].th0[k]-header.comb[icomb].th0[k] )<1e-10)
                if (fabs(header.comb[i].ths[k]-header.comb[icomb].ths[k] )<1e-10)
                if (fabs(header.comb[i].mu1-header.comb[icomb].mu1 )<1e-10)    
                if (fabs(header.comb[i].mu2-header.comb[icomb].mu2 )<1e-10)        
                if (fabs(header.comb[i].off-header.comb[icomb].off )<1e-10){
                    foundk++;
                }
            }
                    
        }
        if (foundk==3){
            ci=i;
            found++;
        }
          
    }
    
    if (found!=1){
        auto c=header.comb[icomb];
        printf("icomb=%d\n",icomb);
        printf("i0 it is=%d  %d  %d\n",c.i0,c.it,c.is);
        printf("mu0  mut  off=%f  %f  %f\n",c.mu1,c.mu2,c.off);
        printf("th0=%f  %f  %f\n",c.th0[0],c.th0[1],c.th0[2]);
        printf("tht=%f  %f  %f\n",c.tht[0],c.tht[1],c.tht[2]);
        printf("ths=%f  %f  %f\n",c.ths[0],c.ths[1],c.ths[2]);
        printf("find_icomb_with_opposite_mu\n");
        printf("Either there is no combination with opposite mu either there are many\n");
        exit(3);
    }
   
        
    return ci;
}  
 

int icomb_2pt_p0k0(struct  header_virph header, int icomb){
    int ci=-1;
    int found=0;
    int foundk;
          
    for (int i =0; i<header.ncomb;i++ ){
        foundk=0;
        for (int k=0; k<3;k++){
            if (fabs(header.comb[i].th0[k]-header.comb[i].tht[k] )<1e-10){
                if (fabs(header.comb[i].th0[k]-0 )<1e-10)
                if (fabs(header.comb[i].ths[k]-0 )<1e-10)
                if (fabs(header.comb[i].mu1-header.comb[icomb].mu1 )<1e-10)    
                if (fabs(header.comb[i].mu2-header.comb[icomb].mu2 )<1e-10)        
                if (fabs(header.comb[i].off-header.comb[icomb].off )<1e-10){
                    foundk++;
                }
            }
                    
        }
        if (foundk==3){
            ci=i;
            found++;
        }
          
    }
    
    if (found!=1){
        auto c=header.comb[icomb];
        printf("icomb=%d\n",icomb);
        printf("i0 it is=%d  %d  %d\n",c.i0,c.it,c.is);
        printf("mu0  mut  off=%f  %f  %f\n",c.mu1,c.mu2,c.off);
        printf("th0=%f  %f  %f\n",c.th0[0],c.th0[1],c.th0[2]);
        printf("tht=%f  %f  %f\n",c.tht[0],c.tht[1],c.tht[2]);
        printf("ths=%f  %f  %f\n",c.ths[0],c.ths[1],c.ths[2]);
        printf("find_icomb_with_opposite_mu\n");
        printf("Either there is no combination with opposite mu either there are many\n");
        exit(3);
    }
   
        
    return ci;
}   
 
 
#endif
