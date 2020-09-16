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
  

  
 
 
 
#endif
