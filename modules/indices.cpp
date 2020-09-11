#define indices_C
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"


  int  index_n_twopt(int si,int ii,int ix0,int ik1,int ik2,int imom1,int imom2)
{
    int nk,nmoms;

    nk=file_head.nk;
    nmoms=file_head.nmoms;

    return ii+si*(ix0+file_head.l0*(imom1+nmoms*(imom2+nmoms*(ik1+nk*ik2))));
}


int  index_mesonic_twopt(int si,int ii,int ix0,int ik1,int ik2,int imom1,int imom2)
{
    int nk,nmoms;

    nk=file_head.nk;
    nmoms=file_head.nmoms;

    return ix0+file_head.l0*(ii+si*(imom1+nmoms*(imom2+nmoms*(ik1+nk*ik2))));
}

  int  index_n_twopt_fit(int ik1,int ik2,int imom1,int imom2)
{
    int nk,nmoms;

    nk=file_head.nk;
    nmoms=file_head.nmoms;

    return (imom1+nmoms*(imom2+nmoms*(ik1+nk*ik2)));
}

  int  index_ensemble_twopt_fit(struct header head,int ik1,int ik2,int imom1,int imom2)
{
    int nk,nmoms;

    nk=head.nk;
    nmoms=head.nmoms;

    return (imom1+nmoms*(imom2+nmoms*(ik1+nk*ik2)));
}


int  index_n_twopt_G_fit(int ikt,int iks, int imom0, int imomt,int imoms)
{
    int nk,nmoms;

    nk=file_head.nk;
    nmoms=file_head.nmoms;

    return (imom0+nmoms*(imomt+nmoms*(imoms+nmoms*(ikt+nk*iks))));
}
int  index_ensemble_twopt_G_fit(struct header head,int ikt,int iks, int imom0, int imomt,int imoms)
{
    int nk,nmoms;

    nk=head.nk;
    nmoms=head.nmoms;

    return (imom0+nmoms*(imomt+nmoms*(imoms+nmoms*(ikt+nk*iks))));
}

  int  index_n_threept(int si,int ii,int ix0,int ik1,int ik2,int ik3,int imom1,int imom2)
{
    int nk,nmoms;

    nk=file_head.nk;
    nmoms=file_head.nmoms;

    return ii+si*(ix0+file_head.l0*(imom1+nmoms*(imom2+nmoms*(ik1+nk*(ik2+nk*ik3)))));
}
  int  index_n_twoptgamma(int si,int ii,int ix0,int ikt,int iks,int imom0,int imomt,int imoms)
{
   int nk,nmoms;

    nk=file_head.nk;
    nmoms=file_head.nmoms;

    return ii+si*(ix0+file_head.l0*(imom0+nmoms*(imomt+nmoms*(imoms+nmoms*(ikt+nk*iks)))));
}


int  index_n_minus_kappa(int ik)
{
    int imk,i;
    double mu;
    
    mu=-file_head.k[ file_head.nk+ik ];
    imk=-1;
    for (i=0;i<file_head.nk;i++){
	if ( file_head.k[ file_head.nk+i ]==mu )
            imk=i;
    }

 //  error(imk==-1,1," index_n_minus_kappa ",  "Unable to find mass=%g",mu);
   return imk; 


}

  int  index_n_minus_theta(int imom)
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

 //  error(imth==-1,1," index_n_minus_theta ",  "Unable to find theta=%g",m3);
   return imth;
}


