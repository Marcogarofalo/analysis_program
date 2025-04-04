#define KandD_C


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "global.hpp"
#include "resampling.hpp"
//#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "non_linear_fit.hpp"
#include "KandD.hpp"
#include "fve.hpp"
 
// #include <omp.h>

 
double two_lines(int n, int Nvar, double *x,int Npar,double  *P){
    double r;
    
    if (n==0)
        r=P[0]+P[1]*x[0];
    if (n==1)
        r=P[2]+P[3]*x[0];
    
    return r;
    
}
double MK_chiral_FVE(int n, int Nvar, double *x,int Npar,double  *P){
    
    double MKw2=0,xi;
    double pi=3.141592653589793;
    
    double mlw=x[0], msw=x[1], w0=x[2], Mpi2=x[3], fpi=x[4], frac_Lw=x[7],  Bw=x[8];
    double fw=x[9],  MK2=x[5], fK=x[6];
    
    double    P0_w=Bw, P1_w=P[0], P3ww=P[1];
    double   Pf1w=P[2],  Pf2w=P[3],  Pf4www=P[4];
    
    
    double KM,Kf;
    
    FVE_K( Bw, fw, frac_Lw,  mlw,  msw ,Mpi2,  fpi,MK2, fK,&KM, &Kf);
    if (n==0){
        MKw2=P0_w*( mlw+msw)*(1+P1_w *mlw+(1/(w0*w0))*P3ww)*KM*KM;
    }
    else if (n==1){
        xi=2*Bw*mlw/(16*pi*pi*fw*fw);
        MKw2=Pf1w*( 1.- (3./2.)* xi*log(xi)+Pf2w*xi+(1/(w0*w0))*Pf4www)*Kf;
    }
        
      
    return MKw2;
    
}

double *MK_linear_chiral_FVE(int n, int Nvar, double *x,int Npar){
    
    double MKw2=0,xi;
    double pi=3.141592653589793;
    
    double mlw=x[0], msw=x[1], w0=x[2], Mpi2=x[3], fpi=x[4], frac_Lw=x[7],  Bw=x[8];
    double fw=x[9],  MK2=x[5], fK=x[6];
    
    double    P0_w=Bw;
    
    double *r;
    r=(double*) malloc(sizeof(double)*Npar);
    
    double KM,Kf;
    
    FVE_K( Bw, fw, frac_Lw,  mlw,  msw ,Mpi2,  fpi,MK2, fK,&KM, &Kf);
    if (n==0){
       // r[0]=( mlw+msw)*KM*KM;//?????????????????????????????
        r[0]=P0_w*( mlw+msw)*mlw*KM*KM;
        r[1]=P0_w*( mlw+msw)*(1./(w0*w0))*KM*KM;
        r[2]=0;
        r[3]=0;
        r[4]=0;
        //MKw2=P0_w*( mlw+msw)*(1+P1_w *mlw+(1/(w0*w0))*P3ww)*KM*KM;
    }
    else if (n==1){
        xi=2.*Bw*mlw/(16.*pi*pi*fw*fw);
        r[0]=0;
        r[1]=0;
        r[2]=(1.- (3./2.)* xi*log(xi))*Kf;
        r[3]=xi*Kf;
        r[4]=(1/(w0*w0))*Kf;
       // MKw2=Pf1w*( 1.- (3./2.)* xi*log(xi)+Pf2w*xi+(1/(w0*w0))*Pf4www)*Kf;
    }
      
    return r;
    
}
double *MK_linear_chiral_FVE_P40(int n, int Nvar, double *x,int Npar){
    
    double MKw2=0,xi;
    double pi=3.141592653589793;
    
    double mlw=x[0], msw=x[1], w0=x[2], Mpi2=x[3], fpi=x[4], frac_Lw=x[7],  Bw=x[8];
    double fw=x[9],  MK2=x[5], fK=x[6];
    
    double    P0_w=Bw;
    
    double *r;
    r=(double*) malloc(sizeof(double)*Npar);
    
    double KM,Kf;
    
    FVE_K( Bw, fw, frac_Lw,  mlw,  msw ,Mpi2,  fpi,MK2, fK,&KM, &Kf);
    if (n==0){
       // r[0]=( mlw+msw)*KM*KM;//?????????????????????????????
        r[0]=P0_w*( mlw+msw)*mlw*KM*KM;
        r[1]=P0_w*( mlw+msw)*(1./(w0*w0))*KM*KM;
        r[2]=0;
        r[3]=0;
        //MKw2=P0_w*( mlw+msw)*(1+P1_w *mlw+(1/(w0*w0))*P3ww)*KM*KM;
    }
    else if (n==1){
        xi=2.*Bw*mlw/(16.*pi*pi*fw*fw);
        r[0]=0;
        r[1]=0;
        r[2]=(1.- (3./2.)* xi*log(xi))*Kf;
        r[3]=xi*Kf;
       // MKw2=Pf1w*( 1.- (3./2.)* xi*log(xi)+Pf2w*xi+(1/(w0*w0))*Pf4www)*Kf;
    }
      
    return r;
    
}


double MD_chiral(int n, int Nvar, double *x,int Npar,double  *P){
    
    double MKw=0,xi;
    double pi=3.141592653589793;
    
    double mlw=x[0], msw=x[1], w0=x[2];
    
    double     P0=P[0], P1w=P[1], P2ww=P[2],   P3ww=P[3];
    double     P0f=P[4], P1fw=P[5], P2fww=P[6],   P3fww=P[7];
    
    
    
    if (n==0)
        MKw=P0+P1w*mlw+P2ww*mlw*mlw+P3ww*(1./(w0*w0));
    else if (n==1){
        MKw=P0f+P1fw*mlw+P2fww*mlw*mlw+P3fww*(1./(w0*w0));
    }
        
      
    return MKw;
    
}
double MD_chiral_P30(int n, int Nvar, double *x,int Npar,double  *P){
    
    double MKw=0,xi;
    double pi=3.141592653589793;
    
    double mlw=x[0], msw=x[1], w0=x[2];
    
    double     P0=P[0], P1w=P[1], P2ww=P[2],   P3ww=P[3];
    double     P0f=P[4], P1fw=P[5], P2fww=P[6],   P3fww=0;
    
    
    
    if (n==0)
        MKw=P0+P1w*mlw+P2ww*mlw*mlw+P3ww*(1./(w0*w0));
    else if (n==1){
        MKw=P0f+P1fw*mlw+P2fww*mlw*mlw+P3fww*(1./(w0*w0));
    }
        
      
    return MKw;
    
}

double MK_chiral_FVE_P40(int n, int Nvar, double *x,int Npar,double  *P){
    
    double MKw2=0,xi;
    double pi=3.141592653589793;
    
    double mlw=x[0], msw=x[1], w0=x[2], Mpi2=x[3], fpi=x[4], frac_Lw=x[7],  Bw=x[8];
    double fw=x[9],  MK2=x[5], fK=x[6];
    
    
    double    P0_w=Bw, P1_w=P[0], P3ww=P[1];
    double   Pf1w=P[2],  Pf2w=P[3],  Pf4www=0;
    
    double KM,Kf;
    
    FVE_K( Bw, fw, frac_Lw,  mlw,  msw ,Mpi2,  fpi,MK2, fK,&KM, &Kf);
    if (n==0)
        MKw2=P0_w*( mlw+msw)*(1+P1_w *mlw+(1/(w0*w0))*P3ww)*KM*KM;
    else if (n==1){
        xi=2*Bw*mlw/(16*pi*pi*fw*fw);
        MKw2=Pf1w*( 1.- (3./2.)* xi*log(xi)+Pf2w*xi+(1/(w0*w0))*Pf4www)*Kf;
    }
        
      
    return MKw2;
    
}


double MK_phys_point(int n, int Nvar, double *x,int Npar,double  *P){
    
    double MKw2=0,xi;
    double pi=3.141592653589793;
    
   
    double mlw=x[0], msw=x[1], w0=x[2], Mpi2=x[3], fpi=x[4], frac_Lw=x[7],  Bw=x[8];
    double fw=x[9],  MK2=x[5], fK=x[6];
    
    
    double    P0_w=Bw, P1_w=P[0], P3ww=P[1];
    double   Pf1w=P[2],  Pf2w=P[3],  Pf4www=P[4];
   
    
    double KM,Kf;
    
    //FVE_K( Bw, fw, frac_Lw,  mlw,  msw ,Mpi2,  fpi,MK2, fK,&KM, &Kf);
    if (n==0)
        MKw2=P0_w*( mlw+msw)*(1+P1_w *mlw);
    else if (n==1){
        xi=2*Bw*mlw/(16*pi*pi*fw*fw);
        MKw2=Pf1w*( 1.- (3./2.)* xi*log(xi)+Pf2w*xi);
    }
        
    
    return MKw2;
    
}
double MD_phys_point(int n, int Nvar, double *x,int Npar,double  *P){
    
    double MKw=0,xi;
    double pi=3.141592653589793;
    
    double mlw=x[0], msw=x[1], w0=x[2];
    
    double     P0=P[0], P1w=P[1], P2ww=P[2],   P3ww=0;
    double     P0f=P[4], P1fw=P[5], P2fww=P[6],   P3fww=0;
    
    
    
    if (n==0)
        MKw=P0+P1w*mlw+P2ww*mlw*mlw;
    else if (n==1){
        MKw=P0f+P1fw*mlw+P2fww*mlw*mlw;
    }
        
      
    return MKw;
    
}
 
 /*
 
double **fit_MK_fK_chiral_FVE(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit; 
   int i,j,e,im;  
   int Npar=6;
   int Nvar=8;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik2_min=1; ik2_max=3;
   int nk=(ik2_max-ik2_min+1);
   
   int n,count,N=2;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=ensembles;
   en[1]=ensembles;
   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
    guess[0]=2.05478;
    guess[1]=0.113021;
    guess[2]=2.54451;
    guess[3]=0.410103;
    guess[4]=4.61247;
    guess[5]=-0.272884;
   
   x=(double**) malloc(sizeof(double*)*(en_tot*nk));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot*nk));

   r=(double**) malloc(sizeof(double*)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot*nk));
        for (n=0;n<N;n++){
            for (ik2=ik2_min;ik2<=ik2_max;ik2++){
                for (i=0;i<en[n];i++){
                    y[j][i+ik2*en[n]+n*en[n]*nk]=(double*) malloc(sizeof(double)*2);
                }
            }
        }
   }
   

   count=0;
   for (n=0;n<N;n++){
       for (ik2=ik2_min;ik2<=ik2_max;ik2++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j];
                        rm[j]*=rm[j];
                    }
                    fit[e+ik2*en[n]+n*en[n]*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[e+ik2*en[n]+n*en[n]*nk]]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+ik2*en[n]+n*en[n]*nk]][0]=rm[j];
                    y[j][e+ik2*en[n]+n*en[n]*nk]][1]=fit[fit[e+ik2*en[n]+n*en[n]*nk]][1];
                }
                
                                
                x[e+ik2*en[n]+n*en[n]*nk]]=(double*) malloc(sizeof(double)*Nvar);
                x[e+ik2*en[n]+n*en[n]*nk]][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[e+ik2*en[n]+n*en[n]*nk]][1]=head[e].k[head[e].nk+ik2]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ms*w0
                
                x[e+ik2*en[n]+n*en[n]*nk]][2]=gJ[e].w0[Njack-1];//w0
                
                x[e+ik2*en[n]+n*en[n]*nk]][3]=gJ[e].M_PS_jack[0][Njack-1]*gJ[e].M_PS_jack[0][Njack-1];//MPi^2
                x[e+ik2*en[n]+n*en[n]*nk]][4]=gJ[e].f_PS_jack[0][Njack-1];//f_Pi
                x[e+ik2*en[n]+n*en[n]*nk]][5]=gJ[e].M_PS_jack[im][Njack-1]*gJ[e].M_PS_jack[im][Njack-1];//MK^2
                x[e+ik2*en[n]+n*en[n]*nk]][6]=gJ[e].f_PS_jack[im][Njack-1];//f_K
            
                x[e+ik2*en[n]+n*en[n]*nk]][7]=double(head[e].l1);//f
             //   printf("%g     %g      %g   %g\n",x[e+count][1],x[e+count][0],fit[e+count][0],fit[e+count][1]);
        }
        count+=en[n];
       }
   }

   
   
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, Mw2_fw_chiral_FVE,guess );

        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, Mw2_fw_chiral_FVE  );
              
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }                
        free(tmp);

   }     
// make_plots_MPi_fPi(en,);
   printf("w0/a[fm]     m*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err\n"); 
   double KM,Kf,K;
   count=0;
   for (n=0;n<N;n++){
        printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                                
               // x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=gJ[e].w0[Njack-1];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][Njack-1]*gJ[e].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
                
                FVE( gJ[e].w0[Njack-1] ,r[2][Njack-1], r[4][Njack-1], r[0][Njack-1], r[1][Njack-1],head[e].l1, x[e+count][0] ,  x[e+count][2], x[e+count][3], &KM, &Kf);
                if (n==0)
                    K=KM*KM;
                else if (n==1)
                    K=Kf;
                printf("%g     %g      %g   %g  \n",x[e+count][1],x[e+count][0],fit[e+count][0]/K,fit[e+count][1]/K);
        }
        count+=en[n];
    }

   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
   printf("$\\chi^2=%f+-%f$\n",chi2m[0],chi2m[1]);
   free(rm);free(chi2m);
   
  
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);   free(fit[e]);
  
   }
   
   free(fit);  

   free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   return r;
    
} 

*/


double **fit_MK_double_chiral_FVE(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1){
   double ***y,**x,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=4;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik2_min=1, ik2_max=3;
   int nk=(ik2_max-ik2_min+1);
   int ms;

   double *mref;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nk);
   mref[0]=0.064;
   mref[1]=0.080;
   mref[2]=0.095;
   int n,count,N=2;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;
   en[1]=nk;

   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
    guess[0]=2.05478;
    guess[1]=0.113021;
    
    
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double***) malloc(sizeof(double**)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double**) malloc(sizeof(double*)*ensembles);
       for(j=0;j<ensembles;j++){
           r[i][j]=(double*) malloc(sizeof(double)*Njack);
       }
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);


   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                    y[j][ms+n*nk]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
   
for (e=0;e<ensembles;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j];
                            rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                
                
            }
       
   } 


   
   for (j=0;j<Njack;j++){
       for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                x[ms+n*nk][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0

            }
       }
 
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
  
        for(i=0;i<Npar;i++){
            r[i][e][j]=tmp[i];
        }                
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   
} 

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   im=mass_index[0][1][0];
   //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
   for (e=0;e<ensembles;e++)
   {im=mass_index[e][1+0][0];
   //printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im][Njack-1]* gJ[e].w0[Njack-1],2) ,r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1],r[2][e][Njack-1]+mref[0]*r[3][e][Njack-1] );
   
   }
    
///////////////////////////////////////////////////////////compute MK at physical point
   en[0]=ensembles;
   en[1]=ensembles;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Nvar=10;
   Npar=5;
   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
//guess[0]=2.040670;  guess[1]=0.428773;  guess[2]=0.410534;  guess[3]=0.126490;   guess[4]=-1.550172;   guess[5]=-0.026200;
  //guess[0]=1.151539;  guess[1]=0.095508;  guess[2]=0.120769;   guess[3]=-2.1775;   guess[4]=0.232919;
guess[0]=2.074829 ;guess[1]=1.636190 ;  guess[2]=0.485904;  guess[3]=0.121129;  guess[4]=-2.204862;

double *xphys=(double*) malloc(sizeof(double)*(Nvar));
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   double *rm1=(double*) malloc(sizeof(double)*(Njack));
   double **y1=(double**) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   y=(double***) malloc(sizeof(double**)*Njack);
   MK=(double**) malloc(sizeof(double*)*(nk*N));
   for(i=0;i<nk*N;i++){
         MK[i]=(double*) malloc(sizeof(double)*Njack);
   }
   

  for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }
    double KM,Kf;
for (ms=0;ms<nk;ms++){
    if(ms==0)   guess[0]=2.111834 ;guess[1]=0.331624 ;  guess[2]=0.275526;  guess[3]=0.125696;  guess[4]=-1.610862; 
    if(ms==1)   guess[0]=2.111834 ;guess[1]=0.802615 ;  guess[2]=0.256295;  guess[3]=0.128061;  guess[4]=-1.639810; 
    if(ms==2)   guess[0]=2.111834 ;guess[1]=1.109423 ;  guess[2]=0.244229;  guess[3]=0.130214;  guess[4]=-1.664364; 

    count=0;
    for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        FVE_K( r1->Bw[j], r1->fw[j], double(head[e].l1)/gJ[e].w0[j],  head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j],  mref[ms] ,gJ[e].M_PS_jack[0][j]*gJ[e].M_PS_jack[0][j],  gJ[e].f_PS_jack[0][j],r[0][e][j]+mref[ms]*r[1][e][j], r[2][e][j]+mref[ms]*r[3][e][j],&KM, &Kf);
                        rm[j]=(r[0][e][j]+mref[ms]*r[1][e][j])  -  KM*KM*r1->Bw[j]*( head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j]+mref[ms]);//-P0_w*( mlw+msw)
                        rm1[j]=(r[0][e][j]+mref[ms]*r[1][e][j])/(KM*KM);
                    }
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[2][e][j]+mref[ms]*r[3][e][j];
                        rm1[j]=rm[j]/Kf;
                    }
                    
                }
                fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                y1[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm1);
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                
                
                
                
        }
        count+=en[n];
    }   
    int ii;
    //#pragma omp parallel for  private(tmp,e,i,xphys,n,count)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)
    for (j=0;j<Njack;j++){
        count=0;
        for (n=0;n<N;n++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                x[e+count][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                x[e+count][1]=mref[ms];//ms*w0
                
                x[e+count][2]=gJ[e].w0[j];//w0
                
                x[e+count][3]=gJ[e].M_PS_jack[0][j]*gJ[e].M_PS_jack[0][j];//MPi^2
                x[e+count][4]=gJ[e].f_PS_jack[0][j];//f_Pi
                x[e+count][5]=r[0][e][j]+mref[ms]*r[1][e][j];//MKw2
                x[e+count][6]=r[2][e][j]+mref[ms]*r[3][e][j];//fkw
            
                x[e+count][7]=double(head[e].l1)/gJ[e].w0[j];//f
                x[e+count][8]=r1->Bw[j];
                x[e+count][9]=r1->fw[j];
                
                for (ii=0;ii<Nvar;ii++)
                    xphys[ii]=x[e+count][ii];
            }
            count+=en[n];
        }

       // tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, MK_chiral_FVE ,guess );
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, MK_chiral_FVE  );
        
        tmp=linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, MK_linear_chiral_FVE );
        chi2[j]=compute_chi_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, MK_linear_chiral_FVE  );
        tmp[3]=tmp[3]/tmp[2];
        tmp[4]=tmp[4]/tmp[2];
        
         if(j==Njack-1){
            printf("\n\n");
            printf("P0_w=%f ;P1_w=%f ;  P3ww=%f;  Pf1w=%f;  Pf2w=%f; Pf4www=%f;  msw=%f;\n",r1->Bw[Njack-1],tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],mref[ms]  );
            printf("chi2=%f\n",chi2[j]);
            
              printf("#mlw    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
              for (e=0;e<ensembles;e++){
                  im=mass_index[e][ik2][ik1];
                  FVE_K( r1->Bw[j], r1->fw[j], double(head[e].l1)/gJ[e].w0[j],  head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j],  mref[ms] ,gJ[e].M_PS_jack[0][j]*gJ[e].M_PS_jack[0][j],  gJ[e].f_PS_jack[0][j],r[0][e][j]+mref[ms]*r[1][e][j], r[2][e][j]+mref[ms]*r[3][e][j],&KM, &Kf);

                  printf("%f    %f   %f     %f       %f       %f  %f\n",x[e][0],y1[e][0],y1[e][1],y1[e+en[0]][0],y1[e+en[0]][1],KM,Kf);
              }
         }

       // printf("guess[0]=%f;  guess[1]=%f;  guess[2]=%f;  guess[3]=%f;   guess[4]=%f;   guess[5]=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]);
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, line  );
         

         xphys[0]=r1->mlw[j];
         MK[ms][j]=MK_phys_point(0,  Nvar, xphys,Npar,tmp);//MK2
       //     printf("MKw02=%f\n",MK[ms][j]);
         MK[ms+1*nk][j]=MK_phys_point(1,  Nvar, xphys,Npar,tmp);//fK
                       
         for (i=0;i<Npar;i++){
              r1->PMK[ms][j][i]=tmp[i];
         }
           

       
         free(tmp);
    } 
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);  free(fit[e]); free(y1[e]);
   }
}
free(y1);free(rm1);

   printf("MKw2(ms1)=%f    MKw2(ms2)=%f     MKw2(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
   printf("fKw(ms1)=%f     fKw(ms2)=%f      fKw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);

 free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);

 ////////////////////////////////////////////////last interpolation  
   en[0]=nk;
   en[1]=nk;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Npar=4;
   Nvar=1;//m_l, w0,M_PS^2,f_PS

   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
   guess[0]=1;guess[1]=1;
   x=(double**) malloc(sizeof(double*)*(en_tot));


   fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   y=(double***) malloc(sizeof(double**)*Njack);
    for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                    y[j][ms+n*nk]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
   out=(double**) malloc(sizeof(double*)*2);
   out[0]=(double*) malloc(sizeof(double)*Njack);
   out[1]=(double*) malloc(sizeof(double)*Njack);

   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms][j];
                          //  rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms+1*nk][j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=mref[ms];//ml*w0
                
            }
       
   } 


   double in;
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines,guess ).P;
      //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        in=r1->MKMeV[j]*r1->w0MeV[j]*r1->MKMeV[j]*r1->w0MeV[j];
        out[0][j]=(in-tmp[0])/tmp[1];
        out[1][j]=tmp[2]+tmp[3]*out[0][j];
                 
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   

   
   
   free(rm);   
   return out;
    
} 













////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double **fit_MK_double_chiral_FVE_P40(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1){
   double ***y,**x,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=4;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik2_min=1, ik2_max=3;
   int nk=(ik2_max-ik2_min+1);
   int ms;

   double *mref;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nk);
   mref[0]=0.064;
   mref[1]=0.080;
   mref[2]=0.095;
   int n,count,N=2;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;
   en[1]=nk;

   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
    guess[0]=2.05478;
    guess[1]=0.113021;
    
   
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double***) malloc(sizeof(double**)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double**) malloc(sizeof(double*)*ensembles);
       for(j=0;j<ensembles;j++){
           r[i][j]=(double*) malloc(sizeof(double)*Njack);
       }
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);


   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                    y[j][ms+n*nk]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
   
for (e=0;e<ensembles;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j];
                            rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                
            }
       
   } 


   
   for (j=0;j<Njack;j++){
       for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                x[ms+n*nk][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
            }
       }
 
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
          
        for(i=0;i<Npar;i++){
            r[i][e][j]=tmp[i];
        }                
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   
}   

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   im=mass_index[0][1][0];
  // printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
   for (e=0;e<ensembles;e++)
   {im=mass_index[e][1+0][0];
   //printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im][Njack-1]* gJ[e].w0[Njack-1],2) ,r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1],r[2][e][Njack-1]+mref[0]*r[3][e][Njack-1] );
   
   }
    
///////////////////////////////////////////////////////////compute MK at physical point
   en[0]=ensembles;
   en[1]=ensembles;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Nvar=10;
   Npar=4;
   guess=(double*) malloc(sizeof(double*)*6);
   for(i=0;i<6;i++)
       guess[i]=1.;
//guess[0]=2.040670;  guess[1]=0.428773;  guess[2]=0.410534;  guess[3]=0.126490;   guess[4]=-1.550172;  
  //guess[0]=1.151539;  guess[1]=0.095508;  guess[2]=0.129585;   guess[3]=-1.802169;   
guess[0]=2.074829 ;guess[1]=1.636190 ;  guess[2]=0.485904;  guess[3]=0.130146;  guess[4]=-1.546381;
double *xphys=(double*) malloc(sizeof(double)*(Nvar));
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   y=(double***) malloc(sizeof(double**)*Njack);
   MK=(double**) malloc(sizeof(double*)*(nk*N));
   for(i=0;i<nk*N;i++){
         MK[i]=(double*) malloc(sizeof(double)*Njack);
   }
   

  for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }
double KM,Kf;
 int ii;
for (ms=0;ms<nk;ms++){
   if(ms==0)   guess[0]=2.111834 ;guess[1]=0.331624 ;  guess[2]=0.275526;  guess[3]=0.125696;  guess[4]=-1.610862; 
   if(ms==1)   guess[0]=2.111834 ;guess[1]=0.802615 ;  guess[2]=0.256295;  guess[3]=0.128061;  guess[4]=-1.639810; 
   if(ms==2)   guess[0]=2.111834 ;guess[1]=1.109423 ;  guess[2]=0.244229;  guess[3]=0.130214;  guess[4]=-1.664364; 
    
   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                         FVE_K( r1->Bw[j], r1->fw[j], double(head[e].l1)/gJ[e].w0[j],  head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j],  mref[ms] ,gJ[e].M_PS_jack[0][j]*gJ[e].M_PS_jack[0][j],  gJ[e].f_PS_jack[0][j],r[0][e][j]+mref[ms]*r[1][e][j], r[2][e][j]+mref[ms]*r[3][e][j],&KM, &Kf);
                        rm[j]=(r[0][e][j]+mref[ms]*r[1][e][j])  -  KM*KM*r1->Bw[j]*( head[e].k[head[e].nk+ik1]* gJ[e].w0[j]/ gJ[e].Zp[j]+ mref[ms]);//-P0_w*( mlw+msw)
                        //rm[j]=r[0][e][j]+mref[ms]*r[1][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[2][e][j]+mref[ms]*r[3][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                
        }
        count+=en[n];
   }    
    for (j=0;j<Njack;j++){
       // tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, MK_chiral_FVE_P40 ,guess );
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, MK_chiral_FVE_P40  );
       count=0;
       for (n=0;n<N;n++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                x[e+count][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                x[e+count][1]=mref[ms];//ms*w0
                
                x[e+count][2]=gJ[e].w0[j];//w0
                
                x[e+count][3]=gJ[e].M_PS_jack[0][j]*gJ[e].M_PS_jack[0][j];//MPi^2
                x[e+count][4]=gJ[e].f_PS_jack[0][j];//f_Pi
                x[e+count][5]=r[0][e][j]+mref[ms]*r[1][e][j];//MKw2
                x[e+count][6]=r[2][e][j]+mref[ms]*r[3][e][j];//fkw
            
                x[e+count][7]=double(head[e].l1)/gJ[e].w0[j];//f
                x[e+count][8]=r1->Bw[j];
                x[e+count][9]=r1->fw[j];
                
                for (ii=0;ii<Nvar;ii++)
                    xphys[ii]=x[e+count][ii];      
            }
            count+=en[n];
       }    

        tmp=linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, MK_linear_chiral_FVE_P40 );
        chi2[j]=compute_chi_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, MK_linear_chiral_FVE_P40  );
        tmp[3]=tmp[3]/tmp[2];
      
         if(j==Njack-1){
            printf("\n\n");
            printf("P0_w=%f ;P1_w=%f ;  P3ww=%f;  Pf1w=%f;  Pf2w=%f; Pf4www=%f;  msw=%f;\n",r1->Bw[Njack-1],tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],mref[ms]  );
            printf("chi2=%f\n",chi2[j]);

              printf("#mlw    MKw2  errr     fkw0    err\n");
              for (e=0;e<ensembles;e++){
                  printf("%f    %f   %f     %f       %f\n",x[e][0],y[Njack-1][e][0],y[Njack-1][e][1],y[Njack-1][e+en[0]][0],y[Njack-1][e+en[0]][1]);
                  
            }
         }

       // printf("guess[0]=%f;  guess[1]=%f;  guess[2]=%f;  guess[3]=%f;   guess[4]=%f;   guess[5]=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]);
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, line  );
          xphys[0]=r1->mlw[j];
            MK[ms][j]=MK_phys_point(0,  Nvar, xphys,Npar,tmp);//MK2
       //     printf("MKw02=%f\n",MK[ms][j]);
            MK[ms+1*nk][j]=MK_phys_point(1,  Nvar, xphys,Npar,tmp);//fK
                       
       

       
         free(tmp);
   } 
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);  free(fit[e]);
   }
}

   printf("MKw2(ms1)=%f    MKw2(ms2)=%f     MKw2(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
   printf("fKw(ms1)=%f     fKw(ms2)=%f      fKw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);

 free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);

 ////////////////////////////////////////////////last interpolation  
   en[0]=nk;
   en[1]=nk;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Npar=4;
   Nvar=1;//m_l, w0,M_PS^2,f_PS

   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
   guess[0]=1;guess[1]=1;
   x=(double**) malloc(sizeof(double*)*(en_tot));


   fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   y=(double***) malloc(sizeof(double**)*Njack);
    for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                    y[j][ms+n*nk]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
   out=(double**) malloc(sizeof(double*)*2);
   out[0]=(double*) malloc(sizeof(double)*Njack);
   out[1]=(double*) malloc(sizeof(double)*Njack);

   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms][j];
                          //  rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms+1*nk][j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=mref[ms];//ml*w0
                
            }
       
   } 


  double in;
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines,guess ).P;
      //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        in=r1->MKMeV[j]*r1->w0MeV[j]*r1->MKMeV[j]*r1->w0MeV[j];
        out[0][j]=(in-tmp[0])/tmp[1];
        out[1][j]=tmp[2]+tmp[3]*out[0][j];
                 
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   
   
   free(rm);   
   return out;
    
} 


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




double **fit_MD_double_chiral(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack r1){
   double ***y,**x,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=4;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik2_min=4, ik2_max=6;
   int nk=(ik2_max-ik2_min+1);
   int ms;

   double *mref;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nk);
   mref[0]=0.74;
   mref[1]=0.84;
   mref[2]=0.93;
   int n,count,N=2;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;
   en[1]=nk;

   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
   guess[0]=2.05478;
   guess[1]=0.113021;
    
   
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double***) malloc(sizeof(double**)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double**) malloc(sizeof(double*)*ensembles);
       for(j=0;j<ensembles;j++){
           r[i][j]=(double*) malloc(sizeof(double)*Njack);
       }
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);


   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                    y[j][ms+n*nk]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
   
for (e=0;e<ensembles;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].M_PS_GEVP_jack[im][j]   *  gJ[e].w0[j];
                            //rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].f_PS_ls_ss_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                
                
            }
       
   } 

   
   for (j=0;j<Njack;j++){
       for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                x[ms+n*nk][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
            }
       }
                
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
          
        for(i=0;i<Npar;i++){
            r[i][e][j]=tmp[i];
        }                
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   
} 

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   im=mass_index[0][1][0];
   //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
   for (e=0;e<ensembles;e++)
   {im=mass_index[e][1+0][0];
   //printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im][Njack-1]* gJ[e].w0[Njack-1],2) ,r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1],r[2][e][Njack-1]+mref[0]*r[3][e][Njack-1] );
   
   }
 printf("HERE1\n");
   
///////////////////////////////////////////////////////////compute MK at physical point
   en[0]=ensembles;
   en[1]=ensembles;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Nvar=3;
   Npar=8;
   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;

  guess[0]=-0.403604;  guess[1]=0.518297;  guess[2]=0.120295;   guess[3]=-2.183413;  
  guess[4]=-0.403604;  guess[5]=0.518297;  guess[6]=0.120295;   guess[7]=-2.183413;   
double *xphys=(double*) malloc(sizeof(double)*(Nvar));
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   y=(double***) malloc(sizeof(double**)*Njack);
   MK=(double**) malloc(sizeof(double*)*(nk*N));
   for(i=0;i<nk*N;i++){
         MK[i]=(double*) malloc(sizeof(double)*Njack);
   }
   

  for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }
    
for (ms=0;ms<nk;ms++){
        
   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[0][e][j]+mref[ms]*r[1][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[2][e][j]+mref[ms]*r[3][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
               
                
               // x[e+count][3]=gJ[e].M_PS_jack[0][Njack-1]*gJ[e].M_PS_jack[0][Njack-1];//MPi^2
                //x[e+count][4]=gJ[e].f_PS_jack[0][Njack-1];//f_Pi
                //x[e+count][5]=r[0][e][Njack-1]+mref[ms]*r[1][e][Njack-1];//MKw2
                //x[e+count][6]=r[2][e][Njack-1]+mref[ms]*r[3][e][Njack-1];//fkw
            
                //x[e+count][7]=double(head[e].l1)/gJ[e].w0[Njack-1];//f
                //x[e+count][8]=r1.Bw[Njack-1];
                //x[e+count][9]=r1.fw[Njack-1];
                
                
        }
        count+=en[n];
   }    
   int ii;
    for (j=0;j<Njack;j++){
        count=0;
        for (n=0;n<N;n++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                x[e+count][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                x[e+count][1]=mref[ms];//ms*w0
                
                x[e+count][2]=gJ[e].w0[j];//w0    
                
                for (ii=0;ii<Nvar;ii++)
                    xphys[ii]=x[e+count][ii];
              }
            count+=en[n];
        }    
        
        
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, MD_chiral ,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, MD_chiral  );
    
         if(j==Njack-1){
            printf("\n\n");
            printf("P0=%f ;P1w=%f ;  P2ww=%f;  P3ww=%f;  P0f=%f ;P1fw=%f ;  P2fww=%f;  P3fww=%f;  mref=%.4f\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],mref[ms]  );
            printf("chi2=%f\n",chi2[j]);

              printf("#mlw    MDw  errr     fDw0    err\n");
              for (e=0;e<ensembles;e++){
                  printf("%f    %f   %f     %f       %f\n",x[e][0],y[Njack-1][e][0],y[Njack-1][e][1],y[Njack-1][e+en[0]][0],y[Njack-1][e+en[0]][1]);
                  
            }
         }

       // printf("guess[0]=%f;  guess[1]=%f;  guess[2]=%f;  guess[3]=%f;   guess[4]=%f;   guess[5]=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]);
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, line  );
            xphys[0]=r1.mlw[j];
            MK[ms][j]=MD_phys_point(0,  Nvar, xphys,Npar,tmp);//MK2
       //     printf("MKw02=%f\n",MK[ms][j]);
            MK[ms+1*nk][j]=MD_phys_point(1,  Nvar, xphys,Npar,tmp);//fK
                       
       

       
         free(tmp);
   } 
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);  free(fit[e]);
   }
}

   printf("MDw(ms1)=%f     MDw(ms2)=%f      MDw(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
   printf("fDw(ms1)=%f     fDw(ms2)=%f      fDw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);

 free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
printf("HERE2\n");

 ////////////////////////////////////////////////last interpolation  
   en[0]=nk;
   en[1]=nk;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Npar=4;
   Nvar=1;//m_l, w0,M_PS^2,f_PS

   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
   guess[0]=1;guess[1]=1;
   x=(double**) malloc(sizeof(double*)*(en_tot));


   fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   y=(double***) malloc(sizeof(double**)*Njack);
    for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                    y[j][ms+n*nk]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
   out=(double**) malloc(sizeof(double*)*2);
   out[0]=(double*) malloc(sizeof(double)*Njack);
   out[1]=(double*) malloc(sizeof(double)*Njack);

   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms][j];
                          //  rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms+1*nk][j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=mref[ms];//ml*w0
                
            }
       
   } 


   double in;
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines,guess ).P;
      //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        in=result.MDMeV[j]*result.w0MeV[j];
        out[0][j]=(in-tmp[0])/tmp[1];
        out[1][j]=tmp[2]+tmp[3]*out[0][j];
                 
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   

   
   
   free(rm);   
   return out;
    
} 


//////////////////////////////////////////////////////////////////////////////////



double **fit_MD_double_chiral_P30(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack r1){
   double ***y,**x,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=4;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik2_min=4, ik2_max=6;
   int nk=(ik2_max-ik2_min+1);
   int ms;

   double *mref;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nk);
   mref[0]=0.74;
   mref[1]=0.84;
   mref[2]=0.93;
   int n,count,N=2;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;
   en[1]=nk;

   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
    guess[0]=2.05478;
    guess[1]=0.113021;
    
   
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double***) malloc(sizeof(double**)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double**) malloc(sizeof(double*)*ensembles);
       for(j=0;j<ensembles;j++){
           r[i][j]=(double*) malloc(sizeof(double)*Njack);
       }
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);


   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                    y[j][ms+n*nk]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
   
for (e=0;e<ensembles;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].M_PS_GEVP_jack[im][j]   *  gJ[e].w0[j];
                            //rm[j]*=rm[j];
                           // printf("MDw0=%f\n",rm[j]);
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].f_PS_ls_ss_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                
            }
       
   } 

   
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
          
        for(i=0;i<Npar;i++){
            r[i][e][j]=tmp[i];
        }                
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   
} 

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   im=mass_index[0][1][0];
   //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
   for (e=0;e<ensembles;e++)
   {im=mass_index[e][0][0];
   //printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im][Njack-1]* gJ[e].w0[Njack-1],2) ,r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1],r[2][e][Njack-1]+mref[0]*r[3][e][Njack-1] );
   
   }
 
///////////////////////////////////////////////////////////compute MK at physical point
   en[0]=ensembles;
   en[1]=ensembles;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Nvar=3;
   Npar=7;
   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;

  guess[0]=-0.403604;  guess[1]=0.518297;  guess[2]=0.120295;   guess[3]=-2.183413;  
  guess[4]=-0.403604;  guess[5]=0.518297;  guess[6]=0.120295;      
double *xphys=(double*) malloc(sizeof(double)*(Nvar));
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   y=(double***) malloc(sizeof(double**)*Njack);
   MK=(double**) malloc(sizeof(double*)*(nk*N));
   for(i=0;i<nk*N;i++){
         MK[i]=(double*) malloc(sizeof(double)*Njack);
   }
   

  for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }
    
for (ms=0;ms<nk;ms++){
        
   count=0;
   for (n=0;n<N;n++){
       // printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[0][e][j]+mref[ms]*r[1][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[2][e][j]+mref[ms]*r[3][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=mref[ms];//ms*w0
                
                x[e+count][2]=gJ[e].w0[Njack-1];//w0
                
               // x[e+count][3]=gJ[e].M_PS_jack[0][Njack-1]*gJ[e].M_PS_jack[0][Njack-1];//MPi^2
                //x[e+count][4]=gJ[e].f_PS_jack[0][Njack-1];//f_Pi
                //x[e+count][5]=r[0][e][Njack-1]+mref[ms]*r[1][e][Njack-1];//MKw2
                //x[e+count][6]=r[2][e][Njack-1]+mref[ms]*r[3][e][Njack-1];//fkw
            
                //x[e+count][7]=double(head[e].l1)/gJ[e].w0[Njack-1];//f
                //x[e+count][8]=r1.Bw[Njack-1];
                //x[e+count][9]=r1.fw[Njack-1];
                
                for (j=0;j<Nvar;j++)
                    xphys[j]=x[e+count][j];
                
                
                
        }
        count+=en[n];
   }    
    for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, MD_chiral_P30 ,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, MD_chiral_P30  );
    
         if(j==Njack-1){
            printf("\n\n");
            printf("P0=%f ;P1w=%f ;  P2ww=%f;  P3ww=%f;  P0f=%f ;P1fw=%f ;  P2fww=%f;  P3fww=0;  mref=%.4f\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],mref[ms]  );
            printf("chi2=%f\n",chi2[j]);

              printf("#mlw    MDw  errr     fDw0    err\n");
              for (e=0;e<ensembles;e++){
                  printf("%f    %f   %f     %f       %f\n",x[e][0],y[Njack-1][e][0],y[Njack-1][e][1],y[Njack-1][e+en[0]][0],y[Njack-1][e+en[0]][1]);
                  
            }
         }

       // printf("guess[0]=%f;  guess[1]=%f;  guess[2]=%f;  guess[3]=%f;   guess[4]=%f;   guess[5]=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]);
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, line  );
            xphys[0]=r1.mlw[j];
            MK[ms][j]=MD_phys_point(0,  Nvar, xphys,Npar,tmp);//MK2
       //     printf("MKw02=%f\n",MK[ms][j]);
            MK[ms+1*nk][j]=MD_phys_point(1,  Nvar, xphys,Npar,tmp);//fK
                       
       

       
         free(tmp);
   } 
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);  free(fit[e]);
   }
}

   printf("MDw(ms1)=%f     MDw(ms2)=%f      MDw(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
   printf("fDw(ms1)=%f     fDw(ms2)=%f      fDw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);

 free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
printf("HERE2\n");

 ////////////////////////////////////////////////last interpolation  
   en[0]=nk;
   en[1]=nk;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Npar=4;
   Nvar=1;//m_l, w0,M_PS^2,f_PS

   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
   guess[0]=1;guess[1]=1;
   x=(double**) malloc(sizeof(double*)*(en_tot));


   fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   y=(double***) malloc(sizeof(double**)*Njack);
    for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                    y[j][ms+n*nk]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
   out=(double**) malloc(sizeof(double*)*2);
   out[0]=(double*) malloc(sizeof(double)*Njack);
   out[1]=(double*) malloc(sizeof(double)*Njack);

   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms][j];
                          //  rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms+1*nk][j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=mref[ms];//ml*w0
                
            }
       
   } 


   double in;
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines,guess ).P;
      //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        in=result.MDMeV[j]*result.w0MeV[j];
        out[0][j]=(in-tmp[0])/tmp[1];
        out[1][j]=tmp[2]+tmp[3]*out[0][j];
                 
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   

   
   
   free(rm);   
   return out;
    
} 



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
/////////////////////////M_DS
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////




double **fit_MDs_double_chiral(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack r1){
   double ***y,**x,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=4;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik1_min=1, ik1_max=3;
   int ik2_min=4, ik2_max=6;
   int nk=(ik2_max-ik2_min+1);
   int ms;

   double ****MDs_mc;
   double *mref;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nk);
   mref[0]=0.064;
   mref[1]=0.080;
   mref[2]=0.095;
   int n,count,N=2;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;
   en[1]=nk;

   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
    guess[0]=2.05478;
    guess[1]=0.113021;
    
   MDs_mc=(double****) malloc(sizeof(double***)*ensembles); 
   for (e=0;e<ensembles;e++){
        MDs_mc[e]=(double***) malloc(sizeof(double**)*nk);
        for (ms=0;ms<nk;ms++){
            MDs_mc[e][ms]=(double**) malloc(sizeof(double*)*N);
            for (i=0;i<N;i++)
                MDs_mc[e][ms][i]=(double*) malloc(sizeof(double)*Njack);
            
        }
   }
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double***) malloc(sizeof(double**)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double**) malloc(sizeof(double*)*ensembles);
       for(j=0;j<ensembles;j++){
           r[i][j]=(double*) malloc(sizeof(double)*Njack);
       }
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

printf("HERE0\n");
   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                    y[j][ms+n*nk]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
for (ik2=ik2_min;ik2<=ik2_max;ik2++){
for (e=0;e<ensembles;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ik2][ms+ik1_min];
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].M_PS_GEVP_jack[im][j]   *  gJ[e].w0[j];
                            rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].f_PS_ls_ss_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=head[e].k[head[e].nk+ms+ik1_min]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                
            }
       
   }


   
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
          
           
        MDs_mc[e][ik2-ik2_min][0][j]=sqrt(tmp[0]+r1.msw[j]*tmp[1]);
        MDs_mc[e][ik2-ik2_min][1][j]=tmp[2]+r1.msw[j]*tmp[3];
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   
} 
}  
  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   
 
   
   free(y); free(guess);
   free(mref);free(en);
   im=mass_index[0][1][0];
   //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
   for (e=0;e<ensembles;e++)
   {im=mass_index[e][1+0][0];
   //printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im][Njack-1]* gJ[e].w0[Njack-1],2) ,r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1],r[2][e][Njack-1]+mref[0]*r[3][e][Njack-1] );
   
   }

 /////////////interpolation m_c
   
   
   
   mref=(double*) malloc(sizeof(double)*nk);
   mref[0]=0.74;
   mref[1]=0.84;
   mref[2]=0.93;
   N=2;
   en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;
   en[1]=nk;

    en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
    guess[0]=2.05478;
    guess[1]=0.113021;
    
   
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double***) malloc(sizeof(double**)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double**) malloc(sizeof(double*)*ensembles);
       for(j=0;j<ensembles;j++){
           r[i][j]=(double*) malloc(sizeof(double)*Njack);
       }
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);


   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                    y[j][ms+n*nk]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
   
for (e=0;e<ensembles;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=MDs_mc[e][ms][0][j];//gJ[e].M_PS_GEVP_jack[im][j]   *  gJ[e].w0[j];
                            //rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]=MDs_mc[e][ms][1][j];//gJ[e].f_PS_ls_ss_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                
            }
       
   } 

   
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
          
        for(i=0;i<Npar;i++){
            r[i][e][j]=tmp[i];
        }                
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   
} 

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   im=mass_index[0][1][0];
   //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
   for (e=0;e<ensembles;e++)
   {im=mass_index[e][1+0][0];
   //printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im][Njack-1]* gJ[e].w0[Njack-1],2) ,r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1],r[2][e][Njack-1]+mref[0]*r[3][e][Njack-1] );
   
   }

///////////////////////////////////////////////////////////compute MK at physical point
   en[0]=ensembles;
   en[1]=ensembles;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Nvar=3;
   Npar=8;
   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;

  guess[0]=-0.403604;  guess[1]=0.518297;  guess[2]=0.120295;   guess[3]=-2.183413;  
  guess[4]=-0.403604;  guess[5]=0.518297;  guess[6]=0.120295;   guess[7]=-2.183413;   
double *xphys=(double*) malloc(sizeof(double)*(Nvar));
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   y=(double***) malloc(sizeof(double**)*Njack);
   MK=(double**) malloc(sizeof(double*)*(nk*N));
   for(i=0;i<nk*N;i++){
         MK[i]=(double*) malloc(sizeof(double)*Njack);
   }
   

  for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
   }
    
for (ms=0;ms<nk;ms++){
        
   count=0;
   for (n=0;n<N;n++){
       // printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[0][e][j]+mref[ms]*r[1][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[2][e][j]+mref[ms]*r[3][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=mref[ms];//ms*w0
                
                x[e+count][2]=gJ[e].w0[Njack-1];//w0
                
               
                for (j=0;j<Nvar;j++)
                    xphys[j]=x[e+count][j];
                
                
        }
        count+=en[n];
   }    
    for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, MD_chiral ,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, MD_chiral  );
    
         if(j==Njack-1){
            printf("\n\n");
            printf("P0=%f ;P1w=%f ;  P2ww=%f;  P3ww=%f;  P0f=%f ;P1fw=%f ;  P2fww=%f;  P3fww=%f;  mref=%.4f\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],mref[ms]  );
            printf("chi2=%f\n",chi2[j]);

              printf("#mlw    MDw  errr     fDw0    err\n");
              for (e=0;e<ensembles;e++){
                  printf("%f    %f   %f     %f       %f\n",x[e][0],y[Njack-1][e][0],y[Njack-1][e][1],y[Njack-1][e+en[0]][0],y[Njack-1][e+en[0]][1]);
                  
            }
         }

       // printf("guess[0]=%f;  guess[1]=%f;  guess[2]=%f;  guess[3]=%f;   guess[4]=%f;   guess[5]=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]);
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, line  );
            xphys[0]=r1.mlw[j];
            MK[ms][j]=MD_phys_point(0,  Nvar, xphys,Npar,tmp);//MK2
       //     printf("MKw02=%f\n",MK[ms][j]);
            MK[ms+1*nk][j]=MD_phys_point(1,  Nvar, xphys,Npar,tmp);//fK
                       
       

       
         free(tmp);
   } 
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);  free(fit[e]);
   }
}

   printf("MDsw(ms1)=%f     MDsw(ms2)=%f      MDsw(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
   printf("fDsw(ms1)=%f     fDsw(ms2)=%f      fDsw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);

 free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);

 ////////////////////////////////////////////////last interpolation  
   en[0]=nk;
   en[1]=nk;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Npar=4;
   Nvar=1;//m_l, w0,M_PS^2,f_PS

   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
   guess[0]=1;guess[1]=1;
   x=(double**) malloc(sizeof(double*)*(en_tot));


   fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   y=(double***) malloc(sizeof(double**)*Njack);
    for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                    y[j][ms+n*nk]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
   out=(double**) malloc(sizeof(double*)*2);
   out[0]=(double*) malloc(sizeof(double)*Njack);
   out[1]=(double*) malloc(sizeof(double)*Njack);

   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms][j];
                          //  rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms+1*nk][j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=mref[ms];//ml*w0
                
            }
       
   } 


   double in;
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines,guess ).P;
      //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        in=result.MDsMeV[j]*result.w0MeV[j];
        out[0][j]=(in-tmp[0])/tmp[1];
        out[1][j]=tmp[2]+tmp[3]*out[0][j];
                 
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   

   
   
   free(rm);   
   return out;
    
} 





