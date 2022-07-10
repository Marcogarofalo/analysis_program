#define KandD_C


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"

#include "resampling.hpp"
#include "read.hpp"
//#include "m_eff.hpp"
#include "gnuplot.hpp"
//#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "rand.hpp"
#include "non_linear_fit.hpp"
#include "KandD.hpp" 
#include "fve.hpp" 
 
#include <omp.h>

 
double two_lines_B0(int n, int Nvar, double *x,int Npar,double  *P){
    double r;
    
    if (n==0)
        r=P[0]+P[1]*x[0];
    if (n==1)
        r=P[2]+P[3]*x[0];
    
    return r;
    
}

double *MK_linear_chiral_B0_FVE(int n, int Nvar, double *x,int Npar){
    
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
        r[0]=( mlw+msw)*KM*KM;//?????????????????????????????
        r[1]=P0_w*( mlw+msw)*mlw*KM*KM;
        r[2]=P0_w*( mlw+msw)*(1./(w0*w0))*KM*KM;
        r[3]=0;
        r[5]=0;
        r[4]=0;
        //MKw2=P0_w*( mlw+msw)*(1+P1_w *mlw+(1/(w0*w0))*P3ww)*KM*KM;
    }
    else if (n==1){
        xi=2.*Bw*mlw/(16.*pi*pi*fw*fw);
        r[0]=0;
        r[1]=0;
        r[2]=0;
        r[3]=(1.- (3./2.)* xi*log(xi))*Kf;
        r[4]=xi*Kf;
        r[5]=(1/(w0*w0))*Kf;
       // MKw2=Pf1w*( 1.- (3./2.)* xi*log(xi)+Pf2w*xi+(1/(w0*w0))*Pf4www)*Kf;
    }
      
    return r;
    
} 


double MK_phys_point_B0(int n, int Nvar, double *x,int Npar,double  *P){
    
    double MKw2=0,xi;
    double pi=3.141592653589793;
    
   
    double mlw=x[0], msw=x[1], w0=x[2], Mpi2=x[3], fpi=x[4], frac_Lw=x[7],  Bw=x[8];
    double fw=x[9],  MK2=x[5], fK=x[6];
    
    
    double    P0_w=P[0], P1_w=P[1], P3ww=P[2];
    double   Pf1w=P[3],  Pf2w=P[4],  Pf4www=P[5];
   
    
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

double **fit_MK_double_chiral_B0_FVE(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1){
   double ***y,**x,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=4;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik2_min=1, ik2_max=3;
   int nk=(ik2_max-ik2_min+1);
   int ms;
   double **B0;
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
   B0=(double**) malloc(sizeof(double*)*nk);
   for (ms=0;ms<nk;ms++)
       B0[ms]=(double*) malloc(sizeof(double)*Njack);
   
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
 
        non_linear_fit_result single_jack_fit=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines_B0,guess );
        tmp=single_jack_fit.P;
        chi2[j]=single_jack_fit.chi2;
          
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
   Npar=6;
   guess=(double*) malloc(sizeof(double*)*6);
   for(i=0;i<6;i++)
       guess[i]=1.;
  
//guess[0]=2.040670;  guess[1]=0.428773;  guess[2]=0.410534;  guess[3]=0.126490;   guess[4]=-1.550172;   guess[5]=-0.026200;
  //guess[0]=1.151539;  guess[1]=0.095508;  guess[2]=0.120769;   guess[3]=-2.1775;   guess[4]=0.232919;
guess[0]=2.074829 ;guess[1]=1.636190 ;  guess[2]=0.485904;  guess[3]=0.121129;  guess[4]=-2.204862; guess[5]=0.372566;

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
    if(ms==0)   guess[0]=2.111834 ;guess[1]=0.331624 ;  guess[2]=0.275526;  guess[3]=0.125696;  guess[4]=-1.610862; guess[5]=0.002318;
    if(ms==1)   guess[0]=2.111834 ;guess[1]=0.802615 ;  guess[2]=0.256295;  guess[3]=0.128061;  guess[4]=-1.639810; guess[5]=0.028341;
    if(ms==2)   guess[0]=2.111834 ;guess[1]=1.109423 ;  guess[2]=0.244229;  guess[3]=0.130214;  guess[4]=-1.664364; guess[5]=0.054130;

    count=0;
    for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        FVE_K( r1->Bw[j], r1->fw[j], double(head[e].l1)/gJ[e].w0[j],  head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j],  mref[ms] ,gJ[e].M_PS_jack[0][j]*gJ[e].M_PS_jack[0][j],  gJ[e].f_PS_jack[0][j],r[0][e][j]+mref[ms]*r[1][e][j], r[2][e][j]+mref[ms]*r[3][e][j],&KM, &Kf);
                        rm[j]=(r[0][e][j]+mref[ms]*r[1][e][j] );//-P0_w*( mlw+msw)
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
        
        tmp=linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, MK_linear_chiral_B0_FVE );
        chi2[j]=compute_chi_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, MK_linear_chiral_B0_FVE  );
        tmp[1]=tmp[1]/tmp[0];
        tmp[2]=tmp[2]/tmp[0];
        
        tmp[4]=tmp[4]/tmp[3];
        tmp[5]=tmp[5]/tmp[3];
         if(j==Njack-1){
            printf("\n\n");
            printf("P0_w=%f ;P1_w=%f ;  P3ww=%f;  Pf1w=%f;  Pf2w=%f; Pf4www=%f;  msw=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],mref[ms]  );
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
         MK[ms][j]=MK_phys_point_B0(0,  Nvar, xphys,Npar,tmp);//MK2
       //     printf("MKw02=%f\n",MK[ms][j]);
         MK[ms+1*nk][j]=MK_phys_point_B0(1,  Nvar, xphys,Npar,tmp);//fK
             
         B0[ms][j]=tmp[0];
        /* for (i=0;i<Npar;i++){
              r1->PMK[ms][j][i]=tmp[i];
         }*/
           

       
         free(tmp);
    } 
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);  free(fit[e]); free(y1[e]);
   }
}
free(y1);free(rm1);
double *B0m;
for (ms=0;ms<nk;ms++){
       B0m=mean_and_error(jack_files[0].sampling,Njack,B0[ms]);
       printf("B0(ms%d)=%g  +-  %2.g\t",ms,B0m[0],B0m[1]);
       free(B0[ms]);free(B0m);
}
free(B0);printf("\n");
   
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
        non_linear_fit_result single_jack_fit=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines_B0,guess );
        tmp=single_jack_fit.P;
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
