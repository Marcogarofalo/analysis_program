#define KandD_form_M_C


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


double one_line(int n, int Nvar, double *x,int Npar,double  *P){
    double r;
    
    
        r=P[0]+P[1]*x[0];
    
    return r;
    
}


double fK_chiral_FVE(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi;
    double pi=3.141592653589793;
    
    double Mpiw=x[0], MKw=x[1], w0=x[2], Mpi2=x[3], fpi=x[4], frac_Lw=x[7],  Bw=x[8];
    double fw=x[9],  MK2=x[5], fK=x[6];
    
    
    //double    P0_w=Bw, P1_w=P[0], P3ww=P[1];
    double   Pf1w=P[0],  Pf2w=P[1],  Pf4www=P[2];
    
    double KM,Kf;
    
    FVE_K( Bw, fw, frac_Lw,  Mpiw*Mpiw/(2.*Bw)/*mlw*/, MKw*MKw-Mpiw*Mpiw/(2.) /*msw*/ ,Mpi2,  fpi,MK2, fK,&KM, &Kf);
    if (n==0){
        xi=Mpiw*Mpiw/(16*pi*pi*fw*fw);
        fKw=Pf1w*( 1.- (3./2.)* xi*log(xi)+Pf2w*xi+(1/(w0*w0))*Pf4www)*Kf;
    }
        
      
    return fKw;
    
}
double *fK_linear_chiral_FVE(int n, int Nvar, double *x,int Npar){
    
    double fKw=0,xi,*r;
    double pi=3.141592653589793;
    
    double Mpiw=x[0], MKw=x[1], w0=x[2], Mpi2=x[3], fpi=x[4], frac_Lw=x[7],  Bw=x[8];
    double fw=x[9],  MK2=x[5], fK=x[6];
     r=(double*) malloc(sizeof(double)*Npar);
    
    //double    P0_w=Bw, P1_w=P[0], P3ww=P[1];
    //double   Pf1w=P[0],  Pf2w=P[1],  Pf4www=P[2];
    
    double KM=1.,Kf=1.;
    
   // FVE_K( Bw, fw, frac_Lw,  Mpiw*Mpiw/(2.*Bw)/*mlw*/, MKw*MKw/Bw-Mpiw*Mpiw/(2.*Bw) /*msw*/ ,Mpi2,  fpi,MK2, fK,&KM, &Kf);
        
    //    fKw=Pf1w*( 1.- (3./2.)* xi*log(xi)+Pf2w*xi+(1/(w0*w0))*Pf4www)*Kf;
    
    xi=Mpiw*Mpiw/(16*pi*pi*fw*fw);
    r[0]=(1.- (3./2.)* xi*log(xi))*Kf;
    r[1]=xi*Kf;
    r[2]=(1/(w0*w0))*Kf;
      
    return r;
    
}


double fK_phys_point(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi;
    double pi=3.141592653589793;
    
   
    double Mpiw=x[0], MKw=x[1], w0=x[2], Mpi2=x[3], fpi=x[4], frac_Lw=x[7],  Bw=x[8];
    double fw=x[9],  MK2=x[5], fK=x[6];
    
    
    double    P0_w=Bw, P1_w=P[0], P3ww=P[1];
    double   Pf1w=P[0],  Pf2w=P[1],  Pf4www=P[2];
   
    
    double KM,Kf;
    
    //FVE_K( Bw, fw, frac_Lw,  mlw,  msw ,Mpi2,  fpi,MK2, fK,&KM, &Kf);
    if (n==0){
        xi=Mpiw*Mpiw/(16*pi*pi*fw*fw);
        fKw=Pf1w*( 1.- (3./2.)* xi*log(xi)+Pf2w*xi);
    }
        
    
    return fKw;
    
}




double fD_chiral(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fDw=0,xi;
    double pi=3.141592653589793;
    
    double Mpiw=x[0], MDw=x[1], w0=x[2];
    
    
    double     P0f=P[0], P1fw=P[1], P2fww=P[2],   P3fww=P[3];
    
    
    fDw=P0f+P1fw*Mpiw+P2fww*Mpiw*Mpiw+P3fww*(1./(w0*w0));
    
        
      
    return fDw;
    
}
double **fit_fK_double_chiral_FVE_from_M(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1){
   double ***y,**x,***r,**fK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=2;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik2_min=1, ik2_max=3;
   int nk=(ik2_max-ik2_min+1);
   int ms;

   double *mref;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nk);
   mref[0]=0.38*0.38;
   mref[1]=0.42*0.42;
   mref[2]=0.46*0.46;
   
   mref[0]=0.14;
   mref[1]=0.16;
   mref[2]=0.19;
   
   
   int n,count,N=1;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;
  // en[1]=nk;

   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
   for(i=0;i<Npar;i++)
       guess[i]=1.;
    guess[0]=1.;
    guess[1]=1.;
    
    
   x=(double**) malloc(sizeof(double*)*(en_tot));

   rm=(double*) malloc(sizeof(double*)*(Njack));
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
                            rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                //x[ms+n*nk][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
               
                
            }
       
   } 


   
   for (j=0;j<Njack;j++){
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                x[ms+n*nk][0]=gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j]*gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j];// M_K*w0
                x[ms+n*nk][0]=x[ms+n*nk][0]-   (gJ[e].M_PS_jack[0][j]   *  gJ[e].w0[j]*gJ[e].M_PS_jack[0][j]   *  gJ[e].w0[j])/2.; //M_K^2 w0^2 - M_pi^2w^2/2
            }
        }
        non_linear_fit_result single_jack_fit=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess );
        tmp= single_jack_fit.P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, one_line  );

        for(i=0;i<Npar;i++){
            r[i][e][j]=tmp[i];
        }                
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
    printf("$\\chi^2/dof=%f+-%f$\n",chi2[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
    free(chi2m);
   
} 
    
 
  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
    
///////////////////////////////////////////////////////////compute MK at physical point
  
   en[0]=ensembles;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Nvar=10;
   Npar=3;
   guess=(double*) malloc(sizeof(double*)*6);
   for(i=0;i<6;i++)
       guess[i]=1.;
  
guess[0]=2.074829 ;guess[1]=1.636190 ;  guess[2]=0.485904;

double *xphys=(double*) malloc(sizeof(double)*(Nvar));
   x=(double**) malloc(sizeof(double*)*(en_tot));

   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   y=(double***) malloc(sizeof(double**)*Njack);
   fK=(double**) malloc(sizeof(double*)*(nk*N));
   for(i=0;i<nk*N;i++){
         fK[i]=(double*) malloc(sizeof(double)*Njack);
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
    if(ms==0)   guess[0]=0.125214 ;guess[1]=-1.614364 ;  guess[2]=0.0023; 
    if(ms==1)   guess[0]=0.128214 ;guess[1]=-1.634364 ;  guess[2]=0.028130; 
    if(ms==2)   guess[0]=0.130214 ;guess[1]=-2.44364  ;  guess[2]=0.094130; 

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
    
                x[e+count][0]=gJ[e].M_PS_jack[0][j]*gJ[e].w0[j];//M_Pi*w0
                //x[e+count][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=sqrt(mref[ms]);//MK*w0
                
                x[e+count][2]=gJ[e].w0[j];//w0
                
                x[e+count][3]=gJ[e].M_PS_jack[0][j]*gJ[e].M_PS_jack[0][j];//MPi^2
                x[e+count][4]=gJ[e].f_PS_jack[0][j];//f_Pi
                x[e+count][5]=mref[ms]/(gJ[e].w0[j]*gJ[e].w0[j]);//MK2
                x[e+count][6]=r[0][e][j]+mref[ms]*r[1][e][j];//fkw
           
                x[e+count][7]=double(head[e].l1)/gJ[e].w0[j];//f
                x[e+count][8]=r1->Bw[j];
                x[e+count][9]=r1->fw_from_M[j];
                 for (ii=0;ii<Nvar;ii++)
                    xphys[ii]=x[e+count][ii];
            }
            count+=en[n];
        } 

        // tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fK_chiral_FVE ,guess );
        //chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fK_chiral_FVE  );
        tmp=linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fK_linear_chiral_FVE  );
        chi2[j]=compute_chi_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fK_linear_chiral_FVE  );
        tmp[1]=tmp[1]/tmp[0];
        tmp[2]=tmp[2]/tmp[0];

         if(j==Njack-1){
            printf("\n\n");
            double KM,Kf;
            
            printf("Pf1w=%f;  Pf2w=%f; Pf4www=%f;  MKw2-Mpi22=%f;\n",tmp[0],tmp[1],tmp[2],mref[ms]  );
            chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
            printf("$\\chi^2/dof=%f+-%f$\n",chi2[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
            free(chi2m);
              printf("#(M_Pi*w0)^2        fkw0    err  Kf  w0\n");
              for (e=0;e<ensembles;e++){
                //  FVE_K( r1->Bw[Njack-1], r1->fw[Njack-1], double(head[e].l1)/gJ[e].w0[Njack-1],  x[e][0]*x[e][0]/(2.*r1->Bw[Njack-1])/*mlw*/, (mref[ms])/r1->Bw[Njack-1]-x[e][0]*x[e][0]/(2.*r1->Bw[Njack-1]) /*msw*/ ,x[e][3],  x[e][4],x[e][5], x[e][6],&KM, &Kf);
                  Kf=1;
                  printf("%f        %f       %f    %f     %f\n",x[e][0]*x[e][0],y[Njack-1][e][0]/Kf,y[Njack-1][e][1]/Kf,Kf,gJ[e].w0[Njack-1]);
              }
         }
         

         xphys[0]=r1->MpiMeV[j]*r1->w0MeV[j];
         fK[ms][j]=fK_phys_point(0,  Nvar, xphys,Npar,tmp);//MK2

       
         free(tmp);
    } 
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);  free(fit[e]);
   }
}

   printf("fKw2(ms1)=%f    fKw(ms2)=%f     fKw(ms3)=%f\n", fK[0][Njack-1],fK[1][Njack-1],fK[2][Njack-1]);


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
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Npar=2;
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
                            rm[j]= fK[ms][j];
                          //  rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
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
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess ).P;
        in=result.MKMeV[j]*result.w0MeV[j];
        in=in*in;
        in=result.MKMeV[j]*result.w0MeV[j]*result.MKMeV[j]*result.w0MeV[j]  -result.MpiMeV[j]*result.w0MeV[j]*result.MpiMeV[j]*result.w0MeV[j]/2.   ;
        out[0][j]=tmp[0]+tmp[1]*in;
                 
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]); free(fK[ms]);
   }
   free(x); free(fit); free(fK);

  
        
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



////////////////////////////////////////////////////////////////////
////D meson
///////////////////////////////////////////////////////////////////////

double **fit_fD_chiral_continuum_from_M(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1){
   double ***y,**x,***r,**fD,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=2;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik2_min=6, ik2_max=9;
   int nk=(ik2_max-ik2_min+1);
   int ms;

   double *mref;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nk);
 /*  mref[0]=1.20;
   mref[1]=1.30;
   mref[2]=1.40;
   mref[3]=1.50;*/
   mref[0]=1.40;
   mref[1]=1.50;
   mref[2]=1.60;
   mref[3]=1.70;
   int n,count,N=1;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;
  // en[1]=nk;

   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
    guess[0]=1.;
    guess[1]=1.;
    
    
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
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
                x[ms+n*nk][0]=gJ[e].M_PS_GEVP_jack[im][j]   *  gJ[e].w0[j];// M_D*w0
            }
             
        }
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, one_line  );

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
   
    
///////////////////////////////////////////////////////////compute MK at physical point
   en[0]=ensembles;
   en[1]=ensembles;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Nvar=3;
   Npar=4;
   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
guess[0]=2.074829 ;guess[1]=1.636190 ;  guess[2]=0.485904; guess[4]=1.;

double *xphys=(double*) malloc(sizeof(double)*(Nvar));
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   y=(double***) malloc(sizeof(double**)*Njack);
   fD=(double**) malloc(sizeof(double*)*(nk*N));
   for(i=0;i<nk*N;i++){
         fD[i]=(double*) malloc(sizeof(double)*Njack);
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
    if(ms==0)   guess[0]=0.125214 ;guess[1]=-1.614364 ;  guess[2]=0.0023; 
    if(ms==1)   guess[0]=0.128214 ;guess[1]=-1.634364 ;  guess[2]=0.028130; 
    if(ms==2)   guess[0]=0.130214 ;guess[1]=-2.44364  ;  guess[2]=0.094130; 

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
                x[e+count][0]=gJ[e].M_PS_jack[0][j]*gJ[e].w0[j];//M_Pi*w0
                x[e+count][0]=x[e+count][0]*x[e+count][0];
                x[e+count][1]=mref[ms];//ms*w0
                x[e+count][2]=gJ[e].w0[j];//w0
              
                for (ii=0;ii<Nvar;ii++)
                    xphys[ii]=x[e+count][ii];
            }
            count+=en[n];
        }   

        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fD_chiral ,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fD_chiral );

         if(j==Njack-1){
            printf("\n\n");
            printf("Pf1w=%f;  Pf2w=%f;  Pf3ww=%f;  Pf4ww=%f;  MDw=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],mref[ms]  );
            chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
            printf("$\\chi^2/dof=%f+-%f$\n",chi2[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
            free(chi2m);

              printf("#w0     M_pi^2w0^2       fDw0    err\n");
              for (e=0;e<ensembles;e++){
                  printf("%f    %f        %f       %f\n",x[e][2],x[e][0],y[Njack-1][e][0],y[Njack-1][e][1]);
              }
         }

         

         xphys[0]=r1->MpiMeV[j]*r1->w0MeV[j];
         xphys[0]=xphys[0]*xphys[0];
         xphys[2]=1e+12;
         fD[ms][j]=fD_chiral(0,  Nvar, xphys,Npar,tmp);//MK2

       
         free(tmp);
    } 
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);  free(fit[e]);
   }
}

   printf("fDw(ms1)=%f    fDw(ms2)=%f     fDw(ms3)=%f\n", fD[0][Njack-1],fD[1][Njack-1],fD[2][Njack-1]);
 //  printf("fKw(ms1)=%f     fKw(ms2)=%f      fKw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);

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
   
   Npar=2;
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
                            rm[j]= fD[ms][j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
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
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess ).P;
   
        in=result.MDMeV[j]*result.w0MeV[j];
   
        out[0][j]=tmp[0]+tmp[1]*in;
                 
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]); free(fD[ms]);
   }
   free(x); free(fit); free(fD);

  
        
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////f_Ds form M
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





double **fit_fDs_chiral_continuum_from_M(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1){
   double ***y,**x,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im,imD,imK;  
   int Npar=2;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik1_min=1, ik1_max=3;
   int ik2_min=6, ik2_max=9;
   int nk=(ik1_max-ik1_min+1);
   int nk2=(ik2_max-ik2_min+1);
   int ms;

   double ****MDs_mc;
   double *mref,in;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nk);
   mref[0]=0.38*0.38;
   mref[1]=0.42*0.42;
   mref[2]=0.46*0.46;

   int n,count,N=1;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;

   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=0.1;
    guess[0]=2.05478;
    guess[1]=0.113021;
    
   MDs_mc=(double****) malloc(sizeof(double***)*ensembles); 
   for (e=0;e<ensembles;e++){
        MDs_mc[e]=(double***) malloc(sizeof(double**)*nk2);
        for (ms=0;ms<nk2;ms++){
            MDs_mc[e][ms]=(double**) malloc(sizeof(double*)*N);
            for (i=0;i<N;i++)
                MDs_mc[e][ms][i]=(double*) malloc(sizeof(double)*Njack);
            
        }
   }
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
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
for (ik2=ik2_min;ik2<=ik2_max;ik2++){
for (e=0;e<ensembles;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ik2][ms+ik1_min];
                imK=mass_index[e][ms+ik1_min][0];
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].f_PS_ls_ss_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
               

            }
       
   }


   
   for (j=0;j<Njack;j++){
          for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ik2][ms+ik1_min];
                imK=mass_index[e][ms+ik1_min][0];
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                //x[ms+n*nk][0]=gJ[e].M_PS_jack[imK][j]   *  gJ[e].w0[j];//MK w0
                //x[ms+n*nk][0]=x[ms+n*nk][0]*x[ms+n*nk][0];
                x[ms+n*nk][0]=pow(gJ[e].M_PS_jack[imK][j]   *  gJ[e].w0[j],2)  -pow(gJ[e].M_PS_jack[0][j]   *  gJ[e].w0[j],2)/2. ;//MK2 w02- M_pi2w02/2
            }       
          }


        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, one_line  );
          
           
        //MDs_mc[e][ik2-ik2_min][0][j]=sqrt(tmp[0]+r1.msw[j]*tmp[1]);
        in=r1->MKMeV[j]*r1->w0MeV[j];
        in=in*in- r1->MpiMeV[j]*r1->w0MeV[j]*r1->MpiMeV[j]*r1->w0MeV[j]/2.;
        MDs_mc[e][ik2-ik2_min][0][j]=tmp[0]+in*tmp[1];
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
   
   
   nk=(ik2_max-ik2_min+1);
   mref=(double*) malloc(sizeof(double)*nk);
   /*mref[0]=1.20;
   mref[1]=1.30;
   mref[2]=1.40;
   mref[3]=1.50;*/
   mref[0]=1.40;
   mref[1]=1.50;
   mref[2]=1.60;
   mref[3]=1.70;

   N=1;
   en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;
   //en[1]=nk;

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
   rm=(double*) malloc(sizeof(double*)*(Njack));
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
               imD=mass_index[e][ms+ik2_min][0];
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=MDs_mc[e][ms][0][j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
            /*    if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]=MDs_mc[e][ms][1][j];//gJ[e].f_PS_ls_ss_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
*/
                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
               
            }
       
   } 

   
   for (j=0;j<Njack;j++){
       for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                 imD=mass_index[e][ms+ik2_min][0];
                 x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                 x[ms+n*nk][0]=gJ[e].M_PS_GEVP_jack[imD][j] *  gJ[e].w0[j];//ml*w0
            }
       } 
  
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, one_line  );
          
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
   //en[1]=ensembles;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Nvar=3;
   Npar=4;
   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
  

  guess[0]=-0.403604;  guess[1]=0.518297;  guess[2]=0.120295;   guess[3]=-2.183413;  
 // guess[4]=-0.403604;  guess[5]=0.518297;  guess[6]=0.120295;   guess[7]=-2.183413;   
double *xphys=(double*) malloc(sizeof(double)*(Nvar));
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
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
    int ii;
for (ms=0;ms<nk;ms++){
        
   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][0][0];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[0][e][j]+mref[ms]*r[1][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                /*if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[2][e][j]+mref[ms]*r[3][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }*/
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
                
                
        }
        count+=en[n];
   }    
    for (j=0;j<Njack;j++){
        count=0;
        for (n=0;n<N;n++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][0][0];
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=gJ[e].M_PS_jack[im][j]*gJ[e].w0[j];
                x[e+count][0]=x[e+count][0]*x[e+count][0];
                x[e+count][1]=mref[ms];//ms*w0
                
                x[e+count][2]=gJ[e].w0[j];//w0
                
                for (ii=0;ii<Nvar;ii++)
                    xphys[ii]=x[e+count][ii];
            }
            count+=en[n];
        }    
  
        if (j==0) guess=guess_for_non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fD_chiral ,guess ); 
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fD_chiral ,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fD_chiral  );
    
         if(j==Njack-1){
            printf("\n\n");
            printf("  P0f=%f ;P1fw=%f ;  P2fww=%f;  P3fww=%f;  mref=%.4f\n",tmp[0],tmp[1],tmp[2],tmp[3],mref[ms]  );
            chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
            printf("$\\chi^2/dof=%f+-%f$\n",chi2[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
            free(chi2m);

              printf("#w0   w0^2M_pi^2    fDsw  errr    \n");
              for (e=0;e<ensembles;e++){
                  printf("%f   %f    %f   %f   \n",x[e][2],x[e][0],y[Njack-1][e][0],y[Njack-1][e][1]);
                  
            }
         }

       // printf("guess[0]=%f;  guess[1]=%f;  guess[2]=%f;  guess[3]=%f;   guess[4]=%f;   guess[5]=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]);
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, line  );
      
         xphys[0]=r1->MpiMeV[j]*r1->w0MeV[j];
         xphys[0]=xphys[0]*xphys[0];
         xphys[2]=1e+12;
         MK[ms][j]=fD_chiral(0,  Nvar, xphys,Npar,tmp);//MK2
       
       

       
         free(tmp);
   } 
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);  free(fit[e]);
   }
}

   printf("fDsw(ms1)=%f     fDsw(ms2)=%f      fDsw(ms3)=%f   fDsw(ms4)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1],MK[3][Njack-1]);
   printf("ms1=%f     ms2=%f      ms3=%f   ms4=%f\n", mref[0],mref[1],mref[2],mref[3]);
  // printf("fDsw(ms1)=%f     fDsw(ms2)=%f      fDsw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);

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
 //  en[1]=nk;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Npar=2;
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
               

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=mref[ms];//ml*w0
                
            }
       
   } 


   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess ).P;
      //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        in=result.MDMeV[j]*result.w0MeV[j];
        out[0][j]=tmp[0]+tmp[1]*in;
                 
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


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////w0 from MDs
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double **w0_from_MDs(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1){
   double ***y,**x,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im,imD,imK;  
   int Npar=2;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik1_min=1, ik1_max=3;
   int ik2_min=4, ik2_max=6;
   int nk=(ik2_max-ik2_min+1);
   int ms;

   double ****MDs_mc;
   double *mref,in;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nk);
   mref[0]=0.38;
   mref[1]=0.42;
   mref[2]=0.46;

   int n,count,N=1;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;

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
   rm=(double*) malloc(sizeof(double*)*(Njack));
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
for (ik2=ik2_min;ik2<=ik2_max;ik2++){
for (e=0;e<ensembles;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ik2][ms+ik1_min];
                imK=mass_index[e][ms+ik1_min][0];
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].M_PS_GEVP_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
            }
    }


   
   for (j=0;j<Njack;j++){
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ik2][ms+ik1_min];
                imK=mass_index[e][ms+ik1_min][0];
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=gJ[e].M_PS_jack[imK][j]   *  gJ[e].w0[j];
            }
        }

        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, one_line  );
          
           
        //MDs_mc[e][ik2-ik2_min][0][j]=sqrt(tmp[0]+r1.msw[j]*tmp[1]);
        in=r1->MKMeV[j]*r1->w0MeV[j];
        MDs_mc[e][ik2-ik2_min][0][j]=tmp[0]+in*tmp[1];
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
   mref[0]=1.20;
   mref[1]=1.30;
   mref[2]=1.40;
   N=1;
   en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;
   //en[1]=nk;

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
   rm=(double*) malloc(sizeof(double*)*(Njack));
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
               imD=mass_index[e][ms+ik2_min][0];
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=MDs_mc[e][ms][0][j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
            /*    if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]=MDs_mc[e][ms][1][j];//gJ[e].f_PS_ls_ss_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
*/
                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
            }
       
   } 

   
   for (j=0;j<Njack;j++){
       for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
               imD=mass_index[e][ms+ik2_min][0];
               x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
               x[ms+n*nk][0]=gJ[e].M_PS_GEVP_jack[imD][j] *  gJ[e].w0[j];//ml*w0
            }
       }
 
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, one_line  );
          
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
   //en[1]=ensembles;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Nvar=3;
   Npar=4;
   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;

  guess[0]=-0.403604;  guess[1]=0.518297;  guess[2]=0.120295;   guess[3]=-2.183413;  
 // guess[4]=-0.403604;  guess[5]=0.518297;  guess[6]=0.120295;   guess[7]=-2.183413;   
double *xphys=(double*) malloc(sizeof(double)*(Nvar));
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
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
    int ii;
for (ms=0;ms<nk;ms++){
        
   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][0][0];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[0][e][j]+mref[ms]*r[1][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                /*if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[2][e][j]+mref[ms]*r[3][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }*/
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
                
                                
                
                for (j=0;j<Nvar;j++)
                    xphys[j]=x[e+count][j];
                
                
        }
        count+=en[n];
   }    
   
   for (j=0;j<Njack;j++){
        count=0;
        for (n=0;n<N;n++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][0][0];
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=gJ[e].M_PS_jack[im][j]*gJ[e].w0[j];
                x[e+count][0]=x[e+count][0]*x[e+count][0];
                x[e+count][1]=mref[ms];//ms*w0
                
                x[e+count][2]=gJ[e].w0[j];//w0
               
                for (ii=0;ii<Nvar;ii++)
                    xphys[ii]=x[e+count][ii];
            }
            count+=en[n];
        }    

        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fD_chiral ,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fD_chiral  );
    
         if(j==Njack-1){
            printf("\n\n");
            printf("  P0f=%f ;P1fw=%f ;  P2fww=%f;  P3fww=%f;  mref=%.4f\n",tmp[0],tmp[1],tmp[2],tmp[3],mref[ms]  );
            printf("chi2=%f\n",chi2[j]);

              printf("#mlw    fDsw  errr    \n");
              for (e=0;e<ensembles;e++){
                  printf("%f    %f   %f   \n",x[e][0],y[Njack-1][e][0],y[Njack-1][e][1]);
                  
            }
         }

       // printf("guess[0]=%f;  guess[1]=%f;  guess[2]=%f;  guess[3]=%f;   guess[4]=%f;   guess[5]=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]);
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, line  );
      
         xphys[0]=r1->MpiMeV[j]*r1->w0MeV[j];
         xphys[0]=xphys[0]*xphys[0];
         xphys[2]=1e+12;
         MK[ms][j]=fD_chiral(0,  Nvar, xphys,Npar,tmp);//MK2
       
       

       
         free(tmp);
   } 
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);  free(fit[e]);
   }
}

   printf("MDsw(ms1)=%f     MDsw(ms2)=%f      MDsw(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
   printf("ms1=%f     ms2=%f      ms3=%f\n", mref[0],mref[1],mref[2]);
  // printf("fDsw(ms1)=%f     fDsw(ms2)=%f      fDsw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);

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
 //  en[1]=nk;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   Npar=2;
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
               

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=mref[ms];//ml*w0
                
            }
       
   } 


   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess ).P;
      //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        in=result.MDMeV[j]*result.w0MeV[j];
        out[0][j]=tmp[0]+tmp[1]*in;
                 
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














