#define Omega_jacob_C


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
#include "tower.hpp" 
 
#include <omp.h>

int ensemble_omega_jacob=7;



static double one_parabola(int n, int Nvar, double *x,int Npar,double  *P){
    double r;
    
    
        r=P[0]+P[1]*x[0]+P[2]*x[0]*x[0];
    
    return r;
    
}


//////////////////////////////////////////////////////////////////////////////////
//OMEGA
double chiral_continuum_M_Omega(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi,r;
    double pi=3.141592653589793;
    
    double Mpiw=x[0],MKw=x[1], w0=x[2];
    
    r=P[0]+P[1]*Mpiw*Mpiw+P[2]/(w0*w0);
    
    return r;
    
}

  


double **fit_Omegaw0_from_M(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1){
   double ***y,**x,***r,**fK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=3, Nparline=2;
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
   for (i=0;i<Npar;i++)
    guess[i]=1.;
    
    
   x=(double**) malloc(sizeof(double*)*(en_tot));

   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double***) malloc(sizeof(double**)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double**) malloc(sizeof(double*)*ensemble_omega_jacob);
       for(j=0;j<ensemble_omega_jacob;j++){
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
   double **strange_jacob=double_malloc_2(ensemble_omega_jacob,nk);
            strange_jacob[6][0]=0.0128;
            strange_jacob[6][1]=0.0161;   
            strange_jacob[6][2]=0.0193;   
            
            strange_jacob[5][0]=0.017;   
            strange_jacob[5][1]=0.0195;   
            strange_jacob[5][2]=0.0220;   
            
            strange_jacob[4][0]=0.0148;   
            strange_jacob[4][1]=0.0185;   
            strange_jacob[4][2]=0.0222;   
            
            strange_jacob[3][0]=0.0182;   
            strange_jacob[3][1]=0.0227;   
            strange_jacob[3][2]=0.0273;   
            
            strange_jacob[2][0]=0.0182;   
            strange_jacob[2][1]=0.0227;   
            strange_jacob[2][2]=0.0273;   
            
            strange_jacob[1][0]=0.0182;   
            strange_jacob[1][1]=0.0227;   
            strange_jacob[1][2]=0.0273;   
            
            strange_jacob[0][0]=0.0182;   
            strange_jacob[0][1]=0.0227;   
            strange_jacob[0][2]=0.0273;   
   
for (e=0;e<ensemble_omega_jacob;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].M_PS_jack[im][j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms+n*nk][1];
                }
                                          
                free(fit[ms+n*nk]);
                //x[ms+n*nk][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
               
            }
       
   } 

   for (j=0;j<Njack;j++){
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=head[e].k[head[e].nk+ms+ik2_min];// M_K*w0
                
            }
        }
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Nparline, one_line,guess );
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Nparline, one_line  );

        for(i=0;i<Nparline;i++){
            r[i][e][j]=tmp[i];
        }                
        free(tmp);
        for (n=0;n<N;n++)
            for (ms=0;ms<nk;ms++)
                    free(x[ms+n*nk]);
        
   }    
}
   
for (e=0;e<ensemble_omega_jacob;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].M_Omega_jack[ms][j]   *  gJ[e].w0[j];
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
                x[ms+n*nk][0]=r[0][e][j]+r[1][e][j]*strange_jacob[e][ms];
                x[ms+n*nk][0]*= gJ[e].w0[j];// M_K*w0
                x[ms+n*nk][0]=x[ms+n*nk][0]*x[ms+n*nk][0]-   (gJ[e].M_PS_jack[0][j]   *  gJ[e].w0[j]*gJ[e].M_PS_jack[0][j]   *  gJ[e].w0[j])/2.; //M_K^2 w0^2 - M_pi^2w^2/2
            }
        }
//        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess );
//        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, one_line  );

       
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_parabola,guess );
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, one_parabola  );

        for(i=0;i<Npar;i++){
            r[i][e][j]=tmp[i];
        }                
        
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        printf("%g    %g    %g\n",x[ms][0],y[Njack-1][ms][0],y[Njack-1][ms][1]);
        free(x[ms]);  free(fit[ms]);
   }
   for(i=0;i<Npar;i++)
       printf("P%d=%g;  \t ",i,r[i][e][Njack-1]);
   printf("\n");
   //printf("P0=%g;  \t P1=%g;\n",r[0][e][Njack-1],r[1][e][Njack-1]);
   
   
   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
    printf("$\\chi^2/dof=%f+-%f$\n",chi2m[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
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
   free_2(ensemble_omega_jacob,strange_jacob);
///////////////////////////////////////////////////////////compute MK at physical point
   
   en_tot=0;
   for (n=0;n<N;n++){
       en[n]=ensemble_omega_jacob;
       en_tot+=en[n];
   }
   
   Nvar=3;
   Npar=3;
   guess=(double*) malloc(sizeof(double*)*6);
  
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
double **par=double_malloc_2(Npar,Njack);

for (ms=0;ms<nk;ms++){
    if(ms==0)   guess[0]=0.125214 ;guess[1]=-1.614364 ;  
    if(ms==1)   guess[0]=0.128214 ;guess[1]=-1.634364 ;  
    if(ms==2)   guess[0]=0.130214 ;guess[1]=-2.44364  ;  

    count=0;
    for (n=0;n<N;n++){
       // printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        //rm[j]=r[0][e][j]+mref[ms]*r[1][e][j];
                        rm[j]=r[0][e][j]+mref[ms]*r[1][e][j]+mref[ms]*mref[ms]*r[2][e][j];
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
    /*            
                x[e+count][3]=gJ[e].M_PS_jack[0][j]*gJ[e].M_PS_jack[0][j];//MPi^2
                x[e+count][4]=gJ[e].f_PS_jack[0][j];//f_Pi
                x[e+count][5]=mref[ms]/(gJ[e].w0[j]*gJ[e].w0[j]);//MK2
                x[e+count][6]=r[0][e][j]+mref[ms]*r[1][e][j];//fkw
           
                x[e+count][7]=double(head[e].l1)/gJ[e].w0[j];//f
                x[e+count][8]=r1->Bw[j];
                x[e+count][9]=r1->fw_from_M[j];*/
                 for (ii=0;ii<Nvar;ii++)
                    xphys[ii]=x[e+count][ii];
            }
            count+=en[n];
        } 

        // tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fK_chiral_FVE ,guess );
        //chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fK_chiral_FVE  );
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, chiral_continuum_M_Omega,guess );
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, chiral_continuum_M_Omega  );
        

         if(j==Njack-1){
            printf("\n\n");
            double KM,Kf;
            
            printf("P0=%f;  P1=%f; P2=%g  MKw2-Mpi22=%f;\n",tmp[0],tmp[1],tmp[2],mref[ms]  );
            chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
            printf("$\\chi^2/dof=%f+-%f$\n",chi2[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
            free(chi2m);
              printf("#(M_Pi*w0)^2        M_Omega w0    err  Kf  w0\n");
              for (e=0;e<ensemble_omega_jacob;e++){
                //  FVE_K( r1->Bw[Njack-1], r1->fw[Njack-1], double(head[e].l1)/gJ[e].w0[Njack-1],  x[e][0]*x[e][0]/(2.*r1->Bw[Njack-1])/*mlw*/, (mref[ms])/r1->Bw[Njack-1]-x[e][0]*x[e][0]/(2.*r1->Bw[Njack-1]) /*msw*/ ,x[e][3],  x[e][4],x[e][5], x[e][6],&KM, &Kf);
                  Kf=1;
                  printf("%f        %f       %f    %f     %f\n",x[e][0]*x[e][0],y[Njack-1][e][0]/Kf,y[Njack-1][e][1]/Kf,Kf,gJ[e].w0[Njack-1]);
              }
         }
         for(i=0;i<Npar;i++)
             par[i][j]=tmp[i];   

       
         free(tmp);
    }
    chi2m=mean_and_error(jack_files[0].sampling,Njack, par[1]);
    if(chi2m[0]<=chi2m[1]*0.5){
        printf("slope is zero\n");
        for (j=0;j<Njack;j++)
            fK[ms][j]=par[0][j]/r1->MOmegaMeV[j] ;
    }
    else {
       for (j=0;j<Njack;j++){
            fK[ms][j]=r1->MOmegaMeV[j]  - sqrt( r1->MOmegaMeV[j]*r1->MOmegaMeV[j]-4.* par[0][j]*par[1][j]*r1->MpiMeV[j]*r1->MpiMeV[j]  );
            fK[ms][j]=fK[ms][j]/(2.*par[1][j]*r1->MpiMeV[j]*r1->MpiMeV[j] );//w0 at that reference ms mass
       }
    }
    free(chi2m);
   
   
   for (e=0;e<en_tot;e++){
        free(x[e]);  free(fit[e]);
   }
}
   
   printf("w0(ms1)=%f    w0(ms2)=%f     w0(ms3)=%f    MeV^-1\n", fK[0][Njack-1],fK[1][Njack-1],fK[2][Njack-1]);
   free_2(Npar,par);

 free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);

 ////////////////////////////////////////////////last interpolation  
   
   en_tot=0;
   for (n=0;n<N;n++){
       en[n]=nk;
       en_tot+=en[n];
   }
   
   Npar=2;
   Nvar=1;//m_l, w0,M_PS^2,f_PS

   par=double_malloc_2(Npar,Njack);
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
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess );
        //w0 at that reference ms mass
        for(i=0;i<Npar;i++)
             par[i][j]=tmp[i]; 
        free(tmp);

   }     
    chi2m=mean_and_error(jack_files[0].sampling,Njack, par[1]);
    if(chi2m[0]<=chi2m[1]*0.5){
        printf("the slope is zero");
        for (j=0;j<Njack;j++)
            out[0][j]=par[0][j];
    }
    else {
        for (j=0;j<Njack;j++){
            in=r1->MKMeV[j]*r1->MKMeV[j]-r1->MpiMeV[j]*r1->MpiMeV[j]/2.;
            out[0][j]=1.-sqrt( 1.-4.* par[0][j]*par[1][j]*in  );
            out[0][j]=out[0][j]/(2.*par[1][j]*in);
        }
    }
    free(chi2m);free_2(Npar,par);
   
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
   
   free(en);
   free(mref);
   free(xphys);
   free(rm);   
   return out;
    
} 

 
