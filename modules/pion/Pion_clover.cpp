#define Pion_clover_C

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
#include "continuum_reph.hpp"
#include "fve.hpp"
#include "indices.hpp"
#include "global_reph.hpp"
#include "tower.hpp"
#include <mutils.hpp>

#include <unistd.h>

#include <omp.h> 
static void init_fit( int N, struct header *head ,int Njack, struct data_jack *gJ,int Nvar,int Npar, int **en,int *en_tot, double ****x, double ***sigmax, double **chi2m, double **rm, double ***r, double ***fit, double ****y,double **chi2,double ****C)
{
    int imoms,imomt,imom0,iG,i,n,e,j;
    int count;
   *en_tot=0;
   
   *en=(int*) calloc(N,sizeof(int));
   
   for (e=0;e<ensembles;e++){
                for (n=0;n<N;n++){
                   
                    (*en)[n]+=1;  
                }
   }
      

 
   for (n=0;n<N;n++)
   {  *en_tot+=(*en)[n];   }
   
   *x=double_malloc_3(Njack,*en_tot,Nvar);//(double**) malloc(sizeof(double*)*(*en_tot));
   *sigmax=double_malloc_2((*en)[0],Nvar);
   
   //*chi2m=(double*) malloc(sizeof(double)*(Npar));
   *rm=(double*) malloc(sizeof(double)*(Njack));

   *fit=(double**) malloc(sizeof(double*)*(*en_tot));

   *r=(double**) malloc(sizeof(double*)*(Npar));
   for(i=0;i<Npar;i++){
       (*r)[i]=(double*) malloc(sizeof(double)*Njack);
   }
   
   *chi2=(double*) malloc(sizeof(double)*Njack);
   (*y)=(double***) malloc(sizeof(double**)*Njack);
   (*C)=(double***) malloc(sizeof(double**)*Njack);
   
   for (j=0;j<Njack;j++){
       (*y)[j]=(double**) malloc(sizeof(double*)*(*en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<(*en)[n];i++){
                (*y)[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=(*en)[n];
        }
   }
}
static struct fit_result close_fit( int N, struct header *head ,int Njack, struct data_jack *gJ,int Npar,int **en,int *en_tot, double ****x, double ***sigmax, double **chi2m, double **rm,double ***r, double ***fit, double ****y,double **chi2, double ****C)
{
    int imoms,imomt,imom0,iG,i,n,e,j;
   int count;
   
   free(*chi2m);
   free(*rm);

  
   count=0;
   for (n=0;n<N;n++){
       for (i=0;i<(*en)[n];i++){
           free((*fit)[i+count]);
           //free((*x)[i+count]);
       }
       count+=(*en)[n];
   }
   free(*fit);     
   free_3(Njack,*en_tot,*x);
   free_2((*en)[0],*sigmax);

   
   for (j=0;j<Njack;j++){
        
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<(*en)[n];i++){
                free((*y)[j][i+count]);
            }
            count+=(*en)[n];
        }
        free((*y)[j]);
   }
   free( (*y));
   free(*en);
   struct fit_result fit_out;
   fit_out.Njack=Njack;
   fit_out.P=(double**) malloc(sizeof(double*)*Npar);
   fit_out.chi2=(double*) malloc(sizeof(double*)*Njack);
   fit_out.C=(double***) malloc(sizeof(double**)*Npar);
   for(i=0;i<Npar;i++){
       fit_out.P[i]=(double*) malloc(sizeof(double*)*Njack);
       for (j=0;j<Njack;j++){
           fit_out.P[i][j]=(*r)[i][j];
       }
       free((*r)[i]);
       fit_out.C[i]=(double**) malloc(sizeof(double*)*Npar);
       for(n=0;n<Npar;n++){     
           fit_out.C[i][n]=(double*) malloc(sizeof(double)*Njack);
           for (j=0;j<Njack;j++){
                fit_out.C[i][n][j]=(*C)[j][i][n];
           }
       }
   }
   for (j=0;j<Njack;j++){
       fit_out.chi2[j]=(*chi2)[j];
   }
   free(*r);
   free(*chi2);
   for (j=0;j<Njack;j++){
       for(n=0;n<Npar;n++)
           free((*C)[j][n]);
       free((*C)[j]);
   }
   free(*C);
   return fit_out;
}
 

 struct fit_result fit_Mpi_fw_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info,char **argv  ,const char *nameout){
   double ***y,***x,**sigmax,**r,*chi2,*tmp,*rm,*chi2m,**fit,***C;
   int i,j,e,im;  
   int Npar=fit_info.Npar;
   int Nvar=5;//fit_info.Nvar;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=0;
   int n,count,N=fit_info.N;
   int *en;
   char fname[NAMESIZE];
   mysprintf(fname,NAMESIZE,"%s/%s",argv[2],nameout);
   FILE *fdat=open_file(fname,"w+");
   
   int en_tot=0;
  
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
   for (i=0;i<Npar;i++)
        guess[i]=rand();

 /*   guess[0]=2.05478;
    guess[1]=0.113021;
    guess[2]=2.54451;
    guess[3]=0.410103;
    guess[4]=4.61247;
    guess[5]=-0.272884;
   */
   init_fit(N,  head , Njack, gJ,Nvar,Npar,&en,&en_tot, &x, &sigmax, &chi2m, &rm, &r, &fit, &y,&chi2,&C);

   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j];
                        rm[j]*=rm[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                    //if (e==0 || e== 1 || e==2  || e==3) y[j][e+count][1]*=100;
                    //if (e==2 ||e==3) y[j][e+count][1]*=100;
                    
                    //if (e==3)  y[j][e+count][1]*=100;
                }
                //x[e+count]=(double*) malloc(sizeof(double)*Nvar);
            
                                
                
             //   printf("%g     %g      %g   %g\n",x[e+count][1],x[e+count][0],fit[e+count][0],fit[e+count][1]);
        }
        count+=en[n];
   }

     double KM,Kf,K;
   
   
     
     
   //#pragma omp parallel for  private(tmp,i,count,n,e,im,x)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2)
   for (j=0;j<Njack;j++){
        count=0;
        //x=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                x[j][e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                x[j][e+count][1]=gJ[e].w0[j];//w0
                x[j][e+count][2]=gJ[e].M_PS_jack[im][j]*gJ[e].M_PS_jack[im][j];//MPS^2
                x[j][e+count][3]=gJ[e].f_PS_jack[im][j];//f_PS
                x[j][e+count][4]=double(head[e].l1);//f_PS
            }
            count+=en[n];
        }
   }
   
   count=0;
   for (n=0;n<1;n++){
      for (e=0;e<en[n];e++){
          for(int v=0 ;v<Nvar;v++){
              for (j=0;j<Njack;j++)
                    rm[j]=x[j][e+count][v];
              tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
             // if (fabs(tmp[1])<1e-6) {printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] ); tmp[1]=tmp[0]/1.0e+8; }
              sigmax[e+count][v]=tmp[1];
              free(tmp);
          }
      }
      count+=en[n];
   }           
         
   double **yy=double_malloc_2(en_tot,Njack);
   for ( i =0;i<en_tot;i++)
        for ( j =0;j<Njack;j++)
            yy[i][j]=y[j][i][0];
    
   double **cov=covariance(jack_files[0].sampling, en_tot, Njack, yy);  
   free_2(en_tot,yy);
   
   double **yx=double_malloc_2(en_tot+Nvar*en[0],Njack);
   for ( i =0;i<en_tot;i++)
        for ( j =0;j<Njack;j++)
            yx[i][j]=y[j][i][0];
   
   for (i=en_tot;i<en_tot+Nvar*en[0];i++){
        for ( j =0;j<Njack;j++){
            yx[i][j]=x[j][ (i-en_tot)/Nvar ][  (i-en_tot)%Nvar  ];
        }
   } 
        
   double **cov_yx=covariance(jack_files[0].sampling, en_tot+Nvar*en[0], Njack, yx);  
   free_2(en_tot+Nvar*en[0],yx);
   int yn=is_it_positive( cov_yx,  en_tot+Nvar*en[0]);
   while(yn==1){
        printf("covariance matrix not positive defined adding eps*cov[0][0]*I \n");
        for(i=0;i<en_tot+Nvar*en[0];i++)
             cov_yx[i][i]+=cov_yx[0][0]*1e-12;
        yn=is_it_positive( cov_yx,  en_tot+Nvar*en[0]);  
        printf("now the matrix is positive defined.  %d\n",yn);
   }  
   double **cov_yx1=symmetric_matrix_inverse(en_tot+Nvar*en[0], cov_yx  );
   
   yn=is_it_positive( cov,  en_tot);
   while(yn==1){
        printf("covariance matrix not positive defined adding eps*cov[0][0]*I \n");
        for(i=0;i<en_tot;i++)
             cov[i][i]+=cov[0][0]*1e-12;
        yn=is_it_positive( cov,  en_tot);  
        printf("now the matrix is positive defined.  %d\n",yn);
   }  
   double **cov1=symmetric_matrix_inverse(en_tot, cov  );
   
   guess=guess_for_non_linear_fit_Nf(N, en,x[Njack-1], y[Njack-1] , Nvar,  Npar, fit_info.function,guess );
  
   
 /*  double *guess1=(double*)  malloc(sizeof(double)*(Npar+Nvar*en[0]));
   for(i=0;i<Npar;i++)
       guess1[i]=guess[i];
   
   for(i=Npar;i<Npar+Nvar*en[0];i++)
       guess1[i]=x[0][ (i-Npar)/Nvar ][  (i-Npar)%Nvar  ];
   free(guess);
   guess=guess1;
   
   
   guess1=non_linear_fit_Nf_sigmax_covariance( N, en ,x[0], sigmax, y[0] , Nvar,  Npar,  fit_info.function , guess ,cov_yx1);
   free(guess);
   guess=guess1;
   */
   #pragma omp parallel for  private(tmp,i,count,n,e,im)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2,C,x,cov,cov_yx1,cov1)
   for (j=0;j<Njack;j++){
        //if (j==0){     }
        tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess );
        //tmp=non_linear_fit_Nf_sigmax( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_iterative( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_covariance( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess ,cov_yx1);
        //tmp=non_linear_fit_Nf_covariance(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess ,cov1);
       
        
        chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        //printf("chi2[%d]=%f",j,chi2[j]);
        //for (i=Npar;i<Npar+en[0]*Nvar;i++)
         //   printf("%.10f   ",fabs(tmp[i]-x[j][(i-Npar)/Nvar ][(i-Npar)%Nvar ] ) );
                //printf("\n");
        //printf("jacknife=%d of %d   chi2=%g\n",j, Njack,chi2[j]);

        C[j]=covariance_non_linear_fit_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function );            
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }         
        if (j==Njack-1){
            printf("w0/a[fm]     mu*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err    KM2/Kf           (Mpi^2/fpi^2)* (Kf^2/KM^2)\n"); 
            fprintf(fdat,"w0/a[fm]     mu*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err    KM2/Kf           (Mpi^2/fpi^2)* (Kf^2/KM^2)\n"); 
            count=0;
            double newline_if_w0=x[j][0][1];
            std::vector<int>  myen={0,1,2,3,   8,4,5,6,   7};
            
            for (n=0;n<N;n++){
                printf("#function %d\n",n);
                for (auto e :myen){
                    KM=fit_info.function(2,Nvar,x[j][e+count],Npar,tmp);
                    Kf=fit_info.function(3,Nvar,x[j][e+count],Npar,tmp);
                    double *tmp1=(double*) malloc(sizeof(double)*Njack);
                    for (int jj=0;jj<Njack;jj++){
                        tmp1[jj]=gJ[e].M_PS_jack[0][jj]/gJ[e].f_PS_jack[0][jj];
                        tmp1[jj]*=Kf/KM;
                        tmp1[jj]*=tmp1[jj];
                    }
                    double *tmp2=mean_and_error(jack_files[0].sampling,Njack,tmp1);
                    if (n==0)
                        K=KM*KM;
                    else if (n==1)
                        K=Kf;
                    
                    if (newline_if_w0!=x[j][e+count][1]){
                        fprintf(fdat,"\n\n");
                        newline_if_w0=x[j][e+count][1] ;
                    }
                    fprintf(fdat,"%.5f     %.5f      %.5f   %.5f     %.5f \t\t %.5f  %.5f\n",x[j][e+count][1],x[j][e+count][0],fit[e+count][0]/K,fit[e+count][1]/K,K, tmp2[0] ,tmp2[1]);
                    
                    printf("%.5f     %.5f      %.5f   %.5f     %.5f \t\t %.5f  %.5f\n",x[j][e+count][1],x[j][e+count][0],y[j][e+count][0]/K,y[j][e+count][1]/K,K, tmp2[0] ,tmp2[1]);
                    free(tmp2);free(tmp1);
                }
                count+=en[n];
            }
            
        }
        free(tmp);
        
   }  
   free_2(en_tot,cov);
   free_2(en_tot,cov1);
   free_2(en_tot+Nvar*en[0],cov_yx);
   free_2(en_tot+Nvar*en[0],cov_yx1);
   
  /* 

   printf("w0/a[fm]     m*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err    KM2/Kf\n"); 
   double KM,Kf,K;
   count=0;
   x=(double**) malloc(sizeof(double*)*(en_tot));

   for (n=0;n<N;n++){
        printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                                
                x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                x[e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=gJ[e].w0[Njack-1];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][Njack-1]*gJ[e].M_PS_jack[im][Njack-1];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][Njack-1];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
                
                //FVE( v_w0GeV, gJ[e].w0[Njack-1] ,r[2][Njack-1], r[4][Njack-1], r[0][Njack-1], r[1][Njack-1],head[e].l1, x[e+count][0] ,  x[e+count][2], x[e+count][3], &KM, &Kf);
                
                double M,L,Rf,RM;
                M=sqrt(x[e+count][2])*gJ[e].w0[Njack-1]/v_w0MeV;
                L=x[e+count][4]*v_w0fm/gJ[e].w0[Njack-1];
                Rf=exp(4.58982-0.0138032*M-2.013*L);
                RM=exp(3.81729-0.0130342*M-2.1714*L);
                Kf=1-Rf;
                KM=RM+1;
                if (n==0)
                    K=KM*KM;
                else if (n==1)
                    K=Kf;
                printf("%g     %g      %g   %g     %g \n",x[e+count][1],x[e+count][0],fit[e+count][0]/K,fit[e+count][1]/K,K);
        }
        count+=en[n];
    }
*/
   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
   printf("$\\chi^2/dof=%f+-%f$  \n",chi2m[0],chi2m[1]);
  
   printf("$\\chi^2=%f+-%f$  \n",chi2m[0]* (en_tot-Npar),chi2m[1]* (en_tot-Npar));

  
  
   struct fit_result fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x,&sigmax, &chi2m, &rm,&r, &fit, &y,&chi2,&C);
   fclose(fdat);
   
   return fit_out;
    
} 








struct fit_result fit_Mpi_fwMpi4_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info, char **argv  ,const char *nameout){
   double ***y,***x,**sigmax,**r,*chi2,*tmp,*rm,*chi2m,**fit,***C;
   int i,j,e,im;  
   int Npar=fit_info.Npar;
   int Nvar=5;//fit_info.Nvar;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=0;
   int n,count,N=fit_info.N;
   int *en;
  
   int en_tot=0;
  
   char fname[NAMESIZE];
   mysprintf(fname,NAMESIZE,"%s/%s",argv[2],nameout);
   FILE *fdat=open_file(fname,"w+");
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
   for (i=0;i<Npar;i++)
        guess[i]=rand();
   init_fit(N,  head , Njack, gJ,Nvar,Npar,&en,&en_tot, &x, &sigmax, &chi2m, &rm, &r, &fit, &y,&chi2,&C);

   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j];
                        rm[j]*=rm[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].f_PS_jack[im][j]   ;
                        rm[j]*=gJ[e].M_PS_jack[im][j]* gJ[e].M_PS_jack[im][j]* gJ[e].M_PS_jack[im][j]* gJ[e].M_PS_jack[im][j];
                        rm[j]*=  gJ[e].w0[j]*  gJ[e].w0[j]*  gJ[e].w0[j]*  gJ[e].w0[j]*  gJ[e].w0[j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                }
               
        }
        count+=en[n];
   }

     double KM,Kf,K;
   
   
     
     
   //#pragma omp parallel for  private(tmp,i,count,n,e,im,x)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2)
   for (j=0;j<Njack;j++){
        count=0;
        //x=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                x[j][e+count][0]=head[e].k[head[e].nk+ik2]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                x[j][e+count][1]=gJ[e].w0[j];//w0
                x[j][e+count][2]=gJ[e].M_PS_jack[im][j]*gJ[e].M_PS_jack[im][j];//MPS^2
                x[j][e+count][3]=gJ[e].f_PS_jack[im][j];//f_PS
                x[j][e+count][4]=double(head[e].l1);//f_PS
            }
            count+=en[n];
        }
   }
   
   count=0;
   for (n=0;n<1;n++){
      for (e=0;e<en[n];e++){
          for(int v=0 ;v<Nvar;v++){
              for (j=0;j<Njack;j++)
                    rm[j]=x[j][e+count][v];
              tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
              //if (fabs(tmp[1])<1e-6) {printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] ); tmp[1]=tmp[0]/1.0e+8; }
              sigmax[e+count][v]=tmp[1];
              free(tmp);
          }
      }
      count+=en[n];
   }           
         
   double **yy=double_malloc_2(en_tot,Njack);
   for ( i =0;i<en_tot;i++)
        for ( j =0;j<Njack;j++)
            yy[i][j]=y[j][i][0];
    
   double **cov=covariance(jack_files[0].sampling, en_tot, Njack, yy);  
   free_2(en_tot,yy);
   
   double **yx=double_malloc_2(en_tot+Nvar*en[0],Njack);
   for ( i =0;i<en_tot;i++)
        for ( j =0;j<Njack;j++)
            yx[i][j]=y[j][i][0];
   
   for (i=en_tot;i<en_tot+Nvar*en[0];i++){
        for ( j =0;j<Njack;j++){
            yx[i][j]=x[j][ (i-en_tot)/Nvar ][  (i-en_tot)%Nvar  ];
        }
   } 
        
   double **cov_yx=covariance(jack_files[0].sampling, en_tot+Nvar*en[0], Njack, yx);  
   free_2(en_tot+Nvar*en[0],yx);
   int yn=is_it_positive( cov_yx,  en_tot+Nvar*en[0]);
   while(yn==1){
        printf("covariance matrix not positive defined adding eps*cov[0][0]*I \n");
        for(i=0;i<en_tot+Nvar*en[0];i++)
             cov_yx[i][i]+=cov_yx[0][0]*1e-12;
        yn=is_it_positive( cov_yx,  en_tot+Nvar*en[0]);  
        printf("now the matrix is positive defined.  %d\n",yn);
   }  
   double **cov_yx1=symmetric_matrix_inverse(en_tot+Nvar*en[0], cov_yx  );
   
   yn=is_it_positive( cov,  en_tot);
   while(yn==1){
        printf("covariance matrix not positive defined adding eps*cov[0][0]*I \n");
        for(i=0;i<en_tot;i++)
             cov[i][i]+=cov[0][0]*1e-12;
        yn=is_it_positive( cov,  en_tot);  
        printf("now the matrix is positive defined.  %d\n",yn);
   }  
   double **cov1=symmetric_matrix_inverse(en_tot, cov  );
   
   guess=guess_for_non_linear_fit_Nf(N, en,x[0], y[0] , Nvar,  Npar, fit_info.function,guess );
  
   
 /*  double *guess1=(double*)  malloc(sizeof(double)*(Npar+Nvar*en[0]));
   for(i=0;i<Npar;i++)
       guess1[i]=guess[i];
   
   for(i=Npar;i<Npar+Nvar*en[0];i++)
       guess1[i]=x[0][ (i-Npar)/Nvar ][  (i-Npar)%Nvar  ];
   free(guess);
   guess=guess1;
   
   
   guess1=non_linear_fit_Nf_sigmax_covariance( N, en ,x[0], sigmax, y[0] , Nvar,  Npar,  fit_info.function , guess ,cov_yx1);
   free(guess);
   guess=guess1;
   */
   //#pragma omp parallel for  private(tmp,i,count,n,e,im)  shared(N, en, y , Nvar,  Npar,guess,Njack,r,chi2,C,x,cov,cov_yx1,cov1)
   for (j=0;j<Njack;j++){
        //if (j==0){     }
        tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess );
        //tmp=non_linear_fit_Nf_sigmax( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_iterative( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess );
        //tmp=non_linear_fit_Nf_sigmax_covariance( N, en ,x[j], sigmax, y[j] , Nvar,  Npar,  fit_info.function , guess ,cov_yx1);
        //tmp=non_linear_fit_Nf_covariance(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess ,cov1);
       
        
        chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        //for (i=Npar;i<Npar+en[0]*Nvar;i++)
         //   printf("%.10f   ",fabs(tmp[i]-x[j][(i-Npar)/Nvar ][(i-Npar)%Nvar ] ) );
                //printf("\n");
        //printf("jacknife=%d of %d   chi2=%g\n",j, Njack,chi2[j]);

        C[j]=covariance_non_linear_fit_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function );            
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }         
        if (j==Njack-1){
            printf("w0/a[fm]     mu*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err    KM2/Kf           (Mpi^2/fpi^2)* (Kf^2/KM^2)\n"); 
            fprintf(fdat,"w0/a[fm]     mu*w0/aZp[]      (M_Pi w0/KM)^2 or fw/Kf   err    KM2/Kf           (Mpi^2/fpi^2)* (Kf^2/KM^2)\n"); 
            
            count=0;
            double newline_if_w0=x[j][0][1];
            std::vector<int>  myen={0,1,2,3,   8,4,5,6,   7};
            
            for (n=0;n<N;n++){
                printf("#function %d\n",n);
                for (auto e :myen){
                    KM=fit_info.function(2,Nvar,x[j][e+count],Npar,tmp);
                    Kf=fit_info.function(3,Nvar,x[j][e+count],Npar,tmp);
                    double *tmp1=(double*) malloc(sizeof(double)*Njack);
                    for (int jj=0;jj<Njack;jj++){
                        tmp1[jj]=gJ[e].M_PS_jack[0][jj]/gJ[e].f_PS_jack[0][jj];
                        tmp1[jj]*=Kf/KM;
                        tmp1[jj]*=tmp1[jj];
                    }
                    double *tmp2=mean_and_error(jack_files[0].sampling,Njack,tmp1);
                    if (n==0)
                        K=KM*KM;
                    else if (n==1)
                        K=Kf;
                    
                    if (newline_if_w0!=x[j][e+count][1]){
                        fprintf(fdat,"\n\n");
                        newline_if_w0=x[j][e+count][1] ;
                    }
                    printf("%.5f     %.5f      %.5f   %.5f     %.5f \t\t %.5f  %.5f\n",x[j][e+count][1],x[j][e+count][0],fit[e+count][0]/K,fit[e+count][1]/K,K, tmp2[0] ,tmp2[1]);
                    fprintf(fdat,"%.5g     %.5g      %.5g   %.5g     %.5g \t\t %.5g  %.5g\n",x[j][e+count][1],x[j][e+count][0],fit[e+count][0]/K,fit[e+count][1]/K,K, tmp2[0] ,tmp2[1]);
                    
                    
                    
                    
                    free(tmp2);free(tmp1);
                }
                count+=en[n];
            }
            
        }
        free(tmp);
        
   }  
   free_2(en_tot,cov);
   free_2(en_tot,cov1);
   free_2(en_tot+Nvar*en[0],cov_yx);
   free_2(en_tot+Nvar*en[0],cov_yx1);
   
  
   chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
   printf("$\\chi^2/dof=%f+-%f$  \n",chi2m[0],chi2m[1]);
   
  
  
   struct fit_result fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x,&sigmax, &chi2m, &rm,&r, &fit, &y,&chi2,&C);
   fclose(fdat);
   return fit_out;
    
} 
