#define global_fit_KandD_C


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
 #include "indices.hpp"

#include <omp.h> 


struct ik_for_n
{
  int ikt;
  int iks;  
};
void init_fit( int N, struct header *head ,int Njack, struct data_jack *gJ,int Npar, double xg_min,struct ik_for_n *ik_comb, int **en,int *en_tot, double ***x, double **chi2m, double **rm, double ***r, double ***fit, double ****y,double **chi2,double ****C)
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
       *en_tot+=(*en)[n];
   
   *x=(double**) malloc(sizeof(double*)*(*en_tot));

   *chi2m=(double*) malloc(sizeof(double)*(Npar));
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
struct fit_result close_fit( int N, struct header *head ,int Njack, struct data_jack *gJ,int Npar,int **en,int *en_tot, double ***x, double **chi2m, double **rm,double ***r, double ***fit, double ****y,double **chi2, double ****C)
{
    int imoms,imomt,imom0,iG,i,n,e,j;
   int count;
   free(*x);
   
   free(*chi2m);
   free(*rm);

   free(*fit);
  
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

struct fit_result fit_fK_of_M(struct database_file_jack  *jack_files,  struct header *head ,int Njack, struct data_jack *gJ,int ***mass_index, struct result_jack *r1 ,struct fit_type fit_info){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit,***C; 
   int i,j,e,im;  
   int Npar=fit_info.Npar;
   int Nvar=fit_info.Nvar;//m_l, w0,M_PS^2,f_PS
   int ikt=0,iks=1;
   int imoms,imomt,imom0,iG;
   int n,count,N=fit_info.N;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=0;
   en[1]=0;
   int en_tot=0, i_m, i_mp;
   double xg_min=0.001;
   struct ik_for_n *ik_comb;
   ik_comb=(struct ik_for_n *) malloc(sizeof(struct ik_for_n )*N);
   ik_comb[0].ikt=1;
   ik_comb[0].iks=0;
   
   ik_comb[1].ikt=2;
   ik_comb[1].iks=0;
   
   ik_comb[2].ikt=3;
   ik_comb[2].iks=0;
   
   init_fit(N,  head , Njack, gJ,Npar,xg_min,ik_comb,&en,&en_tot, &x, &chi2m, &rm, &r, &fit, &y,&chi2,&C);

   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=0.1;

   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<ensembles;e++){
                    im=mass_index[e][ik_comb[n].ikt][ik_comb[n].iks];//index_ensemble_twopt_G_fit(head[e],ikt,iks+n,imom0,imomt,imoms);
              
                
                    for (j=0;j<Njack;j++)                           
                            rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    
                    fit[count]=mean_and_error_jack_biased(Njack, rm);
                
                    for (j=0;j<jack_tot;j++){
                        y[j][count][0]=rm[j];
                        y[j][count][1]=fit[count][1];
                    }
                // printf("HERE n=%d e=%d  en_tot=%d   e+count= %d   y=%f\n",n,e,en_tot,count, y[Njack-1][count ][0]);
                                    
                    x[count]=(double*) malloc(sizeof(double)*Nvar);
                    count++;
            
                
          
                 //  printf("HERE end\n");
          
       }
      // count+=en[n];
   }
       

  // #pragma omp parallel for  private(tmp,i)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)
   for (j=0;j<Njack;j++){
       count=0;
       for (n=0;n<N;n++){
        for (e=0;e<ensembles;e++){
                im=mass_index[e][ik_comb[n].ikt][ik_comb[n].iks];
                
                
                  
                x[count][0]=gJ[e].M_PS_jack[0][j]*gJ[e].w0[j];//M_Pi*w0
                //x[e+count][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[count][1]=gJ[e].M_PS_jack[im][j]*gJ[e].w0[j];//MK*w0
                
                x[count][2]=gJ[e].w0[j];//w0
                
                x[count][3]=gJ[e].M_PS_jack[0][j]*gJ[e].M_PS_jack[0][j];//MPi^2
                x[count][4]=gJ[e].f_PS_jack[0][j];//f_Pi
                x[count][5]=gJ[e].M_PS_jack[im][j]*gJ[e].M_PS_jack[im][j];//MK2
                x[count][6]=gJ[e].f_PS_jack[im][j];//fkw
           
                x[count][7]=double(head[e].l1)/gJ[e].w0[j];//f
                x[count][8]=r1->Bw[j];
                x[count][9]=r1->fw_from_M[j];
                
                
                count++;
                
                
        
        
         //        if (j==Njack-1)
//            printf("mw=%f;   w0=%f;   ZV=%f; ZA=%f; mu=%f; L=%d; Mpi=%f;  a=%f; msw=%f;  MK=%f; \n",x[count-1][0],gJ[e].w0[j],gJ[e].ZV[j],gJ[e].ZA[j],head[e].k[head[e].nk+0],head[e].l1,gJ[e].M_PS[i_mp][j],jack_files[e].a,x[count-1][4],gJ[e].M_PS[i_m][j]);
 
        }
       }
        if (j==0) guess=guess_for_non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );

        non_linear_fit_result single_jack_fit=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );
        tmp= single_jack_fit.P;
        chi2[j]=single_jack_fit.chi2/(en_tot-Npar);
        C[j]=covariance_non_linear_fit_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function );      
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }         

        free(tmp);

   }  
   struct fit_result fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &chi2m, &rm,&r, &fit, &y,&chi2,&C);
   return fit_out;
    
} 



struct fit_result global_fit_Omega_from_Mpi_MK(struct database_file_jack  *jack_files,  struct header *head ,int Njack, struct data_jack *gJ,int ***mass_index, struct result_jack *r1 ,struct fit_type fit_info){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit,***C; 
   int i,j,e,im;  
   int Npar=fit_info.Npar;
   int Nvar=fit_info.Nvar;//m_l, w0,M_PS^2,f_PS
   int ikt=0,iks=1;
   int imoms,imomt,imom0,iG;
   int n,count,N=fit_info.N;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=0;
   en[1]=0;
   int en_tot=0, i_m, i_mp;
   double xg_min=0.001;
   struct ik_for_n *ik_comb;
   ik_comb=(struct ik_for_n *) malloc(sizeof(struct ik_for_n )*N);
   ik_comb[0].ikt=1;
   ik_comb[0].iks=0;
   
   ik_comb[1].ikt=2;
   ik_comb[1].iks=0;
   
   ik_comb[2].ikt=3;
   ik_comb[2].iks=0;
   
   init_fit(N,  head , Njack, gJ,Npar,xg_min,ik_comb,&en,&en_tot, &x, &chi2m, &rm, &r, &fit, &y,&chi2,&C);

   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=0.1;

   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<ensembles;e++){
                    im=mass_index[e][ik_comb[n].ikt][ik_comb[n].iks];//index_ensemble_twopt_G_fit(head[e],ikt,iks+n,imom0,imomt,imoms);
              
                
                    for (j=0;j<Njack;j++)                           
                            rm[j]=gJ[e].M_Omega_jack[n][j]   *  gJ[e].w0[j];
                    
                    fit[count]=mean_and_error_jack_biased(Njack, rm);
                
                    for (j=0;j<jack_tot;j++){
                        y[j][count][0]=rm[j];
                        y[j][count][1]=fit[count][1];
                    }
                // printf("HERE n=%d e=%d  en_tot=%d   e+count= %d   y=%f\n",n,e,en_tot,count, y[Njack-1][count ][0]);
                                    
                    x[count]=(double*) malloc(sizeof(double)*Nvar);
                    count++;
            
                
          
                 //  printf("HERE end\n");
          
       }
      // count+=en[n];
   }
       

  // #pragma omp parallel for  private(tmp,i)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)
   for (j=0;j<Njack;j++){
       count=0;
       for (n=0;n<N;n++){
        for (e=0;e<ensembles;e++){
                im=mass_index[e][ik_comb[n].ikt][ik_comb[n].iks];
                
                
                  
                x[count][0]=gJ[e].M_PS_jack[0][j]*gJ[e].w0[j];//M_Pi*w0
                //x[e+count][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[count][1]=gJ[e].M_PS_jack[im][j]*gJ[e].w0[j];//MK*w0
                x[count][2]=gJ[e].w0[j];//w0
                /*
                x[count][3]=gJ[e].M_PS_jack[0][j]*gJ[e].M_PS_jack[0][j];//MPi^2
                x[count][4]=gJ[e].f_PS_jack[0][j];//f_Pi
                x[count][5]=gJ[e].M_PS_jack[im][j]*gJ[e].M_PS_jack[im][j];//MK2
                x[count][6]=gJ[e].f_PS_jack[im][j];//fkw
           
                x[count][7]=double(head[e].l1)/gJ[e].w0[j];//f
                x[count][8]=r1->Bw[j];
                x[count][9]=r1->fw_from_M[j];
              */  
                
                count++;
                
                
        
        
         //        if (j==Njack-1)
//            printf("mw=%f;   w0=%f;   ZV=%f; ZA=%f; mu=%f; L=%d; Mpi=%f;  a=%f; msw=%f;  MK=%f; \n",x[count-1][0],gJ[e].w0[j],gJ[e].ZV[j],gJ[e].ZA[j],head[e].k[head[e].nk+0],head[e].l1,gJ[e].M_PS[i_mp][j],jack_files[e].a,x[count-1][4],gJ[e].M_PS[i_m][j]);
 
        }
    }
        if (j==0) guess=guess_for_non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );

        non_linear_fit_result single_jack_fit=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );
        tmp=single_jack_fit.P;
        chi2[j]=single_jack_fit.chi2/(en_tot-Npar);
        C[j]=covariance_non_linear_fit_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function );      
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }         

        free(tmp);

   }  
   struct fit_result fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &chi2m, &rm,&r, &fit, &y,&chi2,&C);
   return fit_out;
    
} 
