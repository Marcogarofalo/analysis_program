#define FAV_treshold_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "global.hpp"
#include "resampling.hpp"
#include "non_linear_fit.hpp"
#include "continuum_reph.hpp"
#include "indices.hpp"
#include "global_reph.hpp"
#include "mutils.hpp"
// #include <omp.h> 




static void init_fit( int N, struct header *head ,int Njack, struct reph_jack *gJ,int Npar, double xg_min,struct ik_for_n *ik_comb, int **en,int *en_tot, double ***x, double **chi2m, double **rm, double ***r, double ***fit, double ****y,double **chi2,double ****C , double Mpit, double r0_min, double r0_max)
{
    int imoms,imomt,imom0,iG,i,n,e,j;
    int count;
    int im;
   *en_tot=0;
   
   *en=(int*) calloc(N,sizeof(int));

   for (e=0;e<ensembles_reph;e++){
           for(imoms=0;imoms<head[e].nmoms;imoms++){
           for(imomt=0;imomt<head[e].nmoms;imomt++){
           for(imom0=0;imom0<head[e].nmoms;imom0++){   
                for (n=0;n<N;n++){
                    iG=index_ensemble_twopt_G_fit(head[e],ik_comb[n].ikt,ik_comb[n].iks,imom0,imomt,imoms); 
                    im=index_ensemble_twopt_fit(head[e],0,0,0,0);
                    if (gJ[e].xG[iG][Njack-1]>xg_min  &&   gJ[e].M_PS[im][Njack-1]*gJ[e].w0[Njack-1]  < Mpit    &&   gJ[e].w0[Njack-1] > r0_min && gJ[e].w0[Njack-1] < r0_max){
                        (*en)[n]+=1;
                        if(imoms==0 && imomt==0 && imom0==1)
                            printf("e=%d\n",e);
                    }
                }
           }}}
   }
      
 
   for (n=0;n<N;n++)
   {  *en_tot+=(*en)[n];   }
   
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
static struct fit_result close_fit( int N, struct header *head ,int Njack, struct reph_jack *gJ,int Npar,int **en,int *en_tot, double ***x, double **chi2m, double **rm,double ***r, double ***fit, double ****y,double **chi2, double ****C)
{
    int imoms,imomt,imom0,iG,i,n,e,j;
   int count;
   
   free(*chi2m);
   free(*rm);

  
   count=0;
   for (n=0;n<N;n++){
       for (i=0;i<(*en)[n];i++){
           free((*fit)[i+count]);
           free((*x)[i+count]);
       }
       count+=(*en)[n];
   }
   free(*fit);     
   free(*x);

   
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


struct fit_result fit_FAV_pion_treshold_Mpi(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info, double r0_min, double r0_max){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit,***C; 
   int i,j,e,im;  
   int Npar=fit_info.Npar;
   int Nvar=fit_info.Nvar;//m_l, w0,M_PS^2,f_PS
   int ikt=0,iks=0;
   int imoms,imomt,imom0,iG;
   int n,count,N=fit_info.N;
   int *en;
   int en_tot=0, i_m;
   double xg_min=0.001;
   double Mpit=0.77;
   //Mpiw0   A60=0.92      B55=0.89   A40=0.756 B35=0.72
   
   struct ik_for_n *ik_comb;
   ik_comb=(struct ik_for_n*) malloc(sizeof(struct ik_for_n)*N);
   ik_comb[0].ikt=0;
   ik_comb[0].iks=0;

   init_fit(N,  head , Njack, gJ,Npar,xg_min, ik_comb,&en,&en_tot, &x, &chi2m, &rm, &r, &fit, &y,&chi2,&C,Mpit,r0_min, r0_max  );

   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=0.1;
 
   for (n=0;n<N;n++){
        count=0;
        for (e=0;e<ensembles_reph;e++){
           for(imoms=0;imoms<head[e].nmoms;imoms++){
           for(imomt=0;imomt<head[e].nmoms;imomt++){
           for(imom0=0;imom0<head[e].nmoms;imom0++){       
                iG=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);  
                i_m=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);      
                if (gJ[e].xG[iG][Njack-1]>xg_min  &&   gJ[e].M_PS[i_m][Njack-1]*gJ[e].w0[Njack-1]  < Mpit    &&   gJ[e].w0[Njack-1] > r0_min && gJ[e].w0[Njack-1] < r0_max){
                if(n==0){
                    if(AV==0){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FA[iG][j]   *gJ[e].ZV[j];
                        }         
                    }
                    if(AV==1){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FV[iG][j]   *gJ[e].ZA[j];
                        }         
                    }
                    if(AV==2){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FAp[iG][j]   *gJ[e].ZV[j]*gJ[e].w0[j];
                        }         
                    }
                    if(AV==3){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FA_from_H0[iG][j]   ;
                        }         
                    }
                    if(AV==4){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FV_from_H0[iG][j] *gJ[e].ZA[j]  ;
                        }         
                    }
                    if(AV==5){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FV_from_H0_HA[iG][j] *gJ[e].ZAV[j]  ;
                        }         
                    }
                    fit[count]=mean_and_error_jack_biased(Njack, rm);
                }
                
                for (j=0;j<jack_tot;j++){
                    y[j][count][0]=rm[j];
                    y[j][count][1]=fit[count][1];
                }
                               
                 x[count]=(double*) malloc(sizeof(double)*Nvar);
                 count++;
                }
           }}}
       }
   }

  
  // #pragma omp parallel for  private(tmp,i)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)
   for (j=0;j<Njack;j++){
       for (n=0;n<N;n++){
       count=0;
        for (e=0;e<ensembles_reph;e++){
        for(imoms=0;imoms<head[e].nmoms;imoms++){
        for(imomt=0;imomt<head[e].nmoms;imomt++){
        for(imom0=0;imom0<head[e].nmoms;imom0++){       
                iG=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);  
                i_m=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);      
                if (gJ[e].xG[iG][Njack-1]>xg_min  &&   gJ[e].M_PS[i_m][Njack-1]*gJ[e].w0[Njack-1]  < Mpit    &&   gJ[e].w0[Njack-1] > r0_min && gJ[e].w0[Njack-1] < r0_max){
                x[count][0]=head[e].k[head[e].nk+ikt]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                x[count][1]=gJ[e].w0[j];//w0
                x[count][2]=gJ[e].xG[iG][j];//x_gamma
                x[count][3]=gJ[e].M_PS[i_m][j]*gJ[e].w0[j];//x_gamma
                
                x[count][9]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3]);
                x[count][9]=2*sin(x[count][9]/2.);
                
                x[count][10]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3]);
                x[count][10]=2*sin(x[count][10]/2.);
                x[count][11]=gJ[e].Zf_PS[i_m][j] *  gJ[e].w0[j];

                count++;
                }
        }}}
         i_m=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);      
        }
       }
        if (j==0) guess=guess_for_non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function ,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function )/(en_tot-Npar);
        C[j]=covariance_non_linear_fit_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function );      
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }   
        
        free(tmp);
   }  
   
   free(guess);
   free(ik_comb);
   struct fit_result fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &chi2m, &rm,&r, &fit, &y,&chi2,&C);gdb_hook();
   return fit_out;
    
} 



struct fit_result fit_FAV_Kphys_treshold(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info, double r0_min, double r0_max){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit,***C; 
   int i,j,e,im;  
   int Npar=fit_info.Npar;
   int Nvar=fit_info.Nvar;//m_l, w0,M_PS^2,f_PS
   int ikt=0,iks=1;
   int imoms,imomt,imom0,iG,iG1;
   int n,count,N=1;
   int *en;
   double Mpit=0.77;
   //Mpiw0   A60=0.92      B55=0.89   A40=0.756 B35=0.72
   
   int en_tot=0, i_m, i_mp,i_m1;
   double xg_min=0.001;
   struct ik_for_n *ik_comb;
   ik_comb=(struct ik_for_n *) malloc(sizeof(struct ik_for_n )*N);
   ik_comb[0].ikt=0;
   ik_comb[0].iks=1;
 
   
   
   init_fit(N,  head , Njack, gJ,Npar,xg_min,ik_comb,&en,&en_tot, &x, &chi2m, &rm, &r, &fit, &y,&chi2,&C,Mpit,r0_min, r0_max);

   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++){
       guess[i]=-0.1;
        /*if (AV==0) {
             if(Npar>7) {guess[0]=0.1;guess[1]=-0.1;guess[2]=0.1;guess[3]=-0.1;guess[4]=0.1;guess[5]=-0.1;}
             else guess[i]=0.1;
        }
        else if (AV==1)  guess[i]=-0.1;
        else if (AV==2)  guess[i]=0.1;
        else if (AV==3)  guess[i]=0.1;
        else error(0==0,1,"fit_FAV_D","specify in the fit AV value AV=0 means A AV=1 means V");
        */
   }
   double mc,mc1,F,F1,m,b,p,k,E_p,E_k,kdp,aMD;
   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<ensembles_reph;e++){
           for(imoms=0;imoms<head[e].nmoms;imoms++){
           for(imomt=0;imomt<head[e].nmoms;imomt++){
           for(imom0=0;imom0<head[e].nmoms;imom0++){       
                iG=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);
                iG1=index_ensemble_twopt_G_fit(head[e],ikt,iks+1,imom0,imomt,imoms);
                i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
                if (gJ[e].xG[iG][Njack-1]>xg_min &&   gJ[e].M_PS[i_mp][Njack-1]*gJ[e].w0[Njack-1]  < Mpit    &&   gJ[e].w0[Njack-1] > r0_min && gJ[e].w0[Njack-1] < r0_max){
                
                    if(AV==0){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FA[iG][j]   *gJ[e].ZV[j];
                        }         
                    }
                    if(AV==1){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FV[iG][j]   *gJ[e].ZA[j];
                        }         
                    }   
                    if(AV==2){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FAp[iG][j]   *gJ[e].ZV[j]*gJ[e].w0[j];
                        }         
                    }
                    if(AV==3){
                        for (j=0;j<Njack;j++){       
                            mc=head[e].k[head[e].nk+iks]*gJ[e].w0[j]/gJ[e].Zp[j];//mc*w0
                            mc1=head[e].k[head[e].nk+iks+1]*gJ[e].w0[j]/gJ[e].Zp[j];//mc*w0
                            F=gJ[e].FA_from_H0[iG][j];
                            F1=gJ[e].FA_from_H0[iG1][j];
                            rm[j]=gJ[e].FA_from_H0[iG][j]   ;
                            m=( F-F1 )/(   mc  -  mc1   );
                            b=F-m*mc;
                            rm[j]=m*(result.msw[j])+b;
           
                        }         
                    }
                    if(AV==4){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FV_from_H0[iG][j] *gJ[e].ZA[j]  ;
                        }         
                    }
                    if(AV==5){
                        for (j=0;j<Njack;j++){   
                            
                            rm[j]=inter_2(head[e].k[head[e].nk+iks]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+iks+1]*gJ[e].w0[j]/gJ[e].Zp[j],
                            gJ[e].FV_from_H0_HA[iG][j] *gJ[e].ZAV[j],               gJ[e].FV_from_H0_HA[iG1][j] *gJ[e].ZAV[j], 
                            result.msw[j]);
                        }         
                    }
                    fit[count]=mean_and_error_jack_biased(Njack, rm);
                
                    for (j=0;j<jack_tot;j++){
                        y[j][count][0]=rm[j];
                        y[j][count][1]=fit[count][1];
                    }
                // printf("HERE n=%d e=%d  en_tot=%d   e+count= %d   y=%f\n",n,e,en_tot,count, y[Njack-1][count ][0]);
                                    
                    x[count]=(double*) malloc(sizeof(double)*Nvar);
                    count++;
                }
                
           }}}
                 //  printf("HERE end\n");
          
       }
      // count+=en[n];
   }

  // #pragma omp parallel for  private(tmp,i)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)
   for (j=0;j<Njack;j++){
       count=0;
       for (n=0;n<N;n++){
        for (e=0;e<ensembles_reph;e++){
        for(imoms=0;imoms<head[e].nmoms;imoms++){
        for(imomt=0;imomt<head[e].nmoms;imomt++){
        for(imom0=0;imom0<head[e].nmoms;imom0++){       
                iG=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);  
                iG1=index_ensemble_twopt_G_fit(head[e],ikt,iks+1,imom0,imomt,imoms);  
                i_m=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);
                i_m1=index_ensemble_twopt_fit(head[e],ikt,iks+1,0,0);
                i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
                if (gJ[e].xG[iG][Njack-1]>xg_min &&   gJ[e].M_PS[i_mp][Njack-1]*gJ[e].w0[Njack-1]  < Mpit    &&   gJ[e].w0[Njack-1] > r0_min && gJ[e].w0[Njack-1] < r0_max){
                    x[count][0]=head[e].k[head[e].nk]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                    x[count][1]=gJ[e].w0[j];//w0
                    x[count][2]=gJ[e].xG[iG][j];//x_gamma
                    x[count][3]=gJ[e].M_PS[i_mp][j]*gJ[e].w0[j];//M_\pi r0
                    
                    x[count][4]=head[e].k[head[e].nk+1+n]*gJ[e].w0[j]/gJ[e].Zp[j];//ms*w0
                    x[count][5]=gJ[e].M_PS[i_m][j]*gJ[e].w0[j];//M_K r0
                    
                    x[count][6]=head[e].k[head[e].nk+ikt+n]*gJ[e].w0[j]/gJ[e].Zp[j];//mc*w0
                   
                    
                    x[count][9]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3]);
                    x[count][9]=2*sin(x[count][9]/2.);
                    
                    x[count][10]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3]);
                    x[count][10]=2*sin(x[count][10]/2.);
           
                    p=2.*sin(  (2.*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3])/2.);
                    k=2.*sin(  (2.*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3])/2.);
                    
                    
                    
                    aMD=inter_2(head[e].k[head[e].nk+iks]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+iks+1]*gJ[e].w0[j]/gJ[e].Zp[j],
                            gJ[e].M_PS[i_m][j],                 gJ[e].M_PS[i_m1][j], 
                            result.msw[j]);
                    x[count][5]=aMD*gJ[e].w0[j];//M_D r0
                    E_p=aMD*aMD;
                    E_p+=p*p;
                    E_p=sqrt(E_p);
                    E_k=2.*asinh(   sqrt(k*k)/2.   );
                    kdp=E_p*E_k- p*k; 
                    x[count][2]=2*kdp/(aMD*aMD);//  /M_D^2
 
                    x[count][11]=gJ[e].Zf_PS[i_mp][j] *  gJ[e].w0[j];
                    x[count][12]=inter_2(head[e].k[head[e].nk+iks]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+iks+1]*gJ[e].w0[j]/gJ[e].Zp[j],
                            gJ[e].Zf_PS[i_m][j] *  gJ[e].w0[j],                gJ[e].Zf_PS[i_m1][j] *  gJ[e].w0[j], 
                            result.msw[j]);
                    
                    count++;
                }
                
        }}}
         i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
         i_m=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);
         //printf("imp=%d im=%d  ikt=%d  iks=%d    n=%d\n",i_mp,i_m,ikt,iks,n);
//        if (j==Njack-1)
//            printf("mw=%f;   w0=%f;   ZV=%f; ZA=%f; mu=%f; L=%d; Mpi=%f;  a=%f; mcw=%f;  MD=%f; \n",x[count-1][0],gJ[e].w0[j],gJ[e].ZV[j],gJ[e].ZA[j],head[e].k[head[e].nk+0],head[e].l1,gJ[e].M_PS[i_mp][j],jack_files[e].a,x[count-1][6],gJ[e].M_PS[i_m][j]);
        }
       }
        //tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, FA_FV_chiral_continuum,guess );
        //chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, FA_FV_chiral_continuum  );
        if (j==0) guess=guess_for_non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );  

        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        C[j]=covariance_non_linear_fit_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function );      
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }         

        free(tmp);

   }  
   
   free(guess);
   free(ik_comb);
   struct fit_result fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &chi2m, &rm,&r, &fit, &y,&chi2,&C);
   return fit_out;
    
} 

struct fit_result fit_FAV_Dphys_treshold(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info, double r0_min, double r0_max){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit,***C; 
   int i,j,e,im;  
   int Npar=fit_info.Npar;
   int Nvar=fit_info.Nvar;//m_l, w0,M_PS^2,f_PS
   int ikt=3,iks=0;
   int imoms,imomt,imom0,iG,iG1;
   int n,count,N=1;
   int *en;
   double Mpit=1e+5;
   
   
   int en_tot=0, i_m, i_mp,i_m1;
   double xg_min=0.001;
   struct ik_for_n *ik_comb;
   ik_comb=(struct ik_for_n *) malloc(sizeof(struct ik_for_n )*N);
   ik_comb[0].ikt=3;
   ik_comb[0].iks=0;
   
   
   init_fit(N,  head , Njack, gJ,Npar,xg_min,ik_comb,&en,&en_tot, &x, &chi2m, &rm, &r, &fit, &y,&chi2,&C,Mpit,r0_min, r0_max);

   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++){
       guess[i]=0.1;
        /*if (AV==0) {
             if(Npar>7) {guess[0]=0.1;guess[1]=-0.1;guess[2]=0.1;guess[3]=-0.1;guess[4]=0.1;guess[5]=-0.1;}
             else guess[i]=0.1;
        }
        else if (AV==1)  guess[i]=-0.1;
        else if (AV==2)  guess[i]=0.1;
        else if (AV==3)  guess[i]=0.1;
        else error(0==0,1,"fit_FAV_D","specify in the fit AV value AV=0 means A AV=1 means V");
        */
   }
   double mc,mc1,F,F1,m,b,p,k,E_p,E_k,kdp,aMD;
   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<ensembles_reph;e++){
           for(imoms=0;imoms<head[e].nmoms;imoms++){
           for(imomt=0;imomt<head[e].nmoms;imomt++){
           for(imom0=0;imom0<head[e].nmoms;imom0++){       
                iG=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);
                iG1=index_ensemble_twopt_G_fit(head[e],ikt+1,iks,imom0,imomt,imoms);
                i_m=index_ensemble_twopt_fit(head[e],0,0,0,0);      
                if (gJ[e].xG[iG][Njack-1]>xg_min  &&   gJ[e].M_PS[i_m][Njack-1]*gJ[e].w0[Njack-1]  < Mpit    &&   gJ[e].w0[Njack-1] > r0_min && gJ[e].w0[Njack-1] < r0_max){
                
                    if(AV==0){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FA[iG][j]   *gJ[e].ZV[j];
                        }         
                    }
                    if(AV==1){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FV[iG][j]   *gJ[e].ZA[j];
                        }         
                    }   
                    if(AV==2){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FAp[iG][j]   *gJ[e].ZV[j]*gJ[e].w0[j];
                        }         
                    }
                    if(AV==3){
                        for (j=0;j<Njack;j++){       
                            mc=head[e].k[head[e].nk+ikt]*gJ[e].w0[j]/gJ[e].Zp[j];//mc*w0
                            mc1=head[e].k[head[e].nk+ikt+1]*gJ[e].w0[j]/gJ[e].Zp[j];//mc*w0
                            F=gJ[e].FA_from_H0[iG][j];
                            F1=gJ[e].FA_from_H0[iG1][j];
                            rm[j]=gJ[e].FA_from_H0[iG][j]   ;
                            m=( F-F1 )/(   mc  -  mc1   );
                            b=F-m*mc;
                            rm[j]=m*(result.mcw[j])+b;
           
                        }         
                    }
                    if(AV==4){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FV_from_H0[iG][j] *gJ[e].ZA[j]  ;
                        }         
                    }
                    if(AV==5){
                        for (j=0;j<Njack;j++){   
                            rm[j]=inter_2(head[e].k[head[e].nk+ikt]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+ikt+1]*gJ[e].w0[j]/gJ[e].Zp[j],
                            gJ[e].FV_from_H0_HA[iG][j] *gJ[e].ZAV[j],               gJ[e].FV_from_H0_HA[iG1][j] *gJ[e].ZAV[j], 
                            result.mcw[j]);
                        }         
                    }
                    fit[count]=mean_and_error_jack_biased(Njack, rm);
                
                    for (j=0;j<jack_tot;j++){
                        y[j][count][0]=rm[j];
                        y[j][count][1]=fit[count][1];
                    }
                // printf("HERE n=%d e=%d  en_tot=%d   e+count= %d   y=%f\n",n,e,en_tot,count, y[Njack-1][count ][0]);
                                    
                    x[count]=(double*) malloc(sizeof(double)*Nvar);
                    count++;
                }
                
           }}}
                 //  printf("HERE end\n");
          
       }
      // count+=en[n];
   }

  // #pragma omp parallel for  private(tmp,i)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)
   for (j=0;j<Njack;j++){
       count=0;
       for (n=0;n<N;n++){
        for (e=0;e<ensembles_reph;e++){
        for(imoms=0;imoms<head[e].nmoms;imoms++){
        for(imomt=0;imomt<head[e].nmoms;imomt++){
        for(imom0=0;imom0<head[e].nmoms;imom0++){       
                iG=index_ensemble_twopt_G_fit(head[e],ikt+n,iks,imom0,imomt,imoms);  
                iG1=index_ensemble_twopt_G_fit(head[e],ikt+1,iks,imom0,imomt,imoms);  
                i_m=index_ensemble_twopt_fit(head[e],ikt+n,iks,0,0);
                i_m1=index_ensemble_twopt_fit(head[e],ikt+1,iks,0,0);
                i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
                if (gJ[e].xG[iG][Njack-1]>xg_min  &&   gJ[e].M_PS[i_mp][Njack-1]*gJ[e].w0[Njack-1]  < Mpit    &&   gJ[e].w0[Njack-1] > r0_min && gJ[e].w0[Njack-1] < r0_max){
                    x[count][0]=head[e].k[head[e].nk+ikt]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                    x[count][1]=gJ[e].w0[j];//w0
                    x[count][2]=gJ[e].xG[iG][j];//x_gamma
                    x[count][3]=gJ[e].M_PS[i_mp][j]*gJ[e].w0[j];//M_\pi r0
                    
                    x[count][4]=head[e].k[head[e].nk+1+n]*gJ[e].w0[j]/gJ[e].Zp[j];//ms*w0
                    x[count][5]=gJ[e].M_PS[i_m][j]*gJ[e].w0[j];//M_K r0
                    
                    x[count][6]=head[e].k[head[e].nk+ikt+n]*gJ[e].w0[j]/gJ[e].Zp[j];//mc*w0
                   
                    
                    x[count][9]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3]);
                    x[count][9]=2*sin(x[count][9]/2.);
                    
                    x[count][10]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3]);
                    x[count][10]=2*sin(x[count][10]/2.);
           
                    p=2.*sin(  (2.*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3])/2.);
                    k=2.*sin(  (2.*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3])/2.);
                    
                    
                    
                    aMD=inter_2(head[e].k[head[e].nk+ikt]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+ikt+1]*gJ[e].w0[j]/gJ[e].Zp[j],
                            gJ[e].M_PS[i_m][j],                 gJ[e].M_PS[i_m1][j], 
                            result.mcw[j]);
                    x[count][7]=aMD*gJ[e].w0[j];//M_D r0
                    E_p=aMD*aMD;
                    E_p+=p*p;
                    E_p=sqrt(E_p);
                    E_k=2.*asinh(   sqrt(k*k)/2.   );
                    kdp=E_p*E_k- p*k; 
                    x[count][2]=2*kdp/(aMD*aMD);//  /M_D^2
 
                    x[count][11]=gJ[e].Zf_PS[i_mp][j] *  gJ[e].w0[j];
                    
                    
                    count++;
                }
                
        }}}
         i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
         i_m=index_ensemble_twopt_fit(head[e],ikt+n,iks,0,0);
         //printf("imp=%d im=%d  ikt=%d  iks=%d    n=%d\n",i_mp,i_m,ikt,iks,n);
//        if (j==Njack-1)
//            printf("mw=%f;   w0=%f;   ZV=%f; ZA=%f; mu=%f; L=%d; Mpi=%f;  a=%f; mcw=%f;  MD=%f; \n",x[count-1][0],gJ[e].w0[j],gJ[e].ZV[j],gJ[e].ZA[j],head[e].k[head[e].nk+0],head[e].l1,gJ[e].M_PS[i_mp][j],jack_files[e].a,x[count-1][6],gJ[e].M_PS[i_m][j]);
        }
       }
        //tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, FA_FV_chiral_continuum,guess );
        //chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, FA_FV_chiral_continuum  );
        if (j==0) guess=guess_for_non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );  

        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        C[j]=covariance_non_linear_fit_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function );      
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }         

        free(tmp);

   }  
   
   free(guess);
   free(ik_comb);
   struct fit_result fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &chi2m, &rm,&r, &fit, &y,&chi2,&C);
   return fit_out;
    
}





struct fit_result fit_FAV_Dsphys_treshold(struct database_file_reph_jack  *jack_files,  struct header *head ,int Njack, struct reph_jack *gJ ,int AV,struct fit_type fit_info,double r0_min, double r0_max){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit,***C; 
   int i,j,e,im;  
   int Npar=fit_info.Npar;
   int Nvar=fit_info.Nvar;//m_l, w0,M_PS^2,f_PS
   int ikt=3,iks=1;
   int imoms,imomt,imom0,iG,iG1,iG2,iG3;
   int n,count,N=1;
   int *en;
   double Mpit=1e+5;
   
   
   int en_tot=0, i_m, i_mp,i_m1;
   double xg_min=0.001;
   struct ik_for_n *ik_comb;
   ik_comb=(struct ik_for_n *) malloc(sizeof(struct ik_for_n )*N);
   ik_comb[0].ikt=3;
   ik_comb[0].iks=1;
   
  
   init_fit(N,  head , Njack, gJ,Npar,xg_min,ik_comb,&en,&en_tot, &x, &chi2m, &rm, &r, &fit, &y,&chi2,&C,Mpit, r0_min,r0_max);

   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++){
       guess[i]=-0.1;
        /*if (AV==0) {
             if(Npar>7) {guess[0]=0.1;guess[1]=-0.1;guess[2]=0.1;guess[3]=-0.1;guess[4]=0.1;guess[5]=-0.1;}
             else guess[i]=0.1;
        }
        else if (AV==1)  guess[i]=-0.1;
        else if (AV==2)  guess[i]=0.1;
        else if (AV==3)  guess[i]=0.1;
        else error(0==0,1,"fit_FAV_D","specify in the fit AV value AV=0 means A AV=1 means V");
        */
   }
   double mc,mc1,F,F1,m,b,p,k,E_p,E_k,kdp,aMD;
   count=0;
   for (n=0;n<N;n++){
        for (e=0;e<ensembles_reph;e++){
           for(imoms=0;imoms<head[e].nmoms;imoms++){
           for(imomt=0;imomt<head[e].nmoms;imomt++){
           for(imom0=0;imom0<head[e].nmoms;imom0++){       
                iG=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);
                iG1=index_ensemble_twopt_G_fit(head[e],ikt,iks+1,imom0,imomt,imoms);
                iG2=index_ensemble_twopt_G_fit(head[e],ikt+1,iks,imom0,imomt,imoms);
                iG3=index_ensemble_twopt_G_fit(head[e],ikt+1,iks+1,imom0,imomt,imoms);
                i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);

                if (gJ[e].xG[iG][Njack-1]>xg_min  &&   gJ[e].M_PS[i_mp][Njack-1]*gJ[e].w0[Njack-1]  < Mpit    &&   gJ[e].w0[Njack-1] > r0_min && gJ[e].w0[Njack-1] < r0_max){
                
                    if(AV==0){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FA[iG][j]   *gJ[e].ZV[j];
                        }         
                    }
                    if(AV==1){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FV[iG][j]   *gJ[e].ZA[j];
                        }         
                    }   
                    if(AV==2){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FAp[iG][j]   *gJ[e].ZV[j]*gJ[e].w0[j];
                        }         
                    }
                    if(AV==3){
                        for (j=0;j<Njack;j++){       
                            
                            rm[j]=inter_4( head[e].k[head[e].nk+iks]*gJ[e].w0[j]/gJ[e].Zp[j],  head[e].k[head[e].nk+iks+1]*gJ[e].w0[j]/gJ[e].Zp[j], 
                                           head[e].k[head[e].nk+ikt]*gJ[e].w0[j]/gJ[e].Zp[j], head[e].k[head[e].nk+ikt+1]*gJ[e].w0[j]/gJ[e].Zp[j],
                                           gJ[e].FA_from_H0[iG][j]      ,   gJ[e].FA_from_H0[iG1][j] ,
                                           gJ[e].FA_from_H0[iG2][j]      ,   gJ[e].FA_from_H0[iG3][j] ,
                                           result.msw[j] ,result.mcw[j]   );
           
                        }         
                    }
                    if(AV==4){
                        for (j=0;j<Njack;j++){                            
                            rm[j]=gJ[e].FV_from_H0[iG][j] *gJ[e].ZA[j]  ;
                        }         
                    }
                    if(AV==5){
                        for (j=0;j<Njack;j++){   
                            
                             
                            rm[j]=inter_4( head[e].k[head[e].nk+iks]*gJ[e].w0[j]/gJ[e].Zp[j],  head[e].k[head[e].nk+iks+1]*gJ[e].w0[j]/gJ[e].Zp[j], 
                                           head[e].k[head[e].nk+ikt]*gJ[e].w0[j]/gJ[e].Zp[j], head[e].k[head[e].nk+ikt+1]*gJ[e].w0[j]/gJ[e].Zp[j],
                                           gJ[e].FV_from_H0_HA[iG][j] *gJ[e].ZAV[j]     ,   gJ[e].FV_from_H0_HA[iG1][j] *gJ[e].ZAV[j],
                                           gJ[e].FV_from_H0_HA[iG2][j] *gJ[e].ZAV[j]     ,   gJ[e].FV_from_H0_HA[iG3][j] *gJ[e].ZAV[j],
                                           result.msw[j] ,result.mcw[j]   );
                            
                            
                        }         
                    }
                    fit[count]=mean_and_error_jack_biased(Njack, rm);
                
                    for (j=0;j<jack_tot;j++){
                        y[j][count][0]=rm[j];
                        y[j][count][1]=fit[count][1];
                    }
                // printf("HERE n=%d e=%d  en_tot=%d   e+count= %d   y=%f\n",n,e,en_tot,count, y[Njack-1][count ][0]);
                                    
                    x[count]=(double*) malloc(sizeof(double)*Nvar);
                    count++;
                }
                
           }}}
                 //  printf("HERE end\n");
          
       }
      // count+=en[n];
   }

   int i_ms1,i_ms2,i_mc1,i_mc2;
   int i_mDs1c1,i_mDs2c1,i_mDs1c2,i_mDs2c2;
  // #pragma omp parallel for  private(tmp,i)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)
   for (j=0;j<Njack;j++){
       count=0;
       for (n=0;n<N;n++){
        for (e=0;e<ensembles_reph;e++){
        for(imoms=0;imoms<head[e].nmoms;imoms++){
        for(imomt=0;imomt<head[e].nmoms;imomt++){
        for(imom0=0;imom0<head[e].nmoms;imom0++){       
                iG=index_ensemble_twopt_G_fit(head[e],ikt,iks,imom0,imomt,imoms);  
                iG1=index_ensemble_twopt_G_fit(head[e],ikt,iks+1,imom0,imomt,imoms);  
                iG2=index_ensemble_twopt_G_fit(head[e],ikt+1,iks,imom0,imomt,imoms);
                iG3=index_ensemble_twopt_G_fit(head[e],ikt+1,iks+1,imom0,imomt,imoms);
            
                i_ms1=index_ensemble_twopt_fit(head[e],0,iks,0,0);
                i_ms2=index_ensemble_twopt_fit(head[e],0,iks+1,0,0);
                
                i_mc1=index_ensemble_twopt_fit(head[e],ikt,0,0,0);
                i_mc2=index_ensemble_twopt_fit(head[e],ikt+1,0,0,0);
                
                i_mDs1c1=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);
                i_mDs2c1=index_ensemble_twopt_fit(head[e],ikt,iks+1,0,0);
                i_mDs1c2=index_ensemble_twopt_fit(head[e],ikt+1,iks,0,0);
                i_mDs2c2=index_ensemble_twopt_fit(head[e],ikt+1,iks+1,0,0);
                
                
                
                i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
                if (gJ[e].xG[iG][Njack-1]>xg_min  &&   gJ[e].M_PS[i_mp][Njack-1]*gJ[e].w0[Njack-1]  < Mpit    &&   gJ[e].w0[Njack-1] > r0_min && gJ[e].w0[Njack-1] < r0_max){
                    x[count][0]=head[e].k[head[e].nk]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                    x[count][1]=gJ[e].w0[j];//w0
                    x[count][2]=gJ[e].xG[iG][j];//x_gamma
                    x[count][3]=gJ[e].M_PS[i_mp][j]*gJ[e].w0[j];//M_\pi r0
                    
                    x[count][4]=head[e].k[head[e].nk+1+n]*gJ[e].w0[j]/gJ[e].Zp[j];//ms*w0
                    
                    x[count][6]=head[e].k[head[e].nk+ikt+n]*gJ[e].w0[j]/gJ[e].Zp[j];//mc*w0
                   
                    
                    x[count][9]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3]);
                    x[count][9]=2*sin(x[count][9]/2.);
                    
                    x[count][10]=(2*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3]);
                    x[count][10]=2*sin(x[count][10]/2.);
           
                    p=2.*sin(  (2.*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imoms][3])/2.);
                    k=2.*sin(  (2.*M_PI/head[e].l3)*(head[e].mom[imom0][3]-head[e].mom[imomt][3])/2.);
                    
                    
                    
                    aMD=inter_2(head[e].k[head[e].nk+iks]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+iks+1]*gJ[e].w0[j]/gJ[e].Zp[j],
                            gJ[e].M_PS[i_ms1][j],                 gJ[e].M_PS[i_ms2][j], 
                            result.msw[j]);
                    x[count][5]=aMD*gJ[e].w0[j];//M_K r0
                    aMD=inter_2(head[e].k[head[e].nk+ikt]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+ikt+1]*gJ[e].w0[j]/gJ[e].Zp[j],
                            gJ[e].M_PS[i_mc1][j],                 gJ[e].M_PS[i_mc2][j], 
                            result.mcw[j]);
                    x[count][7]=aMD*gJ[e].w0[j];//M_D r0
                    aMD=inter_4(head[e].k[head[e].nk+iks]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+iks+1]*gJ[e].w0[j]/gJ[e].Zp[j],
                                head[e].k[head[e].nk+ikt]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+ikt+1]*gJ[e].w0[j]/gJ[e].Zp[j],
                                gJ[e].M_PS[i_mDs1c1][j],                 gJ[e].M_PS[i_mDs2c1][j], 
                                gJ[e].M_PS[i_mDs1c2][j],                 gJ[e].M_PS[i_mDs2c2][j], 
                                result.msw[j], result.mcw[j]);
                    x[count][8]=aMD*gJ[e].w0[j];//M_D r0
                    E_p=aMD*aMD;
                    E_p+=p*p;
                    E_p=sqrt(E_p);
                    E_k=2.*asinh(   sqrt(k*k)/2.   );
                    kdp=E_p*E_k- p*k; 
                    x[count][2]=2*kdp/(aMD*aMD);//  /M_D^2
 
                    x[count][11]=gJ[e].Zf_PS[i_mp][j] *  gJ[e].w0[j];
                    x[count][12]=inter_2(head[e].k[head[e].nk+iks]*gJ[e].w0[j]/gJ[e].Zp[j],        head[e].k[head[e].nk+iks+1]*gJ[e].w0[j]/gJ[e].Zp[j],
                            gJ[e].Zf_PS[i_ms1][j] *  gJ[e].w0[j],                gJ[e].Zf_PS[i_ms2][j] *  gJ[e].w0[j], 
                            result.msw[j]);
                    
                    
                    count++;
                }
                
        }}}
         i_mp=index_ensemble_twopt_fit(head[e],0,0,0,0);
         i_m=index_ensemble_twopt_fit(head[e],ikt,iks,0,0);
         //printf("imp=%d im=%d  ikt=%d  iks=%d    n=%d\n",i_mp,i_m,ikt,iks,n);
//        if (j==Njack-1)
//            printf("mw=%f;   w0=%f;   ZV=%f; ZA=%f; mu=%f; L=%d; Mpi=%f;  a=%f; mcw=%f;  MD=%f; \n",x[count-1][0],gJ[e].w0[j],gJ[e].ZV[j],gJ[e].ZA[j],head[e].k[head[e].nk+0],head[e].l1,gJ[e].M_PS[i_mp][j],jack_files[e].a,x[count-1][6],gJ[e].M_PS[i_m][j]);
        }
       }
        //tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, FA_FV_chiral_continuum,guess );
        //chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, FA_FV_chiral_continuum  );
        if (j==0){ guess=guess_for_non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );  }

        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess ).P;
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        C[j]=covariance_non_linear_fit_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function );      
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }         

        free(tmp);

   }  
   
   free(guess);
   free(ik_comb);
   struct fit_result fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &chi2m, &rm,&r, &fit, &y,&chi2,&C);
   return fit_out;
    
} 
