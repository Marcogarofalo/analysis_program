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

int ensemble_omega_jacob=8;


double ***data_omega(char *option, int Njack  , struct header *head){
    
    double ***r=double_malloc_3(ensemble_omega_jacob, 3 ,Njack );
    //double **strange_jacob=double_malloc_2(ensemble_omega_jacob,3);
    double **mean=double_malloc_2(ensemble_omega_jacob,2);
    double ***cov=double_malloc_3(ensemble_omega_jacob,2,2);
    double **w0=double_malloc_2(ensemble_omega_jacob,2);
    int ik2_min=1;
    if (ensemble_omega_jacob>8){
        //B25.32
        int e=8;
        mean[e][0]=0.555225978284128; mean[e][1]= 7.76055858378567;
        cov[e][0][0]=2.438932e-05;     cov[e][0][1]=-0.0007886619;
        cov[e][1][0]=-0.0007886619;     cov[e][1][1]=0.0281728981;
        w0[e][0]=2.09818154701605  ; w0[e][1]= 2.09818154701605  ;
   }
   if (ensemble_omega_jacob>7){
            //C06
       int e=7;
       mean[e][0]=0.458088472927548;    mean[e][1]= 7.89927254874169;
       cov[e][0][0]=2.383969e-05;       cov[e][0][1]=-0.0008416204;
       cov[e][1][0]= -8.416204e-04;     cov[e][1][1]= 0.0330829020 ;
       w0[e][0]=2.50451148219  ;        w0[e][1]= 0.0017199806536  ;
    }
   if (ensemble_omega_jacob>6){
            //B072
       int e=6;
       mean[e][0]= 0.5469;    mean[e][1]=  7.922;
       cov[e][0][0]=1.084535e-05;       cov[e][0][1]= -0.0003198579;
       cov[e][1][0]= -3.198579e-04;     cov[e][1][1]= 0.0108084694 ;
       //cov[e][0][0]=0.0025*0.0025;       cov[e][0][1]= -0.0003198579;
       //cov[e][1][0]= -3.198579e-04;     cov[e][1][1]= 0.081* 0.081 ;
       
       w0[e][0]=2.12715290702  ;        w0[e][1]= 0.00194383531  ;
    }
   if (ensemble_omega_jacob>5){
            //B14
       int e=5;
       mean[e][0]=0.543264208692097;  mean[e][1]= 7.93881124536661;
       cov[e][0][0]=3.780476e-05;     cov[e][0][1]= -0.001208199;
       cov[e][1][0]= -1.208199e-03;   cov[e][1][1]= 0.041925784 ;
       w0[e][0]=2.11760645  ;         w0[e][1]= 0.001313096946  ;
   }
   if (ensemble_omega_jacob>4){
            //B25
       int e=4;
       mean[e][0]=0.555225978284128;  mean[e][1]= 7.76055858378567;
       cov[e][0][0]=2.438932e-05;     cov[e][0][1]= -0.0007886619;
       cov[e][1][0]= -7.886619e-04;   cov[e][1][1]=  0.0281728981;
       w0[e][0]=2.09818154701605  ;   w0[e][1]= 0.0019  ;
    }
   if (ensemble_omega_jacob>3){
       //A12
       int e=3;
       mean[e][0]=0.638555160256861; mean[e][1]= 7.94066602977537;
       cov[e][0][0]=0.0000492634;    cov[e][0][1]= -0.001270452;
       cov[e][1][0]=-0.0012704516;   cov[e][1][1]=  0.036395029;
       w0[e][0]=1.8249134215  ;      w0[e][1]= 0.0032611047  ;
   }
   if (ensemble_omega_jacob>2){
       //A30
       int e=2;
       mean[e][0]=0.645171588364307; mean[e][1]= 7.92242532027645;
       cov[e][0][0]=1.561641e-05;    cov[e][0][1]= -0.0004021516;
       cov[e][1][0]=-4.021516e-04;   cov[e][1][1]= 0.0114912103;       
       w0[e][0]=1.79280474 ;         w0[e][1]= 0.00166108579  ;
   } 
   if (ensemble_omega_jacob>1){
       //A40
       int e=1;
       mean[e][0]=0.675378596527479; mean[e][1]= 7.28044649567117;
       cov[e][0][0]=6.575932e-05;    cov[e][0][1]= -0.001829619;
       cov[e][1][0]=-1.829619e-03;   cov[e][1][1]=  0.054854808;
       w0[e][0]=1.77660378;          w0[e][1]= 0.0032616817  ;
   } 
   if (ensemble_omega_jacob>0){
       //A53
       int e=0;
       mean[e][0]=0.672938611351423; mean[e][1]= 7.50457062778335;
       cov[e][0][0]=3.094976e-05;    cov[e][0][1]= -0.0008455663;
       cov[e][1][0]=-8.455663e-04;   cov[e][1][1]=  0.0248091489;       
       w0[e][0]=1.759718279588;      w0[e][1]= 0.004311795563  ;
   } 
   
   for (int e=0;e<ensemble_omega_jacob;e++){
       //power up the errors
       cov[e][0][0]+=1.2e-3*1.2e-3; // routhly d(M_K^2) *slope
       //cov[e][1][1]*=1.1; 
       if (e==0 || e==1) {cov[e][0][0]+=1  ; cov[e][1][1]+=1  ;}
       //make the matrix symmetric
       cov[e][1][0]=cov[e][0][1];
       double  **P=fake_sampling_covariance(option,  mean[e],  Njack, 2 ,cov[e], e);
       double *Pw=fake_sampling(option,w0[e][0],w0[e][1],Njack,e+ensemble_omega_jacob);
       printf("e=%d  f1= %f  %f  f2=%f  %f\n",e,P[0][Njack-1], error_jackboot(option, Njack,P[0]) ,P[1][Njack-1], error_jackboot(option, Njack,P[1]) );
       for (int j =0; j< Njack; j++)
            for (int i=0 ; i< 3; i++){
              r[e][i][j]=P[0][j]+ P[1][j]*   head[e].k[head[e].nk+i +ik2_min]  ; 
              //r[e][i][j]*=Pw[j];
            }
       //for (int i=0 ; i< 3; i++)
        //    printf("mu=%f   omega=%f\n",head[e].k[head[e].nk+i +ik2_min] , r[e][i][Njack-1] );
       free_2(2, P);
       free(Pw);
   }
   free_2(ensemble_omega_jacob,w0);
   free_2(ensemble_omega_jacob,mean);
   free_3(ensemble_omega_jacob,2,cov);
   return r;
    
}


double ***data_omega_2state(char *option, int Njack  , struct header *head){
    
    double ***r=double_malloc_3(ensemble_omega_jacob, 3 ,Njack );
    //double **strange_jacob=double_malloc_2(ensemble_omega_jacob,3);
    double **mean=double_malloc_2(ensemble_omega_jacob,2);
    double ***cov=double_malloc_3(ensemble_omega_jacob,2,2);
    double **w0=double_malloc_2(ensemble_omega_jacob,2);
    int ik2_min=1;
    if (ensemble_omega_jacob>8){
        //B25.32
        int e=8;
        mean[e][0]=0.560717194370018; mean[e][1]= 7.47964876600531;
        cov[e][0][0]=2.378689e-05;     cov[e][0][1]=-0.0007529395;
        cov[e][1][0]=-0.0007529395;     cov[e][1][1]=0.0284132779;
        w0[e][0]=2.09818154701605  ; w0[e][1]= 2.09818154701605  ;
   }
   if (ensemble_omega_jacob>7){
            //C06
       int e=7;
       mean[e][0]=0.458413516863134; mean[e][1]= 7.93864137084707;
       cov[e][0][0]=5.437673e-06;     cov[e][0][1]=-0.0001625593;
       cov[e][1][0]= -1.625593e-04;     cov[e][1][1]= 0.0065089772 ;
        w0[e][0]=2.50451148219  ; w0[e][1]= 0.0017199806536  ;
    }
   if (ensemble_omega_jacob>6){
            //B072
       int e=6;
       mean[e][0]=0.5520; mean[e][1]= 7.721;
       cov[e][0][0]=5.938507e-06;     cov[e][0][1]= -0.0001559411;
       cov[e][1][0]= -1.559411e-04;     cov[e][1][1]= 0.0052884152 ;
       w0[e][0]=2.12715290702  ; w0[e][1]= 0.00194383531  ;
    }
   if (ensemble_omega_jacob>5){
            //B14
       int e=5;
       mean[e][0]=0.542550880775462; mean[e][1]= 7.88707387097471;
       cov[e][0][0]=0.0000232027;     cov[e][0][1]= -0.0006704852;
       cov[e][1][0]= -0.0006704852;     cov[e][1][1]= 0.0245880891 ;
       w0[e][0]=2.11760645  ; w0[e][1]= 0.001313096946  ;
   }
   if (ensemble_omega_jacob>4){
            //B25
       int e=4;
       mean[e][0]=0.554656884525161; mean[e][1]= 7.67302632107142;
       cov[e][0][0]=4.665973e-06;     cov[e][0][1]= -0.000133651;
       cov[e][1][0]= -1.336510e-04;     cov[e][1][1]= 0.004773261;
        w0[e][0]=2.09818154701605  ; w0[e][1]= 0.0019   ;
    }
   if (ensemble_omega_jacob>3){
       //A12
       int e=3;
       mean[e][0]=0.634678036526639; mean[e][1]= 7.96288367008374;
       cov[e][0][0]=2.505978e-05;     cov[e][0][1]= -0.000574197;
       cov[e][1][0]=-5.741970e-04;     cov[e][1][1]=  0.016208291;
        w0[e][0]=1.8249134215  ; w0[e][1]= 0.0032611047  ;
   }
   if (ensemble_omega_jacob>2){
       //A30
       int e=2;
       mean[e][0]=0.643284221508349; mean[e][1]= 7.94051970360903;
       cov[e][0][0]=1.030913e-05;     cov[e][0][1]= -0.0002451341;
       cov[e][1][0]=-2.451341e-04;     cov[e][1][1]= 0.0068209526;       
        w0[e][0]=1.79280474 ; w0[e][1]= 0.00166108579  ;
   } 
   if (ensemble_omega_jacob>1){
       //A40
       int e=1;
       mean[e][0]=0.670074239998651; mean[e][1]= 7.39070068608978;
       cov[e][0][0]=4.167383e-05;     cov[e][0][1]= -0.001069454;
       cov[e][1][0]=-1.069454e-03;     cov[e][1][1]=  0.031649366;
        w0[e][0]=1.77660378; w0[e][1]= 0.0032616817  ;
   } 
   if (ensemble_omega_jacob>0){
       //A53
       int e=0;
       mean[e][0]=0.675272570053787; mean[e][1]= 7.35567699766334;
       cov[e][0][0]=2.839506e-05;     cov[e][0][1]= -0.0007315468;
       cov[e][1][0]=-7.315468e-04;     cov[e][1][1]= 0.0215777123;
        w0[e][0]=1.759718279588; w0[e][1]= 0.004311795563  ;
   } 
   
   for (int e=0;e<ensemble_omega_jacob;e++){
       //power up the errors
       //if (e==6 ) {cov[e][0][0]+=100; cov[e][1][1]+=100; cov[e][0][1]=0; }
       //cov[e][0][0]*=2; cov[e][1][1]*=2;
       //make the matrix symmetric
       cov[e][1][0]=cov[e][0][1];
       double  **P=fake_sampling_covariance(option,  mean[e],  Njack, 2 ,cov[e], e);
       double *Pw=fake_sampling(option,w0[e][0],w0[e][1],Njack,e+ensemble_omega_jacob);
       printf("e=%d  f1= %f  %f  f2=%f  %f\n",e,P[0][Njack-1], error_jackboot(option, Njack,P[0]) ,P[1][Njack-1], error_jackboot(option, Njack,P[1]) );
       for (int j =0; j< Njack; j++)
            for (int i=0 ; i< 3; i++){
              r[e][i][j]=P[0][j]+ P[1][j]*   head[e].k[head[e].nk+i +ik2_min]  ; 
              //r[e][i][j]*=Pw[j];
            }
      for (int i=0 ; i< 3; i++)
            printf("mu=%f   omega=%f  +- %f\n",head[e].k[head[e].nk+i +ik2_min] , r[e][i][Njack-1], error_jackboot(option, Njack,r[e][i])  );
       free_2(2, P);
       free(Pw);
   }
   free_2(ensemble_omega_jacob,w0);
   free_2(ensemble_omega_jacob,mean);
   free_3(ensemble_omega_jacob,2,cov);
   return r;
    
}

double **w0_ensemble(char *option, int Njack  , struct header *head){
    
    double **r=double_malloc_2(ensemble_omega_jacob,Njack );
    //double **strange_jacob=double_malloc_2(ensemble_omega_jacob,3);
    double **w0=double_malloc_2(ensemble_omega_jacob,2);
    int ik2_min=1;
    if (ensemble_omega_jacob>8){
        //B25.32
        int e=8;
        
        w0[e][0]=2.09818154701605  ; w0[e][1]= 2.09818154701605  ;
   }
   if (ensemble_omega_jacob>7){
            //C06
       int e=7;
       
       w0[e][0]=2.50451148219  ;        w0[e][1]= 0.0017199806536  ;
    }
   if (ensemble_omega_jacob>6){
            //B072
       int e=6;
       
       w0[e][0]=2.12715290702  ;        w0[e][1]= 0.00194383531  ;
    }
   if (ensemble_omega_jacob>5){
            //B14
       int e=5;
       
       w0[e][0]=2.11760645  ;         w0[e][1]= 0.001313096946  ;
   }
   if (ensemble_omega_jacob>4){
            //B25
       int e=4;
       
       w0[e][0]=2.09818154701605  ;   w0[e][1]= 0.0019  ;
    }
   if (ensemble_omega_jacob>3){
       //A12
       int e=3;
       
       w0[e][0]=1.8249134215  ;      w0[e][1]= 0.0032611047  ;
   }
   if (ensemble_omega_jacob>2){
       //A30
       int e=2;
              
       w0[e][0]=1.79280474 ;         w0[e][1]= 0.00166108579  ;
   } 
   if (ensemble_omega_jacob>1){
       //A40
       int e=1;
       
       w0[e][0]=1.77660378;          w0[e][1]= 0.0032616817  ;
   } 
   if (ensemble_omega_jacob>0){
       //A53
       int e=0;
        w0[e][0]=1.759718279588;      w0[e][1]= 0.004311795563  ;
   } 
   
   for (int e=0;e<ensemble_omega_jacob;e++){
       
       double *Pw=fake_sampling(option,w0[e][0],w0[e][1],Njack,e+ensemble_omega_jacob);
       //printf("e=%d  f1= %f  %f  f2=%f  %f\n",e,P[0][Njack-1], error_jackboot(option, Njack,P[0]) ,P[1][Njack-1], error_jackboot(option, Njack,P[1]) );
       for (int j =0; j< Njack; j++)
              r[e][j]=Pw[j];
           
       
       free(Pw);
   }
   free_2(ensemble_omega_jacob,w0);
   
   return r;
    
}




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


double chiral_continuum_M_Omega_no_a(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi,r;
    double pi=3.141592653589793;
    
    double Mpiw=x[0],MKw=x[1], w0=x[2];
    
    r=P[0]+P[1]*Mpiw*Mpiw;
    
    return r;
    
}


double chiral_continuum_M_Omega_const(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi,r;
    double pi=3.141592653589793;
    
    double Mpiw=x[0],MKw=x[1], w0=x[2];
    
    r=P[0];
    
    return r;
    
}

  


double **fit_Omegaw0_from_M(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1){
   double ***y,**x,***r,**fK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=2, Nparline=2;
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
   
   double ***omega=data_omega(jack_files[0].sampling,  Njack  , head);
   double **w0e=w0_ensemble(jack_files[0].sampling,  Njack  , head);
   
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
            if (ensemble_omega_jacob>8){
             //B25.32
            strange_jacob[8][0]=0.0148;   
            strange_jacob[8][1]=0.0185;   
            strange_jacob[8][2]=0.0222; 
            }
            if (ensemble_omega_jacob>7){
            //C06
            strange_jacob[7][0]=0.0128;
            strange_jacob[7][1]=0.0161;   
            strange_jacob[7][2]=0.0193;
            }
            //B072
            strange_jacob[6][0]=0.017;   
            strange_jacob[6][1]=0.0195;   
            strange_jacob[6][2]=0.0220;   
            //B14
            strange_jacob[5][0]=0.0148;   
            strange_jacob[5][1]=0.0185;   
            strange_jacob[5][2]=0.0222;   
            //B25
            strange_jacob[4][0]=0.0148;   
            strange_jacob[4][1]=0.0185;   
            strange_jacob[4][2]=0.0222;   
            //A12
            strange_jacob[3][0]=0.0182;   
            strange_jacob[3][1]=0.0227;   
            strange_jacob[3][2]=0.0273;   
            //A30
            strange_jacob[2][0]=0.0182;   
            strange_jacob[2][1]=0.0227;   
            strange_jacob[2][2]=0.0273;   
            //A40
            strange_jacob[1][0]=0.0182;   
            strange_jacob[1][1]=0.0227;   
            strange_jacob[1][2]=0.0273;   
            //A53
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
                            //rm[j]=gJ[e].M_Omega_jack[ms][j]   *  gJ[e].w0[j]; // data old
                            //rm[j]=omega[e][ms][j]*  gJ[e].w0[j];
                            rm[j]=omega[e][ms][j]*w0e[e][j] ;
                    }
                    //printf("e=%d   ms=%d  n=%d   M=%f  %f  \n",e,ms,n,omega[e][ms][Njack-1],error_jackboot("jack",Njack,omega[e][ms] ) );
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
                //x[ms+n*nk][0]=r[0][e][j]+r[1][e][j]*strange_jacob[e][ms];
                x[ms+n*nk][0]=gJ[e].M_PS_jack[im][j];
                x[ms+n*nk][0]*= gJ[e].w0[j];// M_K*w0
                double xi=gJ[e].M_PS_jack[0][j]/(4*pi_greco*gJ[e].f_PS_jack[0][j]);
                xi=xi*xi;
                double delta=FVE_GL_Mpi(head[e].l1 , xi,   gjack[e].f_PS_jack[0][j]     );
                x[ms+n*nk][0]=x[ms+n*nk][0]*x[ms+n*nk][0]-  
                (gJ[e].M_PS_jack[0][j]   *  gJ[e].w0[j]*gJ[e].M_PS_jack[0][j]   *  gJ[e].w0[j])/(2.*((1-0.25*delta)* (1-0.25*delta) )); //M_K^2 w0^2 - M_pi^2w^2/2
            }
        }
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess );
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, one_line  );

       
   //     tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_parabola,guess );
   //     chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, one_parabola  );

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
                        rm[j]=r[0][e][j]+mref[ms]*r[1][e][j];
                        //rm[j]=r[0][e][j]+mref[ms]*r[1][e][j]+mref[ms]*mref[ms]*r[2][e][j];
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
            double Kf;
            for(i=0;i<Npar;i++)   
                printf("P%d=%f;  ",i,tmp[i]);
            printf("  MKw2-Mpi22=%f;\n",mref[ms]  );
            chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
            printf("$\\chi^2/dof=%f+-%f$\n",chi2[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
            free(chi2m);
              printf("#(M_Pi*w0)^2        M_Omega w0    err  Kf  w0\n");
              for (e=0;e<ensemble_omega_jacob;e++){
                //  FVE_K( r1->Bw[Njack-1], r1->fw[Njack-1], double(head[e].l1)/gJ[e].w0[Njack-1],  x[e][0]*x[e][0]/(2.*r1->Bw[Njack-1])/*mlw*/, (mref[ms])/r1->Bw[Njack-1]-x[e][0]*x[e][0]/(2.*r1->Bw[Njack-1]) /*msw*/ ,x[e][3],  x[e][4],x[e][5], x[e][6],&KM, &Kf);
                  Kf=1;
                  printf("%f        %f       %f    %f     %f\n",x[e][0]*x[e][0],y[Njack-1][e][0] ,y[Njack-1][e][1],Kf,gJ[e].w0[Njack-1]);
              }
         }
         for(i=0;i<Npar;i++)
             par[i][j]=tmp[i];   

       
         free(tmp);
    }
    for(i=0;i<Npar;i++){
        double *ave=mean_and_error(jack_files[0].sampling,Njack, par[i]);
        printf("P%d= %g   %g\t",i,ave[0],ave[1]);
        free(ave);
    }
    printf("\n");
    //extract w0 at given m_ref  after the chiral and continuum fit
    
    chi2m=mean_and_error(jack_files[0].sampling,Njack, par[1]);
    if(chi2m[0]<=chi2m[1]*0.5){
        printf("slope is zero\n");
        for (j=0;j<Njack;j++)
            fK[ms][j]=par[0][j]/r1->MOmegaMeV[j] ;
    }
    //div_jackboot(Njack,fK[ms], par[0],r1->MOmegaMeV );
    
    
    else {
       for (j=0;j<Njack;j++){
            fK[ms][j]=r1->MOmegaMeV[j]  - sqrt( r1->MOmegaMeV[j]*r1->MOmegaMeV[j]-4.* par[0][j]*par[1][j]*r1->MpiMeV[j]*r1->MpiMeV[j]  );
            fK[ms][j]=fK[ms][j]/(2.*par[1][j]*r1->MpiMeV[j]*r1->MpiMeV[j] );//w0 at that reference ms mass
            //compute Omega from w0
            //fK[ms][j]=par[0][j]+par[1][j]*r1->MpiMeV[j]*r1->MpiMeV[j] *(0.1701*0.1701/(197.3*197.3));
            //fK[ms][j]/=(0.1701/(197.3));
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
            // compute Omega from w0
            //out[0][j]=par[0][j]+par[1][j]*(r1->MKMeV[j]*r1->MKMeV[j]-r1->MpiMeV[j]*r1->MpiMeV[j]/2.)*(0.1701*0.1701/(197.3*197.3));
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
   
   free_3(ensemble_omega_jacob, 3  ,omega);
   return out;
    
} 

 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double **fit_Omegaw0_from_M_new(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ , struct result_jack *r1){
   double ***y,**x,***r,**fK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=2, Nparline=2;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik2_min=1, ik2_max=3;
   int nk=(ik2_max-ik2_min+1);
   int ms;

   double *mref;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nk);
   mref[0]=0.38*0.38;
   
   
   double ***omega=data_omega(jack_files[0].sampling,  Njack  , head);
   double **w0e=w0_ensemble(jack_files[0].sampling,  Njack  , head);

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
   
   
for (e=0;e<ensemble_omega_jacob;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                if (n==0){
                    for (j=0;j<Njack;j++){
                            //rm[j]=gJ[e].M_Omega_jack[ms][j]   *  gJ[e].w0[j]; // data old
                            //rm[j]=omega[e][ms][j]*  gJ[e].w0[j];
                            rm[j]=omega[e][ms][j]*w0e[e][j] ;
                            //rm[j]=gJ[e].M_PS_jack[im][j]*  gJ[e].w0[j];
                    }
                    //printf("e=%d   ms=%d  n=%d   M=%f  %f  \n",e,ms,n,omega[e][ms][Njack-1],error_jackboot("jack",Njack,omega[e][ms] ) );
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
                //x[ms+n*nk][0]=r[0][e][j]+r[1][e][j]*strange_jacob[e][ms];
                double xi=gJ[e].M_PS_jack[0][j]/(4*pi_greco*gJ[e].f_PS_jack[0][j]);
                xi=xi*xi;
                double delta=FVE_GL_Mpi(head[e].l1 , xi,   gjack[e].f_PS_jack[0][j]     );
                x[ms+n*nk][0]=gJ[e].M_PS_jack[im][j];
                x[ms+n*nk][0]=x[ms+n*nk][0]*x[ms+n*nk][0]-   (gJ[e].M_PS_jack[0][j]   *gJ[e].M_PS_jack[0][j]/((1-0.25*delta)* (1-0.25*delta))    )/2.; //M_K^2 w0^2 - M_pi^2w^2/2
                x[ms+n*nk][0]/=omega[e][ms][j]*omega[e][ms][j];
            }
        }
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_line,guess );
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, one_line  );

       
   //     tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, one_parabola,guess );
   //     chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, one_parabola  );

        for(i=0;i<Npar;i++){
            r[i][e][j]=tmp[i];
        }                
        
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        //printf("%g    %g    %g\n",x[ms][0],y[Njack-1][ms][0],y[Njack-1][ms][1]);
        free(x[ms]);  free(fit[ms]);
   }
   /*for(i=0;i<Npar;i++)
       printf("P%d=%g;  \t ",i,r[i][e][Njack-1]);
   printf("\n");
   */
   
   
    chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
    //printf("$\\chi^2/dof=%f+-%f$\n",chi2m[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
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
   
   en_tot=0;
   for (n=0;n<N;n++){
       en[n]=ensemble_omega_jacob;
       en_tot+=en[n];
   }
   
   Nvar=3;
   Npar=2;
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

for (ms=0;ms<1;ms++){
    if(ms==0)   guess[0]=0.125214 ;guess[1]=-1.614364 ;  
   

    count=0;
    for (n=0;n<N;n++){
       // printf("#function %d\n",n);
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        double ref_MK_MPI_MOMEGA=r1->MKMeV[j]*r1->MKMeV[j];
                        ref_MK_MPI_MOMEGA-=(r1->MpiMeV[j]*r1->MpiMeV[j]/2);
                        ref_MK_MPI_MOMEGA/=r1->MOmegaMeV[j]*r1->MOmegaMeV[j];
                        
                        rm[j]=r[0][e][j]+ref_MK_MPI_MOMEGA*r[1][e][j];
                        //rm[j]=r[0][e][j]+mref[ms]*r[1][e][j]+mref[ms]*mref[ms]*r[2][e][j];
                    }
                    fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                for (j=0;j<jack_tot;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1] ;
                
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
    
                double xi=gJ[e].M_PS_jack[0][j]/(4*pi_greco*gJ[e].f_PS_jack[0][j]);
                xi=xi*xi;
                double delta=FVE_GL_Mpi(head[e].l1 , xi,   gjack[e].f_PS_jack[0][j]     );
                x[e+count][0]=gJ[e].M_PS_jack[0][j]*gJ[e].w0[j] /(1-0.25*delta);//M_Pi*w0
                //x[e+count][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1];//ml*w0
                x[e+count][1]=sqrt(mref[ms]);//MK*w0
                x[e+count][2]=gJ[e].w0[j];//w0
                //x[e+count][2]=w0e[e][j];//w0
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
        if (Npar==3){
          tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, chiral_continuum_M_Omega,guess );
          chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, chiral_continuum_M_Omega  );
        }
        else if(Npar==2){
            tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, chiral_continuum_M_Omega_no_a,guess );
            chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, chiral_continuum_M_Omega_no_a  );
        }
        else if(Npar==1){
            tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, chiral_continuum_M_Omega_const,guess );
            chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, chiral_continuum_M_Omega_const  );
        }
         if(j==Njack-1){
            printf("\n\n");
            double Kf;
            for(i=0;i<Npar;i++)   
                printf("P%d=%f;  ",i,tmp[i]);
            printf("  MKw2-Mpi22/MOMEGA2=%f ish;\n",0.084068872956 );
            chi2m=mean_and_error(jack_files[0].sampling,Njack, chi2);
            printf("$\\chi^2/dof=%f+-%f$\n",chi2[0]/(en_tot-Npar),chi2m[1]/(en_tot-Npar));
            free(chi2m);
              printf("#(M_Pi*w0)^2        M_Omega w0    err  Kf  w0\n");
              for (e=0;e<ensemble_omega_jacob;e++){
                //  FVE_K( r1->Bw[Njack-1], r1->fw[Njack-1], double(head[e].l1)/gJ[e].w0[Njack-1],  x[e][0]*x[e][0]/(2.*r1->Bw[Njack-1])/*mlw*/, (mref[ms])/r1->Bw[Njack-1]-x[e][0]*x[e][0]/(2.*r1->Bw[Njack-1]) /*msw*/ ,x[e][3],  x[e][4],x[e][5], x[e][6],&KM, &Kf);
                  Kf=1;
                  printf("%f        %f       %f    %f     %f\n",x[e][0]*x[e][0],y[Njack-1][e][0] ,y[Njack-1][e][1],Kf,gJ[e].w0[Njack-1]);
              }
         }
         for(i=0;i<Npar;i++)
             par[i][j]=tmp[i];   

       
         free(tmp);
    }
    for(i=0;i<Npar;i++){
        double *ave=mean_and_error(jack_files[0].sampling,Njack, par[i]);
        printf("P%d= %g   %g\t",i,ave[0],ave[1]);
        free(ave);
    }
    printf("\n");
    
    
    
    //extract w0 at given m_ref  after the chiral and continuum fit
    if (Npar>1){
    chi2m=mean_and_error(jack_files[0].sampling,Njack, par[1]);
    if(fabs(chi2m[0])<=chi2m[1]*0.5   ){
        printf("slope is zero\n");
        for (j=0;j<Njack;j++)
            fK[ms][j]=par[0][j]/r1->MOmegaMeV[j] ;
    }
    //div_jackboot(Njack,fK[ms], par[0],r1->MOmegaMeV );
    else {
       for (j=0;j<Njack;j++){
            fK[ms][j]=r1->MOmegaMeV[j]  - sqrt( r1->MOmegaMeV[j]*r1->MOmegaMeV[j]-4.* par[0][j]*par[1][j]*r1->MpiMeV[j]*r1->MpiMeV[j]  );
            fK[ms][j]=fK[ms][j]/(2.*par[1][j]*r1->MpiMeV[j]*r1->MpiMeV[j] );//w0 at that reference ms mass
            //compute Omega from w0
            //fK[ms][j]=par[0][j]+par[1][j]*r1->MpiMeV[j]*r1->MpiMeV[j] *(0.1701*0.1701/(197.3*197.3));
            //fK[ms][j]/=(0.1701/(197.3));
       }
    }
    free(chi2m);
    }
    else{
        for (j=0;j<Njack;j++)
            fK[ms][j]=par[0][j]/r1->MOmegaMeV[j] ;
    }
    
   
   for (e=0;e<en_tot;e++){
        free(x[e]);  free(fit[e]);
   }
}
   
   printf("w0=%f       MeV^-1\n", fK[0][Njack-1]);
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
   
   free_3(ensemble_omega_jacob, 3  ,omega);
   return fK;
    
} 

 
