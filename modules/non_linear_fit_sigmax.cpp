#define non_linear_fit_sigmax_C


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "mutils.hpp"
#include "tower.hpp"

#define MAXIT 1000000

//#include "jacknife.h"
//#include "bootstrap.h"
#include "linear_fit.hpp"
#include "non_linear_fit.hpp"


    
double fun_for_non_linear_fit_sigmax(int n , int e, int N,int N1 ,int Nvar ,double *x,int Npar, int Npar1,double *P, double fun(int,int,double*,int,double*)){
    double r;
    //printf("n=%d \t e=%d Npar+e=%d\n",n,e,Npar+e*Nvar);
    r=fun(n,Nvar,P+Npar+e*Nvar,Npar,P);
    if (n>=N){
        //printf("Npar=%d\n",Npar);
        r=P[Npar+n-N];
    }
    
    return r;
    
}

double compute_chi_non_linear_Nf1( int N, int N1,int *ensemble1,double **x1, double **y1, double *P ,int Nvar, int Npar, int Npar1,  double f1(int,int,int,int,int,double*,int,int,double*,double (*)(int,int,double*,int,double*)), double fun(int,int,double*,int,double*)) {
    double chi2=0,f;
    int e,n,count;
    
    count=0;
    for (n=0;n<N1;n++){
        for (e=0;e<ensemble1[n];e++){
            f=f1(n,e,N,N1,Nvar,x1[count],Npar,Npar1,P,fun)-y1[count][0];
            f/=y1[count][1];
            chi2+=f*f;
            //printf("chi2[%d]=%g    f1=%g      y1=%g     erry1=%g  P[%d]=%g, Npar=%d n=%d N=%d\n",count,f*f,f1(n,e,N,N1,Nvar,x1[count],Npar,Npar1,P,fun),y1[count][0],y1[count][1],Npar+n-N,P[Npar+n-N],Npar,n,N );
            count++;
        } 
        //count+=ensemble[n];
    }
    return chi2;
}
//https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference 
//derivative 1 accuracy 4
double *der_fun_Nf_h1(int n , int e, int N, int N1, int Nvar, double *x,int Npar, int Npar1 ,double  *P, double f1(int,int,int,int,int,double*,int,int,double*,double (*)(int,int,double*,int,double*)), double fun(int,int,double*,int,double*), double h){
    
    double *Ph=(double*) malloc(sizeof(double)*Npar1);
    double *df=(double*) calloc(Npar1,sizeof(double));
    int i;
    
    for (i=0;i<Npar1;i++)
        Ph[i]=P[i];
    
    for (i=0;i<Npar1;i++){
        Ph[i]=P[i]-2.*h;
        df[i]=f1( n ,  e,  N ,N1, Nvar , x, Npar, Npar1, Ph,  fun);
    
        Ph[i]=P[i]-h;
        df[i]-=8*f1( n ,  e,  N , N1, Nvar , x, Npar, Npar1, Ph,  fun);
        
        Ph[i]=P[i]+h;
        df[i]+=8*f1( n ,  e,  N , N1, Nvar , x, Npar, Npar1, Ph,  fun);
    
        Ph[i]=P[i]+2.*h;
        df[i]-=f1( n ,  e,  N , N1 , Nvar , x, Npar, Npar1, Ph,  fun);
    
        Ph[i]=P[i];//you need to leave the parameter as it was before you move to the next parameter
        df[i]/=(12.*h);
    }
    
    
    free(Ph);
    
    return df;
}

// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(index_function,Nvariables,variables[], Nparameters,parameters[])
//sigmax [ensemble][variable number]
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2 
double  *non_linear_fit_Nf_sigmax(int N, int *ensemble ,double **x, double **sigmax, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess ){
  
    int i,e;
    
    int N1;    int *ensemble1;     double **x1;    double **sigmxax1; double **y1; int Nvar1; int Npar1;  ;double *guess1;
    double *P;
    int en_tot=0;
    double (*f1)(int,int,int,int,int,double*,int,int,double*, double (*)(int,int,double*,int,double*)  );
    
    f1=fun_for_non_linear_fit_sigmax;

    
    N1=N+Nvar*ensemble[0];
    Npar1=Npar+ensemble[0]*Nvar;
    Nvar1=Nvar;
    ensemble1=(int*) malloc(sizeof(int)*N1);
    for (i=0;i<N;i++){
        ensemble1[i]=ensemble[i];
        en_tot+=ensemble[i];
    }
    for (i=N;i<N1;i++){
        ensemble1[i]=1;
        en_tot+=1;
    }
    
    x1=double_malloc_2(en_tot,Nvar1);
    y1=double_malloc_2(en_tot,2);
    
    int count=0;
    for (i=0;i<N;i++){
        for (e=0;e<ensemble1[i];e++){
            y1[count][0]=y[count][0];
            y1[count][1]=y[count][1];
            count++;
        }
    }
    for (i=N;i<N1;i++){
        for (e=0;e<ensemble1[i];e++){
            //printf("en_tot=%d  \tcount=%d     %d     %d\n",en_tot,count, (i-N)/Nvar, (i-N)%Nvar);  //j-N= v+e *Nvar;
            y1[count][0]=x[ (i-N)/Nvar ][  (i-N)%Nvar  ];
            y1[count][1]=sigmax[ (i-N)/Nvar ][  (i-N)%Nvar  ];
            count++;
        }
    }
    guess1=(double*) malloc(sizeof(double)*Npar1);
    for (i=0;i<Npar;i++){
        guess1[i]=guess[i];
    }
    count=0;
    for (i=0;i<Nvar;i++){
        for (e=0;e<ensemble1[i];e++){
            guess1[count+Npar]=x[e][i];
            count++;
        }
    }

    double **alpha,*X,*beta,**a,**C,*sigma;
    int j,k;
    double f,*fk;
    double chi2,chi2_tmp;
    double *P_tmp,lambda,res;
    int n,Niter=0;
    double h=0.00001;
    lambda=0.001;
    res=1;
    
    P=(double*) malloc(Npar1*sizeof(double));
    P_tmp=(double*) malloc(Npar1*sizeof(double));
   
   // printf("Npar %d  -> %d\n",Npar,Npar1);
   // printf("N %d  -> %d\n",N,N1);
    
    for (j=0;j<Npar;j++){
        P[j]=guess[j];
        P_tmp[j]=P[j];
    }
    for (j=Npar;j<Npar1;j++){
        int e=(j-Npar)/Nvar;    //j-Npar= v+e *Nvar;
        int v=j-Npar-e*Nvar;
        //P[j]=x[j-Npar][0];
        //P[j]=x[e][v];
        P[j]=x[ (j-Npar)/Nvar ][  (j-Npar)%Nvar  ];
        P_tmp[j]=P[j];
        //printf("P[%d]=%g\t",j,P[j]);
    }
   //printf("\n");
    
    beta=(double*) calloc(Npar1,sizeof(double));
    alpha=(double**) malloc(sizeof(double*)*Npar1);
    for (j=0;j<Npar1;j++){
        alpha[j]=(double*) calloc(Npar1,sizeof(double));
    }   

    chi2=compute_chi_non_linear_Nf1(N,N1, ensemble1,x1, y1, P_tmp ,Nvar1,Npar , Npar1,  f1,fun);
//printf("chi2 with extra parameter=%g\n",chi2);
    chi2_tmp=chi2+1;
   
    while (res>0.001){
        chi2_tmp=chi2+1;  
        if(Niter>200){ printf("Niter=%d of the Levenberg-Marquardt chi2 minimization: exeeds max number\n",Niter); break;}
        Niter++;
    while (chi2_tmp>chi2){
          
            count=0;
            for (n=0;n<N1;n++){
                for (e=0;e<ensemble1[n];e++){
                    f=f1( n ,  e, N, N1 , Nvar1 , x1[e],Npar, Npar1, P,  fun);
                    fk=der_fun_Nf_h1(n,e,N,N1,  Nvar1, x1[e],Npar, Npar1,P,  f1,fun,  h);
                    /*for (j=0;j<Npar1;j++)
                        printf("f_%d=%g\t",j,fk[j]);
                    printf("\n");*/
                    for (j=0;j<Npar1;j++){
                        beta[j]+=(y1[e+count][0]-f)*fk[j]/(y1[e+count][1]*y1[e+count][1]);
                        for (k=j;k<Npar1;k++){
                            alpha[j][k]+=fk[j]*fk[k]/(y1[e+count][1]*y1[e+count][1]);
                        }
                        
                        
                    }
                    free(fk);
                }
                count+=ensemble1[n];
            }
                
            for (j=0;j<Npar1;j++){
                alpha[j][j]*=(lambda+1.);
                for (k=0;k<j;k++)
                    alpha[j][k]=alpha[k][j];
                
            
            }

            free(P_tmp);
            P_tmp=cholesky_solver_if_possible(Npar1 , alpha , beta);
            for (j=0;j<Npar1;j++)
                P_tmp[j]+=P[j];
           /*for (j=0;j<Npar1;j++)
                    printf("P[%d]=%g\t",j,P_tmp[j]);
            printf("\n");*/
            chi2_tmp=compute_chi_non_linear_Nf1(N, N1, ensemble1,x1, y1, P_tmp , Nvar ,Npar,  Npar1,  f1,fun);
	        if (chi2_tmp!=chi2_tmp) chi2_tmp=chi2+1;		    

            if (chi2_tmp>chi2)
                lambda*=10;
            
            for (j=0;j<Npar1;j++){
                beta[j]=0;
                for (k=0;k<Npar1;k++){
                       alpha[j][k]=0;
                }
            }
            
            //printf("chi2=%g \t chi2_tmp=%g  Npar=%d\n\n",chi2,chi2_tmp,Npar);
            if(lambda>1e+15){
                printf("lambda of the Levenberg-Marquardt chi2 minimization: exeeds 1e+15 lambda=%g\n RESET lambda=0.001\n",lambda); 
                   lambda=0.001;
            }

            
        }
        res=chi2-chi2_tmp;
        chi2=chi2_tmp;
        lambda/=10;
        for (j=0;j<Npar1;j++){
            P[j]=P_tmp[j];
        }
    }

    for (j=0;j<Npar1;j++){
        free(alpha[j]);
    }
    free(P_tmp);
    free(alpha);free(beta);
    return P;  
    
    
    
    
    //end
    free_2(en_tot,x1);
    free_2(en_tot,y1);
    free(guess1);
    free(ensemble1);
    return P;
}

double chi_parameters_square(int Npar, double *P, double *P1, double **C){
    double r=0;
    for (int i=0;i<Npar;i++)
        r+=(P[i]-P1[i])*(P[i]-P1[i])/ (C[i][i]);
    
    return r;
    
}
 
double  *non_linear_fit_Nf_sigmax_iterative(int N, int *ensemble ,double **x, double **sigmax, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess ){ 
 
    int en_tot=0;
    for (int n=0;n<N;n++)
        for (int e=0;e<ensemble[n];e++)
            en_tot++;
    double **y_tmp=double_malloc_2(en_tot,2);
    double *P1;

    double h=0.0001;
    double *P=non_linear_fit_Nf( N, ensemble ,x/*, sigmax*/, y , Nvar,  Npar,   fun , guess ).P;  //first fit 
    double  **C=covariance_non_linear_fit_Nf(N, ensemble,x ,y ,P , Nvar,  Npar, fun );        
    double  chiP2=1;   
    while(chiP2> 0.001) {
        int count=0;
        for(int n=0; n<N; n++){
            for (int e=0;e<ensemble[n];e++){
                double *df=derN_fun_Nf_var_h( n,  Nvar, x[e] ,Npar, P,  fun, h, 1);  //derivative respect all the variable, computed at the value of the ensemble e
                y_tmp[count][0]=y[count][0];  //this can be done onece at the beginning
                
                double dfx=0;
                for(int v=0; v<Nvar; v++){
                    dfx+=df[v]*sigmax[e][v];
                    dfx*=dfx;
                   // printf("df[%d]=%g\t",v,df[v]);
                }
                //printf("\n");
                
                y_tmp[count][1]=y[count][1]*y[count][1]+dfx;
                y_tmp[count][1]=sqrt(y_tmp[count][1]);
                count++;
                free(df);
            }
        }
        P1=non_linear_fit_Nf( N, ensemble ,x/*, sigmax*/, y_tmp , Nvar,  Npar,   fun , guess ).P;
        
        chiP2=chi_parameters_square( Npar, P, P1,C);
       /* printf(" old:\n");
        for(int ip=0;ip<Npar;ip++)
            printf("P[%d]=%.5e\t",ip,P[ip]);
        printf("\n new:\n");
        for(int ip=0;ip<Npar;ip++)
            printf("P[%d]=%.5e\t",ip,P1[ip]);
        printf("\n difference chiP2=%g\n\n",chiP2 );
        */
        free(P);
        P=P1;
     
    }
    
    free_2(Npar,C);
    free_2(en_tot,y_tmp);
    
    return P1;
}
 
double compute_chi_non_linear_Nf1_cov1( int N, int N1,int *ensemble1,double **x1, double **y1, double *P ,int Nvar, int Npar, int Npar1,  double f1(int,int,int,int,int,double*,int,int,double*,double (*)(int,int,double*,int,double*)), double fun(int,int,double*,int,double*) ,double **cov1) {
    double chi2=0,f;
    int e,n,count;
    
    
    int en_tot=0;
    for (n=0;n<N1;n++)
        for (e=0;e<ensemble1[n];e++)
            en_tot++;
    
    double *tmp=(double*) malloc(sizeof(double)*en_tot);     
    count=0;
    for (n=0;n<N1;n++){
        for (e=0;e<ensemble1[n];e++){
            tmp[count]=   f1(n,e,N,N1,Nvar,x1[count],Npar,Npar1,P,fun)-y1[count][0];
            count++;
        }
    }
    
    for (int i=0;i<en_tot;i++)
        chi2+=tmp[i]*cov1[i][i]*tmp[i];
            
    
    for (int i=0;i<en_tot;i++)
        for (int j=i+1;j<en_tot;j++)
            chi2+=2*tmp[i]*cov1[i][j]*tmp[j];
        
    free(tmp);    
    
    
    
    /*
    count=0;
    for (n=0;n<N1;n++){
        for (e=0;e<ensemble1[n];e++){
            f=f1(n,e,N,N1,Nvar,x1[count],Npar,Npar1,P,fun)-y1[count][0];
            int count1=0;
            for (int n1=0;n1<N1;n1++){
                for (int e1=0;e1<ensemble1[n1];e1++){
                                double fr=f1(n1,e1,N,N1,Nvar,x1[count1],Npar,Npar1,P,fun)-y1[count1][0];
                                chi2+=f *cov1[count][count1]*fr;
                                count1++;
                }
            }
            count++;
        } 
    }
    */
    return chi2;
} 
 
 
// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(index_function,Nvariables,variables[], Nparameters,parameters[])
//sigmax [ensemble][variable number]
//cov1 [ent_tot+ Nvar* ensemble][ent_tot+ Nvar* ensemble]= covariance ^-1
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2 
double  *non_linear_fit_Nf_sigmax_covariance(int N, int *ensemble ,double **x, double **sigmax, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess, double **cov1 ){
  
    int i,e;
    
    int N1;    int *ensemble1;     double **x1;    double **sigmxax1; double **y1; int Nvar1; int Npar1;  ;double *guess1;
    double *P;
    int en_tot=0;
    double (*f1)(int,int,int,int,int,double*,int,int,double*, double (*)(int,int,double*,int,double*)  );
    
    f1=fun_for_non_linear_fit_sigmax;

    
    N1=N+Nvar*ensemble[0];
    Npar1=Npar+ensemble[0]*Nvar;
    Nvar1=Nvar;
    ensemble1=(int*) malloc(sizeof(int)*N1);
    for (i=0;i<N;i++){
        ensemble1[i]=ensemble[i];
        en_tot+=ensemble[i];
    }
    for (i=N;i<N1;i++){
        ensemble1[i]=1;
        en_tot+=1;
    }
    
    x1=double_malloc_2(en_tot,Nvar1);
    y1=double_malloc_2(en_tot,2);
    
    int count=0;
    for (i=0;i<N;i++){
        for (e=0;e<ensemble1[i];e++){
            y1[count][0]=y[count][0];
            y1[count][1]=y[count][1];
            count++;
        }
    }
    for (i=N;i<N1;i++){
        for (e=0;e<ensemble1[i];e++){
            //printf("en_tot=%d  \tcount=%d     %d     %d\n",en_tot,count, (i-N)/Nvar, (i-N)%Nvar);  //j-N= v+e *Nvar;//[e=(i-N)/Nvar][v=(i-N)%Nvar]
            y1[count][0]=x[ (i-N)/Nvar ][  (i-N)%Nvar  ];
            y1[count][1]=sigmax[ (i-N)/Nvar ][  (i-N)%Nvar  ];
            count++;
        }
    }
    guess1=(double*) malloc(sizeof(double)*Npar1);
    for (i=0;i<Npar;i++){
        guess1[i]=guess[i];
    }
    count=0;
    for (i=0;i<Nvar;i++){
        for (e=0;e<ensemble1[i];e++){
            guess1[count+Npar]=x[e][i];
            count++;
        }
    }
 /*   
    int yn=is_it_positive( cov,  en_tot);
    while(yn==1){
        printf("covariance matrix not positive defined adding 0.0001*cov[0][0]*I \n");
        for(i=0;i<en_tot;i++)
             cov[i][i]+=cov[0][0]*1e-12;
        yn=is_it_positive( cov,  en_tot);  
        printf("now the matrix is positive defined.  %d\n",yn);
    }  
        
    double **cov1=symmetric_matrix_inverse(en_tot, cov  );
*/
    double **alpha,*X,*beta,**a,**C,*sigma;
    int j,k;
    double f,*fk;
    double chi2,chi2_tmp;
    double *P_tmp,lambda,res;
    int n,Niter=0;
    double h=0.00001;
    lambda=0.001;
    res=1;
    
    P=(double*) malloc(Npar1*sizeof(double));
    P_tmp=(double*) malloc(Npar1*sizeof(double));
   
   // printf("Npar %d  -> %d\n",Npar,Npar1);
   // printf("N %d  -> %d\n",N,N1);
    
    for (j=0;j<Npar;j++){
        P[j]=guess[j];
        P_tmp[j]=P[j];
    }
    
    for (j=Npar;j<Npar1;j++){
        int e=(j-Npar)/Nvar;    //j-Npar= v+e *Nvar;
        int v=j-Npar-e*Nvar;
        //P[j]=x[j-Npar][0];
        //P[j]=x[e][v];
        P[j]=x[ (j-Npar)/Nvar ][  (j-Npar)%Nvar  ];
        P_tmp[j]=P[j];
        //printf("P[%d]=%g\t",j,P[j]);
    }
    
   //printf("\n");
    
    beta=(double*) calloc(Npar1,sizeof(double));
    alpha=(double**) malloc(sizeof(double*)*Npar1);
    for (j=0;j<Npar1;j++){
        alpha[j]=(double*) calloc(Npar1,sizeof(double));
    }   

    chi2=compute_chi_non_linear_Nf1_cov1(N,N1, ensemble1,x1, y1, P_tmp ,Nvar1,Npar , Npar1,  f1,fun,cov1);
//printf("chi2 with extra parameter=%g\n",chi2);
    chi2_tmp=chi2+1;
   
    
    double *f_value=(double*) malloc(sizeof(double)*en_tot);
    double **df_value=(double**) malloc(sizeof(double*)*en_tot);
    
    while (res>0.001){
        chi2_tmp=chi2+1;  
        if(Niter>1000){ printf("Niter=%d of the Levenberg-Marquardt chi2 minimization: exeeds max number\n",Niter); break;}
        Niter++;//printf("Niter=%d chi2=%g\n",Niter,chi2);
    while (chi2_tmp>chi2){
          
            count=0;
            for (n=0;n<N1;n++){
                for (e=0;e<ensemble1[n];e++){
                    f_value[count]=f1( n ,  e, N, N1 , Nvar1 , x1[e],Npar, Npar1, P,  fun);
                    df_value[count]=der_fun_Nf_h1(n,e,N,N1,  Nvar1, x1[e],Npar, Npar1,P,  f1,fun,  h);
                    count++;
                }
            }
            for (j=0;j<Npar1;j++){
                //for(i=0;i<en_tot;i++)   
                //    beta[j]+=(y1[i][0]-f_value[i])*cov1[i][i]*df_value[i][j];
                for(i=0;i<en_tot;i++)   
                    for(int ii=0;ii<en_tot;ii++) 
                        beta[j]+=(y1[i][0]-f_value[i])*cov1[i][ii]*df_value[ii][j];
                
                
                for (k=j;k<Npar1;k++){
                    //for(i=0;i<en_tot;i++)
                    //    alpha[j][k]+=df_value[i][j]*cov1[i][i]*df_value[i][k];
                    for(i=0;i<en_tot;i++)   
                        for(int ii=0;ii<en_tot;ii++) 
                            alpha[j][k]+=df_value[i][j]*cov1[i][ii]*df_value[ii][k];            
                }
            }
            for(i=0;i<en_tot;i++) 
                free(df_value[i]);
               
        
        
             /*
            count=0;
            for (n=0;n<N1;n++){
                for (e=0;e<ensemble1[n];e++){
                    f=f1( n ,  e, N, N1 , Nvar1 , x1[e],Npar, Npar1, P,  fun);
                    fk=der_fun_Nf_h1(n,e,N,N1,  Nvar1, x1[e],Npar, Npar1,P,  f1,fun,  h);
                    int count1=0;
                    for (int n1=0;n1<N1;n1++){
                        for (int e1=0;e1<ensemble1[n1];e1++){
                            double fr=f1( n1 ,  e1, N, N1 , Nvar1 , x1[e1],Npar, Npar1, P,  fun);
                            double *fk1=der_fun_Nf_h1(n,e1,N,N1,  Nvar1, x1[e1],Npar, Npar1,P,  f1,fun,  h);   
                            for (j=0;j<Npar1;j++){
                                beta[j]+=(y1[e1+count1][0]-fr)*cov1[e1+count1][e+count]*fk[j];
                                for (k=j;k<Npar1;k++){
                                    alpha[j][k]+=fk1[j]*cov1[e1+count1][e+count]*fk[k];
                                }
                            }
                            free(fk1);
                        }
                        count1+=ensemble1[n1];
                    }
                    
                    free(fk);
                }
                count+=ensemble1[n];
            }
             */   
            for (j=0;j<Npar1;j++){
                alpha[j][j]*=(lambda+1.);
                for (k=0;k<j;k++)
                    alpha[j][k]=alpha[k][j];
                
            
            }

            free(P_tmp);
            P_tmp=cholesky_solver_if_possible(Npar1 , alpha , beta);
            for (j=0;j<Npar1;j++)
                P_tmp[j]+=P[j];
           /*for (j=0;j<Npar1;j++)
                    printf("P[%d]=%g\t",j,P_tmp[j]);
            printf("\n");*/
            chi2_tmp=compute_chi_non_linear_Nf1_cov1(N, N1, ensemble1,x1, y1, P_tmp , Nvar ,Npar,  Npar1,  f1,fun,cov1);
	        if (chi2_tmp!=chi2_tmp) chi2_tmp=chi2+1;		    

            if (chi2_tmp>chi2)
                lambda*=10;
            
            for (j=0;j<Npar1;j++){
                beta[j]=0;
                for (k=0;k<Npar1;k++){
                       alpha[j][k]=0;
                }
            }
            
            //printf("chi2=%g \t chi2_tmp=%g  Npar=%d\n\n",chi2,chi2_tmp,Npar);
            if(lambda>1e+15){
                printf("lambda of the Levenberg-Marquardt chi2 minimization: exeeds 1e+15 lambda=%g\n RESET lambda=0.001\n",lambda); 
                   lambda=0.001;
                for(j=0;j<Npar1;j++)
                    P[j]=guess1[j];
            }

            
        }
        res=chi2-chi2_tmp;
        chi2=chi2_tmp;
        lambda/=10;
        for (j=0;j<Npar1;j++){
            P[j]=P_tmp[j];
        }
    }

    for (j=0;j<Npar1;j++){
        free(alpha[j]);
    }
    free(P_tmp);
    free(alpha);free(beta);
    //return P;  
    printf("fit converged after %d iterations",Niter);
    
    //free_2(en_tot,cov1);
    
    //end
    free(f_value);free(df_value);
    free_2(en_tot,x1);
    free_2(en_tot,y1);
    free(guess1);
    free(ensemble1);
    return P;
}
 
