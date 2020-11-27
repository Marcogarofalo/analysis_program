#define non_linear_fit_C


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "mutils.hpp"

#define MAXIT 1000000

//#include "jacknife.h"
//#include "bootstrap.h"
#include "linear_fit.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"
 /*
double *LU_decomposition_solver(int N, double **M, double *b){
 double **U,**L,*y,*x;
 int i,j,k;
 
 x=(double*) malloc(sizeof(double)*N);
 y=(double*) malloc(sizeof(double)*N);
 L=(double**) malloc(sizeof(double*)*N);
 U=(double**) malloc(sizeof(double*)*N);
 for (i=0;i<N;i++){
    U[i]=(double*) calloc(N,sizeof(double));
    L[i]=(double*) calloc(N,sizeof(double));
 }
 
 for (i=0;i<N;i++)
     L[i][i]=1;
     
 for (j=0;j<N;j++){
    for (i=0;i<=j;i++){
        U[i][j]=M[i][j];
        for (k=0;k<i;k++)
            U[i][j]-=L[i][k]*U[k][j];
    }
    for (i=j+1;i<N;i++){
        L[i][j]=M[i][j];
        for (k=0;k<j;k++)
            L[i][j]-=L[i][k]*U[k][j];
        L[i][j]/=U[j][j];
    }
 }
 
 y[0]=b[0]/L[0][0];
 for (i=0;i<N;i++){
     y[i]=b[i];
     for (k=0;k<i;k++)
         y[i]-=L[i][k]*y[k];
     y[i]/=L[i][i];
 }
 
 x[N-1]=y[N-1]/U[N-1][N-1];
 for (i=N-2;i>=0;i--){
     x[i]=y[i];
     for (k=i+1;k<N;k++)
         x[i]-=U[i][k]*x[k];
     x[i]/=U[i][i];
 }
 
 free(y);
 for (i=0;i<N;i++){
     free(L[i]);
     free(U[i]);
 }
 free(U);
 free(L);
 
 return x;
  
}

//return the inverse matrix of M
double **matrix_inverse(int N, double **M  ){
    double *b,**r,*a;
    int i,j;
    
    b=(double*) calloc(N,sizeof(double));
    r=(double**) malloc(sizeof(double*)*N);
    for (i=0;i<N;i++){
        r[i]=(double*) malloc(N*sizeof(double));
    }
    
    for (i=0;i<N;i++){
        b[i]=1.;
        a=LU_decomposition_solver(N, M, b);
        for (j=0;j<N;j++)
            r[j][i]=a[j];
        free(a);
        b[i]=0;
    }
    
    free(b);
        
    return r;
}*/


//https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference 
//derivative 1 accuracy 4
double *der_fun_h( int Nvar, double *x,int Npar,double  *P, double fun(int,double*,int,double*), double h){
    
    double *Ph=(double*) malloc(sizeof(double)*Npar);
    double *df=(double*) calloc(Npar,sizeof(double));
    int i;
    
    for (i=0;i<Npar;i++)
        Ph[i]=P[i];
    
    for (i=0;i<Npar;i++){
        Ph[i]=P[i]-2.*h;
        df[i]=fun(Nvar,x,Npar,Ph);
    
        Ph[i]=P[i]-h;
        df[i]-=8*fun(Nvar,x,Npar,Ph);
        
        Ph[i]=P[i]+h;
        df[i]+=8*fun(Nvar,x,Npar,Ph);
    
        Ph[i]=P[i]+2.*h;
        df[i]-=fun(Nvar,x,Npar,Ph);
    
        Ph[i]=P[i];//you need to leave the parameter as it was before you move to the next parameter
        df[i]/=(12.*h);
    }
    
    
    free(Ph);
    
    return df;
}

 
 
double compute_chi_non_linear(int ensemble,double **x, double **y, double *P ,int Nvar, int Npar, double fun(int,double*,int,double*)) {
    double chi2=0,f;
    int e,j;
    
    for (e=0;e<ensemble;e++){
        f=fun(Nvar,x[e],Npar,P);
        chi2+=(y[e][0]-f)*(y[e][0]-f)/ (y[e][1]*y[e][1]);
     
    } 
    
    return chi2;
    
}
 
// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(Nvariables,variables[], Nparameters,parameters[])
//*funk(Nvariables,variables[], Nparameters,parameters[]) must return a vector[Nparameters] that contain the derivatives of fun respect to the parameters
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2 
//double  *non_linear_fit(int ensemble,double **x, double **y ,int Nvar, int Npar, double fun(int,double*,int,double*) , double *funk(int,double*,int,double*) ){
double  *non_linear_fit(int ensemble,double **x, double **y ,int Nvar, int Npar, double fun(int,double*,int,double*)  ){
  
    double **alpha,*X,*beta,**a,**C,*sigma;
    int i,j,k,e;
    double f,*fk;
    double chi2,chi2_tmp;
    double *P,*P_tmp,lambda,res;
    double h=0.00001;
    lambda=0.001;
    res=1;
    
 //   double *fkk;
    P=(double*) malloc(Npar*sizeof(double));
    P_tmp=(double*) malloc(Npar*sizeof(double));
   
   
    for (j=0;j<Npar;j++){
      //  P[j]=0.1;
        P_tmp[j]=P[j];
    }
    
    beta=(double*) calloc(Npar,sizeof(double));
    alpha=(double**) malloc(sizeof(double*)*Npar);
    for (j=0;j<Npar;j++){
        alpha[j]=(double*) calloc(Npar,sizeof(double));
    }   
    chi2=compute_chi_non_linear( ensemble,x, y, P ,Nvar,  Npar,  fun);

   
    // printf("chi2=%f   res=%.10f Bw=%f   fw=%f\n",chi2,res,P[0],P[1]);
    while (res>0.001){
        chi2_tmp=chi2+1;   

        while (chi2_tmp>chi2){
          
            for (e=0;e<ensemble;e++){
                f=fun(Nvar,x[e],Npar,P);
                fk=der_fun_h(  Nvar, x[e], Npar,P,  fun,  h);
             //   fk=funk(Nvar,x[e],Npar,P);
                for (j=0;j<Npar;j++){
                    beta[j]+=(y[e][0]-f)*fk[j]/(y[e][1]*y[e][1]);
              //      printf("|analitic-numeric|=  |%g -%g|   = %g\n",fk[j],fkk[j],fabs(fk[j]-fkk[j]));
                    for (k=0;k<Npar;k++){
                        alpha[j][k]+=fk[j]*fk[k]/(y[e][1]*y[e][1]);
                    }
                }
                free(fk);
            }    

            for (j=0;j<Npar;j++){
                alpha[j][j]*=(lambda+1.);
            }

            if (Npar==1){
                C=(double**) malloc(sizeof(double*)*1);
                C[0]=(double*) malloc(sizeof(double)*1);
                C[0][0]=1./alpha[0][0];
            }
            else
                C=matrix_inverse(Npar, alpha  );
            
            for (j=0;j<Npar;j++)
                P_tmp[j]=P[j];
            for (j=0;j<Npar;j++){
                for (k=0;k<Npar;k++){
                     P_tmp[j]+=C[j][k]*beta[k];
                 }
            }    

            //printf("lambda=%f\n",lambda);
            chi2_tmp=compute_chi_non_linear( ensemble,x, y, P_tmp ,Nvar,  Npar,  fun);
	        if (chi2_tmp!=chi2_tmp) chi2_tmp=chi2+1;
             //       printf("chi2=%f chi2_tmp=%f  res=%f P0=%f P0_tmp=%f  P0=%f P0_tmp=%f\n",chi2,chi2_tmp,res,P[0],P_tmp[0],P[1],P_tmp[1]);
		    

            if (chi2_tmp>chi2)
                lambda*=10;
            for (j=0;j<Npar;j++){
                free(C[j]);
                 beta[j]=0;
                 for (k=0;k<Npar;k++){
                        alpha[j][k]=0;
                 }
            }
            free(C);

            
        }
        res=chi2-chi2_tmp;
        //error(res<0,2,"non_linear_fit","The Levenberg-Marquardt accepted a configuration when the chi2 increased");
        chi2=chi2_tmp;
        lambda/=10;
        for (j=0;j<Npar;j++){
            P[j]=P_tmp[j];
        }
        // printf("chi2=%f   res=%.10f Bw=%f   fw=%f  P1=%f  P2=%f\n",chi2,res,P[0],P[1],P[2],P[3]);
        //printf("check residue=%f\n",res);
    }


    for (j=0;j<Npar;j++){
        free(alpha[j]);
    }
    free(P_tmp);
    free(alpha);free(beta);
    return P;
}



//https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference 
//derivative 1 accuracy 4
double *der_fun_Nf_h(int n, int Nvar, double *x,int Npar,double  *P, double fun(int,int,double*,int,double*), double h){
    
    double *Ph=(double*) malloc(sizeof(double)*Npar);
    double *df=(double*) calloc(Npar,sizeof(double));
    int i;
    
    for (i=0;i<Npar;i++)
        Ph[i]=P[i];
    
    for (i=0;i<Npar;i++){
        Ph[i]=P[i]-2.*h;
        df[i]=fun(n,Nvar,x,Npar,Ph);
    
        Ph[i]=P[i]-h;
        df[i]-=8*fun(n,Nvar,x,Npar,Ph);
        
        Ph[i]=P[i]+h;
        df[i]+=8*fun(n,Nvar,x,Npar,Ph);
    
        Ph[i]=P[i]+2.*h;
        df[i]-=fun(n,Nvar,x,Npar,Ph);
    
        Ph[i]=P[i];//you need to leave the parameter as it was before you move to the next parameter
        df[i]/=(12.*h);
    }
    
    
    free(Ph);
    
    return df;
}
//https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference 
//derivative 1 accuracy 4
double *derN_fun_Nf_var_h(int n, int Nvar, double *x,int Npar,double  *P, double fun(int,int,double*,int,double*), double h,int N){
    
    double *xh=(double*) malloc(sizeof(double)*Nvar);
    double *df=(double*) calloc(Nvar,sizeof(double));
    double *tmp;
    int i;

    for (i=0;i<Nvar;i++)
        xh[i]=x[i];
    
    if (N==0){
        for (i=0;i<Nvar;i++)
            df[i]=fun(n,Nvar,xh,Npar,P);
    }
    else if (N==1){
        for (i=0;i<Nvar;i++){
            xh[i]=x[i]-2.*h;    
            df[i]=fun(n,Nvar,xh,Npar,P);
            
            xh[i]=x[i]-h;
            df[i]-=8*fun(n,Nvar,xh,Npar,P);

            xh[i]=x[i]+h;
            df[i]+=8*fun(n,Nvar,xh,Npar,P);
        
            xh[i]=x[i]+2.*h;
            df[i]-=fun(n,Nvar,xh,Npar,P);
        
            xh[i]=x[i];//you need to leave the parameter as it was before you move to the next parameter
            df[i]/=(12.*h);
        }
    }
    else{
        for (i=0;i<Nvar;i++){
            xh[i]=x[i]-2.*h;
            tmp=derN_fun_Nf_var_h( n,  Nvar, xh, Npar,P, fun,  h, N-1);
            df[i]=tmp[i]; free(tmp);
            
            xh[i]=x[i]-h;
            tmp=derN_fun_Nf_var_h( n,  Nvar, xh, Npar,P, fun,  h, N-1);
            df[i]-=8*tmp[i]; free(tmp);
            
            xh[i]=x[i]+h;
            tmp=derN_fun_Nf_var_h( n,  Nvar, xh, Npar,P, fun,  h, N-1);
            df[i]+=8*tmp[i]; free(tmp);
        
            xh[i]=x[i]+2.*h;
            tmp=derN_fun_Nf_var_h( n,  Nvar, xh, Npar,P, fun,  h, N-1);
            df[i]-=tmp[i]; free(tmp);
        
            xh[i]=x[i];//you need to leave the parameter as it was before you move to the next parameter
            df[i]/=(12.*h);
        }
    }
    
    for (i=1;i<=N;i++)
        df[i]/=((double) i);
    free(xh);
    
    return df;
}



double compute_chi_non_linear_Nf(int N,int *ensemble,double **x, double **y, double *P ,int Nvar, int Npar,  double fun(int,int,double*,int,double*)) {
    double chi2=0,f;
    int e,n,count;
    
    count=0;
    for (n=0;n<N;n++){
        for (e=0;e<ensemble[n];e++){
             f=fun(n,Nvar,x[count],Npar,P)-y[count][0];
            f/=y[count][1];
            chi2+=f*f;
            count++;
        } 
    }
    return chi2;
}

// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(index_function,Nvariables,variables[], Nparameters,parameters[])
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2 
double  **covariance_non_linear_fit_Nf(int N, int *ensemble ,double **x, double **y,double *P ,int Nvar, int Npar,  double fun(int,int,double*,int,double*)  ){
  
    double **alpha,**C;
    int i,j,k,e;
    double f,*fk;
    int n,count;
    double h=0.00001;
    
    alpha=(double**) malloc(sizeof(double*)*Npar);
    for (j=0;j<Npar;j++){
        alpha[j]=(double*) calloc(Npar,sizeof(double));
    }   
    
    count=0;
    for (n=0;n<N;n++){
        for (e=0;e<ensemble[n];e++){
            f=fun(n,Nvar,x[e+count],Npar,P);
            fk=der_fun_Nf_h(n,  Nvar, x[e+count], Npar,P,  fun,  h);
            for (j=0;j<Npar;j++){
                    for (k=j;k<Npar;k++){
                    alpha[j][k]+=fk[j]*fk[k]/(y[e+count][1]*y[e+count][1]);
                }
            }
            free(fk);
        }
        count+=ensemble[n];
    }
    for (j=0;j<Npar;j++){
         for (k=0;k<j;k++)
              alpha[j][k]=alpha[k][j];
            
            }
    if (Npar==1){
                C=(double**) malloc(sizeof(double*)*1);
                C[0]=(double*) malloc(sizeof(double)*1);
                C[0][0]=1./alpha[0][0];
    }
    else
                C=symmetric_matrix_inverse(Npar, alpha  );
     
    
    for (j=0;j<Npar;j++){
        free(alpha[j]);
    }
    
    free(alpha);
    return C;
}



// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(index_function,Nvariables,variables[], Nparameters,parameters[])
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2 
double  *non_linear_fit_Nf(int N, int *ensemble ,double **x, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess ){
  
    double **alpha,*X,*beta,**a,**C,*sigma;
    int i,j,k,e;
    double f,*fk;
    double chi2,chi2_tmp;
    double *P,*P_tmp,lambda,res;
    int n,count,Niter=0;
    double h=0.00001;
    int nerror=0;
    lambda=0.001;
    res=1;
    
 //   double *fkk;
    P=(double*) malloc(Npar*sizeof(double));
    P_tmp=(double*) malloc(Npar*sizeof(double));
   
   
    for (j=0;j<Npar;j++){
        P[j]=guess[j];
        P_tmp[j]=P[j];
    }
   
    
    beta=(double*) calloc(Npar,sizeof(double));
    alpha=(double**) malloc(sizeof(double*)*Npar);
    for (j=0;j<Npar;j++){
        alpha[j]=(double*) calloc(Npar,sizeof(double));
    }   

    chi2=compute_chi_non_linear_Nf(N, ensemble,x, y, P_tmp ,Nvar,  Npar,  fun);//printf("chi2 in fit function=%g\n",chi2/(ensemble[0]*N-Npar));

    chi2_tmp=chi2+1;
   
   // printf("chi2=%f   res=%.10f Bw=%f   fw=%f\n",chi2,res,P[0],P[1]);
    while (res>0.001){
        chi2_tmp=chi2+1;  
        if(Niter>200){ printf("Niter=%d of the Levenberg-Marquardt chi2 minimization: exeeds max number\n",Niter); break;}
        Niter++;
        nerror=0;
    while (chi2_tmp-chi2>0) {  //do {} while()   , at least one time is done. if chi is too big chi_tmp=chi+1 = chi 
        //printf("lambda=%g\n",lambda);
            count=0;
            for (n=0;n<N;n++){
                for (e=0;e<ensemble[n];e++){//printf("e=%d   n=%d   en[%d]=%d\n",e,n,n,ensemble[n]);
                    f=fun(n,Nvar,x[e+count],Npar,P);
                    fk=der_fun_Nf_h(n,  Nvar, x[e+count], Npar,P,  fun,  h);
                //   fk=funk(Nvar,x[e],Npar,P);
                    for (j=0;j<Npar;j++){
                        beta[j]+=(y[e+count][0]-f)*fk[j]/(y[e+count][1]*y[e+count][1]);
                //      printf("|analitic-numeric|=  |%g -%g|   = %g\n",fk[j],fkk[j],fabs(fk[j]-fkk[j]));
                        for (k=j;k<Npar;k++){
                            alpha[j][k]+=fk[j]*fk[k]/(y[e+count][1]*y[e+count][1]);
                        }
                        
                        
                    }
                    free(fk);
                }
                count+=ensemble[n];
            }
                
            for (j=0;j<Npar;j++){
                alpha[j][j]*=(lambda+1.);
                for (k=0;k<j;k++)
                    alpha[j][k]=alpha[k][j];
            
            }
/*
            if (Npar==1){
                C=(double**) malloc(sizeof(double*)*1);
                C[0]=(double*) malloc(sizeof(double)*1);
                C[0][0]=1./alpha[0][0];
            }
            else
                C=matrix_inverse(Npar, alpha  );
  
            for (j=0;j<Npar;j++){
                P_tmp[j]=P[j];
                for (k=0;k<Npar;k++){
                     P_tmp[j]+=C[j][k]*beta[k];
                 }
            }    
*/
            free(P_tmp);
            P_tmp=cholesky_solver_if_possible(Npar , alpha , beta);
            //P_tmp=LU_decomposition_solver(Npar , alpha , beta);
            for (j=0;j<Npar;j++)
                P_tmp[j]+=P[j];
            //printf("lambda=%f\n",lambda);
            chi2_tmp=compute_chi_non_linear_Nf(N, ensemble,x, y, P_tmp ,Nvar,  Npar,  fun);
	        if (chi2_tmp!=chi2_tmp) chi2_tmp=chi2+1;
            //       printf("chi2=%f chi2_tmp=%f  res=%f P0=%f P0_tmp=%f  P0=%f P0_tmp=%f\n",chi2,chi2_tmp,res,P[0],P_tmp[0],P[1],P_tmp[1]);
		    

            if (chi2_tmp>chi2)
                lambda*=10;
            
            for (j=0;j<Npar;j++){
                //free(C[j]);
                beta[j]=0;
                for (k=0;k<Npar;k++){
                       alpha[j][k]=0;
                }
            }
           // free(C);
            //error(lambda>1e+15,1,"non_linear_fit_Nf","lambda of the Levenberg-Marquardt chi2 minimization: exeeds 1e+15 lambda=%g",lambda);
            if(lambda>1e+15){
                printf("lambda of the Levenberg-Marquardt chi2 minimization: exeeds 1e+15 lambda=%g\n RESET lambda=0.001\n",lambda); 
                /*free(P_tmp);
                free(alpha);free(beta);
                return P;*/
                lambda=0.001;
                nerror++;
                if (nerror>20){
                    printf("\n !!!!!!!!!!!!!!!! error:  Impossible to minimise the chi2 with Levenberg-Marquardt for this starting point   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
                    for (j=0;j<Npar;j++){
                        free(alpha[j]);
                    }
                    free(P_tmp);
                    free(alpha);free(beta);
                    return P;
                }
                    
            }

            
        } 
        res=chi2-chi2_tmp;
        //error(res<0,2,"non_linear_fit","The Levenberg-Marquardt accepted a configuration when the chi2 increased");
        chi2=chi2_tmp;
        lambda/=10;
        for (j=0;j<Npar;j++){
            P[j]=P_tmp[j];
          //  printf("P[%d]= %f    tot_par=%d\t",j,P[j],Npar);
        }
         //printf("chi2=%f   res=%.10f Bw=%f   fw=%f  P1=%f  P2=%f \n",chi2,res,P[0],P[1],P[2],P[3]);
        //printf("check residue=%f\n",res);
    }

    for (j=0;j<Npar;j++){
        free(alpha[j]);
    }
    free(P_tmp);
    free(alpha);free(beta);
    return P;
}

// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(index_function,Nvariables,variables[], Nparameters,parameters[])
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2 
double  *guess_for_non_linear_fit_Nf(int N, int *ensemble ,double **x, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess ){
  
    int i,j,jmax,k,e,n;
    double en_tot=0;
    double chi2,chi2_tmp,chi2_tmp1;
    double *P,*P_tmp;
    double r,rm,gmax=2;
    rm=(double) RAND_MAX;
    
    for (n=0;n<N;n++)
       en_tot+=ensemble[n];
    //for(i=0;i<Npar;i++)
      //      guess[i]=((r/rm)-0.5)*2;
    P=non_linear_fit_Nf(N, ensemble ,x, y , Nvar,  Npar,   fun ,guess );
    chi2=compute_chi_non_linear_Nf(N, ensemble,x, y, P ,Nvar,  Npar,  fun)/(en_tot-Npar);
    jmax=3+((int) (chi2*2));
    
    if (jmax>15  || jmax <=3) jmax=15;
    
    chi2_tmp1=chi2;

    for (j=0;j<jmax;j++){
        for(i=0;i<Npar;i++){
            r=(double) rand();
            guess[i]=((r/rm)-0.5)*(j+1);
          //  printf("%f\t",guess[i]);
        }
       // printf("\n");
        P_tmp=non_linear_fit_Nf(N, ensemble ,x, y , Nvar,  Npar,   fun ,guess );
        chi2_tmp=compute_chi_non_linear_Nf(N, ensemble,x, y, P_tmp ,Nvar,  Npar,  fun)/(en_tot-Npar);
        printf("chi2=%.10f \tchi2_tmp=%.10f\n",chi2,chi2_tmp);
        /*for(i=0;i<Npar;i++)
            printf("P[%d]=%g\t",i,P[i]);
        printf("\n");*/
        if (fabs(chi2-chi2_tmp)<1e-3 ){
            gmax=gmax*10;//printf("the chi2 didn't change\n\n");
            free(P_tmp);
        }
        else if (chi2_tmp-chi2<-1e-3 ){
            free(P); P=P_tmp;
            chi2=chi2_tmp;
            //printf("chi2 smaller founded\n\n");
        }
        else{
            free(P_tmp);//printf("chi2 LARGER\n\n");
            chi2_tmp1=chi2_tmp;
        }
    }
    free(guess);
    return P;
}



double compute_chi_non_linear_Nf_cov1(int N,int *ensemble,double **x, double **y, double *P ,int Nvar, int Npar,  double fun(int,int,double*,int,double*), double **cov1) {
    double chi2=0,f,f1;
    int e,n,count,e1,n1,count1;
    
    
     
    int en_tot=0;
    for (n=0;n<N;n++)
        for (e=0;e<ensemble[n];e++)
            en_tot++;
    
    double *tmp=(double*) malloc(sizeof(double)*en_tot);     
    count=0;
    for (n=0;n<N;n++){
        for (e=0;e<ensemble[n];e++){
            tmp[count]=   fun(n,Nvar,x[count],Npar,P)-y[count][0];// f1(n,e,N,N1,Nvar,x1[count],Npar,Npar1,P,fun)-y1[count][0];
            count++;
        }
    }
    
    for (int i=0;i<en_tot;i++)
        chi2+=tmp[i]*cov1[i][i]*tmp[i];
            
    
    for (int i=0;i<en_tot;i++)
        for (int j=i+1;j<en_tot;j++)
            chi2+=2.*tmp[i]*cov1[i][j]*tmp[j];
        
    free(tmp);    
    /*double chi2o=0;
    count=0;
    for (n=0;n<N;n++){
        for (e=0;e<ensemble[n];e++){
            count1=0;
            f=fun(n,Nvar,x[count],Npar,P)-y[count][0];
            for (n1=0;n1<N;n1++){
                for (e1=0;e1<ensemble[n1];e1++){
                    f1=fun(n1,Nvar,x[count1],Npar,P)-y[count1][0];
                    chi2o+=f *cov1[count][count1]*f1;
                    count1++;
                }
            }
            count++;
        } 
    }
    error(fabs(chi2-chi2o)>e1-6,1,"chi2 diverso","");
    */
    return chi2;
}

double compute_chi_non_linear_Nf_cov(int N,int *ensemble,double **x, double **y, double *P ,int Nvar, int Npar,  double fun(int,int,double*,int,double*), double **cov) {
    int n,e;
    int en_tot=0;
    for (n=0;n<N;n++)
        for (e=0;e<ensemble[n];e++)
            en_tot+=1;
        
    double **cov1=symmetric_matrix_inverse(en_tot, cov  );
    double chi2=compute_chi_non_linear_Nf_cov1( N, ensemble, x, y,  P , Nvar,  Npar,  fun, cov1) ;
    free_2(en_tot,cov1);
    return chi2;
}
// x[ensemble][variable number] ,   y[ensemble][0=mean,1=error], fun(index_function,Nvariables,variables[], Nparameters,parameters[]), 
//cov[en_tot][en_tot]  is the covariance matrix, with en_tot=sum_i^N ensemble[i],
//the function return an array[Nparameter]  with the value of the parameters that minimise the chi2 
double  *non_linear_fit_Nf_covariance(int N, int *ensemble ,double **x, double **y ,int Nvar, int Npar,  double fun(int,int,double*,int,double*) ,double *guess, double **cov1 ){
  
    double **alpha,*X,*beta,**a,**C,*sigma;
    int i,j,k,e;
    double f,*fk,f1,*fk1;
    double chi2,chi2_tmp;
    double *P,*P_tmp,lambda,res;
    int n,count,n1,count1,e1,Niter=0;
    double h=1.0e-5;
    lambda=0.001;
    res=1;
    int en_tot=0;
    for (n=0;n<N;n++)
        for (e=0;e<ensemble[n];e++)
            en_tot+=1;
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
   
 /*   printf("covariance\n");
    for (j=0;j<en_tot;j++){
        for (i=0;i<en_tot;i++){
        printf("%g\t",cov[j][i]);
        }
        printf("\n");
    }
    printf("inverse\n");
    for (j=0;j<en_tot;j++){
        for (i=0;i<en_tot;i++){
        printf("%g\t",cov1[j][i]);
        }
        printf("\n");
    }*/
        
 //   double *fkk;
    P=(double*) malloc(Npar*sizeof(double));
    P_tmp=(double*) malloc(Npar*sizeof(double));
   
   
    for (j=0;j<Npar;j++){
        P[j]=guess[j];
        P_tmp[j]=P[j];
    }
   
    
    beta=(double*) calloc(Npar,sizeof(double));
    alpha=(double**) malloc(sizeof(double*)*Npar);
    for (j=0;j<Npar;j++){
        alpha[j]=(double*) calloc(Npar,sizeof(double));
    }   

    chi2=compute_chi_non_linear_Nf_cov1(N, ensemble,x, y, P_tmp ,Nvar,  Npar,  fun,cov1);

    chi2_tmp=chi2+1;
    double *f_value=(double*) malloc(sizeof(double)*en_tot);
    double **df_value=(double**) malloc(sizeof(double*)*en_tot);
   
    while (res>0.001){
        chi2_tmp=chi2+1;  
        if(Niter>200){ printf("Niter=%d of the Levenberg-Marquardt chi2 minimization: exeeds max number\n",Niter); break;}
        Niter++;
    while (chi2_tmp>chi2){
        
          count=0;
          for (n=0;n<N;n++){
                for (e=0;e<ensemble[n];e++){
                    f_value[count]=fun(n,Nvar,x[e],Npar,P);
                    df_value[count]=der_fun_Nf_h(n,  Nvar, x[e], Npar,P,  fun,  h);
                    count++;
                }
          }
          for (j=0;j<Npar;j++){
                /*for(i=0;i<en_tot;i++)   
                    beta[j]+=(y[i][0]-f_value[i])*cov1[i][i]*df_value[i][j];*/
                for(i=0;i<en_tot;i++)   
                    for(int ii=0;ii<en_tot;ii++) 
                        beta[j]+=(y[i][0]-f_value[i])*cov1[i][ii]*df_value[ii][j];
                
                
                for (k=j;k<Npar;k++){
                    //for(i=0;i<en_tot;i++)
                    //    alpha[j][k]+=df_value[i][j]*cov1[i][i]*df_value[i][k];
                    for(i=0;i<en_tot;i++)   
                        for(int ii=0;ii<en_tot;ii++) 
                            alpha[j][k]+=df_value[i][j]*cov1[i][ii]*df_value[ii][k];            
                }
          }
          for(i=0;i<en_tot;i++)
              free(df_value[i]);
            
            double *beta1=(double*) calloc(Npar,sizeof(double));
            double **alpha1=(double**) malloc(Npar*sizeof(double));
            for (j=0;j<Npar;j++)
                alpha1[j]=(double*) calloc(Npar,sizeof(double));
                 
            count=0;
            for (n=0;n<N;n++){
                for (e=0;e<ensemble[n];e++){
                    
                    f=fun(n,Nvar,x[e+count],Npar,P);
                    fk=der_fun_Nf_h(n,  Nvar, x[e+count], Npar,P,  fun,  h);
                    count1=0;
                    for (n1=0;n1<N;n1++){
                        for (e1=0;e1<ensemble[n1];e1++){
                            f1=fun(n1,Nvar,x[e1+count1],Npar,P);
                            fk1=der_fun_Nf_h(n1,  Nvar, x[e1+count1], Npar,P,  fun,  h);
                            
                            for (j=0;j<Npar;j++){
                                beta1[j]+=(y[e1+count1][0]-f1)*cov1[e1+count1][e+count]*fk[j];
                                for (k=j;k<Npar;k++){
                                    alpha1[j][k]+=fk1[j]*cov1[e1+count1][e+count]*fk[k];
                                }
                                
                                
                            }
                            free(fk1);
                        }
                        count1+=ensemble[n];
                    }
                    free(fk);
                }
                count+=ensemble[n];
            }
            
            for (j=0;j<Npar;j++)
                error( fabs(beta[j]-beta1[j])>1e-6 ,1,"fit", "beta differs:   %d    %g    %g",j,beta[j],beta1[j]);
            for (j=0;j<Npar;j++)
                for (k=j;k<Npar;k++)
                error( fabs(alpha[j][k]-alpha1[j][k])>1e-6 ,1,"fit", "alpha differs:   %d  %d   %g    %g",j,k,alpha[j][k],alpha1[j][k]);
            
            free(beta1);
                
                
            for (j=0;j<Npar;j++){
                alpha[j][j]*=(lambda+1.);
                for (k=0;k<j;k++)
                    alpha[j][k]=alpha[k][j];
            
            }
            free(P_tmp);
            P_tmp=cholesky_solver_if_possible(Npar , alpha , beta);

            for (j=0;j<Npar;j++)
                P_tmp[j]+=P[j];

            chi2_tmp=compute_chi_non_linear_Nf_cov1(N, ensemble,x, y, P_tmp ,Nvar,  Npar,  fun,cov1);
	        if (chi2_tmp!=chi2_tmp) chi2_tmp=chi2+1;

		    

            if (chi2_tmp>chi2)
                lambda*=10;
            
            for (j=0;j<Npar;j++){

                beta[j]=0;
                for (k=0;k<Npar;k++){
                       alpha[j][k]=0;
                }
            }
            if(lambda>1e+15){
                printf("lambda of the Levenberg-Marquardt chi2 minimization: exeeds 1e+15 lambda=%g\n RESET lambda=0.001\n",lambda); 
                lambda=0.001;
            }

            
        }
        res=chi2-chi2_tmp;

        chi2=chi2_tmp;
        lambda/=10;
        for (j=0;j<Npar;j++){
            P[j]=P_tmp[j];

        }
    }

    //free_2(en_tot,cov1);
    for (j=0;j<Npar;j++){
        free(alpha[j]);
    }
    free(f_value);free(df_value);

    free(P_tmp);
    free(alpha);free(beta);
    return P;
}



double rtsafe(void (*funcd)(double,int,double*, double *, double *),int Npar,double *P , double x1, double x2,
double xacc)
//Using a combination of Newton-Raphson and bisection, find the root of a function bracketed
//between x1 and x2. The root, returned as the function value rtsafe, will be refined until
//its accuracy is known within ±xacc. funcd is a user-supplied routine that returns both the
//function value and the first derivative of the function.
{
int j;
double df,dx,dxold,f,fh,fl;
double temp,xh,xl,rts;

(*funcd)(x1,Npar,P,&fl,&df);
(*funcd)(x2,Npar,P,&fh,&df);
error(((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)),1,"rtsafe","Root must be bracketed in rtsafe");
if (fl == 0.0) return x1;
if (fh == 0.0) return x2;
if (fl < 0.0) {// Orient the search so that f (xl) < 0.
xl=x1;
xh=x2;
} else {
xh=x1;
xl=x2;
}
rts=0.5*(x1+x2);// Initialize the guess for root,
dxold=fabs(x2-x1);  //the “stepsize before last,”
dx=dxold;// and the last step.

(*funcd)(rts,Npar,P,&f,&df); 
for (j=1;j<=MAXIT;j++) { // Loop over allowed iterations.
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)//Bisect if Newton out of range,
        || (fabs(2.0*f) > fabs(dxold*df))) {//or not decreasing fast enough.
        dxold=dx;
        dx=0.5*(xh-xl);
        rts=xl+dx;
        if (xl == rts) return rts;// Change in root is negligible.
    } 
    else {//Newton step acceptable. Take it.
        dxold=dx;
        dx=f/df;
        temp=rts;
        rts -= dx;
        if (temp == rts) return rts;
    }
    if (fabs(dx) < xacc) return rts; //Convergence criterion.

    (*funcd)(rts,Npar,P,&f,&df);
    if (f < 0.0)
        xl=rts;
    else
        xh=rts;
}
error(1>0,1,"rtsafe","Maximum number of iterations exceeded in rtsafe");
return 0.0;
}

double rtbis(double (*func)(double , double,int,double*),double input,int Npar, double *P, double x1, double x2, double xacc)
//Using bisection, find the root of a function func known to lie between x1 and x2. The root,
//returned as rtbis, will be refined until its accuracy is ±xacc.
{

    int j;
    double dx,f,fmid,xmid,rtb;

    f=(*func)(input,x1,Npar,P);
    fmid=(*func)(input,x2,Npar,P);

    error(f*fmid >= 0.0,1,"rtbis","Root must be bracketed for bisection in rtbis f(x1)=%f   f(x2)=%f",f,fmid);
    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);// Orient the search so that f>0
    for (j=1;j<=MAXIT;j++) {// lies at x+dx.
        fmid=(*func)(input,xmid=rtb+(dx *= 0.5),Npar,P); //Bisection loop.
        if (fmid <= 0.0) rtb=xmid;
        if (fabs(dx) < xacc || fmid == 0.0) return rtb;
    }
    error(1>0,1,"rtbis","Too many bisections in rtbis");
    return 0.0;
}


double rtbis_func_eq_input(double (*func)(int , int , double*,int,double*),int n, int Nvar, double *x,int Npar, double *P, int ivar,double input, double x1, double x2, double xacc)
//Using bisection, find the root of a function func-input known to lie between x1 and x2. The root,
//returned as rtbis, will be refined until its accuracy is ±xacc.
//func return different values for different n
//it solves function=input 
{

    double *xt = (double*) malloc(sizeof(double)*Nvar);
    for (int i=0;i<Nvar; i++)
        xt[i]=x[i];
    
    
    int j;
    double dx,f,fmid,xmid,rtb;
    xt[ivar]=x1;
    f=(*func)(n,Nvar,xt,Npar,P)-input;
    xt[ivar]=x2;
    fmid=(*func)(n,Nvar,xt,Npar,P)-input;

    error(f*fmid >= 0.0,1,"rtbis","Root must be bracketed for bisection in rtbis f(x1)=%f   f(x2)=%f",f,fmid);
    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);// Orient the search so that f>0
    for (j=1;j<=MAXIT;j++) {// lies at x+dx.
        //printf(" x=%f \n",xt[ivar]   );
        xt[ivar]=rtb+(dx *= 0.5);
        fmid=(*func)(n,Nvar,xt,Npar,P)-input; //Bisection loop.
        //printf(" x+dx=%f fmind=%f   %d\n",xt[ivar], fmid ,MAXIT  );
        if (fmid <= 0.0) rtb=xt[ivar];
        if (fabs(dx) < xacc || fmid == 0.0){
            free(xt);
            return rtb;
        }
    }
    free(xt);
    error(1>0,1,"rtbis","Too many bisections in rtbis");
    return 0.0;
}

/*

double f1( int Nvar, double *x,int Npar,double  *P)
{
    return P[0]*log(P[0]*x[0])+P[0]*log(P[1]*x[1]);
}
double *f1k( int Nvar, double *x,int Npar,double  *P)
{
    double *r;
    r=(double*) malloc(sizeof(double)*Npar);
    r[0]=log(P[0]*x[0])+1+log(P[1]*x[1]);
    r[1]=(P[0]/P[1]);
    return r;
    
}

int main(){
    
    double *P;
    double **x,**y;
    int i;
    
    x=(double**) malloc(sizeof(double*)*4);
    y=(double**) malloc(sizeof(double*)*4);
    for (i=0;i<4;i++){
       x[i]=(double*) malloc(sizeof(double)*2);   
       y[i]=(double*) malloc(sizeof(double)*2);
    }
   x[0][0]=1;x[1][0]=4;x[2][0]=7;x[3][0]=10;
   x[0][1]=2;x[1][1]=3;x[2][1]=5;x[3][1]=9;
   
   y[0][0]=2;y[1][0]=5;y[2][0]=8;y[3][0]=11;
   y[0][1]=0.5;y[1][1]=2;y[2][1]=2;y[3][1]=2;
    int ensemble=4;
    int Nvar=2;
    int Npar=2;
        

//P=non_linear_fit( ensemble,x, y , Nvar,  Npar,  f1 , f1k);
P=non_linear_fit( ensemble,x, y , Nvar,  Npar,  f1 );
printf("P[0]=%f\n",P[0]);
printf("P[1]=%f\n",P[1]);

    return 0;
}
*/

/*Gnuplot result
gnuplot> f(x)=A*log(A*x)
gnuplot> fit f(x) 'tmp' i 0 u 1:2:3 yerr via A
iter      chisq       delta/lim  lambda   A            
   0 1.0935357306e+02   0.00e+00  1.35e+01    4.343868e+00
   1 9.0384567941e+00  -1.11e+06  1.35e+00    3.012440e+00
   2 2.7981042656e+00  -2.23e+05  1.35e-01    2.553096e+00
   3 2.7944989111e+00  -1.29e+02  1.35e-02    2.541241e+00
   4 2.7944983222e+00  -2.11e-02  1.35e-03    2.541377e+00
iter      chisq       delta/lim  lambda   A            

After 4 iterations the fit converged.
final sum of squares of residuals : 2.7945
rel. change during last iteration : -2.10747e-07

degrees of freedom    (FIT_NDF)                        : 3
rms of residuals      (FIT_STDFIT) = sqrt(WSSR/ndf)    : 0.965142
variance of residuals (reduced chisquare) = WSSR/ndf   : 0.931499
p-value of the Chisq distribution (FIT_P)              : 0.424406

Final set of parameters            Asymptotic Standard Error
=======================            ==========================
A               = 2.54138          +/- 0.1895       (7.455%)
*/
