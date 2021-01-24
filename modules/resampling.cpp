#define resampling_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "mutils.hpp"
#include "resampling.hpp"
#include "rand.hpp"
#include "global.hpp"
#include "tower.hpp"

void free_jack(int N,int var , int t, double ****in){
    int i,j,k,l;
    
    
    for (i=0;i<N;i++){
        for (j=0;j<var;j++){
            for(k=0;k<t;k++)
                free(in[i][j][k]);
            free(in[i][j]);    
        }
        free(in[i]);
    }
    free(in);
}

void write_jack(int N, double *jack,  char *outname){
    int j;
    FILE *f;
    
    f=fopen(outname,"ab");
    error(f==NULL,1,"write_jack ",
         "Unable to open output file %s",outname);
    
    for (j=0;j<N;j++)
        fprintf(f,"%.15g\n",jack[j]);
    
    fclose(f);
}



void write_jack_bin(int N, double *jack,  char *outname){
    int j;
    FILE *f;
    
    f=fopen(outname,"ab");
    error(f==NULL,1,"write_jack_bin ",
         "Unable to open output file %s",outname);
    
    fwrite(jack,sizeof(double),N,  f);
    
    fclose(f);
}



double ***read_jack(int N, int var, int T){
    int i,j,k;
    FILE **f;
    double ***jack;
    char **datafile;
    
    f=(FILE**) malloc(sizeof(FILE*)*var);
    datafile=(char**) malloc(sizeof(char*)*var);
    for (i=0;i<var;i++){
        datafile[i]=(char*) malloc(sizeof(char)*var);
        sprintf(datafile[i],"in%d",i);
        f[i]=fopen(datafile[i],"r");
        error(f==NULL,1,"read_jack ",
         "Unable to open input file %s",datafile[i]);
    }

    jack=(double***) malloc(sizeof(double**)*N);
    for (j=0;j<N;j++){
        jack[j]=(double**) malloc(sizeof(double*)*var);
        for (i=0;i<var;i++){
            jack[j][i]=(double*) malloc(sizeof(double)*T);
	    for (k=0;k<T;k++){
                 fscanf(f[i],"%lf\n",&jack[j][i][k]);
            } 
        }
    }
    
    
    for (i=0;i<var;i++){
    	fclose(f[i]);
        free(datafile[i]);
    }
    free(datafile);
    return jack;
}

void write_jack_corr(int N, int t,double **jack, char *outname){
    int i,j;
    FILE *f;
    
    f=fopen(outname,"w+");
    error(f==NULL,1,"write_jack ",
         "Unable to open output file %s",outname);
    
    for (j=0;j<N;j++)
        for (i=0;i<t;i++)
            fprintf(f,"%.15g\n",jack[i][j]);
    
    fclose(f);
}

//create_jack
//in[#conf.][#variable][#time_cordinate][#re or im]
//returns the jacknife configuration from the data ****in
//the last entry of [#conf] is the average
double ****create_jack(int  N, int var, int t, double ****in){
    double ****jack;
    int i,j,k,l;
    
    jack=(double****) malloc(sizeof(double***)*(N+1));
    for (i=0;i<=N;i++){
        jack[i]=(double***) malloc(sizeof(double**)*var);
        for (j=0;j<var;j++){
            jack[i][j]=(double**) malloc(sizeof(double*)*t);
            for(k=0;k<t;k++)
                jack[i][j][k]=(double*) calloc(2,sizeof(double));
        }
    }
    
    
    for(i=0;i<N;i++){
        for(j=0;j<var;j++){
            for(k=0;k<t;k++){
                for(l=0;l<2;l++){
                    jack[N][j][k][l]+=in[i][j][k][l];
                }
            }
        }
    }
    
    for(i=0;i<N;i++){
        for(j=0;j<var;j++){
            for(k=0;k<t;k++){
                for(l=0;l<2;l++){
                    jack[i][j][k][l]=(jack[N][j][k][l]-in[i][j][k][l])/( (double) (N-1) );
                }
            }
        }
    }
    
    for(j=0;j<var;j++){
            for(k=0;k<t;k++){
                for(l=0;l<2;l++){
                    jack[N][j][k][l]/=( (double) (N) );
                }
            }
        }
    

return  jack;
}



double *mean_jack(int N,int var,int t, int call, double ****jack, double function_jack(int,int,int,double ***) ){
    double *r;
    int i;
    
    r=(double*) malloc(sizeof(double)*N);
    for (i=0;i<N;i++)
        r[i]=function_jack(var,t,call,jack[i]);
    
    return r;
}

//mean_and_error_jack
//returns the mean and error from set of N  jacknife called *in  and the average stored in in[N]
double *mean_and_error_jack(int Np1, double *in){
    double *r;
    int i,N;
    
    N=Np1-1;
    r=(double*) calloc(2,sizeof(double));
    
    for(i=0;i<N;i++)
        r[0]+=in[i];
    
    r[0]/=((double) N);
    
    
    for(i=0;i<N;i++)
        r[1]+=(in[i]-r[0])*(in[i]-r[0]);
    
    r[1]*=(N-1.)/((double) N);
    r[1]=sqrt(r[1]);
    
    r[0]=((double)N)*in[N]-((double) (N-1))*r[0];
    
    return r;
}


double *mean_and_error_jack_biased1(int Np1, double *in){
    double *r;
    int i,N;
    
    N=Np1-1;
    r=(double*) calloc(2,sizeof(double));
    
    for(i=0;i<N;i++)
        r[0]+=in[i];
    
    r[0]/=((double) N);
    
    
    for(i=0;i<N;i++)
        r[1]+=(in[i]-r[0])*(in[i]-r[0]);
    
    r[1]*=(N-1.)/((double) N);
    r[1]=sqrt(r[1]);
    
    
    return r;
}

double *mean_and_error_jack_biased(int Np1, double *in){
    double *r;
    int i,N;
    
    N=Np1-1;
    r=(double*) calloc(2,sizeof(double));
    
    for(i=0;i<N;i++)
        r[0]+=in[i];
    
    r[0]/=((double) N);
    
    
    for(i=0;i<N;i++){
        r[1]+=(in[i]-r[0])*(in[i]-r[0]);
        //error(in[i]!=in[i],1,"mean_and_error_jack_biased","errore jack=%d is nan",i);
    }
    r[1]*=(N-1.)/((double) N);
    r[1]=sqrt(r[1]);
    
    r[0]=in[N];
    
    return r;
}

double *fake_jack(double mean,double error, int Njack,int seed){
    int i;
    double *r;
    
    r=(double*) malloc(sizeof(double)*Njack);
    srand(seed);   
    for (i=0;i<Njack-1;i++){
        r[i]=generateGaussianNoise(mean,  error/sqrt(Njack-2));
    }
    r[Njack-1]=mean;
    
    return r;
}
/////////////////////////////boot


double ****create_boot(int  N, int Nboot, int var, int t, double ****in){
    double ****boot;
    int i,j,k,l,b,ib;
    
    srand(123);
    
    
    boot=(double****) malloc(sizeof(double***)*(Nboot+1));
    for (i=0;i<=Nboot;i++){
        boot[i]=(double***) malloc(sizeof(double**)*var);
        for (j=0;j<var;j++){
            boot[i][j]=(double**) malloc(sizeof(double*)*t);
            for(k=0;k<t;k++)
                boot[i][j][k]=(double*) calloc(2,sizeof(double));
        }
    }
    
    
    for(ib=0;ib<Nboot;ib++){
       for(i=0;i<N;i++){
        b=(int )rand()%N;  //have to stay here, because I want that all the var, t , l to have the same random number 
        for(j=0;j<var;j++){
            for(k=0;k<t;k++){
                for(l=0;l<2;l++){
                    boot[ib][j][k][l]+=in[b][j][k][l];
                }
            }
        }
      }
    }
    
    for(i=0;i<N;i++){
        for(j=0;j<var;j++){
            for(k=0;k<t;k++){
                for(l=0;l<2;l++){
                    boot[Nboot][j][k][l]+=in[i][j][k][l];
                }
            }
        }
    }
    
    for(i=0;i<=Nboot;i++){
    for(j=0;j<var;j++){
            for(k=0;k<t;k++){
                for(l=0;l<2;l++){
                    boot[i][j][k][l]/=( (double) (N) );
                }
            }
        }
     }

return  boot;
}

double *fake_boot(double mean,double error, int Njack,int seed){
    int i;
    double *r;
    
    r=(double*) malloc(sizeof(double)*Njack);
    srand(seed);   
    for (i=0;i<Njack-1;i++){
        r[i]=generateGaussianNoise(mean,  error);
    }
    r[Njack-1]=mean;
    
    return r;
}


//mean_and_error_jack
//returns the mean and error from set of N  jacknife called *in  and the average stored in in[N]
double *mean_and_error_boot(int Np1, double *in){
    double *r;
    int i,N;
    
    N=Np1-1;
    r=(double*) calloc(2,sizeof(double));
    
    for(i=0;i<N;i++)
        r[0]+=in[i];
    
    r[0]/=((double) N);
    
    
    for(i=0;i<N;i++)
        r[1]+=(in[i]-r[0])*(in[i]-r[0]);
    
    
    r[1]/=((double) (N));
    r[1]=sqrt(r[1]);
    
    r[0]=2*in[N]-r[0];
    
    return r;
}
////////////////////////////////

double *mean_and_error(const char *option , int Np1, double *in){
    double *r;
    if( strcmp(option,"jack")==0)
              r=mean_and_error_jack_biased(Np1, in);
    else if( strcmp(option,"boot")==0)
              r=mean_and_error_boot(Np1, in);
    else 
        error(0==0,1,"mea_and_error call","mea_and_error called with %s while the only options supported are jack or boot",option);
    return r;
}


// give the error the jackknife series in[Np1=Njack]
double error_jackboot(const char *option , int Np1, double *in){
    double *r,r1;
    if( strcmp(option,"jack")==0)
              r=mean_and_error_jack_biased(Np1, in);
    else if( strcmp(option,"boot")==0)
              r=mean_and_error_boot(Np1, in);
    else 
        error(0==0,1,"mea_and_error call","mea_and_error called with %s while the only options supported are jack or boot",option);
    r1=r[1];
    free(r);
    return r1;
}

const char *smean_and_error(const char *option , int Np1, double *in){
    double *r;
    char *s;
    if( strcmp(option,"jack")==0)
              r=mean_and_error_jack_biased(Np1, in);
    else if( strcmp(option,"boot")==0)
              r=mean_and_error_boot(Np1, in);
    else 
        error(0==0,1,"mea_and_error call","mea_and_error called with %s while the only options supported are jack or boot",option);
    
    //error(r[0]!=r[0],0,"smean_and_error","error  mean value is  nan");
    //error(r[1]!=r[1],0,"smean_and_error","error  error value is  nan");
    //error(r[1]<0,0,    "smean_and_error","error  is nevative");
    if (r[0]!=r[0]) {s=(char*) malloc(sizeof(char)*(4)); mysprintf(s,4,"nan");return s;}
    if (r[1]!=r[1]) {s=(char*) malloc(sizeof(char)*(4)); mysprintf(s,4,"nan"); return s;}
    if (r[1]<=0)    {s=(char*) malloc(sizeof(char)*(4)); mysprintf(s,4,"nan"); return s;}
    int m10=(int) log10(fabs(r[0]));
    if (m10<-19) m10=-19;
    int e10=(int) log10(r[1]);
    int e10_10;
    if (e10<=1e-6) e10_10=1;
    else e10_10=abs((int) log10( abs(e10)))+1;
    int d=abs(m10-e10)+1;
    int s_size=(d+2+2+3+2+e10_10)+2+130;
    if (r[0]<0) s_size++; 
    double wm=( r[0]/pow(10,e10-1));
    double we=(r[1]/pow(10,e10-1));
    s=(char*) malloc(sizeof(char)*(s_size));
   /* printf("%g  %g\n",r[0],r[1]);
    printf("%d  %d\n",m10,e10);
    printf("%.1f(%.1f)e%+-d     %d\n",wm,we,e10-1,s_size);*/
    mysprintf(s,s_size,"%.1f(%.1f)10^{%d}",wm,we,e10-1);
    free(r);
    return s;
}


double ****create_resampling(const char *option, int  N, int var, int t, double ****in){
    double ****r;
    int Nboot=Nbootstrap;
    if( strcmp(option,"jack")==0)
              r=create_jack(  N,  var,  t,  in);
    else if( strcmp(option,"boot")==0)
              r=create_boot(  N,  Nboot,  var,  t, in);
    else 
        error(0==0,1,"create_resampling call","create_resampling called with %s while the only options supported are jack or boot",option);    
    return r;

}

double *fake_sampling(const char *option,double mean,double delta_mean, int Njack,int seed){
    double *r;
    int Nboot=Nbootstrap;
    if( strcmp(option,"jack")==0)
              r=fake_jack(  mean, delta_mean,  Njack,seed);
    else if( strcmp(option,"boot")==0)
              r=fake_boot(  mean, delta_mean,  Njack,seed);
    else    
        error(0==0,1,"fake_sampling call","fake_sampling called with %s while the only options supported are jack or boot",option);    

    return r;
}


double **fake_jack_covariance(double *mean, int Njack,int seed,int N, double **cov){
    int i,k;
    double **r,**r1;
    
    double **covj=(double**) malloc(sizeof(double*)*N);
    for (i=0;i<N;i++){
        covj[i]=(double*) malloc(sizeof(double)*N);
        for (k=0;k<N;k++){
            covj[i][k]=cov[i][k]/(Njack-2);
        }
    }
    
    r=(double**) malloc(sizeof(double*)*Njack);
    srand(seed);   
    for (i=0;i<Njack-1;i++){
        r[i]=generate_correlatedNoise(N,mean,  covj);//generateGaussianNoise(mean,  error/sqrt(Njack-2));
    }
    r[Njack-1]=(double*) malloc(sizeof(double)*N);
    for (i=0;i<N;i++){
        r[Njack-1][i]=mean[i];
        free(covj[i]);
    }
    free(covj);
    
    r1=swap_indices(Njack,N,r);
    free_2(Njack,r);
    return r1;
}

double **fake_boot_covariance(double *mean, int Njack,int seed,int N, double **cov){
    int i;
    double **r,**r1;
    
    r=(double**) malloc(sizeof(double*)*Njack);
    srand(seed);   
    for (i=0;i<Njack-1;i++){
        r[i]=generate_correlatedNoise(N,mean,  cov);
    }
    r[Njack-1]=(double*) malloc(sizeof(double)*N);
    for (i=0;i<N;i++)
        r[Njack-1][i]=mean[i];
    
    r1=swap_indices(Njack,N,r);
    free_2(Njack,r);
    return r1;
}
/*
double **reduce_covariance_matrix1(int N, double **cov, int *N1){
    int Nuncor=0;
    int i,j,ii,jj;
    double **r;
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            if (  fabs(cov[i][j]/sqrt(cov[i][i]*cov[j][j])-1)<0.001  ){
                Nuncor++;
                r=double_malloc_2(N-1,N-1);
                for(ii=0;ii<j;ii++){//remove j row and column
                    for(jj=0;jj<j;jj++)
                        r[ii][jj]=cov[ii][jj];
                    for(jj=j+1;jj<N;jj++)
                        r[ii][jj]=cov[ii][jj+1];
                }
                for(ii=j+1;ii<N;ii++){
                    for(jj=0;jj<j;jj++)
                        r[ii][jj]=cov[ii+1][jj];
                    for(jj=j+1;jj<N;jj++)
                        r[ii][jj]=cov[ii+1][jj+1];
                }// end remove j row and column
                *N1=N-1;
                return r;
            }
        }
    }
    if (Nuncor==0){
        r=double_malloc_2(N,N);
        for(i=0;i<N;i++)
            for(j=i+1;j<N;j++)
                r[i][j]=cov[i][j];
    }
    N1=N;
    return r;
}

double **reduce_covariance_matrix(int N, double **cov){
    int N1=N;
    int N2=N,i,j;
    int stop=1;
    
    double **tmp,**tmp1;
    tmp=double_malloc_2(N,N);
    for(i=0;i<N;i++)
            for(j=i+1;j<N;j++)
                tmp[i][j]=cov[i][j];
            
    while(stop>0){
        tmp1=reduce_covariance_matrix1(N2,tmp,N1);
        free_2(tmp);
        tmp=tmp1;
        stop=N2-N1;
        
    }
    
    return tmp;
    
}
*/

//given a covariance matrix  cov[N][N]
//returns a double array of r[N][jacknife] of jacknife or bootstrap
double **fake_sampling_covariance(const char *option,double *mean, int Njack,int N, double **cov,int seed){
    double **r;
    int Nboot=Nbootstrap;
    /*
    int i,j;
    int Nuncor=0    ;
    int *sametable;
    
    if (N==1) error(0==0,1,"you can not call fake_sampling_covariance for one signle variable","");
    
    yn=0;
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            if (  cov[i][j]fabs(cov[i][j]/sqrt(cov[i][i]*cov[j][j])-1)<0.001  )
                Nuncor++;
        }
    }
    Nuncor=(N*(N-1)/2) ;
    */
    if( strcmp(option,"jack")==0)
              r=fake_jack_covariance(  mean,  Njack,seed,N,cov);
    else if( strcmp(option,"boot")==0)
              r=fake_boot_covariance(  mean,   Njack,seed,N,cov);
    else    
        error(0==0,1,"fake_sampling call","fake_sampling called with %s while the only options supported are jack or boot",option);    

    return r;
}


double **covariance_jack(int Nobs, int Np1, double **in){
    double **r,*ave;
    int i,k,l,N;
    
    N=Np1-1;
    ave=(double*) calloc(Nobs,sizeof(double));
    r=(double**) malloc(Nobs*sizeof(double*));
    
    for(k=0;k<Nobs;k++){
        for(i=0;i<N;i++)
            ave[k]+=in[k][i];
        ave[k]/=((double) N);
        r[k]=(double*) calloc(Nobs,sizeof(double));
    }
    for(k=0;k<Nobs;k++){
        for(l=k;l<Nobs;l++){    
            for(i=0;i<N;i++){
                r[k][l]+=(in[k][i]-ave[k])*(in[l][i]-ave[l]);
                //error(in[i]!=in[i],1,"mean_and_error_jack_biased","errore jack=%d is nan",i);
            }
        r[k][l]*=(N-1.)/((double) N);
        //r[k][l]=sqrt(r[k][l]);
        }
        for(l=0;l<k;l++){
            //r[l][k]/=sqrt(r[k][k]*r[l][l]);
            r[k][l]=r[l][k];
        }
    }
    free(ave);
    return r;
    
}

double **covariance_boot(int Nobs, int Np1, double **in){
    double **r,*ave;
    int i,k,l,N;
    
    N=Np1-1;
    ave=(double*) calloc(Nobs,sizeof(double));
    r=(double**) malloc(Nobs*sizeof(double*));
    
    for(k=0;k<Nobs;k++){
         for(i=0;i<N;i++)
                ave[k]+=in[k][i];
            ave[k]/=((double) N);
            r[k]=(double*) calloc(Nobs,sizeof(double));
    }
    for(k=0;k<Nobs;k++){
        for(l=k;l<Nobs;l++){    
            for(i=0;i<N;i++){
                r[k][l]+=(in[k][i]-ave[k])*(in[l][i]-ave[l]);
                //error(in[i]!=in[i],1,"mean_and_error_jack_biased","errore jack=%d is nan",i);
            }
        r[k][l]/=((double) N);
        //r[k][l]=sqrt(r[k][l]);
        }
        for(l=0;l<k;l++){
            //r[l][k]/=sqrt(r[k][k]*r[l][l]);
            r[k][l]=r[l][k];
        }
        r[k][k]=1;
    }
    free(ave);
    return r;
    
}

//conpute the covariance of in[variables=0,...,Nobs-1][jacknifes]
//return Nobs x Nobs matrix
double **covariance(const char *option , int Nobs, int Np1, double **in){
    double **r;
    
    if( strcmp(option,"jack")==0)
              r=covariance_jack(  Nobs,  Np1, in);
    else if( strcmp(option,"boot")==0)
              r=covariance_boot(   Nobs,  Np1, in);
    else    
        error(0==0,1,"covariance call","covariance called with %s while the only options supported are jack or boot",option);    

    return r;
}
