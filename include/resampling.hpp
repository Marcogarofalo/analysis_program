#ifndef resampling_H
#define resampling_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <iostream>


class jackknife{
public:
    int Njack;
    double  *jack;
    std::string resampling="";
    
     //constructor
    jackknife(){}
        
    //constructor
    jackknife(const int Np1){
//         std::cout<< "constructing jack"<< std::endl;
        jack=(double*) malloc(sizeof(double) *Np1); 
        Njack=Np1;
    };
    jackknife(std::string s, const int Np1 ){
        jack=(double*) malloc(sizeof(double) *Np1); 
        Njack=Np1;
        resampling=s;
        if (s!="jack" && s!="boot"){ printf("jackknife construction:\n resampling must be jack or boot\n"); exit(-10);}
    };
    jackknife(std::string s, const int Np1, double *p){
        jack=(double*) malloc(sizeof(double) *Np1); 
        Njack=Np1;
        resampling=s;
        if (s!="jack" && s!="boot"){ printf("jackknife construction:\n resampling must be jack or boot\n"); exit(-10);}
        for(int j=0; j<Njack; j++)
            jack[j] = p[j];
    };
    
    
    //copy constructor
    jackknife (const jackknife &rhs)  
    {  
        std::cout<<"Copy constructor called "<<std::endl; 
        Njack=rhs.Njack;
        resampling=rhs.resampling;
        jack=(double*) malloc(sizeof(double)*Njack);
        for(int j=0; j<Njack; j++)
            jack[j] = rhs.jack[j];
        
    }  
    // copy assignment
    jackknife& operator = (const jackknife &rhs)    { 
        std::cout<<"Copy assignment called "<<std::endl; 
        if(Njack!=rhs.Njack) {printf("impossible to assign jack of different sizes"), exit(1);};
        for(int j=0; j<Njack; j++)
            jack[j] = rhs.jack[j];
        return *this; 
    } 
    jackknife& operator = ( double *rhs)    { 
        jack = rhs;
        return *this; 
    }
    jackknife& operator += (const jackknife &rhs)    { 

        Njack=rhs.Njack;
        for(int j=0; j<Njack; j++)
            jack[j] += rhs.jack[j];
        return *this; 
    }
    jackknife Add(const jackknife &rhs)     {
        jackknife tmp(rhs.Njack);
        for(int j=0; j<Njack; j++)
            tmp.jack[j] += rhs.jack[j];
        return tmp;
    }
    jackknife operator + (const jackknife &rhs)     { 
        jackknife tmp(rhs.Njack);
        for(int j=0; j<Njack; j++)
            tmp.jack[j] =(*this).jack[j]+ rhs.jack[j];
        
        return tmp; 
    }
    jackknife operator - (const jackknife &rhs)     { 
        jackknife tmp(rhs.Njack);
        for(int j=0; j<Njack; j++)
            tmp.jack[j] =(*this).jack[j]- rhs.jack[j];
        
        return tmp; 
    }
    jackknife operator * (const jackknife &rhs)     { 
        jackknife tmp(rhs.Njack);
        for(int j=0; j<Njack; j++)
            tmp.jack[j] =(*this).jack[j]* rhs.jack[j];
        
        return tmp; 
    }
    jackknife operator / (const jackknife &rhs) { 
        jackknife tmp(rhs.Njack);
        for(int j=0; j<Njack; j++)
            tmp.jack[j] =(*this).jack[j]/ rhs.jack[j];
        
        return tmp; 
    }
    
    //destructor
    ~jackknife(){
        std::cout<< "destructing jack of size="<< Njack << std::endl;
        free(jack);
    };
    // conversion to A (type-cast operator)
    operator double*() {return jack;};
    
//     //move contructor 
//     jackknife(jackknife&& rhs) noexcept{
//         Njack=rhs.Njack;
//         resampling=rhs.resampling;
//         jack = rhs.jack;
//         
//         rhs.Njack=0;
//         rhs.jack=nullptr;
//         rhs.resampling=nullptr;
//     }
//     // move assignment
//     jackknife& operator=(jackknife&& rhs) noexcept 
//     {
//         //if(this==&rhs); {printf("jackknife move assignment by itself");exit(-10);}
//         Njack=rhs.Njack;
//         resampling=rhs.resampling;
//         jack = rhs.jack;
//         
//         rhs.Njack=0;
//         rhs.jack=nullptr;
//         rhs.resampling=nullptr;
//         return *this;
//     }
//  
};




void free_jack(int N,int var , int t, double ****in);
void write_jack(int N, double *jack, char *outname);
void write_jack_bin(int N, double *jack, char *outname);

double ***read_jack(int N, int var, int T);
void write_jack_corr(int N, int t,double **jack, char *outname);
//create_jack
//in[#conf.][#variable][#time_cordinate][#re or im]
//returns the jacknife configuration from the data ****in
//the last entry of [#conf] is the average
double ****create_jack(int  N, int var, int t, double ****in);
void symmetrise_jackboot(int  Np1, int var, int T, double ****in, int sym=1 );
double *mean_jack(int N,int var,int t, int call, double ****jack, double function_jack(int,int,int,double ***) );
//mean_and_error_jack
//returns the mean and error from set of N  jacknife called *in  and the average stored in in[N]
double *mean_and_error_jack(int Np1, double *in);
double *mean_and_error_jack_biased(int Np1, double *in);
double *mean_and_error_jack_biased1(int Np1, double *in);

double *fake_jack(double mean,double error, int Njack,int seed);
double** covariance_jack(int Nobs, int Np1, double** in);




/////////////boot
double ****create_boot(int  N, int Nboot, int var, int t, double ****in);
double *mean_and_error_boot(int Np1, double *in);
double *fake_boot(double mean,double error, int Njack,int seed);
double** covariance_boot(int Nobs, int Np1, double** in);

/////////////////////
double *mean_and_error( const char *option , int Np1, double *in);  
double error_jackboot( const char *option , int Np1, double *in);  

const char  *smean_and_error( const char *option , int Np1, double *in);  
double ****create_resampling(const char *option, int  N, int var, int t, double ****in,int seed=123);
double *fake_sampling(const char *option,double mean,double error, int Njack,int seed);

double **covariance(const char *option , int Nobs, int Np1, double **in);
double **error_covariance(const char *option , int Nobs, int Np1, double **in);




double **fake_sampling_covariance(const char *option,double *mean, int Njack,int N, double **cov,int seed);
double**** create_boot_of_boot_from_jack(int  Njack, int Nboot, int en_tot, double*** in, int seed = 123);

///////////////////////////////////////////////////////
double* malloc_copy_jackboot(int Np1,  double *a);
    
void sum_jackboot(int Np1,  double *r, double *a, double *b);
void sub_jackboot(int Np1,  double *r, double *a, double *b);
void mult_jackboot(int Np1,  double *r, double *a, double *b);
void div_jackboot(int Np1,  double *r, double *a, double *b);
void invert_jackboot(int Np1,  double *r, double *a);

void scalar_times_jackboot(int Np1,  double *r, double *a, double s);


#endif
 
