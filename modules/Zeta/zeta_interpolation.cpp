#define zeta_interpolation_C

#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>

#include <global.hpp>
#include <linear_fit.hpp>
#include <tower.hpp>
#include "header_phi4.hpp"
#include "zeta_interpolation.hpp"
#include "resampling.hpp"
#include "mutils.hpp"
//#include "header_phi4.hpp"

extern "C" { 
    #include "../external/rzeta/src/dzeta_function.h"
}

/*
void zeta_interpolation::Init(char *resampling, std::vector<int>  myen,  std::vector<cluster::IO_params> paramsj, std::vector<data_phi> gjack ){
        double a=timestamp();
        double z[2];
        int tmp_vec[3]={0,0,0};
        // init the workspace of the dzeta_function
        dzeta_function(z,  1 ,0 , 0, tmp_vec, 1, 1. , 1.e-3, 1.e6 ,1);
        etot=myen.size();
        std::vector< std::vector<int> > momenta(5,std::vector<int>(3));
        momenta[0][0]=0; momenta[0][1]=0; momenta[0][2]=0; 
        momenta[1][0]=1; momenta[1][1]=0; momenta[1][2]=0; 
        momenta[2][0]=1; momenta[2][1]=1; momenta[2][2]=0; 
        momenta[3][0]=1; momenta[3][1]=1; momenta[3][2]=1; 
        momenta[4][0]=0; momenta[4][1]=0; momenta[4][2]=0; 
        
        mom=momenta;
        krange=double_malloc_4(etot,ntot,Nm,2);
        m=double_malloc_2(etot,Nm);
        printf("initializing zeta grid for interpolation \n");
        int Njack=gjack[0].Njack;
        grid=(double****) malloc(sizeof(double***)*myen.size());
        kint=(double****) malloc(sizeof(double***)*myen.size());
        for (int e=0; e<etot;e++){
            L.emplace_back(paramsj[myen[e]].data.L[1] );
            double err=error_jackboot(resampling,Njack,gjack[myen[e]].jack[1]);
            m[e][0]=gjack[myen[e]].jack[1][Njack-1]-err; 
            for (int im=0; im<Nm;im++){
                m[e][im]=m[e][0] +2.*err*im/((double) Nm) ;
            }
            
            double twopiL2=(2.*pi_greco/L[e]);
            twopiL2*=twopiL2;
            grid[e]=(double***) malloc(sizeof(double**)*momenta.size());
            kint[e]=(double***) malloc(sizeof(double**)*momenta.size());
            for (int d=0; d<ntot;d++){
                grid[e][d]=(double**) malloc(sizeof(double*)*Nm);
                kint[e][d]=(double**) malloc(sizeof(double*)*Nm);
                for (int im=0; im<Nm;im++){
                    //int dvec[3]={momenta[d][0],momenta[d][1],momenta[d][2]};
                    int dvec[3],dvec1[3],dvec2[3];
                    dvec[0]=mom[d][0]; dvec[1]=mom[d][1]; dvec[2]=mom[d][2];
                    
                    if(d==0){//E2_0
                        dvec1[0]=0; dvec1[1]=0; dvec1[2]=0;
                        dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
                    }
                    else if(d==1){//E2_0_p1
                        dvec1[0]=1; dvec1[1]=0; dvec1[2]=0;
                        dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
                    }
                    else if(d==2){//E2_0_p11
                        dvec1[0]=1; dvec1[1]=1; dvec1[2]=0;
                        dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
                    }
                    else if(d==3){//E2_0_p111
                        dvec1[0]=1; dvec1[1]=1; dvec1[2]=1;
                        dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
                    }
                    else if(d==4){//E2_0_A1
                        dvec1[0]=1; dvec1[1]=0; dvec1[2]=0;
                        dvec2[0]=-1; dvec2[1]=0; dvec2[2]=0;
                    }
                    else {
                        exit(1);
                    }
                    
                    //init the range of k
                    double E1f=sqrt(m[e][im]*m[e][im]+twopiL2*(dvec1[0]*dvec1[0]+dvec1[1]*dvec1[1]+dvec1[2]*dvec1[2])   );
                    double E2f=sqrt(m[e][im]*m[e][im]+twopiL2*(dvec2[0]*dvec2[0]+dvec2[1]*dvec2[1]+dvec2[2]*dvec2[2])   );
                    double Ef=E1f+E2f;
                    double ECMfsq=Ef*Ef-twopiL2*(dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
                    double kf=sqrt(ECMfsq/4. -m[e][im]*m[e][im]);
                    krange[e][d][im][0]=kf+1e-10;
                    E1f=sqrt(m[e][im]*m[e][im]+twopiL2*((dvec1[0])*(dvec1[0])+dvec1[1]*dvec1[1]+(dvec1[2]+1)*(dvec1[2]+1))   );
                    E2f=sqrt(m[e][im]*m[e][im]+twopiL2*((dvec2[0])*(dvec2[0])+dvec2[1]*dvec2[1]+(dvec2[2]-1)*(dvec2[2]-1))   );
                    Ef=E1f+E2f;
                    ECMfsq=Ef*Ef-twopiL2*((dvec[0])*(dvec[0])+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
                    double kf1=sqrt(ECMfsq/4.-m[e][im]*m[e][im]);
                    krange[e][d][im][1]=kf1-1e-10;
                    int iter=(krange[e][d][im][1]-krange[e][d][im][0]) /h+1 ;
                    error(iter<3,1,"zeta_interpolation::Init","points to evaluate the zeta in k are %d to little to interpolate try decreasing h",iter);
//                     printf("e=%d  mom=%d=(%d,%d,%d)  im=%d  kiter=%d\n",e,d,dvec[0],dvec[1],dvec[2],im,iter );
                    
                    grid[e][d][im]=(double*) malloc(sizeof(double)*iter);
                    kint[e][d][im]=(double*) malloc(sizeof(double)*iter);
                    #pragma omp parallel for  private(z)  shared(h,iter,e,d,im,twopiL2,dvec,m,krange) 
                    for (int ik=0;ik< iter;ik++){
                        double k=krange[e][d][im][0]+h*ik;
                        if(ik==iter-1) k=krange[e][d][im][1];
                        double ECM2=4*(k*k+m[e][im]*m[e][im]);
                        double E2=ECM2+twopiL2*(dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
                        double gamma=sqrt(E2/ECM2);
    //                     printf("k=%g   gamma=%g  dvec=(%d,%d,%d) mass=%g  e=%d\n",k,gamma,dvec[0],dvec[1],dvec[2], mass[ens[e]],ens[e]);
                        dzeta_function(z,  k*k/(twopiL2) ,0 , 0, dvec, gamma, 1. , 1.e-3, 1.e6 ,4);
                        grid[e][d][im][ik]=z[0]; //real part
                        kint[e][d][im][ik]=k; //real part
                    }
                    for (int ik=iter-10;ik< iter;ik++){
                        error(grid[e][d][im][0]*grid[e][d][im][ik]>0,1,"zeta_interpolation::Init","over the pole ");
                    }
                }
            }
        }
        double b=timestamp();
        printf("time to Init zeta_interpolation=%g s\n",b-a);
}

*/
void zeta_interpolation::Init(char *resampling, std::vector<int>  myen,  std::vector<cluster::IO_params> paramsj, std::vector<data_phi> gjack ){
    double a=timestamp();
    double z[2];
    int tmp_vec[3]={0,0,0};
    // init the workspace of the dzeta_function
    dzeta_function(z,  1 ,0 , 0, tmp_vec, 1, 1. , 1.e-3, 1.e6 ,1);
    etot=myen.size();
    std::vector< std::vector<int> > momenta(ntot,std::vector<int>(3));
    momenta[0][0]=0; momenta[0][1]=0; momenta[0][2]=0; 
    momenta[1][0]=1; momenta[1][1]=0; momenta[1][2]=0; 
    momenta[2][0]=1; momenta[2][1]=1; momenta[2][2]=0; 
    momenta[4][0]=1; momenta[4][1]=1; momenta[4][2]=1; 
    momenta[3][0]=0; momenta[3][1]=0; momenta[3][2]=0; 
    
    mom=momenta;
    krange=double_malloc_4(etot,ntot,Nm,2);
    m=double_malloc_2(etot,Nm);
    printf("initializing zeta grid for interpolation \n");
    int Njack=gjack[0].Njack;
    grid=(double****) malloc(sizeof(double***)*myen.size());
    kint=(double****) malloc(sizeof(double***)*myen.size());
    for (int e=0; e<etot;e++){
        L.emplace_back(paramsj[myen[e]].data.L[1] );
        double err=error_jackboot(resampling,Njack,gjack[myen[e]].jack[1]);
        double twopiL=(2.*pi_greco/L[e]);
        double twopiL2=twopiL*twopiL;
        m[e][0]=(gjack[myen[e]].jack[1][Njack-1]-3.*err)/twopiL; 
        for (int im=1; im<Nm;im++){
            m[e][im]=m[e][0] +(2.*(3.*err)*im/((double) Nm))/twopiL ;
        }
        
        
        grid[e]=(double***) malloc(sizeof(double**)*momenta.size());
        kint[e]=(double***) malloc(sizeof(double**)*momenta.size());
        for (int d=0; d<ntot;d++){
            grid[e][d]=(double**) malloc(sizeof(double*)*Nm);
            kint[e][d]=(double**) malloc(sizeof(double*)*Nm);
            for (int im=0; im<Nm;im++){
                //int dvec[3]={momenta[d][0],momenta[d][1],momenta[d][2]};
                int dvec[3],dvec1[3],dvec2[3];
                dvec[0]=mom[d][0]; dvec[1]=mom[d][1]; dvec[2]=mom[d][2];
                
                if(d==0){//E2_0
                    dvec1[0]=0; dvec1[1]=0; dvec1[2]=0;
                    dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
                }
                else if(d==1){//E2_0_p1
                    dvec1[0]=1; dvec1[1]=0; dvec1[2]=0;
                    dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
                }
                else if(d==2){//E2_0_p11
                    dvec1[0]=1; dvec1[1]=1; dvec1[2]=0;
                    dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
                }
                else if(d==4){//E2_0_p111
                    dvec1[0]=1; dvec1[1]=1; dvec1[2]=1;
                    dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
                }
                else if(d==3){//E2_0_A1
                    dvec1[0]=1; dvec1[1]=0; dvec1[2]=0;
                    dvec2[0]=-1; dvec2[1]=0; dvec2[2]=0;
                }
                else {
                    exit(1);
                }
                
                //init the range of k
                double E1f=sqrt(m[e][im]*m[e][im]+(dvec1[0]*dvec1[0]+dvec1[1]*dvec1[1]+dvec1[2]*dvec1[2])   );
                double E2f=sqrt(m[e][im]*m[e][im]+(dvec2[0]*dvec2[0]+dvec2[1]*dvec2[1]+dvec2[2]*dvec2[2])   );
                double Ef=E1f+E2f;
                double ECMfsq=Ef*Ef-(dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
                double kf=sqrt(ECMfsq/4. -m[e][im]*m[e][im]);
                krange[e][d][im][0]=kf*kf;//+1e-10;
                E1f=sqrt(m[e][im]*m[e][im]+((dvec1[0])*(dvec1[0])+dvec1[1]*dvec1[1]+(dvec1[2]+1)*(dvec1[2]+1))   );
                E2f=sqrt(m[e][im]*m[e][im]+((dvec2[0])*(dvec2[0])+dvec2[1]*dvec2[1]+(dvec2[2]-1)*(dvec2[2]-1))   );
                Ef=E1f+E2f;
                ECMfsq=Ef*Ef-((dvec[0])*(dvec[0])+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
                double kf1=sqrt(ECMfsq/4.-m[e][im]*m[e][im]);
                krange[e][d][im][1]=kf1*kf1;//-1e-10;
//                 int iter=(krange[e][d][im][1]-krange[e][d][im][0]) /h+1 ;
                int iter=Nh;
                double h=(krange[e][d][im][1]-krange[e][d][im][0])/(Nh+2.);
                krange[e][d][im][0]+=h;
                krange[e][d][im][1]-=h;
//                 double h=(krange[e][d][im][1]-krange[e][d][im][0])/(Nh+2);
                error(iter<3,1,"zeta_interpolation::Init","points to evaluate the zeta in k are %d to little to interpolate try decreasing h",iter);
                //                     printf("e=%d  mom=%d=(%d,%d,%d)  im=%d  kiter=%d\n",e,d,dvec[0],dvec[1],dvec[2],im,iter );
                
                grid[e][d][im]=(double*) malloc(sizeof(double)*(iter));
                kint[e][d][im]=(double*) malloc(sizeof(double)*(iter));
                #pragma omp parallel for  private(z)  shared(h,iter,e,d,im,twopiL2,dvec,m,krange) 
                for (int ik=0;ik< iter;ik++){
                    double k=krange[e][d][im][0]+h*ik;
                    if(ik==iter-1) k=krange[e][d][im][1];
                    double ECM2=4*(k+m[e][im]*m[e][im]);
                    double E2=ECM2+(dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
                    double gamma=sqrt(E2/ECM2);
                    //                     printf("k=%g   gamma=%g  dvec=(%d,%d,%d) mass=%g  e=%d\n",k,gamma,dvec[0],dvec[1],dvec[2], mass[ens[e]],ens[e]);
                    dzeta_function(z,  k ,0 , 0, dvec, gamma, 1. , 1.e-3, 1.e6 ,3);
                    grid[e][d][im][ik]=z[0]; //real part
                    kint[e][d][im][ik]=k; //real part
//                     printf("e=%d  n=%d  q2=%g h=%g  z=%g\n",e,d,k,h,z[0]);
                }
                for (int ik=iter-2;ik< iter;ik++){
                    error(grid[e][d][im][0]*grid[e][d][im][ik]>0,1,"zeta_interpolation::Init","over the pole e=%d n=%d  mass(im=%d)=%g  q2=[%g,%g]  h=%g",e,d,im,m[e][im],krange[e][d][im][0],krange[e][d][im][1],h );
                }
            }
        }
    }
    double b=timestamp();
    allocated=0;
    printf("time to Init zeta_interpolation=%g s\n",b-a);
}

void zeta_interpolation::Init_Lmq( std::vector<int>  Ls,  std::vector<double> masses, std::vector<double> err_mass ){
    double a=timestamp();
    error(Ls.size()!=masses.size(),1," Init_Lmq", "Ls vector length %d not the same ad masses %d",Ls.size(),masses.size());
    error(err_mass.size()!=masses.size(),1," Init_Lmq", "error vector length  %d not the same ad masses %d",err_mass.size(),masses.size());
    double z[2];
    int tmp_vec[3]={0,0,0};
    // init the workspace of the dzeta_function
    dzeta_function(z,  1 ,0 , 0, tmp_vec, 1, 1. , 1.e-3, 1.e6 ,1);
    etot=Ls.size();
    std::vector< std::vector<int> > momenta(ntot,std::vector<int>(3));
    momenta[0][0]=0; momenta[0][1]=0; momenta[0][2]=0; 
    momenta[1][0]=1; momenta[1][1]=0; momenta[1][2]=0; 
    momenta[2][0]=1; momenta[2][1]=1; momenta[2][2]=0; 
    momenta[4][0]=1; momenta[4][1]=1; momenta[4][2]=1; 
    momenta[3][0]=0; momenta[3][1]=0; momenta[3][2]=0; 
    
    mom=momenta;
    krange=double_malloc_4(etot,ntot,Nm,2);
    m=double_malloc_2(etot,Nm);
    printf("initializing zeta grid for interpolation \n");
    
    grid=(double****) malloc(sizeof(double***)*etot);
    kint=(double****) malloc(sizeof(double***)*etot);
    for (int e=0; e<etot;e++){
        L.emplace_back(Ls[e] );
        printf(" init L=%d\n",Ls[e]);
        double twopiL=(2.*pi_greco/L[e]);
        double twopiL2=twopiL*twopiL;
        m[e][0]=(masses[e]-5.*err_mass[e]); 
        for (int im=1; im<Nm;im++){
            m[e][im]=m[e][0]+(2.*(5.*err_mass[e])*im/((double) Nm)) ;
        }
        
        grid[e]=(double***) malloc(sizeof(double**)*momenta.size());
        kint[e]=(double***) malloc(sizeof(double**)*momenta.size());
        for (int d=0; d<ntot;d++){
            grid[e][d]=(double**) malloc(sizeof(double*)*Nm);
            kint[e][d]=(double**) malloc(sizeof(double*)*Nm);
            for (int im=0; im<Nm;im++){
                //int dvec[3]={momenta[d][0],momenta[d][1],momenta[d][2]};
                int dvec[3],dvec1[3],dvec2[3],dmax1[3],dmax2[3];
                dvec[0]=mom[d][0]; dvec[1]=mom[d][1]; dvec[2]=mom[d][2]; 
                
                if(d==0){//E2_0
                    dvec1[0]=0; dvec1[1]=0; dvec1[2]=0;
                    dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
                    dmax1[0]=1; dmax1[1]=0; dmax1[2]=0;
                    dmax2[0]=-1; dmax2[1]=0; dmax2[2]=0;
                }
                else if(d==1){//E2_0_p1
                    dvec1[0]=1; dvec1[1]=0; dvec1[2]=0;
                    dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
                    dmax1[0]=1; dmax1[1]=1; dmax1[2]=0;
                    dmax2[0]=0; dmax2[1]=-1; dmax2[2]=0;
                }
                else if(d==2){//E2_0_p11
                    dvec1[0]=1; dvec1[1]=1; dvec1[2]=0;
                    dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
                    dmax1[0]=1; dmax1[1]=0; dmax1[2]=0;
                    dmax2[0]=0; dmax2[1]=1; dmax2[2]=0;
                }
                else if(d==4){//E2_0_p111
                    dvec1[0]=1; dvec1[1]=1; dvec1[2]=1;
                    dvec2[0]=0; dvec2[1]=0; dvec2[2]=0;
                    dmax1[0]=1; dmax1[1]=1; dmax1[2]=0;
                    dmax2[0]=0; dmax2[1]=0; dmax2[2]=1;
                }
                else if(d==3){//E2_0_A1
                    dvec1[0]=1; dvec1[1]=0; dvec1[2]=0;
                    dvec2[0]=-1; dvec2[1]=0; dvec2[2]=0;
                    dmax1[0]=1; dmax1[1]=0; dmax1[2]=1;
                    dmax2[0]=-1; dmax2[1]=0; dmax2[2]=-1;
                }
                else {
                    exit(1);
                }
                
                //init the range of k
                double E1f=sqrt(m[e][im]*m[e][im]/twopiL2+(dvec1[0]*dvec1[0]+dvec1[1]*dvec1[1]+dvec1[2]*dvec1[2])   );
                double E2f=sqrt(m[e][im]*m[e][im]/twopiL2+(dvec2[0]*dvec2[0]+dvec2[1]*dvec2[1]+dvec2[2]*dvec2[2])   );
                double Ef=E1f+E2f;
                double ECMfsq=Ef*Ef-(dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
                double kf=(ECMfsq/4. -m[e][im]*m[e][im]/twopiL2);
                krange[e][d][im][0]=kf;//+1e-10;
                E1f=sqrt(m[e][im]*m[e][im]/twopiL2+((dmax1[0])*(dmax1[0])+dmax1[1]*dmax1[1]+(dmax1[2])*(dmax1[2]))   );
                E2f=sqrt(m[e][im]*m[e][im]/twopiL2+((dmax2[0])*(dmax2[0])+dmax2[1]*dmax2[1]+(dmax2[2])*(dmax2[2]))   );
                Ef=E1f+E2f;
                ECMfsq=Ef*Ef-((dvec[0])*(dvec[0])+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
                double kf1=(ECMfsq/4.-m[e][im]*m[e][im]/twopiL2);
                krange[e][d][im][1]=kf1;//-1e-10;
                //                 int iter=(krange[e][d][im][1]-krange[e][d][im][0]) /h+1 ;
                int iter=Nh;
                double h=(krange[e][d][im][1]-krange[e][d][im][0])/(Nh+2.);
                krange[e][d][im][0]+=h;
                krange[e][d][im][1]-=h;
                //                 double h=(krange[e][d][im][1]-krange[e][d][im][0])/(Nh+2);
                error(iter<3,1,"zeta_interpolation::Init","points to evaluate the zeta in k are %d to little to interpolate try decreasing h",iter);
                //                     printf("e=%d  mom=%d=(%d,%d,%d)  im=%d  kiter=%d\n",e,d,dvec[0],dvec[1],dvec[2],im,iter );
                
                grid[e][d][im]=(double*) malloc(sizeof(double)*(iter));
                kint[e][d][im]=(double*) malloc(sizeof(double)*(iter));
                 #pragma omp parallel for  private(z)  shared(h,iter,e,d,im,twopiL2,dvec,m,krange) 
                for (int ik=0;ik< iter;ik++){
                    double k=krange[e][d][im][0]+h*ik;
                    if(ik==iter-1) k=krange[e][d][im][1];
                    double ECM2=4*(k+m[e][im]*m[e][im]/twopiL2);
                    double E2=ECM2+(dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
                    double gamma=sqrt(E2/ECM2); 
//                     printf("k=%g   gamma=%g  dvec=(%d,%d,%d) mass=%g  e=%d  kmin=%g kmax=%g h=%g, n=%d kf=%g kf1=%g\n",k,gamma,dvec[0],dvec[1],dvec[2],m[e][im],e, krange[e][d][im][0],krange[e][d][im][1],h,d, kf,kf1);
                    dzeta_function(z,  k ,0 , 0, dvec, gamma, 1. , 1.e-3, 1.e6 ,4);
                    grid[e][d][im][ik]=z[0]; //real part
                    kint[e][d][im][ik]=k; //real part
                    //                     printf("e=%d  n=%d  q2=%g h=%g  z=%g\n",e,d,k,h,z[0]);
                }
                for (int ik=iter-2;ik< iter;ik++){
                    error(grid[e][d][im][0]*grid[e][d][im][ik]>0,1,"zeta_interpolation::Init","over the pole e=%d n=%d  mass(im=%d)=%g  q2=[%g,%g]  h=%g",e,d,im,m[e][im],krange[e][d][im][0],krange[e][d][im][1],h );
                }
            }
        }
    }
    double b=timestamp();
    allocated=0;
    printf("time to Init zeta_interpolation=%g s\n",b-a);
}

double zeta_interpolation::compute(double inL, int n, double mass,  double k ){
    
    int e=-1;
    double z; // to be returned
    
    for (int e1=0; e1<L.size(); e1++){
        if (fabs(inL-L[e1])<0.1)
            e=e1;
    }
    //error(e==-1,3, "zeta_interpolation::compute","no ensemble with L=%g ",inL);
//     printf("e=%d    L=%g\n",e,inL);
    
    double r[3];
    double **M, *y;
    M=double_malloc_2(3,3);
    y=(double*) malloc(sizeof(double)*3);
    
    if (e!=-1){
        //select masses 
        int  m_min=0;
        for(int im=1; im<(Nm-2);im++) {
            if (mass> m[e][im] )
                m_min++;
        }
        for(int  im=m_min;im<m_min+3;im++){
            double kmin=krange[e][n][im][0];
            double kmax=krange[e][n][im][1];
//             printf("qsq: e=%d  n=%d  im=%d min=%.12g max=%.12g\n",e,n,im,kmin,kmax);
            //error(k<kmin || k>kmax  , 3,"zeta_interpolation::compute"," k=%g out of range [%g,%g ] ",k,kmin,kmax);
            
            // select k for interpolation
            int Niter=Nh;//(kmax-kmin) /h+1 ;
            double h=(kmax-kmin)/Nh;
            int x0=(k-kmin)/h ;
            if (x0<0){// if on the edege use the first 3 point
                x0=0; 
    //             printf("Warning extrapolating out of range  k=%g  kmin=%g\n",k,kmin);
            }
            int x1=x0+1;
            int x2=x1+1;
            double k2=x2*h+kmin;

            if (x2>=Niter){ 
                x2=Niter-1;
                x1=x2-1;
                x0=x1-1;
                k2=kmax;
            }
            double k0=x0*h+kmin;
            double k1=x1*h+kmin;
            
            
            
            M[0][0]=k0*k0;  M[0][1]=k0;   M[0][2]=1.;
            M[1][0]=k1*k1;  M[1][1]=k1;   M[1][2]=1.;
            M[2][0]=k2*k2;  M[2][1]=k2;   M[2][2]=1.;
            y[0]=grid[e][n][im][x0]; y[1]=grid[e][n][im][x1];  y[2]=grid[e][n][im][x2];
    //         printf("im=%d  mass=%g L=%d e=%d n=%d x=(%d,%d,%d)   k=(%g,%g,%g)    Z= %g   %g   %g    k=(%g,%g,%g)\n", im, m[e][im],inL,e,n,x0,x1,x2,k0,k1,k2,y[0],y[1],y[2],  kint[e][n][im][x0],kint[e][n][im][x1],kint[e][n][im][x2]);
            double *P=LU_decomposition_solver(3, M, y  );
            
            r[im-m_min]=P[0]*k*k+P[1]*k+P[2];
        
            free(P);
            
        }
        M[0][0]=m[e][m_min]*m[e][m_min];      M[0][1]=m[e][m_min];     M[0][2]=1.;
        M[1][0]=m[e][m_min+1]*m[e][m_min+1];  M[1][1]=m[e][m_min+1];   M[1][2]=1.;
        M[2][0]=m[e][m_min+2]*m[e][m_min+2];  M[2][1]=m[e][m_min+2];   M[2][2]=1.;
        y[0]=r[0]; y[1]=r[1];  y[2]=r[2];
        double *P=LU_decomposition_solver(3, M, y  );
        
        z=P[0]*mass*mass+P[1]*mass+P[2];
        free(P);
        
        
    }
    else { //e==-1 interpolate in L // inL is not a double
//         double z1[2];
//         int dvec[3]={mom[n][0],mom[n][1],mom[n][2]};
//         double ECM2=4*(k+mass*mass);
//         double E2=ECM2+(dvec[0]*dvec[0]+dvec[1]*dvec[1]+dvec[2]*dvec[2]);
//         double gamma=sqrt(E2/ECM2);
//         dzeta_function(z1,  k ,0 , 0, dvec, gamma, 1. , 1.e-3, 1.e6 ,2);
//         return z1[0];
        
        int Lmin=0;
        for(int iL=1; iL<(etot-2);iL++) {
            if (inL> L[iL] )
                Lmin++;
        }
        double zL[3];
        int Lmax=Lmin+3;
        for(int iL=Lmin;iL<Lmax;iL++){
            //select masses 
            int  m_min=0;
            for(int im=1; im<(Nm-2);im++) {
                if (mass> m[iL][im] )
                    m_min++;
            }
            for(int  im=m_min;im<m_min+3;im++){
                double kmin=krange[iL][n][im][0];
                double kmax=krange[iL][n][im][1];
                //error(k<kmin || k>kmax  , 3,"zeta_interpolation::compute"," k=%g out of range [%g,%g ] ",k,kmin,kmax);
                
                // select k for interpolation
                int Niter=Nh;//(kmax-kmin) /h+1 ;
                double h=(kmax-kmin)/Nh;
                
                int x0=(k-kmin)/h ;
                if (x0<0){// if on the edege use the first 3 point
                    x0=0; 
                    //             printf("Warning extrapolating out of range  k=%g  kmin=%g\n",k,kmin);
                }
                int x1=x0+1;
                int x2=x1+1;
                double k2=x2*h+kmin;
                
                if (x2>=Niter){ 
                    x2=Niter-1;
                    x1=x2-1;
                    x0=x1-1;
                    k2=kmax;
                }
                double k0=x0*h+kmin;
                double k1=x1*h+kmin;
                
                
                M[0][0]=k0*k0;  M[0][1]=k0;   M[0][2]=1.;
                M[1][0]=k1*k1;  M[1][1]=k1;   M[1][2]=1.;
                M[2][0]=k2*k2;  M[2][1]=k2;   M[2][2]=1.;
                y[0]=grid[iL][n][im][x0]; y[1]=grid[iL][n][im][x1];  y[2]=grid[iL][n][im][x2];
//                 printf("im=%d  mass=%g L=%d  n=%d x=(%d,%d,%d)   k=(%g,%g,%g)    Z= %g   %g   %g   \n", im, m[iL][im],L[iL],n,x0,x1,x2,k0,k1,k2,y[0],y[1],y[2]);
                double *P=LU_decomposition_solver(3, M, y  );
                
                r[im-m_min]=P[0]*k*k+P[1]*k+P[2];
//                 printf("interpl k=%g -> %g\n",k,r[im-m_min]);
                free(P);
                
            }// loop in mass
            
            
            M[0][0]=m[iL][m_min]*m[iL][m_min];      M[0][1]=m[iL][m_min];     M[0][2]=1.;
            M[1][0]=m[iL][m_min+1]*m[iL][m_min+1];  M[1][1]=m[iL][m_min+1];   M[1][2]=1.;
            M[2][0]=m[iL][m_min+2]*m[iL][m_min+2];  M[2][1]=m[iL][m_min+2];   M[2][2]=1.;
            y[0]=r[0]; y[1]=r[1];  y[2]=r[2];
            double *P=LU_decomposition_solver(3, M, y  );
            
            zL[iL-Lmin]=P[0]*mass*mass+P[1]*mass+P[2];
            free(P);
//             printf("interpl m=%g -> %g\n",mass,zL[iL-Lmin]);
            
        } // loop in L
        M[0][0]=L[Lmin]*L[Lmin];      M[0][1]=L[Lmin];     M[0][2]=1.;
        M[1][0]=L[Lmin+1]*L[Lmin+1];  M[1][1]=L[Lmin+1];   M[1][2]=1.;
        M[2][0]=L[Lmin+2]*L[Lmin+2];  M[2][1]=L[Lmin+2];   M[2][2]=1.;
        y[0]=zL[0]; y[1]=zL[1];  y[2]=zL[2];
        double *P=LU_decomposition_solver(3, M, y  );
        
        z=P[0]*inL*inL+P[1]*inL+P[2];
//         printf("final L=%g -> Z=%g   \n",inL,z);
        free(P);
        
    }// end if 
    
    free(y);
    free_2(3,M);
    return z;
    
}
zeta_interpolation::~zeta_interpolation(){
    if(allocated==0){
        free_4(etot,ntot,3,grid);
        free_2(etot,m);
        free_4(etot,ntot,3,krange);
//         free_4(etot,ntot,3,kint);
    }
}
    /*
} ;*/

    
void zeta_interpolation::write(){
    FILE *f=open_file("zeta_interpolation.dat","w+");
    double a=timestamp();
    size_t ir;
    double z[2];
    ir=fwrite(&etot,sizeof(int),1,f );
    ir=fwrite(&ntot,sizeof(int),1,f );
    ir=fwrite(&Nh,sizeof(double),1,f );
    ir=fwrite(&Nm,sizeof(int),1,f );
    
    for (int i=0;i<ntot;i++){
        ir=fwrite(&mom[i][0],sizeof(double),1,f );
        ir=fwrite(&mom[i][2],sizeof(double),1,f );
        ir=fwrite(&mom[i][3],sizeof(double),1,f );
    }
    printf("writing zeta grid for interpolation \n");
    for (int e=0; e<etot;e++){
        ir=fwrite(&L[e],sizeof(int),1,f );        
             
        for (int im=0; im<Nm;im++)
            ir=fwrite(&m[e][im],sizeof(double),1,f );
            
        for (int d=0; d<ntot;d++){
            for (int im=0; im<Nm;im++){
                ir=fwrite(&krange[e][d][im][0],sizeof(double),1,f );
                ir=fwrite(&krange[e][d][im][1],sizeof(double),1,f );
                    
                int iter=Nh;//(krange[e][d][im][1]-krange[e][d][im][0]) /h+1 ;
                for (int ik=0;ik< iter;ik++){
                    ir=fwrite(&grid[e][d][im][ik],sizeof(double),1,f );
//                     ir=fwrite(&kint[e][d][im][ik],sizeof(double),1,f );
                }
                
                
            }
        }
        
    }
    double b=timestamp();
    printf("time to write zeta_interpolation=%g s\n",b-a);
}
    

    
void zeta_interpolation::read(){
    error(allocated==0,1,"zeta_interpolation::read", "zeta grid already allocated");
    FILE *f=open_file("zeta_interpolation.dat","r+");
    double a=timestamp();
    size_t ir;
    ir=fread(&etot,sizeof(int),1,f );
    ir=fread(&ntot,sizeof(int),1,f );
    ir=fread(&Nh,sizeof(double),1,f );
    ir=fread(&Nm,sizeof(int),1,f );
    
    std::vector< std::vector<int> > momenta(ntot,std::vector<int>(3));
    for (int i=0;i<ntot;i++){
        ir=fread(&momenta[i][0],sizeof(double),1,f );
        ir=fread(&momenta[i][2],sizeof(double),1,f );
        ir=fread(&momenta[i][3],sizeof(double),1,f );
    }
    mom=momenta;
    
    krange=double_malloc_4(etot,ntot,Nm,2);
    m=double_malloc_2(etot,Nm);
    grid=(double****)  malloc(sizeof(double***)*etot);
    printf("reading zeta grid for interpolation \n");
    for (int e=0; e<etot;e++){
        int Ltmp;
        ir=fread(&Ltmp,sizeof(int),1,f );        
        L.emplace_back(Ltmp );
        
        for (int im=0; im<Nm;im++)
            ir=fread(&m[e][im],sizeof(double),1,f );
        grid[e]=(double***)  malloc(sizeof(double**)*ntot);
        for (int d=0; d<ntot;d++){
            grid[e][d]=(double**)  malloc(sizeof(double*)*Nm);
            for (int im=0; im<Nm;im++){
                ir=fread(&krange[e][d][im][0],sizeof(double),1,f );
                ir=fread(&krange[e][d][im][1],sizeof(double),1,f );
                    
                int iter=Nh;//(krange[e][d][im][1]-krange[e][d][im][0]) /h+1 ;
                grid[e][d][im]=(double*)  malloc(sizeof(double)*iter);
                for (int ik=0;ik< iter;ik++){
                    ir=fread(&grid[e][d][im][ik],sizeof(double),1,f );
//                     ir=fread(&kint[e][d][im][ik],sizeof(double),1,f );
                }
            }
        }
    }
    double b=timestamp();
    printf("time to read zeta_interpolation=%g s\n",b-a);
    allocated=0;
}
    
    
    
////////////////////////////////////////////////////////////////////////////////////////////////////
void zeta_interpolation_qsqg::Init(){
    double a=timestamp();
    double z[2];
    int tmp_vec[3]={0,0,0};
    // init the workspace of the dzeta_function
    dzeta_function(z,  1 ,0 , 0, tmp_vec, 1, 1. , 1.e-3, 1.e6 ,1);
    
    std::vector< std::vector<int> > momenta(5,std::vector<int>(3));
    momenta[0][0]=0; momenta[0][1]=0; momenta[0][2]=0; 
    momenta[1][0]=1; momenta[1][1]=0; momenta[1][2]=0; 
    momenta[2][0]=1; momenta[2][1]=1; momenta[2][2]=0; 
    momenta[4][0]=1; momenta[4][1]=1; momenta[4][2]=1; 
    momenta[3][0]=0; momenta[3][1]=0; momenta[3][2]=0; 
    
    mom=momenta;
    
    Ng=maxg/hg;
    qsqrange=double_malloc_3(ntot,Ng,2);
    grid=(double***) malloc(sizeof(double**)*ntot);
    for (int n=0; n<ntot;n++){
        grid[n]=(double**) malloc(sizeof(double*)*Ng);
        double normd=mom[n][0]*mom[n][0]/4. + mom[n][1]*mom[n][1]/4. + mom[n][2]*mom[n][2]/4.;
        int dvec[3];
        dvec[0]=mom[n][0]; dvec[1]=mom[n][1]; dvec[2]=mom[n][2];
        for (int ig=0; ig<Ng;ig++){
            double gamma=1+ig*hg;
            qsqrange[n][ig][0]=normd/(gamma*gamma)+1e-8;
            qsqrange[n][ig][1]=1.+normd/(gamma*gamma)-1e-8;
            if(n==4){
                qsqrange[n][ig][0]=1+normd/(gamma*gamma)+1e-8;
                qsqrange[n][ig][1]=2.+normd/(gamma*gamma)-1e-8; 
            }
            
            int iter=(qsqrange[n][ig][1]-qsqrange[n][ig][0]) /hqsq+1 ;
            error(iter<3,1,"zeta_interpolation::Init","points to evaluate the zeta in qsq are %d to little to interpolate try decreasing hqsq",iter);            
            grid[n][ig]=(double*) malloc(sizeof(double)*iter);
            #pragma omp parallel for  private(z)  shared(hqsq,iter,dvec,n,ig,qsqrange,grid) 
            for (int iqsq=0;iqsq< iter;iqsq++){
                double qsq=qsqrange[n][ig][0]+hqsq*iqsq;
                if(iqsq==iter-1) qsq=qsqrange[n][ig][1];
                //                     printf("k=%g   gamma=%g  dvec=(%d,%d,%d) mass=%g  e=%d\n",k,gamma,dvec[0],dvec[1],dvec[2], mass[ens[e]],ens[e]);
                dzeta_function(z,  qsq ,0 , 0, dvec, gamma, 1. , 1.e-3, 1.e6 ,3);
                grid[n][ig][iqsq]=z[0]; //real part
            }
            for (int iqsq=iter-10;iqsq< iter;iqsq++){
                error(grid[n][ig][0]*grid[n][ig][iqsq]>0,1,"zeta_interpolation::Init","over the pole    n=%d |(%d,%d,%d)|^2=%g   ig=%d   gamma=%g  qsq=[%g,%g]\n",n,mom[n][0],mom[n][1],mom[n][2],normd,ig,gamma,qsqrange[n][ig][0],qsqrange[n][ig][1]);
            }
        }
    }
    double b=timestamp();
    allocated=0;
    printf("time to Init zeta_interpolation=%g s\n",b-a);
    
}
  
  
double zeta_interpolation_qsqg::compute( double qsq, int n,  double gamma ){
    double r[3];
    double **M, *y;
    M=double_malloc_2(3,3);
    y=(double*) malloc(sizeof(double)*3);
    
    double gmax=maxg;
    double gmin=1.;
    int ig[3];
    double g[3];
    
    // select g for interpolation
    ig[0]=(gamma-gmin)/hg ;
    error(ig[0]<0 || ig[0]>=Ng ,1,"zeta_interpolation_qsqg::compute", "gamma=%g out of range [%g,%g] \n",gamma,gmin,gmax);
    ig[1]=ig[0]+1;
    ig[2]=ig[1]+1;
    g[2]=ig[2]*hg+gmin;
    
    if (ig[2]>=Ng){ 
        ig[2]=Ng-1;
        ig[1]=ig[2]-1;
        ig[0]=ig[1]-1;
        g[2]=gmax;
    }
    g[0]=ig[0]*hg+gmin;
    g[1]=ig[1]*hg+gmin;
    for (int i=0; i<3 ;i++){
        double qsqmin=qsqrange[n][ig[i]][0];
        double qsqmax=qsqrange[n][ig[i]][1];
        //error(k<kmin || k>kmax  , 3,"zeta_interpolation::compute"," k=%g out of range [%g,%g ] ",k,kmin,kmax);
        
        // select k for interpolation
        int Niter=(qsqmax-qsqmin) /hqsq+1 ;
        int x0=(qsq-qsqmin)/hqsq ;
        if (x0<0){// if on the edege use the first 3 point
            x0=0; 
        }
        int x1=x0+1;
        int x2=x1+1;
        double qsq2=x2*hqsq+qsqmin;
        
        if (x2>=Niter){ 
            x2=Niter-1;
            x1=x2-1;
            x0=x1-1;
            qsq2=qsqmax;
        }
        double qsq0=x0*hqsq+qsqmin;
        double qsq1=x1*hqsq+qsqmin;
        
        M[0][0]=qsq0*qsq0;  M[0][1]=qsq0;   M[0][2]=1.;
        M[1][0]=qsq1*qsq1;  M[1][1]=qsq1;   M[1][2]=1.;
        M[2][0]=qsq2*qsq2;  M[2][1]=qsq2;   M[2][2]=1.;
        y[0]=grid[n][ig[i]][x0]; y[1]=grid[n][ig[i]][x1];  y[2]=grid[n][ig[i]][x2];
        //         printf("im=%d  mass=%g L=%d e=%d n=%d x=(%d,%d,%d)   k=(%g,%g,%g)    Z= %g   %g   %g    k=(%g,%g,%g)\n", im, m[e][im],inL,e,n,x0,x1,x2,k0,k1,k2,y[0],y[1],y[2],  kint[e][n][im][x0],kint[e][n][im][x1],kint[e][n][im][x2]);
        double *P=LU_decomposition_solver(3, M, y  );
        
        r[i]=P[0]*qsq*qsq+P[1]*qsq+P[2];
        
        free(P);
        
    }
    M[0][0]=g[0]*g[0];  M[0][1]=g[0];     M[0][2]=1.;
    M[1][0]=g[1]*g[1];  M[1][1]=g[1];   M[1][2]=1.;
    M[2][0]=g[2]*g[2];  M[2][1]=g[2];   M[2][2]=1.;
    y[0]=r[0]; y[1]=r[1];  y[2]=r[2];
    double *P=LU_decomposition_solver(3, M, y  );
    
    double z=P[0]*gamma*gamma+P[1]*gamma+P[2];
    free(P);
        
    free(y);
    free_2(3,M);
    return z;
}      

zeta_interpolation_qsqg::~zeta_interpolation_qsqg(){
    if(allocated==0){
    free_3(ntot,Ng,grid);
    free_3(ntot,Ng,qsqrange);
    }
}
    
    
