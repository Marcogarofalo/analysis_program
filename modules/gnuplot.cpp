#define gnuplot_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "linear_fit.hpp"
#include "mutils.hpp"
#include "resampling.hpp"
#include "gnuplot.hpp"
#include "global.hpp"

 


 
void plotting(int T, double **corr,int *tmin,int *tmax  , int *sep){
    FILE * temp = fopen("data.temp", "w");
    int i,j=0;
    char c[NAMESIZE];
    FILE * gnuplotPipe = popen ("gnuplot ", "w"); //-p if you want the graph to stay after you close gnuplot
    
    for(i=1;i<T/2;i++){
        if (isinf(corr[i][0])!=1 && corr[i][0]==corr[i][0])
            fprintf(temp,"%d \t %.15g \t %.15g\n",i,corr[i][0],corr[i][1]);
        else fprintf(temp,"%d \t %.15g \t %.15g\n",i,0.0,0.0);
    }
    fclose(temp);
    fprintf(gnuplotPipe, "plot 'data.temp' i 0 u 1:2:3 w e \n");
    fflush(gnuplotPipe);

    //system("gnuplot\n load 'gnuplot.gp'\n");  // system("shell command");
    
    fprintf(stderr,"Please enter the plateau interval and separation:\n"); 
    
    myscanf(3, (char*) "%d %d %d",tmin,tmax,sep);
//printf("%d  %d  %d\n", *tmin,*tmax,*sep);
/*    while (j!=3){
        printf("the entry must be numbers\n");
        j=scanf("%d   %d   %d",tmin, tmax,sep);
        if (j!=3)
            for (i=0;i<3-j;i++) scanf("%s",c);
        printf("%d\n",j);
    }*/
    while (*tmin>*tmax && *tmax>= T/2 &&  *tmin <0){
            fprintf(stderr,"please enter a valid number for tmin and tmax\n");
            myscanf(3,(char*)"%d %d %d",tmin,tmax,sep);
    }
    
    
    
    fprintf(gnuplotPipe,"exit \n");
  
    pclose(gnuplotPipe);

}

void plotting_fit_p(int T, double **corr,int tmin,int tmax, double *fit ){
    FILE * temp = fopen("data.temp", "w");
    int i,r=-1;
    char answer[NAMESIZE];
    FILE * gnuplotPipe = popen ("gnuplot -p", "w"); 
    
    for(i=1;i<T/2;i++){
        fprintf(temp,"%d \t %.15g \t %.15g\n",i,corr[i][0],corr[i][1]);
    }
    fclose(temp);
    fprintf(gnuplotPipe,"fm(x)= x<(%d-0.1) ? 1/0 : x>(%d+0.1) ? 1/0 : %g - %g \n",tmin,tmax,fit[0],fit[1]);
    fprintf(gnuplotPipe,"f(x)= x<(%d-0.1) ? 1/0 : x>(%d+0.1) ? 1/0 : %g  \n",tmin,tmax,fit[0]);
    fprintf(gnuplotPipe,"fp(x)= x<(%d-0.1) ? 1/0 : x>(%d+0.1) ? 1/0 : %g + %g \n",tmin,tmax,fit[0],fit[1]);

    fprintf(gnuplotPipe,"set sample 200\n");
    fprintf(gnuplotPipe, "plot 'data.temp' i 0 u 1:2:3 w e, fm(x) t '' lc 7, f(x) t '' lc 7, fp(x) t '' lc 7\n");
    fflush(gnuplotPipe);

    
    fprintf(gnuplotPipe,"exit \n");
  
    pclose(gnuplotPipe);
    
    

}


int plotting_fit(int T, double **corr,int tmin,int tmax, double *fit,double *chi2 ){
    FILE * temp = fopen("data.temp", "w");
    int i,r=-1;
    char answer[NAMESIZE],p[NAMESIZE];
    FILE * gnuplotPipe = popen ("gnuplot ", "w"); 
    
    for(i=1;i<T/2;i++){
        fprintf(temp,"%d \t %.15g \t %.15g\n",i,corr[i][0],corr[i][1]);
    }
    fclose(temp);
    fprintf(gnuplotPipe,"fm(x)= x<(%d-0.1) ? 1/0 : x>(%d+0.1) ? 1/0 : %g - %g \n",tmin,tmax,fit[0],fit[1]);
    fprintf(gnuplotPipe,"f(x)= x<(%d-0.1) ? 1/0 : x>(%d+0.1) ? 1/0 : %g  \n",tmin,tmax,fit[0]);
    fprintf(gnuplotPipe,"fp(x)= x<(%d-0.1) ? 1/0 : x>(%d+0.1) ? 1/0 : %g + %g \n",tmin,tmax,fit[0],fit[1]);

    fprintf(gnuplotPipe,"set title 'chi2=%g'\n",chi2[0]);
    fprintf(gnuplotPipe,"set sample 200\n");
    fprintf(gnuplotPipe, "plot 'data.temp' i 0 u 1:2:3 w e, fm(x) t '' lc 7, f(x) t '' lc 7, fp(x) t '' lc 7\n");
    fflush(gnuplotPipe);

    
    fprintf(stderr,"is the fit ok? [yes or no]\n");
    
    scanf(" %s",answer);
    
    while (r<0){
        if ( strcmp(answer,"yes")==0 ){
            r=0;
            system("pkill gnuplot_qt");
        }
        else if ( strcmp(answer,"y")==0 ){
            r=0;
            system("pkill gnuplot_qt");
        }
        else if ( strcmp(answer,"n")==0 ){
            r=1;
            plotting_fit_p( T, corr , tmin,tmax,fit);
        }
        else if ( strcmp(answer,"no")==0 ){
            r=1;
            plotting_fit_p( T, corr , tmin,tmax,fit);
        }
        else{
            printf("please enter a valid answer:\n");
            scanf("%s",answer);
        }
    }
    
    fprintf(gnuplotPipe,"exit \n");
  
    pclose(gnuplotPipe);
    
    return r;

}



void contraction_name(int ii,char *r);

void plotting_fit_pdf(char **argv, const char *string,int L, int T, double beta, double k_sea,double mu_sea, double **corr,int tmin,int tmax, double *fit,const char *contraction ,struct kinematic kinematic_2pt){
    FILE * temp = fopen("data.temp", "w");
    int i;
    char answer[NAMESIZE],p[NAMESIZE];
    FILE * gnuplotPipe = popen ("gnuplot ", "w"); 
    char name[NAMESIZE],instruction[NAMESIZE];
    double t2p=2*3.1415926536*sqrt(3)/((double)T);
    double k2=kinematic_2pt.k2;
    int    r2=kinematic_2pt.r2;
    double k1=kinematic_2pt.k1;
    int    r1=kinematic_2pt.r1;
    double mom2=kinematic_2pt.mom2;
    double mom1=kinematic_2pt.mom1;
    int ik2=kinematic_2pt.ik2;
    int ik1=kinematic_2pt.ik1;
    
    for(i=1;i<T/2;i++){
        fprintf(temp,"%d \t %.15g \t %.15g\n",i,corr[i][0],corr[i][1]);
    }
    fclose(temp);
    fprintf(gnuplotPipe,"fm(x)= x<(%d-0.1) ? 1/0 : x>(%d+0.1) ? 1/0 : %g - %g \n",tmin,tmax,fit[0],fit[1]);
    fprintf(gnuplotPipe,"f(x)= x<(%d-0.1) ? 1/0 : x>(%d+0.1) ? 1/0 : %g  \n",tmin,tmax,fit[0]);
    fprintf(gnuplotPipe,"fp(x)= x<(%d-0.1) ? 1/0 : x>(%d+0.1) ? 1/0 : %g + %g \n",tmin,tmax,fit[0],fit[1]);


    fprintf(gnuplotPipe,"set sample 200\n");
    fprintf(gnuplotPipe,"set term epslatex  color colortext standalone size 16.5cm,11.5cm linewidth 2 \n");
    fprintf(gnuplotPipe,"set linetype 1 lw 1.5 pt 5 ps 2 lc 'blue'\n");
    fprintf(gnuplotPipe,"set linetype 2 lw 1.5 pt 7 ps 2 lc 'red'\n");
    fprintf(gnuplotPipe,"set mytics 2\n");
    fprintf(gnuplotPipe,"set mxtics 2\n");
    fprintf(gnuplotPipe,"set grid\n");
    fprintf(gnuplotPipe,"set key spacing 1.2\n");

    mysprintf(name,NAMESIZE,"%s_%s_L%dT%d_b%.4f_ks%.7f_mus%.4f_mu2_%.4f_r2_%d_mu1_%.4f_r1_%d_theta1_%.4f_theta2_%.4f",string,contraction,L,T,beta,k_sea,mu_sea,k2,r2,k1,r1,mom2,mom1);
    fprintf(gnuplotPipe,"set output '%s.tex'\n",name);
    
    fprintf(gnuplotPipe,"set title '%d-%d, L%dT%d  %s $\\beta=%.4f$ $\\kappa_{sea}=%.7f$  $\\mu_{sea}=%.4f$'\n",ik2,ik1,L,T,contraction,beta,k_sea,mu_sea);
    int tmax_plot,tmin_plot;
    if (tmin-T/8 <0) tmin_plot=0;
    else tmin_plot=tmin-T/8;
    if (tmax+T/12 >T/2) tmax_plot=T/2;
    else tmax_plot=tmax+T/12;
       
    double ymax_plot,ymin_plot;
    ymax_plot=fit[0]+35*fit[1];
    ymin_plot=fit[0]-25*fit[1];
    fprintf(gnuplotPipe,"set xrange [%d-0.1:%d+0.1]\n",tmin_plot,tmax_plot);
    fprintf(gnuplotPipe,"set yrange [%f:%f]\n",ymin_plot,ymax_plot);
    fprintf(gnuplotPipe, "plot 'data.temp' i 0 u 1:2:3 w e t '$a\\mu_2=%.4f$,  $a\\mu_1=%.4f$, $a\\theta_2=%.4f$, $a\\theta_1=%.4f$ ', fm(x) t ''  lc 7, f(x) t '$%s=%g\\pm%.2g$' lc 7, fp(x) t '' lc 7\n",k2,k1,mom2,mom1,string,fit[0],fit[1]);
    fflush(gnuplotPipe);

   
    
    fprintf(gnuplotPipe,"exit \n");
  
    pclose(gnuplotPipe);
    mysprintf(instruction,NAMESIZE,"pdflatex -output-directory=%s/out %s.tex   >> /dev/null  ",argv[3],name);
    system(instruction);
    
    mysprintf(instruction,NAMESIZE,"rm %s.tex %s/out/%s.log %s/out/%s.aux  %s-inc.eps %s-inc-eps-converted-to.pdf ",name,argv[3],name,argv[3],name,name,name);
    system(instruction);
    
} 


void plotting_fit_pdf_G(char **argv, const char *string ,double **corr,int tmin,int tmax, double *fit,const char *contraction , struct kinematic_G kinematic_2pt_G ){
    
    int L=file_head.l1, T=file_head.l0;
    double beta=file_head.beta, k_sea=file_head.ksea, mu_sea=file_head.musea;
  
    FILE * temp = fopen("data.temp", "w");
    int i;
    char answer[NAMESIZE],p[NAMESIZE];
    FILE * gnuplotPipe = popen ("gnuplot ", "w"); 
    char name[NAMESIZE],instruction[NAMESIZE];
    double t2p=2*3.1415926536/((double)T);
    
    
     
    
    for(i=1;i<T/2;i++){
        if (corr[i][0]!=corr[i][0]) return ;
        fprintf(temp,"%d \t %.15g \t %.15g\n",i,corr[i][0],corr[i][1]);
    }
    fclose(temp);
    fprintf(gnuplotPipe,"fm(x)= x<(%d-0.1) ? 1/0 : x>(%d+0.1) ? 1/0 : %g - %g \n",tmin,tmax,fit[0],fit[1]);
    fprintf(gnuplotPipe,"f(x)= x<(%d-0.1) ? 1/0 : x>(%d+0.1) ? 1/0 : %g  \n",tmin,tmax,fit[0]);
    fprintf(gnuplotPipe,"fp(x)= x<(%d-0.1) ? 1/0 : x>(%d+0.1) ? 1/0 : %g + %g \n",tmin,tmax,fit[0],fit[1]);


    fprintf(gnuplotPipe,"set sample 200\n");
    fprintf(gnuplotPipe,"set term epslatex  color colortext standalone size 16.5cm,11.5cm linewidth 2 \n");
    fprintf(gnuplotPipe,"set linetype 1 lw 1.5 pt 5 ps 2 lc 'blue'\n");
    fprintf(gnuplotPipe,"set linetype 2 lw 1.5 pt 7 ps 2 lc 'red'\n");
    fprintf(gnuplotPipe,"set mytics 2\n");
    fprintf(gnuplotPipe,"set key spacing 1.2\n");

    mysprintf(name,NAMESIZE,"%s_%s_L%dT%d_b%.4f_ks%.7f_mus%.4f_mut_%.4f_rt_%d_mus_%.4f_rs_%d_theta0_%.4f_thetat_%.4f_thetas_%.4f",string,contraction,L,T,beta,k_sea,mu_sea,
            kinematic_2pt_G.kt,kinematic_2pt_G.rt,kinematic_2pt_G.ks,kinematic_2pt_G.rs,
            kinematic_2pt_G.Mom0[3],kinematic_2pt_G.Momt[3],kinematic_2pt_G.Moms[3]);
    fprintf(gnuplotPipe,"set output '%s.tex'\n",name);
    
    fprintf(gnuplotPipe,"set title '$L%dT%d$  $%s$ $\\beta=%.4f$ $\\kappa_{sea}=%.7f$  $\\mu_{sea}=%.4f$'\n",L,T,contraction,beta,k_sea,mu_sea);
    int tmax_plot,tmin_plot;
    if (tmin-8 <0) tmin_plot=0;
    else tmin_plot=tmin-8;
    if (tmax+4 >T/2) tmax_plot=T/2;
    else tmax_plot=tmax+4;
    
    double ymax_plot,ymin_plot;
    ymax_plot=fit[0]+20*fit[1];
    ymin_plot=fit[0]-20*fit[1];
    fprintf(gnuplotPipe,"set xrange [%d-0.1:%d+0.1]\n",tmin_plot,tmax_plot);
    fprintf(gnuplotPipe,"set yrange [%f:%f]\n",ymin_plot,ymax_plot);
    fprintf(gnuplotPipe, "plot 'data.temp' i 0 u 1:2:3 w e t '$a\\mu_t=%.4f$,  $a\\mu_s=%.4f$, $ap_0=%.4f$, $ap_t=%.4f$ $ap_s=%.4f$ ', fm(x) t ''  lc 7, f(x) t '$%s=%g\\pm%.2g$' lc 7, fp(x) t '' lc 7\n",
              kinematic_2pt_G.kt, kinematic_2pt_G.ks,
            kinematic_2pt_G.Mom0[3]*t2p,kinematic_2pt_G.Momt[3]*t2p,kinematic_2pt_G.Moms[3]*t2p,
            string,fit[0],fit[1]);
    fflush(gnuplotPipe);

   
    
    fprintf(gnuplotPipe,"exit \n");
  
    pclose(gnuplotPipe);
    mysprintf(instruction,NAMESIZE,"pdflatex -output-directory=%s/out %s.tex  > /dev/null &",argv[3],name);
    system(instruction);
    
    mysprintf(instruction,NAMESIZE,"rm %s.tex %s/out/%s.log %s/out/%s.aux  %s-inc.eps %s-inc-eps-converted-to.pdf",name,argv[3],name,argv[3],name,name,name);
    system(instruction);
    
} 

void plotting_fit_deg_n_pdf(char **argv, const char *string,int npoints , int deg ,double *x, double **y, double **fit,double *chi2 ){
    FILE * temp = fopen("data.temp", "w");
    int i;
    char answer[NAMESIZE],p[NAMESIZE];
    FILE * gnuplotPipe = popen ("gnuplot ", "w"); 
    char name[NAMESIZE],instruction[NAMESIZE];
    char polynomial[NAMESIZE];
   
    for(i=0;i<npoints;i++){
        fprintf(temp,"%g \t %.15g \t %.15g\n",x[i],y[i][0],y[i][1]);
    }
    fclose(temp);
    fprintf(gnuplotPipe,"f(x)= ");
    for(i=0;i<=deg;i++)
        fprintf(gnuplotPipe,"+( %g )*(x**%d) ",fit[i][0],i);
    fprintf(gnuplotPipe,"\n");


    fprintf(gnuplotPipe,"set sample 200\n");
    fprintf(gnuplotPipe,"set key left \n");
    fprintf(gnuplotPipe,"set term epslatex  color colortext standalone size 16.5cm,11.5cm linewidth 2 \n");
    fprintf(gnuplotPipe,"set linetype 1 lw 1.5 pt 5 ps 1 lc 'blue'\n");
    fprintf(gnuplotPipe,"set linetype 2 lw 1.5 pt 7 ps 1 lc 'red'\n");
    fprintf(gnuplotPipe,"set mytics 2\n");
    fprintf(gnuplotPipe,"set key spacing 1.2\n");

    mysprintf(name,NAMESIZE,"%s_fit_deg%d_L%dT%d_b%.4f_ks%.7f_mus%.4f",string,deg,file_head.l1,file_head.l0,file_head.beta,file_head.ksea,file_head.musea);
    fprintf(gnuplotPipe,"set output '%s.tex'\n",name);
    
    mysprintf(polynomial,NAMESIZE," $(%g\\pm %.2g)",fit[0][0],fit[0][1]);
    for(i=1;i<=deg;i++)
        mysprintf(polynomial,NAMESIZE,"%s+( %g\\pm %.2g )x^{%d} ",polynomial,fit[i][0],fit[i][1],i);
    mysprintf(polynomial,NAMESIZE,"%s \t$",polynomial);

    fprintf(gnuplotPipe,"set title 'L%dT%d  $%s$ $\\beta=%.4f$ $\\kappa_{sea}=%.7f$  $\\mu_{sea}=%.4f$'\n",file_head.l1,file_head.l0,string,file_head.beta,file_head.ksea,file_head.musea);
    printf("%s",polynomial);
    fprintf(gnuplotPipe, "plot 'data.temp' i 0 u 1:2:3 w e t '$%s$', f(x) t '%s' \n",string,polynomial);
    fflush(gnuplotPipe);

    fprintf(gnuplotPipe,"exit \n");
    pclose(gnuplotPipe);
    mysprintf(instruction,NAMESIZE,"pdflatex -output-directory=%s/out %s.tex > null ",argv[3],name);
    system(instruction);
    
    mysprintf(instruction,NAMESIZE,"rm %s.tex %s/out/%s.log %s/out/%s.aux  %s-inc.eps %s-inc-eps-converted-to.pdf",name,argv[3],name,argv[3],name,name,name);
    system(instruction);
  
} 


void plotting_fit_poles_deg_n_pdf(char **argv, const char *string,int npoints , int deg ,double *x, double **y, double **fit,double *chi2 ){
    FILE * temp = fopen("data.temp", "w");
    int i;
    char answer[NAMESIZE],p[NAMESIZE];
    FILE * gnuplotPipe = popen ("gnuplot ", "w"); 
    char name[NAMESIZE],instruction[NAMESIZE];
    char polynomial[NAMESIZE];
    
    for(i=0;i<npoints;i++){
        fprintf(temp,"%g \t %.15g \t %.15g\n",x[i],y[i][0],y[i][1]);
    }
    fclose(temp);
    fprintf(gnuplotPipe,"f(x)= ");
    for(i=0;i<=deg;i++)
        fprintf(gnuplotPipe,"+( %g )*(x**(%d)) ",fit[i][0],-i);
    fprintf(gnuplotPipe,"\n");


    fprintf(gnuplotPipe,"set sample 200\n");
    fprintf(gnuplotPipe,"set key left \n");
    fprintf(gnuplotPipe,"set term epslatex  color colortext standalone size 16.5cm,11.5cm linewidth 2 \n");
    fprintf(gnuplotPipe,"set linetype 1 lw 1.5 pt 5 ps 1 lc 'blue'\n");
    fprintf(gnuplotPipe,"set linetype 2 lw 1.5 pt 7 ps 1 lc 'red'\n");
    fprintf(gnuplotPipe,"set mytics 2\n");
    fprintf(gnuplotPipe,"set key spacing 1.2\n");

    mysprintf(name,NAMESIZE,"%s_fit_deg%d_L%dT%d_b%.4f_ks%.7f_mus%.4f",string,deg,file_head.l1,file_head.l0,file_head.beta,file_head.ksea,file_head.musea);
    fprintf(gnuplotPipe,"set output '%s.tex'\n",name);
    
    mysprintf(polynomial,NAMESIZE," $(%g\\pm %.2g)",fit[0][0],fit[0][1]);
    for(i=1;i<=deg;i++)
        mysprintf(polynomial,NAMESIZE,"%s+( %g\\pm %.2g )x^{-%d} ",polynomial,fit[i][0],fit[i][1],i);
    mysprintf(polynomial,NAMESIZE,"%s$",polynomial);
    
    fprintf(gnuplotPipe,"set title 'L%dT%d  $%s$ $\\beta=%.4f$ $\\kappa_{sea}=%.7f$  $\\mu_{sea}=%.4f$'\n",file_head.l1,file_head.l0,string,file_head.beta,file_head.ksea,file_head.musea);
    fprintf(gnuplotPipe, "plot 'data.temp' i 0 u 1:2:3 w e t '$%s$', f(x) t '%s' \n",string,polynomial);
    fflush(gnuplotPipe);

   
    fprintf(gnuplotPipe,"exit \n");
  
    pclose(gnuplotPipe);
    mysprintf(instruction,NAMESIZE,"pdflatex -output-directory=%s/out %s.tex   > null",argv[3],name);
    system(instruction);
    
    mysprintf(instruction,NAMESIZE,"rm %s.tex %s/out/%s.log %s/out/%s.aux  %s-inc.eps %s-inc-eps-converted-to.pdf",name,argv[3],name,argv[3],name,name,name);
    system(instruction);
    
} 





void contraction_index(int *ii ,const char*r ){

if(  strcmp(r,"V0P5")==0 )   *ii=0;
else if(  strcmp(r,"P5P5")==0 )   *ii=1;

}
 
 
void contraction_name(int ii, char *r){
   
    
if (ii==0)  mysprintf(r,NAMESIZE,"V0P5");
else if (ii==1) mysprintf(r,NAMESIZE,"P5P5"); 
/*
if (ii==0)  mysprintf(r,NAMESIZE,"P5S0");
else if (ii==1)  mysprintf(r,NAMESIZE,"P5V1");
else if (ii==2)  mysprintf(r,NAMESIZE,"P5V2");
else if (ii==3)  mysprintf(r,NAMESIZE,"P5V3");
else if (ii==4)  mysprintf(r,NAMESIZE,"P5V0");
else if (ii==5)  mysprintf(r,NAMESIZE,"P5P5");
else if (ii==6)  mysprintf(r,NAMESIZE,"P5A1");
else if (ii==7)  mysprintf(r,NAMESIZE,"P5A2");
else if (ii==8)  mysprintf(r,NAMESIZE,"P5A3");
else if (ii==9)  mysprintf(r,NAMESIZE,"P5A0");
else if (ii==10)  mysprintf(r,NAMESIZE,"P5T1");
else if (ii==11)  mysprintf(r,NAMESIZE,"P5T2");
else if (ii==12)  mysprintf(r,NAMESIZE,"P5T3");
else if (ii==13)  mysprintf(r,NAMESIZE,"P5B1");
else if (ii==14)  mysprintf(r,NAMESIZE,"P5B2");
else if (ii==15)  mysprintf(r,NAMESIZE,"P5B3");
else if (ii==16)  mysprintf(r,NAMESIZE,"S0P5");
else if (ii==17)  mysprintf(r,NAMESIZE,"V1P5");
else if (ii==18)  mysprintf(r,NAMESIZE,"V2P5");
else if (ii==19)  mysprintf(r,NAMESIZE,"V3P5");
else if (ii==20)  mysprintf(r,NAMESIZE,"V0P5");
else if (ii==21)  mysprintf(r,NAMESIZE,"P5P5");
else if (ii==22)  mysprintf(r,NAMESIZE,"A1P5");
else if (ii==23)  mysprintf(r,NAMESIZE,"A2P5");
else if (ii==24)  mysprintf(r,NAMESIZE,"A3P5");
else if (ii==25)  mysprintf(r,NAMESIZE,"A0P5");
else if (ii==26)  mysprintf(r,NAMESIZE,"T1P5");
else if (ii==27)  mysprintf(r,NAMESIZE,"T2P5");
else if (ii==28)  mysprintf(r,NAMESIZE,"T3P5");
else if (ii==29)  mysprintf(r,NAMESIZE,"B1P5");
else if (ii==30)  mysprintf(r,NAMESIZE,"B2P5");
else if (ii==31)  mysprintf(r,NAMESIZE,"B3P5");
else if (ii==32)  mysprintf(r,NAMESIZE,"S0S0");
else if (ii==33)  mysprintf(r,NAMESIZE,"V0V0");
else if (ii==34)  mysprintf(r,NAMESIZE,"A0A0");
else if (ii==35)  mysprintf(r,NAMESIZE,"V1V1");
else if (ii==36)  mysprintf(r,NAMESIZE,"V2V2");
else if (ii==37)  mysprintf(r,NAMESIZE,"V3V3");
else if (ii==38)  mysprintf(r,NAMESIZE,"A1A1");
else if (ii==39)  mysprintf(r,NAMESIZE,"A2A2");
else if (ii==40)  mysprintf(r,NAMESIZE,"A3A3");
else if (ii==41)  mysprintf(r,NAMESIZE,"T1T1");
else if (ii==42)  mysprintf(r,NAMESIZE,"T2T2");
else if (ii==43)  mysprintf(r,NAMESIZE,"T3T3");
else if (ii==44)  mysprintf(r,NAMESIZE,"V1T1");
else if (ii==45)  mysprintf(r,NAMESIZE,"V2T2");
else if (ii==46)  mysprintf(r,NAMESIZE,"V3T3");
else if (ii==47)  mysprintf(r,NAMESIZE,"T1V1");
else if (ii==48)  mysprintf(r,NAMESIZE,"T2V2");
else if (ii==49)  mysprintf(r,NAMESIZE,"T3V3");
else if (ii==50)  mysprintf(r,NAMESIZE,"B1B1");
else if (ii==51)  mysprintf(r,NAMESIZE,"B2B2");
else if (ii==52)  mysprintf(r,NAMESIZE,"B3B3");
else if (ii==53)  mysprintf(r,NAMESIZE,"A1B1");
else if (ii==54)  mysprintf(r,NAMESIZE,"A2B2");
else if (ii==55)  mysprintf(r,NAMESIZE,"A3B3");
else if (ii==56)  mysprintf(r,NAMESIZE,"B1A1");
else if (ii==57)  mysprintf(r,NAMESIZE,"B2A2");
else if (ii==58)  mysprintf(r,NAMESIZE,"B3A3");
else if (ii==59)  mysprintf(r,NAMESIZE,"V0S0");
else if (ii==60)  mysprintf(r,NAMESIZE,"S0V0");
*/
}



