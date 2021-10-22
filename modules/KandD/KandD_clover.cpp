#define KandD_clover_C

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
#include "fve.hpp"
#include "indices.hpp"
#include "KandD.hpp"
#include "tower.hpp"
#include "mutils.hpp"

#include <unistd.h>

#include <omp.h> 


static double two_parabolas(int n, int Nvar, double *x,int Npar,double  *P){
    double r;
    
    if (n==0)
        r=P[0]+P[1]*x[0]+P[2]*x[0]*x[0];
    if (n==1)
        r=P[3]+P[4]*x[0]+P[5]*x[0]*x[0];
    
    return r;
    
}

static void print_chiral_continuum_fit(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, const char *AV,const char *namefile, struct header *head ,struct data_jack *gJ){
    
    int i,j;
    char name[NAMESIZE];
    FILE *fc=NULL,*fcA=NULL,*fcB=NULL,*fcC=NULL;
    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum.txt",argv[2],namefile);
    fc=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_A.txt",argv[2],namefile);
    fcA=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_B.txt",argv[2],namefile);
    fcB=open_file(name,"w+");
    mysprintf(name,NAMESIZE,"%s/%s_chiral_continuum_C.txt",argv[2],namefile);
    fcC=open_file(name,"w+");
    int N=fit_info.N,n;
    
    double h=0.0003;
    double R=12.;
    //double h=0.003;
    double **tif=swap_indices(fit_info.Npar,jack_tot,fit_out.P);   
    double *r=(double*) malloc(sizeof(double)*jack_tot);
    double *r1=(double*) malloc(sizeof(double)*jack_tot);
    double *m;
    double **x=(double**) malloc(sizeof(double*)*jack_tot);
    for(j=0;j<jack_tot;j++){
        x[j]=(double*) malloc(sizeof(double)*fit_info.Nvar);
        //mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
        x[j][0]=0;//mlw
        x[j][1]=1e+6;//r0
        x[j][2]=1e+8;//Mpi2
        x[j][3]=1e+8;//fpi
        x[j][4]=1e+9;//L
        x[j][5]=phys_point[j][5];
        x[j][6]=phys_point[j][6];
        x[j][7]=phys_point[j][7];
        x[j][8]=phys_point[j][8];
        x[j][9]=phys_point[j][9]; 
        x[j][10]=phys_point[j][10];
        x[j][11]=phys_point[j][11];
        
        
        
    }
    for (i=2;i<100;i++){
        fprintf(fc,"%g \t",((double)i)*h);
        for (n=0;n<N;n++){
            for(j=0;j<jack_tot;j++){
                x[j][0]=((double)i)*h;//=1e+10;//xG
                /*if (j==jack_tot-1){
                    printf("KM=%g  Kf=%g\n",fit_info.function(2,fit_info.Nvar,x[j],fit_info.Npar,tif[j]),fit_info.function(3,fit_info.Nvar,x[j],fit_info.Npar,tif[j]));
                    printf("L_w=%f, mw=%f, fw=%f, Bw=%f\n",x[j][4], x[j][0], x[j][8], x[j][9] );
                }*/
                
                r[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
                x[j][10]=((double)i)*h*R;//=1e+10;//xG
                r1[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
            }
            m=mean_and_error_jack_biased(jack_tot,r);
            fprintf(fc," %g  %g \t",m[0],m[1]);
            free(m);
            m=mean_and_error_jack_biased(jack_tot,r1);
            fprintf(fc,"%g %g  %g \t",x[jack_tot-1][10],m[0],m[1]);
            free(m);
        }
        fprintf(fc,"\n");
    }
    for (i=2;i<100;i++){
        fprintf(fcA,"%g \t",((double)i)*h);
        for (n=0;n<N;n++){
            for(j=0;j<jack_tot;j++){
                x[j][1]=gJ[0].w0[j];
                x[j][0]=((double)i)*h;//=1e+10;//xG
                r[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
                
                x[j][10]=((double)i)*h*R;//=1e+10;//xG
                r1[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
            }
            m=mean_and_error_jack_biased(jack_tot,r);
            fprintf(fcA," %g  %g \t",m[0],m[1]);
            free(m);
            m=mean_and_error_jack_biased(jack_tot,r1);
            fprintf(fcA,"%g %g  %g \t",x[jack_tot-1][10],m[0],m[1]);
            free(m);
        }
       
        fprintf(fcA,"%g \n",gJ[0].w0[jack_tot-1]);
        
        fprintf(fcB,"%g \t",((double)i)*h);
        for (n=0;n<N;n++){
            for(j=0;j<jack_tot;j++){
                x[j][1]=gJ[4].w0[j];
                x[j][0]=((double)i)*h;//=1e+10;//xG
                r[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
                
                x[j][10]=((double)i)*h*R;//=1e+10;//xG
                r1[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
            }
            m=mean_and_error_jack_biased(jack_tot,r);
            fprintf(fcB," %g  %g \t",m[0],m[1]);
            free(m);
            m=mean_and_error_jack_biased(jack_tot,r1);
            fprintf(fcB,"%g %g  %g \t",x[jack_tot-1][10],m[0],m[1]);
            free(m);
        }
        fprintf(fcB,"\n");

       fprintf(fcC,"%g \t",((double)i)*h);
        for (n=0;n<N;n++){
            for(j=0;j<jack_tot;j++){
                x[j][1]=gJ[7].w0[j];
                x[j][2]=((double)i)*h;//=1e+10;//xG
                r[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
                
                x[j][10]=((double)i)*h*R;//=1e+10;//xG
                r1[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
                
            }
            m=mean_and_error_jack_biased(jack_tot,r);
            fprintf(fcC," %g  %g \t",m[0],m[1]);
            free(m);
            m=mean_and_error_jack_biased(jack_tot,r1);
            fprintf(fcC,"%g %g  %g \t",x[jack_tot-1][10],m[0],m[1]);
            free(m);
        }
        fprintf(fcC,"\n");

    }
    fclose(fc);
   

    free(r);
    free(r1);
    free_2(jack_tot,tif);
    free_2(jack_tot,x);
    fclose(fcA);    fclose(fcB);    fclose(fcC);
    
}


static void  print_fit_info(char **argv,int jack_tot,struct fit_result fit_out, struct fit_type fit_info, double **phys_point, struct data_jack *grephJ, struct header *head , const char *AV,const char *namefile){
    int i,j,k;
    
    double **fit=fit_out.P,*chi2m;
    
    double **tif,*tmp,**tmp2,*fk;
    double **Ci=(double**) malloc(sizeof(double*)*fit_info.Npar);
    char nametex[NAMESIZE],namegp[NAMESIZE];
    mysprintf(nametex,NAMESIZE,"%s/%s.tex",argv[2],namefile);
    mysprintf(namegp,NAMESIZE,"%s/%s.gp",argv[2],namefile);
    FILE *ftex=NULL, *fgp=NULL;
    ftex=open_file(nametex,"w+");
    fgp=open_file(namegp,"w+");
    const char *s;
    tmp=(double*)  malloc(sizeof(double) *jack_tot);
    tmp2=(double**) malloc(sizeof(double*)*fit_info.Npar);
    for (i=0;i<fit_info.Npar;i++)
        tmp2[i]=(double*)  malloc(sizeof(double) *jack_tot);

    tif=swap_indices(fit_info.Npar,jack_tot,fit);
    
    printf("chi2=%g\n",fit_out.chi2[jack_tot-1] );
    chi2m=mean_and_error_jack_biased(jack_tot,fit_out.chi2);
    for (i=0;i<fit_info.Npar;i++){
        Ci[i]=mean_and_error_jack_biased(jack_tot,fit[i]);
    }
    fprintf(ftex,"\\begin{align}\n");
    fprintf(ftex,"& \\chi^2/d.o.f.= %+.5f \\pm \t%.2g \\\\ \n",chi2m[0],chi2m[1]);
    for (i=0;i<fit_info.Npar;i++){
         if (strcmp(namefile,"fit_Mpi_Fpi")==0 ||  strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M1")==0  ||  strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M2a")==0 ||  strcmp(namefile,"fit_Mpi_Fpi_GL_w0_M2b")==0){
             if(i==0)     fprintf(ftex,"& Bw_{0}= %+.5f \\pm \t%.2g   \\\\ \n",Ci[i][0],Ci[i][1]);
             if(i==1)     fprintf(ftex,"& fw_{0}= %+.5f \\pm \t%.2g   \\\\ \n",Ci[i][0],Ci[i][1]);
             if(i==2)     fprintf(ftex,"& \\bar{\\ell_3}= %+.5f \\pm \t%.2g   \\\\ \n",Ci[i][0],Ci[i][1]);
             if(i==3)     fprintf(ftex,"& P_2= %+.5f \\pm \t%.2g   \\\\ \n",Ci[i][0],Ci[i][1]);
             if(i==4)     fprintf(ftex,"& \\bar{\\ell_4}= %+.5f \\pm \t%.2g   \\\\ \n",Ci[i][0],Ci[i][1]);
             if(i==5)     fprintf(ftex,"& P_4= %+.5f \\pm \t%.2g   \\\\ \n",Ci[i][0],Ci[i][1]);


        }
        else
                fprintf(ftex,"& P_{%d}= %+.5f \\pm \t%.2g   \\\\ \n",i,Ci[i][0],Ci[i][1]);
        
       
    }
    fprintf(ftex,"\\end{align}\n");
    
    for (i=0;i<fit_info.Npar;i++){
        fprintf(fgp,"P%d= %+.5f;\t",i,Ci[i][0]);
    }
    fprintf(fgp,"\n");
    for (i=0;i<fit_info.Npar;i++){       
        free(Ci[i]);
    }
    /*
        fprintf(ftex,"{\\tiny\\begin{gather}\n C=\\begin{pmatrix}\n");
    for(j=0;j<jack_tot;j++){
        for (i=0;i<fit_info.Npar;i++){
            for (k=0;k<i;k++)
                    jack_tot,fit_out.C[i][k][j]/=sqrt(fit_out.C[i][i][j]*fit_out.C[k][k][j]);
            for (k=i+1;k<fit_info.Npar;k++)
                    jack_tot,fit_out.C[i][k][j]/=sqrt(fit_out.C[i][i][j]*fit_out.C[k][k][j]);
        }
    }

     for (i=0;i<fit_info.Npar;i++){
        for (j=0;j<fit_info.Npar;j++){
            
            s=smean_and_error("jack",jack_tot,fit_out.C[i][j]);
            if (j==0)  fprintf(ftex,"%s",  s );
            else       fprintf(ftex,"& %s", s );
            free((void*)s);
        }
        if (i!=fit_info.Npar) fprintf(ftex,"\\\\ \n");
        else fprintf(ftex,"\n");
    }
    fprintf(ftex,"\\end{pmatrix}\n\\end{gather}}\n");
    */
    print_chiral_continuum_fit(argv, jack_tot,  fit_out,   fit_info,  phys_point,  AV, namefile,  head , grephJ);
   
    for (i=0;i<fit_info.Npar;i++){
        free(fit_out.P[i]);
        for (j=0;j<fit_info.Npar;j++){
            free(fit_out.C[i][j]);
        }
        free(fit_out.C[i]);
        free(tmp2[i]);
    }     
    free(tmp2);
    free(fit_out.chi2);
    free(fit_out.P);
    free(fit_out.C);
    
    free_2(jack_tot,tif);
    //free_tif(fit_info.Npar,fit);
    fclose(ftex);fclose(fgp);free(chi2m);free(tmp); free(Ci);
    
}

static void init_fit( int N, struct header *head ,int Njack, struct data_jack *gJ,int Nvar,int Npar, int **en,int *en_tot, double ****x, double ***sigmax, double **chi2m, double **rm, double ***r, double ***fit, double ****y,double **chi2,double ****C, int ensembles_diff=ensembles)
{
    int imoms,imomt,imom0,iG,i,n,e,j;
    int count;
   *en_tot=0;
   
   *en=(int*) calloc(N,sizeof(int));
   
   for (e=0;e<ensembles_diff;e++){
                for (n=0;n<N;n++){
                   
                    (*en)[n]+=1;  
                }
   }
      

 
   for (n=0;n<N;n++)
   {  *en_tot+=(*en)[n];   }
   
   *x=double_malloc_3(Njack,*en_tot,Nvar);//(double**) malloc(sizeof(double*)*(*en_tot));
   *sigmax=double_malloc_2(*en_tot,Nvar);

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
static struct fit_result close_fit( int N, struct header *head ,int Njack, struct data_jack *gJ,int Npar,int **en,int *en_tot, double ****x,  double ***sigmax, double **chi2m, double **rm,double ***r, double ***fit, double ****y,double **chi2, double ****C)
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
   //free(*x);
   free_3(Njack,*en_tot,*x);
   free_2(*en_tot,*sigmax);
   
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
               // fit_out.C[i][n][j]=(*C)[j][i][n];
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
 










 
 double **fit_MK_Mpi_fK_fpi_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info , struct result_jack *r1, const char *prefix,char **argv, store_fit_clover mud,  std::vector<int> myen){
     
     int ensembles_here=myen.size();  
     double ***y,***x,**sigmax,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit;   double **out;
     int i,j,im;  
     int Npar;
     int Nvar=1;//m_l, w0,M_PS^2,f_PS
     int ik1=0,ik2=1;
     int ik2_min=1, ik2_max=3;
     int nk=(ik2_max-ik2_min+1);
     int ms;
     double *tmp_chi=(double*) calloc(Njack,sizeof(double));
     
     double *mref;//[Nms]={0.52,0.68,0.81};
     mref=(double*) malloc(sizeof(double)*nk);
     mref[0]=0.064;
     mref[1]=0.080;
     mref[2]=0.095;
     int n,count,N=fit_info.N;
     int *en=(int*) malloc(sizeof(int)*N);
     
     for (i=0;i<N;i++)
         en[i]=nk;
     
     if (N==1)
         Npar=2;
     else if (N==2)
         Npar=4;
     
     
     
     int en_tot=0;
     
     for (n=0;n<N;n++)
         en_tot+=en[n];
     
     double *guess=(double*) malloc(sizeof(double)*Npar);
     for(i=0;i<Npar;i++)
         guess[i]=1.;
     
     
     //x=(double**) malloc(sizeof(double*)*(en_tot));
     
     chi2m=(double*) malloc(sizeof(double)*(Npar));
     rm=(double*) malloc(sizeof(double)*(Njack));
     fit=(double**) malloc(sizeof(double*)*(en_tot));
     
     r=(double***) malloc(sizeof(double**)*(Npar));
     for(i=0;i<Npar;i++){
         r[i]=(double**) malloc(sizeof(double*)*ensembles_here);
         for(j=0;j<ensembles_here;j++){
             r[i][j]=(double*) malloc(sizeof(double)*Njack);
         }
     }
     
     chi2=(double*) malloc(sizeof(double)*Njack);
     y=(double***) malloc(sizeof(double**)*Njack);
     y=double_malloc_3(Njack,en_tot,2);
     
     
     count=0;
     for (int e :myen){     
         for (n=0;n<N;n++){
             for (ms=0;ms<nk;ms++){
                 im=mass_index[e][ms+ik2_min][ik1];
                 int imp=mass_index[e][0][0];
                 if (n==0){
                     for (j=0;j<Njack;j++){
                         rm[j]=gJ[e].M_PS_jack[im][j]   /gJ[e].M_PS_jack[0][j]  ;
                         rm[j]*=rm[j];
                     }
                     fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                 }
                 if (n==1){
                     for (j=0;j<Njack;j++){
                         rm[j]=gJ[e].f_PS_jack[im][j]  *  gJ[e].w0[j];                           
                     }
                     fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                 }
                 
                 for (j=0;j<Njack;j++){
                     y[j][ms+n*nk][0]=rm[j];
                     y[j][ms+n*nk][1]=fit[ms+n*nk][1];
                     
                 }
                 
                 //x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                 
                 
             }
             
         } 
         
         x=double_malloc_3(Njack,en_tot,Nvar);
         sigmax=double_malloc_2(en_tot,Nvar);
         
         
         for (j=0;j<Njack;j++){
             for (n=0;n<N;n++){
                 for (ms=0;ms<nk;ms++){
                     im=mass_index[e][ms+ik2_min][ik1];
                     x[j][ms+n*nk][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                 }
             }
         }
         
         for (n=0;n<N;n++){
             for (ms=0;ms<nk;ms++){
                 for(int v=0 ;v<Nvar;v++){
                     for (j=0;j<Njack;j++)
                         rm[j]=x[j][ms+n*nk][v];
                     tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                     if (fabs(tmp[1])<1e-6) {
                         //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] );
                         tmp[1]=tmp[0]/1.0e+8;
                         
                     }
                     sigmax[ms+n*nk][v]=tmp[1];
                     free(tmp);
                 }
             }
         }
         
         
         
         for (j=0;j<Njack;j++){
             tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_lines,guess );
             //tmp=non_linear_fit_Nf_sigmax(N, en,x[j], sigmax, y[j] , Nvar,  Npar, two_lines,guess );
             //tmp=non_linear_fit_Nf_sigmax_iterative(N, en,x[j], sigmax, y[j] , Nvar,  Npar, two_lines,guess );
             chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, two_lines  );
             
             
             for(i=0;i<Npar;i++){
                 r[i][count][j]=tmp[i];
             }                
             free(tmp);
             
         }     
         for (ms=0;ms<en_tot;ms++){
             free(fit[ms]);
         }
         count++; 
     } 
     free_3(Njack,en_tot,x);
     free_2(en_tot,sigmax);
     
     
     free(fit);     
     
     free_3(Njack,en_tot,y);
     free(guess);
     im=mass_index[0][1][0];
     //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
     for (int e:myen)
     {im=mass_index[e][1+0][0];
         //    printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im]//[Njack-1]* gJ[e].w0[Njack-1],2) ,r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1],r[2][e][Njack-1]+mref[0]*r[3][e][Njack-1] );
         
     }
     free(chi2);
     ///////////////////////////////////////////////////////////compute MK at physical point
     free(en);
     Nvar=fit_info.Nvar;
     Npar=fit_info.Npar;
     guess=(double*) malloc(sizeof(double*)*Npar);
     for(i=0;i<Npar;i++)
         guess[i]=1.;
     
     
     double *rm1=(double*) malloc(sizeof(double)*(Njack));
     double **y1=(double**) malloc(sizeof(double*)*(Njack));
     // fit=(double**) malloc(sizeof(double*)*(en_tot));
     
     MK=(double**) malloc(sizeof(double*)*(nk*N));
     for(i=0;i<nk*N;i++){
         MK[i]=(double*) malloc(sizeof(double)*Njack);
     }
     
     double **xphys=double_malloc_2(Njack,Nvar);
     
     
     double KM,Kf;
     for (ms=0;ms<nk;ms++){
         
         char fname[NAMESIZE];
         mysprintf(fname,NAMESIZE,"%s/%s_ms%d",argv[2],prefix,ms);
         printf("writing data to :%s\n",fname);
         FILE *fdat=open_file(fname,"w+");
         double ***C,**tmp1;
         init_fit(  N, head , Njack, gJ,Nvar, Npar, &en, &en_tot, &x, &sigmax, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C, ensembles_here);
         
         count=0;
         for (n=0;n<N;n++){
             for (int e=0 ; e< myen.size();  e++){
                 im=mass_index[myen[e]][ik2][ik1];
                 if(n==0){
                     for (j=0;j<Njack;j++){
                         rm[j]=(r[0][e][j]+mref[ms]*r[1][e][j]);//(KM*KM);
                         rm1[j]=(r[0][e][j]+mref[ms]*r[1][e][j]);//(KM*KM);
                     }
                 }
                 if(n==1){
                     for (j=0;j<Njack;j++){
                         rm[j]=r[2][e][j]+mref[ms]*r[3][e][j];
                         //rm1[j]=rm[j]/Kf;
                     }
                     
                 }
                 fit[count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                 y1[count]=mean_and_error(jack_files[0].sampling,Njack, rm1);
                 
                 for (j=0;j<Njack;j++){
                     y[j][count][0]=rm[j];
                     y[j][count][1]=fit[count][1];
                     //if (e==3) {y[j][e+count][1]*=100.; }
                 }
                 count++;
                 
                 //x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                 
             }
             
         }   
         int ii;
         //#pragma omp parallel for  private(tmp,e,i,xphys,n,count)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)
         for (j=0;j<Njack;j++){
             count=0;
             for (n=0;n<N;n++){
                 for (int e=0 ; e< myen.size();  e++){
                     
                     
                     im=mass_index[myen[e]][ik2][ik1];
                     x[j][count][0]=head[myen[e]].k[head[myen[e]].nk+ik1]*gJ[myen[e]].w0[j]/gJ[myen[e]].Zp[j];//ml*w0
                     x[j][count][1]=gJ[myen[e]].w0[j];//w0
                     x[j][count][2]=gJ[myen[e]].M_PS_jack[im][j]*gJ[myen[e]].M_PS_jack[im][j];//MPS^2
                     x[j][count][3]=gJ[myen[e]].f_PS_jack[im][j];//f_PS
                     x[j][count][4]=double(head[myen[e]].l1);//f_PS
                     x[j][count][5]=mref[ms];//ms*w0
                     x[j][count][6]=r[0][e][j]+mref[ms]*r[1][e][j];//MKw2
                     if (N>1)
                         x[j][count][7]=r[2][e][j]+mref[ms]*r[3][e][j];//fkw
                         x[j][count][8]=mud.jack_B[j]*mud.w0[j]/197.326963;
                     x[j][count][9]=mud.jack_f[j]*mud.w0[j]/197.326963;
                     count++;
                 }
                 
             }
         }
         /*        count=0;
          *      for (n=0;n<1;n++){
          *          for (int e :myen){
          *              for(int v=0 ;v<Nvar;v++){
          *                  for (j=0;j<Njack;j++)
          *                      rm[j]=x[j][count][v];
          *                  tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
          *                  if (fabs(tmp[0])<1e-6) {
          *                      tmp[1]=tmp[0]/1.0e+8; }
          *                      sigmax[count][v]=tmp[1];
          *                      free(tmp);
     }
     count++;
     }
     count+=en[n];
     }
     */      
         guess=guess_for_non_linear_fit_Nf(N, en,x[0], y[0] , Nvar,  Npar, fit_info.function,guess ); 
         for (j=0;j<Njack;j++){
             
             
             //if (j==0){    }
             tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess );
             //tmp=non_linear_fit_Nf_sigmax(N, en,x[j], sigmax, y[j] , Nvar,  Npar, fit_info.function,guess );
             //tmp=non_linear_fit_Nf_sigmax_iterative(N, en,x[j], sigmax, y[j] , Nvar,  Npar, fit_info.function,guess );
             chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
             tmp_chi[j]+=chi2[j];
             C[j]=covariance_non_linear_fit_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function );  
             
             if(j==Njack-1){
                 //printf("#P0_w=%f ;P1_w=%f ;  P3ww=%f;  Pf1w=%f;  Pf2w=%f; Pf4www=%f;  msw=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],mref[ms]  );
                 printf("#chi2=%f    N=%d\n",chi2[j]/(en_tot-Npar),N);
                 
                 printf("# w0    mlw    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
                 fprintf(fdat,"# w0    mlw    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
                 double newline_if_w0=x[j][0][1];
                 
                 for (int e=0 ; e<myen.size();  e++){
                     im=mass_index[myen[e]][ik2][ik1];
                     KM=fit_info.function(2,Nvar,x[j][e],Npar,tmp);
                     Kf=fit_info.function(3,Nvar,x[j][e],Npar,tmp);
                     
                     if (newline_if_w0!=x[j][e][1]){
                         fprintf(fdat,"\n\n");
                         newline_if_w0=x[j][e][1] ;
                     }
                     if(N==1){
                         fprintf(fdat,"%f   %f    %f   %f           %g  %g\n",x[j][e][1],x[j][e][0], y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),KM,Kf);
                         printf("%f   %f    %f   %f           %g  %g\n",x[j][e][1],x[j][e][0], y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),KM,Kf);
                         
                     }
                     else if (N==2){
                         fprintf(fdat,"%f   %f    %f   %f     %f       %f       %g  %g\n",x[j][e][1],x[j][e][0],y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),y[j][e+en[0]][0]/(Kf),y[j][e+en[0]][1]/(Kf),KM,Kf);
                         printf("%f   %f    %f   %f     %f       %f       %g  %g\n",x[j][e][1],x[j][e][0],y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),y[j][e+en[0]][0]/(Kf),y[j][e+en[0]][1]/(Kf),KM,Kf);
                     }
                 }
             }
             
             
             xphys[j][0]=mud.jack_m[j] * mud.w0[j]/197.326963;
             xphys[j][1]=1e+6;  //w0
             xphys[j][2]=r1->MpiMeV[j]*r1->MpiMeV[j]*(mud.w0[j]/197.326963)*(mud.w0[j]/197.326963);
             xphys[j][3]=r1->fpiw[j];
             xphys[j][4]=1e+10;  // L such that L/w0=1e+6
             xphys[j][5]=mref[ms];
             xphys[j][6]=r1->MKMeV[j]*(mud.w0[j]/197.326963)*r1->MKMeV[j]*(mud.w0[j]/197.326963);//MKw2
             xphys[j][7]=0;//fkw
             xphys[j][8]=mud.jack_B[j]*mud.w0[j]/197.326963;
             xphys[j][9]=mud.jack_f[j]*mud.w0[j]/197.326963;
             
             /*if(j==Njack-1){
              *         for (i=0;i<10;i++)
              *            printf("xphis[%d]=%g\n",i,xphys[j][i]);
              *            printf("delta=%g  %g\n",fit_info.function(2,  Nvar, xphys[j],Npar,tmp),fit_info.function(3,  Nvar, xphys[j],Npar,tmp));       
         }*/
             //add w0 -inft
             MK[ms][j]=fit_info.function(0,  Nvar, xphys[j],Npar,tmp);//MK2
             //     printf("MKw02=%f\n",MK[ms][j]);
             if (N>1)
                 MK[ms+1*nk][j]=fit_info.function(1,  Nvar, xphys[j],Npar,tmp);//fK
                 
            for (i=0;i<Npar;i++){
                 tmp1[i][j]=tmp[i];
            }
                 
                 
            free(tmp);
         } 
         struct fit_result  fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &sigmax, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C);
         
         
         char namefile[NAMESIZE];
         mysprintf(namefile,NAMESIZE,"%s_ms%d",prefix,ms );
         print_fit_info(argv, Njack, fit_out,  fit_info, xphys, gJ, head , "K",namefile);
         
         
         for (int e=0;e<myen.size();e++){
             /*free(x[e]);  free(fit[e]);*/   //moved to close fit;
             free(y1[e]);
         }
         fclose(fdat);
     }//end loop ms
     free(y1);free(rm1);
     
     printf("MKw2(ms1)=%f    MKw2(ms2)=%f     MKw2(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
     if (N>1)
         printf("fKw(ms1)=%f     fKw(ms2)=%f      fKw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);
     free(guess);
     
     ////////////////////////////////////////////////last interpolation  
     en=(int*) malloc(sizeof(int)*N);
     
     en_tot=0;
     for (n=0;n<N;n++){
         en[n]=nk;
         en_tot+=en[n];
     }
     if (N==1)
         Npar=2;
     else if (N==2)    
         Npar=4;
     Nvar=1;//m_l, w0,M_PS^2,f_PS
     
     guess=(double*) malloc(sizeof(double*)*Npar);
     for(i=0;i<Npar;i++)
         guess[i]=1.;
     
     //x=(double**) malloc(sizeof(double*)*(en_tot));
     rm=(double*) malloc(sizeof(double)*Njack);
     
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
     
     y=double_malloc_3(Njack,en_tot,2);
     
     x=double_malloc_3(Njack,en_tot,Nvar);
     sigmax=double_malloc_2(en_tot,Nvar);
     
     //the last one is the chi2
     out=double_malloc_2(N+1,Njack);
     
     for (n=0;n<N;n++){
         for (ms=0;ms<nk;ms++){
             
             if (n==0){
                 for (j=0;j<Njack;j++){
                     rm[j]= MK[ms][j];
                     //  rm[j]*=rm[j];
                 }
                 fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
             }
             if (n==1){
                 for (j=0;j<Njack;j++){
                     rm[j]= MK[ms+1*nk][j];
                 }
                 fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
             }
             
             for (j=0;j<Njack;j++){
                 y[j][ms+n*nk][0]=rm[j];
                 y[j][ms+n*nk][1]=fit[ms+n*nk][1];
             }
             
             //x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
             //x[ms+n*nk][0]=mref[ms];//ml*w0
             
         }
         
     } 
     for (j=0;j<Njack;j++){
         for (n=0;n<N;n++){
             for (ms=0;ms<nk;ms++){
                 x[j][ms+n*nk][0]=mref[ms];//ml*w0
             }
         }
     }
     
     for (n=0;n<N;n++){
         for (ms=0;ms<nk;ms++){
             for(int v=0 ;v<Nvar;v++){
                 for (j=0;j<Njack;j++)
                     rm[j]=x[j][ms+n*nk][v];
                 tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                 if (fabs(tmp[1])<1e-6) {
                     //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] );
                     tmp[1]=tmp[0]/1.0e+8; }
                     sigmax[ms+n*nk][v]=tmp[1];
                     free(tmp);
             }
         }
     }
     
     
     double in;
     for (j=0;j<Njack;j++){
         tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_lines,guess );
         
         //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
         //in=r1->MKMeV[j]*(mud.w0[j]/197.326963)*r1->MKMeV[j]*mud.w0[j]/197.326963;
         in=r1->MKMeV[j]/r1->MpiMeV[j];
         in=in*in;
         out[0][j]=(in-tmp[0])/tmp[1];
         //printf("%g  %g  %g  %g  \n",tmp[0],tmp[1],out[0][j],in);
         if (N>1)
             out[1][j]=tmp[2]+tmp[3]*out[0][j];
         
         out[N][j]=tmp_chi[j]/nk;//average chi2 of the three fit
         
         free(tmp);
         
     }     
     for (ms=0;ms<en_tot;ms++){
         free(fit[ms]);
     }
     
     free_3(Njack,en_tot,x);
     free_2(en_tot,sigmax);
     
     
     free(fit);     
     
     free_3(Njack,en_tot,y);
     free(guess);
     free(tmp_chi);
     
     
     
     free(rm);   
     return out;
     
} 

void interpolate_at_fixed_a(fit_type fit_info, int nk, int Npar, int Nvar, int Njack, double **MK, const char *resampling, double *mref, struct result_jack *r1, store_fit_clover mud, double *w0, double *Zp){
    
    ////////////////////////////////////////////////last interpolation 
    int N=fit_info.N;
    
    int *en=(int*) malloc(sizeof(int)*N);
    
    
    int en_tot=0;
    for (int n=0;n<N;n++){
        en[n]=nk;
        en_tot+=en[n];
    }
    
    
    double *guess=(double*) malloc(sizeof(double*)*Npar);
    for(int i=0;i<Npar;i++)
        guess[i]=1.;
   
    //x=(double**) malloc(sizeof(double*)*(en_tot));
    double *rm=(double*) malloc(sizeof(double)*Njack);
    
    double **fit=(double**) malloc(sizeof(double*)*(en_tot));
    
    double ***y=(double***) malloc(sizeof(double**)*Njack);
    for (int j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (int n=0;n<N;n++){
            for (int ms=0;ms<nk;ms++){
                y[j][ms+n*nk]=(double*) malloc(sizeof(double)*2);
                
            }
        }
    }
    
    y=double_malloc_3(Njack,en_tot,2);
    
    double ***x=double_malloc_3(Njack,en_tot,Nvar);
    double **sigmax=double_malloc_2(en_tot,Nvar);
    
    //the last one is the chi2
    double **out=double_malloc_2(N+1,Njack);
    
    for (int n=0;n<N;n++){
        for (int ms=0;ms<nk;ms++){
            
            if (n==0){
                for (int j=0;j<Njack;j++){
                    rm[j]= MK[ms][j];
                    //  rm[j]*=rm[j];
                }
                fit[ms]=mean_and_error(resampling,Njack, rm);
            }
            if (n==1){
                for (int j=0;j<Njack;j++){
                    rm[j]= MK[ms+1*nk][j];
                }
                fit[ms+n*nk]=mean_and_error(resampling,Njack, rm);
            }
            
            for (int j=0;j<Njack;j++){
                y[j][ms+n*nk][0]=rm[j];
                y[j][ms+n*nk][1]=fit[ms+n*nk][1];
            }
            
            //x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
            //x[ms+n*nk][0]=mref[ms];//ml*w0
            
        }
        
    } 
    for (int j=0;j<Njack;j++){
        for (int n=0;n<N;n++){
            for (int ms=0;ms<nk;ms++){
                x[j][ms+n*nk][0]=mref[ms];//ml*w0
            }
        }
    }
    
    for (int n=0;n<N;n++){
        for (int ms=0;ms<nk;ms++){
            for(int v=0 ;v<Nvar;v++){
                for (int j=0;j<Njack;j++)
                    rm[j]=x[j][ms+n*nk][v];
                double *tmp=mean_and_error(resampling,Njack, rm);
                if (fabs(tmp[1])<1e-6) {
                    //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] );
                    tmp[1]=tmp[0]/1.0e+8;
                    
                }
                sigmax[ms+n*nk][v]=tmp[1];
                free(tmp);
            }
        }
    }
    
    
    double in;
    for (int j=0;j<Njack;j++){
        double *tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_lines,guess );
        
        //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        //in=r1->MKMeV[j]*(mud.w0[j]/197.326963)*r1->MKMeV[j]*mud.w0[j]/197.326963;
        in=r1->MKMeV[j]*(mud.w0[j]/197.326963);
        in=in*in;
        out[0][j]=((in-tmp[0])/tmp[1])*Zp[j]/w0[j];
        //printf("%g  %g  %g  %g  \n",tmp[0],tmp[1],out[0][j],in);
        if (N>1)
            out[1][j]=(tmp[2]+tmp[3]*out[0][j]);
        
        
        free(tmp);
        
    }  
    printf("ms at fixed a=%g   +-  %g\n",out[0][Njack-1],error_jackboot(resampling,Njack,out[0]) );
    
    for (int ms=0;ms<en_tot;ms++){
        free(fit[ms]);
    }
    
    free_3(Njack,en_tot,x);
    free_2(en_tot,sigmax);
    
    
    free(fit);     
                
    free_3(Njack,en_tot,y);
    free(guess);
    
}



double **fit_MK_fK_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info , struct result_jack *r1, const char *prefix,char **argv, store_fit_clover mud,  std::vector<int> myen){
    
    int ensembles_here=myen.size();  
    double ***y,***x,**sigmax,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit;   double **out;
    int i,j,im;  
    int Npar;
    int Nvar=1;//m_l, w0,M_PS^2,f_PS
    int ik1=0,ik2=1;
    int ik2_min=1, ik2_max=3;
    int nk=(ik2_max-ik2_min+1);
    int ms;
    double *tmp_chi=(double*) calloc(Njack,sizeof(double));
    
    double *mref;//[Nms]={0.52,0.68,0.81};
    mref=(double*) malloc(sizeof(double)*nk);
    mref[0]=0.064;
    mref[1]=0.080;
    mref[2]=0.095;
    int n,count,N=fit_info.N;
    int *en=(int*) malloc(sizeof(int)*N);
    
    for (i=0;i<N;i++)
        en[i]=nk;
    
    if (N==1)
        Npar=2;
    else if (N==2)
        Npar=4;
    
    
    
    int en_tot=0;
    
    for (n=0;n<N;n++)
        en_tot+=en[n];
    
    double *guess=(double*) malloc(sizeof(double)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
    
    
    //x=(double**) malloc(sizeof(double*)*(en_tot));
    
    chi2m=(double*) malloc(sizeof(double)*(Npar));
    rm=(double*) malloc(sizeof(double)*(Njack));
    fit=(double**) malloc(sizeof(double*)*(en_tot));
    
    r=(double***) malloc(sizeof(double**)*(Npar));
    for(i=0;i<Npar;i++){
        r[i]=(double**) malloc(sizeof(double*)*ensembles_here);
        for(j=0;j<ensembles_here;j++){
            r[i][j]=(double*) malloc(sizeof(double)*Njack);
        }
    }
    
    chi2=(double*) malloc(sizeof(double)*Njack);
    y=(double***) malloc(sizeof(double**)*Njack);
    y=double_malloc_3(Njack,en_tot,2);
    
    
    count=0;
    for (int e :myen){     
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                int imp=mass_index[e][0][0];
                if (n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j];
                        rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].f_PS_jack[im][j]  *  gJ[e].w0[j];                           
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<Njack;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms+n*nk][1];
                    
                }
                
                //x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                
                
            }
            
        } 
        
        x=double_malloc_3(Njack,en_tot,Nvar);
        sigmax=double_malloc_2(en_tot,Nvar);
        
        
        for (j=0;j<Njack;j++){
            for (n=0;n<N;n++){
                for (ms=0;ms<nk;ms++){
                    im=mass_index[e][ms+ik2_min][ik1];
                    x[j][ms+n*nk][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                }
            }
        }
        
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                for(int v=0 ;v<Nvar;v++){
                    for (j=0;j<Njack;j++)
                        rm[j]=x[j][ms+n*nk][v];
                    tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                    if (fabs(tmp[1])<1e-6) {
                        //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] );
                        tmp[1]=tmp[0]/1.0e+8;
                        
                    }
                    sigmax[ms+n*nk][v]=tmp[1];
                    free(tmp);
                }
            }
        }
        
        
        
        for (j=0;j<Njack;j++){
            tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_lines,guess );
            //tmp=non_linear_fit_Nf_sigmax(N, en,x[j], sigmax, y[j] , Nvar,  Npar, two_lines,guess );
            //tmp=non_linear_fit_Nf_sigmax_iterative(N, en,x[j], sigmax, y[j] , Nvar,  Npar, two_lines,guess );
            chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, two_lines  );
            
            
            for(i=0;i<Npar;i++){
                r[i][count][j]=tmp[i];
            }                
            free(tmp);
            
        }     
        for (ms=0;ms<en_tot;ms++){
            free(fit[ms]);
        }
        count++; 
    } 
    free_3(Njack,en_tot,x);
    free_2(en_tot,sigmax);
    
    
    free(fit);     
    
    free_3(Njack,en_tot,y);
    free(guess);
    im=mass_index[0][1][0];
    //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
    for (int e:myen)
    {im=mass_index[e][1+0][0];
    //    printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im]//[Njack-1]* gJ[e].w0[Njack-1],2) ,r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1],r[2][e][Njack-1]+mref[0]*r[3][e][Njack-1] );
        
    }
    free(chi2);
    ///////////////////////////////////////////////////////////compute MK at physical point
    free(en);
    Nvar=fit_info.Nvar;
    Npar=fit_info.Npar;
    guess=(double*) malloc(sizeof(double*)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
    
    
    double *rm1=(double*) malloc(sizeof(double)*(Njack));
    double **y1=(double**) malloc(sizeof(double*)*(Njack));
    // fit=(double**) malloc(sizeof(double*)*(en_tot));
    
    MK=(double**) malloc(sizeof(double*)*(nk*N));
    for(i=0;i<nk*N;i++){
        MK[i]=(double*) malloc(sizeof(double)*Njack);
    }
    double ***MK_a=double_malloc_3(3,nk*N, Njack);
    
    double **xphys=double_malloc_2(Njack,Nvar);
    
    
    double KM,Kf;
    for (ms=0;ms<nk;ms++){
        
        char fname[NAMESIZE];
        mysprintf(fname,NAMESIZE,"%s/%s_ms%d",argv[2],prefix,ms);
        printf("writing data to :%s\n",fname);
        FILE *fdat=open_file(fname,"w+");
        double ***C,**tmp1;
        init_fit(  N, head , Njack, gJ,Nvar, Npar, &en, &en_tot, &x, &sigmax, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C, ensembles_here);
        
        count=0;
        for (n=0;n<N;n++){
            for (int e=0 ; e< myen.size();  e++){
                im=mass_index[myen[e]][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=(r[0][e][j]+mref[ms]*r[1][e][j]);//(KM*KM);
                        rm1[j]=(r[0][e][j]+mref[ms]*r[1][e][j]);//(KM*KM);
                    }
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[2][e][j]+mref[ms]*r[3][e][j];
                        //rm1[j]=rm[j]/Kf;
                    }
                    
                }
                fit[count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                y1[count]=mean_and_error(jack_files[0].sampling,Njack, rm1);
                
                for (j=0;j<Njack;j++){
                    y[j][count][0]=rm[j];
                    y[j][count][1]=fit[count][1];
                    //if (e==3) {y[j][e+count][1]*=100.; }
                }
                count++;
                
                //x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                
            }
            
        }   
        int ii;
        //#pragma omp parallel for  private(tmp,e,i,xphys,n,count)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)
        for (j=0;j<Njack;j++){
            count=0;
            for (n=0;n<N;n++){
                for (int e=0 ; e< myen.size();  e++){
                    
                    
                    im=mass_index[myen[e]][ik2][ik1];
                    x[j][count][0]=head[myen[e]].k[head[myen[e]].nk+ik1]*gJ[myen[e]].w0[j]/gJ[myen[e]].Zp[j];//ml*w0
                    x[j][count][1]=gJ[myen[e]].w0[j];//w0
                    x[j][count][2]=gJ[myen[e]].M_PS_jack[im][j]*gJ[myen[e]].M_PS_jack[im][j];//MPS^2
                    x[j][count][3]=gJ[myen[e]].f_PS_jack[im][j];//f_PS
                    x[j][count][4]=double(head[myen[e]].l1);//f_PS
                    x[j][count][5]=mref[ms];//ms*w0
                    x[j][count][6]=r[0][e][j]+mref[ms]*r[1][e][j];//MKw2
                    if (N>1)
                        x[j][count][7]=r[2][e][j]+mref[ms]*r[3][e][j];//fkw
                    x[j][count][8]=mud.jack_B[j]*mud.w0[j]/197.326963;
                    x[j][count][9]=mud.jack_f[j]*mud.w0[j]/197.326963;
                    count++;
                }
                
            }
        }
/*        count=0;
        for (n=0;n<1;n++){
            for (int e :myen){
                for(int v=0 ;v<Nvar;v++){
                    for (j=0;j<Njack;j++)
                        rm[j]=x[j][count][v];
                    tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                    if (fabs(tmp[0])<1e-6) {
                        tmp[1]=tmp[0]/1.0e+8; }
                        sigmax[count][v]=tmp[1];
                        free(tmp);
                }
                count++;
            }
            count+=en[n];
        }
  */      
        guess=guess_for_non_linear_fit_Nf(N, en,x[0], y[0] , Nvar,  Npar, fit_info.function,guess ); 
        for (j=0;j<Njack;j++){
            
            
            //if (j==0){    }
            tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess );
            //tmp=non_linear_fit_Nf_sigmax(N, en,x[j], sigmax, y[j] , Nvar,  Npar, fit_info.function,guess );
            //tmp=non_linear_fit_Nf_sigmax_iterative(N, en,x[j], sigmax, y[j] , Nvar,  Npar, fit_info.function,guess );
            chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
            tmp_chi[j]+=chi2[j];
            C[j]=covariance_non_linear_fit_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function );  
            
            if(j==Njack-1){
                //printf("#P0_w=%f ;P1_w=%f ;  P3ww=%f;  Pf1w=%f;  Pf2w=%f; Pf4www=%f;  msw=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],mref[ms]  );
                printf("#chi2=%f    N=%d\n",chi2[j]/(en_tot-Npar),N);
                
                printf("# w0    mlw    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
                fprintf(fdat,"# w0    mlw    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
                double newline_if_w0=x[j][0][1];
               
                for (int e=0 ; e<myen.size();  e++){
                    im=mass_index[myen[e]][ik2][ik1];
                    KM=fit_info.function(2,Nvar,x[j][e],Npar,tmp);
                    Kf=fit_info.function(3,Nvar,x[j][e],Npar,tmp);
                    
                    if (newline_if_w0!=x[j][e][1]){
                        fprintf(fdat,"\n\n");
                        newline_if_w0=x[j][e][1] ;
                    }
                    if(N==1){
                        fprintf(fdat,"%f   %f    %f   %f           %g  %g\n",x[j][e][1],x[j][e][0], y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),KM,Kf);
                        printf("%f   %f    %f   %f           %g  %g\n",x[j][e][1],x[j][e][0], y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),KM,Kf);
                        
                    }
                    else if (N==2){
                        fprintf(fdat,"%f   %f    %f   %f     %f       %f       %g  %g\n",x[j][e][1],x[j][e][0],y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),y[j][e+en[0]][0]/(Kf),y[j][e+en[0]][1]/(Kf),KM,Kf);
                        printf("%f   %f    %f   %f     %f       %f       %g  %g\n",x[j][e][1],x[j][e][0],y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),y[j][e+en[0]][0]/(Kf),y[j][e+en[0]][1]/(Kf),KM,Kf);
                    }
                }
            }
            
            
            xphys[j][0]=mud.jack_m[j] * mud.w0[j]/197.326963;
            xphys[j][1]=1e+6;  //w0
            xphys[j][2]=r1->MpiMeV[j]*r1->MpiMeV[j]*(mud.w0[j]/197.326963)*(mud.w0[j]/197.326963);
            xphys[j][3]=r1->fpiw[j];
            xphys[j][4]=1e+10;  // L such that L/w0=1e+6
            xphys[j][5]=mref[ms];
            xphys[j][6]=r1->MKMeV[j]*(mud.w0[j]/197.326963)*r1->MKMeV[j]*(mud.w0[j]/197.326963);//MKw2
            xphys[j][7]=0;//fkw
            xphys[j][8]=mud.jack_B[j]*mud.w0[j]/197.326963;
            xphys[j][9]=mud.jack_f[j]*mud.w0[j]/197.326963;
            
            /*if(j==Njack-1){
             *         for (i=0;i<10;i++)
             *            printf("xphis[%d]=%g\n",i,xphys[j][i]);
             *            printf("delta=%g  %g\n",fit_info.function(2,  Nvar, xphys[j],Npar,tmp),fit_info.function(3,  Nvar, xphys[j],Npar,tmp));       
        }*/
            //add w0 -inft
            MK[ms][j]=fit_info.function(0,  Nvar, xphys[j],Npar,tmp);//MK2
            //     printf("MKw02=%f\n",MK[ms][j]);
            if (N>1)
                MK[ms+1*nk][j]=fit_info.function(1,  Nvar, xphys[j],Npar,tmp);//fK
            
            for (i=0;i<Npar;i++){
                tmp1[i][j]=tmp[i];
            }
            // interpolate MK at given lattice spacing to the physical value
            for (int a =0; a<3;a++){
                double *w0;
                if (a==0)       w0=gJ[0].w0;
                if (a==1)       w0=gJ[4].w0;
                if (a==2)       w0=gJ[7].w0;
                
                xphys[j][0]=mud.jack_m[j] *( mud.w0[j]/197.326963);//*w0[j];
                xphys[j][1]=w0[j];  //w0
                xphys[j][2]=r1->MpiMeV[j]*r1->MpiMeV[j]*(mud.w0[j]/197.326963)*(mud.w0[j]/197.326963);//*w0[j]*w0[j];
                xphys[j][3]=r1->fpiw[j];
                xphys[j][4]=1e+10;  // L such that L/w0=1e+6
                xphys[j][5]=mref[ms];
                xphys[j][6]=r1->MKMeV[j]*(mud.w0[j]/197.326963)*r1->MKMeV[j]*(mud.w0[j]/197.326963);//*w0[j]*w0[j];//MKw2
                xphys[j][7]=0;//fkw
                xphys[j][8]=mud.jack_B[j]*(mud.w0[j]/197.326963);//*w0[j];
                xphys[j][9]=mud.jack_f[j]*(mud.w0[j]/197.326963);//*w0[j];
                MK_a[a][ms][j]=fit_info.function(0,  Nvar, xphys[j],Npar,tmp);
            }
            free(tmp);
        } 
        struct fit_result  fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &sigmax, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C);
        
        
        char namefile[NAMESIZE];
        mysprintf(namefile,NAMESIZE,"%s_ms%d",prefix,ms );
        print_fit_info(argv, Njack, fit_out,  fit_info, xphys, gJ, head , "K",namefile);
        
        
        for (int e=0;e<myen.size();e++){
            /*free(x[e]);  free(fit[e]);*/   //moved to close fit;
            free(y1[e]);
        }
        fclose(fdat);
    }//end loop ms
    free(y1);free(rm1);
    
    printf("MKw2(ms1)=%f    MKw2(ms2)=%f     MKw2(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
    if (N>1)
        printf("fKw(ms1)=%f     fKw(ms2)=%f      fKw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);
    free(guess);
    
    ////////////////////////////////////////////////last interpolation  
    en=(int*) malloc(sizeof(int)*N);
    
    en_tot=0;
    for (n=0;n<N;n++){
        en[n]=nk;
        en_tot+=en[n];
    }
    if (N==1)
        Npar=2;
    else if (N==2)    
        Npar=4;
    Nvar=1;//m_l, w0,M_PS^2,f_PS
    
    interpolate_at_fixed_a( fit_info,  nk,  Npar,  Nvar,  Njack, MK_a[0], jack_files[0].sampling, mref, r1, mud, gJ[0].w0, gJ[0].Zp);
    interpolate_at_fixed_a( fit_info,  nk,  Npar,  Nvar,  Njack, MK_a[1], jack_files[0].sampling, mref, r1, mud, gJ[4].w0, gJ[4].Zp);
    interpolate_at_fixed_a( fit_info,  nk,  Npar,  Nvar,  Njack, MK_a[2], jack_files[0].sampling, mref, r1, mud, gJ[7].w0, gJ[7].Zp);
    
    guess=(double*) malloc(sizeof(double*)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
   
    //x=(double**) malloc(sizeof(double*)*(en_tot));
    rm=(double*) malloc(sizeof(double)*Njack);
    
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
    
    y=double_malloc_3(Njack,en_tot,2);
    
    x=double_malloc_3(Njack,en_tot,Nvar);
    sigmax=double_malloc_2(en_tot,Nvar);
    
    //the last one is the chi2
    out=double_malloc_2(N+1,Njack);
    
    for (n=0;n<N;n++){
        for (ms=0;ms<nk;ms++){
            
            if (n==0){
                for (j=0;j<Njack;j++){
                    rm[j]= MK[ms][j];
                    //  rm[j]*=rm[j];
                }
                fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
            }
            if (n==1){
                for (j=0;j<Njack;j++){
                    rm[j]= MK[ms+1*nk][j];
                }
                fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
            }
            
            for (j=0;j<Njack;j++){
                y[j][ms+n*nk][0]=rm[j];
                y[j][ms+n*nk][1]=fit[ms+n*nk][1];
            }
            
            //x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
            //x[ms+n*nk][0]=mref[ms];//ml*w0
            
        }
        
    } 
    for (j=0;j<Njack;j++){
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                x[j][ms+n*nk][0]=mref[ms];//ml*w0
            }
        }
    }
    
    for (n=0;n<N;n++){
        for (ms=0;ms<nk;ms++){
            for(int v=0 ;v<Nvar;v++){
                for (j=0;j<Njack;j++)
                    rm[j]=x[j][ms+n*nk][v];
                tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                if (fabs(tmp[1])<1e-6) {
                    //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] );
                    tmp[1]=tmp[0]/1.0e+8; }
                    sigmax[ms+n*nk][v]=tmp[1];
                    free(tmp);
            }
        }
    }
    
    
    double in;
    for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_lines,guess );
        
        //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        //in=r1->MKMeV[j]*(mud.w0[j]/197.326963)*r1->MKMeV[j]*mud.w0[j]/197.326963;
        in=r1->MKMeV[j]*mud.w0[j]/197.326963;
        in=in*in;
        out[0][j]=(in-tmp[0])/tmp[1];
        //printf("%g  %g  %g  %g  \n",tmp[0],tmp[1],out[0][j],in);
        if (N>1)
            out[1][j]=tmp[2]+tmp[3]*out[0][j];
        
        out[N][j]=tmp_chi[j]/nk;//average chi2 of the three fit
        
        free(tmp);
        
    }     
    for (ms=0;ms<en_tot;ms++){
        free(fit[ms]);
    }
    
    free_3(Njack,en_tot,x);
    free_2(en_tot,sigmax);
    
    
    free(fit);     
                
    free_3(Njack,en_tot,y);
    free(guess);
    free(tmp_chi);
    
    free_2(nk*N, MK);
    free_3(3, nk*N, MK_a);
    
    
    free(rm);   
    return out;
    
} 





double **fit_Mpi_MK_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info , struct result_jack *r1, const char *prefix,char **argv, store_fit_clover mud,  std::vector<int> myen){
    
    int ensembles_here=myen.size();  
    double ***y,***x,**sigmax,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit;   double **out;
    int i,j,im;  
    int Npar;
    int Nvar=1;//m_l, w0,M_PS^2,f_PS
    int ik1=0,ik2=1;
    int ik2_min=1, ik2_max=3;
    int nk=(ik2_max-ik2_min+1);
    int ms;
    double *tmp_chi=(double*) calloc(Njack,sizeof(double));
    
    double *mref;//[Nms]={0.52,0.68,0.81};
    mref=(double*) malloc(sizeof(double)*nk);
    mref[0]=0.064;
    mref[1]=0.080;
    mref[2]=0.095;
    int n,count,N=fit_info.N;
    int *en=(int*) malloc(sizeof(int)*N);
    
    for (i=0;i<N;i++)
        en[i]=nk;
    
    if (N==1)
        Npar=2;
    else if (N==2)
        Npar=4;
    
    
    
    int en_tot=0;
    
    for (n=0;n<N;n++)
        en_tot+=en[n];
    
    double *guess=(double*) malloc(sizeof(double)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
    
    
    //x=(double**) malloc(sizeof(double*)*(en_tot));
    
    chi2m=(double*) malloc(sizeof(double)*(Npar));
    rm=(double*) malloc(sizeof(double)*(Njack));
    fit=(double**) malloc(sizeof(double*)*(en_tot));
    
    r=(double***) malloc(sizeof(double**)*(Npar));
    for(i=0;i<Npar;i++){
        r[i]=(double**) malloc(sizeof(double*)*ensembles_here);
        for(j=0;j<ensembles_here;j++){
            r[i][j]=(double*) malloc(sizeof(double)*Njack);
        }
    }
    
    chi2=(double*) malloc(sizeof(double)*Njack);
    y=(double***) malloc(sizeof(double**)*Njack);
    y=double_malloc_3(Njack,en_tot,2);
    
    
    count=0;
    for (int e :myen){     
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                int imp=mass_index[e][0][0];
                if (n==0){
                    for (j=0;j<Njack;j++){
                        //rm[j]=gJ[e].M_PS_jack[imp][j]/gJ[e].M_PS_jack[im][j]   ;
                        //rm[j]*=rm[j];
                        rm[j]=gJ[e].M_PS_jack[im][j]   ;
                        rm[j]*=rm[j];
                        
                        
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[e].f_PS_jack[im][j]  *  gJ[e].w0[j];                           
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<Njack;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms+n*nk][1];
                    
                }
                
                //x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                
                
            }
            
        } 
        
        x=double_malloc_3(Njack,en_tot,Nvar);
        sigmax=double_malloc_2(en_tot,Nvar);
        
        
        for (j=0;j<Njack;j++){
            for (n=0;n<N;n++){
                for (ms=0;ms<nk;ms++){
                    im=mass_index[e][ms+ik2_min][ik1];
                    x[j][ms+n*nk][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                }
            }
        }
        
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                for(int v=0 ;v<Nvar;v++){
                    for (j=0;j<Njack;j++)
                        rm[j]=x[j][ms+n*nk][v];
                    tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                    if (fabs(tmp[1])<1e-6) {
                        //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] );
                        tmp[1]=tmp[0]/1.0e+8;
                        
                    }
                    sigmax[ms+n*nk][v]=tmp[1];
                    free(tmp);
                }
            }
        }
        
        
        
        for (j=0;j<Njack;j++){
            tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_lines,guess );
            //tmp=non_linear_fit_Nf_sigmax(N, en,x[j], sigmax, y[j] , Nvar,  Npar, two_lines,guess );
            //tmp=non_linear_fit_Nf_sigmax_iterative(N, en,x[j], sigmax, y[j] , Nvar,  Npar, two_lines,guess );
            chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, two_lines  );
            
            
            for(i=0;i<Npar;i++){
                r[i][count][j]=tmp[i];
            }                
            free(tmp);
            
        }     
        for (ms=0;ms<en_tot;ms++){
            free(fit[ms]);
        }
        count++; 
    } 
    free_3(Njack,en_tot,x);
    free_2(en_tot,sigmax);
    
    
    free(fit);     
    
    free_3(Njack,en_tot,y);
    free(guess);
    im=mass_index[0][1][0];
    //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
    for (int e:myen)
    {im=mass_index[e][1+0][0];
        //    printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im]//[Njack-1]* gJ[e].w0[Njack-1],2) ,r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1],r[2][e][Njack-1]+mref[0]*r[3][e][Njack-1] );
        
    }
    free(chi2);
    ///////////////////////////////////////////////////////////compute MK at physical point
    free(en);
    Nvar=fit_info.Nvar;
    Npar=fit_info.Npar;
    guess=(double*) malloc(sizeof(double*)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
    
    
    double *rm1=(double*) malloc(sizeof(double)*(Njack));
    double **y1=(double**) malloc(sizeof(double*)*(Njack));
    // fit=(double**) malloc(sizeof(double*)*(en_tot));
    
    MK=(double**) malloc(sizeof(double*)*(nk*N));
    for(i=0;i<nk*N;i++){
        MK[i]=(double*) malloc(sizeof(double)*Njack);
    }
    
    double **xphys=double_malloc_2(Njack,Nvar);
    
    
    double KM,Kf;
    for (ms=0;ms<nk;ms++){
        
        char fname[NAMESIZE];
        mysprintf(fname,NAMESIZE,"%s/%s_ms%d",argv[2],prefix,ms);
        printf("writing data to :%s\n",fname);
        FILE *fdat=open_file(fname,"w+");
        double ***C,**tmp1;
        init_fit(  N, head , Njack, gJ,Nvar, Npar, &en, &en_tot, &x, &sigmax, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C, ensembles_here);
        
        count=0;
        for (n=0;n<N;n++){
            for (int e=0 ; e< myen.size();  e++){
                im=mass_index[myen[e]][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[myen[e]].M_PS_jack[0][j]*gJ[myen[e]].M_PS_jack[0][j]/(r[0][e][j]+mref[ms]*r[1][e][j]);//(KM*KM);
                        rm1[j]=gJ[myen[e]].M_PS_jack[0][j]*gJ[myen[e]].M_PS_jack[0][j]/(r[0][e][j]+mref[ms]*r[1][e][j]);//(KM*KM);
                    }
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[2][e][j]+mref[ms]*r[3][e][j];
                        //rm1[j]=rm[j]/Kf;
                    }
                    
                }
                fit[count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                y1[count]=mean_and_error(jack_files[0].sampling,Njack, rm1);
                
                for (j=0;j<Njack;j++){
                    y[j][count][0]=rm[j];
                    y[j][count][1]=fit[count][1];
                    //if (e==3) {y[j][e+count][1]*=100.; }
                }
                count++;
                
                //x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                
            }
            
        }   
        int ii;
        //#pragma omp parallel for  private(tmp,e,i,xphys,n,count)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)
        for (j=0;j<Njack;j++){
            count=0;
            for (n=0;n<N;n++){
                for (int e=0 ; e< myen.size();  e++){
                    
                    
                    im=mass_index[myen[e]][ik2][ik1];
                    x[j][count][0]=head[myen[e]].k[head[myen[e]].nk+ik1]*gJ[myen[e]].w0[j]/gJ[myen[e]].Zp[j];//ml*w0
                    
                    
                    x[j][count][1]=gJ[myen[e]].w0[j];//w0
                    x[j][count][2]=gJ[myen[e]].M_PS_jack[0][j]*gJ[myen[e]].M_PS_jack[im][j];//MPS^2
                    x[j][count][3]=gJ[myen[e]].f_PS_jack[0][j];//f_PS
                    x[j][count][4]=double(head[myen[e]].l1);//f_PS
                    x[j][count][5]=mref[ms];//ms*w0
                    x[j][count][6]=r[0][e][j]+mref[ms]*r[1][e][j];//MKw2
                    if (N>1)
                        x[j][count][7]=r[2][e][j]+mref[ms]*r[3][e][j];//fkw
                    x[j][count][8]=mud.jack_B[j]*mud.w0[j]/197.326963;
                    x[j][count][9]=mud.jack_f[j]*mud.w0[j]/197.326963;
                    x[j][count][10]=gJ[myen[e]].M_PS_jack[0][j]*gJ[myen[e]].w0[j];
                    //x[j][count][11]=gJ[myen[e]].f_PS_jack[0][j]*gJ[myen[e]].w0[j];//f_PS
                    x[j][count][11]=124.0*gJ[0].scalefm[j]/197.326963;
                    
                    
                    count++;
                }
                
            }
        }
        /*        count=0;
         *      for (n=0;n<1;n++){
         *          for (int e :myen){
         *              for(int v=0 ;v<Nvar;v++){
         *                  for (j=0;j<Njack;j++)
         *                      rm[j]=x[j][count][v];
         *                  tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
         *                  if (fabs(tmp[0])<1e-6) {
         *                      tmp[1]=tmp[0]/1.0e+8; }
         *                      sigmax[count][v]=tmp[1];
         *                      free(tmp);
    }
    count++;
    }
    count+=en[n];
    }
    */      
        guess=guess_for_non_linear_fit_Nf(N, en,x[0], y[0] , Nvar,  Npar, fit_info.function,guess ); 
        for (j=0;j<Njack;j++){
            
            
            //if (j==0){    }
            tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess );
            //tmp=non_linear_fit_Nf_sigmax(N, en,x[j], sigmax, y[j] , Nvar,  Npar, fit_info.function,guess );
            //tmp=non_linear_fit_Nf_sigmax_iterative(N, en,x[j], sigmax, y[j] , Nvar,  Npar, fit_info.function,guess );
            chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
            tmp_chi[j]+=chi2[j];
            C[j]=covariance_non_linear_fit_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function );  
            
            if(j==Njack-1){
                //printf("#P0_w=%f ;P1_w=%f ;  P3ww=%f;  Pf1w=%f;  Pf2w=%f; Pf4www=%f;  msw=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],mref[ms]  );
                printf("#chi2=%f    N=%d\n",chi2[j]/(en_tot-Npar),N);
                
                printf("# w0    xi    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
                fprintf(fdat,"# w0    xi    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
                double newline_if_w0=x[j][0][1];
                
                for (int e=0 ; e<myen.size();  e++){
                    im=mass_index[myen[e]][ik2][ik1];
                    KM=fit_info.function(2,Nvar,x[j][e],Npar,tmp);
                    Kf=fit_info.function(3,Nvar,x[j][e],Npar,tmp);
                    
                    if (newline_if_w0!=x[j][e][1]){
                        fprintf(fdat,"\n\n");
                        newline_if_w0=x[j][e][1] ;
                    }
                    
                    double xi=x[j][e][10]/(4*pi_greco*x[j][e][11]);
                    xi=xi*xi;
                    double L_w=double(head[myen[e]].l1)/gJ[myen[e]].w0[j];
                    double Delta=FVE_GL_Mpi( L_w /* L/w0 */,  xi,gJ[myen[e]].f_PS_jack[0][j]  *gJ[myen[e]].w0[j]);
                    //xi*= ((1+ Delta)*(1+ Delta)) /((1-0.25 *Delta)*(1-0.25 *Delta));
                    
                    if(N==1){
                        fprintf(fdat,"%f   %f    %f   %f           %g  %g\n",x[j][e][1],xi, y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),KM,Kf);
                        printf("%f   %f    %f   %f           %g  %g\n",x[j][e][1],xi, y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),KM,Kf);
                    }
                    else if (N==2){
                        fprintf(fdat,"%f   %f    %f   %f     %f       %f       %g  %g\n",x[j][e][1],x[j][e][10],y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),y[j][e+en[0]][0]/(Kf),y[j][e+en[0]][1]/(Kf),KM,Kf);
                        printf("%f   %f    %f   %f     %f       %f       %g  %g\n",x[j][e][1],x[j][e][10],y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),y[j][e+en[0]][0]/(Kf),y[j][e+en[0]][1]/(Kf),KM,Kf);
                    }
                }
            }
            
            
            //xphys[j][0]=mud.jack_m[j] * mud.w0[j]/197.326963;
            xphys[j][10]=r1->MpiMeV[j]*mud.w0[j]/197.326963;
            
    
            xphys[j][1]=1e+6;  //w0
            xphys[j][2]=r1->MpiMeV[j]*r1->MpiMeV[j]*(gJ[0].scalefm[j]/197.326963)*(gJ[0].scalefm[j]/197.326963);
            xphys[j][3]=r1->fpiw[j];
            xphys[j][4]=1e+10;  // L such that L/w0=1e+6
            xphys[j][5]=mref[ms];
            xphys[j][6]=r1->MKMeV[j]*(gJ[0].scalefm[j]/197.326963)*r1->MKMeV[j]*(gJ[0].scalefm[j]/197.326963);//MKw2
            xphys[j][7]=0;//fkw
            xphys[j][8]=mud.jack_B[j]*gJ[0].scalefm[j]/197.326963;
            xphys[j][9]=mud.jack_f[j]*gJ[0].scalefm[j]/197.326963;
            xphys[j][10]=r1->MpiMeV[j]*gJ[0].scalefm[j]/197.326963;
            //xphys[j][11]=r1->fpiMeV_exp[j]*gJ[0].scalefm[j]/197.326963;
            xphys[j][11]=124.0*gJ[0].scalefm[j]/197.326963;
            
            /*if(j==Njack-1){
             *         for (i=0;i<10;i++)
             *            printf("xphis[%d]=%g\n",i,xphys[j][i]);
             *            printf("delta=%g  %g\n",fit_info.function(2,  Nvar, xphys[j],Npar,tmp),fit_info.function(3,  Nvar, xphys[j],Npar,tmp));       
        }*/
            //add w0 -inft
            MK[ms][j]=fit_info.function(0,  Nvar, xphys[j],Npar,tmp);//MK2
            //     printf("MKw02=%f\n",MK[ms][j]);
            if (N>1)
                MK[ms+1*nk][j]=fit_info.function(1,  Nvar, xphys[j],Npar,tmp);//fK
                
            for (i=0;i<Npar;i++){
                tmp1[i][j]=tmp[i];
            }
                
                
            free(tmp);
        } 
        struct fit_result  fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &sigmax, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C);
        
        
        char namefile[NAMESIZE];
        mysprintf(namefile,NAMESIZE,"%s_ms%d",prefix,ms );
        print_fit_info(argv, Njack, fit_out,  fit_info, xphys, gJ, head , "K",namefile);
        
        
        for (int e=0;e<myen.size();e++){
            /*free(x[e]);  free(fit[e]);*/   //moved to close fit;
            free(y1[e]);
        }
        fclose(fdat);
    }//end loop ms
    free(y1);free(rm1);
    
    printf("MKw2(ms1)=%f    MKw2(ms2)=%f     MKw2(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
    if (N>1)
        printf("fKw(ms1)=%f     fKw(ms2)=%f      fKw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);
    free(guess);
    
    ////////////////////////////////////////////////last interpolation  
    en=(int*) malloc(sizeof(int)*N);
    
    en_tot=0;
    for (n=0;n<N;n++){
        en[n]=nk;
        en_tot+=en[n];
    }
    if (N==1)
        Npar=2;
    else if (N==2)    
        Npar=4;
    Nvar=1;//m_l, w0,M_PS^2,f_PS
    
    guess=(double*) malloc(sizeof(double*)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
    
    //x=(double**) malloc(sizeof(double*)*(en_tot));
    rm=(double*) malloc(sizeof(double)*Njack);
    
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
    
    y=double_malloc_3(Njack,en_tot,2);
    
    x=double_malloc_3(Njack,en_tot,Nvar);
    sigmax=double_malloc_2(en_tot,Nvar);
    
    //the last one is the chi2
    out=double_malloc_2(N+1,Njack);
    
    for (n=0;n<N;n++){
        for (ms=0;ms<nk;ms++){
            
            if (n==0){
                for (j=0;j<Njack;j++){
                    rm[j]= 1./MK[ms][j];
                    //  rm[j]*=rm[j];
                }
                fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
            }
            if (n==1){
                for (j=0;j<Njack;j++){
                    rm[j]= MK[ms+1*nk][j];
                }
                fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
            }
            
            for (j=0;j<Njack;j++){
                y[j][ms+n*nk][0]=rm[j];
                y[j][ms+n*nk][1]=fit[ms+n*nk][1];
            }
            
            //x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
            //x[ms+n*nk][0]=mref[ms];//ml*w0
            
        }
        
    } 
    for (j=0;j<Njack;j++){
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                x[j][ms+n*nk][0]=mref[ms];//ml*w0
            }
        }
    }
    
    for (n=0;n<N;n++){
        for (ms=0;ms<nk;ms++){
            for(int v=0 ;v<Nvar;v++){
                for (j=0;j<Njack;j++)
                    rm[j]=x[j][ms+n*nk][v];
                tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                if (fabs(tmp[1])<1e-6) {
                    //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] );
                    tmp[1]=tmp[0]/1.0e+8; }
                    sigmax[ms+n*nk][v]=tmp[1];
                    free(tmp);
            }
        }
    }
    
    
    double in;
    for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_lines,guess );
        
        //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        //in=r1->MKMeV[j]*(mud.w0[j]/197.326963)*r1->MKMeV[j]*mud.w0[j]/197.326963;
        //in=r1->MKMeV[j]*mud.w0[j]/197.326963;
        in=r1->MpiMeV[j]/r1->MKMeV[j];
        in=1./(in*in);
        out[0][j]=(in-tmp[0])/tmp[1];
        //printf("%g  %g  %g  %g  \n",tmp[0],tmp[1],out[0][j],in);
        if (N>1)
            out[1][j]=tmp[2]+tmp[3]*out[0][j];
        
        out[N][j]=tmp_chi[j]/nk;//average chi2 of the three fit
        
        free(tmp);
        
    }     
    for (ms=0;ms<en_tot;ms++){
        free(fit[ms]);
    }
    
    free_3(Njack,en_tot,x);
    free_2(en_tot,sigmax);
    
    
    free(fit);     
    
    free_3(Njack,en_tot,y);
    free(guess);
    free(tmp_chi);
    
    
    
    free(rm);   
    return out;
    
} 
//////////////////////////////////////////////////////// spline in ms///
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*

double **fit_MK_fK_chiral_spline_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info , struct result_jack *r1, const char *prefix,char **argv){
double ***y,**x,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=6;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik2_min=1, ik2_max=3;
   int nk=(ik2_max-ik2_min+1);
   int ms;

   double *mref;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nk);
   mref[0]=0.064;
   mref[1]=0.080;
   mref[2]=0.095;
   int n,count,N=2;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;
   en[1]=nk;

   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
    guess[0]=2.05478;
    guess[1]=0.113021;
    
    
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double***) malloc(sizeof(double**)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double**) malloc(sizeof(double*)*ensembles);
       for(j=0;j<ensembles;j++){
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
   
for (e=0;e<ensembles;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].M_PS_jack[im][j]   *  gJ[e].w0[j];
                            rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<Njack;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms+n*nk][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                
                
            }
       
   } 


   
   for (j=0;j<Njack;j++){
       for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                x[ms+n*nk][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0

            }
       }
 
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_parabolas,guess );
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_parabolas  );
  
        for(i=0;i<Npar;i++){
            r[i][e][j]=tmp[i];
        }                
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   
} 

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   im=mass_index[0][1][0];
   //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
   for (e=0;e<ensembles;e++)
   {im=mass_index[e][1+0][0];
   printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im][Njack-1]* gJ[e].w0[Njack-1],2) ,   r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1]+mref[0]*mref[0]*r[2][e][Njack-1],  r[3][e][Njack-1]+mref[0]*r[4][e][Njack-1] +mref[0]*mref[0]*r[5][e][Njack-1]);
   
   }
    free(chi2);
///////////////////////////////////////////////////////////compute MK at physical point
   free(en);
   Nvar=fit_info.Nvar;
   Npar=fit_info.Npar;
   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
//guess[0]=2.040670;  guess[1]=0.428773;  guess[2]=0.410534;  guess[3]=0.126490;   guess[4]=-1.550172;   guess[5]=-0.026200;
  //guess[0]=1.151539;  guess[1]=0.095508;  guess[2]=0.120769;   guess[3]=-2.1775;   guess[4]=0.232919;
//guess[0]=2.074829 ;guess[1]=1.636190 ;  guess[2]=0.485904;  guess[3]=0.121129;  guess[4]=-2.204862;

  //double *xphys=(double*) malloc(sizeof(double)*(Nvar));
   

   
   
   double *rm1=(double*) malloc(sizeof(double)*(Njack));
   double **y1=(double**) malloc(sizeof(double*)*(Njack));
  // fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   MK=(double**) malloc(sizeof(double*)*(nk*N));
   for(i=0;i<nk*N;i++){
         MK[i]=(double*) malloc(sizeof(double)*Njack);
   }

   double **xphys=double_malloc_2(Njack,Nvar);
   
  
    double KM,Kf;
for (ms=0;ms<nk;ms++){
    double ***C,**tmp1;
    init_fit(  N, head , Njack, gJ, Npar, &en, &en_tot, &x, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C);
    
    count=0;
    for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        //FVE_K( r1->Bw[j], r1->fw[j], double(head[e].l1)/gJ[e].w0[j],  head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j],  mref[ms] ,gJ[e].M_PS_jack[0][j]*gJ[e].M_PS_jack[0][j],  gJ[e].f_PS_jack[0][j],r[0][e][j]+mref[ms]*r[1][e][j], r[2][e][j]+mref[ms]*r[3][e][j],&KM, &Kf);
                        //rm[j]=(r[0][e][j]+mref[ms]*r[1][e][j])  -  KM*KM*r1->Bw[j]*( head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j]+mref[ms]);//-P0_w*( mlw+msw)
                        rm[j]=(r[0][e][j]+mref[ms]*r[1][e][j] +  mref[ms]* mref[ms]*r[2][e][j]);//(KM*KM);
                        rm1[j]=(r[0][e][j]+mref[ms]*r[1][e][j] +  mref[ms]* mref[ms]*r[2][e][j]);//(KM*KM);
                    }
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[3][e][j]+mref[ms]*r[4][e][j] +  mref[ms]* mref[ms]*r[5][e][j];
                        //rm1[j]=rm[j]/Kf;
                    }
                    
                }
                fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                y1[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm1);
                
                for (j=0;j<Njack;j++){
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
               
                
                im=mass_index[e][ik2][ik1];
                x[e+count][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                x[e+count][1]=gJ[e].w0[j];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][j]*gJ[e].M_PS_jack[im][j];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][j];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
                x[e+count][5]=mref[ms];//ms*w0
                x[e+count][6]=r[0][e][j]+mref[ms]*r[1][e][j] +  mref[ms]* mref[ms]*r[2][e][j];//MKw2
                x[e+count][7]=r[3][e][j]+mref[ms]*r[4][e][j] +  mref[ms]* mref[ms]*r[5][e][j];//fkw
                x[e+count][8]=r1->Bw[j];
                x[e+count][9]=r1->fw[j];
               
            
            
               
            }
            count+=en[n];
        }

       // tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, MK_chiral_FVE ,guess );
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, MK_chiral_FVE  );
        if (j==0){ guess=guess_for_non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );    }
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        C[j]=covariance_non_linear_fit_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function );  
         
         if(j==Njack-1){
            printf("#P0_w=%f ;P1_w=%f ;  P3ww=%f;  Pf1w=%f;  Pf2w=%f; Pf4www=%f;  msw=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],mref[ms]  );
            printf("#chi2=%f\n",chi2[j]/(en_tot-Npar));
            
              printf("# w0    mlw    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
              for (e=0;e<ensembles;e++){
                  im=mass_index[e][ik2][ik1];
                  KM=fit_info.function(2,Nvar,x[e],Npar,tmp);
                  Kf=fit_info.function(3,Nvar,x[e],Npar,tmp);
         
                  printf("%f   %f    %f   %f     %f       %f       %g  %g\n",x[e][1],x[e][0],y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),y[j][e+en[0]][0]/(Kf),y[j][e+en[0]][1]/(Kf),KM,Kf);
              }
         }

       // printf("guess[0]=%f;  guess[1]=%f;  guess[2]=%f;  guess[3]=%f;   guess[4]=%f;   guess[5]=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]);
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, line  );
         

         xphys[j][0]=r1->mlw[j];
         xphys[j][1]=1e+6;  //w0
         xphys[j][2]=r1->MpiMeV[j]*r1->MpiMeV[j]*r1->w0MeV[j]*r1->w0MeV[j];
         xphys[j][3]=r1->fpiw[j];
         xphys[j][4]=1e+10;  // L such that L/w0=1e+6
         xphys[j][5]=mref[ms];
         xphys[j][6]=r1->MKMeV[j]*r1->w0MeV[j]*r1->MKMeV[j]*r1->w0MeV[j];//MKw2
         xphys[j][7]=0;//fkw
         xphys[j][8]=r1->Bw[j];
         xphys[j][9]=r1->fw[j];
         
         //add w0 -inft
         MK[ms][j]=fit_info.function(0,  Nvar, xphys[j],Npar,tmp);//MK2
       //     printf("MKw02=%f\n",MK[ms][j]);
         MK[ms+1*nk][j]=fit_info.function(1,  Nvar, xphys[j],Npar,tmp);//fK
         for (i=0;i<Npar;i++){
              tmp1[i][j]=tmp[i];
          }
       
         
         free(tmp);
    } 
    struct fit_result  fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C);
    
    
    char namefile[NAMESIZE];
    mysprintf(namefile,NAMESIZE,"%s_ms%d",prefix,ms );
    print_fit_info(argv, Njack, fit_out,  fit_info, xphys, gJ, head , "K",namefile);

   
   for (e=0;e<en_tot;e++){
           //moved to close fit;
        free(y1[e]);
   }
}
free(y1);free(rm1);

   printf("MKw2(ms1)=%f    MKw2(ms2)=%f     MKw2(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
   printf("fKw(ms1)=%f     fKw(ms2)=%f      fKw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);
   free(guess);

 ////////////////////////////////////////////////last interpolation  
   en=(int*) malloc(sizeof(int)*N);
   
   en_tot=0;
   for (n=0;n<N;n++){
       en[n]=nk;
       en_tot+=en[n];
   }
   
   Npar=6;
   Nvar=1;//m_l, w0,M_PS^2,f_PS

   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
   x=(double**) malloc(sizeof(double*)*(en_tot));
   rm=(double*) malloc(sizeof(double)*Njack);

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
   chi2=(double*) malloc(sizeof(double)*Njack);
   out=(double**) malloc(sizeof(double*)*2);
   out[0]=(double*) malloc(sizeof(double)*Njack);
   out[1]=(double*) malloc(sizeof(double)*Njack);

   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms][j];
                          //  rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms+1*nk][j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<Njack;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms+n*nk][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=mref[ms];//ml*w0
                
            }
       
   } 


   double in, msroot[1];
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_parabolas ,guess );
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_parabolas  );
        in=r1->MKMeV[j]*r1->w0MeV[j]*r1->MKMeV[j]*r1->w0MeV[j];
        out[0][j]=rtbis_func_eq_input(two_parabolas,0 ,  Nvar,x[0], Npar, tmp, 0,in ,  0.01, 0.1, 1e-10);
        if(j==Njack-1){ 
            printf("%g  vs   %g or %g \n",out[0][j],  (-tmp[1] + sqrt(tmp[1]*tmp[1]-4*tmp[2]*(tmp[0]-in)) )/(2 *tmp[2]),  (-tmp[1] - sqrt(tmp[1]*tmp[1]-4*tmp[2]*(tmp[0]-in)) )/(2 *tmp[2]) );
            printf("chi2=%g\n", chi2[Njack-1]);
        }
        //out[0][j]=(in-tmp[0])/tmp[1];
        //out[1][j]=tmp[2]+tmp[3]*out[0][j];
        msroot[0]=out[0][j];
        out[1][j]=two_parabolas(1 ,  Nvar,msroot, Npar, tmp);
                 
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   free(chi2);

   
   
   free(rm);   
   return out;
   
} 


//////////////////////////////////////////////////////////////////////////////





double **fit_MK_Mpi_fK_fpi_chiral_spline_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info , struct result_jack *r1, const char *prefix,char **argv){
double ***y,**x,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=6;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik2_min=1, ik2_max=3;
   int nk=(ik2_max-ik2_min+1);
   int ms;

   double *mref;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nk);
   mref[0]=0.064;
   mref[1]=0.080;
   mref[2]=0.095;
   int n,count,N=2;
   int *en=(int*) malloc(sizeof(int)*N);
   en[0]=nk;
   en[1]=nk;

   int en_tot=0;
   
   for (n=0;n<N;n++)
       en_tot+=en[n];
   
   double *guess=(double*) malloc(sizeof(double)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
   
    
   x=(double**) malloc(sizeof(double*)*(en_tot));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(en_tot));

   r=(double***) malloc(sizeof(double**)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double**) malloc(sizeof(double*)*ensembles);
       for(j=0;j<ensembles;j++){
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
   
for (e=0;e<ensembles;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                int imp=mass_index[e][0][0];
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].M_PS_jack[im][j]   /  gJ[e].M_PS_jack[imp][j];
                            rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].f_PS_jack[im][j]  /gJ[e].f_PS_jack[imp][j]
                            ;
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<Njack;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms+n*nk][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                
                
            }
       
   } 


   
   for (j=0;j<Njack;j++){
       for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                x[ms+n*nk][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0

            }
       }
 
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_parabolas,guess );
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_parabolas  );
  
        for(i=0;i<Npar;i++){
            r[i][e][j]=tmp[i];
        }                
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   
} 

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   im=mass_index[0][1][0];
   //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
    for (e=0;e<ensembles;e++)
   {im=mass_index[e][1+0][0];
   printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im][Njack-1]* gJ[e].w0[Njack-1],2) ,   r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1]+mref[0]*mref[0]*r[2][e][Njack-1],  r[3][e][Njack-1]+mref[0]*r[4][e][Njack-1] +mref[0]*mref[0]*r[5][e][Njack-1]);
   
   }
    free(chi2);
///////////////////////////////////////////////////////////compute MK at physical point
   free(en);
   Nvar=fit_info.Nvar;
   Npar=fit_info.Npar;
   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
   
   
   double *rm1=(double*) malloc(sizeof(double)*(Njack));
   double **y1=(double**) malloc(sizeof(double*)*(Njack));
  // fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   MK=(double**) malloc(sizeof(double*)*(nk*N));
   for(i=0;i<nk*N;i++){
         MK[i]=(double*) malloc(sizeof(double)*Njack);
   }

   double **xphys=double_malloc_2(Njack,Nvar);
   
  
    double KM,Kf;
for (ms=0;ms<nk;ms++){
 
    double ***C,**tmp1;
    init_fit(  N, head , Njack, gJ, Npar, &en, &en_tot, &x, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C);
    
    count=0;
    for (n=0;n<N;n++){
        for (e=0;e<en[n];e++){
                im=mass_index[e][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=(r[0][e][j]+mref[ms]*r[1][e][j]+mref[ms]*mref[ms]*r[2][e][j]);//(KM*KM);
                        rm1[j]=(r[0][e][j]+mref[ms]*r[1][e][j]+mref[ms]*mref[ms]*r[2][e][j]);//(KM*KM);
                    }
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[3][e][j]+mref[ms]*r[4][e][j]+mref[ms]*mref[ms]*r[5][e][j];
                        //rm1[j]=rm[j]/Kf;
                    }
                    
                }
                fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                y1[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm1);
                
                for (j=0;j<Njack;j++){
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
               
                
                im=mass_index[e][ik2][ik1];
                x[e+count][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                x[e+count][1]=gJ[e].w0[j];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][j]*gJ[e].M_PS_jack[im][j];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][j];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
                x[e+count][5]=mref[ms];//ms*w0
                x[e+count][6]=r[0][e][j]+mref[ms]*r[1][e][j];//MKw2
                x[e+count][7]=r[2][e][j]+mref[ms]*r[3][e][j];//fkw
                x[e+count][8]=r1->Bw[j];
                x[e+count][9]=r1->fw[j];
               
            }
            count+=en[n];
        }

        if (j==0){ guess=guess_for_non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );    }
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );
        chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
        C[j]=covariance_non_linear_fit_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function );  
         
         if(j==Njack-1){
            printf("#P0_w=%f ;P1_w=%f ;  P3ww=%f;  Pf1w=%f;  Pf2w=%f; Pf4www=%f;  msw=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],mref[ms]  );
            printf("#chi2=%f\n",chi2[j]/(en_tot-Npar));
            
              printf("# w0    mlw    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
              for (e=0;e<ensembles;e++){
                  im=mass_index[e][ik2][ik1];
                  KM=fit_info.function(2,Nvar,x[e],Npar,tmp);
                  Kf=fit_info.function(3,Nvar,x[e],Npar,tmp);
                 
                  printf("%f   %f    %f   %f     %f       %f       %g  %g\n",x[e][1],x[e][0],y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),y[j][e+en[0]][0]/(Kf),y[j][e+en[0]][1]/(Kf),KM,Kf);
              }
         }

 
         xphys[j][0]=r1->mlw[j];
         xphys[j][1]=1e+6;  //w0
         xphys[j][2]=r1->MpiMeV[j]*r1->MpiMeV[j]*r1->w0MeV[j]*r1->w0MeV[j];
         xphys[j][3]=r1->fpiw[j];
         xphys[j][4]=1e+10;  // L such that L/w0=1e+6
         xphys[j][5]=mref[ms];
         xphys[j][6]=r1->MKMeV[j]*r1->w0MeV[j]*r1->MKMeV[j]*r1->w0MeV[j];//MKw2
         xphys[j][7]=0;//fkw
         xphys[j][8]=r1->Bw[j];
         xphys[j][9]=r1->fw[j];
         
         //add w0 -inft
         MK[ms][j]=fit_info.function(0,  Nvar, xphys[j],Npar,tmp);//MK2
       //     printf("MKw02=%f\n",MK[ms][j]);
         MK[ms+1*nk][j]=fit_info.function(1,  Nvar, xphys[j],Npar,tmp);//fK
         for (i=0;i<Npar;i++){
              tmp1[i][j]=tmp[i];
          }
       
         
         free(tmp);
    } 
    struct fit_result  fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C);
    
    
    char namefile[NAMESIZE];
    mysprintf(namefile,NAMESIZE,"%s_ms%d",prefix,ms );
    print_fit_info(argv, Njack, fit_out,  fit_info, xphys, gJ, head , "K",namefile);

   
   for (e=0;e<en_tot;e++){
           //moved to close fit;
        free(y1[e]);
   }
}
free(y1);free(rm1);

   printf("MKw2(ms1)=%f    MKw2(ms2)=%f     MKw2(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
   printf("fKw(ms1)=%f     fKw(ms2)=%f      fKw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);
   free(guess);

 ////////////////////////////////////////////////last interpolation  
   en=(int*) malloc(sizeof(int)*N);
   
   en_tot=0;
   for (n=0;n<N;n++){
       en[n]=nk;
       en_tot+=en[n];
   }
   
   Npar=6;
   Nvar=1;//m_l, w0,M_PS^2,f_PS

   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=1.;
   guess[0]=1;guess[1]=1;
   x=(double**) malloc(sizeof(double*)*(en_tot));
   rm=(double*) malloc(sizeof(double)*Njack);

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
                            rm[j]= MK[ms][j];
                          //  rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms+1*nk][j];
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<Njack;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms+n*nk][1];
                }
                                          
                x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nk][0]=mref[ms];//ml*w0
                
            }
       
   } 


   double in,msroot[1];
   for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_parabolas,guess );
      //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        in=r1->MKMeV[j]*r1->w0MeV[j]*r1->MKMeV[j]*r1->w0MeV[j];
        in=r1->MKMeV[j]/r1->MpiMeV[j];
        in=in*in;
        out[0][j]=rtbis_func_eq_input(two_parabolas,0 ,  Nvar,x[0], Npar, tmp, 0,in ,  0.01, 0.1, 1e-10);
        if(j==Njack-1){ 
            printf("%g  vs   %g or %g \n",out[0][j],  (-tmp[1] + sqrt(tmp[1]*tmp[1]-4*tmp[2]*(tmp[0]-in)) )/(2 *tmp[2]),  (-tmp[1] - sqrt(tmp[1]*tmp[1]-4*tmp[2]*(tmp[0]-in)) )/(2 *tmp[2]) );
            printf("chi2=%g\n", chi2[Njack-1]);
        }
        //out[0][j]=(in-tmp[0])/tmp[1];
        //out[1][j]=tmp[2]+tmp[3]*out[0][j];
        msroot[0]=out[0][j];
        out[1][j]=two_parabolas(1 ,  Nvar,msroot, Npar, tmp);
                 
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nk*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   

   
   
   free(rm);   
   return out;
   
} 
*/




double **fit_MD_fD_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info , struct result_jack *r1, const char *prefix,char **argv, store_fit_clover mud ,std::vector<int> myen){
    
    double ***y,***x,**sigmax,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit;   double **out;
    
    int ensembles_D=myen.size();
    int i,j,im;  
    int n,count,N=fit_info.N;
    int Npar=4;
    double *tmp_chi2=(double*) calloc(Njack,sizeof(double));
    int Npar_inter=2;
    if (N==1)
        Npar=Npar_inter;
    else if (N==2)
        Npar=Npar_inter*2;
    
    int Nvar=1;//m_l, w0,M_PS^2,f_PS
    int ik1=0,ik2=4;
    int ik2_min=4, ik2_max=7;
    int nk=(ik2_max-ik2_min+1);
    int ms;
    int refs=3;
    
    double *mref;//[Nms]={0.52,0.68,0.81};
    mref=(double*) malloc(sizeof(double)*refs);
    mref[0]=0.94;
    mref[1]=1.04;
    mref[2]=1.08;
    //mref[3]=1.09;

    int *en=(int*) malloc(sizeof(int)*N);
    for (n=0;n<N;n++)
        en[n]=nk;
    
    
    
    int en_tot=0;
    
    for (n=0;n<N;n++)
        en_tot+=en[n];
    
    double *guess=(double*) malloc(sizeof(double)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
    
    
    //x=(double**) malloc(sizeof(double*)*(en_tot));
    
    chi2m=(double*) malloc(sizeof(double)*(Npar));
    rm=(double*) malloc(sizeof(double)*(Njack));
    fit=(double**) malloc(sizeof(double*)*(en_tot));
    
    r=(double***) malloc(sizeof(double**)*(Npar));
    for(i=0;i<Npar;i++){
        r[i]=(double**) malloc(sizeof(double*)*ensembles_D);
        for(j=0;j<ensembles_D;j++){
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
    
    for (int e=0 ; e<myen.size(); e++ ){     
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                im=mass_index[myen[e]][ms+ik2_min][ik1];
                int imp=mass_index[myen[e]][0][0];
                if (n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[myen[e]].M_PS_GEVP_jack[im][j] *gJ[myen[e]].w0[j]  ;// /  gJ[e].M_PS_jack[imp][j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=gJ[myen[e]].f_PS_ls_ss_jack[im][j]*gJ[myen[e]].w0[j];//  /gJ[e].f_PS_jack[imp][j];                            ;
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<Njack;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms+n*nk][1];
                    
                }
                
                //x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                
                
            }
            
        } 
        
        x=double_malloc_3(Njack,en_tot,Nvar);
        sigmax=double_malloc_2(en_tot,Nvar);
        
        
        for (j=0;j<Njack;j++){
            for (n=0;n<N;n++){
                for (ms=0;ms<nk;ms++){
                    im=mass_index[myen[e]][ms+ik2_min][ik1];
                    x[j][ms+n*nk][0]=head[myen[e]].k[head[myen[e]].nk+ik2_min+ms]*gJ[myen[e]].w0[j]/gJ[myen[e]].Zp[j];//ml*w0
                }
            }
        }
        
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                for(int v=0 ;v<Nvar;v++){
                    for (j=0;j<Njack;j++)
                        rm[j]=x[j][ms+n*nk][v];
                    tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                    if (fabs(tmp[1])<1e-6) {
                        //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] );
                        tmp[1]=tmp[0]/1.0e+8;
                        
                    }
                    sigmax[ms+n*nk][v]=tmp[1];
                    free(tmp);
                }
            }
        }
        
        
        
        for (j=0;j<Njack;j++){
            if (Npar==2 || Npar==4 ){
                tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_lines,guess );
                chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, two_lines  );
            }
            else if (Npar==3 || Npar==6 ){
                tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_parabolas,guess );
                chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, two_parabolas  );
            }
            if (j==Njack-1){printf("myen=%d   chi2=%g\n",myen[e],chi2[j]);  }
            //tmp=non_linear_fit_Nf_sigmax(N, en,x[j], sigmax, y[j] , Nvar,  Npar, two_lines,guess );
            //tmp=non_linear_fit_Nf_sigmax_iterative(N, en,x[j], sigmax, y[j] , Nvar,  Npar, two_lines,guess );
            
            
            for(i=0;i<Npar;i++){
                r[i][e][j]=tmp[i];
            }                
            free(tmp);
            
        }     
        for (ms=0;ms<en_tot;ms++){
            free(fit[ms]);
        }
        
    } 
    free_3(Njack,en_tot,x);
    free_2(en_tot,sigmax);
    
    
    free(fit);     
    for (j=0;j<Njack;j++){
        for (int e=0;e<nk*N;e++){
            free(y[j][e]);
        }
        free(y[j]);
    }
    free(y); free(guess);
    im=mass_index[0][1][0];
    //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
    for (int e=0; e< myen.size();e++){
       printf("myen[e]=%d\n",myen[e]);
       for (ms=0;ms<nk;ms++){
            im=mass_index[myen[e]][ms+ik2_min][ik1];
            double *M=(double*) malloc(sizeof(double)*Njack);
            double *mw=(double*) malloc(sizeof(double)*Njack);
            for (j=0;j<Njack;j++){
                M[j]=gJ[myen[e]].M_PS_GEVP_jack[im][j]*gJ[myen[e]].w0[j];
                mw[j]=head[myen[e]].k[head[myen[e]].nk+ik2_min+ms]*gJ[myen[e]].w0[j]/gJ[myen[e]].Zp[j];//ml*w0
            }
            printf("%g  %g    %g  %g \t", mw[Njack-1],  error_jackboot(jack_files[0].sampling,Njack  ,mw)  ,M[Njack-1], error_jackboot(jack_files[0].sampling,Njack  ,M) );
            printf("  %g    %g",r[0][e][Njack-1], r[1][e][Njack-1]);
            if (Npar==2 || Npar==4 ){
                printf(" \n");
            }
            else if (Npar==3 || Npar==6 ){
                printf("  %g    \n",r[2][e][Njack-1]);
            }
            free(M);
            
        }
        
    }
    free(chi2);
    ///////////////////////////////////////////////////////////compute MK at physical point
    free(en);
    Nvar=fit_info.Nvar;
    Npar=fit_info.Npar;
    guess=(double*) malloc(sizeof(double*)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
    
    
    double *rm1=(double*) malloc(sizeof(double)*(Njack));
    double **y1=(double**) malloc(sizeof(double*)*(Njack));
    // fit=(double**) malloc(sizeof(double*)*(en_tot));
    
    MK=(double**) malloc(sizeof(double*)*(refs*N));
    for(i=0;i<refs*N;i++){
        MK[i]=(double*) malloc(sizeof(double)*Njack);
    }
    
    double **xphys=double_malloc_2(Njack,Nvar);
    
    
    double KM,Kf;
    for (ms=0;ms<refs;ms++){
        
        char fname[NAMESIZE];
        mysprintf(fname,NAMESIZE,"%s/%s_ms%d",argv[2],prefix,ms);
        printf("writing data to :%s\n",fname);
        FILE *fdat=open_file(fname,"w+");
        double ***C,**tmp1;
        init_fit(  N, head , Njack, gJ,Nvar, Npar, &en, &en_tot, &x, &sigmax, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C, ensembles_D);
        
        count=0;
        for (n=0;n<N;n++){
            for (int e=0;e<en[n];e++){
                im=mass_index[myen[e]][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        if (Npar_inter==2  ){
                        rm[j]=(r[0][e][j]+mref[ms]*r[1][e][j]);//(KM*KM);
                        rm1[j]=(r[0][e][j]+mref[ms]*r[1][e][j]);//(KM*KM);
                        }
                        else if (Npar_inter==3 ){
                            rm[j]=(r[0][e][j]+mref[ms]*r[1][e][j] +mref[ms]*mref[ms]*r[2][e][j] );//(KM*KM);
                            rm1[j]=(r[0][e][j]+mref[ms]*r[1][e][j]+mref[ms]*mref[ms]*r[2][e][j]);//(KM*KM);
                        }
                    }
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        if (Npar_inter==2 ){
                            rm[j]=r[2][e][j]+mref[ms]*r[3][e][j];
                        }
                        else if (Npar_inter==3  ){
                            rm[j]=r[3][e][j]+mref[ms]*r[4][e][j]+mref[ms]*mref[ms]*r[5][e][j] ;
                        }
                        //rm1[j]=rm[j]/Kf;
                    }
                    
                }
                fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                y1[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm1);
                
                for (j=0;j<Njack;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                    //if (e==3) {y[j][e+count][1]*=100.; }
                }
                
                
                //x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                
            }
            count+=en[n];
        }   
        int ii;
        //#pragma omp parallel for  private(tmp,e,i,xphys,n,count)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)
        for (j=0;j<Njack;j++){
            count=0;
            for (n=0;n<N;n++){
                for (int e=0;e<en[n];e++){
                    
                    
                    im=mass_index[myen[e]][ik2][ik1];
                    x[j][e+count][0]=head[myen[e]].k[head[myen[e]].nk+ik1]*gJ[myen[e]].w0[j]/gJ[myen[e]].Zp[j];//ml*w0
                    x[j][e+count][1]=gJ[myen[e]].w0[j];//w0
                    x[j][e+count][2]=gJ[myen[e]].M_PS_jack[im][j]*gJ[myen[e]].M_PS_jack[im][j];//MPS^2
                    x[j][e+count][3]=gJ[myen[e]].f_PS_jack[im][j];//f_PS
                    x[j][e+count][4]=double(head[e].l1);//f_PS
                    x[j][e+count][5]=mref[ms];//ms*w0
                    x[j][e+count][6]=r[0][e][j]+mref[ms]*r[1][e][j];//MKw2
                    if (N>1)
                        x[j][e+count][7]=r[2][e][j]+mref[ms]*r[3][e][j];//fkw
                    x[j][e+count][8]=mud.jack_B[j]*mud.w0[j]/197.326963;
                    x[j][e+count][9]=mud.jack_f[j]*mud.w0[j]/197.326963;
                    
                }
                count+=en[n];
            }
        }
        /*
        count=0;
        for (n=0;n<1;n++){
            for (int e=0;e<en[n];e++){
                for(int v=0 ;v<Nvar;v++){
                    for (j=0;j<Njack;j++)
                        rm[j]=x[j][e+count][v];
                    tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                    if (fabs(tmp[0])<1e-6) {
                        //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] ); 
                        tmp[1]=tmp[0]/1.0e+8; }
                        sigmax[e+count][v]=tmp[1];
                        free(tmp);
                }
            }
            count+=en[n];
        }*/
        
        guess=guess_for_non_linear_fit_Nf(N, en,x[0], y[0] , Nvar,  Npar, fit_info.function,guess ); 
        for (j=0;j<Njack;j++){
            
            
            //if (j==0){    }
            tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess );
            //tmp=non_linear_fit_Nf_sigmax(N, en,x[j], sigmax, y[j] , Nvar,  Npar, fit_info.function,guess );
            //tmp=non_linear_fit_Nf_sigmax_iterative(N, en,x[j], sigmax, y[j] , Nvar,  Npar, fit_info.function,guess );
            chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
            C[j]=covariance_non_linear_fit_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function );  
            tmp_chi2[j]+=chi2[j];
            
            if(j==Njack-1){
                //printf("#P0_w=%f ;P1_w=%f ;  P3ww=%f;  Pf1w=%f;  Pf2w=%f; Pf4www=%f;  msw=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],mref[ms]  );
                printf("#chi2=%f\n",chi2[j]/(en_tot-Npar));
                
                printf("# w0    mlw    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
                fprintf(fdat,"# w0    mlw    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
                double newline_if_w0=x[j][0][1];
                //std::vector<int>  myen={0,1,2,3,   4,5,6,   7};
                for (int e=0 ;e<myen.size();e++){
                    im=mass_index[myen[e]][ik2][ik1];
                    KM=fit_info.function(2,Nvar,x[j][e],Npar,tmp);
                    Kf=fit_info.function(3,Nvar,x[j][e],Npar,tmp);
                    
                    if (newline_if_w0!=x[j][e][1]){
                        fprintf(fdat,"\n\n");
                        newline_if_w0=x[j][e][1] ;
                    }
                    if (N==1){
                        fprintf(fdat,"%f   %f    %f   %f           %g  %g\n",x[j][e][1], x[j][e][0], y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),KM,Kf);
                        printf("%f   %f    %f   %f           %g  %g\n",x[j][e][1], x[j][e][0], y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),KM,Kf);
                        
                    }
                    else if(N==2){
                        fprintf(fdat,"%f   %f    %f   %f     %f       %f       %g  %g\n",x[j][e][1], x[j][e][0], y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),y[j][e+en[0]][0]/(Kf),y[j][e+en[0]][1]/(Kf),KM,Kf);
                        printf("%f   %f    %f   %f     %f       %f       %g  %g\n",x[j][e][1], x[j][e][0], y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),y[j][e+en[0]][0]/(Kf),y[j][e+en[0]][1]/(Kf),KM,Kf);
                    }
                }
            }
            
            
            xphys[j][0]=mud.jack_m[j]*mud.w0[j]/197.326963;
            xphys[j][1]=1e+6;  //w0
            xphys[j][2]=r1->MpiMeV[j]*r1->MpiMeV[j]*mud.w0[j]*mud.w0[j]/(197.326963*197.326963);
            xphys[j][3]=r1->fpiw[j];
            xphys[j][4]=1e+10;  // L such that L/w0=1e+6
            xphys[j][5]=mref[ms];
            xphys[j][6]=r1->MKMeV[j]*mud.w0[j]*r1->MKMeV[j]*mud.w0[j]/(197.326963*197.326963);//MKw2
            xphys[j][7]=0;//fkw
            xphys[j][8]=mud.jack_B[j]*mud.w0[j]/197.326963;
            xphys[j][9]=mud.jack_f[j]*mud.w0[j]/197.326963;
            
            /*if(j==Njack-1){
             *         for (i=0;i<10;i++)
             *            printf("xphis[%d]=%g\n",i,xphys[j][i]);
             *            printf("delta=%g  %g\n",fit_info.function(2,  Nvar, xphys[j],Npar,tmp),fit_info.function(3,  Nvar, xphys[j],Npar,tmp));       
        }*/
            //add w0 -inft
            MK[ms][j]=fit_info.function(0,  Nvar, xphys[j],Npar,tmp);//MK2
            //     printf("MKw02=%f\n",MK[ms][j]);
            if (N>1)
                MK[ms+1*refs][j]=fit_info.function(1,  Nvar, xphys[j],Npar,tmp);//fK
            for (i=0;i<Npar;i++){
                tmp1[i][j]=tmp[i];
            }
            
            
            free(tmp);
        } 
        struct fit_result  fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &sigmax, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C);
        
        
        char namefile[NAMESIZE];
        mysprintf(namefile,NAMESIZE,"%s_ms%d",prefix,ms );
        print_fit_info(argv, Njack, fit_out,  fit_info, xphys, gJ, head , "D",namefile);
        
        
        for (int e=0;e<en_tot;e++){
            /*free(x[e]);  free(fit[e]);*/   //moved to close fit;
            free(y1[e]);
        }
        fclose(fdat);
    }//end loop ms
    free(y1);free(rm1);
    
    printf("MKw2(ms1)=%f    MKw2(ms2)=%f     MKw2(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
    if (N>1)
        printf("fKw(ms1)=%f     fKw(ms2)=%f      fKw(ms3)=%f\n", MK[0+refs][Njack-1],MK[1+refs][Njack-1],MK[2+refs][Njack-1]);
    free(guess);
    
    ////////////////////////////////////////////////last interpolation  
    en=(int*) malloc(sizeof(int)*N);
    
    en_tot=0;
    for (n=0;n<N;n++){
        en[n]=refs;
        en_tot+=en[n];
    }
    if (N==1)
        Npar=2;
    else if (N==2)
        Npar=2*2;    
        
    
    Nvar=1;//m_l, w0,M_PS^2,f_PS
    
    guess=(double*) malloc(sizeof(double*)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
    guess[0]=1;guess[1]=1;
    //x=(double**) malloc(sizeof(double*)*(en_tot));
    rm=(double*) malloc(sizeof(double)*Njack);
    
    fit=(double**) malloc(sizeof(double*)*(en_tot));
    
    y=(double***) malloc(sizeof(double**)*Njack);
    for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        for (n=0;n<N;n++){
            for (ms=0;ms<refs;ms++){
                y[j][ms+n*refs]=(double*) malloc(sizeof(double)*2);
                
            }
        }
    }
    
    x=double_malloc_3(Njack,en_tot,Nvar);
    sigmax=double_malloc_2(en_tot,Nvar);
    
    out=double_malloc_2(N+1,Njack);
    
    for (n=0;n<N;n++){
        for (ms=0;ms<refs;ms++){
            
            if (n==0){
                for (j=0;j<Njack;j++){
                    rm[j]= MK[ms][j];
                    //  rm[j]*=rm[j];
                }
                fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
            }
            if (n==1){
                for (j=0;j<Njack;j++){
                    rm[j]= MK[ms+1*refs][j];
                }
                fit[ms+n*refs]=mean_and_error(jack_files[0].sampling,Njack, rm);
            }
            
            for (j=0;j<Njack;j++){
                y[j][ms+n*refs][0]=rm[j];
                y[j][ms+n*refs][1]=fit[ms+n*refs][1];
            }
            
            //x[ms+n*refs]=(double*) malloc(sizeof(double)*Nvar);
            //x[ms+n*refs][0]=mref[ms];//ml*w0
            
        }
        
    } 
    for (j=0;j<Njack;j++){
        for (n=0;n<N;n++){
            for (ms=0;ms<refs;ms++){
                x[j][ms+n*refs][0]=mref[ms];//ml*w0
            }
        }
    }
    
    for (n=0;n<N;n++){
        for (ms=0;ms<refs;ms++){
            for(int v=0 ;v<Nvar;v++){
                for (j=0;j<Njack;j++)
                    rm[j]=x[j][ms+n*refs][v];
                tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                if (fabs(tmp[1])<1e-6) {
                    //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] );
                    tmp[1]=tmp[0]/1.0e+8; }
                    sigmax[ms+n*refs][v]=tmp[1];
                    free(tmp);
            }
        }
    }
    
    
    double in;
    for (j=0;j<Njack;j++){
        in=r1->MDMeV[j]*mud.w0[j]/197.326963;
        
        if (Npar==2 || Npar==4){
            tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_lines,guess );
            out[0][j]=(in-tmp[0])/tmp[1];
            if (N>1)
                out[1][j]=tmp[2]+tmp[3]*out[0][j];
            
        }
        else if (Npar==3 || Npar==6){
            tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_parabolas,guess );
            double c=tmp[0]-in;
            double b=tmp[1], a=tmp[2];
            out[0][j]=-b+sqrt(b*b-4.*a*c);            
            out[0][j]/=(2.*a);
            if (N>1)
                out[1][j]=tmp[3]+tmp[4]*out[0][j]+tmp[5]*out[0][j]*out[0][j];
            
            
        }
        
        out[N][j]=tmp_chi2[j]/ refs;
        
        free(tmp);
        
    }     
    for (ms=0;ms<en_tot;ms++){
        free(fit[ms]);
    }
    
    free_3(Njack,en_tot,x);
    free_2(en_tot,sigmax);
    
    
    free(fit);     
    for (j=0;j<Njack;j++){
        for (int e=0;e<refs*N;e++){
            free(y[j][e]);
        }
        free(y[j]);
    }
    free(y); free(guess);
    free(tmp_chi2);
    
    
    
    free(rm);   
    return out;
    
} 






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


double **fit_MDs_fDs_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info , struct result_jack *r1, const char *prefix,char **argv, store_fit_clover mud, store_fit_clover result_ms, std::vector<int> myen){
    double ***y,***x,**sigmax,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit;   double **out;
    
    int ensembles_D=myen.size();//myen.size();
    int i,j,im;  
    
    int Nvar=1;//m_l, w0,M_PS^2,f_PS
    int ik1=0,ik2=4;
    int ik2_min=4, ik2_max=7;
    int ik1_min=1, ik1_max=3;
    int nk1=(ik1_max-ik1_min+1);
    int nk=(ik2_max-ik2_min+1);
    int ms;
    int refs=3;
    double *tmp_chi2=(double *) calloc(Njack,sizeof(double));
    double *mref;//[Nms]={0.52,0.68,0.81};
    mref=(double*) malloc(sizeof(double)*refs);
    mref[0]=0.94;
    mref[1]=1.04;
    mref[2]=1.08;
    //mref[3]=1.09;
    int n,count,N=fit_info.N;
    int Npar=2;
    if(N>1)
        Npar=4;    
    
    int *en=(int*) malloc(sizeof(int)*N);
    for (n=0;n<N;n++)
        en[n]=nk1;
    
    
    
    int en_tot=0;
    
    for (n=0;n<N;n++)
        en_tot+=en[n];
    
    double *guess=(double*) malloc(sizeof(double)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
    guess[0]=2.05478;
    guess[1]=0.113021;
    double ****MDs_mc;
    MDs_mc=(double****) malloc(sizeof(double***)*ensembles_D); 
    for (int e=0;e<ensembles_D;e++){
        MDs_mc[e]=(double***) malloc(sizeof(double**)*nk);
        for (ms=0;ms<nk;ms++){
            MDs_mc[e][ms]=(double**) malloc(sizeof(double*)*N);
            for (i=0;i<N;i++)
                MDs_mc[e][ms][i]=(double*) malloc(sizeof(double)*Njack);
            
        }
    }
    x=double_malloc_3(Njack,en_tot,Nvar);
    
    chi2m=(double*) malloc(sizeof(double)*(Npar));
    rm=(double*) malloc(sizeof(double)*(Njack));
    fit=(double**) malloc(sizeof(double*)*(en_tot));
    
 
    r=double_malloc_3(Npar,ensembles_D,Njack);
    
    chi2=(double*) malloc(sizeof(double)*Njack);
    y=double_malloc_3(Njack,en_tot,2);
    
    for (ik2=ik2_min;ik2<=ik2_max;ik2++){
        for (int e=0;e<ensembles_D;e++){     
            for (n=0;n<N;n++){
                for (ms=0;ms<nk1;ms++){
                    im=mass_index[myen[e]][ik2][ms+ik1_min];
                    if (n==0){
                        for (j=0;j<Njack;j++){
                            rm[j]=gJ[myen[e]].M_PS_GEVP_jack[im][j]   *  gJ[myen[e]].w0[j];
                        }
                        fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                    }
                    if (n==1){
                        for (j=0;j<Njack;j++){
                            rm[j]=gJ[myen[e]].f_PS_ls_ss_jack[im][j]   *  gJ[myen[e]].w0[j];
                        }
                        fit[ms+n*nk1]=mean_and_error(jack_files[0].sampling,Njack, rm);
                    }
                    for (j=0;j<Njack;j++){
                        y[j][ms+n*nk1][0]=rm[j];
                        y[j][ms+n*nk][1]=fit[ms][1];
                        
                        x[j][ms+n*nk1][0]=head[myen[e]].k[head[myen[e]].nk+ms+ik1_min]*gJ[myen[e]].w0[j]/gJ[myen[e]].Zp[j];//ml*w0
                    }
                    
                    
                    
                    
                }
                
            }
            
            
            for (j=0;j<Njack;j++){
                
                tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_lines,guess );
                chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, two_lines  );
                
                
                MDs_mc[e][ik2-ik2_min][0][j]=tmp[0]+result_ms.ms[j]*(mud.w0[j]/197.326963)*tmp[1];
                if (N>1)
                    MDs_mc[e][ik2-ik2_min][1][j]=tmp[2]+result_ms.ms[j]*(mud.w0[j]/197.326963)*tmp[3];
                
                free(tmp);
                
            }     
            
            
        } 
    }  
    
    free(fit);     free_3(Njack,en_tot,x);

    free_3(Njack,en_tot,y);
    
    
     free(guess);
    free(mref);free(en);
    
    
    
    /////////////interpolation m_c
    
    
    
    mref=(double*) malloc(sizeof(double)*nk);
    mref[0]=0.74;
    mref[1]=0.84;
    mref[2]=0.93;
    
    
    en=(int*) malloc(sizeof(int)*N);
    
    en_tot=0;
    
    for (n=0;n<N;n++){
        en[n]=nk;
        en_tot+=en[n];
    }
    
    guess=(double*) malloc(sizeof(double)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
    guess[0]=2.05478;
    guess[1]=0.113021;
    
    
    //x=(double**) malloc(sizeof(double*)*(en_tot));
    
    chi2m=(double*) malloc(sizeof(double)*(Npar));
    rm=(double*) malloc(sizeof(double)*(Njack));
    fit=(double**) malloc(sizeof(double*)*(en_tot));
    
   
    r=double_malloc_3(Npar, ensembles_D,Njack);
    
    chi2=(double*) malloc(sizeof(double)*Njack);
    y=double_malloc_3(Njack,en_tot,2);
    
    for (int e=0;e<ensembles_D;e++){     
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                
                if (n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=MDs_mc[e][ms][0][j]  ;// /  gJ[e].M_PS_jack[imp][j];
                        
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=MDs_mc[e][ms][1][j];//  /gJ[e].f_PS_jack[imp][j];                            ;
                    }
                    fit[ms+n*nk]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                
                for (j=0;j<Njack;j++){
                    y[j][ms+n*nk][0]=rm[j];
                    y[j][ms+n*nk][1]=fit[ms+n*nk][1];
                    
                }
                
                //x[ms+n*nk]=(double*) malloc(sizeof(double)*Nvar);
                
                
            }
            
        } 
        
        x=double_malloc_3(Njack,en_tot,Nvar);
        sigmax=double_malloc_2(en_tot,Nvar);
        
        
        for (j=0;j<Njack;j++){
            for (n=0;n<N;n++){
                for (ms=0;ms<nk;ms++){
                    im=mass_index[myen[e]][ms+ik2_min][ik1];
                    x[j][ms+n*nk][0]=head[myen[e]].k[head[myen[e]].nk+ik2_min+ms]*gJ[myen[e]].w0[j]/gJ[myen[e]].Zp[j];//ml*w0
                }
            }
        }
        
        for (n=0;n<N;n++){
            for (ms=0;ms<nk;ms++){
                for(int v=0 ;v<Nvar;v++){
                    for (j=0;j<Njack;j++)
                        rm[j]=x[j][ms+n*nk][v];
                    tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                    if (fabs(tmp[1])<1e-6) {
                        //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] );
                        tmp[1]=tmp[0]/1.0e+8;
                        
                    }
                    sigmax[ms+n*nk][v]=tmp[1];
                    free(tmp);
                }
            }
        }
        
        
        
        for (j=0;j<Njack;j++){
            tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_lines,guess );
            //tmp=non_linear_fit_Nf_sigmax(N, en,x[j], sigmax, y[j] , Nvar,  Npar, two_lines,guess );
            //tmp=non_linear_fit_Nf_sigmax_iterative(N, en,x[j], sigmax, y[j] , Nvar,  Npar, two_lines,guess );
            chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, two_lines  );
            
            
            for(i=0;i<Npar;i++){
                r[i][e][j]=tmp[i];
            }                
            free(tmp);
            
        }     
        for (ms=0;ms<en_tot;ms++){
            free(fit[ms]);
        }
        
    } 
    free_3(Njack,en_tot,x);
    free_2(en_tot,sigmax);
    
    
    free(fit);     
    
    free_3(Njack,nk*N,y);
    free(guess);
    im=mass_index[0][1][0];
    //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
    for (int e=0;e<ensembles_D;e++)
    {im=mass_index[myen[e]][ik2_min+0][0];
       // printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im][Njack-1]* gJ[e].w0[Njack-1],2) ,r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1],r[2][e][Njack-1]+mref[0]*r[3][e][Njack-1] );
        
    }
    free(chi2);
    ///////////////////////////////////////////////////////////compute MK at physical point
    free(en);
    Nvar=fit_info.Nvar;
    Npar=fit_info.Npar;
    guess=(double*) malloc(sizeof(double*)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
    
    
    double *rm1=(double*) malloc(sizeof(double)*(Njack));
    double **y1=(double**) malloc(sizeof(double*)*(Njack));
    // fit=(double**) malloc(sizeof(double*)*(en_tot));
    
    MK=(double**) malloc(sizeof(double*)*(refs*N));
    for(i=0;i<refs*N;i++){
        MK[i]=(double*) malloc(sizeof(double)*Njack);
    }
    
    double **xphys=double_malloc_2(Njack,Nvar);
    
    
    double KM,Kf;
    for (ms=0;ms<refs;ms++){
        
        char fname[NAMESIZE];
        mysprintf(fname,NAMESIZE,"%s/%s_ms%d",argv[2],prefix,ms);
        printf("writing data to :%s\n",fname);
        FILE *fdat=open_file(fname,"w+");
        double ***C,**tmp1;
        init_fit(  N, head , Njack, gJ,Nvar, Npar, &en, &en_tot, &x, &sigmax, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C, ensembles_D);
        
        count=0;
        for (n=0;n<N;n++){
            for (int e=0;e<en[n];e++){
                im=mass_index[myen[e]][ik2][ik1];
                if(n==0){
                    for (j=0;j<Njack;j++){
                        rm[j]=(r[0][e][j]+mref[ms]*r[1][e][j]);//(KM*KM);
                        rm1[j]=(r[0][e][j]+mref[ms]*r[1][e][j]);//(KM*KM);
                    }
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[2][e][j]+mref[ms]*r[3][e][j];
                        //rm1[j]=rm[j]/Kf;
                    }
                    
                }
                fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                y1[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm1);
                
                for (j=0;j<Njack;j++){
                    y[j][e+count][0]=rm[j];
                    y[j][e+count][1]=fit[e+count][1];
                    //if (e==3) {y[j][e+count][1]*=100.; }
                }
                
                
                //x[e+count]=(double*) malloc(sizeof(double)*Nvar);
                
            }
            count+=en[n];
        }   
        int ii;
        //#pragma omp parallel for  private(tmp,e,i,xphys,n,count)  shared(N, en,x, y , Nvar,  Npar,guess,Njack,r,chi2)
        for (j=0;j<Njack;j++){
            count=0;
            for (n=0;n<N;n++){
                for (int e=0;e<en[n];e++){
                    
                    
                    im=mass_index[myen[e]][ik2][ik1];
                    x[j][e+count][0]=head[myen[e]].k[head[myen[e]].nk+ik1]*gJ[myen[e]].w0[j]/gJ[myen[e]].Zp[j];//ml*w0
                    x[j][e+count][1]=gJ[myen[e]].w0[j];//w0
                    x[j][e+count][2]=gJ[myen[e]].M_PS_jack[im][j]*gJ[myen[e]].M_PS_jack[im][j];//MPS^2
                    x[j][e+count][3]=gJ[myen[e]].f_PS_jack[im][j];//f_PS
                    x[j][e+count][4]=double(head[e].l1);//f_PS
                    x[j][e+count][5]=mref[ms];//ms*w0
                    x[j][e+count][6]=r[0][e][j]+mref[ms]*r[1][e][j];//MKw2
                    if (N>1)
                        x[j][e+count][7]=r[2][e][j]+mref[ms]*r[3][e][j];//fkw
                    x[j][e+count][8]=mud.jack_B[j]*mud.w0[j]/197.326963;
                    x[j][e+count][9]=mud.jack_f[j]*mud.w0[j]/197.326963;
                    
                }
                count+=en[n];
            }
        }
        /*  for sigmax
        count=0;
        for (n=0;n<1;n++){
            for (int e=0;e<en[n];e++){
                for(int v=0 ;v<Nvar;v++){
                    for (j=0;j<Njack;j++)
                        rm[j]=x[j][e+count][v];
                    tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                    if (fabs(tmp[1])<1e-6) {
                        //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] ); 
                        tmp[1]=tmp[0]/1.0e+8; }
                        sigmax[e+count][v]=tmp[1];
                        free(tmp);
                }
            }
            count+=en[n];
        }*/
        
        guess=guess_for_non_linear_fit_Nf(N, en,x[0], y[0] , Nvar,  Npar, fit_info.function,guess ); 
        for (j=0;j<Njack;j++){
            
            
            //if (j==0){    }
            tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, fit_info.function,guess );
            //tmp=non_linear_fit_Nf_sigmax(N, en,x[j], sigmax, y[j] , Nvar,  Npar, fit_info.function,guess );
            //tmp=non_linear_fit_Nf_sigmax_iterative(N, en,x[j], sigmax, y[j] , Nvar,  Npar, fit_info.function,guess );
            chi2[j]=compute_chi_non_linear_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function  )/(en_tot-Npar);
            C[j]=covariance_non_linear_fit_Nf(N, en,x[j], y[j],tmp , Nvar,  Npar, fit_info.function );  
            tmp_chi2[j]+=chi2[j];
            
            if(j==Njack-1){
                //printf("#P0_w=%f ;P1_w=%f ;  P3ww=%f;  Pf1w=%f;  Pf2w=%f; Pf4www=%f;  msw=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],mref[ms]  );
                printf("#chi2=%f\n",chi2[j]/(en_tot-Npar));
                
                printf("# w0    mlw    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
                fprintf(fdat,"# w0    mlw    MKw2/KM2  errr     fkw0/Kf    err  KM  Kf\n");
                double newline_if_w0=x[j][0][1];

                for (int e=0 ;e< ensembles_D; e++){
                    im=mass_index[myen[e]][ik2][ik1];
                    KM=fit_info.function(2,Nvar,x[j][e],Npar,tmp);
                    Kf=fit_info.function(3,Nvar,x[j][e],Npar,tmp);
                    
                    if (newline_if_w0!=x[j][e][1]){
                        fprintf(fdat,"\n\n");
                        newline_if_w0=x[j][e][1] ;
                    }
                    if (N==1){
                        fprintf(fdat,"%f   %f    %f   %f          %g  %g\n",x[j][e][1],x[j][e][0],y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),KM,Kf);
                        printf("%f   %f    %f   %f            %g  %g\n",x[j][e][1],x[j][e][0],y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),KM,Kf);
                    }
                    else if (N>1){
                        fprintf(fdat,"%f   %f    %f   %f     %f       %f       %g  %g\n",x[j][e][1],x[j][e][0],y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),y[j][e+en[0]][0]/(Kf),y[j][e+en[0]][1]/(Kf),KM,Kf);
                        printf("%f   %f    %f   %f     %f       %f       %g  %g\n",x[j][e][1],x[j][e][0],y[j][e][0]/(KM*KM),y[j][e][1]/(KM*KM),y[j][e+en[0]][0]/(Kf),y[j][e+en[0]][1]/(Kf),KM,Kf);
                    }
                    
                }
            }
            
            
            xphys[j][0]=mud.jack_m[j]*mud.w0[j]/197.326963;
            xphys[j][1]=1e+6;  //w0
            xphys[j][2]=r1->MpiMeV[j]*r1->MpiMeV[j]*mud.w0[j]*mud.w0[j]/(197.326963*197.326963);
            xphys[j][3]=r1->fpiw[j];
            xphys[j][4]=1e+10;  // L such that L/w0=1e+6
            xphys[j][5]=mref[ms];
            xphys[j][6]=r1->MKMeV[j]*mud.w0[j]*r1->MKMeV[j]*mud.w0[j]/(197.326963*197.326963);//MKw2
            xphys[j][7]=0;//fkw
            xphys[j][8]=mud.jack_B[j]*mud.w0[j]/197.326963;
            xphys[j][9]=mud.jack_f[j]*mud.w0[j]/197.326963;
            
            /*if(j==Njack-1){
             *         for (i=0;i<10;i++)
             *            printf("xphis[%d]=%g\n",i,xphys[j][i]);
             *            printf("delta=%g  %g\n",fit_info.function(2,  Nvar, xphys[j],Npar,tmp),fit_info.function(3,  Nvar, xphys[j],Npar,tmp));       
        }*/
            //add w0 -inft
            MK[ms][j]=fit_info.function(0,  Nvar, xphys[j],Npar,tmp);//MK2
            //     printf("MKw02=%f\n",MK[ms][j]);
            if (N>1)
                MK[ms+1*refs][j]=fit_info.function(1,  Nvar, xphys[j],Npar,tmp);//fK
            for (i=0;i<Npar;i++){
                tmp1[i][j]=tmp[i];
            }
            
            
            free(tmp);
        } 
        struct fit_result  fit_out=close_fit(N,  head , Njack, gJ,Npar,&en,&en_tot, &x, &sigmax, &chi2m, &rm,&tmp1, &fit, &y,&chi2,&C);
        
        
        char namefile[NAMESIZE];
        mysprintf(namefile,NAMESIZE,"%s_ms%d",prefix,ms );
        print_fit_info(argv, Njack, fit_out,  fit_info, xphys, gJ, head , "D",namefile);
        
        
        for (int e=0;e<en_tot;e++){
            /*free(x[e]);  free(fit[e]);*/   //moved to close fit;
            free(y1[e]);
        }
        fclose(fdat);
    }//end loop ms
    free(y1);free(rm1);
    
    printf("MKw2(ms1)=%f    MKw2(ms2)=%f     MKw2(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
    if(N>1)
        printf("fKw(ms1)=%f     fKw(ms2)=%f      fKw(ms3)=%f\n", MK[0+refs][Njack-1],MK[1+refs][Njack-1],MK[2+refs][Njack-1]);
    free(guess);
    
    ////////////////////////////////////////////////last interpolation  
    en=(int*) malloc(sizeof(int)*N);
    
    en_tot=0;
    for (n=0;n<N;n++){
        en[n]=refs;
        en_tot+=en[n];
    }
    if(N==1)
        Npar=2;
    if (N>1)
        Npar=4;
    Nvar=1;//m_l, w0,M_PS^2,f_PS
    
    guess=(double*) malloc(sizeof(double*)*Npar);
    for(i=0;i<Npar;i++)
        guess[i]=1.;
    guess[0]=1;guess[1]=1;
    //x=(double**) malloc(sizeof(double*)*(en_tot));
    rm=(double*) malloc(sizeof(double)*Njack);
    
    fit=(double**) malloc(sizeof(double*)*(en_tot));
    
    
    y=double_malloc_3(Njack,en_tot,2);
    
    x=double_malloc_3(Njack,en_tot,Nvar);
    sigmax=double_malloc_2(en_tot,Nvar);
    
    out=double_malloc_2(N+1,Njack);
    for (n=0;n<N;n++){
        for (ms=0;ms<refs;ms++){
            
            if (n==0){
                for (j=0;j<Njack;j++){
                    rm[j]= MK[ms][j];
                    //  rm[j]*=rm[j];
                }
                fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
            }
            if (n==1){
                for (j=0;j<Njack;j++){
                    rm[j]= MK[ms+1*refs][j];
                }
                fit[ms+n*refs]=mean_and_error(jack_files[0].sampling,Njack, rm);
            }
            
            for (j=0;j<Njack;j++){
                y[j][ms+n*refs][0]=rm[j];
                y[j][ms+n*refs][1]=fit[ms+n*refs][1];
            }
            
            //x[ms+n*refs]=(double*) malloc(sizeof(double)*Nvar);
            //x[ms+n*refs][0]=mref[ms];//ml*w0
            
        }
        
    } 
    for (j=0;j<Njack;j++){
        for (n=0;n<N;n++){
            for (ms=0;ms<refs;ms++){
                x[j][ms+n*refs][0]=mref[ms];//ml*w0
            }
        }
    }
    /*
    for (n=0;n<N;n++){
        for (ms=0;ms<refs;ms++){
            for(int v=0 ;v<Nvar;v++){
                for (j=0;j<Njack;j++)
                    rm[j]=x[j][ms+n*refs][v];
                tmp=mean_and_error(jack_files[0].sampling,Njack, rm);
                if (fabs(tmp[1])<1e-6) {
                    //printf("e=%d    v=%d   %g +- %g\n", e,v,tmp[0],tmp[1] );
                    tmp[1]=tmp[0]/1.0e+8; }
                    sigmax[ms+n*refs][v]=tmp[1];
                    free(tmp);
            }
        }
    }
    */
    
    double in;
    for (j=0;j<Njack;j++){
        tmp=non_linear_fit_Nf(N, en,x[j], y[j] , Nvar,  Npar, two_lines,guess );
        
        //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        in=r1->MDsMeV[j]*mud.w0[j]/197.326963;
        //in=r1->MKMeV[j]/r1->MpiMeV[j];
        //in=in*in;
        out[0][j]=(in-tmp[0])/tmp[1];
        if (N>1)
            out[1][j]=tmp[2]+tmp[3]*out[0][j];
        out[N][j]=tmp_chi2[j]/refs;
        
        free(tmp);
        
    }     
    for (ms=0;ms<en_tot;ms++){
        free(fit[ms]);
    }
    
    free_3(Njack,en_tot,x);
    free_2(en_tot,sigmax);
    
    
    free(fit);     
    
    free_3(Njack,en_tot,y);
    free(guess);
    free(tmp_chi2);
    
    
    
    free(rm);   
    return out;
    
} 


