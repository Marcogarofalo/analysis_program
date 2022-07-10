#define fKoverfpi_clover_C

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

 
double const_fit(int n, int Nvar, double *x,int Npar,double  *P){
    return P[0];
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
    double **tif=swap_indices(fit_info.Npar,jack_tot,fit_out.P);   
    double *r=(double*) malloc(sizeof(double)*jack_tot);
    double *m;
    double **x=(double**) malloc(sizeof(double*)*jack_tot);
    for(j=0;j<jack_tot;j++){
        x[j]=(double*) malloc(sizeof(double)*fit_info.Nvar);
        //mw=x[0], w0=x[1], dmpi2=x[2], dfpi=x[3];
        x[j][0]=0;//mlw
        x[j][1]=1e+6;//r0
        x[j][2]=1e+8;//Mpi2
        x[j][3]=1e+8;//fpi
        x[j][4]=1e+8;//L
        x[j][5]=phys_point[j][5];
        x[j][6]=phys_point[j][6];
        x[j][7]=phys_point[j][7];
        x[j][8]=phys_point[j][8];
        x[j][9]=phys_point[j][9];
        
        
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
            }
            m=mean_and_error_jack_biased(jack_tot,r);
            fprintf(fc," %g  %g \t",m[0],m[1]);
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
            }
            m=mean_and_error_jack_biased(jack_tot,r);
            fprintf(fcA," %g  %g \t",m[0],m[1]);
            free(m);
        }
        fprintf(fcA,"%g \n",gJ[0].w0[jack_tot-1]);

        fprintf(fcB,"%g \t",((double)i)*h);
        for (n=0;n<N;n++){
            for(j=0;j<jack_tot;j++){
                x[j][1]=gJ[4].w0[j];
                x[j][0]=((double)i)*h;//=1e+10;//xG
                r[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
            }
            m=mean_and_error_jack_biased(jack_tot,r);
            fprintf(fcB," %g  %g \t",m[0],m[1]);
            free(m);
        }
        fprintf(fcB,"\n");

       fprintf(fcC,"%g \t",((double)i)*h);
        for (n=0;n<N;n++){
            for(j=0;j<jack_tot;j++){
                x[j][1]=gJ[7].w0[j];
                x[j][2]=((double)i)*h;//=1e+10;//xG
                r[j]=fit_info.function(n,fit_info.Nvar,x[j],fit_info.Npar,tif[j]);
            }
            m=mean_and_error_jack_biased(jack_tot,r);
            fprintf(fcC," %g  %g \t",m[0],m[1]);
            free(m);
        }
        fprintf(fcC,"\n");

    }
    fclose(fc);
   

    free(r);
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

static void init_fit( int N, struct header *head ,int Njack, struct data_jack *gJ,int Npar, int **en,int *en_tot, double ***x, double **chi2m, double **rm, double ***r, double ***fit, double ****y,double **chi2,double ****C)
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
static struct fit_result close_fit( int N, struct header *head ,int Njack, struct data_jack *gJ,int Npar,int **en,int *en_tot, double ***x, double **chi2m, double **rm,double ***r, double ***fit, double ****y,double **chi2, double ****C)
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
 

double **fit_fKoverfpi_chiral_FVE_clover(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ,struct fit_type fit_info , struct result_jack *r1, const char *prefix,char **argv){
double ***y,**x,***r,**MK,*chi2,*tmp,*rm,*chi2m,**fit; 
   double **out;
   int i,j,e,im;  
   int Npar=2;
   int Nvar=1;//m_l, w0,M_PS^2,f_PS
   int ik1=0,ik2=1;
   int ik2_min=1, ik2_max=3;
   
   int nks=(ik2_max-ik2_min+1), nks1=1;  
   int ms;

   double *mref;//[Nms]={0.52,0.68,0.81};
   mref=(double*) malloc(sizeof(double)*nks1);
   mref[0]=r1->msw[Njack-1];
   //mref[1]=0.080;
   //mref[2]=0.095;
   int n,count,N=fit_info.N;
   int *en=(int*) malloc(sizeof(int)*N);
   
   int en_tot=0;
   
   for (n=0;n<N;n++){
       en[0]=nks;
       en_tot+=en[n];
   }
   
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
            for (ms=0;ms<(nks);ms++){
                    y[j][ms+n*(nks)]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
   
for (e=0;e<ensembles;e++){     
   for (n=0;n<N;n++){
            for (ms=0;ms<(nks);ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                int imp=mass_index[e][0][0];
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].f_PS_jack[im][j] /gJ[e].f_PS_jack[imp][j]    ;
                            
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]=gJ[e].f_PS_jack[im][j]   *  gJ[e].w0[j];
                    }
                    fit[ms+n*nks]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*(nks)][0]=rm[j];
                    y[j][ms+n*(nks)][1]=fit[ms][1];
                }
                                          
                x[ms+n*nks]=(double*) malloc(sizeof(double)*Nvar);
                
             // printf("e=%d  ms=%d   x=%g    fk/fpi=%g  ms=%g\n", e,ms,   head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,fit[ms][0],r1->msw[Njack-1]);  
            }
       
   } 


   
   for (j=0;j<Njack;j++){
       for (n=0;n<N;n++){
            for (ms=0;ms<(nks);ms++){
                im=mass_index[e][ms+ik2_min][ik1];
                x[ms+n*nks][0]=head[e].k[head[e].nk+ik2_min+ms]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0

            }
       }
        non_linear_fit_result single_jack_fit =non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, two_lines,guess );
        tmp=single_jack_fit.P;
        chi2[j]=single_jack_fit.chi2;
  
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
        for (e=0;e<(nks)*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   im=mass_index[0][1][0];
   //printf("A53: Mk(ms1)=%f   ms1=%f\n",gJ[0].M_PS_jack[im][Njack-1],head[0].k[head[0].nk+ik2_min+0]*gJ[0].w0[Njack-1] );
 /*  for (e=0;e<ensembles;e++)
   {im=mass_index[e][1+0][0];
   printf("%d   MKw2(ms=%f)=%f    MKw2=%f      fk=%f\n",e   ,  head[e].k[head[e].nk+ik2_min+0]*gJ[e].w0[Njack-1]/gJ[e].Zp[Njack-1]   ,pow(gJ[e].M_PS_jack[im][Njack-1]* gJ[e].w0[Njack-1],2) ,r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1],r[2][e][Njack-1]+mref[0]*r[3][e][Njack-1] );
   
   }*/
    for (e=0;e<ensembles;e++){
        printf("fK/fpi[e=%d]=%g  \n",e,(r[0][e][Njack-1]+mref[0]*r[1][e][Njack-1]) );
    }
    free(chi2);
///////////////////////////////////////////////////////////compute MK at physical point
 /*  en[0]=ensembles;
   en[1]=ensembles;
   en_tot=0;
   for (n=0;n<N;n++)
       en_tot+=en[n];
   */free(en);
   Nvar=fit_info.Nvar;
   Npar=fit_info.Npar;
   guess=(double*) malloc(sizeof(double*)*Npar);
   for(i=0;i<Npar;i++)
       guess[i]=rand();
//guess[0]=2.040670;  guess[1]=0.428773;  guess[2]=0.410534;  guess[3]=0.126490;   guess[4]=-1.550172;   guess[5]=-0.026200;
  //guess[0]=1.151539;  guess[1]=0.095508;  guess[2]=0.120769;   guess[3]=-2.1775;   guess[4]=0.232919;
//guess[0]=2.074829 ;guess[1]=1.636190 ;  guess[2]=0.485904;  guess[3]=0.121129;  guess[4]=-2.204862;

  //double *xphys=(double*) malloc(sizeof(double)*(Nvar));
   

   
   
   double *rm1=(double*) malloc(sizeof(double)*(Njack));
   double **y1=(double**) malloc(sizeof(double*)*(Njack));
  // fit=(double**) malloc(sizeof(double*)*(en_tot));
 
   MK=(double**) malloc(sizeof(double*)*(nks1*N));
   for(i=0;i<nks1*N;i++){
         MK[i]=(double*) malloc(sizeof(double)*Njack);
   }

   double **xphys=double_malloc_2(Njack,Nvar);
   
  
    double KM,Kf;
for (ms=0;ms<nks1;ms++){
  /*  if(ms==0)   guess[0]=2.111834 ;guess[1]=0.331624 ;  guess[2]=0.275526;  guess[3]=0.125696;  guess[4]=-1.610862; 
    if(ms==1)   guess[0]=2.111834 ;guess[1]=0.802615 ;  guess[2]=0.256295;  guess[3]=0.128061;  guess[4]=-1.639810; 
    if(ms==2)   guess[0]=2.111834 ;guess[1]=1.109423 ;  guess[2]=0.244229;  guess[3]=0.130214;  guess[4]=-1.664364; 
*/
   /* double **tmp1=double_malloc_2(Npar,Njack);
    chi2m=(double*) malloc(sizeof(double)*(Npar)); 
    chi2=(double*) malloc(sizeof(double)*Njack);
    rm=(double*) malloc(sizeof(double)*(Njack));
    x=(double**) malloc(sizeof(double*)*(en_tot));
    y=(double***) malloc(sizeof(double**)*Njack);
    for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(en_tot));
        count=0;
        for (n=0;n<N;n++){
            for (i=0;i<en[n];i++){
                y[j][i+count]=(double*) malloc(sizeof(double)*2);
            }
            count+=en[n];
        }
    }*/
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
                        rm[j]=(r[0][e][j]+r1->msw[j]*r[1][e][j]);//(KM*KM);
                        rm1[j]=(r[0][e][j]+r1->msw[j]*r[1][e][j]);//(KM*KM);
                    }
                }
                if(n==1){
                    for (j=0;j<Njack;j++){
                        rm[j]=r[2][e][j]+r1->msw[j]*r[3][e][j];
                        //rm1[j]=rm[j]/Kf;
                    }
                    
                }
                fit[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm);
                y1[e+count]=mean_and_error(jack_files[0].sampling,Njack, rm1);
                
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
               
                
                im=mass_index[e][ik2][ik1];
                x[e+count][0]=head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j];//ml*w0
                x[e+count][1]=gJ[e].w0[j];//w0
                x[e+count][2]=gJ[e].M_PS_jack[im][j]*gJ[e].M_PS_jack[im][j];//MPS^2
                x[e+count][3]=gJ[e].f_PS_jack[im][j];//f_PS
                x[e+count][4]=double(head[e].l1);//f_PS
                x[e+count][5]=r1->msw[j];//ms*w0
                x[e+count][6]=r[0][e][j]+r1->msw[j]*r[1][e][j];//MKw2
                x[e+count][7]=r1->fkw[j];//fkw
                x[e+count][8]=r1->Bw[j];
                x[e+count][9]=r1->fw[j];
               
            
            
                /*for (ii=0;ii<Nvar;ii++)
                    xphys[ii]=x[e+count][ii];*/
            }
            count+=en[n];
        }

       // tmp=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, MK_chiral_FVE ,guess );
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, MK_chiral_FVE  );
        if (j==0){ guess=guess_for_non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );    }
        non_linear_fit_result single_jack_fit =non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, fit_info.function,guess );
        tmp=single_jack_fit.P;
        chi2[j]=single_jack_fit.chi2/(en_tot-Npar);
        C[j]=covariance_non_linear_fit_Nf(N, en,x, y[j],tmp , Nvar,  Npar, fit_info.function );  
         
         if(j==Njack-1){
            //printf("#P0_w=%f ;P1_w=%f ;  P3ww=%f;  Pf1w=%f;  Pf2w=%f; Pf4www=%f;  msw=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],mref[ms]  );
            printf("#chi2=%f\n",chi2[j]);
            
              printf("# w0    mlw        fkw0/Kf    err  K\n");
              for (e=0;e<ensembles;e++){
                  im=mass_index[e][ik2][ik1];
                  //KM=fit_info.function(2,Nvar,x[e],Npar,tmp);
                  Kf=fit_info.function(1,Nvar,x[e],Npar,tmp);
                  /*if (j==jack_tot-1){
                    printf("KM=%g  Kf=%g\n",fit_info.function(2,fit_info.Nvar,x[e],fit_info.Npar,tmp),fit_info.function(3,fit_info.Nvar,x[e],fit_info.Npar,tmp));
                    printf("L_w=%f, mw=%f, fw=%f, Bw=%f\n",x[e][4], x[e][0], x[e][8], x[e][9] );
                  }*/
                 // FVE_K( r1->Bw[j], r1->fw[j], double(head[e].l1)/gJ[e].w0[j],  head[e].k[head[e].nk+ik1]*gJ[e].w0[j]/gJ[e].Zp[j],  mref[ms] ,gJ[e].M_PS_jack[0][j]*gJ[e].M_PS_jack[0][j],  gJ[e].f_PS_jack[0][j],r[0][e][j]+mref[ms]*r[1][e][j], r[2][e][j]+mref[ms]*r[3][e][j],&KM, &Kf);

                  printf("%f   %f    %f   %f     %f       \n",x[e][1],x[e][0],y[j][e][0]/(Kf),y[j][e][1]/(Kf),(Kf));
              }
         }

       // printf("guess[0]=%f;  guess[1]=%f;  guess[2]=%f;  guess[3]=%f;   guess[4]=%f;   guess[5]=%f;\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]);
       // chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, line  );
         

         xphys[j][0]=r1->mlw[j];
         xphys[j][1]=1e+6;  //w0
         xphys[j][2]=r1->MpiMeV[j]*r1->MpiMeV[j]*r1->w0MeV[j]*r1->w0MeV[j];
         xphys[j][3]=r1->fpiw[j];
         xphys[j][4]=1e+10;  // L such that L/w0=1e+6
         xphys[j][5]=r1->msw[j];
         xphys[j][6]=r1->MKMeV[j]*r1->w0MeV[j]*r1->MKMeV[j]*r1->w0MeV[j];//MKw2
         xphys[j][7]=0;//fkw
         xphys[j][8]=r1->Bw[j];
         xphys[j][9]=r1->fw[j];
     
         
         //add w0 -inft
         MK[ms][j]=fit_info.function(0,  Nvar, xphys[j],Npar,tmp);//MK2
       //     printf("MKw02=%f\n",MK[ms][j]);
       //  MK[ms+1*nk][j]=fit_info.function(1,  Nvar, xphys[j],Npar,tmp);//fK
                       
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
        /*free(x[e]);  free(fit[e]);*/   //moved to close fit;
        free(y1[e]);
   }
}
free(y1);free(rm1);

 //  printf("MKw2(ms1)=%f    MKw2(ms2)=%f     MKw2(ms3)=%f\n", MK[0][Njack-1],MK[1][Njack-1],MK[2][Njack-1]);
 //  printf("fKw(ms1)=%f     fKw(ms2)=%f      fKw(ms3)=%f\n", MK[0+nk][Njack-1],MK[1+nk][Njack-1],MK[2+nk][Njack-1]);
   printf("fKw/fpi(ms1)=%f\n",MK[0][Njack-1]);
 /*free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<en_tot;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); */free(guess);

 ////////////////////////////////////////////////last interpolation  
   en=(int*) malloc(sizeof(int)*N);
   
   en_tot=0;
   for (n=0;n<N;n++){
       en[n]=nks1;
       en_tot+=en[n];
   }
   
   Npar=1;
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
            for (ms=0;ms<nks1;ms++){
                    y[j][ms+n*nks1]=(double*) malloc(sizeof(double)*2);
                
            }
        }
   }
   out=(double**) malloc(sizeof(double*)*2);
   out[0]=(double*) malloc(sizeof(double)*Njack);
   out[1]=(double*) malloc(sizeof(double)*Njack);

   for (n=0;n<N;n++){
            for (ms=0;ms<nks1;ms++){
               
                if (n==0){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms][j];
                          //  rm[j]*=rm[j];
                    }
                    fit[ms]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }
                if (n==1){
                    for (j=0;j<Njack;j++){
                            rm[j]= MK[ms+1*nks1][j];
                    }
                    fit[ms+n*nks1]=mean_and_error(jack_files[0].sampling,Njack, rm);
                }

                for (j=0;j<jack_tot;j++){
                    y[j][ms+n*nks1][0]=rm[j];
                    y[j][ms+n*nks1][1]=fit[ms][1];
                }
                                          
                x[ms+n*nks1]=(double*) malloc(sizeof(double)*Nvar);
                x[ms+n*nks1][0]=r1->msw[Njack-1];//ml*w0
                
            }
       
   } 


   double in;
   for (j=0;j<Njack;j++){
        non_linear_fit_result single_jack_fit=non_linear_fit_Nf(N, en,x, y[j] , Nvar,  Npar, const_fit,guess );
        tmp = single_jack_fit.P;
       // printf("j=%d    tmp=%g     MK=%g\n",j,tmp[0], y[j][0][0]);
      //  chi2[j]=compute_chi_non_linear_Nf(N, en,x, y[j],tmp , Nvar,  Npar, two_lines  );
        /*in=r1->MKMeV[j]*r1->w0MeV[j]*r1->MKMeV[j]*r1->w0MeV[j];
        out[0][j]=(in-tmp[0])/tmp[1];
        out[1][j]=tmp[2]+tmp[3]*out[0][j];
         */
        out[0][j]=tmp[0];
        free(tmp);

   }     
   for (ms=0;ms<en_tot;ms++){
        free(x[ms]);  free(fit[ms]);
   }
   

  
   free(fit);     free(x);
   for (j=0;j<Njack;j++){
        for (e=0;e<nks1*N;e++){
             free(y[j][e]);
         }
         free(y[j]);
   }
   free(y); free(guess);
   

   
   
   free(rm);   
   return out;
   
} 

