#define CONTROL

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
#include "pion.hpp"
#include "KandD.hpp"
#include "global_fit_KandD.hpp"
#include "tower.hpp"
#include "mutils.hpp"

#include "fve.hpp"


int ***mass_index;

double ***init_omega_jacob(){
    int e,ik;
    double ***r;
    
    r=double_malloc_3(ensembles,3, 2);
    if(ensembles>6){
        //C60
        r[6][0][0]=0.558; r[6][0][1]=0.007;  
        r[6][1][0]=0.586; r[6][1][1]=0.005;    
        r[6][2][0]=0.609; r[6][2][1]=0.005;  
    }if(ensembles>5){
        //B072
        r[5][0][0]=0.682; r[5][0][1]=0.006;  
        r[5][1][0]=0.701; r[5][1][1]=0.005;  
        r[5][2][0]=0.719; r[5][2][1]=0.004;  
    }if(ensembles>4){
        //B25
        r[4][0][0]=0.662; r[4][0][1]=0.006;  
        r[4][1][0]=0.699; r[4][1][1]=0.004;  
        r[4][2][0]=0.723; r[4][2][1]=0.004;  
    }if(ensembles>3){
        //A12
        r[3][0][0]=0.787; r[3][0][1]=0.004;  
        r[3][1][0]=0.818; r[3][1][1]=0.007;  
        r[3][2][0]=0.853; r[3][2][1]=0.006;  
    }if(ensembles>2){
        //A30
        r[2][0][0]=0.804; r[2][0][1]=0.005;  
        r[2][1][0]=0.840; r[2][1][1]=0.003;  
        r[2][2][0]=0.873; r[2][2][1]=0.003;  
    }if(ensembles>1){
        //A40
        r[1][0][0]=0.821; r[1][0][1]=0.006;  
        r[1][1][0]=0.852; r[1][1][1]=0.005;  
        r[1][2][0]=0.882; r[1][2][1]=0.006;  
    }if(ensembles>0){
        //A53
        r[0][0][0]=0.811; r[0][0][1]=0.007;  
        r[0][1][0]=0.845; r[0][1][1]=0.005;  
        r[0][2][0]=0.879; r[0][2][1]=0.004;  
    }
    return r;
    
}

int ***init_mass_index_ave_r(struct header *head)
{
     int k1, k2,i;
     int nk,e;
     int ***mass_index;
     
     mass_index=(int ***)  malloc(sizeof(int**)*ensembles);
     for (e=0;e<ensembles;e++){
        nk=head[e].nk;
        mass_index[e]=(int**) malloc(sizeof(int*)*nk);
        for (k2=0;k2<nk;k2++){
            mass_index[e][k2]=(int*) malloc(sizeof(int)*(k2+1));
        }

        i=0;
        for (k2=0;k2<nk;k2++)
        for (k1=0;k1<=k2;k1++)
        {
            mass_index[e][k2][k1]=i;
            i++;
        }

     }
     return mass_index;
}

double global_fK_from_M(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi,r;
    double pi=3.141592653589793;
    
    double Mpiw=x[0], MKw=x[1], w0=x[2], Mpi2=x[3], fpi=x[4], frac_Lw=x[7],  Bw=x[8];
    double fw=x[9],  MK2=x[5], fK=x[6];
    
    //double    P0_w=Bw, P1_w=P[0], P3ww=P[1];
    //double   Pf1w=P[0],  Pf2w=P[1],  Pf4www=P[2];
    
    double KM=1.,Kf=1.0;
    
   //FVE_K( Bw, fw, frac_Lw,  Mpiw*Mpiw/(2.*Bw)/*mlw*/, MKw*MKw/Bw-Mpiw*Mpiw/(2.*Bw) /*msw*/ ,Mpi2,  fpi,MK2, fK,&KM, &Kf);
        
    //    fKw=Pf1w*( 1.- (3./2.)* xi*log(xi)+Pf2w*xi+(1/(w0*w0))*Pf4www)*Kf;
    
    xi=Mpiw*Mpiw/(16*pi*pi*fw*fw);
    double P1=P[0]+P[3]*MKw*MKw;
    double P2=P[1]+P[4]*MKw*MKw;
    double P4=P[2]+P[5]*MKw*MKw;
    
    r=P1*(1.- (3./2.)* xi*log(xi)+ P2*xi +  P4 *(1/(w0*w0))    )*Kf;
    
    return r;
    
}

double  global_Omega_Mpi_MK(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi,r;
    double pi=3.141592653589793;
    
    double Mpiw=x[0], MKw=x[1], w0=x[2];
    
    r=P[0]+P[1]*Mpiw*Mpiw+P[2]*MKw*MKw+P[3]/(w0*w0);
    
    return r;
    
}
double  global_Omega_MK(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi,r;
    double pi=3.141592653589793;
    
    double Mpiw=x[0], MKw=x[1], w0=x[2];
    
    r=P[0]+P[1]*MKw*MKw+P[2]/(w0*w0);
    
    return r;
    
}
double  global_Omega_propMK(int n, int Nvar, double *x,int Npar,double  *P){
    
    double fKw=0,xi,r;
    double pi=3.141592653589793;
    
    double Mpiw=x[0], MKw=x[1], w0=x[2];
    
    r=(P[0])*MKw*MKw+P[1]/(w0*w0);
    
    return r;
    
}

static void  read_file_head_jack(FILE *stream,struct header *head)
{
    int i;
    
    fread(&(head->twist),sizeof(int),1,stream);
    fread(&(head->nf),sizeof(int),1,stream);
    fread(&(head->nsrc),sizeof(int),1,stream);
    fread(&(head->l0),sizeof(int),1,stream);
    fread(&(head->l1),sizeof(int),1,stream);
    fread(&(head->l2),sizeof(int),1,stream);
    fread(&(head->l3),sizeof(int),1,stream);
    fread(&(head->nk),sizeof(int),1,stream);
    fread(&(head->nmoms),sizeof(int),1,stream);
    
    fread(&(head->beta),sizeof(double),1,stream);
    fread(&(head->ksea),sizeof(double),1,stream);
    fread(&(head->musea),sizeof(double),1,stream);
    fread(&(head->csw),sizeof(double),1,stream);
   
    head->k=(double*) malloc(sizeof(double)*2*head->nk);
    for(i=0;i<2*head->nk;++i)
    	fread(&(head->k[i]),sizeof(double),1,stream);
    
    head->mom=(double**) malloc(sizeof(double*)*head->nmoms);
    for(i=0;i<head->nmoms;i++) {
    	head->mom[i]=(double*) malloc(sizeof(double)*4);
        fread(&(head->mom[i][0]),sizeof(double),1,stream);
        fread(&(head->mom[i][1]),sizeof(double),1,stream);
        fread(&(head->mom[i][2]),sizeof(double),1,stream);
        fread(&(head->mom[i][3]),sizeof(double),1,stream);

    }
}

void setup_reading_single_jack( struct  database_file_jack *jack_files, struct header *head){
     int N;
     
     jack_files->f_M_PS=fopen(jack_files->M_PS,"r");
     error(jack_files->f_M_PS==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->M_PS);
     read_file_head_jack(jack_files->f_M_PS,head);
     fread(&(jack_files->Njack),sizeof(int),1,jack_files->f_M_PS);
     
     jack_files->f_M_PS_GEVP=fopen(jack_files->M_PS_GEVP,"r");
     error(jack_files->f_M_PS_GEVP==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->M_PS_GEVP);
     read_file_head_jack(jack_files->f_M_PS_GEVP,head);
     fread(&(jack_files->Njack),sizeof(int),1,jack_files->f_M_PS_GEVP);
     
     jack_files->f_f_PS=fopen(jack_files->f_PS,"r");
     error(jack_files->f_f_PS==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->f_PS);
     read_file_head_jack(jack_files->f_f_PS,head);
     fread(&N,sizeof(int),1,jack_files->f_f_PS);
     
     jack_files->f_f_PS_ls_ss=fopen(jack_files->f_PS_ls_ss,"r");
     error(jack_files->f_f_PS_ls_ss==NULL,1,"setup_reading_single_jack",
         "Unable to open output file %s",jack_files->f_PS_ls_ss);
     read_file_head_jack(jack_files->f_f_PS_ls_ss,head);
     fread(&N,sizeof(int),1,jack_files->f_f_PS_ls_ss);
     
     error(jack_files->Njack!=N,1,"setup_reading_single_jack", " files \n %s has %d elements \n %s has %d elements\n ",
         jack_files->M_PS_GEVP, jack_files->Njack,jack_files->f_PS,N); 
     
}

void  setup_reading_jack(char **argv,struct  database_file_jack *jack_files, struct header *head,const char  *name)  {

    
        mysprintf(jack_files->M_PS,NAMESIZE,"%s/M_{PS}_%s",name,argv[1]);
        mysprintf(jack_files->f_PS,NAMESIZE,"%s/Zf_{PS}_%s",name,argv[1]);  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  I rename for simplicity  Zf_PS->f_PS !!!!!!!!!!!!!!!!!!!!!!!!!1
        

        mysprintf(jack_files->M_PS_GEVP,NAMESIZE,"%s/M_{PS}^{GEVP}_%s",name,argv[1]);
        mysprintf(jack_files->f_PS_ls_ss,NAMESIZE,"%s/f_{PS}_ls_ss_%s",name,argv[1]);         
        setup_reading_single_jack(jack_files, head);

}



void  files_declarations(char **argv,struct database_file_jack *jack_files,struct header *head){
  
        

    setup_reading_jack( argv,&jack_files[0],&head[0],"../../beta1.726/cA211ab.53.24/analysis/main/jackknife");  
    jack_files[0].a=0.096;
    if (ensembles>1){
        setup_reading_jack(argv, &jack_files[1],&head[1],"../../beta1.726/cA211ab.40.24/analysis/main/jackknife");  
        jack_files[1].a=0.096;
    }
    if (ensembles>2){
        setup_reading_jack(argv, &jack_files[2],&head[2],"../../beta1.726/cA211ab.30.32/analysis/main/jackknife");  
        jack_files[2].a=0.096;
    }
    if (ensembles>3){
     //   setup_reading_jack(argv, &jack_files[3],&head[3],"../../beta1.726/cA211ab.12.48/analysis/main/jackknife");  
        setup_reading_jack(argv, &jack_files[3],&head[3],"../../beta1.726/cA211ab.12.48_no_rew/analysis/main/jackknife");  
        jack_files[3].a=0.096;
    }
    if (ensembles>4){
        setup_reading_jack(argv, &jack_files[4],&head[4],"../../beta1.778/cB211ab.25.48/analysis/main/jackknife");  
        jack_files[4].a=0.081;
    }
    if (ensembles>5){
        setup_reading_jack(argv, &jack_files[5],&head[5],"../../beta1.778/cB211ab.072.64/analysis/main/jackknife");  
        jack_files[5].a=0.081;
    }
    
    if (ensembles>6){
        setup_reading_jack(argv, &jack_files[6],&head[6],"../../beta1.836/cC211ab.06.80/analysis/main/jackknife");  
        jack_files[6].a=0.070;
    }
    if (ensembles>7){
        setup_reading_jack(argv, &jack_files[7],&head[7],"../../beta1.778/cB211ab.14.64/analysis/main/jackknife");  
        jack_files[7].a=0.070;
    }
    
    
    //take the header form ensemble 4
   /* if (ensembles>5){
        jack_files[5].a=0.096;jack_files[5].Njack=100;     
        fseek(jack_files[4].f_M_PS_GEVP,0,SEEK_SET);
        read_file_head_jack(jack_files[4].f_M_PS_GEVP,&head[5]);
        head[5].k[head[5].nk]=0.0012;
        head[5].musea=0.0012;
        head[5].l0=96;head[5].l1=48;head[5].l2=48;head[5].l3=48;
    }*/
    
 /*   if (ensembles>6){
        jack_files[6].a=0.070;jack_files[6].Njack=100;
        fseek(jack_files[4].f_M_PS_GEVP,0,SEEK_SET);
        read_file_head_jack(jack_files[4].f_M_PS_GEVP,&head[6]);
        head[6].k[head[6].nk]=0.0006;
        head[5].musea=0.0006;
        head[6].l0=160;head[6].l1=80;head[6].l2=80;head[6].l3=80;
    }
  
    if (ensembles>6){
        fseek(jack_files[4].f_M_PS_GEVP,sizeof(int),SEEK_CUR);
    }*/
}
/*
void  read_files_jack( struct database_file_jack *jack_files, struct header *head,int ***mass_index, double ***M_PS_GEVP_jack, double ***f_PS_jack){
      
      int i,ik1,ik2;
      
      for(i=0;i<ensembles;i++){
            M_PS_GEVP_jack[i]=(double**) malloc(sizeof(double*)*head[i].nk*head[i].nk);
            f_PS_jack[i]=(double**) malloc(sizeof(double*)*head[i].nk*head[i].nk);
            for(ik1=0;ik1<2;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[i].nk;ik2++){
   
                M_PS_GEVP_jack[i][mass_index[i][ik2][ik1]]=(double*) malloc(sizeof(double)*jack_files[i].Njack);
                f_PS_jack[i][mass_index[i][ik2][ik1]]=(double*) malloc(sizeof(double)*jack_files[i].Njack);

                fread(M_PS_GEVP_jack[i][mass_index[i][ik2][ik1]],   sizeof(double),   jack_files[i].Njack,   jack_files[i].f_M_PS_GEVP );
                fread(f_PS_jack[i][mass_index[i][ik2][ik1]],        sizeof(double),   jack_files[i].Njack,  jack_files[i].f_f_PS );
            }
            }
            
      }
    
}
*/

double *create_fake_jack_Zp_M12GEv(int ensemble, struct database_file_jack *jack_files)
{
    double *tmp;
    
    
     if (fabs(jack_files[ensemble].a-0.096)<0.000001){
         //tmp=fake_sampling(jack_files[0].sampling,0.483,0.004,jack_files[ensemble].Njack);// M1 2GeV report 16/2/2019  a^2g*^2
        // tmp=fake_sampling(jack_files[0].sampling,0.474,0.005,jack_files[ensemble].Njack);// M1 2GeV report 15/06/2018 tilde p
        // tmp=fake_sampling(jack_files[0].sampling,0.501,0.003,jack_files[ensemble].Njack);//M2 2Gev  16/2/2019  a^2g*^2
         
       //  tmp=fake_sampling(jack_files[0].sampling,0.459,0.005,jack_files[ensemble].Njack);//M1 2Gev  uncorrected
       //  tmp=fake_sampling(jack_files[0].sampling,0.527,0.004,jack_files[ensemble].Njack);//M2 2Gev  uncorrected
          //tmp=fake_sampling(jack_files[0].sampling, 0.485,0.005,jack_files[ensemble].Njack);//M1 2Gev  29/3/2019  a^2g*^2
         tmp=fake_sampling(jack_files[0].sampling,0.502,0.004,jack_files[ensemble].Njack,rand());//M2 2Gev  29/3/2019  a^2g*^2
          //tmp=fake_sampling(jack_files[0].sampling, 0.482,0.005,jack_files[ensemble].Njack);//M1 2Gev  29/3/2019  a^inf g_0^2
         // tmp=fake_sampling(jack_files[0].sampling, 0.530,0.004,jack_files[ensemble].Njack);//M2 2Gev  29/3/2019  a^inf g_0^2
          
         //         tmp=fake_sampling(jack_files[0].sampling,0.495,0.004,jack_files[ensemble].Njack);

     }
     if (fabs(jack_files[ensemble].a-0.081)<0.000001){
         //tmp=fake_sampling(jack_files[0].sampling, 0.478,0.002,jack_files[ensemble].Njack);// M1 2GeV report 16/2/2019  a^2g*^2
         //tmp=fake_sampling(jack_files[0].sampling,0.469,0.004,jack_files[ensemble].Njack);// M1 2GeV report 15/06/2018  tilde p
    //tmp=fake_sampling(jack_files[0].sampling,0.496,0.002,jack_files[ensemble].Njack);//M2 2Gev  16/2/2019  a^2g*^2

      //  tmp=fake_sampling(jack_files[0].sampling,0.471,0.007,jack_files[ensemble].Njack);//M1 2Gev  uncorrected
      //  tmp=fake_sampling(jack_files[0].sampling,0.502,0.005,jack_files[ensemble].Njack);//M2 2Gev  uncorrected
         //tmp=fake_sampling(jack_files[0].sampling,0.484,0.007,jack_files[ensemble].Njack);//M1 2Gev  29/3/2019  a^2g*^2
        tmp=fake_sampling(jack_files[0].sampling,0.491,0.005,jack_files[ensemble].Njack,rand());//M2 2Gev  29/3/2019  a^2g*^2
       // tmp=fake_sampling(jack_files[0].sampling, 0.485,0.007,jack_files[ensemble].Njack);//M1 2Gev  29/3/2019  a^inf g_0^2
        //tmp=fake_sampling(jack_files[0].sampling, 0.509,0.005,jack_files[ensemble].Njack);//M2 2Gev  29/3/2019  a^inf g_0^2

          //     tmp=fake_sampling(jack_files[0].sampling,0.496,0.006,jack_files[ensemble].Njack);

     }
     if (fabs(jack_files[ensemble].a-0.070)<0.000001){
      //  tmp=fake_sampling(jack_files[0].sampling,0.492333,0.003,jack_files[ensemble].Njack);// extrapolated from M2 2Gev  16/2/2019  linear
        tmp=fake_sampling(jack_files[0].sampling,0.497,0.003,jack_files[ensemble].Njack,rand());// extrapolated from M2 2Gev  16/2/2019  RF

     }
     
     return tmp;
}
double *create_fake_jack_w0(int ensemble, struct database_file_jack *jack_files)
{
    double *tmp;
    
    
     if (fabs(jack_files[ensemble].a-0.096)<0.000001)
        tmp=fake_sampling(jack_files[0].sampling,1.8346689, 0.005178046,jack_files[ensemble].Njack,rand());
     
     if (fabs(jack_files[ensemble].a-0.081)<0.000001)
        tmp=fake_sampling(jack_files[0].sampling,2.1330729,0.00468807,jack_files[ensemble].Njack,rand());
     
     if (fabs(jack_files[ensemble].a-0.070)<0.000001)
        tmp=fake_sampling(jack_files[0].sampling,2.49879971,0.0034,jack_files[ensemble].Njack,rand());
     
     return tmp;
}


void  read_files_jack( struct database_file_jack *jack_files, struct header *head,int ***mass_index, struct data_jack *dataJ){
      
      int i,ik1,ik2;
      
      for(i=0;i<ensembles;i++){
            
          
            dataJ[i].M_PS_jack=(double**) malloc(sizeof(double*)*head[i].nk*head[i].nk);
            dataJ[i].f_PS_jack=(double**) malloc(sizeof(double*)*head[i].nk*head[i].nk);

            dataJ[i].M_PS_GEVP_jack=(double**) malloc(sizeof(double*)*head[i].nk*head[i].nk);
            dataJ[i].f_PS_ls_ss_jack=(double**) malloc(sizeof(double*)*head[i].nk*head[i].nk);


            for(ik1=0;ik1<4;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[i].nk;ik2++){
               
                dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=(double*) malloc(sizeof(double)*jack_files[i].Njack);
                dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=(double*) malloc(sizeof(double)*jack_files[i].Njack);
                dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=(double*) malloc(sizeof(double)*jack_files[i].Njack);
                dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=(double*) malloc(sizeof(double)*jack_files[i].Njack);
                
                
                fread(dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]],   sizeof(double),   jack_files[i].Njack,   jack_files[i].f_M_PS );
                fread(dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]],        sizeof(double),   jack_files[i].Njack,  jack_files[i].f_f_PS );
                fread(dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]],   sizeof(double),   jack_files[i].Njack,   jack_files[i].f_M_PS_GEVP );
                fread(dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]],        sizeof(double),   jack_files[i].Njack,  jack_files[i].f_f_PS_ls_ss );
                /*
                if (ik2==0 && ik1==0){
                
               
                 if (i==0){
                  // dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.166249064969,4.35388233333e-4,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0713600049579,1.84703451511e-4,jack_files[i].Njack);
                   
                   dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.166249064969,4.35388233333e-4,jack_files[i].Njack);
                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0713600049579,1.84703451511e-4,jack_files[i].Njack);
                  
                }
               else if (i==1){
                   //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.145705225353, 1.19948022076e-3,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0678279826442,2.25713230671e-4,jack_files[i].Njack);
                   
                   dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.145705225353, 1.19948022076e-3,jack_files[i].Njack);
                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0678279826442,2.25713230671e-4,jack_files[i].Njack);
                }
                else if (i==2){
                  //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.125775584839, 2.94951334969e-4,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0666157700528, 1.34015838616e-4,jack_files[i].Njack);
                     dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.125775584839, 2.94951334969e-4,jack_files[i].Njack);
                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0666157700528, 1.34015838616e-4,jack_files[i].Njack);
                    
                }
                else if (i==3){
                   //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0799,0.0002,jack_files[i].Njack);
                   dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0799,0.0002,jack_files[i].Njack);

                
                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0622,0.0005,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0622,0.0005,jack_files[i].Njack);

                }
               else if (i==4){
                   //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.10429304904, 9.99583163124e-5,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0576087558838, 7.3009468564e-5,jack_files[i].Njack);
                   
                    dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.10429304904, 9.99583163124e-5,jack_files[i].Njack);
                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0576087558838, 7.3009468564e-5,jack_files[i].Njack);
                }
                
                else if (i==5){
                   //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0567528, 0.0000691,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0526446, 0.0000852,jack_files[i].Njack);
                   
                   dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0567528, 0.0000691,jack_files[i].Njack);
                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0526446, 0.0000852,jack_files[i].Njack);
                }
                else if (i==6){
                   //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0799,0.0002,jack_files[i].Njack);
                   dataJ[i].M_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,4.723748e-02,8.6e-5,jack_files[i].Njack);

                   //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0782,0.0002,jack_files[i].Njack);//lower by supposed FSE of pion mass splitting
                  //dataJ[i].M_PS_GEVP_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0799,0.02,jack_files[i].Njack);//larger error to eliminate it from the fit

                   //dataJ[i].f_PS_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,0.0622,0.0005,jack_files[i].Njack);
                   dataJ[i].f_PS_ls_ss_jack[mass_index[i][ik2][ik1]]=fake_sampling(jack_files[0].sampling,4.496516e-2,1.3653e-4,jack_files[i].Njack);

                }
                
                }*/
            }
            }
      }
 
}

double *lin_a_lin_mu(int np1, int var,double in){
    double *r;
    int i,j;
    
    r=(double*) malloc(sizeof(double)*(np1));

    for (i=0;i<np1;i++){
        r[i]=1.;
        for (j=0;j<i;j++){
            r[i]*=in;
            
        }
    }

    return r;
}
double **fit_nth_mass_la_lmu(int ik2,int ik1,struct database_file_jack  *jack_files,  struct header *head ,int Njack, double ***M_PS_GEVP_jack_tot  ){
       double ***y,**x,**r,*chi2,**tmp,**rm,*chi2m,**fit; 
   int i,j,im;  
   int npara=2;
   
   
   x=(double**) malloc(sizeof(double*)*(ensembles));

   chi2m=(double*) malloc(sizeof(double)*(npara));
   rm=(double**) malloc(sizeof(double*)*(npara));
   fit=(double**) malloc(sizeof(double*)*(ensembles));

   r=(double**) malloc(sizeof(double*)*(npara));
   for(i=0;i<npara;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(ensembles));
        for (i=0;i<ensembles;i++){
            y[j][i]=(double*) malloc(sizeof(double)*2);
        }
   }
   printf("a[fm]     mu[fm^-1]     M_PS^2 [MeV] err\n");
   for (i=0;i<ensembles;i++){
        im=mass_index[i][ik2][ik1];
        fit[i]=mean_and_error_jack_biased(Njack, M_PS_GEVP_jack_tot[i][im]);
        for (j=0;j<Njack;j++){ 
            y[j][i][0]=M_PS_GEVP_jack_tot[i][im][j]*197.326963/jack_files[i].a;
            y[j][i][0]*=y[j][i][0];
            y[j][i][1]=2*fit[i][1]*M_PS_GEVP_jack_tot[i][im][Njack-1];
            y[j][i][1]*=197.326963*197.326963/(jack_files[i].a*jack_files[i].a);
        }
        x[i]=(double*) malloc(sizeof(double)*npara);
        x[i][0]=jack_files[i].a*jack_files[i].a;
        x[i][1]=head[i].k[head[i].nk+ik2]/(jack_files[i].a);
         printf("%g     %g  %g   %g\n",jack_files[i].a,x[i][1],y[Njack-1][i][0],fit[i][1]);
   }   

   for (j=0;j<Njack;j++){
        tmp=global_linear_fit( ensembles, x, y[j],  npara );

        chi2[j]=compute_chisqr_global_fit( ensembles, x, y[j],npara, tmp  );
         
        for(i=0;i<npara;i++){
            r[i][j]=tmp[i][0];
            free(tmp[i]);
        }                

   }   
           free(tmp);
           /*
    y[0][0][0]=2.95996;  y[0][0][1]=0.0004388 ;
    y[0][1][0]=2.2355 ;   y[0][1][1]=0.000653378;
    y[0][2][0]=1.68077;   y[0][2][1]=0.000341078;
    y[0][3][0]=1.6609  ;  y[0][3][1]=0.000251249 ; 
    y[0][4][0]=0.481672;  y[0][4][1]=0.000141917;
           tmp=global_linear_fit( ensembles, x, y[0],  npara );
           
    printf("compare with gnuplot=(%g+-%g)  +  (%g+-%g)  a^2 +  (%g+-%g) \\mu\n",tmp[0][0],tmp[0][1],tmp[1][0],tmp[1][1],tmp[2][0],tmp[2][1] );
        
*/           
   chi2m=mean_and_error_jack_biased(Njack, chi2);
  
   //////free////
  /* for(i=0;i<npara;i++){
        free(rm[i]);       
    }*/
    free(rm);free(chi2m);
    
    for (i=0;i<ensembles;i++){
        free(fit[i]);
        free(x[i]);   

    }
    free(fit);

    free(x);
    for (j=0;j<Njack;j++){
            for (i=0;i<ensembles;i++){
                free(y[j][i]);
            }
            free(y[j]);
    }
    free(y);
    return r;

}
//P[0]=B,  P[1]=f, P[2]=P1, P[3]=P2
//x[0]=m_l, x[1]=w0,  x[2]=KM
double fun_Mw2( int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double Bw=P[0], fw=P[1], P1=P[2], P2=P[3];
    double ml=x[0], w0=x[1], KM=x[2];
    
    xi=2*Bw*ml*w0/(16.*pi*pi*fw*fw);
    
    Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2;
    Mw2*=2*Bw*w0*ml*KM*KM;
    
    return Mw2;
}

//P[0]=B,  P[1]=f, P[2]=P1, P[3]=P2
//x[0]=m_l, x[1]=w0,  x[2]=KM
double *fun_Mw2_k( int Nvar, double *x,int Npar,double  *P){
    
    double Mw2=0,xi;
    double pi=3.141592653589793;
    double *Mw2_k=(double*) calloc(Npar,sizeof(double));
    double Bw=P[0], fw=P[1], P1=P[2], P2=P[3];
    double ml=x[0], w0=x[1], KM=x[2];
    
    xi=2*Bw*ml*w0/(16.*pi*pi*fw*fw);
    
    Mw2=1+xi*log(xi)+P1*xi+ (1./(w0*w0))*P2;
    Mw2*=2*Bw*w0*ml*KM;
    
    Mw2_k[0]=Mw2/Bw+ 2*Bw*w0*ml*KM* ( log(xi)+ 1+ P1  )*(xi/Bw);
    Mw2_k[1]= 2*Bw*w0*ml*KM* ( log(xi)+ 1+ P1  )*(-2*xi/fw);
    Mw2_k[2]=2*Bw*w0*ml*KM*(xi);
    Mw2_k[3]=2*Bw*w0*ml*KM*(1/(w0*w0));
    
    return Mw2_k;
    
}




void mud(double xi,int Npar,double *P, double *f, double *df){
    
    double h=0.001;
  
    
    *f=m_over_f_xi( xi, Npar,P);
    
    xi=xi-2.*h;
    *df=m_over_f_xi( xi, Npar,P);
    
    xi=xi+h;
    *df-=8*m_over_f_xi( xi, Npar,P);
    
    xi=xi+2*h;
    *df+=8*m_over_f_xi( xi, Npar,P);
    
    xi=xi+h;
    *df-=m_over_f_xi( xi, Npar,P);
    
    xi=xi-2*h;
    *df/=(12.*h);
    
}




double **fit_Mpi(struct database_file_jack  *jack_files,  struct header *head ,int Njack,int ***mass_index, struct data_jack *gJ ){
   double ***y,**x,**r,*chi2,*tmp,*rm,*chi2m,**fit; 
   int i,j,e,im;  
   int Npar=4;
   int Nvar=3;//m_l, w0,KM
   int ik1=0,ik2=0;
   
   x=(double**) malloc(sizeof(double*)*(ensembles));

   chi2m=(double*) malloc(sizeof(double)*(Npar));
   rm=(double*) malloc(sizeof(double*)*(Njack));
   fit=(double**) malloc(sizeof(double*)*(ensembles));

   r=(double**) malloc(sizeof(double*)*(Npar));
   for(i=0;i<Npar;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(ensembles));
        for (i=0;i<ensembles;i++){
             y[j][i]=(double*) malloc(sizeof(double)*2);
        }
   }
   
   printf("w0/a[fm]     am/aZp[]   KM   (M_Pi w0)^2 [MeV] err\n");
   for (e=0;e<ensembles;e++){
        im=mass_index[e][ik2][ik1];
        for (j=0;j<Njack;j++){
             rm[j]=gJ[e].M_PS_GEVP_jack[im][j]   *  gJ[e].w0[j];
             rm[j]*=rm[j];
        }
        fit[e]=mean_and_error_jack_biased(Njack, rm);
        for (j=0;j<jack_tot;j++){
             y[j][e][0]=rm[j];
             y[j][e][1]=fit[e][1];
        }
        
                        
        x[e]=(double*) malloc(sizeof(double)*Nvar);
        x[e][0]=head[e].k[head[e].nk+ik2]/gJ[e].Zp[Njack-1];//ml
        x[e][1]=gJ[e].w0[Njack-1];//w0
        x[e][2]=gJ[e].KM[im];//KM
        printf("%g     %g    %g   %g   %g\n",x[e][1],x[e][0],x[e][2],fit[e][0],fit[e][1]);
        free(fit[e]);
   }   

   for (j=0;j<Njack;j++){
        tmp=non_linear_fit( ensembles,x, y[j] , Nvar,  Npar, fun_Mw2 );

        chi2[j]=compute_chi_non_linear( ensembles,x, y[j],tmp , Nvar,  Npar, fun_Mw2  );
              
        for(i=0;i<Npar;i++){
            r[i][j]=tmp[i];
        }                
        free(tmp);

   }         

   chi2m=mean_and_error_jack_biased(Njack, chi2);

   free(rm);free(chi2m);
   for (e=0;e<ensembles;e++){
        free(x[e]);   

   }
   free(fit);  

   free(x);
   for (j=0;j<Njack;j++){
           for (e=0;e<ensembles;e++){
               free(y[j][e]);
           }
           free(y[j]);
   }
   free(y); 
   return r;
    
}



double **fit_nth_f_PS_la_lmu(int ik2,int ik1,struct database_file_jack  *jack_files,  struct header *head ,int Njack, double ***M_PS_GEVP_jack_tot  ){
       double ***y,**x,**r,*chi2,**tmp,**rm,*chi2m,**fit; 
   int i,j,im;  
   int npara=3;
   
   
   x=(double**) malloc(sizeof(double*)*(ensembles));

   chi2m=(double*) malloc(sizeof(double)*(npara));
   rm=(double**) malloc(sizeof(double*)*(npara));
   fit=(double**) malloc(sizeof(double*)*(ensembles));

   r=(double**) malloc(sizeof(double*)*(npara));
   for(i=0;i<npara;i++){
       r[i]=(double*) malloc(sizeof(double)*Njack);
   }
   
   chi2=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);

   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(ensembles));
        for (i=0;i<ensembles;i++){
            y[j][i]=(double*) malloc(sizeof(double)*2);
        }
   }
   printf("a[fm]     mu[fm^-1]     f_PS [MeV] err\n");
   for (i=0;i<ensembles;i++){
        im=mass_index[i][ik2][ik1];
        fit[i]=mean_and_error_jack_biased(Njack, M_PS_GEVP_jack_tot[i][im]);
        for (j=0;j<Njack;j++){ 
            y[j][i][0]=M_PS_GEVP_jack_tot[i][im][j]*197.326963/jack_files[i].a;
            y[j][i][1]=fit[i][1];
            y[j][i][1]/=197.326963/(jack_files[i].a);
        }
        x[i]=(double*) malloc(sizeof(double)*npara);
        x[i][0]=1.;
        x[i][1]=jack_files[i].a*jack_files[i].a;
        x[i][2]=head[i].k[head[i].nk+ik2]/(jack_files[i].a);
         printf("%g     %g  %g   %g\n",jack_files[i].a,x[i][2],y[Njack-1][i][0],fit[i][1]);
   }   

   for (j=0;j<Njack;j++){
        tmp=global_linear_fit( ensembles, x, y[j],  npara );

        chi2[j]=compute_chisqr_global_fit( ensembles, x, y[j],npara, tmp  );
         
        for(i=0;i<npara;i++){
            r[i][j]=tmp[i][0];
            free(tmp[i]);
        }                

   }   
           free(tmp);
           /*
    y[0][0][0]=2.95996;  y[0][0][1]=0.0004388 ;
    y[0][1][0]=2.2355 ;   y[0][1][1]=0.000653378;
    y[0][2][0]=1.68077;   y[0][2][1]=0.000341078;
    y[0][3][0]=1.6609  ;  y[0][3][1]=0.000251249 ; 
    y[0][4][0]=0.481672;  y[0][4][1]=0.000141917;
           tmp=global_linear_fit( ensembles, x, y[0],  npara );
           
    printf("compare with gnuplot=(%g+-%g)  +  (%g+-%g)  a^2 +  (%g+-%g) \\mu\n",tmp[0][0],tmp[0][1],tmp[1][0],tmp[1][1],tmp[2][0],tmp[2][1] );
        
*/           
   chi2m=mean_and_error_jack_biased(Njack, chi2);
  
   //////free////
  /* for(i=0;i<npara;i++){
        free(rm[i]);       
    }*/
    free(rm);free(chi2m);
    
    for (i=0;i<ensembles;i++){
        free(fit[i]);
        free(x[i]);   

    }
    free(fit);

    free(x);
    for (j=0;j<Njack;j++){
            for (i=0;i<ensembles;i++){
                free(y[j][i]);
            }
            free(y[j]);
    }
    free(y);
    return r;

}
void  create_fake_distribution(const char *jackboot,double **w0A,double **w0B,double **w0C,double **ZpA,double **ZpB,double **ZpC,int jack_tot){
      //w0A=fake_sampling(jackboot,1.8346689, 0.005178046,*jack_tot);
      //w0A=fake_sampling(jackboot,1.83005,3.48173757101e-3,*jack_tot,rand());// MG fit in M_PS^2
      //(*w0A)=fake_sampling(jackboot,1.8355,0.0042,jack_tot,123); // PD fit
      (*w0A)=fake_sampling(jackboot,1.83616,0.0050,jack_tot,123);  // MG fit 26/11/2119  in M_PS^2/f_PS^2
      //(*w0A)=fake_sampling(jackboot,1.33582,0.001038,jack_tot,123); //t/w0
      
      //w0B=fake_sampling(jackboot,2.1330729,0.00468807,*jack_tot,rand());// MG fit in M_PS^2
      //(*w0B)=fake_sampling(jackboot,2.1347,0.0047,jack_tot,1234);// MG fit in M_PS^2
      //(*w0B)=fake_sampling(jackboot,2.12650,0.0023,jack_tot,1234);// MG fit 26/11/2119  in M_PS^2/f_PS^2
      (*w0B)=fake_sampling(jackboot,2.12842 ,0.002013,jack_tot,1234);// MG fit 26/03/2020  in M_PS^2/f_PS^2
      //(*w0B)=fake_sampling(jackboot,1.52764 ,0.000342,jack_tot,1234);// t/w0
      
      //(*w0C)=fake_sampling(jackboot,2.49879971,0.0034,jack_tot,12345);// MG fit in M_PS^2
      (*w0C)=fake_sampling(jackboot,2.50310,0.0018,jack_tot,12345);//  MG fit 26/11/2119  in M_PS^2/f_PS^2
      //(*w0C)=fake_sampling(jackboot,1.77671,0.00048,jack_tot,12345);//  t/w0

      //  ZpA=fake_sampling(jackboot,0.459,0.005,*jack_tot);//M1 2Gev  uncorrected
      //  ZpA=fake_sampling(jackboot,0.527,0.004,*jack_tot);//M2 2Gev  uncorrected
      //  ZpA=fake_sampling(jackboot, 0.485,0.005,*jack_tot);//M1 2Gev  29/3/2019  a^2g*^2
      //ZpA=fake_sampling(jackboot,0.502,0.004,*jack_tot);//M2 2Gev  29/3/2019  a^2g*^2
      //  ZpA=fake_sampling(jackboot, 0.482,0.005,*jack_tot);//M1 2Gev  29/3/2019  a^inf g_0^2
      //  ZpA=fake_sampling(jackboot, 0.530,0.004,*jack_tot);//M2 2Gev  29/3/2019  a^inf g_0^2
      //ZpA=fake_sampling(jackboot,0.471,0.005,*jack_tot);//M1 2Gev  fiorenza  
      //ZpA=fake_sampling(jackboot,0.491,0.004,*jack_tot);//M2a 2Gev  fiorenza  
      //ZpA=fake_sampling(jackboot, 0.508,0.003,*jack_tot,rand());//M2b 2Gev  fiorenza  
      // (*ZpA)=fake_sampling(jackboot, 0.4770,0.0045,jack_tot,321);//M1 2Gev  Petros  2/9/2019
       //(*ZpA)=fake_sampling(jackboot, 0.4628,0.0052,jack_tot,321);//M1 2Gev  Petros  15/12/2019  quadratic
       //(*ZpA)=fake_sampling(jackboot, 0.478,0.004,jack_tot,321);//M1 2Gev  Petros  15/12/2019   constant
       (*ZpA)=fake_sampling(jackboot, 0.474,0.002,jack_tot,321);//M1 2Gev  Enrico-Matteo  19/05/2020  constant

      
      //  ZpB=fake_sampling(jackboot,0.471,0.007,*jack_tot);//M1 2Gev  uncorrected
      //  ZpB=fake_sampling(jackboot,0.502,0.005,*jack_tot);//M2 2Gev  uncorrected
      //  ZpB=fake_sampling(jackboot,0.484,0.007,*jack_tot);//M1 2Gev  29/3/2019  a^2g*^2
      //ZpB=fake_sampling(jackboot,0.491,0.005,*jack_tot);//M2 2Gev  29/3/2019  a^2g*^2
      //  ZpB=fake_sampling( jackboot,0.485,0.007,*jack_tot);//M1 2Gev  29/3/2019  a^inf g_0^2
      //  ZpB=fake_sampling( jackboot,0.509,0.005,*jack_tot);//M2 2Gev  29/3/2019  a^inf g_0^2
     // ZpB=fake_sampling(jackboot,0.476,0.008,*jack_tot);//M1 2Gev  fiorenza
      //ZpB=fake_sampling(jackboot,0.486,0.005,*jack_tot);//M2a 2Gev  fiorenza
      //ZpB=fake_sampling(jackboot,0.500,0.003,*jack_tot,rand());//M2b 2Gev  fiorenza
      //(*ZpB)=fake_sampling(jackboot, 0.4860,0.0070,jack_tot,3214);//M1 2Gev  Petros  2/9/2019
      //(*ZpB)=fake_sampling(jackboot, 0.4780,0.0070,jack_tot,3214);//M1 2Gev  Petros  15/12/2019 quadratic 
      //(*ZpB)=fake_sampling(jackboot, 0.487,0.0040,jack_tot,3214);//M1 2Gev  Petros  15/12/2019  constant
      (*ZpB)=fake_sampling(jackboot, 0.482,0.0030,jack_tot,3214);//M1 2Gev Enrico-Matteo  19/05/2020  constant

     //  ZpC=fake_sampling(jackboot,0.492333,0.003,*jack_tot);// extrapolated from M2 2Gev  16/2/2019  linear
      //ZpC=fake_sampling(jackboot,0.497,0.003,*jack_tot);// extrapolated from M2 2Gev  16/2/2019  RF
      //(*ZpC)=fake_sampling(jackboot, 0.4860,0.0030,jack_tot,32145);//M1 2Gev  Petros  15/12/2019 quadratic
      //(*ZpC)=fake_sampling(jackboot, 0.484,0.0060,jack_tot,32145);//M1 2Gev  Petros  15/12/2019 constant
      (*ZpC)=fake_sampling(jackboot, 0.493,0.0030,jack_tot,32145);//M1 2Gev  Enrico-Matteo  19/05/2020  constant

      int seedw0=213;
      result.w0fm=fake_sampling(jackboot,v_w0fm,err_w0fm,jack_tot,seedw0);
      result.w0MeV=fake_sampling(jackboot,v_w0MeV,err_w0fm/197.326963 ,jack_tot,seedw0);
      result.MpiMeV=fake_sampling(jackboot,v_MpiMeV,err_MpiMeV,jack_tot,2134);
      result.MKMeV=fake_sampling(jackboot,v_MKMeV,err_MKMeV,jack_tot,21345);
      result.MDMeV=fake_sampling(jackboot,v_MDMeV,err_MDMeV,jack_tot,321);
      result.MDsMeV=fake_sampling(jackboot,v_MDsMeV,err_MDsMeV,jack_tot,3124);
      result.fpiMeV_exp=fake_sampling(jackboot,v_fpiMeV_exp,err_fpiMeV_exp,jack_tot,31245);
      result.MOmegaMeV=fake_sampling(jackboot,v_MOmegaMeV,err_MOmegaMeV,jack_tot,111);
    
}


struct data_jack *create_generalised_boot( struct database_file_jack *jack_files, struct header *head,int *jack_tot,int ***mass_index, struct data_jack *dJ){
      int j,e,e1,ik1,ik2,counter;
      int im;
      double ***M_PS_GEVP_jack_tot;
      struct data_jack *gJ;
      double ***omega;
      omega=init_omega_jacob();
      
      gJ=(struct data_jack *) malloc (sizeof(struct data_jack )*ensembles);
      
      *jack_tot=0;
      for(e=0;e<ensembles-1;e++){
          error(jack_files[e].Njack!=jack_files[e+1].Njack,1,"create_generalised_boot","bootstrap of the file %d has different number",e+1 );
      }
      *jack_tot=jack_files[0].Njack;

      for(e=0;e<ensembles;e++){
          gJ[e].M_PS_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          gJ[e].f_PS_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          
          gJ[e].M_PS_GEVP_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          gJ[e].f_PS_ls_ss_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          
          gJ[e].KM=(double*) malloc(sizeof(double)*head[e].nk*head[e].nk);
          gJ[e].Kf=(double*) malloc(sizeof(double)*head[e].nk*head[e].nk);
          
            for(ik1=0;ik1<4;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[e].nk;ik2++){   
                 im=mass_index[e][ik2][ik1];
                 gJ[e].M_PS_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 gJ[e].f_PS_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 
                 gJ[e].M_PS_GEVP_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 gJ[e].f_PS_ls_ss_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 
                 
                 
                 gJ[e].KM[im]=dJ[e].KM[im];
                 gJ[e].Kf[im]=dJ[e].Kf[im];
                 
                 
                 
                 for(j=0;j<(*jack_tot);j++){
                           
                           gJ[e].M_PS_jack[im][j]=dJ[e].M_PS_jack[  im ][j];
                           gJ[e].f_PS_jack[im][j]=dJ[e].f_PS_jack[  im ][j];
                           
                           gJ[e].M_PS_GEVP_jack[im][j]=dJ[e].M_PS_GEVP_jack[  im ][j];
                           gJ[e].f_PS_ls_ss_jack[im][j]=dJ[e].f_PS_ls_ss_jack[  im ][j];
                           
                 }
                 
            }}
            gJ[e].M_Omega_jack=(double**) malloc(sizeof(double*)*3);
            gJ[e].M_Omega_jack[0]= fake_sampling(jack_files[0].sampling, omega[e][0][0],omega[e][0][1],*jack_tot,e);        
            gJ[e].M_Omega_jack[1]= fake_sampling(jack_files[0].sampling, omega[e][1][0],omega[e][1][1],*jack_tot,e);        
            gJ[e].M_Omega_jack[2]= fake_sampling(jack_files[0].sampling, omega[e][2][0],omega[e][2][1],*jack_tot,e);        
            
      }
      double *w0A,*w0B,*w0C, *ZpA,*ZpB,*ZpC;
      create_fake_distribution(jack_files[0].sampling, &w0A, &w0B, &w0C, &ZpA, &ZpB, &ZpC,*jack_tot);
      
      
      if (ensembles>0){
        gJ[0].w0=w0A;
        gJ[0].Zp=ZpA;
      }
      if (ensembles>1){
        gJ[1].w0=w0A;
        gJ[1].Zp=ZpA;        
      }
      if (ensembles>2){
        gJ[2].w0=w0A;
        gJ[2].Zp=ZpA;        
      }
      if (ensembles>3){
        gJ[3].w0=w0A;
        gJ[3].Zp=ZpA;
      }
      if (ensembles>4){
        gJ[4].w0=w0B;
        gJ[4].Zp=ZpB;
      }
      if (ensembles>5){
        gJ[5].w0=w0B;
        gJ[5].Zp=ZpB;
      }
      
      if (ensembles>6){
        gJ[6].w0=w0C;
        gJ[6].Zp=ZpC;
      }
      if (ensembles>7){
        gJ[7].w0=w0B;
        gJ[7].Zp=ZpB;
      }
      
     // free(w0A); free(w0B); free(w0C); free(ZpA); free(ZpB); free(ZpC);
      
      
      for(e=0;e<ensembles;e++){
            for(ik1=0;ik1<2;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[e].nk;ik2++){
               free(dJ[e].M_PS_jack[ mass_index[e][ik2][ik1] ]);
               free(dJ[e].f_PS_jack[ mass_index[e][ik2][ik1] ]);
               free(dJ[e].M_PS_GEVP_jack[ mass_index[e][ik2][ik1] ]);
               free(dJ[e].f_PS_ls_ss_jack[ mass_index[e][ik2][ik1] ]);

            }
            }
            free(dJ[e].KM);free(dJ[e].Kf);
            free(dJ[e].M_PS_jack); free(dJ[e].M_PS_GEVP_jack); free(dJ[e].f_PS_jack); free(dJ[e].f_PS_ls_ss_jack);
            free(omega[e][0]);free(omega[e][1]);free(omega[e][2]);
      }
      free(dJ);          free(omega);  
      
     return gJ;
}
struct data_jack *create_generalised_jack( struct database_file_jack *jack_files, struct header *head,int *jack_tot,int ***mass_index, struct data_jack *dJ){
      int j,e,e1,ik1,ik2,counter;
      int im;
      double ***M_PS_GEVP_jack_tot;
      struct data_jack *gJ;
      double ***omega;
      omega=init_omega_jacob();
 
      gJ=(struct data_jack *) malloc (sizeof(struct data_jack )*ensembles);
      
      *jack_tot=0;
      for(e=0;e<ensembles;e++){
          *jack_tot+=jack_files[e].Njack;
      }
      *jack_tot=*jack_tot-ensembles+1;

      for(e=0;e<ensembles;e++){
          gJ[e].M_PS_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          gJ[e].f_PS_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          
          gJ[e].M_PS_GEVP_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          gJ[e].f_PS_ls_ss_jack=(double**) malloc(sizeof(double*)*head[e].nk*head[e].nk);
          
          
          gJ[e].KM=(double*) malloc(sizeof(double)*head[e].nk*head[e].nk);
          gJ[e].Kf=(double*) malloc(sizeof(double)*head[e].nk*head[e].nk);
            for(ik1=0;ik1<4;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[e].nk;ik2++){   
                 im=mass_index[e][ik2][ik1];
                 gJ[e].M_PS_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 gJ[e].f_PS_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 
                 gJ[e].M_PS_GEVP_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 gJ[e].f_PS_ls_ss_jack[mass_index[e][ik2][ik1]]=(double*) calloc(*jack_tot,sizeof(double));
                 
                 gJ[e].KM[im]=dJ[e].KM[im];
                 gJ[e].Kf[im]=dJ[e].Kf[im];
                 
                 counter=0;
                 for(e1=0;e1<ensembles;e1++){
                      for(j=0;j<(jack_files[e1].Njack-1);j++){
                           if (e==e1){
                           gJ[e].M_PS_jack[im][j+counter]=dJ[e].M_PS_jack[  im ][j];
                           gJ[e].f_PS_jack[im][j+counter]=dJ[e].f_PS_jack[  im ][j];
                           
                           gJ[e].M_PS_GEVP_jack[im][j+counter]=dJ[e].M_PS_GEVP_jack[  im ][j];
                           gJ[e].f_PS_ls_ss_jack[im][j+counter]=dJ[e].f_PS_ls_ss_jack[  im ][j];
                           }
                           else{ 
                           gJ[e].M_PS_jack[im][j+counter]=dJ[e].M_PS_jack[im][ jack_files[e].Njack-1 ];
                           gJ[e].f_PS_jack[im][j+counter]=dJ[e].f_PS_jack[im][ jack_files[e].Njack-1 ];   
                           
                           gJ[e].M_PS_GEVP_jack[im][j+counter]=dJ[e].M_PS_GEVP_jack[im][ jack_files[e].Njack-1 ];
                           gJ[e].f_PS_ls_ss_jack[im][j+counter]=dJ[e].f_PS_ls_ss_jack[im][ jack_files[e].Njack-1 ];    
                           }
                      }
                      counter+=jack_files[e1].Njack-1;
                 }
                 gJ[e].M_PS_jack[im][*jack_tot-1]=dJ[e].M_PS_jack[im][jack_files[e].Njack-1];
                 gJ[e].f_PS_jack[im][*jack_tot-1]=dJ[e].f_PS_jack[im][jack_files[e].Njack-1];
                 gJ[e].M_PS_GEVP_jack[im][*jack_tot-1]=dJ[e].M_PS_GEVP_jack[im][jack_files[e].Njack-1];
                 gJ[e].f_PS_ls_ss_jack[im][*jack_tot-1]=dJ[e].f_PS_ls_ss_jack[im][jack_files[e].Njack-1];
                 
            }}

            counter=0;
            gJ[e].M_Omega_jack=(double**) malloc(sizeof(double*)*3);
            gJ[e].M_Omega_jack[0]= fake_sampling(jack_files[0].sampling, omega[e][0][0],omega[e][0][1],*jack_tot,e);        
            gJ[e].M_Omega_jack[1]= fake_sampling(jack_files[0].sampling, omega[e][1][0],omega[e][1][1],*jack_tot,e);        
            gJ[e].M_Omega_jack[2]= fake_sampling(jack_files[0].sampling, omega[e][2][0],omega[e][2][1],*jack_tot,e);        

            /*gJ[e].Zp=(double*) malloc(sizeof(double)*(*jack_tot));//jack_tot is a pointer here
            gJ[e].w0=(double*) malloc(sizeof(double)*(*jack_tot));//jack_tot is a pointer here
            for(e1=0;e1<ensembles;e1++){
                for(j=0;j<(jack_files[e1].Njack-1);j++){
                    if (e==e1){
                        gJ[e].Zp[j+counter]=dJ[e].Zp[j];
                        gJ[e].w0[j+counter]=dJ[e].w0[j];

                    }
                    else{ 
                        gJ[e].Zp[j+counter]=dJ[e].Zp[ jack_files[e].Njack-1 ];//jack_tot is a pointer here
                        gJ[e].w0[j+counter]=dJ[e].w0[ jack_files[e].Njack-1 ];
                           
                    }
                }
                counter+=jack_files[e1].Njack-1;
            } 
            gJ[e].Zp[*jack_tot-1]=dJ[e].Zp[ jack_files[e].Njack-1 ];
            gJ[e].w0[*jack_tot-1]=dJ[e].w0[ jack_files[e].Njack-1 ];*/
            
      }
      double *w0A,*w0B,*w0C, *ZpA,*ZpB,*ZpC;
      create_fake_distribution(jack_files[0].sampling, &w0A, &w0B, &w0C, &ZpA, &ZpB, &ZpC,*jack_tot);
      
      if (ensembles>0){
        gJ[0].w0=w0A;
        gJ[0].Zp=ZpA;
      }
      if (ensembles>1){
        gJ[1].w0=w0A;
        gJ[1].Zp=ZpA;        
      }
      if (ensembles>2){
        gJ[2].w0=w0A;
        gJ[2].Zp=ZpA;        
      }
      if (ensembles>3){
        gJ[3].w0=w0A;
        gJ[3].Zp=ZpA;
      }
      if (ensembles>4){
        gJ[4].w0=w0B;
        gJ[4].Zp=ZpB;
      }
      if (ensembles>5){
        gJ[5].w0=w0B;
        gJ[5].Zp=ZpB;
      }
      
      if (ensembles>6){
        gJ[6].w0=w0C;
        gJ[6].Zp=ZpC;
      }
      if (ensembles>7){
        gJ[7].w0=w0B;
        gJ[7].Zp=ZpB;
      }
      for(e=0;e<ensembles;e++){
            for(ik1=0;ik1<2;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[e].nk;ik2++){
               free(dJ[e].M_PS_jack[ mass_index[e][ik2][ik1] ]);
               free(dJ[e].f_PS_jack[ mass_index[e][ik2][ik1] ]);
               free(dJ[e].M_PS_GEVP_jack[ mass_index[e][ik2][ik1] ]);
               free(dJ[e].f_PS_ls_ss_jack[ mass_index[e][ik2][ik1] ]);

            }
            }
            free(dJ[e].KM);free(dJ[e].Kf);
            free(dJ[e].M_PS_jack); free(dJ[e].M_PS_GEVP_jack); free(dJ[e].f_PS_jack); free(dJ[e].f_PS_ls_ss_jack);
            free(omega[e][0]);free(omega[e][1]);free(omega[e][2]);
      }
      free(dJ);            free(omega);
      
     return gJ;
}
//table 3 of arXiv:hep-lat/0503014  fit: log R= A+B*M+C*L // R=(M_L-M_inf)/M_inf //  K*M_inf=M_L
void KM_FSE(struct database_file_jack *jack_files, struct header *head, struct data_jack *dJ){
    double RM,L,M;
    int ik1,ik2,im,e;
    
    for (e=0;e<ensembles;e++){
            dJ[e].KM=(double*) malloc(sizeof(double)*head[e].nk*head[e].nk);
            for(ik1=0;ik1<2;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[e].nk;ik2++){   
                im=mass_index[e][ik2][ik1];
                L=jack_files[e].a*head[e].l1;
                M=dJ[e].M_PS_GEVP_jack[im][jack_files[e].Njack-1]*197.326963/jack_files[e].a;
                
                RM=exp(3.81729-0.0130342*M-2.1714*L);
                dJ[e].KM[im]=RM+1;
                if (ik1==0 && ik2==0)printf("%f   %f   im=%d\n",dJ[e].M_PS_GEVP_jack[im][jack_files[e].Njack-1],dJ[e].KM[im], im );
            }}
            
    }
    
}

void Kf_FSE(struct database_file_jack *jack_files, struct header *head, struct data_jack *dJ){
    double Rf,L,M;
    int ik1,ik2,im,e;
    
    for (e=0;e<ensembles;e++){
            dJ[e].Kf=(double*) malloc(sizeof(double)*head[e].nk*head[e].nk);
            for(ik1=0;ik1<2;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
            for(ik2=ik1;ik2<head[e].nk;ik2++){   
                im=mass_index[e][ik2][ik1];
                L=jack_files[e].a*head[e].l1;
                M=dJ[e].M_PS_GEVP_jack[im][jack_files[e].Njack-1]*197.326963/jack_files[e].a;
                
                Rf=exp(4.58982-0.0138032*M-2.013*L);
                dJ[e].Kf[im]=-Rf+1;
            }}
    }
    
}


int main(int argc, char **argv){
    
    int i,j;
    struct header  *head;
    struct database_file_jack  *jack_files;
    //double ***M_PS_GEVP_jack,***f_PS_jack;
    double ***M_PS_GEVP_jack_tot,***f_PS_jack_tot;
   
    double *tmp,**fit,*tmp1;
    
   error(argc!=2,1,"main ",
         "usage: ./fit_all_beta   jack/boot");

   error(strcmp(argv[1],"jack")!=0 && strcmp(argv[1],"boot")!=0 ,2,"main ",
         "choose jack or boot \n usage: ./fit_all_beta   jack/boot");
 
    head=(struct header*) malloc(sizeof(struct header)*ensembles);
    jack_files=(struct database_file_jack *) malloc (sizeof(struct database_file_jack )*ensembles);
    if( strcmp(argv[1],"jack")==0){
                mysprintf(jack_files[0].sampling,NAMESIZE,"jack");
    }
    if( strcmp(argv[1],"boot")==0){
                mysprintf(jack_files[0].sampling,NAMESIZE,"boot");
    }
    

//    M_PS_GEVP_jack=(double***) malloc(sizeof(double**)*ensembles);
 //   f_PS_jack=(double***) malloc(sizeof(double**)*ensembles);
//    alloca_data_jack(dataJ);
    dataJ=(struct data_jack *) malloc(sizeof(struct data_jack)*ensembles);

    files_declarations(argv,jack_files,head);
    mass_index=init_mass_index_ave_r(head);
   // read_files_jack(jack_files,head,mass_index,M_PS_GEVP_jack,f_PS_jack);
    read_files_jack(jack_files,head,mass_index,dataJ);
    KM_FSE(jack_files, head, dataJ);
    Kf_FSE(jack_files, head, dataJ);
    
    int atimesObs,e;
    atimesObs=2*(2+1);
    //lattice_spacings*(N )
   //N=obs+ prior
   //obs=2, i.e. M and f
    index_a=(int**) malloc(sizeof(int*)*(atimesObs)); 
    for(i=0;i<atimesObs;i++)
       index_a[i]=(int*) malloc(sizeof(int)*ensembles); 
   
   for (e=0;e<ensembles;e++){
        index_a[0][e]=e;  
        index_a[1][e]=e+3;
        index_a[2][e]=e;
        index_a[3][e]=e+3;
        index_a[4][e]=0;
        index_a[5][e]=3;
   }

    
    
   printf("ensembles\n");
    for (i=0;i<ensembles;i++){
        
        printf("L%dT%d N=%d  musea=%.5f  KM=%f  Kf=%f ",head[i].l1,head[i].l0, jack_files[i].Njack-1 , head[i].musea, dataJ[i].KM[0], dataJ[i].Kf[0]); 
        tmp=mean_and_error(argv[1],jack_files[i].Njack,   dataJ[i].M_PS_jack[0] );
        printf("M_PS=%g  +-  %g  ",  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_files[i].Njack,   dataJ[i].f_PS_jack[0] );
        printf("  f_PS=%g  +-  %g  \t ",  tmp[0],tmp[1] );
        free(tmp);
        
        tmp1=(double*) malloc(sizeof(double)* jack_files[i].Njack);
        for(j=0;j<jack_files[i].Njack;j++)
            tmp1[j]= dataJ[i].M_PS_jack[0][j]*dataJ[i].M_PS_jack[0][j]/( dataJ[i].f_PS_jack[0][j]* dataJ[i].f_PS_jack[0][j]);
        tmp=mean_and_error(argv[1],jack_files[i].Njack,   tmp1 );
        printf("  M_PS^2/f_PS^2=%g  +-  %g  \n ",  tmp[0],tmp[1] );
        free(tmp);free(tmp1);
        
    }
  /*      printf("ensembles\n");

    for (i=0;i<ensembles;i++){
      
        tmp=mean_and_error_jack_biased1(jack_files[i].Njack,   dataJ[i].M_PS_GEVP_jack[0] );
        printf("L%dT%d N=%d  musea=%.5f  M_PS=%g  +-  %g  ",head[i].l1,head[i].l0, jack_files[i].Njack-1 , head[i].musea,  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error_jack_biased1(jack_files[i].Njack,   dataJ[i].f_PS_jack[0] );
        printf("  f_PS=%g  +-  %g   \t",  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error_jack_biased1(jack_files[i].Njack,   dataJ[i].Zp );
        printf("  Zp=%g  +-  %g    %s\n",  tmp[0],tmp[1]  ,jack_files[i].f_PS);
        free(tmp);
        
    }
       printf("ensembles\n");

    for (i=0;i<ensembles;i++){
      
        tmp=mean_and_error_jack_biased(jack_files[i].Njack,   dataJ[i].M_PS_GEVP_jack[0] );
        printf("L%dT%d N=%d  musea=%.5f  M_PS=%g  +-  %g  ",head[i].l1,head[i].l0, jack_files[i].Njack-1 , head[i].musea,  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error_jack_biased(jack_files[i].Njack,   dataJ[i].f_PS_jack[0] );
        printf("  f_PS=%g  +-  %g   \t",  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error_jack_biased(jack_files[i].Njack,   dataJ[i].Zp );
        printf("  Zp=%g  +-  %g    %s\n",  tmp[0],tmp[1]  ,jack_files[i].f_PS);
        free(tmp);
        
    }*/

    
    
    //M_PS_GEVP_jack_tot=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, dataJ.M_PS_GEVP_jack);//every time jack_tot is initialised
    //f_PS_jack_tot=create_generalised_jack( jack_files, head, &jack_tot,mass_index, dataJ.f_PS_jack);
    if( strcmp(argv[1],"jack")==0){
                gjack=create_generalised_jack( jack_files, head, &jack_tot ,mass_index, dataJ);
                mysprintf(jack_files[0].sampling,NAMESIZE,"jack");
    }
    if( strcmp(argv[1],"boot")==0){
                gjack=create_generalised_boot( jack_files, head, &jack_tot ,mass_index, dataJ);
                mysprintf(jack_files[0].sampling,NAMESIZE,"boot");
    }

    
    int im,ik1,ik2;
   for(ik1=0;ik1<1;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
   for(ik2=ik1;ik2<4;ik2++){
    printf("ensambles after generalised jack ik2=%d  ik1=%d\n",ik2,ik1);

    for (i=0;i<ensembles;i++){
        im=mass_index[i][ik2][ik1];
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].M_PS_jack[im] );
        printf("L%dT%d  mu2=%.5f mu1=%.5f  M_PS=%g  +-  %g  ",head[i].l1,head[i].l0, head[i].k[head[i].nk+ik2], head[i].k[head[i].nk+ik1] ,  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].f_PS_jack[im] );
        printf("  f_PS=%g  +-  %g   \t",  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].Zp );
        printf("  Zp=%g  +-  %g    %s\n",  tmp[0],tmp[1]  ,jack_files[i].f_PS);
        free(tmp);
        
    }
    printf("\n");
    }}
    for(ik1=0;ik1<1;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
   for(ik2=4;ik2<7;ik2++){
    printf("ensambles after generalised jack ik2=%d  ik1=%d\n",ik2,ik1);
    for (i=0;i<ensembles;i++){
        im=mass_index[i][ik2][ik1];
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].M_PS_GEVP_jack[im] );
        printf("L%dT%d  mu2=%.5f mu1=%.5f  M_PS=%g  +-  %g  ",head[i].l1,head[i].l0, head[i].k[head[i].nk+ik2], head[i].k[head[i].nk+ik1] ,  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].f_PS_ls_ss_jack[im] );
        printf("  f_PS=%g  +-  %g   \t",  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].Zp );
        printf("  Zp=%g  +-  %g    %s\n",  tmp[0],tmp[1]  ,jack_files[i].f_PS);
        free(tmp);
        
    }
     printf("\n");
    }}
    for(ik1=1;ik1<2;ik1++){     //for(ik1=0;ik1<=ik2;ik1++){
   for(ik2=4;ik2<7;ik2++){
    printf("ensambles after generalised jack ik2=%d  ik1=%d\n",ik2,ik1);
    for (i=0;i<ensembles;i++){
        im=mass_index[i][ik2][ik1];
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].M_PS_GEVP_jack[im] );
        printf("L%dT%d  mu2=%.5f mu1=%.5f  M_PS=%g  +-  %g  ",head[i].l1,head[i].l0, head[i].k[head[i].nk+ik2], head[i].k[head[i].nk+ik1] ,  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].f_PS_ls_ss_jack[im] );
        printf("  f_PS=%g  +-  %g   \t",  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error(argv[1],jack_tot,   gjack[i].Zp );
        printf("  Zp=%g  +-  %g    %s\n",  tmp[0],tmp[1]  ,jack_files[i].f_PS);
        free(tmp);
        
    }
    printf("\n");
    }}
    
  /*   printf("ensambles\n");
    for (i=0;i<ensembles;i++){
        tmp=mean_and_error_jack(jack_tot,   gjack[i].M_PS_GEVP_jack[0] );
        printf("L%dT%d  L=%f musea=%.5f  M_PS=%g  +-  %g  ",head[i].l1,head[i].l0, head[i].l1*jack_files[i].a ,head[i].musea,  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error_jack(jack_tot,   gjack[i].f_PS_jack[0] );
        printf("  f_PS=%g  +-  %g   \t",  tmp[0],tmp[1] );
        free(tmp);
        tmp=mean_and_error_jack(jack_tot,   gjack[i].Zp );
        printf("  Zp=%g  +-  %g    %s\n",  tmp[0],tmp[1]  ,jack_files[i].f_PS);
        free(tmp);
        
    }
    */
  


 //  fit=fit_Mpi(jack_files, head , jack_tot, mass_index, gjack );
/*    Ci=(double**) malloc(sizeof(double*)*4);
    Ci[0]=mean_and_error_jack_biased(jack_tot,fit[0]);
    Ci[1]=mean_and_error_jack_biased(jack_tot,fit[1]);
    Ci[2]=mean_and_error_jack_biased(jack_tot,fit[2]);
    Ci[3]=mean_and_error_jack_biased(jack_tot,fit[3]);
    printf("Bw=(%g+-%g)  fw=(%g+-%g)  P1=(%g+-%g)  P2=(%g+-%g) \n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1],Ci[2][0],Ci[2][1],Ci[3][0],Ci[3][1]);

    for (i=0;i<4;i++)
        free(Ci[i]);
    free(Ci);
    
    Ci=(double**) malloc(sizeof(double*)*6);

    
    fit=fit_Mpi_fw(jack_files, head , jack_tot, mass_index, gjack );
    Ci[0]=mean_and_error_jack_biased(jack_tot,fit[0]);
    Ci[1]=mean_and_error_jack_biased(jack_tot,fit[1]);
    Ci[2]=mean_and_error_jack_biased(jack_tot,fit[2]);
    Ci[3]=mean_and_error_jack_biased(jack_tot,fit[3]);
    Ci[4]=mean_and_error_jack_biased(jack_tot,fit[4]);
    Ci[5]=mean_and_error_jack_biased(jack_tot,fit[5]);
    printf("Bw=(%g) ;  fw=(%g) ; P1=(%g);  P2=(%g);  P3=(%g);  P4=(%g); \n",Ci[0][0],Ci[1][0],Ci[2][0],Ci[3][0],Ci[4][0],Ci[5][0]);
    printf("\\begin{gather}\\nn B_0w_0=(%g\\pm%.2g) \\quad fw_0=(%g\\pm%.2g) \\\\\\nn  P1=(%g\\pm%.2g) \\quad P2=(%g\\pm%.2g) \\\\\\nn P3=(%g\\pm%.2g) \\quad P4=(%g\\pm%.2g)  \\end{gather}\n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1],Ci[2][0],Ci[2][1],Ci[3][0],Ci[3][1],Ci[4][0],Ci[4][1],Ci[5][0],Ci[5][1]);
*/
    double *tif=(double*) malloc(sizeof(double)*6),*xi=(double*) malloc(sizeof(double)*jack_tot);
    double *fw_phys=(double*) malloc(sizeof(double)*jack_tot);
    double *phys_point=(double*) malloc(sizeof(double)*2);
    
    double *Mw2=(double*) malloc(sizeof(double)*jack_tot);
    double *B_point=(double*) malloc(sizeof(double)*2);
        double **Ci;

    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////    
/*    printf("\n\nw0 from MDs(Mpi, MK, MD)\n");

    tif=(double*) malloc(sizeof(double)*1);
    Ci=(double**) malloc(sizeof(double*)*1);
    
    
    fit=w0_from_MDs(jack_files, head , jack_tot, mass_index, gjack,  &result );
    for (j=0;j<jack_tot;j++){
        fit[0][j]=fit[0][j]/result.MDsMeV[j];
    }
    Ci[0]=mean_and_error_jack_biased(jack_tot,fit[0]);

    printf("Imposing $M_\\Ds =%.2f \\pm %.2g$ \n",v_MDsMeV,err_MDsMeV);
    
    printf("\\begin{gather}\n   w_{0}=(%g\\pm%.2g) MeV^{-1}              =(%g\\pm%.2g) fm     \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1],Ci[0][0]*197.326963,Ci[0][1]*197.326963);
    
     
    free(tif);
    for (i=0;i<1;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);    
 */
    
   
    
    
  /*
    for (j=0;j<jack_tot;j++){
        for (i=0;i<6;i++)
            tif[i]=fit[i][j];
       xi[j]=rtbis( m_over_f,6,tif, 0.0001, 0.01, 1e-10);
       phys_point[0]=xi[j];
       phys_point[1]=1e+16;
       fw_phys[j]=Mw2_fw(1, 2, phys_point,6,tif)/130.41;
      
       B_point[0]=head[4].k[head[4].nk]*gjack[4].w0[j]/gjack[4].Zp[j];
       B_point[1]=gjack[4].w0[j];
       Mw2[j]=Mw2_fw(0, 2, B_point,6,tif);
      
    }
    printf("Imposing M_pi/f_pi and f_pi\n");
    Ci[0]=mean_and_error_jack_biased(jack_tot,xi);
    printf("m_ud*w0=(%g+-%.2g)\n",Ci[0][0],Ci[0][1]);
    Ci[0]=mean_and_error_jack_biased(jack_tot,fw_phys);
    printf("w0=(%g+-%.2g) Mev^-1   \t  (%g+-%.2g) fm \n",Ci[0][0],Ci[0][1],Ci[0][0]*197.326963,Ci[0][1]*197.326963);
    
     
    for (j=0;j<jack_tot;j++){
        for (i=0;i<6;i++)
            tif[i]=fit[i][j];
       xi[j]=rtbis( Mw2_fw_a0_minus,6,tif, 0.0001, 0.01, 1e-10);
       phys_point[0]=xi[j];
       phys_point[1]=1e+16;
       fw_phys[j]=Mw2_fw(1, 2, phys_point,6,tif)/(w0fm/197.326963);
    }
    printf("Imposing $M_\\pi w_0=0.117244859$ and $w_0=%f$ fm\n",w0fm);
    Ci[0]=mean_and_error_jack_biased(jack_tot,xi);
    printf("\\begin{gather} m_{ud}w_0=(%g\\pm%.2g) \t  m_{ud}=(%g\\pm%.2g) MeV \\nn \\\\ \n",Ci[0][0],Ci[0][1],Ci[0][0]/(w0fm/197.326963),Ci[0][1]/(w0fm/197.326963));
    Ci[0]=mean_and_error_jack_biased(jack_tot,fw_phys);
    printf("f_\\pi=(%g\\pm%.2g) Mev  \\nn \\end{gather}\n",Ci[0][0],Ci[0][1]);
   
    
    
    printf("\n\nchecks\n");
    Ci[0]=mean_and_error_jack_biased(jack_tot,Mw2);
    printf("M*w0(B.072)=(%g+-%g)\n",Ci[0][0],Ci[0][1]);

    for (j=0;j<jack_tot;j++){  
         for (i=0;i<6;i++)
            tif[i]=fit[i][j];
       B_point[0]=0.00322426 ;
       B_point[1]=2.4988;
       Mw2[j]=Mw2_fw(0, 2, B_point,6,tif);
    }
    Ci[0]=mean_and_error_jack_biased(jack_tot,Mw2);
    printf("(M*w0)^2(C.06)=(%g+-%g)\n",Ci[0][0],Ci[0][1]);

    for (j=0;j<jack_tot;j++){  
         for (i=0;i<6;i++)
            tif[i]=fit[i][j];
       B_point[0]=0.00322426 ;
       B_point[1]=2.4988;
       Mw2[j]=Mw2_fw(1, 2, B_point,6,tif);
    }
    Ci[0]=mean_and_error_jack_biased(jack_tot,Mw2);
    printf("f*w0(C.06)=(%g+-%g)\n",Ci[0][0],Ci[0][1]);
    
    
    free(tif);
    for (i=0;i<6;i++)
        free(Ci[i]);
    free(Ci);*/
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2
  
    /*printf("\n\nChiral c2\n");

    tif=(double*) malloc(sizeof(double)*7);
    Ci=(double**) malloc(sizeof(double*)*7);
    
    fit=fit_Mpi_fw_chiral_c2(jack_files, head , jack_tot, mass_index, gjack );
    Ci[0]=mean_and_error_jack_biased(jack_tot,fit[0]);
    Ci[1]=mean_and_error_jack_biased(jack_tot,fit[1]);
    Ci[2]=mean_and_error_jack_biased(jack_tot,fit[2]);
    Ci[3]=mean_and_error_jack_biased(jack_tot,fit[3]);
    Ci[4]=mean_and_error_jack_biased(jack_tot,fit[4]);
    Ci[5]=mean_and_error_jack_biased(jack_tot,fit[5]);
   Ci[6]=mean_and_error_jack_biased(jack_tot,fit[6]);
    printf("Bwc=(%g) ;  fwc=(%g) ; P1c=(%g);  P2c=(%g);  P3c=(%g);  P4c=(%g);  c2wwc=(%g);\n",Ci[0][0],Ci[1][0],Ci[2][0],Ci[3][0],Ci[4][0],Ci[5][0],Ci[6][0]);
    printf("\\begin{gather}\\nn B_0w_0=(%g\\pm%.2g) \\quad fw_0=(%g\\pm%.2g) \\\\\\nn  P1=(%g\\pm%.2g) \\quad P2=(%g\\pm%.2g) \\\\\\nn P3=(%g\\pm%.2g) \\quad P4=(%g\\pm%.2g) \\\\\\nn  c_2w_0^2=(%g\\pm%0.2g) \\end{gather}\n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1],Ci[2][0],Ci[2][1],Ci[3][0],Ci[3][1],Ci[4][0],Ci[4][1],Ci[5][0],Ci[5][1],Ci[6][0],Ci[6][1]);


       for (j=0;j<jack_tot;j++){
        for (i=0;i<7;i++)
            tif[i]=fit[i][j];
       xi[j]=rtbis( m_over_f,6,tif, 0.0001, 0.01, 1e-10);
       phys_point[0]=xi[j];
       phys_point[1]=1e+16;
       fw_phys[j]=Mw2_fw_chiral_c2(1, 2, phys_point,7,tif)/130.41;
      
       B_point[0]=head[4].k[head[4].nk]*gjack[4].w0[j]/gjack[4].Zp[j];
       B_point[1]=gjack[4].w0[j];
       Mw2[j]=Mw2_fw_chiral_c2(0, 2, B_point,7,tif);
      
    }
       printf("Imposing M_pi/f_pi and f_pi\n");
    Ci[0]=mean_and_error_jack_biased(jack_tot,xi);
    printf("m_ud*w0=(%g+-%.2g)\n",Ci[0][0],Ci[0][1]);
    Ci[0]=mean_and_error_jack_biased(jack_tot,fw_phys);
    printf("w0=(%g+-%.2g) Mev^-1   \t  (%g+-%.2g) fm \n",Ci[0][0],Ci[0][1],Ci[0][0]*197.326963,Ci[0][1]*197.326963);

    
    
    
     for (j=0;j<jack_tot;j++){
        for (i=0;i<7;i++)
            tif[i]=fit[i][j];
       xi[j]=rtbis( Mw2_fw_chiral_c2_a0_minus,7,tif, 0.0001, 0.01, 1e-10);
       phys_point[0]=xi[j];
       phys_point[1]=1e+16;
       fw_phys[j]=Mw2_fw(1, 2, phys_point,6,tif)/(w0fm/197.326963);
    }
    printf("Imposing $M_\\pi w_0=0.117244859$ and $w_0=%f$ fm\n",w0fm);
    Ci[0]=mean_and_error_jack_biased(jack_tot,xi);
    printf("\\begin{gather} m_{ud}w_0=(%g\\pm%.2g) \t  m_{ud}=(%g\\pm%.2g) MeV \\nn \\\\ \n",Ci[0][0],Ci[0][1],Ci[0][0]/(w0fm/197.326963),Ci[0][1]/(w0fm/197.326963));
    Ci[0]=mean_and_error_jack_biased(jack_tot,fw_phys);
    printf("f_\\pi=(%g\\pm%.2g) Mev  \\nn \\end{gather}\n",Ci[0][0],Ci[0][1]);
   
     printf("\n\nchecks\n");
     
    Ci[0]=mean_and_error_jack_biased(jack_tot,Mw2);
    printf("M*w0(B.072)=(%g+-%g)\n",Ci[0][0],Ci[0][1]);

    for (j=0;j<jack_tot;j++){  
         for (i=0;i<6;i++)
            tif[i]=fit[i][j];
       B_point[0]=0.00322426 ;
       B_point[1]=2.4988;
       Mw2[j]=Mw2_fw_chiral_c2(0, 2, B_point,7,tif);
    }
    Ci[0]=mean_and_error_jack_biased(jack_tot,Mw2);
    printf("(M*w0)^2(C.06)=(%g+-%g)\n",Ci[0][0],Ci[0][1]);

    for (j=0;j<jack_tot;j++){  
         for (i=0;i<6;i++)
            tif[i]=fit[i][j];
       B_point[0]=0.00322426 ;
       B_point[1]=2.4988;
       Mw2[j]=Mw2_fw_chiral_c2(1, 2, B_point,7,tif);
    }
    Ci[0]=mean_and_error_jack_biased(jack_tot,Mw2);
    printf("f*w0(C.06)=(%g+-%g)\n",Ci[0][0],Ci[0][1]);

    
    
    
    free(tif);
    for (i=0;i<7;i++)
        free(Ci[i]);
    free(Ci);
    */
    
     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2
    printf("\n\nChiral with FVE\n");

    tif=(double*) malloc(sizeof(double)*6);
    double **Ci1=(double**) malloc(sizeof(double*)*6);
    double *w0_estimate=(double*) malloc(sizeof(double)*jack_tot);
    
    fit=fit_Mpi_fw_chiral_FVE(jack_files, head , jack_tot, mass_index, gjack );
    Ci1[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci1[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    Ci1[2]=mean_and_error(argv[1],jack_tot,fit[2]);
    Ci1[3]=mean_and_error(argv[1],jack_tot,fit[3]);
    Ci1[4]=mean_and_error(argv[1],jack_tot,fit[4]);
    Ci1[5]=mean_and_error(argv[1],jack_tot,fit[5]);
    
    result.Bw=(double*) malloc(sizeof(double)*jack_tot);
    result.fw=(double*) malloc(sizeof(double)*jack_tot);
    result.l3b=(double*) malloc(sizeof(double)*jack_tot);
    result.l4b=(double*) malloc(sizeof(double)*jack_tot);
    result.mlw=(double*) malloc(sizeof(double)*jack_tot);
    result.fpiw=(double*) malloc(sizeof(double)*jack_tot);
    for (j=0;j<jack_tot;j++){
        result.Bw[j]=fit[0][j];
        result.fw[j]=fit[1][j];
        result.l3b[j]=fit[2][j];
        result.l4b[j]=fit[4][j];
        
    }
    tmp1=(double*) malloc(sizeof(double)*jack_tot);
    
    printf("Bw=(%g) ;  fw=(%g) ; l3b=(%g);  P2=(%g);  l4b=(%g);  P4=(%g);  \n",Ci1[0][0],Ci1[1][0],Ci1[2][0],Ci1[3][0],Ci1[4][0],Ci1[5][0]);
    printf("\\begin{gather}\\nn B_0w_0=(%g\\pm%.2g) \\quad fw_0=(%g\\pm%.2g) \\\\\\nn  \\bar\\ell_3=(%g\\pm%.2g) \\quad P2=(%g\\pm%.2g) \\\\\\nn \\bar\\ell_4=(%g\\pm%.2g) \\quad P4=(%g\\pm%.2g) \\nn   \\end{gather}\n",Ci1[0][0],Ci1[0][1],Ci1[1][0],Ci1[1][1],Ci1[2][0],Ci1[2][1],Ci1[3][0],Ci1[3][1],Ci1[4][0],Ci1[4][1],Ci1[5][0],Ci1[5][1]);

    
    
    double in;
     for (j=0;j<jack_tot;j++){
        for (i=0;i<6;i++)
            tif[i]=fit[i][j];
       in=result.MpiMeV[j]*result.w0MeV[j]*result.MpiMeV[j]*result.w0MeV[j];
       xi[j]=rtbis( Mw2_chiral_FVE_a0_minus,in,6,tif, 0.0001, 0.01, 1e-10);
       phys_point[0]=xi[j];
       phys_point[1]=1e+16;
       result.mlw[j]=xi[j];
       //fw_phys[j]=Mw2_fw(1, 2, phys_point,6,tif)/(w0fm/197.326963);
       result.fpiw[j]=fPSw_chiral_FVE(  phys_point[0],6,tif);
       fw_phys[j]=result.fpiw[j]/(result.w0MeV[j]);
       xi[j]=xi[j]/result.w0MeV[j];
       
       in=result.MpiMeV[j]*result.MpiMeV[j]/(result.fpiMeV_exp[j]*result.fpiMeV_exp[j]);
       xi[j]=rtbis(Mw2_over_fw2_chiral_FVE_a0_minus,in,6,tif, 0.0001, 0.01, 1e-10);
       phys_point[0]=xi[j];
       result.fpiw[j]=fPSw_chiral_FVE(  phys_point[0],6,tif);
       w0_estimate[j]=result.fpiw[j]/(v_fpiMeV_exp/197.326963);
       xi[j]=xi[j]/(w0_estimate[j]/197.326963);
       
    }
    double **C1=(double**) malloc(sizeof(double*)*2);
   
    printf("Imposing $M_\\pi =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %.2g$ fm\n",v_MpiMeV,err_MpiMeV,v_w0fm,err_w0fm);
    
    C1[0]=mean_and_error(argv[1],jack_tot,result.mlw);
    printf("\\begin{gather}   m_{ud}=(%g\\pm%.2g) MeV \\nn \\\\ \n",C1[0][0],C1[0][1]);
    free(C1[0]);
    
     C1[1]=mean_and_error(argv[1],jack_tot,fw_phys);
    printf("f_\\pi=(%g\\pm%.2g) Mev  \\nn \\end{gather}\n",C1[1][0],C1[1][1]);
    free(tif);free(C1[1]);
    
    printf("Imposing $M_\\pi =%.2f \\pm %.2g$ and $f_pi=%.4f\\pm %.2g$ MeV\n",v_MpiMeV,err_MpiMeV,v_fpiMeV_exp,err_fpiMeV_exp);

    C1[0]=mean_and_error(argv[1],jack_tot,xi);
    printf("\\begin{gather}\n   m_{ud}=(%g\\pm%.2g) MeV \\nn \\\\ \n",C1[0][0],C1[0][1]);
    free(C1[0]);

    C1[1]=mean_and_error(argv[1],jack_tot,w0_estimate);
    printf("w_0=(%g\\pm%.2g) fm (from f_pi)  \\nn \n\\end{gather}\n",C1[1][0],C1[1][1]);
    free(C1[1]);
    
    for (j=0;j<jack_tot;j++)
        xi[j]=result.Bw[j]/(w0_estimate[j]/197.326963);
    C1[0]=mean_and_error(argv[1],jack_tot,xi);
    printf("   B=(%g\\pm%.2g) MeV \\nn \\\\ \n",C1[0][0],C1[0][1]);
    free(C1[0]);
    
    for (j=0;j<jack_tot;j++)
        xi[j]=result.fw[j]/(w0_estimate[j]/197.326963);
     C1[0]=mean_and_error(argv[1],jack_tot,xi);
    printf("   f=(%g\\pm%.2g) MeV \\nn \\\\ \n",C1[0][0],C1[0][1]);
    free(C1[0]);
    
       
        
    
   
    
    
    
    ////////////////////////////////////////////////////////////////////////////
     printf("\n\nChiral with FVE P1 and P3\n");
    tif=(double*) malloc(sizeof(double)*6);
    Ci=(double**) malloc(sizeof(double*)*6);

    fit=fit_Mpi_fw_chiral_FVE_P1_P3(jack_files, head , jack_tot, mass_index, gjack );
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    Ci[2]=mean_and_error(argv[1],jack_tot,fit[2]);
    Ci[3]=mean_and_error(argv[1],jack_tot,fit[3]);
    Ci[4]=mean_and_error(argv[1],jack_tot,fit[4]);
    Ci[5]=mean_and_error(argv[1],jack_tot,fit[5]);
    
    printf("Bw=(%g) ;  fw=(%g) ; P1=(%g);  P2=(%g);  P3=(%g);  P4=(%g);  \n",Ci[0][0],Ci[1][0],Ci[2][0],Ci[3][0],Ci[4][0],Ci[5][0]);
    printf("\\begin{gather}\\nn\n B_0w_0=(%g\\pm%.2g) \\quad\n fw_0=(%g\\pm%.2g) \\\\\\nn\n  P_1=(%g\\pm%.2g) \\quad\n P2=(%g\\pm%.2g) \\\\\\nn\n P_3=(%g\\pm%.2g) \\quad\n P4=(%g\\pm%.2g)\n \\nn\n   \\end{gather}\n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1],Ci[2][0],Ci[2][1],Ci[3][0],Ci[3][1],Ci[4][0],Ci[4][1],Ci[5][0],Ci[5][1]);

    
    
     for (j=0;j<jack_tot;j++){
        for (i=0;i<6;i++)
            tif[i]=fit[i][j];
       in=result.MpiMeV[j]*result.w0MeV[j]*result.MpiMeV[j]*result.w0MeV[j];
       xi[j]=rtbis( Mw2_chiral_FVE_a0_minus_P1_P3,in,6,tif, 0.0001, 0.01, 1e-10);
       phys_point[0]=xi[j];
       phys_point[1]=1e+16;
       //result.mlw[j]=xi[j];
       //fw_phys[j]=Mw2_fw(1, 2, phys_point,6,tif)/(w0fm/197.326963);
       fw_phys[j]=fPSw_chiral_FVE_P1_P3(  phys_point[0],6,tif)/(result.w0MeV[j]);
       xi[j]=xi[j]/result.w0MeV[j];
       
    }
   free(Ci[0]); free(Ci[1]);
    printf("Imposing $M_\\pi =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %.2g$ fm\n",v_MpiMeV,err_MpiMeV,v_w0fm,err_w0fm);
    Ci[0]=mean_and_error(argv[1],jack_tot,xi);
    printf("\\begin{gather}   m_{ud}=(%g\\pm%.2g) MeV \\nn \\\\ \n",Ci[0][0],Ci[0][1]);
    
    Ci[1]=mean_and_error(argv[1],jack_tot,fw_phys);
    printf("f_\\pi=(%g\\pm%.2g) Mev  \\nn \\end{gather}\n",Ci[1][0],Ci[1][1]);
    free(tif);
    
    
    for (i=0;i<6;i++)
    {    free(Ci[i]); free(fit[i]);}
    free(Ci); free(fit);
    
    
    
    
    
    
 /*   for (i=0;i<6;i++)
        free(Ci1[i]);
    free(Ci1);
   */ 
 
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2
 /*   printf("\n\nChiral with FVE prior\n");

    tif=(double*) malloc(sizeof(double)*8);
    Ci=(double**) malloc(sizeof(double*)*8);
    
    fit=fit_Mpi_fw_chiral_FVE_prior(jack_files, head , jack_tot, mass_index, gjack );
    printf("HERE\n");
    printf("fit %f\n",fit[6][0]);
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    Ci[2]=mean_and_error(argv[1],jack_tot,fit[2]);
    Ci[3]=mean_and_error(argv[1],jack_tot,fit[3]);
    Ci[4]=mean_and_error(argv[1],jack_tot,fit[4]);
    Ci[5]=mean_and_error(argv[1],jack_tot,fit[5]);
    Ci[6]=mean_and_error(argv[1],jack_tot,fit[6]);
    Ci[7]=mean_and_error(argv[1],jack_tot,fit[7]);
    
    printf("Bw=(%g) ;  fw=(%g) ; l3b=(%g);  P2=(%g);  l4b=(%g);  P4=(%g);  ZPA=%g  ZPB=%g\n",Ci[0][0],Ci[1][0],Ci[2][0],Ci[3][0],Ci[4][0],Ci[5][0],Ci[6][0],Ci[7][0]);
    printf("\\begin{gather}\\nn B_0w_0=(%g\\pm%.2g) \\quad fw_0=(%g\\pm%.2g) \\\\\\nn  \\bar\\ell_3=(%g\\pm%.2g) \\quad P2=(%g\\pm%.2g) \\\\\\nn \\bar\\ell_4=(%g\\pm%.2g) \\quad P4=(%g\\pm%.2g)      ZPA=(%g\\pm%.2g)   \\ ZPA=(%g\\pm%.2g)\\nn   \\end{gather}\n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1],Ci[2][0],Ci[2][1],Ci[3][0],Ci[3][1],Ci[4][0],Ci[4][1],Ci[5][0],Ci[5][1],Ci[6][0],Ci[6][1],Ci[7][0],Ci[7][1]);

    
    
     for (j=0;j<jack_tot;j++){
        for (i=0;i<6;i++)
            tif[i]=fit[i][j];
       in=result.MpiMeV[j]*result.w0MeV[j]*result.MpiMeV[j]*result.w0MeV[j];
       xi[j]=rtbis( Mw2_chiral_FVE_a0_minus,in,6,tif, 0.0001, 0.01, 1e-10);
       phys_point[0]=xi[j];
       phys_point[1]=1e+16;
       result.mlw[j]=xi[j];
       //fw_phys[j]=Mw2_fw(1, 2, phys_point,6,tif)/(w0fm/197.326963);
       result.fpiw[j]=fPSw_chiral_FVE(  phys_point[0],6,tif);
       fw_phys[j]=result.fpiw[j]/(result.w0MeV[j]);
       xi[j]=xi[j]/result.w0MeV[j];
       
    }
    double **C1p=(double**) malloc(sizeof(double*)*2);
   
    printf("Imposing $M_\\pi =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %.2g$ fm\n",v_MpiMeV,err_MpiMeV,v_w0fm,err_w0fm);
    C1p[0]=mean_and_error(argv[1],jack_tot,xi);
    printf("\\begin{gather}   m_{ud}=(%g\\pm%.2g) MeV \\nn \\\\ \n",C1p[0][0],C1p[0][1]);
    
    C1p[1]=mean_and_error(argv[1],jack_tot,fw_phys);
    printf("f_\\pi=(%g\\pm%.2g) Mev  \\nn \\end{gather}\n",C1p[1][0],C1p[1][1]);
    free(tif);
   */ 
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2
       printf("\n\nChiral with FVE l3b and l4b from flag\n");

    tif=(double*) malloc(sizeof(double)*4);
    Ci=(double**) malloc(sizeof(double*)*4);
    
    fit=fit_Mpi_fw_chiral_FVE_flag(jack_files, head , jack_tot, mass_index, gjack );
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    Ci[2]=mean_and_error(argv[1],jack_tot,fit[2]);
    Ci[3]=mean_and_error(argv[1],jack_tot,fit[3]);
   
    printf("Bw=(%g) ;  fw=(%g) ;   P2=(%g);    P4=(%g);  \n",Ci[0][0],Ci[1][0],Ci[2][0],Ci[3][0]);
    printf("\\begin{gather}\\nn B_0w_0=(%g\\pm%.2g) \\quad fw_0=(%g\\pm%.2g) \\\\\\nn   \\quad P2=(%g\\pm%.2g)  \\quad P4=(%g\\pm%.2g) \\nn   \\end{gather}\n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1],Ci[2][0],Ci[2][1],Ci[3][0],Ci[3][1]);

    
    
    
     for (j=0;j<jack_tot;j++){
        for (i=0;i<4;i++)
            tif[i]=fit[i][j];
       in=result.MpiMeV[j]*result.w0MeV[j]*result.MpiMeV[j]*result.w0MeV[j];
       xi[j]=rtbis( Mw2_chiral_FVE_a0_minus,in,4,tif, 0.0001, 0.01, 1e-10);
       phys_point[0]=xi[j];
       phys_point[1]=1e+16;
       //fw_phys[j]=Mw2_fw(1, 2, phys_point,6,tif)/(w0fm/197.326963);
       fw_phys[j]=fPSw_chiral_FVE(  phys_point[0],4,tif)/(result.w0MeV[j]);
       xi[j]=xi[j]/result.w0MeV[j];
    }
    printf("Imposing $M_\\pi =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %.2g$ fm\n",v_MpiMeV,err_MpiMeV,v_w0fm,err_w0fm);
    Ci[0]=mean_and_error(argv[1],jack_tot,xi);
    printf("\\begin{gather}  m_{ud}=(%g\\pm%.2g) MeV \\nn \\\\ \n",Ci[0][0],Ci[0][1]);
    Ci[0]=mean_and_error(argv[1],jack_tot,fw_phys);
    printf("f_\\pi=(%g\\pm%.2g) Mev  \\nn \\end{gather}\n",Ci[0][0],Ci[0][1]);
    free(tif);
      for (i=0;i<4;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    
          
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2
       printf("\n\nf from M^2\n");
    int Npar=3;
    tif=(double*) malloc(sizeof(double)*Npar);
    Ci=(double**) malloc(sizeof(double*)*Npar);
    
    fit=fit_fw_of_Mw_chiral_FVE(jack_files, head , jack_tot, mass_index, gjack );
   
    for(i=0;i<Npar;i++){
            Ci[i]=mean_and_error(argv[1],jack_tot,fit[i]);
    }
    /*
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    Ci[2]=mean_and_error(argv[1],jack_tot,fit[2]);
    Ci[3]=mean_and_error(argv[1],jack_tot,fit[3]);
    Ci[4]=mean_and_error(argv[1],jack_tot,fit[4]);
   */
   
        
    printf(" fw=(%g) ;  l4b=(%g);  P4=(%g);   \n",Ci[0][0],Ci[1][0],Ci[2][0]);
    printf("\\begin{gather}\\nn fw_0=(%g\\pm%.2g)  \\quad    \\bar\\ell_4=(%g\\pm%.2g) \\quad P4=(%g\\pm%.2g)    \\end{gather}\n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1],Ci[2][0],Ci[2][1]); 
  //  printf(" fw=(%g) ;  l4b=(%g);  P4=(%g); c2ww=(%g);   \n",Ci[0][0],Ci[1][0],Ci[2][0],Ci[3][0]);
   // printf("\\begin{gather}\\nn fw_0=(%g\\pm%.2g)  \\quad    \\bar\\ell_4=(%g\\pm%.2g) \\quad P4=(%g\\pm%.2g)  \\quad c_2w_0^2=(%g\\pm%.2g) \\nn   \\end{gather}\n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1],Ci[2][0],Ci[2][1],Ci[3][0],Ci[3][1]);
    //printf(" fw=(%g) ;    P4=(%g); c2ww=(%g);  c3=%g; \n",Ci[0][0],Ci[1][0],Ci[2][0],Ci[3][0]);
//printf("\\begin{gather}\\nn fw_0=(%g\\pm%.2g)  \\quad     P4=(%g\\pm%.2g)  \\quad c_2w_0^2=(%g\\pm%.2g) \\nn c3=(%g\\pm%.2g)  \\end{gather}\n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1],Ci[2][0],Ci[2][1],Ci[3][0],Ci[3][1]);

     result.fpiw_from_M=(double*) malloc(sizeof(double)*jack_tot);
     
     
    result.fw_from_M=(double*) malloc(sizeof(double)*jack_tot);
    result.l4b_from_M=(double*) malloc(sizeof(double)*jack_tot);
    for (j=0;j<jack_tot;j++){
        result.fw_from_M[j]=fit[0][j];
        result.l4b_from_M[j]=fit[1][j];
       
        

    }
     
     for (j=0;j<jack_tot;j++){
         for (i=0;i<Npar;i++)
            tif[i]=fit[i][j];
   // computing fpi from w0
        in=result.MpiMeV[j]*result.w0MeV[j]*result.MpiMeV[j]*result.w0MeV[j];
        fw_phys[j]=fw_of_Mw2_physical_point(in,3,tif)/result.w0MeV[j];
        result.fpiw_from_M[j]=fw_of_Mw2_physical_point(in,3,tif);
        
        // computing w0 from fpi
       
        xi[j]=rtbis( fw_over_Mw_of_Mw_minus_exp,  result.fpiMeV_exp[j]/result.MpiMeV[j]   ,Npar,tif, 0.01, 0.2, 1e-10);
        in=xi[j]*xi[j];
        w0_estimate[j]=fw_of_Mw2_physical_point(in,3,tif)/(result.fpiMeV_exp[j]/197.326963);
    }
    printf("Imposing $M_\\pi =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %.2g$ fm\n",v_MpiMeV,err_MpiMeV,v_w0fm,err_w0fm);
   // Ci[0]=mean_and_error(argv[1],jack_tot,xi);
   // printf("\\begin{gather} m_{ud}w_0=(%g\\pm%.2g) \t  m_{ud}=(%g\\pm%.2g) MeV \\nn \\\\ \n",Ci[0][0],Ci[0][1],Ci[0][0]/(w0fm/197.326963),Ci[0][1]/(w0fm/197.326963));
    double **C2=(double**)  malloc(sizeof(double*)*1);
    
    for (j=0;j<jack_tot;j++){
        tmp1[j]=fit[0][j]/result.w0MeV[j];
    }
    C2[0]=mean_and_error(argv[1],jack_tot,tmp1);
    printf("f=(%g\\pm%.2g) MeV  \\nn \\end{gather}\n",C2[0][0],C2[0][1]);
    free(C2[0]);
    
    C2[0]=mean_and_error(argv[1],jack_tot,fw_phys);
    printf("f_\\pi=(%g\\pm%.2g) Mev  \\nn \\end{gather}\n",C2[0][0],C2[0][1]);
    free(C2[0]);
    
    printf("Imposing $M_\\pi =%.2f \\pm %.2g$ and $f_{pi}=%.4f\\pm %.2g$ fm\n",v_MpiMeV,err_MpiMeV,v_fpiMeV_exp,err_fpiMeV_exp);
    C2[0]=mean_and_error(argv[1],jack_tot,w0_estimate);
    printf("w_0=(%g\\pm%.2g) fm  \\nn \\end{gather}\n",C2[0][0],C2[0][1]);
    free(C2[0]);
    
    for (j=0;j<jack_tot;j++){
        tmp1[j]=fit[0][j]*197.326963/w0_estimate[j];
    }
    C2[0]=mean_and_error(argv[1],jack_tot,tmp1);
    printf("f=(%g\\pm%.2g) MeV  \\nn \\end{gather}\n",C2[0][0],C2[0][1]);
    free(C2[0]);
    
    
      printf("\n\\hline\n");
      printf("$B_0w_0$ & $%g\\pm%.2g$  &  -\\\\\n",Ci1[0][0],Ci1[0][1]);
      printf("$f_0w_0$ & $%g\\pm%.2g$  &  $%g\\pm%.2g$\\\\\n",Ci1[1][0],Ci1[1][1],Ci[0][0],Ci[0][1]);
      printf("$\\bar\\ell_3$ & $%g\\pm%.2g$  & - \\\\ \n",Ci1[2][0],Ci1[2][1]);
      printf("$P_2$ & $%g\\pm%.2g$  &  -\\\\\n",Ci1[3][0],Ci1[3][1]);
      printf("$\\bar\\ell_4$ & $%g\\pm%.2g$ &  $%g\\pm%.2g$\\\\\n",Ci1[4][0],Ci1[4][1],Ci[1][0],Ci[1][1]);
      printf("$P_4$ & $%g\\pm%.2g$  &  $%g\\pm%.2g$\\\\\n",Ci1[5][0],Ci1[5][1],Ci[2][0],Ci[2][1]);

      printf("\n\n");
      //printf("$m_\\ell $ & $%g\\pm%.2g$  &  -  \\\\\n",C1[0][0],C1[0][1]);
      //printf("$f_{\\pi}$ & $%g\\pm%.2g$  &  $%g\\pm%.2g$ \\\\\n",C1[1][0],C1[1][1],C2[0][0],C2[0][1]);
      free(tif);
      for (i=0;i<Npar;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2
      printf("\n\nChiral with FVE  P4=00\n");
       
    for (i=0;i<6;i++)
        free(Ci1[i]);
    free(Ci1);
    tif=(double*) malloc(sizeof(double)*5);
     Ci1=(double**) malloc(sizeof(double*)*5);
    
    fit=fit_Mpi_fw_chiral_FVE_P40(jack_files, head , jack_tot, mass_index, gjack );
    Ci1[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci1[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    Ci1[2]=mean_and_error(argv[1],jack_tot,fit[2]);
    Ci1[3]=mean_and_error(argv[1],jack_tot,fit[3]);
    Ci1[4]=mean_and_error(argv[1],jack_tot,fit[4]);
    
    /*
    for (j=0;j<jack_tot;j++){
        result.Bw[j]=fit[0][j];
        result.fw[j]=fit[1][j];
        result.l3b[j]=fit[2][j];
        result.l4b[j]=fit[4][j];
        
    }*/
    
    printf("Bw=(%g) ;  fw=(%g) ; l3b=(%g);  P2=(%g);  l4b=(%g); P4=0   \n",Ci1[0][0],Ci1[1][0],Ci1[2][0],Ci1[3][0],Ci1[4][0]);
    printf("\\begin{gather}\\nn B_0w_0=(%g\\pm%.2g) \\quad fw_0=(%g\\pm%.2g) \\\\\\nn  \\bar\\ell_3=(%g\\pm%.2g) \\quad P2=(%g\\pm%.2g) \\\\\\nn \\bar\\ell_4=(%g\\pm%.2g)  \\nn   \\end{gather}\n",Ci1[0][0],Ci1[0][1],Ci1[1][0],Ci1[1][1],Ci1[2][0],Ci1[2][1],Ci1[3][0],Ci1[3][1],Ci1[4][0],Ci1[4][1]);

    
    
    
     for (j=0;j<jack_tot;j++){
       for (i=0;i<5;i++)
            tif[i]=fit[i][j];
       in=result.MpiMeV[j]*result.w0MeV[j]*result.MpiMeV[j]*result.w0MeV[j];
       xi[j]=rtbis( Mw2_chiral_FVE_a0_minus,in,6,tif, 0.0001, 0.01, 1e-10);//keep the number of parameter 6 or 4 if you want l34b from flag
       phys_point[0]=xi[j];
       phys_point[1]=1e+16;
       //result.mlw[j]=xi[j];
       //fw_phys[j]=Mw2_fw(1, 2, phys_point,6,tif)/(w0fm/197.326963);
       //result.fpiw[j]=fPSw_chiral_FVE(  phys_point[0],6,tif);//keep the number of parameter 6 or 4 if you want l34b from flag
       fw_phys[j]=fPSw_chiral_FVE(  phys_point[0],6,tif)/(result.w0MeV[j]);
       xi[j]=xi[j]/result.w0MeV[j];
       
    }
   
    printf("Imposing $M_\\pi =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %.2g$ fm\n",v_MpiMeV,err_MpiMeV,v_w0fm,err_w0fm);
    C1[0]=mean_and_error(argv[1],jack_tot,xi);
    printf("\\begin{gather}   m_{ud}=(%g\\pm%.2g) MeV \\nn \\\\ \n",C1[0][0],C1[0][1]);
  
    C1[1]=mean_and_error(argv[1],jack_tot,fw_phys);
    printf("f_\\pi=(%g\\pm%.2g) Mev  \\nn \\end{gather}\n",C1[1][0],C1[1][1]);
    free(tif);
 /*   for (i=0;i<6;i++)
        free(Ci1[i]);
    free(Ci1);

   */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2
 /*   printf("\n\nChiral with FVE P4=0 prior\n");

    tif=(double*) malloc(sizeof(double)*7);
    Ci=(double**) malloc(sizeof(double*)*7);
    
    fit=fit_Mpi_fw_chiral_FVE_P40_prior(jack_files, head , jack_tot, mass_index, gjack );
    printf("HERE\n");
    printf("fit %f\n",fit[6][0]);
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    Ci[2]=mean_and_error(argv[1],jack_tot,fit[2]);
    Ci[3]=mean_and_error(argv[1],jack_tot,fit[3]);
    Ci[4]=mean_and_error(argv[1],jack_tot,fit[4]);
    Ci[5]=mean_and_error(argv[1],jack_tot,fit[5]);
    Ci[6]=mean_and_error(argv[1],jack_tot,fit[6]);
    
    
    printf("Bw=(%g) ;  fw=(%g) ; l3b=(%g);  P2=(%g);  l4b=(%g);  P4=(0);  ZPA=%g  ZPB=%g\n",Ci[0][0],Ci[1][0],Ci[2][0],Ci[3][0],Ci[4][0],Ci[5][0],Ci[6][0]);
    printf("\\begin{gather}\\nn B_0w_0=(%g\\pm%.2g) \\quad fw_0=(%g\\pm%.2g) \\\\\\nn  \\bar\\ell_3=(%g\\pm%.2g) \\quad P2=(%g\\pm%.2g) \\\\\\nn \\bar\\ell_4=(%g\\pm%.2g) \\quad P4=(0)      ZPA=(%g\\pm%.2g)   \\ ZPA=(%g\\pm%.2g)\\nn   \\end{gather}\n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1],Ci[2][0],Ci[2][1],Ci[3][0],Ci[3][1],Ci[4][0],Ci[4][1],Ci[5][0],Ci[5][1],Ci[6][0],Ci[6][1]);

    
    
     
     for (j=0;j<jack_tot;j++){
       for (i=0;i<7;i++)
            tif[i]=fit[i][j];
       in=result.MpiMeV[j]*result.w0MeV[j]*result.MpiMeV[j]*result.w0MeV[j];
       xi[j]=rtbis( Mw2_chiral_FVE_a0_minus,in,6,tif, 0.0001, 0.01, 1e-10);//keep the number of parameter 6 or 4 if you want l34b from flag
       phys_point[0]=xi[j];
       phys_point[1]=1e+16;
       //result.mlw[j]=xi[j];   add these two line if you want to pass to the K analysis this values
       //result.fpiw[j]=fPSw_chiral_FVE(  phys_point[0],6,tif);//keep the number of parameter 6 or 4 if you want l34b from flag
       fw_phys[j]=fPSw_chiral_FVE(  phys_point[0],6,tif)/(result.w0MeV[j]);
       xi[j]=xi[j]/result.w0MeV[j];
       
    }
   
    double **C1p0=(double**) malloc(sizeof(double*)*2);
   
    printf("Imposing $M_\\pi =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %.2g$ fm\n",v_MpiMeV,err_MpiMeV,v_w0fm,err_w0fm);
    C1p0[0]=mean_and_error(argv[1],jack_tot,xi);
    printf("\\begin{gather}   m_{ud}=(%g\\pm%.2g) MeV \\nn \\\\ \n",C1p0[0][0],C1p0[0][1]);
    
    C1p0[1]=mean_and_error(argv[1],jack_tot,fw_phys);
    printf("f_\\pi=(%g\\pm%.2g) Mev  \\nn \\end{gather}\n",C1p0[1][0],C1p0[1][1]);
    free(tif);
    
    */
     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2
       printf("\n\nf of M^2 P40\n");

    tif=(double*) malloc(sizeof(double)*2);
    Ci=(double**) malloc(sizeof(double*)*2);
    
    fit=fit_fw_of_Mw_chiral_FVE_P40(jack_files, head , jack_tot, mass_index, gjack );
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
   
  
    printf(" fw=(%g) ;  l4b=(%g);  P4=(0);  \n",Ci[0][0],Ci[1][0]);
    printf("\\begin{gather}\\nn fw_0=(%g\\pm%.2g)  \\quad    \\bar\\ell_4=(%g\\pm%.2g) \\quad P4=0 \\nn   \\end{gather}\n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1]);

    
    
    
     for (j=0;j<jack_tot;j++){
        for (i=0;i<2;i++)
            tif[i]=fit[i][j];
      in=result.MpiMeV[j]*result.w0MeV[j]*result.MpiMeV[j]*result.w0MeV[j];

       fw_phys[j]=fw_of_Mw2_physical_point(in,3,tif)/(result.w0MeV[j]);

    }
    printf("Imposing $M_\\pi =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %.2g$ fm\n",v_MpiMeV,err_MpiMeV,v_w0fm,err_w0fm);
   // Ci[0]=mean_and_error(argv[1],jack_tot,xi);
   // printf("\\begin{gather} m_{ud}w_0=(%g\\pm%.2g) \t  m_{ud}=(%g\\pm%.2g) MeV \\nn \\\\ \n",Ci[0][0],Ci[0][1],Ci[0][0]/(w0fm/197.326963),Ci[0][1]/(w0fm/197.326963));
    C2=(double**)  malloc(sizeof(double*)*1);
    C2[0]=mean_and_error(argv[1],jack_tot,fw_phys);
    printf("f_\\pi=(%g\\pm%.2g) Mev  \\nn \\end{gather}\n",C2[0][0],C2[0][1]);
  
      printf("\n\\hline\n");
      printf("$B_0w_0$ & $%g\\pm%.2g$  &  -\\\\\n",Ci1[0][0],Ci1[0][1]);
      printf("$f_0w_0$ & $%g\\pm%.2g$  &  $%g\\pm%.2g$\\\\\n",Ci1[1][0],Ci1[1][1],Ci[0][0],Ci[0][1]);
      printf("$\\bar\\ell_3$ & $%g\\pm%.2g$  & - \\\\ \n",Ci1[2][0],Ci1[2][1]);
      printf("$P_2$ & $%g\\pm%.2g$  &  -\\\\\n",Ci1[3][0],Ci1[3][1]);
      printf("$\\bar\\ell_4$ & $%g\\pm%.2g$ &  $%g\\pm%.2g$\\\\\n",Ci1[4][0],Ci1[4][1],Ci[1][0],Ci[1][1]);
      printf("$P_4$ & 0  &  0\\\\\n");

      printf("\n\n");
      printf("$m_\\ell $ & $%g\\pm%.2g$  &  -  \\\\\n",C1[0][0],C1[0][1]);
      printf("$f_{\\pi}$ & $%g\\pm%.2g$  &  $%g\\pm%.2g$ \\\\\n",C1[1][0],C1[1][1],C2[0][0],C2[0][1]);
      free(tif);
      for (i=0;i<2;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    
    //////////////////////////
 /*   printf("\n\nPolynomial\n");
   
   
    tif=(double*) malloc(sizeof(double)*8);
    Ci=(double**) malloc(sizeof(double*)*8);

    fit=fit_Mpi_fw_polynomial(jack_files, head , jack_tot, mass_index, gjack );
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    Ci[2]=mean_and_error(argv[1],jack_tot,fit[2]);
    Ci[3]=mean_and_error(argv[1],jack_tot,fit[3]);
    Ci[4]=mean_and_error(argv[1],jack_tot,fit[4]);
    Ci[5]=mean_and_error(argv[1],jack_tot,fit[5]);
    Ci[6]=mean_and_error(argv[1],jack_tot,fit[6]);
    Ci[7]=mean_and_error(argv[1],jack_tot,fit[7]);
    
    printf("Bw=(%g+-%g)  fw=(%g+-%g)  P1=(%g+-%g)  P2=(%g+-%g)  P3=(%g+-%g)  P4=(%g+-%g)  P5=(%g+-%g) P6=(%g+-%g) \n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1],Ci[2][0],Ci[2][1],Ci[3][0],Ci[3][1],Ci[4][0],Ci[4][1],Ci[5][0],Ci[5][1],Ci[6][0],Ci[6][1],Ci[7][0],Ci[7][1]);

    
    for (j=0;j<jack_tot;j++){
        for (i=0;i<8;i++)
            tif[i]=fit[i][j];
        xi[j]=rtbis( m_over_f_pol,8,tif,  0.0001, 0.2, 1e-10);
        phys_point[0]=xi[j];
        phys_point[1]=1e+16;
        fw_phys[j]=Mw2_fw_polynomial(1, 2, phys_point,6,tif)/130.41;
    
    }
   Ci[0]=mean_and_error(argv[1],jack_tot,xi);
   printf("m_ud*w0=(%g+-%g)\n",Ci[0][0],Ci[0][1]);
   Ci[0]=mean_and_error(argv[1],jack_tot,fw_phys);
   printf("w0=(%g+-%g) Mev^-1   \t  (%g+-%g) fm \n",Ci[0][0],Ci[0][1],Ci[0][0]*197.326963,Ci[0][1]*197.326963);


    
    printf("\n\nPolynomial 6\n");
    free(tif);
    for (i=0;i<8;i++)
        free(Ci[i]);
    free(Ci);
    
    tif=(double*) malloc(sizeof(double)*6);
    Ci=(double**) malloc(sizeof(double*)*6);

    fit=fit_Mpi_fw_polynomial_6(jack_files, head , jack_tot, mass_index, gjack );
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    Ci[2]=mean_and_error(argv[1],jack_tot,fit[2]);
    Ci[3]=mean_and_error(argv[1],jack_tot,fit[3]);
    Ci[4]=mean_and_error(argv[1],jack_tot,fit[4]);
    Ci[5]=mean_and_error(argv[1],jack_tot,fit[5]);
    
    printf("Bw=(%g+-%g)  fw=(%g+-%g)  P1=(%g+-%g)  P2=(%g+-%g)   P4=(%g+-%g)  P5=(%g+-%g) \n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1],Ci[2][0],Ci[2][1],Ci[3][0],Ci[3][1],Ci[4][0],Ci[4][1],Ci[5][0],Ci[5][1]);

    double tif_pol[8];
    
    for (j=0;j<jack_tot;j++){
        for (i=0;i<6;i++)
            tif[i]=fit[i][j];
        xi[j]=rtbis( m_over_f_pol_6,6,tif, 0.001, 0.2, 1e-10);
        phys_point[0]=xi[j];
        phys_point[1]=1e+16;
        fw_phys[j]=Mw2_fw_polynomial_6(1, 2, phys_point,6,tif)/130.41;
         
        double test,dtest;
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,xi);
    printf("m_ud*w0=(%g+-%g)\n",Ci[0][0],Ci[0][1]);
    Ci[0]=mean_and_error(argv[1],jack_tot,fw_phys);
    printf("w0=(%g+-%g) Mev^-1   \t  (%g+-%g) fm \n",Ci[0][0],Ci[0][1],Ci[0][0]*197.326963,Ci[0][1]*197.326963);

    
        free(tif);
    for (i=0;i<6;i++)
        free(Ci[i]);
    free(Ci);*/
    /////////////////////////////////////////// K meson
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2

       printf("\n\nfK and MK^2\n");

    tif=(double*) malloc(sizeof(double)*3);
    Ci=(double**) malloc(sizeof(double*)*2);
    result.fkw=(double*) malloc(sizeof(double)*jack_tot);
    result.msw=(double*) malloc(sizeof(double)*jack_tot);
    
    result.fk_fpi=(double*) malloc(sizeof(double)*jack_tot);
    result.ms_mud=(double*) malloc(sizeof(double)*jack_tot);
    
    result.PMK=(double***) malloc(sizeof(double**)*3);
    for (i=0;i<3;i++){
        result.PMK[i]=(double**) malloc(sizeof(double*)*jack_tot);
        for (j=0;j<jack_tot;j++){
               result.PMK[i][j]=(double*) malloc(sizeof(double)*5);
        }
    }

    fit=fit_MK_double_chiral_FVE(jack_files, head , jack_tot, mass_index, gjack,  &result );
    //printf("\n\ntest:\n P[0]=%g P[0]=%g P[0]=%g P[0]=%g P[0]=%g \n", result.PMK[0][jack_tot-1][0],result.PMK[0][jack_tot-1][1],result.PMK[0][jack_tot-1][2],result.PMK[0][jack_tot-1][3],result.PMK[0][jack_tot-1][4]);
    for (j=0;j<jack_tot;j++){
        result.msw[j]=fit[0][j];
        result.fkw[j]=fit[1][j];
        
        result.ms_mud[j]=result.msw[j]/result.mlw[j];
       // printf("%f\n",result.ms_mud[j]);
        result.fk_fpi[j]=result.fkw[j]/result.fpiw[j];
        fit[0][j]=fit[0][j]/result.w0MeV[j];
        fit[1][j]=fit[1][j]/result.w0MeV[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    
    printf("Imposing $M_\\K =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %2g$ fm\n",v_MKMeV,err_MKMeV,v_w0fm,err_w0fm);
    printf("\\begin{gather}\n  m_{s}=(%g\\pm%.2g) MeV \\nn  \\\\ \n",Ci[0][0],Ci[0][1]);
    printf("   f_{K}=(%g\\pm%.2g) MeV \\nn \n\\end{gather} \n",Ci[1][0],Ci[1][1]);
    
     for (i=0;i<2;i++)
        free(Ci[i]);
     
    Ci[0]=mean_and_error(argv[1],jack_tot,result.ms_mud);
    Ci[1]=mean_and_error(argv[1],jack_tot,result.fk_fpi);
    
    printf("\\begin{gather}\n m_{s}/m_{ub}=(%g\\pm%.2g),\\quad \t  f_{K}/f_{\\pi}=(%g\\pm%.2g)  \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1]);

    
    
    
    free(tif);
    for (i=0;i<2;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2

       printf("\n\nfK and MK^2 fitting B0\n");

    tif=(double*) malloc(sizeof(double)*3);
    Ci=(double**) malloc(sizeof(double*)*2);
    
    fit=fit_MK_double_chiral_B0_FVE(jack_files, head , jack_tot, mass_index, gjack,  &result );
    for (j=0;j<jack_tot;j++){
        result.msw[j]=fit[0][j];
        result.fkw[j]=fit[1][j];
        
        result.ms_mud[j]=result.msw[j]/result.mlw[j];
       // printf("%f\n",result.ms_mud[j]);
        result.fk_fpi[j]=result.fkw[j]/result.fpiw[j];
        fit[0][j]=fit[0][j]/result.w0MeV[j];
        fit[1][j]=fit[1][j]/result.w0MeV[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    
    printf("Imposing $M_\\K =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %2g$ fm\n",v_MKMeV,err_MKMeV,v_w0fm,err_w0fm);
    printf("\\begin{gather}\n  m_{s}=(%g\\pm%.2g) MeV \\nn  \\\\ \n",Ci[0][0],Ci[0][1]);
    printf("   f_{K}=(%g\\pm%.2g) MeV \\nn \n\\end{gather} \n",Ci[1][0],Ci[1][1]);
    
     for (i=0;i<2;i++)
        free(Ci[i]);
     
    Ci[0]=mean_and_error(argv[1],jack_tot,result.ms_mud);
    Ci[1]=mean_and_error(argv[1],jack_tot,result.fk_fpi);
    
    printf("\\begin{gather}\n m_{s}/m_{ub}=(%g\\pm%.2g),\\quad \t  f_{K}/f_{\\pi}=(%g\\pm%.2g)  \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1]);

    
    
    
    free(tif);
    for (i=0;i<2;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2

       printf("\n\nfK from MK\n");

    tif=(double*) malloc(sizeof(double)*1);
    Ci=(double**) malloc(sizeof(double*)*2);
    
    fit=fit_fK_double_chiral_FVE_from_M(jack_files, head , jack_tot, mass_index, gjack,  &result );
    result.fKw_from_M=(double*) malloc(sizeof(double)*jack_tot);
    result.fk_fpi_from_M=(double*) malloc(sizeof(double)*jack_tot);
    
    for (j=0;j<jack_tot;j++){
        result.fKw_from_M[j]=fit[0][j];
        tmp1[j]=fit[0][j]/result.w0MeV[j];
        result.fk_fpi_from_M[j]=result.fKw_from_M[j]/result.fpiw_from_M[j];
    }

    Ci[0]=mean_and_error(argv[1],jack_tot,tmp1);
    Ci[1]=mean_and_error(argv[1],jack_tot,result.fk_fpi_from_M);
    printf("Imposing $M_\\K =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %2g$ fm\n",v_MKMeV,err_MKMeV,v_w0fm,err_w0fm);
    printf("   f_{K}=(%g\\pm%.2g) MeV \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1]);
    printf("\\begin{gather}\n   f_{K}/f_{\\pi}=(%g\\pm%.2g) MeV \\nn \n\\end{gather} \n",Ci[1][0],Ci[1][1]);
     
    
    free(Ci[0]);free(Ci[1]); free(fit);
    
    for (j=0;j<jack_tot;j++){  //putting in result.w0MeV the w0_estimate from f_pi_exp
        tmp1[j]=result.w0MeV[j];
        result.w0MeV[j]=w0_estimate[j]/197.326963;
    }
    fit=fit_fK_double_chiral_FVE_from_M(jack_files, head , jack_tot, mass_index, gjack,  &result );
    for (j=0;j<jack_tot;j++){  //restoring result.w0MeV 
        result.w0MeV[j]=tmp1[j];
    }
    
    printf("Imposing $M_\\K =%.2f \\pm %.2g$ and $f_pi=%.4f\\pm %2g$ fm\n",v_MKMeV,err_MKMeV,v_fpiMeV_exp,err_fpiMeV_exp);
    for (j=0;j<jack_tot;j++){
        tmp1[j]=fit[0][j]*197.326963/w0_estimate[j];
        result.fk_fpi_from_M[j]=tmp1[j]/result.fpiMeV_exp[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,tmp1);
    Ci[1]=mean_and_error(argv[1],jack_tot,result.fk_fpi_from_M);
    printf("   f_{K}=(%g\\pm%.2g) MeV \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1]);
    printf("\\begin{gather}\n   f_{K}/f_{\\pi}=(%g\\pm%.2g) MeV \\nn \n\\end{gather} \n",Ci[1][0],Ci[1][1]);
     
  
    
    free(tif);
    for (i=0;i<1;i++)
    {    free(Ci[i]);   free(fit[i]);}free(Ci[1]);
    free(Ci);free(fit);    
    
    printf("/////////////////////////////////////////////////////////// global fit ///////////////////////////////////////////////////////////\n ");
    struct fit_result fit_out;
    struct fit_type fit_info;

    fit_info.Npar=6;
    fit_info.Nvar=10;
    fit_info.N=3;
    fit_info.function=global_fK_from_M;
    fit_out=fit_fK_of_M(jack_files,  head ,jack_tot , gjack,mass_index, &result , fit_info);

    for(i=0;i<fit_info.Npar;i++){
        tmp=mean_and_error(argv[1],jack_tot,fit_out.P[i]);
        printf("P%d=%g   %g\n",i,tmp[0],tmp[1]);
    }
    printf("/////////////////////////////////////////////////////////// global fit OMEGA  Mpi MK  ///////////////////////////////////////////////////////////\n ");
/* 
    fit_info.Npar=4;
    fit_info.Nvar=3;
    fit_info.N=3;
    fit_info.function=global_Omega_Mpi_MK;
    fit_out=global_fit_Omega_from_Mpi_MK(jack_files,  head ,jack_tot , gjack,mass_index, &result , fit_info);

    for(i=0;i<fit_info.Npar;i++){
        tmp=mean_and_error(argv[1],jack_tot,fit_out.P[i]);
        printf("P%d=%g   %g\n",i,tmp[0],tmp[1]);
    }
    printf("/////////////////////////////////////////////////////////// global fit OMEGA  MK ///////////////////////////////////////////////////////////\n ");
 
    fit_info.Npar=4;
    fit_info.Nvar=3;
    fit_info.N=3;
    fit_info.function=global_Omega_MK;
    fit_out=global_fit_Omega_from_Mpi_MK(jack_files,  head ,jack_tot , gjack,mass_index, &result , fit_info);

    for(i=0;i<fit_info.Npar;i++){
        tmp=mean_and_error(argv[1],jack_tot,fit_out.P[i]);
        printf("P%d=%g   %g\n",i,tmp[0],tmp[1]);
    }
    printf("/////////////////////////////////////////////////////////// global fit OMEGA propMK ///////////////////////////////////////////////////////////\n ");
 
    fit_info.Npar=2;
    fit_info.Nvar=3;
    fit_info.N=3;
    fit_info.function=global_Omega_propMK;
    fit_out=global_fit_Omega_from_Mpi_MK(jack_files,  head ,jack_tot , gjack,mass_index, &result , fit_info);

    for(i=0;i<fit_info.Npar;i++){
        tmp=mean_and_error(argv[1],jack_tot,fit_out.P[i]);
        printf("P%d=%g   %g\n",i,tmp[0],tmp[1]);
    }
    
    */
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2

    printf("\n\nOMEGA MK\n");

    tif=(double*) malloc(sizeof(double)*1);
    Ci=(double**) malloc(sizeof(double*)*2);
    
    fit=fit_Omegaw0_from_M(jack_files, head , jack_tot, mass_index, gjack,  &result );
    
    for (j=0;j<jack_tot;j++){
        tmp1[j]=fit[0][j]*197.326963;
    }

    Ci[0]=mean_and_error(argv[1],jack_tot,tmp1);
    printf("Imposing $M_\\K =%.2f \\pm %.2g$ , $M_\\pi=%.4f\\pm %2g$ fm  and M_\\Omega=%g \\pm %g\n",v_MKMeV,err_MKMeV,v_MpiMeV,err_MpiMeV,v_MOmegaMeV,err_MOmegaMeV);
    printf("   w_0=(%g\\pm%.2g) fm \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1]);
     
    
    free(tif);
    for (i=0;i<1;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);    
    
    
/////////////////////////////////////////////////////////////////////////////////////////

    printf("\n\nfK and MK^2   P4=0\n");


    tif=(double*) malloc(sizeof(double)*3);
    Ci=(double**) malloc(sizeof(double*)*2);
    //result.fkw=(double*) malloc(sizeof(double)*jack_tot);
    //result.msw=(double*) malloc(sizeof(double)*jack_tot);
    
    result.fk_fpi=(double*) malloc(sizeof(double)*jack_tot);
    result.ms_mud=(double*) malloc(sizeof(double)*jack_tot);
    
    fit=fit_MK_double_chiral_FVE_P40(jack_files, head , jack_tot, mass_index, gjack,  &result );
    for (j=0;j<jack_tot;j++){
        result.msw[j]=fit[0][j];
        result.fkw[j]=fit[1][j];
        
        result.ms_mud[j]=result.msw[j]/result.mlw[j];
       // printf("%f\n",result.ms_mud[j]);
        result.fk_fpi[j]=result.fkw[j]/result.fpiw[j];
        fit[0][j]=fit[0][j]/result.w0MeV[j];
        fit[1][j]=fit[1][j]/result.w0MeV[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);

    
    
    printf("Imposing $M_\\K =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %2g$ fm\n",v_MKMeV,err_MKMeV,v_w0fm,err_w0fm);
    printf("\\begin{gather}\n   m_{s}=(%g\\pm%.2g) MeV \\nn  \\\\ \n",Ci[0][0],Ci[0][1]);
    printf("   f_{K}=(%g\\pm%.2g) MeV \\nn \n\\end{gather} \n",Ci[1][0],Ci[1][1]);
    
     for (i=0;i<2;i++)
        free(Ci[i]);
     
    Ci[0]=mean_and_error(argv[1],jack_tot,result.ms_mud);
    Ci[1]=mean_and_error(argv[1],jack_tot,result.fk_fpi);
    
    printf("\\begin{gather}\n m_{s}/m_{ub}=(%g\\pm%.2g),\\quad \t  f_{K}/f_{\\pi}=(%g\\pm%.2g)  \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1],Ci[1][0],Ci[1][1]);
    
    
    
    free(tif);
    for (i=0;i<2;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    
  
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2
       printf("\n\nfD and MD\n");

    tif=(double*) malloc(sizeof(double)*3);
    Ci=(double**) malloc(sizeof(double*)*2);
    result.fDw=(double*) malloc(sizeof(double)*jack_tot);
    result.mcw=(double*) malloc(sizeof(double)*jack_tot);
    
   
    result.mc_ms=(double*) malloc(sizeof(double)*jack_tot);
    
    fit=fit_MD_double_chiral(jack_files, head , jack_tot, mass_index, gjack,  result );
    
    for (j=0;j<jack_tot;j++){
        result.mcw[j]=fit[0][j];
        result.fDw[j]=fit[1][j];
        
        result.mc_ms[j]=result.mcw[j]/result.msw[j];
       // printf("%f\n",result.ms_mud[j]);
        fit[0][j]=fit[0][j]/result.w0MeV[j];
        fit[1][j]=fit[1][j]/result.w0MeV[j];
    }

    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    
    printf("Imposing $M_\\D =%.2f \\pm %.2g$ MeV and $w_0=%.4f\\pm %2g$ fm\n",v_MDMeV,err_MDMeV,v_w0fm,err_w0fm);
    printf("\\begin{gather}\n  m_{c}=(%g\\pm%.2g) MeV \\nn  \\\\ \n",Ci[0][0],Ci[0][1]);
    printf("   f_{D}=(%g\\pm%.2g) MeV \\nn \n\\end{gather} \n",Ci[1][0],Ci[1][1]);
    
     for (i=0;i<2;i++)
        free(Ci[i]);
     
    Ci[0]=mean_and_error(argv[1],jack_tot,result.ms_mud);
    
    printf("\\begin{gather}\n m_{c}/m_{s}=(%g\\pm%.2g) \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1]);
    
    
    
    free(tif);
    for (i=0;i<1;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2
       printf("\n\nfD and MD P30\n");

    tif=(double*) malloc(sizeof(double)*3);
    Ci=(double**) malloc(sizeof(double*)*2);
    result.fDw=(double*) malloc(sizeof(double)*jack_tot);
    result.mcw=(double*) malloc(sizeof(double)*jack_tot);
    
   
    result.mc_ms=(double*) malloc(sizeof(double)*jack_tot);
    
    fit=fit_MD_double_chiral_P30(jack_files, head , jack_tot, mass_index, gjack,  result );
    
    for (j=0;j<jack_tot;j++){
        result.mcw[j]=fit[0][j];
        result.fDw[j]=fit[1][j];
        
        result.mc_ms[j]=result.mcw[j]/result.msw[j];
       // printf("%f\n",result.ms_mud[j]);
        fit[0][j]=fit[0][j]/result.w0MeV[j];
        fit[1][j]=fit[1][j]/result.w0MeV[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    
    printf("Imposing $M_D =%.2f \\pm %.2g$ MeV and $w_0=%.4f\\pm %2g$ fm\n",v_MDMeV,err_MDMeV,v_w0fm,err_w0fm);
    printf("\\begin{gather}\n  m_{c}=(%g\\pm%.2g) MeV \\nn  \\\\ \n",Ci[0][0],Ci[0][1]);
    printf("   f_{D}=(%g\\pm%.2g) MeV \\nn \n\\end{gather} \n",Ci[1][0],Ci[1][1]);

    
     for (i=0;i<2;i++)
     {free(Ci[i]);free(fit[i]);}
     
    Ci[0]=mean_and_error(argv[1],jack_tot,result.mc_ms);
    
    printf("\\begin{gather}\n m_{c}/m_{s}=(%g\\pm%.2g) \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1]);
    
    
    
    free(tif);

    for (i=0;i<1;i++)
    {    free(Ci[i]);   }
    free(Ci);free(fit);
//////////////////////////////////////////////////////////////////////////////////////////////////////////    
    printf("\n\nfD from MD\n");

    tif=(double*) malloc(sizeof(double)*1);
    Ci=(double**) malloc(sizeof(double*)*1);
    result.fDw_from_M=(double*) malloc(sizeof(double)*jack_tot);

    for (j=0;j<jack_tot;j++){  //putting in result.w0MeV the w0_estimate from f_pi_exp
        tmp1[j]=result.w0MeV[j];
        result.w0MeV[j]=w0_estimate[j]/197.326963;
    }
    fit=fit_fD_chiral_continuum_from_M(jack_files, head , jack_tot, mass_index, gjack,  &result );
    for (j=0;j<jack_tot;j++){  //restoring result.w0MeV 
        result.w0MeV[j]=tmp1[j];
    }
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    
    printf("Imposing $M_\\D =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %2g$ fm\n",v_MDMeV,err_MDMeV,v_w0fm,err_w0fm);
    
    printf("   f_{D}=(%g\\pm%.2g) MeV \\nn \n\\end{gather} \n",Ci[0][0]/v_w0MeV,Ci[0][1]/v_w0MeV);

    for (j=0;j<jack_tot;j++)
                result.fDw_from_M[j]=fit[0][j];

     
     
    free(tif);
    for (i=0;i<1;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);    
    
  
    
    
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////c2
       printf("\n\nfDs and MDs^2\n");

    tif=(double*) malloc(sizeof(double)*3);
    Ci=(double**) malloc(sizeof(double*)*2);
    result.fDsw=(double*) malloc(sizeof(double)*jack_tot);
    result.mcw=(double*) malloc(sizeof(double)*jack_tot);
    
   
    result.mc_ms=(double*) malloc(sizeof(double)*jack_tot);
    
    fit=fit_MDs_double_chiral(jack_files, head , jack_tot, mass_index, gjack,  result );
    
    for (j=0;j<jack_tot;j++){
        result.mcw[j]=fit[0][j];
        result.fDsw[j]=fit[1][j];
        
        result.mc_ms[j]=result.mcw[j]/result.msw[j];
       // printf("%f\n",result.ms_mud[j]);
        fit[0][j]=fit[0][j]/result.w0MeV[j];
        fit[1][j]=fit[1][j]/result.w0MeV[j];
    }

    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    Ci[1]=mean_and_error(argv[1],jack_tot,fit[1]);
    
    printf("Imposing $M_\\Ds =%.2f \\pm %.2g$ MeV and $w_0=%.4f\\pm %2g$ fm\n",v_MDsMeV,err_MDsMeV,v_w0fm,err_w0fm);
    printf("\\begin{gather}\n  m_{c}=(%g\\pm%.2g) MeV \\nn  \\\\ \n",Ci[0][0],Ci[0][1]);
    printf("   f_{Ds}=(%g\\pm%.2g) MeV \\nn \n\\end{gather} \n",Ci[1][0],Ci[1][1]);
    
     for (i=0;i<2;i++)
        free(Ci[i]);
     
    Ci[0]=mean_and_error(argv[1],jack_tot,result.ms_mud);
    
    printf("\\begin{gather}\n m_{c}/m_{s}=(%g\\pm%.2g) \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1]);
    
    
    
    free(tif);
    for (i=0;i<1;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////    
    printf("\n\nfDs from Mpi, MK, MD\n");

    tif=(double*) malloc(sizeof(double)*1);
    Ci=(double**) malloc(sizeof(double*)*1);
    result.fDsw_from_M=(double*) malloc(sizeof(double)*jack_tot);
    result.fDs_fD_from_M=(double*) malloc(sizeof(double)*jack_tot);

    for (j=0;j<jack_tot;j++){  //putting in result.w0MeV the w0_estimate from f_pi_exp
        tmp1[j]=result.w0MeV[j];
        result.w0MeV[j]=w0_estimate[j]/197.326963;
    }
    fit=fit_fDs_chiral_continuum_from_M(jack_files, head , jack_tot, mass_index, gjack,  &result );
    for (j=0;j<jack_tot;j++){  //restoring result.w0MeV 
        result.w0MeV[j]=tmp1[j];
    }
    
    Ci[0]=mean_and_error(argv[1],jack_tot,fit[0]);
    
    printf("Imposing $M_\\D =%.2f \\pm %.2g$ and $w_0=%.4f\\pm %2g$ fm\n",v_MDMeV,err_MDMeV,v_w0fm,err_w0fm);
    
    printf("   f_{Ds}=(%g\\pm%.2g) MeV \\nn \n\\end{gather} \n",Ci[0][0]/v_w0MeV,Ci[0][1]/v_w0MeV);
    for (j=0;j<jack_tot;j++){
        result.fDsw_from_M[j]=fit[0][j];
        result.fDs_fD_from_M[j]=result.fDsw_from_M[j]/result.fDw_from_M[j];
        
    }
    free(Ci[0]);
    Ci[0]=mean_and_error(argv[1],jack_tot,result.fDs_fD_from_M);
    printf("   f_{Ds}/f_{D}=(%g\\pm%.2g)  \\nn \n\\end{gather} \n",Ci[0][0],Ci[0][1]);

     
    free(tif);
    for (i=0;i<1;i++)
    {    free(Ci[i]);   free(fit[i]);}
    free(Ci);free(fit);    
  
       
  
    
    
    
    

    return 0;
}
