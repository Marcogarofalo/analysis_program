#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>


#include "resampling.hpp"
#include "read.hpp"
#include "linear_fit.hpp"
#include "global.hpp"
#include "mutils.hpp"

 


static void  print_file_head(FILE *stream)
{
    int i,dsize;
    double *dstd;
    
    fprintf(stream,"%d\n",file_head.twist);
    fprintf(stream,"%d\n",file_head.nf);
    fprintf(stream,"%d\n",file_head.nsrc);
    fprintf(stream,"%d\n",file_head.l0);
    fprintf(stream,"%d\n",file_head.l1);
    fprintf(stream,"%d\n",file_head.l2);
    fprintf(stream,"%d\n",file_head.l3);
    fprintf(stream,"%d\n",file_head.nk);
    fprintf(stream,"%d\n",file_head.nmoms);
    
    fprintf(stream,"%f\n",file_head.beta);
    fprintf(stream,"%f\n",file_head.ksea);
    fprintf(stream,"%f\n",file_head.musea);
    fprintf(stream,"%f\n",file_head.csw);
   
    for(i=0;i<4*file_head.nk;++i)
            fprintf(stream,"%f\n",file_head.k[i]);

    
    for(i=0;i<file_head.nmoms;++i)
           fprintf(stream,"%f  %f   %f   %f\n",file_head.mom[i][0],file_head.mom[i][1],file_head.mom[i][2],file_head.mom[i][3]);
}

static void  read_file_head_bin(FILE *stream)
{
    int i,dsize;
    double *dstd;
    
    fread(&file_head.twist,sizeof(int),1,stream);
    fread(&file_head.nf,sizeof(int),1,stream);
    fread(&file_head.nsrc,sizeof(int),1,stream);
    fread(&file_head.l0,sizeof(int),1,stream);
    fread(&file_head.l1,sizeof(int),1,stream);
    fread(&file_head.l2,sizeof(int),1,stream);
    fread(&file_head.l3,sizeof(int),1,stream);
    fread(&file_head.nk,sizeof(int),1,stream);
    fread(&file_head.nmoms,sizeof(int),1,stream);
    
    fread(&file_head.beta,sizeof(double),1,stream);
    fread(&file_head.ksea,sizeof(double),1,stream);
    fread(&file_head.musea,sizeof(double),1,stream);
    fread(&file_head.csw,sizeof(double),1,stream);
   
    file_head.k=(double*) malloc(sizeof(double)*2*file_head.nk);
    for(i=0;i<2*file_head.nk;++i)
    	fread(&file_head.k[i],sizeof(double),1,stream);
    
    file_head.mom=(double**) malloc(sizeof(double*)*file_head.nmoms);
    for(i=0;i<file_head.nmoms;i++) {
    	file_head.mom[i]=(double*) malloc(sizeof(double)*4);
        fread(&file_head.mom[i][0],sizeof(double),1,stream);
        fread(&file_head.mom[i][1],sizeof(double),1,stream);
        fread(&file_head.mom[i][2],sizeof(double),1,stream);
        fread(&file_head.mom[i][3],sizeof(double),1,stream);

    }
}

static void  write_file_head(FILE *stream)
{
    int i,dsize;
    double *dstd;
    
    fwrite(&file_head.twist,sizeof(int),1,stream);
    fwrite(&file_head.nf,sizeof(int),1,stream);
    fwrite(&file_head.nsrc,sizeof(int),1,stream);
    fwrite(&file_head.l0,sizeof(int),1,stream);
    fwrite(&file_head.l1,sizeof(int),1,stream);
    fwrite(&file_head.l2,sizeof(int),1,stream);
    fwrite(&file_head.l3,sizeof(int),1,stream);
    fwrite(&file_head.nk,sizeof(int),1,stream);
    fwrite(&file_head.nmoms,sizeof(int),1,stream);
    
    fwrite(&file_head.beta,sizeof(double),1,stream);
    fwrite(&file_head.ksea,sizeof(double),1,stream);
    fwrite(&file_head.musea,sizeof(double),1,stream);
    fwrite(&file_head.csw,sizeof(double),1,stream);
   
    fwrite(file_head.k,sizeof(double),2*file_head.nk,stream);

    for(i=0;i<file_head.nmoms;i++)  
        fwrite(file_head.mom[i],sizeof(double),4,stream);
}

void read_nconfs(int *s, int *c, FILE *stream){

   FILE *f1;
   long int tmp;

   fread(s,sizeof(int),1,stream);
   f1=stream;
   
   fseek(stream, 0, SEEK_END);
   tmp = ftell(stream);
   tmp-= sizeof(double)* (file_head.nmoms*4 + file_head.nk*2+4 )+ sizeof(int)*9 ;
   tmp-= 2*sizeof(int);
   (*c)= tmp/ (sizeof(int)+ (*s)*sizeof(double) );
 
  rewind(stream);
  read_file_head_bin(stream);
  fread(s,sizeof(int),1,stream);

printf("size=%d  confs=%d\n",*s,*c);
}


double *constant_fit(int M, double in){
    double *r;
    
    r=(double*) malloc(sizeof(double)*M);
    r[0]=1.;
    
    return r;
}
double M_eff( int t, double ***in){
    double mass;
 
    mass=acosh( (in[0][t+1][0]+ in[0][t-1][0])/(2*in[0][t][0]) );
    return mass;
}

double ratio( int t, double ***in){
    double mass;
    double a,b,c,d;    
     
    a=in[0][t][0];
    b=in[1][t][0];
    c=in[2][t][0];
    d=in[3][t][0];
    
    mass=( a*b/(c*d) ) ;
   // mass=in[0][t][0];
    return mass;
}

static int index_twopt(int si,int ii,int ix0,int imom1,int imom2,int ik1,int ik2)
{
    int nk,nmoms;

    nk=file_head.nk;
    nmoms=file_head.nmoms;

    return ii+si*(ix0+file_head.l0*(imom1+nmoms*(imom2+nmoms*(ik1+nk*ik2))));
}
static int index_threept(int si,int ii,int ix0,int imom1,int imom2,int ik1,int ik2,int ik3)
{
    int nk,nmoms;

    nk=file_head.nk;
    nmoms=file_head.nmoms;

    return ii+si*(ix0+file_head.l0*(imom1+nmoms*(imom2+nmoms*(ik1+nk*(ik2+nk*ik3)))));
} 
int main(){
   int size;
   int i,j,t;
   FILE  *f=NULL;
   int *iconf,confs;
   double ****data, **out,**tmp; 
   char c;
   clock_t t1,t2;
   double *in;
 
   double *fit,***y,*x,*m,*me;
   FILE *outfile=NULL;

   double ****conf_jack,**r,**mt,**met;
   int Ncorr=1;
   outfile=fopen("out_E0_temp","w+");
   error(outfile==NULL,1,"main ",
         "Unable to open output file");


   double E_B,E_Pi, x_SCHET,q2,vec_pB,vec_pPi;
   int Neff,Njack;
   t1=clock();

   f=fopen("./meas_2pts.dat","r");
   f=fopen("./to_read_bin.dat","r");
   if (f==NULL) {printf("to_read file not found\n"); exit(0);}

     
   read_file_head_bin(f);
   print_file_head(outfile);
   fflush(outfile);

   read_nconfs(&size,&confs,f);   
   int aaa;
   fread(&aaa,sizeof(int),1,f);
   printf("size=%d  confs=%d   aaa=%d\n",size,confs,aaa);
   
   iconf=(int*) malloc(sizeof(int)*confs);
   out=(double**) malloc(sizeof(double*)*confs);
   
   
   printf("size=%d  confs=%d\n",size,confs);
   printf("reading confs:\n");
   for(i=0;i<confs;i++){
       fread(&iconf[i],sizeof(int),1,f);
       out[i]=(double*) malloc(sizeof(double)*size);
       fread(out[i],sizeof(double),size,f);
       printf("%d\t",iconf[i]);
   }
   printf("\n");
   
   fclose(f);
////
   int bin=10;
   FILE *fbin=fopen("meas_2pts_bin10.dat","w+");
   int conf_out=aaa/bin;
   write_file_head(fbin);
   fwrite(&size,sizeof(int),1,fbin);
   fwrite(&conf_out,sizeof(int),1,fbin);
 
   
   int rconf;
   rconf=(confs/bin)*bin;
   for (i=0;i<rconf;i+=bin){
     for (t=1;t<bin;t++){
        for (j=0;j<size;j++){
     	   out[i][j]+=out[i+t][j];
        }
     }
     for (j=0;j<size;j++){
        out[i][j]/=(double) bin;
     }
       fwrite(&iconf[i],sizeof(int),1,fbin);
       fwrite(out[i],sizeof(double),size,fbin);
     
   }

   fclose(fbin);
/////
   //copyng the input
   int var=1,si;
    

   data=(double****) malloc(sizeof(double***)*confs);
   for (i=0;i<confs;i++){
        data[i]=(double***) malloc(sizeof(double**)*var);
        for (j=0;j<var;j++){
            data[i][j]=(double**) malloc(sizeof(double*)*file_head.l0);
            for(t=0;t<file_head.l0;t++)
                data[i][j][t]=(double*) malloc(2*sizeof(double));
        }
   }
   
   int ii, imom1, imom2, ik1,ik2,ik3;
   ii=0;
   imom1=0;
   imom2=0;
   ik1=0;
   ik2=0;
   ik3=0;
   int tmin=12, tmax=15;
  
   int index;
   double **E0;
   E0=(double**) malloc(sizeof(double*)*size);
   for (i=0;i<size ;i++)
	E0[i]=(double*) malloc(sizeof(double)*2); 

   si=size/(4*(file_head.nk*(file_head.nk+1)/2)*file_head.nmoms*file_head.nmoms*file_head.l0*2);

   for(ik2=0;ik2<file_head.nk;ik2++){
   for(ik1=0;ik1<file_head.nk;ik1++){
   for(imom2=0;imom2<file_head.nmoms;imom2++){
   for(imom1=0;imom1<file_head.nmoms;imom1++){


   for (i=0;i<confs;i++){
   index=2*index_twopt(si,ii,0,imom1,imom2,ik1,ik2);
//      index =(ii+si*(0 +file_head.l0 *(imom1+file_head.nmoms*(imom2+file_head.nmoms*(ik1+file_head.nk*(ik2))))))*2;   
       for(t=0;t<file_head.l0;t++){
           data[i][0][t][0]=out[i][index];
           data[i][0][t][1]=out[i][index+1];
	   index+=si*2;	
       }
   }
   Neff=confs;
   symmetrise_corr(Neff, 0, file_head.l0,data);
   conf_jack=create_jack(Neff, 1, file_head.l0, data);
  
   Njack=Neff+1;

////////////////////allocation
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
        r[i]=(double*) malloc(sizeof(double)*Njack);
    
   mt=(double**) malloc(sizeof(double*)*file_head.l0);
   fit=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);
   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(tmax-tmin));
        for (i=tmin;i<tmax;i++){
	       y[j][i-tmin]=(double*) malloc(sizeof(double)*2);
        }
    }
    x=(double*) malloc(sizeof(double)*(tmax-tmin));
////////////////end 

/////////////////M_PS
   fprintf(outfile,"#m_eff(t) from  propagators:1) mu %g theta %g 2) mu %g  theta %g\n",file_head.k[ik1+file_head.nk],file_head.mom[imom1][1],file_head.k[ik2+file_head.nk], file_head.mom[imom2][1] );
   
   for(i=1;i<file_head.l0/2;i++){    
        for (j=0;j<Njack;j++){
            r[i][j]=M_eff(i,conf_jack[j]);
            if (i>=tmin && i<tmax){
                y[j][i-tmin][0]=r[i][j];
            }

        }
        mt[i]=mean_and_error_jack(Neff, r[i]);
        
        if (i>=tmin && i<tmax){
            for (j=0;j<Njack;j++){
                y[j][i-tmin][1]=mt[i][1];
            }
        }
    fprintf(outfile,"%d   %.15g    %.15g\n",i,mt[i][0],mt[i][1]);
    free(mt[i]);
    }

    for (j=0;j<Njack;j++){
        tmp=linear_fit( tmax-tmin, x, y[j],  1,constant_fit );
        fit[j]=tmp[0][0];
        free(tmp[0]);free(tmp);
    }
    m=mean_and_error_jack(Njack, fit);
    index=2*index_twopt(si,ii,0,imom1,imom2,ik1,ik2);
    E0[index][0]=m[0];
    E0[index][1]=m[1];	
    fprintf(outfile,"\n\n #M_PS fit in [%d,%d]\n  %.15g    %.15g\n\n\n",tmin,tmax-1,m[0],m[1]);

    fflush(outfile);
/////////////////////free memory
    for(i=0;i<file_head.l0;i++)
   	  free(r[i]);
    free(r); free(fit);
    free(m);free(mt);
    for (j=0;j<Njack;j++){
        for (i=tmin;i<tmax;i++){
	       free(y[j][i-tmin]);
        }
        free(y[j]);
    }
    free(y);
    free(x);
//////////////////////////end
   }}}} // end loop imom1 imom2 ik1 ik2
   fclose(outfile);
   for (i=0;i<confs;i++){
        for (j=0;j<var;j++){
            for(t=0;t<file_head.l0;t++)
               free(data[i][j][t]);
            free(data[i][j]);
        }
        free(data[i]);
   }
   for (i=0;i<confs;i++)
       free(out[i]);
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //3pt
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   tmin=12; tmax=16;
   
   printf("starting 3pt");
   outfile=fopen("out_3pt","w+");
   error(outfile==NULL,1,"main ",
         "Unable to open output file");

   f=fopen("./meas_3pt.dat","r");
   if (f==NULL) {printf("to_read file not found\n"); exit(0);}

   read_file_head_bin(f);
   print_file_head(outfile);
   fflush(outfile);

   read_nconfs(&size,&confs,f);   
   printf("reading confs:\n");
   for(i=0;i<confs;i++){
       fread(&iconf[i],sizeof(int),1,f);
       out[i]=(double*) malloc(sizeof(double)*size);
       fread(out[i],sizeof(double),size,f);
       printf("%d\t",iconf[i]);
       fflush(stdout);
   }
   printf("\n");
   fclose(f);



////
   fbin=fopen("oPVmuPo-sss_bin10.dat","w");
   write_file_head(fbin);
   fwrite(&size,sizeof(int),1,fbin);
 
//   confs=168;
   rconf=(confs/bin)*bin;
   for (i=0;i<rconf;i+=bin){
     for (t=1;t<bin;t++){
        for (j=0;j<size;j++){
     	   out[i][j]+=out[i+t][j];
        }
     }
     for (j=0;j<size;j++){
        out[i][j]/=(double) bin;
     }
       fwrite(&iconf[i],sizeof(int),1,fbin);
       fwrite(out[i],sizeof(double),size,fbin);
     
   }

   fclose(fbin);

/////



   si=size/(file_head.nk*file_head.nk*file_head.nk*file_head.nmoms*file_head.nmoms*file_head.l0*2);
   if (si==0)
       si=size/(file_head.nk*file_head.nk*file_head.nmoms*file_head.nmoms*file_head.l0*2);

   var=4;

   data=(double****) malloc(sizeof(double***)*confs);
   for (i=0;i<confs;i++){
        data[i]=(double***) malloc(sizeof(double**)*var);
        for (j=0;j<var;j++){
            data[i][j]=(double**) malloc(sizeof(double*)*file_head.l0);
            for(t=0;t<file_head.l0;t++)
                data[i][j][t]=(double*) malloc(2*sizeof(double));
        }
   }


   for(ik3=file_head.nk-2;ik3<file_head.nk;++ik3){
   for(ik2=0;ik2<file_head.nk;++ik2){
   for(ik1=file_head.nk-2;ik1<file_head.nk;++ik1){
   for(imom2=0;imom2<file_head.nmoms;imom2++){
   for(imom1=0;imom1<file_head.nmoms;imom1++){


   for (i=0;i<confs;i++){
       index=2* index_threept(si,ii,0,imom1,imom2,ik1,ik2,ik3);
       for(t=0;t<file_head.l0;t++){
           data[i][0][t][0]=out[i][index];
           data[i][0][t][1]=out[i][index+1];
	   index+=si*2;	
       }

       index=2* index_threept(si,ii,0,imom2,imom1,ik2,ik1,ik3);
       for(t=0;t<file_head.l0;t++){
           data[i][1][t][0]=out[i][index];
           data[i][1][t][1]=out[i][index+1];
	   index+=si*2;	
       }

       index=2* index_threept(si,ii,0,imom2,imom2,ik2,ik2,ik3);
       for(t=0;t<file_head.l0;t++){
           data[i][2][t][0]=out[i][index];
           data[i][2][t][1]=out[i][index+1];
	   index+=si*2;	
       }
       
       index=2* index_threept(si,ii,0,imom1,imom1,ik1,ik1,ik3);
       for(t=0;t<file_head.l0;t++){
           data[i][3][t][0]=out[i][index];
           data[i][3][t][1]=out[i][index+1];
	   index+=si*2;	
       }

   }

   if (imom1 >0  && imom2>0){

	   for (i=0;i<confs;i++){
	       index=2* index_threept(si,ii,0,imom2,imom1,ik1,ik2,ik3);
	       for(t=0;t<file_head.l0;t++){
		   data[i][0][t][0]+=out[i][index];
		   data[i][0][t][0]/=2;
		   data[i][0][t][1]+=out[i][index+1];
		   data[i][0][t][1]/=2;
		   index+=si*2;	
	       }

	       index=2* index_threept(si,ii,0,imom1,imom2,ik2,ik1,ik3);
	       for(t=0;t<file_head.l0;t++){
		   data[i][1][t][0]+=out[i][index];
		   data[i][0][t][0]/=2;
		   data[i][1][t][1]+=out[i][index+1];
		   data[i][0][t][1]/=2;
		   index+=si*2;	
	       }
            }

   }


   Neff=confs;
   symmetrise_corr(Neff, 0, file_head.l0,data);
   symmetrise_corr(Neff, 1, file_head.l0,data);
   symmetrise_corr(Neff, 2, file_head.l0,data);
   symmetrise_corr(Neff, 3, file_head.l0,data);
   conf_jack=create_jack(Neff, var, file_head.l0, data);
  
   Njack=Neff+1;
////////////////////allocation
   r=(double**) malloc(sizeof(double*)*file_head.l0);
   for(i=0;i<file_head.l0;i++)
        r[i]=(double*) malloc(sizeof(double)*Njack);
    
   met=(double**) malloc(sizeof(double*)*file_head.l0);
   fit=(double*) malloc(sizeof(double)*Njack);
   y=(double***) malloc(sizeof(double**)*Njack);
   for (j=0;j<Njack;j++){
        y[j]=(double**) malloc(sizeof(double*)*(tmax-tmin));
        for (i=tmin;i<tmax;i++){
	       y[j][i-tmin]=(double*) malloc(sizeof(double)*2);
        }
    }
    x=(double*) malloc(sizeof(double)*(tmax-tmin));
////////////////end 

///////////////////////kynematic
   
   fprintf(outfile,"#matrix elements(t) from  propagators:1) mu %g theta %g 2) mu %g theta %g 3) mu %g theta 0\n",file_head.k[ik1+file_head.nk],file_head.mom[imom1][1],file_head.k[ik2+file_head.nk], file_head.mom[imom2][1], file_head.k[ik3+file_head.nk] );   
   index=2*index_twopt(1,0,0,imom2,0,ik2,ik3);
   E_B=E0[index][0];
   fprintf(outfile,"E_B=%g ",E0[index][0]);
   index=2*index_twopt(1,0,0,imom1,0,ik3,ik3);
   E_Pi=E0[index][0];
   fprintf(outfile,"E_Pi=%g\n",E0[index][0]);
   vec_pB=sqrt(3)*2*3.14159265358979 *file_head.mom[imom2][1]/file_head.l1 ;
   vec_pPi=sqrt(3)*2*3.14159265358979 *file_head.mom[imom1][1]/file_head.l1 ;

   fprintf(outfile,"p_B=%g p_Pi=%g\n",vec_pB,vec_pPi);
   x_SCHET=2*(   E_B*E_Pi- vec_pB*vec_pPi )/ ( (E_B*E_B- vec_pB*vec_pB) );
   q2=(E_B-E_Pi)*(E_B-E_Pi)  -   (vec_pB-vec_pPi)*(vec_pB-vec_pPi);
   fprintf(outfile,"x_SCHET=%g q^2=%g\n",x_SCHET,q2);
   
   
/////////////////matrix element
   for(i=1;i<file_head.l0/2;i++){    
        for (j=0;j<Njack;j++){
            r[i][j]=ratio(i,conf_jack[j]);
////            r[i][j]=r[i][j]*4.*E_B*E_Pi;
           // r[i][j]=sqrt(r[i][j]*4.*E_B*E_Pi);
 
            if (i>=tmin && i<tmax){
                y[j][i-tmin][0]=r[i][j];
            }

        }
        met[i]=mean_and_error_jack(Neff, r[i]);
        
        if (i>=tmin && i<tmax){
            for (j=0;j<Njack;j++){
                y[j][i-tmin][1]=met[i][1];
            }
        }
    fprintf(outfile,"%d   %.15g    %.15g\n",i,met[i][0],met[i][1]);
    free(met[i]);
    }

    for (j=0;j<Njack;j++){
        tmp=linear_fit( tmax-tmin, x, y[j],  1,constant_fit );
        fit[j]=tmp[0][0];
        free(tmp[0]);free(tmp);
    }
    me=mean_and_error_jack(Njack, fit);
    fprintf(outfile,"\n\n #matrix elemnt fit in [%d,%d]\n  %.15g    %.15g\n\n\n",tmin,tmax-1,me[0],me[1]);

    fflush(outfile);
/////////////////////free memory
    for(i=0;i<file_head.l0;i++)
   	  free(r[i]);
    free(r); free(fit);
    free(me);free(met);
    for (j=0;j<Njack;j++){
        for (i=tmin;i<tmax;i++){
	       free(y[j][i-tmin]);
        }
        free(y[j]);
    }
    free(y);
    free(x);
//////////////////////////end


   }}}}}///end loop imom1 imom2 ik1 ik2 ik3
return 0;
}
