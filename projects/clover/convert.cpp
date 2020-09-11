#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include "read.hpp"
#include "global.hpp"
#include "mutils.hpp"



static void  read_file_head(FILE *stream)
{
    int i,dsize;
    double *dstd;
    
    fscanf(stream,"%d\n",&file_head.twist);
    fscanf(stream,"%d\n",&file_head.nf);
    fscanf(stream,"%d\n",&file_head.nsrc);
    fscanf(stream,"%d\n",&file_head.l0);
    fscanf(stream,"%d\n",&file_head.l1);
    fscanf(stream,"%d\n",&file_head.l2);
    fscanf(stream,"%d\n",&file_head.l3);
    fscanf(stream,"%d\n",&file_head.nk);
    fscanf(stream,"%d\n",&file_head.nmoms);
    
    fscanf(stream,"%lf\n",&file_head.beta);
    fscanf(stream,"%lf\n",&file_head.ksea);
    fscanf(stream,"%lf\n",&file_head.musea);
    fscanf(stream,"%lf\n",&file_head.csw);
   
    file_head.k=(double*) malloc(sizeof(double)*2*file_head.nk);
    for(i=0;i<2*file_head.nk;++i)
            fscanf(stream,"%lf\n",&file_head.k[i]);

//    file_head.mom=malloc(sizeof(file_head.mom[4])*file_head.nmoms);
    file_head.mom=(double**)  malloc(sizeof(double*)*file_head.nmoms);
    
    for(i=0;i<file_head.nmoms;++i){
       file_head.mom[i]=(double*) malloc(sizeof(double)*4);    
       fscanf(stream,"%lf  %lf   %lf   %lf\n",&file_head.mom[i][0],&file_head.mom[i][1],&file_head.mom[i][2],&file_head.mom[i][3]);
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


 
int main(int argc, char **argv){
   int size;
   int tmp,j;
   FILE  *f=NULL;
   int *iconf,confs;
   double ***data, **out; 
   char c;
   clock_t t1,t2;
   
   error(argc!=2,1,"main ",
         "usage: ./convert ASCI_file");

   t1=clock();

   f=fopen(argv[1],"r");
   error(f==NULL,1,"main ",
         "file %s not found",argv[1]);
   
  printf("HERE\n") ;
   read_file_head(f);
   
   fscanf(f,"%d\n",&size);
int   aaa ;
   fscanf(f,"%d\n",&aaa);
   
   long int i=0  ;
   // Extract characters from file and store in character c
   for (c = getc(f); c != EOF; c = getc(f))
       if (c == '\n') // Increment count if this character is newline   
        i = i + 1;
// i=1761177775-14-2*file_head.nk-file_head.nmoms;
   confs=i/( (size/2) +1 );
   printf("number of configuration are %d,  %ld  size=%d,  confs=%d  \n",confs, i , size,aaa);    
       
   fclose(f);
   f=fopen(argv[1],"r");
   error(f==NULL,1,"main ",
         "file %s not found",argv[1]);
   read_file_head(f);
   
   fscanf(f,"%d\n",&size);
   fscanf(f,"%d\n",&aaa);
   iconf=(int*) malloc(sizeof(int)*confs);
 
   FILE *f1;
   char sout[100];
   printf("HERE1\n");
   sprintf(sout,"%s_bin.dat",argv[1]);
   printf("HERE2\n");
   f1=fopen(sout,"w+");
   if (f1==NULL) {printf("unable to open binary file\n"); exit(0);}
   write_file_head(f1);
   printf("HERE1\n");
   
   fwrite(&size,sizeof(int),1,f1);
   fwrite(&aaa,sizeof(int),1,f1);
   fflush(f1);
long   int ii,si,t,index;

   data=(double***) malloc(sizeof(double**)*confs);       
   out=(double**) malloc(sizeof(double*)*confs);

 

   fflush(stdout);
       
   si=size/( 4*(file_head.nk*(file_head.nk+1)/2)*file_head.nmoms*file_head.nmoms*file_head.l0*2);
printf("si=%ld\n",si);
long  int count=0;
   for(i=0;i<confs;i++){
       data[i]=(double**) malloc(sizeof(double*)*2);       
       data[i][0]=(double*) malloc( sizeof(double)*size/2);
       data[i][1]=(double*) malloc(sizeof(double)*size/2); 
       out[i]=(double*) malloc(sizeof(double)*(size)); 
       fscanf(f,"%d\n",&iconf[i]);
       printf("conf %d read %ld\n",iconf[i],i); fflush(stdout);
   

      count=0;
       for(j=0;j<size/2;j++){
           fscanf(f,"%lf    %lf\n",&data[i][0][count],&data[i][1][count]);
          // printf("%.15g    %.15g\n",data[i][0][count],data[i][1][count]);
           count++;
              
       }
printf("readed\n");
       count=0;
       for(j=0;j<(size/2);j+=si*file_head.l0)
           for(t=0;t<file_head.l0;t++)
               for(ii=0;ii<si;ii++){
		   index=t+(ii*file_head.l0)+j;
                   out[i][count*2]=data[i][0][index];
                   out[i][count*2+1]=data[i][1][index];
                   count++;
               }
printf("converted\n");
       free(data[i][0]);
       free(data[i][1]);
       free(data[i]);

       fwrite(&iconf[i],sizeof(int),1,f1);
       fwrite(out[i],sizeof(double),size,f1);
printf("written\n");
if (i==0)
for(j=0;j<62;j++)
printf("%d   %.15g\n",j, out[0][j*2] );   
       free(out[i]);

   }   

fclose(f1);
   fclose(f);

 

  
   t2=clock();
   double dt=(double)(t2-t1)/(double)(CLOCKS_PER_SEC);
   printf("./convert took %f to convert from ASCI to binary \n",dt);
return 0;        
  
}
