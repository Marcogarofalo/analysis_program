#ifndef header_BSM_H
#define header_BSM_H


#include <array>
#include <cstring> 
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>


struct  header_BSM
{
    int T;
    int L;
    double rho;
    double eta;
    double csw;
    double mu03;
    double m0;
    int confs;
    double Njack;
    double Nobs;
    
} ;

void read_header_BSM(header_BSM &head, FILE *stream){
    
    int count_lines = 0;
    char chr = getc(stream);
    while (chr != EOF)
    {
        //Count whenever new line is encountered
        if (chr == '\n')
        {
            count_lines = count_lines + 1;
        }
        //take next character from file.
        chr = getc(stream);
    }
    rewind(stream);
    
    fscanf(stream,"T  %d\n",&head.T);
    fscanf(stream,"L  %d\n",&head.L);
    fscanf(stream,"rho  %lf\n",&head.rho);
    fscanf(stream,"eta  %lf\n",&head.eta);
    fscanf(stream,"csw  %lf\n",&head.csw);
    fscanf(stream,"mu03  %lf\n",&head.mu03);
    fscanf(stream,"m0  %lf\n",&head.m0);
    head.confs=(count_lines-7)/head.T;
    error((count_lines-7)%head.T!=0,1,"read_header_BSM","numer of line= %d - header_lines=%d is not a multiple of T=%d",count_lines,7,head.T);
    
}
void check_header_BSM(header_BSM head, FILE *stream, std::string namefile){
    int count_lines = 0;
    char chr = getc(stream);
    while (chr != EOF)
    {
        //Count whenever new line is encountered
        if (chr == '\n')
        {
            count_lines = count_lines + 1;
        }
        //take next character from file.
        chr = getc(stream);
    }
    rewind(stream);
    error((count_lines-7)%head.T!=0,1,"header","numer of line= %d - header_lines=%d is not a multiple of T=%d\n file:%s",count_lines, 7,head.T,namefile.c_str());
    if (head.confs>(count_lines-7)/head.T){
        printf("file: %s \n confs does not match ref=%d read=%d  getting the minimum of the two\n",namefile.c_str(),head.confs,(count_lines-7)/head.T);
        head.confs=(count_lines-7)/head.T;
    }
    
//     error(head.confs!=(count_lines-7)/head.T,1,"header","file: %s \n confs does not match ref=%d read=%d ",namefile.c_str(),head.confs,(count_lines-7)/head.T);
    
    int tmp;
    fscanf(stream,"T  %d\n",&tmp); error(tmp!=head.T,1,"check header ","file: %s \n T does not match ref=%d read=%d ",namefile.c_str(),head.T,tmp);
    fscanf(stream,"L  %d\n",&tmp); error(tmp!=head.L,1,"check header ","file: %s \n L does not match ref=%d read=%d ",namefile.c_str(),head.L,tmp);
    double tmp1;
    fscanf(stream,"rho  %lf\n",&tmp1); error(fabs(tmp1-head.rho)>1e-6,1,"check header ","file: %s \n rho does not match ref=%f read=%f ",namefile.c_str(),head.rho,tmp1);
    fscanf(stream,"eta  %lf\n",&tmp1); error(fabs(tmp1-head.eta)>1e-6,1,"check header ","file: %s \n rho does not match ref=%f read=%f ",namefile.c_str(),head.eta,tmp1);
    fscanf(stream,"csw  %lf\n",&tmp1); error(fabs(tmp1-head.csw)>1e-6,1,"check header ","file: %s \n rho does not match ref=%f read=%f ",namefile.c_str(),head.csw,tmp1);
    fscanf(stream,"mu03  %lf\n",&tmp1); error(fabs(tmp1-head.mu03)>1e-6,1,"check header ","file: %s \n rho does not match ref=%f read=%f ",namefile.c_str(),head.mu03,tmp1);
    fscanf(stream,"m0  %lf\n",&tmp1); error(fabs(tmp1-head.m0)>1e-6,1,"check header ","file: %s \n rho does not match ref=%f read=%f ",namefile.c_str(),head.m0,tmp1);
}

void write_header_BSM_bin(header_BSM head, FILE *stream){
    fwrite(&head.T,sizeof(int),1,stream);
    fwrite(&head.L,sizeof(int),1,stream);
    fwrite(&head.rho,sizeof(double),1,stream);
    fwrite(&head.eta,sizeof(double),1,stream);
    fwrite(&head.csw,sizeof(double),1,stream);
    fwrite(&head.mu03,sizeof(double),1,stream);
    fwrite(&head.m0,sizeof(double),1,stream);
    fwrite(&head.Njack,sizeof(double),1,stream);
    
}


void read_header_BSM_bin(header_BSM &head, FILE *stream){
    fread(&head.T,sizeof(int),1,stream);
    fread(&head.L,sizeof(int),1,stream);
    fread(&head.rho,sizeof(double),1,stream);
    fread(&head.eta,sizeof(double),1,stream);
    fread(&head.csw,sizeof(double),1,stream);
    fread(&head.mu03,sizeof(double),1,stream);
    fread(&head.m0,sizeof(double),1,stream);
    fread(&head.Njack,sizeof(double),1,stream);
    
}

#endif
