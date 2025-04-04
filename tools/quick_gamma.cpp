#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "mutils.hpp"
#include "global.hpp"
#include "gamma_analysis.hpp"

using namespace std;

double *function(int var, int order,int flow ,double *ah){
    double *r=(double*) calloc((1),sizeof(double)); 
    r[0]=ah[1]-ah[0]*ah[0] ;
    return r;
}

int main(int argc, char **argv){
    
    error(argc!=2,1,"main ", "usage:./quinck_mass  file");
    char namefile[NAMESIZE];
    mysprintf(namefile,NAMESIZE,"%s",argv[1]); 
    FILE *infile=open_file(namefile,"r");     
    
    int count=0;
    string line;
    ifstream file(argv[1]);
    while (getline(file, line))
        count++;
 
    cout << "Numbers of lines in the file : " << count << endl;
    int confs=count;
    double *data;
    data=(double *) malloc(sizeof(double) * confs * 2);
    
    for (int iconf=0; iconf < confs ;iconf++){
        fscanf(infile,"%lf %lf ", &data[iconf*2], &data[ 1+iconf*2]);
        data[ 1+iconf*2]=data[iconf*2]*data[iconf*2];
    }
    
    double *r=analysis_gamma(  2, 1, confs, 0, data  , function);
    for (int i =0; i< 5; i++)
        printf("%g\t", r[i] );
    cout << endl;
    
}
