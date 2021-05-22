#define header_phi4_C

#include "header_phi4.hpp"

#include <array>
#include <cstring> 
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>

void read_header_phi4(FILE *stream, cluster::IO_params &params ){
    size_t i;  
    i=fread(&params.data.L[0], sizeof(int), 1, stream); 
    i=fread(&params.data.L[1], sizeof(int), 1, stream); 
    i=fread(&params.data.L[2], sizeof(int), 1, stream); 
    i=fread(&params.data.L[3], sizeof(int), 1, stream); 
    printf("%d %d %d %d\n",params.data.L[0],params.data.L[1],params.data.L[2],params.data.L[3]);
    char string[100];
    i=fread(string, sizeof(char)*100, 1, stream); 
    std::string tmp(string);
    params.data.formulation=tmp;
    
    i=fread(&params.data.msq0, sizeof(double), 1, stream); 
    i=fread(&params.data.msq1, sizeof(double), 1, stream); 
    i=fread(&params.data.lambdaC0, sizeof(double), 1, stream); 
    i=fread(&params.data.lambdaC1, sizeof(double), 1, stream); 
    i=fread(&params.data.muC, sizeof(double), 1, stream); 
    i=fread(&params.data.gC, sizeof(double), 1, stream); 
    
    i=fread(&params.data.metropolis_local_hits, sizeof(int), 1, stream); 
    i=fread(&params.data.metropolis_global_hits, sizeof(int), 1, stream); 
    i=fread(&params.data.metropolis_delta, sizeof(double), 1, stream); 
    
    i=fread(&params.data.cluster_hits, sizeof(int), 1, stream); 
    i=fread(&params.data.cluster_min_size, sizeof(double), 1, stream); 
    
    i=fread(&params.data.seed, sizeof(int), 1, stream); 
    i=fread(&params.data.replica, sizeof(int), 1, stream); 
    
    
    i=fread(&params.data.ncorr, sizeof(int), 1, stream); 
    //printf("correlators=%d\n",params.data.ncorr);
    
    i=fread(&params.data.size, sizeof(size_t), 1, stream); 
    //printf("size=%ld\n",params.data.size);
    
    params.data.header_size=ftell(stream);
    //printf("header size=%d\n",params.data.header_size);
}


void write_header_phi4(FILE *stream, cluster::IO_params params ){
    
    fwrite(&params.data.L[0], sizeof(int), 1, stream); 
    fwrite(&params.data.L[1], sizeof(int), 1, stream); 
    fwrite(&params.data.L[2], sizeof(int), 1, stream); 
    fwrite(&params.data.L[3], sizeof(int), 1, stream); 
    char string[100];
    mysprintf(string,100,"%s0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",params.data.formulation.c_str());
    fwrite(string, sizeof(char)*100, 1, stream); 
    //params.data.formulation = std::to_string(string);
    
    fwrite(&params.data.msq0, sizeof(double), 1, stream); 
    fwrite(&params.data.msq1, sizeof(double), 1, stream); 
    fwrite(&params.data.lambdaC0, sizeof(double), 1, stream); 
    fwrite(&params.data.lambdaC1, sizeof(double), 1, stream); 
    fwrite(&params.data.muC, sizeof(double), 1, stream); 
    fwrite(&params.data.gC, sizeof(double), 1, stream); 
    
    fwrite(&params.data.metropolis_local_hits, sizeof(int), 1, stream); 
    fwrite(&params.data.metropolis_global_hits, sizeof(int), 1, stream); 
    fwrite(&params.data.metropolis_delta, sizeof(double), 1, stream); 
    
    fwrite(&params.data.cluster_hits, sizeof(int), 1, stream); 
    fwrite(&params.data.cluster_min_size, sizeof(double), 1, stream); 
    
    fwrite(&params.data.seed, sizeof(int), 1, stream); 
    fwrite(&params.data.replica, sizeof(int), 1, stream); 
    
    
    fwrite(&params.data.ncorr, sizeof(int), 1, stream); 
    
    fwrite(&params.data.size, sizeof(size_t), 1, stream);     
}



void read_Njack_Nobs( FILE *stream, cluster::IO_params params, int &Njack, int &Nobs ){
    
    long int tmp;
    int s=params.data.header_size;
    
    size_t i=fread(&Njack, sizeof(int), 1, stream );
    
    
    fseek(stream, 0, SEEK_END);
    tmp = ftell(stream);
    tmp-= params.data.header_size+sizeof(int) ;
    
    s=Njack;
    
    Nobs= (tmp)/ ((s)*sizeof(double) );
    
    fseek(stream, params.data.header_size+sizeof(int), SEEK_SET);
    
    
    
}

void read_dataj(FILE *stream,cluster::IO_params params, data_phi &dj){
    read_Njack_Nobs(stream, params, dj.Njack, dj.Nobs );
    //printf("Nobs=%d   Njack=%d\n",dj.Nobs,dj.Njack)
    dj.jack=double_malloc_2( dj.Nobs, dj.Njack);
    size_t i=0;
    for (int obs=0; obs<dj.Nobs; obs++ ){
        i+=fread(dj.jack[obs], sizeof(double ), dj.Njack, stream );
    }
    
}

void emplace_back_par_data( char *namefile , vector<cluster::IO_params> &paramsj, vector<data_phi> &dataj){
    cluster::IO_params params;
    data_phi  data;
    FILE *f=open_file(namefile,"r");
    read_header_phi4( f  , params);
    read_dataj(f,params,data );
    fclose(f);
    
    paramsj.emplace_back(params);
    dataj.emplace_back(data);
    //printf("E1=%g    %g\n",dataj[0].jack[1][   dataj[0].Njack-1 ],    data.jack[1][data.Njack-1]);
    printf("ending the function: emplace_back_par_data\n");
}


vector<data_phi> create_generalised_resampling(  vector<data_phi> &dataj ){
    // if the length is the same return dataj
    int same=0;
    for( auto &d :dataj){
        printf("jacks=%d\n",d.Njack);
        if (d.Njack==dataj[0].Njack)
            same++;
    }
    if (same==dataj.size()){
        cout << "all the files have the same number of jack/boot , do nothing"<<endl;
        return dataj;
    }
    else{
        cout << "creating generalised jack"<<endl;
        vector<data_phi> gjack;
        //gjack.resize(dataj.size());
        //jac_tot is the summ of all jackknife +1 
        //remember alle the dataj have one extra entry for the mean
        int jack_tot=0;
        for( int e=0 ; e<dataj.size();e++)
            jack_tot+=dataj[e].Njack;
        jack_tot=jack_tot-dataj.size()+1;
        cout << "ensembles "<< dataj.size() << endl;
        cout<< "jack tot= "<< jack_tot<< endl;
        
        //get Nobs the minimum number of observable between the diles
        int Nobs=1000;
        for( auto &d :dataj)
            if(Nobs> d.Nobs )
                Nobs=d.Nobs;
            
            
            
            
            for(int e=0;e<dataj.size();e++){
                data_phi tmp;
                tmp.Njack=jack_tot;
                tmp.Nobs=Nobs;
                cout<< Nobs<<" " <<jack_tot<< endl;
                tmp.jack=double_malloc_2( Nobs, jack_tot);
                int counter=0;
                for(int e1=0;e1<dataj.size();e1++){
                    for(int j=0;j<(dataj[e1].Njack-1);j++){
                        for(int o=0;o<Nobs;o++){ 
                            if (e==e1){
                                tmp.jack[o][j+counter]=dataj[e].jack[o][j];
                            }
                            else{
                                tmp.jack[o][j+counter]=dataj[e].jack[o][ dataj[e].Njack-1  ];
                            }
                        }
                    }
                    counter+=dataj[e1].Njack-1 ;
                }
                for(int o=0;o<Nobs;o++)
                    tmp.jack[o][ jack_tot-1 ]=dataj[e].jack[o][ dataj[e].Njack-1  ];
                
                free_2(dataj[e].Nobs, dataj[e].jack);
                gjack.emplace_back(tmp);
            }
            vector<data_phi>().swap(dataj);
            
            return gjack;
            
    }
    
}
