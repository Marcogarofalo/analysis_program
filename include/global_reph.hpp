#ifndef global_reph_H
#define global_reph_H



struct  database_file_reph_jack
{
    char FA[NAMESIZE],M_PS[NAMESIZE];
    char FAp[NAMESIZE];
    char FV[NAMESIZE];
    char xG[NAMESIZE];
    char f_PS[NAMESIZE];
    char f_PS_ls_ss[NAMESIZE];
    char FAp1[NAMESIZE];
    char FV_exclude[NAMESIZE];
    char FAp1_exclude[NAMESIZE];
    char kp[NAMESIZE];
    int Njack;    
    double a;
    FILE *f_FA,*f_M_PS;
    FILE *f_FAp;
    FILE *f_FV;
    FILE *f_xG;
    
    char FA_from_H0[NAMESIZE];
    FILE *f_FA_from_H0;
    
    char FV_from_H0[NAMESIZE];
    FILE *f_FV_from_H0;
    
    char FV_from_H0_HA[NAMESIZE];
    FILE *f_FV_from_H0_HA;
    
    char Zf_PS[NAMESIZE];
    FILE *f_Zf_PS;
    
   

};




struct  reph_jack
{
    ///rephJ[ensemble].observable[mass_index][xG][jack]
    double **FAp,**FV,**xG,**FA,**M_PS,**FA_from_H0,**Zf_PS,**FV_from_H0,**FV_from_H0_HA;
    
     ///rephJ[ensemble].observable[jack]
    double *w0,*Zp,*ZA,*ZV,*ZAV;
    
    ///rephJ[ensemble].observable[mass_index]
    double *KM,*Kf;
} ;


struct ik_for_n
{
  int ikt;
  int iks;  
};
#endif
 
