#ifndef global_H
#define global_H



#ifdef CONTROL 
#define EXTERN 
#else
#define EXTERN extern
#endif

#define Nbootstrap 200
#define NAMESIZE  1000

//#define ensembles_reph 11
#define ensembles_reph 11
#define ensembles 9
#define pi_greco 3.141592653589793
#define pdiffmax  3

#define v_w0fm 0.1714  // MILC
#define err_w0fm 0.0015  // MILC
//#define w0fm 0.1755  //BMW 

//#define v_w0fm 0.172  // ETMC
//#define err_w0fm 0.000001  // ETMC

#define v_w0GeV ((v_w0fm/197.326963)*1000)
#define v_w0MeV (v_w0fm/197.326963)

#define v_MpiMeV 134.80
#define err_MpiMeV 0.004
#define v_Mpiw0  (v_MpiMeV*v_w0MeV)

#define v_MKMeV 494.2
#define err_MKMeV 0.4
#define v_MKw0  (v_MKMeV*v_w0MeV)

#define v_MDMeV  1867.0
#define err_MDMeV  0.4

#define v_MDsMeV  1969.0
#define err_MDsMeV  0.4


#define v_fpiMeV_exp  130.41
#define err_fpiMeV_exp  0.0001

#define v_fKMeV_exp  155.0 
#define err_fKMeV_exp  1.9


#define v_MOmegaMeV 1672.45 
#define err_MOmegaMeV   0.29

EXTERN int jack_tot;


struct  header
{
    int twist;
    int nf;
    int nsrc,nsrcd;
    int l0,l1,l2,l3;
    int nk,nmoms;
    double beta,ksea,musea,csw;
    double *k,*mu,**mom;
    
    int allocated=0;
} ;
 

struct  kinematic
{
      double k2,k1,mom2,mom1,mom02,mom01;
      int ik1,ik2;
      double Mom1[4],Mom2[4];
      int r2,r1;
      char   plateau_m_ll[NAMESIZE];
      char   plateau_m_sl[NAMESIZE];
      char   plateau_f[NAMESIZE];
      char   plateau_m_GEVP[NAMESIZE];
} ;

struct  kinematic_G
{
    double k0,kt,ks;
    int ik0,ikt,iks;
    double Mom0[4],Momt[4],Moms[4];
    double p[4],k[4];
    int rt,rs,i;
    double  E_g,E_gT;
    double kp;
    const double eps1[4]={0,-1./sqrt(2.),-1./sqrt(2.),0};
    const double eps2[4]={0,1./sqrt(2.),-1./sqrt(2.),0};
    double eps1_curl_p[4],eps2_curl_p[4];
    double eps1_curl_k[4],eps2_curl_k[4];
    char plateau_RA[NAMESIZE];
    char plateau_RV[NAMESIZE];
    char plateau_H_H0_A[NAMESIZE];
    
} ;


struct  database_file_jack
{
    char M_PS_GEVP[NAMESIZE];
    char M_PS[NAMESIZE];
    char f_PS[NAMESIZE];
    char Zf_PS[NAMESIZE];
    char f_PS_ls_ss[NAMESIZE];
    char FV[NAMESIZE],FV_autoplateaux[NAMESIZE];
    char FAp1[NAMESIZE];
    char FA[NAMESIZE],FA_autoplateaux[NAMESIZE];
    char FAp[NAMESIZE];
    char xG[NAMESIZE];
    char FV_exclude[NAMESIZE];
    char FAp1_exclude[NAMESIZE];
    char kp[NAMESIZE];
    char sampling[10];
    int Njack;    
    double a;
    FILE *f_M_PS;
    FILE *f_f_PS,*f_Zf_PS;
    FILE *f_M_PS_GEVP;
    FILE *f_f_PS_ls_ss;
    FILE *f_FAp, *f_FA, *f_FA_autoplateaux;
    FILE *f_FV,*f_FV_autoplateaux;
    FILE *f_xG;
    
    char FA_from_H0[NAMESIZE];
    char FA_from_H0_autoplateaux[NAMESIZE];
    char FV_from_H0[NAMESIZE];
    char FV_from_H0_autoplateaux[NAMESIZE];
    char FV_from_H0_HA[NAMESIZE];

};


struct  data_jack
{
    ///dataJ[ensemble].observable[mass_index][jack]
    double **M_PS_GEVP_jack,**f_PS_jack;
    double **M_PS_jack,**f_PS_ls_ss_jack;
    double **M_Omega_jack;
    
     ///dataJ[ensemble].observable[jack]
    double *w0,*Zp,*Zv;
    
     ///dataJ[ensemble].observable[mass_index]
    double *KM,*Kf;
};


struct  result_jack
{
    ///dataJ[ensemble].observable[mass_index][jack]
    double *Bw,*fw,*l3b,*l4b;
    double *Bw_from_M,*fw_from_M,*l4b_from_M;
    
     ///dataJ[ensemble].observable[jack]
    double *mlw,*fpiw ;
    double *fpiw_from_M,*fKw_from_M,*fDw_from_M,*fDsw_from_M,*fDs_fD_from_M;
    
    double  *msw,*fkw,*ms_mud,*fk_fpi,*fk_fpi_from_M;
    double *mcw,*fDw,*fDsw,*mc_ms;
    
    double  *w0fm, *w0MeV,*MpiMeV,*MKMeV,*MDMeV,*MDsMeV;
    double **PMPi,***PMK,***PMD;
    
    double *fpiMeV_exp;
    double *fKMeV_exp;
    
    double *MOmegaMeV;
    
};




struct  jack_fit
{
  double **m,**oPp;
  double **RA,**RV,**HA,**HV,**RA_autoplateaux,**RV_autoplateaux;
  double **f_PS,**Zf_PS;  
    
};

struct  fit_type
{
  double (*function)(int,int,double*,int,double*);
  int N,Npar,Nvar;  
  double (*f1)(int,int,double*,int,double*);
  double (*f2)(int,int,double*,int,double*);
};

struct fit_result
{
  int Njack;
  double **P;
  double *chi2;
  double ***C;   
    
};

struct fit_all
{
    int Nfits;
    struct fit_type  *info;
    struct fit_result  *out;
};

struct observable
{
   char name_out[NAMESIZE]; 
   FILE *f_out;
   
   char name_plateaux[NAMESIZE];
   FILE *f_plateaux;
   
   char name_jack[NAMESIZE];
   FILE *f_jack;
   
    
};

EXTERN struct header file_head;
EXTERN struct database_file_jack file_jack;//used to save the jack
EXTERN struct data_jack *dataJ,*gjack;
EXTERN struct result_jack result;
EXTERN int **index_a;
#endif
