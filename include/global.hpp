#ifndef global_H
#define global_H

#include <vector>

#ifdef CONTROL 
#define EXTERN 
#else
#define EXTERN extern
#endif

#define Nbootstrap 400
#define NAMESIZE  1000

//#define ensembles_reph 11
#define ensembles_reph 11
#define ensembles 10
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
#define err_MpiMeV 0.2
//#define v_MpiMeV 134.98
//#define err_MpiMeV 0.01
#define v_Mpiw0  (v_MpiMeV*v_w0MeV)

#define v_MKMeV 494.2
#define err_MKMeV 0.4
#define v_MKw0  (v_MKMeV*v_w0MeV)

#define v_MDMeV  1867.0
#define err_MDMeV  0.4

#define v_MDsMeV  1969.0
#define err_MDsMeV  0.4


#define v_fpiMeV_exp  130.4
#define err_fpiMeV_exp  0.2
//#define err_fpiMeV_exp  0.0001

#define v_fKMeV_exp  155.0 
#define err_fKMeV_exp  1.9


#define v_MOmegaMeV 1672.45 
#define err_MOmegaMeV   0.29

#define m_from2to3GEV 0.901112318


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
      int line;
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
    int Twall;
    const double eps1[4]={0,-1./1.41421356237,-1./1.41421356237,0};//1.41421356237=sqrt(2.)
    const double eps2[4]={0,1./1.41421356237,-1./1.41421356237,0};
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
    double *scalefm;
    double *s0_w0;
    double *s0_w0_cont;
    
    
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
  double (*function)(int,int,double*,int,double*);//N, Nvar, x ,Npar,P
  int N,Npar,Nvar,Njack;  
  int n;// an index to be passed to (*function)( iN ,int,double*,int,double*)
  double (*f1)(int,int,double*,int,double*);
  double (*f2)(int,int,double*,int,double*);
  int n_ext_P; //number of external parameter that will no be fitted
  double **ext_P;  //parameter to not fit, they will be add to the Nvar,  Nvar=variables+(prameter to not fit)
  int custom=false; // false=0 means default fit , 1 custom fit options
  double lambda=0.001;
  double acc=0.001; 
  double h=1e-5; 
  std::vector<int> corr_id={};
  std::vector<double>  Prange={};
  std::vector<double>  guess={};
  int t0_GEVP=3;
  int devorder=4;  // 2 , 4 , -2 adaptive h=h*parameter
  int repeat_start=1;
  
  bool plateaux_scan=false;
  FILE *f_plateaux_scan=NULL;
  std::string name_plateaux_scan;
  
  std::string resampling;
  double chi2_gap_jackboot=1;
  int guess_per_jack=3;
  int precision_sum=0;// 0 float, 1 float kahan,  >1 long double
  int verbosity=0;
  
  void restore_default(){
     custom=false; // 1 means default fit , 0 custom fit options
     lambda=0.001;
     acc=0.001; 
     h=1e-5; 
     Prange=std::vector<double>();
     guess=std::vector<double>();
     corr_id=std::vector<int>();
     t0_GEVP=3;
     devorder=4; 
     repeat_start=1;
     
     plateaux_scan=false;
     guess_per_jack=0;
     chi2_gap_jackboot=1;
     
     precision_sum=0;
     verbosity=0;
//      f_plateaux_scan=NULL;
//      name_plateaux_scan="\0";
  };
  
  
};

struct fit_result
{
  int Njack;
  double **P;// [par][jack]
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


struct store_fit_clover{
    std::string M;
    std::string GF;
    double *jack_m;
    double *jack_f;
    double *jack_B;
    double *chi2;
    double *Sigma_13;
    
    double *ms;
    double *ms_mud;
    double *fk;
    double *fk_fpi;
    
    double *mc;
    double *mc_ms;
    double *fD;
    double *fD_fk;
    
    double *fDs;
    double *fDs_fD;
    
    double *w0;
    
    int Njack;
    std::string name;
};

EXTERN struct header file_head;
EXTERN struct database_file_jack file_jack;//used to save the jack
EXTERN struct data_jack *dataJ,*gjack;
EXTERN struct result_jack result;
EXTERN int **index_a;
EXTERN int corr_counter;
EXTERN fit_type default_fit_info;
#endif
