#ifndef global_barion_QED_H
#define global_barion_QED_H
 

struct  database_file_barion_QED_jack
{
    int Njack;    
    double a;
    
    char M_PS[NAMESIZE];
    FILE *f_M_PS;
    
    char Zf_PS[NAMESIZE];
    FILE *f_Zf_PS;
    
    char m_PCAC[NAMESIZE];
    FILE *f_m_PCAC;
    

};


#endif
