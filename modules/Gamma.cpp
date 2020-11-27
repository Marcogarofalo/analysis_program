#define GAMMA_C

#include <iostream>
#include <cmath>

static double order=0;
static double obs[1],dobs[1],ddobs[1],taubb_intF[1],dtau[1];
static double *abb;
static double CbbF[1],w[1];
static double **gammaFbb;

double **Gamma(int t, int var, int order, int rep,int nconf, double *a, double *bba)
{
  double **r;
  int i0,i1,i2,i3,alpha,N;
  
  alpha=(order+1)*var;
  N=alpha*rep;
  r=(double**) malloc(sizeof(double*)*(alpha));
  for(i0=0;i0<alpha;i0++)
    r[i0]=(double*) calloc(alpha,sizeof(double));
  
  for(i0=0;i0<alpha;i0++)
    for(i1=0;i1<alpha;i1++)
        for(i2=0;i2<rep;i2++)
            for(i3=0;i3<(nconf-t);i3++)
                r[i0][i1]+= (a[i0+i2*alpha+i3*N]-bba[i0])*(a[i1+i2*alpha+(i3+t)*N]-bba[i1]);
    // if(t==0) printf("r[0][0]=%f\n",r[0][0]);  
          
    
          /*printf("a-bba=%f\n",a[0]-bba[0]);*/
for(i0=0;i0<alpha;i0++)
  for(i1=0;i1<alpha;i1++)
	    r[i0][i1]/=((double)(rep*nconf-rep*t));
     
return r;
}




double **barf( int var, int order, int rep,int nconf,int flow, double *bba,double **ga, double *function(int var, int order,int flow ,double *ah))
{
double **r,*tmp,*tmp1;
double *h,*ah;
int i,j,N,alpha;

N=rep*nconf;
alpha=(order+1)*var;
ah=(double*) malloc(alpha*sizeof(double));
h=(double*) malloc(alpha*sizeof(double));
r=(double**) malloc(alpha*sizeof(double*));
for(i=0;i<alpha;i++)
{
  r[i]=(double*) calloc(order+1,sizeof(double));
  h[i]=sqrt(   ga[i][i]/((double)N*4.)  );
  ah[i]=bba[i];
}

for(i=0;i<alpha;i++)
{
  ah[i]+=h[i];
  tmp1=function(var,order,flow,ah);
  ah[i]-=h[i];ah[i]-=h[i]; 
  tmp=function(var,order,flow,ah);
  //sub_pseries(order,tmp1,tmp,r[i]);
  for(j=0;j<=order;j++){
      r[i][j]=tmp1[j]-tmp[j];
      r[i][j]/=2.*h[i];
  }
  free(tmp);free(tmp1);ah[i]+=h[i];
  //scale_pseries(order,0,1./(2.*h[i]),r[i]  );
}

free(h);free(ah);
return r;
}

// *a is an array of data a[ a+ rep *alpha + conf* Nv ] 
// a =0 ... (alpha-1)
//alpha =(order+1)*var
void mean_value(int var, int order,int rep, int nconf,int flow,double *a, double *function(int var, int order,int flow ,double *ah))
{
    int i0,i1,i2,j,alpha,imax;
    double **ab,N;
    double *Fbb,*Fb,*tmp;

    alpha=(order+1)*var;
    N=nconf*rep;
    imax=alpha*rep;    

    for(i0=0;i0<alpha;i0++)
       abb[i0]=0;
        int i,Nfit=1;
// printf("a=%f\n",a[1+38*2]); 
    ab=(double**) malloc(rep*sizeof(double*));
    for(i1=0;i1<rep;i1++)
        ab[i1]=(double*) calloc(alpha,sizeof(double));
    
    Fb=(double*) calloc(order+1,sizeof(double));
    
    for(i0=0;i0<alpha;i0++)
        for(i1=0;i1<rep;i1++)
            for(i2=0;i2<nconf;i2++)
	     ab[i1][i0]+=a[i0+i1*alpha+i2*imax];   
    //printf("nconf=%d\t replicas=%d\t alpha=%d\n",nconf,rep,alpha);
   
    for(i0=0;i0<alpha;i0++)
        for(i1=0;i1<rep;i1++)
            abb[i0]+=ab[i1][i0];

    for(i0=0;i0<alpha;i0++)
        for(i1=0;i1<rep;i1++)
            ab[i1][i0]/=(double)nconf;

    for(i0=0;i0<alpha;i0++)        
	abb[i0]/=(double) (((double) nconf)*rep); 
    
   
    Fbb=function(var,order,flow,abb);    
   /* printf("Fbb=%f\n",Fbb[0]);*/
    if(rep==1)   for(i0=0;i0<=order;i0++) obs[i0]=Fbb[i0];
    else
    {    
        for(i1=0;i1<rep;i1++)
        {
            printf(" ok untill here  replica %d\n",i1);
            tmp=function(var,order,flow,ab[i1]);
            for(j=0;j<=order;j++){
                tmp[j]/=(double)nconf;
                Fb[j]+=tmp[j];
            }
            //scale_pseries(order,0,nconf,tmp);
            //add_pseries(order,tmp,Fb,Fb);
        }
        for(j=0;j<=order;j++)
            Fb[j]/=(double) N;
        //scale_pseries(order,0,1./((double)N),Fb);
        
        for(i0=0;i0<=order;i0++)
            obs[i0]=(((double)rep)*Fbb[i0]-Fb[i0])/(((double)rep)-1.);
    }
        
    //free_dpseries(rep-1,ab);
    for(j=0;j<=order;j++){
    free(ab[j]);
    }
}

double *Gammaf( int var, int order,double **ga,double **fa)
{
int i,j,k,alpha;
double *r;
r=(double*) calloc(order+1,sizeof(double));
alpha=var*(order+1);

for(i=0;i<=order;i++)
  for(j=0;j<alpha;j++)
    for(k=0;k<alpha;k++)
      r[i]+=fa[j][i]*fa[k][i]*ga[j][k];
   
return r;
}

void windowing(int var,int order, int rep, int nconf, int flow,double *a, double *bba ,double *function(int var, int order,int flow ,double *ah))
{
    double **fbba,**tmp,*g,*tau,Caa=0;
    int count=0, i,j,i1,N,alpha;
    double S=1.5;
    
    FILE  *file_tau=NULL;
    file_tau=fopen("tau_int","w+");
    if(file_tau==NULL){
	printf("unable to open analysis file"); exit(0);
    }
    
    
    alpha=(order+1)*var;
    
    g=(double*) calloc(order+1,sizeof(double));
    tau=(double*) calloc(order+1,sizeof(double));
    
    N=rep*nconf;
      
     for(i=0;i<=order;i++)
    {
        CbbF[i]=0;
        w[i]=-1;
    }

    tmp=Gamma(0,  var,  order,  rep, nconf, a,  bba);
    fbba=barf(  var,  order,  rep, nconf, flow, bba, tmp ,function);
    gammaFbb[0]=Gammaf(var,order,tmp,fbba);

    for(i1=0;i1<=order;i1++)
        CbbF[i1]+=gammaFbb[0][i1];
    Caa+=tmp[0][0];
 
    for(i1=0;i1<alpha;i1++)
        free(tmp[i1]);

    fprintf(file_tau,"%d   \t %d %g\n",flow,0,0.5);
    
    for(i=1;i<nconf;i++)
    {   
      
        tmp=Gamma(i,  var,  order,  rep, nconf, a,  bba);
        gammaFbb[i]=Gammaf(var,order,tmp,fbba);
        if(i==1) for(j=0;j<=order;j++)
        {
            if(gammaFbb[1][j]<0){
                /*printf("there is no autocorrelation\n");*/
                w[j]=0;count++;
            }
        }  
        if(w[0]==-1)  Caa+=2*tmp[0][0];
        //free_dpseries(alpha-1,tmp);
        for(i1=0;i1<alpha;i1++)
            free(tmp[i1]);
        for(j=0;j<=order;j++)
        {
            
            if(w[j]==-1)
            {
                CbbF[j]+=2.*gammaFbb[i][j];
                taubb_intF[j]=CbbF[j]/(2.*gammaFbb[0][j]);
 		        fprintf(file_tau,"%d   \t %d %g\n",flow,i,taubb_intF[j]);
        	    tau[j]=0.6;
	            if(taubb_intF[j]>0.5)
                        tau[j]=S/(  log( (2.*taubb_intF[j]+1.)/(2.*taubb_intF[j]-1.)  ));
                g[j]=exp(-((double)i)/tau[j])- (tau[j]/ (sqrt((double)(i*N) ))  );
                /*printf("g=%f",g[j]);*/
                if(g[j]<0)
                {  count++;  w[j]=i; }
                
            }
            /*if(j==0) printf("gammaFbb[%d]=%0.10f\n",i,gammaFbb[i][0]);*/
             if(count==order+1) break;
        }
        free(gammaFbb[i]);
        if(count==order+1) break;
    }
    
            

    //free_dpseries(alpha-1,fbba);
    for(i1=0;i1<alpha;i1++)
        free(fbba[i1]);

    free(g);free(tau);

        for(j=0;j<=order;j++)
        {   
            if(w[j]==-1)
            {
		printf("Windowing condition order %d failed up to W = %d\n",j,nconf-1);
                w[j]=nconf-1;
	    }
        }


    for(j=0;j<=order;j++)
    {   
        
        gammaFbb[0][j]+=CbbF[j]/((double)N);
        CbbF[j]+=CbbF[j]*(2.*w[j]+1)/((double)N);
        taubb_intF[j]=CbbF[j]/(2.*gammaFbb[0][j]);
    }
    free(gammaFbb[0]);
    fprintf(file_tau,"\n");
    fclose(file_tau);
    
}


void return_answer( int var, int order ,int rep, int nconf)
{
    int i,N;
    
    
    N=rep*nconf;
    for(i=0;i<=order;i++)
    {
        dobs[i]=CbbF[i]/((double)N);
        dobs[i]=sqrt(dobs[i]);
       
        ddobs[i]=dobs[i]*sqrt((w[i]+0.5)/N);
        dtau[i]=sqrt( (w[i]+0.5-taubb_intF[i])/ ((double)N) )*2.* taubb_intF[i] ;
    }
    
}

//replicas =rep 
// flow is an extra parameter to pass to the function
double  *analysis_gamma (  int var, int rep, int nconf,int flow, double *a , double *function(int var, int order,int flow ,double *ah)){
    double *r=(double*) malloc(sizeof(double) *5 );
    
    abb=(double*) calloc(var,sizeof(double));
    gammaFbb=(double**) malloc(sizeof(double*)*nconf);
    
    
    mean_value(var,order,rep, nconf,flow, a, function);
    windowing(var,order,rep, nconf,flow,a,abb, function) ;
    return_answer( var,order, rep, nconf);//printf("HERE  %d %d %f\n",k,i+L0,abb[2*39]);
    //fprintf(file_out,"%d \t %.15f  %.15f %.15f %.15f %.15f \n",i,obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);
    
    
    r[0]=obs[0];
    r[1]=dobs[0]; 
    r[2]=ddobs[0];
    r[3]=taubb_intF[0];
    r[4]=dtau[0];
    
    free(abb);
    free(gammaFbb);
    return r;
}
