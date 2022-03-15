#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"
#include "resampling.hpp"
#include "read.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"
#include <string>
#include "correlators_analysis.hpp"
#include "eigensystem.hpp"
#include "header_phi4.hpp"
#include "zeta_interpolation.hpp"
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>


//local folder
#include "mass_phi4.hpp"
#include "fit_function.hpp"

#include <QC3_interface.hpp>



using namespace std;
int Ne=0;


void print_fit_band_L_M(char **argv,vector<data_phi> gjack ,struct fit_type fit_info,struct fit_type fit_info_m0 , const char* label, struct fit_result fit_out, struct fit_result fit_out_m0,    vector<cluster::IO_params> params, std::vector<int> myen, std::vector<int> Lrange={16,50}){
	int Npar=fit_info.Npar;
	int Nvar=fit_info.Nvar+fit_info.n_ext_P;
	int Njack=gjack[0].Njack;
	int N=fit_info.N;
	char namefile[NAMESIZE];
	FILE *f;
	
	double **tif=swap_indices(fit_info.Npar,Njack,fit_out.P);
	double **tif_m0=swap_indices(fit_info_m0.Npar,Njack,fit_out_m0.P);
	
	for (int n=0;n< N; n++){
		
		mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_n%d_L.txt",argv[3], label,n);
		f=open_file(namefile,"w+");
		double *tmpx=(double*) malloc(sizeof(double)* Nvar);
		double *tmpy=(double*) malloc(sizeof(double)* Njack);
		printf("writing: %s\n",namefile);
		
		for (int i=Lrange[0] ; i<Lrange[1]; i++){
			double finalL=i;
			tmpx[0]=finalL;
			double *E3_m=(double*) malloc(sizeof(double) *Njack);
			for (int j=0;j<Njack;j++){
				E3_m[j]=lhs_E3_m(n,0,j,params,gjack,fit_info);
			}
			double E3_m_err=error_jackboot(argv[1], Njack, E3_m );
			for (int j=0;j<Njack;j++){
				tmpx[1]=fit_info_m0.function(0,fit_info_m0.Nvar,tmpx,fit_info_m0.Npar,tif_m0[j]); //m0   put for each n the mass of the last ensemble
				
				tmpx[2]=gjack[0].jack[2][j];//m1  //
				tmpx[3]=gjack[0].jack[4][j];//E20
				tmpx[4]=gjack[0].jack[5][j];//E21
				tmpx[5]=(double) params[0].data.L[0];//T
				tmpx[6]=1;//k
				tmpx[7]=1;//MpiL
				
				double *x=(double*) malloc(sizeof(double)*myen.size());
				double *y=(double*) malloc(sizeof(double)*myen.size());
				for (int e=0; e<myen.size();e++){
					y[e]= lhs_E3_m(n,myen[e],j,params,gjack,fit_info);// E3( \vec{n} )/mass
					if (params[0].data.gC>0){
						y[e]= lhs_E3_m_g(n,myen[e],j,params,gjack,fit_info);
					}
					x[e]=params[myen[e]].data.L[1];
				}
				tmpx[8]=inter_spline( tmpx[0], myen.size(), x, y   ); // tmpx[0]=L 

				if (params[0].data.gC>0){
					tmpx[1]=gjack[0].jack[443][j];//m0
					tmpx[10]=compute_k_m_g(n,0,j,params,gjack,fit_info);//k_m
					tmpx[6]=tmpx[10]*tmpx[1];//k
					tmpx[3]=gjack[0].jack[594][j];//E20
					tmpx[7]=gjack[0].jack[443][j]*(double) params[0].data.L[1]/(2.*pi_greco);//mL_2pi
				}    

				free(x);free(y);
				//tmpx[8]=E3_m[j];// E3( \vec{n} )/mass
				
				// as error on E3/m  we take the error on the ensemble=0
				tmpx[9]=E3_m_err;
				
				
				for(int i=fit_info.Nvar ; i<fit_info.Nvar+ fit_info.n_ext_P; i++)
					tmpx[i]=fit_info.ext_P[i-fit_info.Nvar][j];
				
				
				
				tmpy[j]=fit_info.function(n,Nvar,tmpx,Npar,tif[j]);
//                 tmpy[j]=tmpx[8];
//                 printf("L=%g  j=%d  E=%g\n",finalL,j,tmpx[8]);
//                 if (fabs(finalL-36)<1e-5  && n==2)
//                     printf("%g   j=%d \t k=%g  m=%g P0=%g  P1=%g\n",finalL,j,tmpy[j], tmpx[1], tif[j][0] ,tif[j][1]  );
//                 
			}
			/*if (fabs(finalL-36)<1e-5  && n==2){
				printf("%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
			}*/
			fprintf(f,"%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
			
			free(E3_m);
		}
		
		free(tmpy);free(tmpx);
		fclose(f); 
	}
	free_2(Njack,tif);
	free_2(Njack,tif_m0);
	
}





void print_fit_band_E3_vs_L(char **argv,vector<data_phi> gjack , struct fit_type fit_info, struct fit_type fit_info_m0 , const char* label, struct fit_result fit_out, struct fit_result fit_out_m0,    vector<cluster::IO_params> params, std::vector<int> myen,  struct fit_type fit_info_E3_poly ,fit_result fit_E3_poly, std::vector<int> Lrange={16,50}){
	int Npar=fit_info.Npar;
	int Nvar=fit_info.Nvar+fit_info.n_ext_P;
	int Njack=gjack[0].Njack;
	int N=fit_info.N;
	char namefile[NAMESIZE];
	FILE *f;
	
	double **tif=swap_indices(fit_info.Npar,Njack,fit_out.P);
	double **tif_m0=swap_indices(fit_info_m0.Npar,Njack,fit_out_m0.P);
	double **tif_E3_poly=swap_indices(fit_info_E3_poly.Npar, Njack, fit_E3_poly.P);

	
	for (int n=0;n< N; n++){
		
		mysprintf(namefile,NAMESIZE,"%s/%s_fit_out_n%d_L.txt",argv[3], label,n);
		f=open_file(namefile,"w+");
		double *tmpx=(double*) malloc(sizeof(double)* Nvar);
		double *tmpy=(double*) malloc(sizeof(double)* Njack);
		printf("writing: %s\n",namefile);
		
		for (int i=Lrange[0] ; i<Lrange[1]; i++){
			double finalL=i;
			tmpx[0]=finalL;
			double *E3_m=(double*) malloc(sizeof(double) *Njack);
			for (int j=0;j<Njack;j++){
				E3_m[j]=fit_info_E3_poly.function(n,fit_info_E3_poly.Nvar , tmpx, fit_info_E3_poly.Npar, tif_E3_poly[j]);
			}
			double E3_m_err=error_jackboot(argv[1], Njack, E3_m );
			
			
			
			for (int j=0;j<Njack;j++){
				tmpx[1]=fit_info_m0.function(0,fit_info_m0.Nvar,tmpx,fit_info_m0.Npar,tif_m0[j]); //m0   put for each n the mass of the last ensemble
				
				tmpx[2]=gjack[0].jack[2][j];//m1  //
				tmpx[3]=gjack[0].jack[4][j];//E20
				tmpx[4]=gjack[0].jack[5][j];//E21
				tmpx[5]=(double) params[0].data.L[0];//T
				tmpx[6]=1;//k
				tmpx[7]=1;//MpiL
				
//                 double *x=(double*) malloc(sizeof(double)*myen.size());
//                 double *y=(double*) malloc(sizeof(double)*myen.size());
//                 for (int e=0; e<myen.size();e++){
//                     y[e]= lhs_E3_m(n,myen[e],j,params,gjack,fit_info);// E3( \vec{n} )/mass
//                     x[e]=params[myen[e]].data.L[1];
//                 }
//                 tmpx[8]=inter_spline( tmpx[0], myen.size(), x, y   ); // tmpx[0]=L 
//                 free(x);free(y);
				
				tmpx[8]=E3_m[Njack-1];// E3( \vec{n} )/mass
				// as error on E3/m  we take the error on the ensemble=0
				tmpx[9]=E3_m_err;
				
				tmpx[11]=E3_m[Njack-1];// E3( \vec{n} )/mass
				// as error on E3/m  we take the error on the ensemble=0
				tmpx[12]=E3_m_err;
				
				if (params[0].data.gC>0){
					tmpx[1]=gjack[0].jack[443][j];//m0
					tmpx[10]=compute_k_m_g(n,0,j,params,gjack,fit_info);//k_m
					tmpx[6]=tmpx[10]*tmpx[1];//k
					tmpx[3]=gjack[0].jack[594][j];//E20
					tmpx[7]=gjack[0].jack[443][j]*(double) params[0].data.L[1]/(2.*pi_greco);//mL_2pi
				}    
				
				for(int i=fit_info.Nvar ; i<fit_info.Nvar+ fit_info.n_ext_P; i++)
					tmpx[i]=fit_info.ext_P[i-fit_info.Nvar][j];
				
				
				
				tmpy[j]=fit_info.function(n,Nvar,tmpx,Npar,tif[j]);
				//                 tmpy[j]=tmpx[8];
				//                 printf("L=%g  j=%d  E=%g\n",finalL,j,tmpx[8]);
				//                 if (fabs(finalL-36)<1e-5  && n==2)
				//                     printf("%g   j=%d \t k=%g  m=%g P0=%g  P1=%g\n",finalL,j,tmpy[j], tmpx[1], tif[j][0] ,tif[j][1]  );
				//                 
			}
			/*if (fabs(finalL-36)<1e-5  && n==2){
			 *                printf("%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
		}*/
			fprintf(f,"%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
			
			free(E3_m);
		}
		
		free(tmpy);free(tmpx);
		fclose(f); 
	}
	free_2(Njack,tif);
	free_2(Njack,tif_m0);
	free_2(Njack,tif_E3_poly);
}





void print_kiso_P0_inf_L_M(char **argv,vector<data_phi> gjack , struct fit_type fit_info, struct fit_type fit_info_m0 , const char* label, struct fit_result fit_out, struct fit_result fit_out_m0,    vector<cluster::IO_params> params, std::vector<int> myen,  struct fit_type fit_info_E3_poly ,fit_result fit_E3_poly, std::vector<int> Lrange={16,50}){
	int Npar=fit_info.Npar;
	int Nvar=fit_info.Nvar+fit_info.n_ext_P;
	int Njack=gjack[0].Njack;
	int N=fit_info.N;
	char namefile[NAMESIZE];
	FILE *f;
	
	double **tif=double_malloc_2( Njack,Npar);
	for (int i=0;i<Npar; i++)
		for (int j=0;j<Njack;j++)
			tif[j][i]=0;
	double **tif_m0=swap_indices(fit_info_m0.Npar,Njack,fit_out_m0.P);
	double **tif_E3_poly=swap_indices(fit_info_E3_poly.Npar, Njack, fit_E3_poly.P);
	
	for (int n=0;n< N; n++){
		
		mysprintf(namefile,NAMESIZE,"%s/kiso_P0_n%d_L.txt",argv[3],n);
		f=open_file(namefile,"w+");
		double *tmpx=(double*) malloc(sizeof(double)* Nvar);
		double *tmpy=(double*) malloc(sizeof(double)* Njack);
		printf("writing: %s\n",namefile);
		
		for (int i=Lrange[0] ; i<Lrange[1]; i++){
			double finalL=i;
			tmpx[0]=finalL;
			double *E3_m=(double*) malloc(sizeof(double) *Njack);
			for (int j=0;j<Njack;j++){
				E3_m[j]=fit_info_E3_poly.function(n,fit_info_E3_poly.Nvar , tmpx, fit_info_E3_poly.Npar, tif_E3_poly[j]);
			}
			double E3_m_err=error_jackboot(argv[1], Njack, E3_m );
			for (int j=0;j<Njack;j++){
				tmpx[1]=fit_info_m0.function(0,fit_info_m0.Nvar,tmpx,fit_info_m0.Npar,tif_m0[j]); //m0   put for each n the mass of the last ensemble
				
				tmpx[2]=gjack[0].jack[2][j];//m1  //
				tmpx[3]=gjack[0].jack[4][j];//E20
				tmpx[4]=gjack[0].jack[5][j];//E21
				tmpx[5]=(double) params[0].data.L[0];//T
				tmpx[6]=1;//k
				tmpx[7]=1;//MpiL
				
					
				tmpx[8]=E3_m[Njack-1];// E3( \vec{n} )/mass
				
				tmpx[9]=E3_m_err;
				
				
				for(int i=fit_info.Nvar ; i<fit_info.Nvar+ fit_info.n_ext_P; i++)
					tmpx[i]=fit_info.ext_P[i-fit_info.Nvar][j];
				
				
				
				tmpy[j]=fit_info.function(n,Nvar,tmpx,Npar,tif[j]);
							  
			}
			
			fprintf(f,"%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
			
			free(E3_m);
		}
		
		free(tmpy);free(tmpx);
		fclose(f); 
	}
	for (int i=0;i<Npar; i++)
		for (int j=0;j<Njack;j++)
			tif[j][i]=-1e+3;
	for (int n=0;n< N; n++){
		
		mysprintf(namefile,NAMESIZE,"%s/kiso_P-1e+3_n%d_L.txt",argv[3],n);
		f=open_file(namefile,"w+");
		double *tmpx=(double*) malloc(sizeof(double)* Nvar);
		double *tmpy=(double*) malloc(sizeof(double)* Njack);
		printf("writing: %s\n",namefile);
		
		for (int i=Lrange[0] ; i<Lrange[1]; i++){
			double finalL=i;
			tmpx[0]=finalL;
			double *E3_m=(double*) malloc(sizeof(double) *Njack);
			for (int j=0;j<Njack;j++){
				E3_m[j]=fit_info_E3_poly.function(n,fit_info_E3_poly.Nvar , tmpx, fit_info_E3_poly.Npar, tif_E3_poly[j]);
			}
			double E3_m_err=error_jackboot(argv[1], Njack, E3_m );
			for (int j=0;j<Njack;j++){
				tmpx[1]=fit_info_m0.function(0,fit_info_m0.Nvar,tmpx,fit_info_m0.Npar,tif_m0[j]); //m0   put for each n the mass of the last ensemble
				
				tmpx[2]=gjack[0].jack[2][j];//m1  //
				tmpx[3]=gjack[0].jack[4][j];//E20
				tmpx[4]=gjack[0].jack[5][j];//E21
				tmpx[5]=(double) params[0].data.L[0];//T
				tmpx[6]=1;//k
				tmpx[7]=1;//MpiL
				
				
				tmpx[8]=E3_m[Njack-1];// E3( \vec{n} )/mass
				
				tmpx[9]=E3_m_err;
				
				
				for(int i=fit_info.Nvar ; i<fit_info.Nvar+ fit_info.n_ext_P; i++)
					tmpx[i]=fit_info.ext_P[i-fit_info.Nvar][j];
				
				
				
				tmpy[j]=fit_info.function(n,Nvar,tmpx,Npar,tif[j]);
				
			}
			
			fprintf(f,"%g  \t %g  %g\n",finalL,tmpy[Njack-1], error_jackboot(argv[1],Njack, tmpy ) );
			
			free(E3_m);
		}
		
		free(tmpy);free(tmpx);
		fclose(f); 
	}
	
	free_2(Njack,tif);
	free_2(Njack,tif_m0);
	free_2(Njack,tif_E3_poly);
}

void print_phase_shift(char **argv,vector<data_phi> gjack , struct fit_type fit_info , const char* label, struct fit_result fit_out){
	
	char namefile[NAMESIZE];
	mysprintf(namefile,NAMESIZE,"%s/%s_fit_phase_shift.txt",argv[3], label);
	FILE *f=open_file(namefile,"w+");
	int Npar=fit_info.Npar;
	int Njack=gjack[0].Njack;
	
	int ikmax=100;
	double dk=0.03;
	double *delta=(double* ) malloc(sizeof(double)*Njack);
	for (int ik=0;ik<ikmax;ik++){
	
		double k_m=ik*dk;
		for (int j=0; j<Njack;j++){
		
			double a0m0=fit_out.P[0][j];
			double r0m0=fit_out.P[1][j];
			
			double kcotdelta_m=1.0/a0m0+ + r0m0*k_m*k_m/2.;  //   (k cot(d) )/ mass
			if (Npar>=3){
				kcotdelta_m+=fit_out.P[2][j]*r0m0*r0m0*r0m0*k_m*k_m*k_m*k_m;
			}
			delta[j]=std::atan( k_m/kcotdelta_m);
			
		}
		fprintf(f,"%g    %g   %g\n",k_m,delta[Njack-1], error_jackboot( argv[1] , Njack, delta) );
	}
	free(delta);
	fclose(f);
}


int main(int argc, char **argv){
	 error(argc!=4,1,"main ",
		 "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");   
	
	 //int Ne=0;   
	 /*cluster::IO_params *params=(cluster::IO_params*) malloc(sizeof(cluster::IO_params)*Ne);
	 
	 
	 data_phi data;
	 data_phi *dataj=(data_phi*) malloc(sizeof(data_phi*)*Ne);
	 
	 
	 char namefile[NAMESIZE];
	 mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L32_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
	 FILE *f=open_file(namefile,"r");
	 read_header_phi4( f  , params[0]);
	 read_dataj(f,params[0],dataj[0] );
	 fclose(f);
	 */
	
	 vector<cluster::IO_params> paramsj;
	 vector<data_phi> dataj;
	 
	 int Ne=0; 
	 cluster::IO_params params;
	 
	 char namefile[NAMESIZE];
	 char jackboot[NAMESIZE];
	 char resampling[NAMESIZE];
	 mysprintf(jackboot,NAMESIZE,"%s" ,argv[1] );
	 mysprintf(resampling,NAMESIZE,"%s" ,argv[1] );
	 
	 //mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T32_L16_msq0-4.925000_msq1-4.850000_l02.500000_l12.500000_mu5.000000_g0.000000_rep0",argv[2],argv[1]);
	 //emplace_back_par_data(namefile,paramsj,dataj);
	
	 //0
	 mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L15_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0",argv[2],argv[1]);
	 emplace_back_par_data(namefile,paramsj,dataj);
	 mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L16_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0",argv[2],argv[1]);
	 emplace_back_par_data(namefile,paramsj,dataj);
	 mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L17_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0",argv[2],argv[1]);
	 emplace_back_par_data(namefile,paramsj,dataj);
	 mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L18_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0",argv[2],argv[1]);
	 emplace_back_par_data(namefile,paramsj,dataj);
	 mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L19_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0",argv[2],argv[1]);
	 emplace_back_par_data(namefile,paramsj,dataj);
	 mysprintf(namefile,NAMESIZE,"%s/%s_G2t_T64_L20_msq0-1.267000_msq1-0.550000_l00.000000_l10.000000_mu0.000000_g5.000000_rep0",argv[2],argv[1]);
	 emplace_back_par_data(namefile,paramsj,dataj);
	 
	 
	 printf("E1_0 =%f   %f\n", dataj[0].jack[1][dataj[0].Njack-1], error_jackboot(argv[1],dataj[0].Njack, dataj[0].jack[1]  ) );
	
	 
	 vector<data_phi> gjack= create_generalised_resampling( dataj );
	 printf("GEVP_E2_01 =%f   %f\n", gjack[1].jack[19][gjack[1].Njack-1], error_jackboot(argv[1],gjack[1].Njack, gjack[1].jack[19]  ) );
	 Ne=gjack.size();
	 printf("number of ensembles = %d\n",Ne);
	 
	 std::vector<int> myen(Ne);
	 for(int i=0;i< Ne; i++)  myen[i]=i;
	
	 int Njack=gjack[0].Njack;
	 ///////////////////////////////////////////////////////////////////////////////////////////////////
	 // init zeta func
	 //////////////////////////////////////////////////////////////////////////////////////////////////
	 


	 std::vector<int> Ls;
	 for (int  i=10;i<42;i++){Ls.emplace_back(i);}

	 std::vector<double> masses;
	 std::vector<double> err_mass;
	 int e1;
	 for (int L: Ls){
		 e1=Ne-1;
		 for(int e=0; e< Ne;e++){
			 if (L==paramsj[e].data.L[1]){
				 e1=e;
				 printf("L=%d found   mass=%g\n",L,gjack[e1].jack[443][Njack-1]);
			 }
		 }
		 masses.emplace_back(gjack[e1].jack[443][Njack-1]);// 443 ->GEVP_0_3_1_001_011_meffl0
		 err_mass.emplace_back(  error_jackboot(argv[1],Njack, gjack[e1].jack[443] )  );
	 }
	 
	// zeta.Init_Lmq_g(Ls, masses, err_mass  );
	//zeta.write("zeta_complex_g5.dat");
	zeta.read("zeta_complex_g5.dat");
	///////////////////////////////////////////////////////////////////////////////////////////////////
	// printing files
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	printing_file_for_maxim_and_fernando_complex(argv, paramsj , gjack, myen);

	 ///////////////////////////////////////////////////////////////////////////////////////////////////
	 // start fitting
	 //////////////////////////////////////////////////////////////////////////////////////////////////
	 struct fit_type fit_info, fit_info_m0;
	 
	 struct fit_result  fit_m1, fit_m0;
	 fit_info.Nvar=14;                   fit_info_m0.Nvar=14;            
	 fit_info.Npar=2;                   fit_info_m0.Npar=2;           
	 fit_info.N=1;                      fit_info_m0.N=1;
	 fit_info.Njack=gjack[0].Njack;     fit_info_m0.Njack=gjack[0].Njack;
	 fit_info.n_ext_P=0;                fit_info_m0.n_ext_P=0;
	 fit_info.function=M_finite_volume; fit_info_m0.function=M_finite_volume;

	 
	 fit_m0=fit_data(argv,  paramsj ,gjack, M0_g_finite_volume_lhs ,fit_info, "M0_finite_vol" ,myen);
	
	 
	 fit_m1=fit_data(argv,  paramsj ,gjack, M1_finite_volume_lhs ,fit_info, "M1_finite_vol" ,myen);
	 
	 // corr 75 =E1_1_px
	 fit_result fit_m1_p1=fit_data(argv,  paramsj ,gjack, M1_p_finite_volume_lhs<75> ,fit_info, "M1_p1_finite_vol" ,myen);
	 free_fit_result(fit_info,fit_m1_p1);
	 printf("\n/////////////////////////////////     E2_0//////////////////\n");
	 ///////////////////////////////////////////////////////////////////////////////////////////////////
	 // E20
	 //////////////////////////////////////////////////////////////////////////////////////////////////
	
	 ///////////////////////////////////////////////////////////////////////////////////////////////////
	 // start fitting
	 
	 fit_info.Npar=1;
	 fit_info.N=1;
	 fit_info.Njack=gjack[0].Njack;// E1_0
	 fit_info.n_ext_P=0;
	 //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
	 fit_info.function=muDE_rhs;
	 
	 struct fit_result fit_a_00=fit_data(argv,  paramsj ,gjack, muDE_00_g_lhs ,fit_info, "a_00_luscher",myen );
	 
	 printf("\n/////////////////////////////////     k cot delta    //////////////////\n");
	 ///////////////////////////////////////////////////////////////////////////////////////////////////
	 // kcot
	 //////////////////////////////////////////////////////////////////////////////////////////////////
	 
	
	 fit_info.Npar=2;
	 fit_info.N=3;
	 fit_info.Njack=gjack[0].Njack;
	 fit_info.n_ext_P=0;
	 //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
	 fit_info.function=rhs_kcotd_m;
	 
	 struct fit_result fit_kcotd=fit_data(argv,  paramsj ,gjack, lhs_kcotd_m_g ,fit_info, "kcotd_m",myen );
	  
	///////////////////////////////////////////////////////////////////////////////////////////////////
	 // kcot  deltaE
	 //////////////////////////////////////////////////////////////////////////////////////////////////
		  printf("\n/////////////////////////////////     k cot delta deltaE    //////////////////\n");

	
	 fit_info.Npar=2;
	 fit_info.N=3;
	 fit_info.Njack=gjack[0].Njack;
	 fit_info.n_ext_P=0;
	 //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
	 fit_info.function=rhs_kcotd_m;
	 
	 struct fit_result fit_kcotd_DeltaE=fit_data(argv,  paramsj ,gjack, lhs_kcotd_m_deltaE_g ,fit_info, "kcotd_m_deltaE",myen );
	  
	

	 ///////////////////////////////////////////////////////////////////////////////////////////////////
	 printf("\n/////////////////////////////////   fit  k  form from_phase_shift    //////////////////\n");
	 //////////////////////////////////////////////////////////////////////////////////////////////////
	 
	 
	//  fit_info.Npar=2;
	//  fit_info.N=3;
	//  fit_info.Njack=gjack[0].Njack;
	//  fit_info.n_ext_P=0;
	//  //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
	//  fit_info.function=rhs_k_from_phase_shift_g;
	//  fit_info.lambda=0.001;
	//  fit_info.acc=0.01;
	//  fit_info.h=1e-3;
	//  fit_info.Prange={10,10};
	//  fit_info.devorder=2;
	//  fit_info.guess= {-0.121902,-4.9};
	 
	//  struct fit_result k_from_phase_shift=fit_data(argv,  paramsj ,gjack, lhs_k_g ,fit_info, "k_from_phase_shift_g",myen  );// {-0.948817,-114.788,0.0003987}
	//  print_fit_band_L_M( argv, gjack , fit_info,fit_info_m0 ,  "k_from_phase_shift_g",   k_from_phase_shift ,fit_m0,    paramsj,  myen, {10,20});
	//  fit_info.guess=std::vector<double>();
	 
	 
	 
	 ///////////////////////////////////////////////////////////////////////////////////////////////////
	 printf("\n/////////////////////////////////   fit  deltaE2_m_quant_cond  //////////////////\n");
	 //////////////////////////////////////////////////////////////////////////////////////////////////
	 fit_info.restore_default();
	 
	 fit_info.Npar=2;
	 fit_info.N=3;
	 fit_info.Njack=gjack[0].Njack;
	 fit_info.n_ext_P=0;
	 //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
	 fit_info.function=rhs_deltaE2_m_quant_cond_g;
	 
	 fit_info.lambda=0.001;
	 fit_info.acc=0.01;
	 fit_info.h=1e-3;
	 fit_info.Prange={100,10};
	 fit_info.devorder=2;
	 fit_info.guess= {-0.123, -4.67};
	 
	 //      the zeta is computed analytically, use the interpolated one for faster result!!!!!!!!!
	 //      struct fit_result k_from_phase_shift_3par=fit_data(argv,  paramsj ,gjack, lhs_k ,fit_info, "k_from_phase_shift_n5_3par",myen ,  {-0.11,-950, 6.4e-6} );// {-0.948817,-114.788,0.0003987}
	 fit_result deltaE2_m_quant_cond=fit_data(argv,  paramsj ,gjack, lhs_deltaE2_m_latt_g ,fit_info, "deltaE2_m_quant_cond",myen );
	 print_fit_band_L_M( argv, gjack , fit_info,fit_info_m0 ,  "deltaE2_m_quant_cond",   deltaE2_m_quant_cond ,fit_m0,    paramsj,  myen, {10,20});
	 
	 print_phase_shift(argv, gjack ,  fit_info , "deltaE2_m_quant_cond", deltaE2_m_quant_cond);
	 fit_info.restore_default();
	 
	 

	printf("\n/////////////////////////////////     delta 2par   //////////////////\n");
	 ///////////////////////////////////////////////////////////////////////////////////////////////////
	 // kcot
	 //////////////////////////////////////////////////////////////////////////////////////////////////
	 
	 
	 fit_info.Npar=2;
	 fit_info.N=1;
	 fit_info.Njack=gjack[0].Njack;
	 fit_info.n_ext_P=0;
	 //fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
	 fit_info.function=rhs_delta_g;
	 
	 struct fit_result fit_delta_2par=fit_data(argv,  paramsj ,gjack, lhs_delta_g ,fit_info, "delta_2par_g",myen );
	 free_fit_result(fit_info,fit_delta_2par);

	


	 ///////////////////////////////////////////////////////////////////////////////////////////////////
	 printf("\n/////////////////////////////////   fit  E3 quant cond  //////////////////\n");
	 //////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef PYTHON
	 //// we need python
	 wchar_t *program = Py_DecodeLocale(argv[0], NULL);
	 if (program == NULL) {
		 fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
		 exit(1);
	 }  
	 Py_SetProgramName(program);  /* optional but recommended */
	 
	 Py_Initialize();
	 ///////////// end python init
	 
	 
	 printf("//////////////////// poly fit E3   ////////////////////////////////////\n");
	 fit_type fit_info_E3_poly;
	 fit_info_E3_poly.N=2;
	 fit_info_E3_poly.Nvar=fit_info.Nvar;
	 
	 fit_info_E3_poly.Njack=gjack[0].Njack;
	 fit_info_E3_poly.n_ext_P=2;
	 fit_info_E3_poly.ext_P=(double**) malloc(sizeof(double*)*fit_info_E3_poly.n_ext_P);
	 
	 fit_info_E3_poly.ext_P[0]=deltaE2_m_quant_cond.P[0];
	 fit_info_E3_poly.ext_P[1]=deltaE2_m_quant_cond.P[1];
	 fit_info_E3_poly.Npar=fit_info_E3_poly.N*4;
	 fit_info_E3_poly.function=rhs_poly_order_E3_m<4>;
	 
	 fit_info_E3_poly.lambda=0.001;
	 fit_info_E3_poly.acc=0.01;
	 fit_info_E3_poly.h=1e-3;
	 fit_info_E3_poly.Prange={1000,10000};
	 fit_info_E3_poly.devorder=2;
	 
	 mysprintf(namefile,NAMESIZE,"poly_E3andE1_N%d",fit_info_E3_poly.N );
	// struct fit_result fit_QC3_poly=fit_data(argv,  paramsj ,gjack, lhs_E3_m ,fit_info_E3_poly, namefile,myen   );
	 
	 struct fit_result fit_QC3_poly=fit_data(argv,  paramsj ,gjack, lhs_E3orE1_g_m ,fit_info_E3_poly, namefile,myen   );
	 print_fit_band_L_M( argv, gjack , fit_info_E3_poly,fit_info_m0 ,  namefile,   fit_QC3_poly ,fit_m0,    paramsj,  myen, {23,41});
//      free_fit_result(fit_info,fit_QC3_poly);
//      fit_info_E3_poly.restore_default();
	 

printf("////////////////////  kiso const fit   ////////////////////////////////////\n");
	 init_python_detQC();
	 init_python_detQC_kcot_kiso("kcot_2par", "kiso_pole_fix", "find_2sol");

	 
	 fit_info.N=2;
	 fit_info.Njack=gjack[0].Njack;
	 
	 fit_info.n_ext_P=2;
	 fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
	 fit_info.ext_P[0]=deltaE2_m_quant_cond.P[0];
	 fit_info.ext_P[1]=deltaE2_m_quant_cond.P[1];
	 
	 
	 
	 fit_info.function=rhs_E3_m_QC3_2sol;
	 fit_info.Npar=2;
	 fit_info.guess={1,10};
	 fit_info.lambda=0.01;
	 fit_info.acc=0.01;
	 fit_info.h=1e-1;
	 fit_info.devorder=-2;
	 fit_info.verbosity=100;
	 fit_info.repeat_start=1;
	 
	 fit_info.mean_only=true;

	 mysprintf(namefile,NAMESIZE,"QC3_N%d_pole_fix",fit_info.N, fit_info.Npar);
	 struct fit_result fit_QC3_const=fit_data(argv,  paramsj ,gjack, lhs_E3orE1_g_m ,fit_info, namefile,
	 /*myen*/ {0,1,2}   );
	//  print_fit_band_E3_vs_L( argv, gjack , fit_info,fit_info_m0 ,  namefile,   fit_QC3_const ,fit_m0,    paramsj,  myen,  fit_info_E3_poly, fit_QC3_poly, {23,41});
	 
	 fit_info.restore_default();
exit(1);

	 printf("////////////////////  kiso pole fit   ////////////////////////////////////\n");
	 init_python_detQC();
	 init_python_detQC_kcot_kiso("kcot_2par", "kiso_pole");
//     init_python_detQC_kcot_kiso("kcot_2par", "kiso_2par");
//      init_python_detQC_kcot_kiso("kcot_2par", "kiso_1par");
	 fit_info.Npar=3;
	 fit_info.N=2;
	 fit_info.Njack=gjack[0].Njack;
	 
	 
	 
	 fit_info.function=rhs_E3_m_QC3_pole;
	 
	 fit_info.n_ext_P=2;
	 fit_info.ext_P=(double**) malloc(sizeof(double*)*fit_info.n_ext_P);
	 fit_info.ext_P[0]=deltaE2_m_quant_cond.P[0];
	 fit_info.ext_P[1]=deltaE2_m_quant_cond.P[1];
	 
	 
	 
	 fit_info.lambda=0.001;
	 fit_info.acc=0.01;
	 fit_info.h=1e-3;
	 //fit_info.Prange={1000,10000};
	 fit_info.devorder=2;
	 fit_info.verbosity=100;
	 fit_info.repeat_start=10;
	 fit_info.guess={-1e+6, 1,1e+6};
	 fit_info.mean_only=true;

	 mysprintf(namefile,NAMESIZE,"QC3_N%d_%dpar_pole",fit_info.N, fit_info.Npar);
	 struct fit_result fit_QC3_1par=fit_data(argv,  paramsj ,gjack, lhs_E3orE1_g_m ,fit_info, namefile,
	 /*myen*/ {0,1,2}   );
	 print_fit_band_E3_vs_L( argv, gjack , fit_info,fit_info_m0 ,  namefile,   fit_QC3_1par ,fit_m0,    paramsj,  myen,  fit_info_E3_poly, fit_QC3_poly, {23,41});
	 
	 fit_info.restore_default();
	
	
	 ///// close python
	 python_detQC_write_database();
	 python_detQC_free();
	 if (Py_FinalizeEx() < 0) {
		 exit(120);
	 }
	 PyMem_RawFree(program);
	 
#endif
	 return 0;
}
