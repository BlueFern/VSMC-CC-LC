#include <iostream>
#include <fstream>
//#include <cstdlib>
#include <math.h>
#include<sys/stat.h>
#include "include/macros.h"
#include "include/Model_constants.h"

/**
 * @file Initializing.cpp
 * \brief This file contains "void initialize(FILE *fs)" function
 */
#define P2(x) ((x)*(x))

extern double *y, *par, *y_converge;
extern double *L_cell, *L_rest;
extern double *G_gj, *IP3_P, *G_gj_varied, *IP3_P_varied;
extern int *num_par, *num_var, *num_cell;
extern double *Q_ip3r, *Q_serca, *Q_ryr;

/**************************************************************************************/
/**/void initialize(FILE *fs, FILE *fL, double tnow)/**/
/**************************************************************************************/
{
	fprintf(fL,"%f\t",tnow);
for (int i=0; i<num_cell[0]; i++)
{
fprintf(fs,"%d\t",i);

				Q_ip3r[i] 						= 	Q_ip3r_base;
				Q_serca[i]						=	Q_serca_base;
				Q_ryr[i]						=	Q_ryr_base;

				y[i * num_var[0] + smc_vm]		= 	vm_ini;
				y[i * num_var[0] + smc_ca]		=	Ca_ini;
				y[i * num_var[0] + smc_SR_ca]	=	SR_Ca_ini;

				par[i * num_par[0] + dl_bar]	= 	1/(1+exp(-y[i * num_var[0] + smc_vm]/8.3));
				par[i * num_par[0] + t_dl]		= 	((2.5*exp(-P2((y[i * num_var[0] + smc_vm]+40)/30)))+1.15)*1e-3; // s
				par[i * num_par[0] + fl_bar]	= 	1/(1+exp((y[i * num_var[0] + smc_vm]+42)/9.1));
				par[i * num_par[0] + t_fl]		= 	((65*exp(-P2((y[i * num_var[0] + smc_vm]+35)/25)))+45)*1e-3; // s
				y[i * num_var[0] + smc_dl]		=	1/(1+exp(-y[i * num_var[0] + smc_vm]/8.3));
				y[i * num_var[0] + smc_fl]		=	1/(1+exp((y[i * num_var[0] + smc_vm]+42)/9.1));

				par[i * num_par[0] + v_kca]		=	-41.7*log10(y[i * num_var[0] + smc_ca]*1e-6)-128.2;
				par[i * num_par[0] + po_bar]	=	1/(1+exp(-(y[i * num_var[0] + smc_vm]-par[i * num_par[0] + v_kca])/18.25));
				y[i * num_var[0] + smc_pf]		=	1/(1+exp(-(y[i * num_var[0] + smc_vm]-par[i * num_par[0] + v_kca])/18.25));
				y[i * num_var[0] + smc_ps]		=	y[i * num_var[0] + smc_pf];


				G_gj[i]= G_gj_varied[0];
				IP3_P[i] = IP3_P_varied[0];

				L_cell[i] = nonstim_NE;
				fprintf(fL,"%f\t",L_cell[i]);

				par[i * num_par[0] + rho_rg] 	=	L_cell[i]/(1.0*(k_1+L_cell[i]));
				y[i * num_var[0] + smc_G]		=	K_aG*(delta_G + par[i * num_par[0] + rho_rg])*G_TG/(K_aG*(delta_G + par[i * num_par[0] + rho_rg]) + K_dG);
				par[i * num_par[0] + r_hg]		=	eta_G*y[i * num_var[0] + smc_G]*y[i * num_var[0] + smc_ca]/(y[i * num_var[0] + smc_ca]+K_cG);
				y[i * num_var[0] + smc_ip3]		=	(par[i * num_par[0] + r_hg]*PIP_2T/gamma_G)/K_degG;
				y[i * num_var[0] + smc_DAG]		=	(par[i * num_par[0] + r_hg]*PIP_2T/gamma_G)/K_DAG;

				y[i * num_var[0] + smc_R01]		=	0.001004;
				y[i * num_var[0] + smc_R10]		=	0.002128;
				y[i * num_var[0] + smc_R00]		=	0.996866;
				y[i * num_var[0] + smc_R11]		=	1-y[i * num_var[0] + smc_R01]-y[i * num_var[0] + smc_R00]-y[i * num_var[0] + smc_R10];

				y[i * num_var[0] + smc_x00]		=	0.789487;
				y[i * num_var[0] + smc_x10]		=	0.002412;
				y[i * num_var[0] + smc_x01]		=	0.206240;
				y[i * num_var[0] + smc_x11]		=	1-y[i * num_var[0] + smc_x00]-y[i * num_var[0] + smc_x10]-y[i * num_var[0] + smc_x01];

				par[i * num_par[0] + E_Ca]			=	(R*T/(z_Ca*F))* log(Ca_e/y[i * num_var[0] + smc_ca]);
				par[i * num_par[0] + E_K]			=	(R*T/(F))* log(K_e/K_cyt);
				par[i * num_par[0] + E_Cl]			=	(R*T/(F*z_Cl))* log(Cl_e/Cl_cyt);
				par[i * num_par[0] + E_Na]			=	(R*T/(F))* log(Na_e/Na_cyt);

				for (int j=0;j<num_var[0];j++)
							y_converge[i * num_var[0] + j] = y[i * num_var[0] + j];
}
fprintf(fL,"\n");
fclose(fs);
}

