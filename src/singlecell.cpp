#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include<sys/stat.h>
#include "include/macros.h"
#include "include/Model_constants.h"

/**
 * @file singlecell.cpp
 * \brief This file contains "void singlecell(double tnow, double interval)" function
 */

extern double *y, *par, *y_temp, *y_converge;
extern double *L_cell;
extern int *num_var, *num_par, *num_cell;
extern double *Q_ip3r, *Q_serca, *Q_ryr;

#define P2(x) ((x)*(x))
#define P3(x) ((x)*(x)*(x))
#define P4(x) ((x)*(x)*(x)*(x))

/**************************************************************************************/
/**/void singlecell(double tnow, double interval) /**/
/**************************************************************************************/
{

	for (int i=0; i<num_cell[0]; i++)
	{

		/* Adrenoceptor cascade */
				par[i * num_par[0] + rho_rg]		=	L_cell[i] / (k_1 + L_cell[i]);
				par[i * num_par[0] + r_hg]			=	eta_G * y_converge[i * num_var[0] + smc_G] / (1 + K_cG / y_converge[i * num_var[0] + smc_ca]);

		/* VOCC channel */
				par[i * num_par[0] + E_Ca]			=	(R * T / (z_Ca * F)) * log(Ca_e / y_converge[i * num_var[0] + smc_ca]);
				par[i * num_par[0] + dl_bar]		= 	1 / (1 + exp(-y_converge[i * num_var[0] + smc_vm] / 8.3));
				par[i * num_par[0] + t_dl]			= 	((2.5 * exp(-P2((y_converge[i * num_var[0] + smc_vm] + 40) / 30))) + 1.15) * 1e-3; // s
				par[i * num_par[0] + fl_bar]		= 	1 /(1 + exp((y_converge[i * num_var[0] + smc_vm] + 42) / 9.1));
				par[i * num_par[0] + t_fl]			= 	((65 * exp(-P2((y_converge[i * num_var[0] + smc_vm] + 35) / 25))) + 45) * 1e-3; // s
				par[i * num_par[0] + I_vocc]		=	G_vocc * y_converge[i * num_var[0] + smc_dl] *  y_converge[i * num_var[0] + smc_fl] * ( y_converge[i * num_var[0] + smc_vm] - par[i * num_par[0] + E_Ca]);


		/* BKCa channel */
				par[i * num_par[0] + v_kca]			=	-41.7 * log10(y_converge[i * num_var[0] + smc_ca] * 1e-6) - 128.2;
				par[i * num_par[0] + p_kca]			=	0.17 * y_converge[i * num_var[0] + smc_pf] + 0.83 * y_converge[i * num_var[0] + smc_ps];
				par[i * num_par[0] + po_bar]		=	1 / (1 + exp(-(y_converge[i * num_var[0] + smc_vm]-par[i * num_par[0] + v_kca]) / 18.25));
				par[i * num_par[0] + I_bkca]		= 	G_bkca * par[i * num_par[0] + p_kca] * (y_converge[i * num_var[0] + smc_vm] - par[i * num_par[0] + E_K]);

		/* CACC channel */
				par[i * num_par[0] + p_cacc]		=	(1 / (1 + pow(k_cacc / y_converge[i * num_var[0] + smc_ca], n_cacc)));
				par[i * num_par[0] + I_cacc]		=	G_cacc * par[i * num_par[0] + p_cacc] * (y_converge[i * num_var[0] + smc_vm] - par[i * num_par[0] + E_Cl]);

		/* NSC channel */
				par[i * num_par[0] + p_nscvm]		=	1 / (1+exp(-(y_converge[i * num_var[0] + smc_vm] + 78) / 45));
				par[i * num_par[0] + p_nscdag]		=	1 / (1 + k_nsc / y_converge[i * num_var[0] + smc_DAG]);
				par[i * num_par[0] + I_nsc]			=	G_nsc * (par[i * num_par[0] + p_nscdag]) * par[i * num_par[0] + p_nscvm] * (y_converge[i * num_var[0] + smc_vm] - E_nsc);

		/* NCX channel */
				par[i * num_par[0] + phiF]			=	exp(gamma_ncx * y_converge[i * num_var[0] + smc_vm] * f_rt);
				par[i * num_par[0] + phiR]			=	exp((gamma_ncx - 1) * y_converge[i * num_var[0] + smc_vm] * f_rt);
				par[i * num_par[0] + p_ncxallo]		=   1 / (1 + P2(k_ncx / y_converge[i * num_var[0] + smc_ca]));
				par[i * num_par[0] + p_ncxelec]		=	((P3(Na_cyt) * Ca_e * par[i * num_par[0] + phiF]) - (P3(Na_e) * y_converge[i * num_var[0] + smc_ca] * par[i * num_par[0] + phiR])) / (1 + d_ncx * (P3(Na_cyt) * Ca_e + P3(Na_e) * y_converge[i * num_var[0] + smc_ca]));
				par[i * num_par[0] + I_ncx]			=	G_ncx * par[i * num_par[0] + p_ncxallo] * par[i * num_par[0] + p_ncxelec];

		/* PMCA pump */
				par[i * num_par[0] + p_pmca]		=	1 / (1 + pow(k_pmca / y_converge[i * num_var[0] + smc_ca], n_pmca));
				par[i * num_par[0] + I_pmca] 		=	Q_pmca * alpha * par[i * num_par[0] + p_pmca];

		/* IP3R channel */
				par[i * num_par[0] + k1_ip3]		=	b1_ip3 / (a1_ip3 * y_converge[i * num_var[0] + smc_ip3]);
				par[i * num_par[0] + k2_ip3]		=	b2_ip3 / (a2_ip3 * y_converge[i * num_var[0] + smc_ip3]);
				par[i * num_par[0] + k3_ip3]		=	b3_ip3 / (a3_ip3 * y_converge[i * num_var[0] + smc_ip3]);
				par[i * num_par[0] + k4_ip3]		=	b4_ip3 / (a4_ip3 * y_converge[i * num_var[0] + smc_ip3]);
				par[i * num_par[0] + k5_ip3]		=	b5_ip3 / (a5_ip3 * y_converge[i * num_var[0] + smc_ip3]);
				par[i * num_par[0] + p_ip3r]		=	P4(y_converge[i * num_var[0] + smc_x10]) + 4 * P3(y_converge[i * num_var[0] + smc_x10]) * (1 - y_converge[i * num_var[0] + smc_x10]);
				par[i * num_par[0] + I_ip3r]		=	beta * Q_ip3r[i]* par[i * num_par[0] + p_ip3r]  * (y_converge[i * num_var[0] + smc_SR_ca] - y_converge[i * num_var[0] + smc_ca]);

		/* SERCA pump */
				par[i * num_par[0] + p_serca]		=	1 /(1 + pow(k_serca / y_converge[i * num_var[0] + smc_ca], n_serca));
				par[i * num_par[0] + I_serca]		=	Q_serca[i] * alpha * par[i * num_par[0] + p_serca];

		/* SRleak channel */
				par[i * num_par[0] + p_srleak]		=	exp(y_converge[i * num_var[0] + smc_SR_ca] / k_leak);
				par[i * num_par[0] + I_srleak]		=	beta * Q_leak * par[i * num_par[0] + p_srleak];

		/* RyR channel */
				par[i * num_par[0] + I_ryr]			=	beta * Q_ryr[i] * P2(y_converge[i * num_var[0] + smc_R10]) * (y_converge[i * num_var[0] + smc_SR_ca] - y_converge[i * num_var[0] + smc_ca]);

		/* Ca buffering */
				par[i * num_par[0] + rhocell]		=	(1 + (scm_bar*k_d / P2(k_d + y_converge[i * num_var[0] + smc_ca])) + (bf_bar * k_db / P2(k_db + y_converge[i * num_var[0] + smc_ca])));

/* State variables */

		/* Adrenoceptor cascade */
				y_temp[i * num_var[0] + smc_G]		=	((y[i * num_var[0] + smc_G] / interval) + (K_aG * (delta_G+par[i * num_par[0] + rho_rg]) * G_TG)) / (K_dG + K_aG * (delta_G + par[i * num_par[0] + rho_rg]) + 1 / interval);

		/* VOCC channel */
				y_temp[i * num_var[0] + smc_dl]		=	((y[i * num_var[0] + smc_dl] / interval) + (par[i * num_par[0] + dl_bar]) / par[i * num_par[0] + t_dl]) / (1 / interval + 1 / par[i * num_par[0] + t_dl]);
				y_temp[i * num_var[0] + smc_fl]		=	((y[i * num_var[0] + smc_fl] / interval) + (par[i * num_par[0] + fl_bar]) / par[i * num_par[0] + t_fl]) / (1 / interval + 1 / par[i * num_par[0] + t_fl]);

		/* BKCa channel */
				y_temp[i * num_var[0] + smc_pf]		=	((y[i * num_var[0] + smc_pf] / interval) + (par[i * num_par[0] + po_bar]) / t_pf) / (1 / interval + 1 / t_pf);
				y_temp[i * num_var[0] + smc_ps]		=	((y[i * num_var[0] + smc_ps] / interval) + (par[i * num_par[0] + po_bar]) / t_ps) / (1 / interval + 1 / t_ps);

		/* RyR channel */
				y_temp[i * num_var[0] + smc_R10]	=	((y[i * num_var[0] + smc_R10] / interval) + SR_kr1 * P2(y_converge[i * num_var[0] + smc_ca]) * y_converge[i * num_var[0] + smc_R00] + SR_kmr2 * y_converge[i * num_var[0] + smc_R11]) / (SR_kmr1 + SR_kr2 * y_converge[i * num_var[0] + smc_ca] + 1 / interval);
				y_temp[i * num_var[0] + smc_R11]	=	((y[i * num_var[0] + smc_R11] / interval) + SR_kr2 * (y_converge[i * num_var[0] + smc_ca]) * y_converge[i * num_var[0] + smc_R10] + SR_kr1 * P2(y_converge[i * num_var[0] + smc_ca]) * y_converge[i * num_var[0] + smc_R01]) / ((SR_kmr1 + SR_kmr2) + 1 / interval);
				y_temp[i * num_var[0] + smc_R01]	=	((y[i * num_var[0] + smc_R01] / interval) + SR_kr2 * (y_converge[i * num_var[0] + smc_ca]) * y_converge[i * num_var[0] + smc_R00] + SR_kmr1 * y_converge[i * num_var[0] + smc_R11]) / ((SR_kmr2 + SR_kr1 * P2(y_converge[i * num_var[0] + smc_ca])) + 1 / interval);
				y_temp[i * num_var[0] + smc_R00]	=	1 - y_temp[i * num_var[0] + smc_R01] - y_temp[i * num_var[0] + smc_R11] - y_temp[i * num_var[0] + smc_R10];

		/* IP3R channel */
				y_temp[i * num_var[0] + smc_x00] 	= 	((y[i * num_var[0] + smc_x00]/interval) + (b4_ip3*par[i * num_par[0] + k3_ip3]+b2_ip3)*(y_converge[i * num_var[0] + smc_x01]/(1+par[i * num_par[0] + k3_ip3]))+b5_ip3*y_converge[i * num_var[0] + smc_x10])/(((a4_ip3*par[i * num_par[0] + k1_ip3]+a5_ip3+a2_ip3)*y_converge[i * num_var[0] + smc_ca]/(1+par[i * num_par[0] + k1_ip3])) + 1/interval);
				y_temp[i * num_var[0] + smc_x10]	= 	((y[i * num_var[0] + smc_x10]/interval) + b2_ip3*y_converge[i * num_var[0] + smc_x11]+(a5_ip3*y_converge[i * num_var[0] + smc_ca])*(y_converge[i * num_var[0] + smc_x00]/(1+par[i * num_par[0] + k1_ip3])))/((a2_ip3*y_converge[i * num_var[0] + smc_ca]+b5_ip3) + 1/interval);
				y_temp[i * num_var[0] + smc_x01]	=	((y[i * num_var[0] + smc_x01]/interval) + (a4_ip3*par[i * num_par[0] + k1_ip3]+a2_ip3)*y_converge[i * num_var[0] + smc_ca]*(y_converge[i * num_var[0] + smc_x00]/(1+par[i * num_par[0] + k1_ip3]))+b5_ip3*y_converge[i * num_var[0] + smc_x11])/(((b4_ip3*par[i * num_par[0] + k3_ip3]+b2_ip3+a5_ip3*y_converge[i * num_var[0] + smc_ca])/(1+par[i * num_par[0] + k3_ip3])) + 1/interval);
				y_temp[i * num_var[0] + smc_x11]	=	1-y_temp[i * num_var[0] + smc_x00]-y_temp[i * num_var[0] + smc_x10]-y_temp[i * num_var[0] + smc_x01];

		/* IP3 and DAG concentration */
				y_temp[i * num_var[0] + smc_ip3]	=	((y[i * num_var[0] + smc_ip3] / interval) + (par[i * num_par[0] + r_hg] * PIP_2T / gamma_G)) / (K_degG + 1 / interval);
				y_temp[i * num_var[0] + smc_DAG]	=	((y[i * num_var[0] + smc_DAG] / interval) + (par[i * num_par[0] + r_hg] * PIP_2T / gamma_G)) / (K_DAG + 1 / interval);

		/* IntraSR calcium */
				y_temp[i * num_var[0] + smc_SR_ca]	=	((y[i * num_var[0] + smc_SR_ca] / interval) + par[i * num_par[0] + I_serca] / beta - par[i * num_par[0] + I_srleak] / beta	+ Q_ip3r[i] * par[i * num_par[0] + p_ip3r] * y_converge[i * num_var[0] + smc_ca] +  Q_ryr[i] * P2(y_converge[i * num_var[0] + smc_R10]) * y_converge[i * num_var[0] + smc_ca]) / ( Q_ip3r[i] * par[i * num_par[0] + p_ip3r] +  Q_ryr[i] * P2(y_converge[i * num_var[0] + smc_R10]) + 1 / interval);

		/* Intracellular calcium */
				y_temp[i * num_var[0] + smc_ca]		=	((y[i * num_var[0] + smc_ca] * alpha * par[i * num_par[0] + rhocell] / (interval)) - par[i * num_par[0] + I_serca] + par[i * num_par[0] + I_srleak] + Q_ip3r[i] * beta * par[i * num_par[0] + p_ip3r] * y_converge[i * num_var[0] + smc_SR_ca] +  beta * Q_ryr[i] * P2(y_converge[i * num_var[0] + smc_R10]) * y_converge[i * num_var[0] + smc_SR_ca] +  2 * par[i * num_par[0] + I_ncx] -par[i * num_par[0] + I_vocc] - par[i * num_par[0] + I_pmca]) / (Q_ip3r[i] * beta * par[i * num_par[0] + p_ip3r] + beta * Q_ryr[i] * P2(y_converge[i * num_var[0] + smc_R10]) + (alpha * par[i * num_par[0] + rhocell] / interval));

		/* Membrane potential */
				y_temp[i * num_var[0] + smc_vm]		=	((y[i * num_var[0] + smc_vm] * C_m*1e-3 / interval) + par[i * num_par[0] + E_K] * G_bkca * par[i * num_par[0] + p_kca] + G_vocc * y_converge[i * num_var[0] + smc_dl] * y_converge[i * num_var[0] + smc_fl] * par[i * num_par[0] + E_Ca] + G_cacc * par[i * num_par[0] + p_cacc] * par[i * num_par[0] + E_Cl] - par[i * num_par[0] + I_ncx] + G_nsc * (par[i * num_par[0] + p_nscdag]) * par[i * num_par[0] + p_nscvm] * E_nsc - par[i * num_par[0] + I_pmca]) / ((C_m*1e-3/interval)  +  G_bkca * par[i * num_par[0] + p_kca] + G_vocc * y_converge[i * num_var[0] + smc_dl] *  y_converge[i * num_var[0] + smc_fl] + G_cacc * par[i * num_par[0] + p_cacc] + G_nsc * (par[i * num_par[0] + p_nscdag]) * par[i * num_par[0] + p_nscvm]);

	}


}
