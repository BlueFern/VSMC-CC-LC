/**
 * @file Model_constants.h
 * @brief Model Constants:
 * Contains Standard VSMC Parameters,  Initial values of state variables,
 * Parameters for Membrane Channels, Pumps and exchangers, Parameters for Sarcoplasmic Reticulum channels and pumps
 * Parameters for Calcium buffering in cytosol, Parameters for alpha-adrenoceptor cascade, IP3 and DAG formation and
 * Boundary values of agonist concentration in the sigmoidal disribution */

/* Standard VSMC Parameters */

const double	C_m		= 	20;  ///< Cell Capacitance, pC V-1

const double	z_K		= 	1.0; ///< Potassium Valency
const double	z_Na	= 	1.0; ///< Sodium Valency
const double	z_Ca	= 	2.0; ///< Calcium Valency
const double	z_Cl	= 	-1.0; ///< Chloride Valency

const double    N_AV	= 	6.022e23; ///< Avagadro's Number
const double    F 		= 	96487.0; ///< Faraday's Constant, C mol-1
const double   	R 		= 	8.341e3; ///< Universal Gas constant, mJ mol-1 K-1
const double   	T 		= 	293.0; ///< Temperature, K

const double   	Vol_smc	= 	1;  ///< Cell Volume, pL
const double   	Vol_cyt	= 	0.7;  ///< Cytosol Volume, pL
const double   	Vol_SR	= 	Vol_smc*0.05;  ///< SR Volume, pL

const double   	beta			= z_Ca*F*Vol_SR*1e-9;  ///< ion concentration to ion current coonversion factor for intraSR, pC nM-1
const double   	alpha		= z_Ca*F*Vol_cyt*1e-9;  ///< ion concentration to ion current coonversion factor for cytosol, pC nM-1

const double  	K_cyt 	= 	140e6;  ///< Potassium Concentration in cytosol, nM
const double  	Na_cyt 	= 	8.4e6;  ///< Sodium Concentration in cytosol, nM
const double  	Cl_cyt 	= 	59.4e6;  ///< Chloride Concentration in cytosol, nM

const double  	Ca_e 	=  	2e6;///< nM
const double 	K_e 	= 	5e6;///< Potassium Concentration in extracellular space, nM
const double 	Na_e 	= 	140e6;///< Sodium Concentration in extracellular space, nM
const double 	Cl_e 	= 	129e6;///< Chloride Concentration in extracellular space, nM

const double 	f2_rt	=	F*F/(R*T);
const double 	fca2_rt	=	z_Ca*z_Ca*F*F/(R*T);
const double 	f_rt	=	F/(R*T);

/* Initial values of state variables */
const double   	SR_Ca_ini 	=	1.2149e6; ///< Initial Concentration of intraSR calcium, nM
const double  	vm_ini 	= 	-54.03;  ///< Initial value of Membrane Potential, mV
const double   	Ca_ini 	= 	80.55;  ///< Initial concentration of intracellular calcium, nM


/* Parameters for Membrane Channels, Pumps and exchangers */
const double  	G_vocc			= 	26.3e-1; ///< Whole cell conductance of VOCC, nS

const double  	G_bkca			= 	111e-1;///< Whole cell conductance of BKCa channel, nS
const double 	t_pf	=	0.84e-3; ///< Fast activation time constant, s
const double 	t_ps	=	35.9e-3; ///< Slow activation time constant, s

const double   	G_cacc	= 	0.112 * C_m; ///< Whole cell conductance of CACC, nS
const double   	n_cacc	= 	2; ///< order of the function shows depenedncy on intracellular calcium
const double   	k_cacc	= 	587; ///< Dissociation constant of CACC, nM

const double    G_ncx	= 	0.25e-3;///< Whole cell conductance of NCX channel, pS
const double   	d_ncx	= 	0.0003;
const double  	gamma_ncx	= 	0.45; ///< Position of the energy barrier of NCX channel
const double  	k_ncx	= 	125; ///< Dissociation constant for allosteric component, nM

const double   	E_nsc	= 	9; ///< Reverse potential of NSC channel, mV
const double   	k_nsc	= 	1e3;  ///< Dissociation constant of NSC, nM
const double   	G_nsc	= 	0.61; ///< Whole cell conductance of NSC, nS

const double  	k_pmca	= 	150; ///< Dissociation constant of PMCA, nM
const double  	Q_pmca	= 	9.993e3;///< Maximum rate of calcium flux of PMCA, nM s-1
const double  	n_pmca	= 	1.0; ///<  order of function

/* Parameters for Sarcoplasmic Reticulum channels and pumps */
/* Binding Constants for IP3R channel */
const double   	a1_ip3	= 	167.6e-3; ///< IP3, nM-1 s-1
const double   	a2_ip3	= 	3.81e-3; ///< Calcium inhibition, nM-1 s-1
const double   	a3_ip3	= 	413.4e-3; ///< IP3, nM-1 s-1
const double   	a4_ip3	= 	0.3101e-3; ///< Calcium inhibition, nM-1 s-1
const double   	a5_ip3	= 	53.9e-3; ///< Calcium activation, nM-1 s-1v
/* Dissociation Constants for IP3R channel */
const double   	b1_ip3	= 	228; ///< IP3, s-1
const double   	b2_ip3	= 	0.409; ///< Calcium inhibition, s-1
const double   	b3_ip3	= 	188.5; ///< IP3, s-1
const double   	b4_ip3	= 	0.096; ///< Calcium inhibition, s-1
const double   	b5_ip3	= 	4.52; ///< Calcium activation, s-1

const double   	Q_leak	= 	12;
const double  	k_leak	= 	145e3; ///< Exponential growth factor, nM

const double   	k_serca	= 	500;  ///< Dissociation constant for SERCA pump, nM
const double   	n_serca	= 	2; ///< order of function

const double   	Q_ip3r_base		= 2000; ///< Base value of Q_ip3r, s-1
const double   	Q_ryr_base		= 4000; ///< Base value of Q_ryr, s-1
const double   	Q_serca_base		= 210e3; // Base value of Q_serca, nM/s

/* Binding Constants for RyR channel */
const double   	SR_kr1			= 2.5e-6;  ///< Activation rate constant, nM-2 s-1
const double   	SR_kmr1			= 7.6;  ///< Unbinding rate constant from activation, s-1
const double   	SR_kr2			= 1.05e-3; ///< Inactivation rate constant, nM-1 s-1
const double   	SR_kmr2			= 84; ///< Unbinding rate constant from inactivation, s-1


/* Parameters for Calcium buffering in cytosol */
const double   	scm_bar 	= 	0.1e6;  ///< Concentration of Calmodulin, nM
const double   	bf_bar 	= 	0.1e6;  ///< Concentration of other buffers, nM
const double   	k_d 		= 	260;  ///< Dissociation constant of Calmodulin, nM
const double   	k_db 	= 	530;  ///< Dissociation constant of other buffers, nM

/* Parameters for alpha-adrenoceptor cascade and IP3 formation */
const double   	k_1		= 	0.06e6;  ///< nM 0.01e6
const double   	G_TG	= 	1e5;  ///< Total number of G-protein molecules
const double   	K_degG	= 	1.25;  ///< IP3 degradation rate, s-1
const double   	K_aG	= 	0.17;  ///< G-protein activation rate, s-1
const double   	K_dG	= 	1.5;  ///< G-protein inactivation rate, s-1
const double   	PIP_2T	= 	5e7;  ///< total number of PIP2 molecules
const double   	K_cG	= 	0.4e3;  ///< Dissociation constant for calcium binding to PLC, nM
const double   	eta_G	= 	2.781e-5;  ///< Effective signal gain parameter, s-1
const double   	gamma_G	= 	N_AV*Vol_smc*1e-21;
const double   	K_pG		= 	0.01;  ///< s-1
const double   	K_DAG	= 	1.0;  ///< DAG degradation rate, s-1
const double	delta_G	=	1.235e-3; ///< G-protein intrinsic activity parameter

/* Boundary values of agonist concentration in the sigmoidal disribution */
const double   stim_NE 			= 400;
const double   nonstim_NE 		= 100;
