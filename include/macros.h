/* VSMC model macros */

/**
 * @file macros.h
 * \brief This file contains macros required for the state variables and model parameters
 */

/* State variables */

#define    		smc_vm			0 ///< Membrane potential
#define    		smc_ca			1 ///< Intracellular calcium
#define			smc_ip3			2 ///< IP3 concentration
#define			smc_SR_ca		3 ///< IntraSR concentration
#define			smc_DAG			4 ///< DAG concentration
#define			smc_x00			5 ///< X00 state in IP3R model
#define			smc_x01			6 ///< X01 state in IP3R model
#define			smc_x10			7 ///< X10 state in IP3R model
#define			smc_x11			8 ///< X11 state in IP3R model
#define			smc_dl			9 ///< VOCC activation
#define			smc_fl			10 ///< VOCC inactivation
#define			smc_pf			11 ///< BKCa fast activation
#define			smc_ps			12 ///< BKCa slow activation
#define			smc_R10			13 ///< R10 state in RyR model
#define			smc_R11			14 ///< R11 state in RyR model
#define			smc_R01			15 ///< R01 state in RyR model
#define			smc_R00			16 ///< R00 state in RyR model
#define			smc_G			17 ///< G-protein concentration

/* Model parameters */

#define			rhocell		0 ///< Buffering parameter
#define			E_K			1 ///< Nernst potential of potassium
#define			E_Ca		2 ///< Nernst potential of Calcium
#define 		E_Cl		3 ///< Nernst potential of Chloride
#define 		E_Na		4 ///< Nernst potential of Sodium
#define			I_bkca		5 ///< BKCa channel current
#define			I_vocc		6 ///< VOCC channel current
#define 		I_ip3r		7 ///< IP3R channel current
#define 		I_serca		8 ///< SERCA pump current
#define 		I_srleak	9 ///< SRleak current
#define 		rho_rg		10 ///< Ratio of agonist binding
#define 		r_hg		11 ///< Rate coefficient of PIP2 hydrolysis
#define 		k1_ip3		12 ///< Ratio of binding and dissociation constants for IP3 with IP3 concentration
#define 		k2_ip3		13 ///< Ratio of binding and dissociation constants for Ca inhibition with IP3 concentration
#define 		k3_ip3		14 ///< Ratio of binding and dissociation constants for IP3 with IP3 concentration
#define 		k4_ip3		15 ///< Ratio of binding and dissociation constants for Ca inhibition with IP3 concentration
#define 		k5_ip3		16 ///< Ratio of binding and dissociation constants for Ca activation with IP3 concentration
#define 		I_cacc		17 ///< CACC channel current
#define 		p_cacc		18 ///< CACC channel open probability
#define 		v_kca		19 ///< V0.5BKCa
#define 		dl_bar		20 ///< Steady state activation of VOCC
#define 		t_dl		21 ///< Activation time constant
#define 		fl_bar		22 ///< Steady state inactivation of VOCC
#define 		t_fl		23 ///< Inactivation time constant
#define 		po_bar		24 ///< Steady state open probability of BKCa channel
#define 		p_kca		25 ///< Open probability of BKCa channel
#define 		I_ncx		26 ///< NCX exchanger current
#define 		p_ncxallo	27 ///< NCX allosteric factor
#define 		p_ncxelec	28 ///< NCX electrochemical part
#define 		phiF		29 ///< NCX dependency on membrane potential
#define 		phiR		30 ///< NCX dependency on membrane potential
#define 		I_pmca		31 ///< PMCA current
#define 		p_pmca		32 ///< PMCA open probability
#define			p_nscvm		33 ///< NSC channel open probability dependency to membrane potential
#define			p_nscdag	34 ///< NSC channel open probability dependency to DAG concentration
#define			I_nsc		35 ///< NSC channel current
#define			I_ryr		36 ///< RyR channel current
#define 		p_serca		37 ///< SERCA pump open probability
#define 		p_srleak	38 ///< SRleak channel open probability
#define 		p_ip3r		39 ///< IP3R channel open probability
#define			I_gj		40 ///< Net gap junctional current
