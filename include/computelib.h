#include <iostream>
#include <fstream>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros.h"
#include "Model_constants.h"

#define P2(x) ((x)*(x))
#define P3(x) ((x)*(x)*(x))
#define P4(x) ((x)*(x)*(x)*(x))

double *par; ///< par - model parameters in macros file
double *y; ///< y - model state variables in macros file
double *y_temp; ///< y_temp - temporary storage of y
double *y_converge; ///< y_converge - temporary storage of converged results

int *num_var; ///< num_var - total number of state variables
int *num_par; ///< num_par - total number of model parameters in macros file
int *num_cell; ///< num_cell - total number of cells coupled

double *TDMB; ///< tridiagonal matrix- main diagonal elements for solving V_m using TDM algorithm
double *TDMD; ///< tridiagonal matrix - constant vector for solving V_m using TDM algorithm
double *TDMCC; ///< tridiagonal matrix - elements for solving V_m using TDM algorithm
double *TDMDD; ///< tridiagonal matrix - elements for solving V_m using TDM algorithm
double *TDMB_IP3; ///< tridiagonal matrix- main diagonal elements for solving ip3 using TDM algorithm
double *TDMD_IP3; ///< tridiagonal matrix - constant vector for solving ip3 using TDM algorithm
double *TDMCC_IP3; ///< tridiagonal matrix elements for solving ip3 using TDM algorithm
double *TDMDD_IP3; ///< tridiagonal matrix elements for solving ip3 using TDM algorithm

double *G_gj; ///< Gap junctional conductance
double *G_gj_base; ///< Gap junctional conductance - base value
double *G_gj_varied; ///< Gap junctional conductance - new value in the loop
double *IP3_P; ///< IP3 coupling coefficient
double *IP3_P_base; ///< IP3 coupling coefficient - base value
double *IP3_P_varied; ///< IP3 coupling coefficient - new value in the loop
double *L_cell; ///< Agonist concentration in the cells
double *L_rest; ///< Agonist cocentration before the stimulation starts

double *Q_ip3r; ///< Maximum rate of ip3r channel
double *Q_serca; ///< Maximum rate of serca pump
double *Q_ryr; ///< Maximum rate of ryr channel

void initialize(FILE *fs); ///< This function initializes the variables
void allocatememory(int var_tot, int par_tot); ///< This function allocates memory for all the variables defined above
void singlecell(double tnow, double interval); ///< This function solves single cell dynamics
void multicell(double tnow, double interval); ///< This function solves coupled cells dynamics
void dump_data(FILE *ft, FILE *fvm, FILE *fca, FILE *fip3, FILE *fncx, FILE *fvocc,  FILE *fbkca, FILE *fip3r, FILE *fsrca, FILE *fryr, FILE *fdag, double tnow); ///< This function writes output as .txt files
