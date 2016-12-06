#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include<sys/stat.h>
#include "include/macros.h"
#include "include/Model_constants.h"

/**
 * @file Outputdata.cpp
 * \brief This file contains "void dump_data(FILE *ft, FILE *fvm, FILE *fca,FILE *fpmca, FILE *fsrleak, FILE *fip3, FILE *fserca, FILE *fncx, FILE *fvocc,  FILE *fbkca, FILE *fip3r, FILE *fcacc,  FILE *fnsc, FILE *fsrca, FILE *fryr, FILE *fdag, double tnow)" function and
 * "void dump_data_timeseries(FILE *fvtk, double tnow)" function
 */

extern double *y, *par;
extern int *num_cell, *num_par, *num_var;


void dump_data(FILE *ft, FILE *fvm, FILE *fca, FILE *fip3, FILE *fncx, FILE *fvocc,  FILE *fbkca, FILE *fip3r, FILE *fsrca, FILE *fryr, FILE *fdag, double tnow)
	{
	fprintf(ft,"%f\t\n",tnow);
	for (int i=0; i<num_cell[0]; i++)
	{
			fprintf(fvm,"%f\t",y[i * num_var[0] + smc_vm]);
			fprintf(fca,"%f\t",y[i * num_var[0] + smc_ca]);
			fprintf(fip3,"%f\t",y[i * num_var[0] + smc_ip3]);
			fprintf(fdag,"%f\t",y[i * num_var[0] + smc_DAG]);
			fprintf(fsrca,"%f\t",y[i * num_var[0] + smc_SR_ca]);
			fprintf(fncx,"%f\t",par[i * num_par[0] + I_ncx]*1e3);
			fprintf(fvocc,"%f\t",par[i * num_par[0] + I_vocc]*1e3);
			fprintf(fbkca,"%f\t",par[i * num_par[0] + I_bkca]*1e3);
			fprintf(fip3r,"%f\t",par[i * num_par[0] + I_ip3r]*1e3);
			fprintf(fryr,"%f\t",par[i * num_par[0] + I_ryr]*1e3);
	}

		fprintf(fvm,"\n");
		fprintf(fca,"\n");
		fprintf(fvocc,"\n");
		fprintf(fbkca,"\n");
		fprintf(fip3,"\n");
		fprintf(fdag,"\n");
		fprintf(fip3r,"\n");
		fprintf(fsrca,"\n");
		fprintf(fncx,"\n");
		fprintf(fryr,"\n");
	}

