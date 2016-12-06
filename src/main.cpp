/**
 * \mainpage Vascular Smooth Muscle Cell (VSMC): Coupled cells Model
 * The coupled cells model is used to simulate calcium dynamics in a population of VSMCs. 
 * The gap junctional intercellular communication (GJIC) is modelled by adding coupling equations for IP$_3$ and V$_m$.
 * A linear coupling equation is assigned for both IP$_3$ and V$_m$.
 * Calcium diffusion between the cells is neglected.
 * A sigmoidal distriubtion of agonist concentration is applied over the 100 coupled cells.
 * @author Jaijus
 */


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include<sys/stat.h>

#include "include/computelib.h"

int main(int argc, char* argv[]){

	/* SMCs coupling parameters */
	    double Exp_G_gj = 178.57; // Gap junction conductance - base value, nS
		double unitary_G_gj = 65; // Unitary conductance - average value, pS
		double number_Gj_base = Exp_G_gj*1e3/unitary_G_gj; // Total number of gap junctions - base value
		double Exp_IP3 = 2;  // IP3 coupling coeffient -  base value assumed, s-1
		double unitary_IP3 = Exp_IP3/number_Gj_base; // Unitary IP3 coupling coefficient -  base value calculated, s-1
		double number_G_gj = number_Gj_base; // Total number of gap junctions


	/* Sigmoidal function */
		double xmean = 0.5; // sigmoidal distribution mean position
		double xstretch_var = 0.1; // sigmoidal function: slope

	/* Simulation controlling parameters */
		double tnow = 0.0; // Initial time , s
		double tfinal = 201; // Final time, s
		double interval = 1e-4; // Time interval, s
		double tol_state = 1e-4; // Tolerance limit for the convergence in the iteration
		double delay = 50; // Time allowed for variables to converge
		double stim_time = tfinal; // Total stimulation time, s
		double stim = delay+stim_time; // Time when stimulation ends, s
		double count_stim = 0; // control variable for the loop
		double count_delay = 0; // control variable for the loop
		int tol_count = 0; // Control varialbe for fixed point iteration
		double file_write_freq = 0.01; // File writing frequency, s
		int folder_status; // variable to check the successful creartion of folder


		char main_folder[300];
		char sub_folder[300];
		char temp_name[300];

		/* Variables for coupled cells */
		num_var = (int *) malloc(sizeof(int)); // Total number of varibles
		num_par = (int *) malloc(sizeof(int)); // Total number of parameters
		num_cell = (int *) malloc(sizeof(int)); // Total number of cells

		num_var[0] 	= 18; // Total number of state variables in a cell
		num_par[0] 	= 41; // Total number of parameters in a cell

		int var_tot; // total number state variables in all the cells
		int par_tot; // total number of parameters in all the cells

		num_cell[0]	= 100; // Total number of SMCs

		var_tot	= num_cell[0] * num_var[0];
	 	par_tot	= num_cell[0] * num_par[0];

	 	allocatememory(var_tot, par_tot); // Allocating memory

		G_gj_base[0] = number_G_gj * unitary_G_gj/1e3; // Base value of gap junctional conductance
		IP3_P_base[0] = unitary_IP3 * number_G_gj; // Base value of IP3 coupling coefficient

		G_gj_varied[0] = number_G_gj * unitary_G_gj/1e3; // value of gap junctional conductance
		IP3_P_varied[0] = unitary_IP3 * number_G_gj;  // value of IP3 coupling coefficient
		L_rest[0] = 0.0; // initial agonist concentration- before assigning sigmoidal distribution


		/* Creating main folder */
		sprintf(main_folder,"%3.1fs_%1.4fs",tfinal,interval);
		folder_status = mkdir(main_folder,0777);
			if(folder_status == -1)
				{
				perror("Couldn't create sub-directory");
				return -1;
				}


			/* Creating subfolder to save all the files for each loop */
			sprintf(sub_folder,"%s/%d_%3.2f_%3.1fnM_%3.1fnM_%4.1fpS_%3.1fs-1",main_folder,num_cell[0],xstretch_var,nonstim_NE,stim_NE,G_gj_varied[0]/1e-3,IP3_P_varied[0]);
			printf("%s\n",sub_folder);
			folder_status = mkdir(sub_folder,0777);
				if(folder_status == -1)
					{
					perror("Couldn't create sub-directory");
					return -1;
					}


			int count = 1; // counter to produce output files
			// OIpening and creating FILE names for output files
			FILE *fca; // Time series of intracellular calcium
			FILE *fvm; // Time series of membrane potential
			FILE *ft; // Time data
			FILE *fbkca; // Time series of BKCa channel flux
			FILE *fvocc; // Time series of BKCa channel flux
			FILE *fdag; // Time series of DAG concentration
			FILE *fip3; // Time series of IP3 concentration
			FILE *fip3r; // Time series of IP3R channel flux
			FILE *fsrca; // Time series of intraSR calcium
			FILE *fncx; // Time series of NCX channel flux
			FILE *fryr; // Time series of RyR channel flux
			FILE *fs; // Cell spatial position
			FILE *fL; // Spatial distribution of agonist concentration
			FILE *fGgj; // Spatial distribution of gap junctional conductance
			FILE *fIP3P; // Spatial distribution of IP3 coupling coefficient


			sprintf(temp_name,"%s/Ca_timeseries.txt",sub_folder);
			fca = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/Vm_timeseries.txt",sub_folder);
			fvm = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/space.txt",sub_folder);
			fs = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/time.txt",sub_folder);
			ft = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_BKCa.txt",sub_folder);
			fbkca = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_VOCC.txt",sub_folder);
			fvocc = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/DAG_timeseries.txt",sub_folder);
			fdag = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/IP3_timeseries.txt",sub_folder);
			fip3 = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_IP3R.txt",sub_folder);
			fip3r = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/SR_Ca_timeseries.txt",sub_folder);
			fsrca = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/flux_NCX.txt",sub_folder);
			fncx = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/flux_RyR.txt",sub_folder);
			fryr = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/NE.txt",sub_folder);
			fL = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/G_gj.txt",sub_folder);
			fGgj = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/IP3_P.txt",sub_folder);
			fIP3P = fopen(temp_name,"w+");


			initialize(fs); // Initialize the variables
			dump_data(ft, fvm, fca, fip3, fncx, fvocc, fbkca, fip3r, fsrca, fryr,fdag, tnow); // Creating output files

			// Simulation starts here
			while (tnow <= tfinal)
			{

				// Assign agonist concentration of the cells to the value before the stimulation
				if (tnow >= stim && count_stim < 1)
				{
					count_stim = 2;
					for (int i = 0; i < num_cell[0]; i++)
						{
						L_cell[i] = L_rest[0];
						}
				}

				// Assigning spatial distribution of agonist concentration to the cells - sigmoidal distribution
				else if (tnow >= (delay) && count_delay < 1)
				{
					count_delay = 2;
					L_cell[0] = ((1+tanh((double(0)/(num_cell[0]-1)-xmean)/xstretch_var))/2); //
					L_cell[num_cell[0]-1] = ((1+tanh((1-xmean)/xstretch_var))/2)-L_cell[0]; //
					for (int i = 0; i < num_cell[0]; i++)
						{
						fprintf(fL,"%d\t%f\t",i,L_cell[i]);
						fprintf(fGgj,"%d\t%f\n",i,G_gj[i]);
						fprintf(fIP3P,"%d\t%f\n",i,IP3_P[i]);
						if(i==0)
							fprintf(fL,"%f\n",stim_NE);
						else if(i == (num_cell[0]-1))
							fprintf(fL,"%f\n",nonstim_NE);
						else
							{
							L_cell[i] = stim_NE-(stim_NE-nonstim_NE)*(((1+tanh((double(i)/(num_cell[0]-1)-xmean)/xstretch_var))/2)-L_cell[0])/L_cell[num_cell[0]-1]; //
							fprintf(fL,"%f\n",L_cell[i]);
							}
						}
					L_cell[0] = stim_NE;
					L_cell[num_cell[0]-1] = nonstim_NE;

					fclose(fL);
					fclose(fGgj);
					fclose(fIP3P);
				}


				while (tol_count<1)
				{
					tol_count = 2;
					if (tnow<delay)
					{
						singlecell(tnow, interval); // Single cell simulation to all cells to converge variables
					}
					else
					{
						multicell(tnow, interval); // Coupled cells simulation
					}

				/* Fixed point iteration checking*/
					for (int i=0;i<num_cell[0];i++)
					for (int j=0;j<num_var[0];j++)
						if (fabs(y_temp[i * num_var[0] + j] - y_converge[i * num_var[0] + j])>=tol_state)
						{
							tol_count = 0;
							/* Not converged */
							for (int k=0;k<num_cell[0];k++)
								for (int m=0;m<num_var[0];m++)
										y_converge[k * num_var[0] + m] = y_temp[k * num_var[0] + m];
							break;
						}
				}

				/* Converged*/
				/* Assigning converged results as the values at the new time step */
				for (int i=0;i<num_cell[0];i++)
					for (int j=0;j<num_var[0];j++)
						{
							y[i * num_var[0] + j] = y_temp[i * num_var[0] + j];
							y_converge[i * num_var[0] + j] = y_temp[i * num_var[0] + j];
						}

				tol_count = 0;

				tnow += interval; // new time step

				//  Checking file writing condition
					if (int(tnow/file_write_freq) == count)
						{
						dump_data(ft, fvm, fca, fip3, fncx, fvocc, fbkca, fip3r,  fsrca, fryr, fdag, tnow);
						count++;
						}

			}

			/* Closing all the files opened for writing */
			fclose(ft);
			fclose(fvm);
			fclose(fca);
			fclose(fip3);
			fclose(fsrca);
			fclose(fncx);
			fclose(fvocc);
			fclose(fbkca);
			fclose(fip3r);
			fclose(fdag);


	// free up memory
			free(y);
			free(y_temp);
			free(y_converge);
			free(par);
			free(num_var);
			free(num_par);
			free(num_cell);
			free(TDMB);
			free(TDMD);
			free(TDMDD);
			free(TDMCC);
			free(TDMB_IP3);
			free(TDMD_IP3);
			free(TDMDD_IP3);
			free(TDMCC_IP3);
			free(G_gj);
			free(G_gj_base);
			free(G_gj_varied);
			free(IP3_P_base);
			free(IP3_P_varied);
			free(IP3_P);
			free(L_cell);
			free(L_rest);
			free(Q_ip3r);
			free(Q_serca);
			free(Q_ryr);

}
