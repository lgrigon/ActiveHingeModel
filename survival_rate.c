/*
survival_rate.c

Function to calculate the survival probability of the rod,
which will be used to estimate the MFPT.

Input: measurement files (obtained from running single_hinge.c)
Output: survival rate files
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define T_MAX 	20000000
#define N_SKIP 	1000
#define N_MEA 	T_MAX/N_SKIP
#define DT 		0.00005
#define KBT 	10.0

#define FILE_NAME_MEA	"data/mea_L%.1f-pf%.2f-F%.1f-ang%.3f-mu%.2f-T%.1f.dat"
#define FILE_NAME_SUR 	"mfpt/surv_L%.1f-pf%.2f-F%.1f-ang%.3f-mu%.2f-T%.1f.dat"

/** Calculating the survival rate **/
void surv_rate(double R, double rho, double PE, double ang, double mu) {

	char name[100];
    sprintf(name,FILE_NAME_MEA,R,rho,PE,ang,mu,KBT);

	if (access(name, F_OK) == 0) {
		int i, j, n;

		double t1, t2, tb;
		int t1i, t2i, tbi, nabs;

		double surv_rate[2*N_MEA][3]; // It has to be 2*N_MEA to account for both rod's absorptions dynamics, the 3 comes from the 3 FPT (first rod, second rod, both rods)
		for (i=0;i<2*N_MEA;i++)
		for (j=0;j<3;j++)
    	surv_rate[i][j] = 0.;

		FILE *inFile;
		inFile = fopen(name, "r");
		int abs[2] = {0,0};
		double meanb[2] = {0.,0.};
		int norm = 0;
		while (fscanf(inFile, "%*d %lf %*lf %*lf %lf %*lf %*lf %lf %*lf %*lf %d",&t1,&t2,&tb,&nabs) != EOF) {

			// Normalization based on nabs
			norm++;
			if (nabs>=1) { // If at least one rod was absorbed
				abs[0]++; // Increase counter for 1st rod
				meanb[0] += t1;
				if (nabs==2) { // If second was also absorbed
					abs[1]++;
					meanb[1] += t2;
				}
			}

			t1i = (int) (t1/DT);
			t1i = (int) round(t1i/N_SKIP); // FPT first rod

			t2i = (int) (t2/DT);
			t2i = (int) round(t2i/N_SKIP); // FPT second rod

			tbi = (int) (tb/DT);
			tbi = (int) round(tbi/N_SKIP); // FPT both rods

			// Account for the t=0 (i.e. if one/both particle was not absorbed)
			if (nabs==0) t1i = N_MEA;
			if (nabs==1) t2i = N_MEA;

			// Calculating the survival rate
			for (i=0;i<t1i;i++) surv_rate[i][0] += 1.;
			for (i=0;i<t2i;i++) surv_rate[i][1] += 1.;
			for (i=0;i<tbi;i++) surv_rate[i][2] += 1.;

		}
		fclose(inFile);

		// calculate error using jacknife
		double srj,err_sr[2*N_MEA];
		sprintf(name,FILE_NAME_MEA,R,rho,PE,ang,mu,KBT);
		inFile = fopen(name, "r");
		for (i=0;i<2*N_MEA;i++) err_sr[i] = 0.;
		while (fscanf(inFile, "%*d %lf %*lf %*lf %lf %*lf %*lf %lf %*lf %*lf %d",&t1,&t2,&tb,&nabs) != EOF) {
			t2i = (int) (t2/DT);
			t2i = (int) round(t2i/N_SKIP);
			for (i=0;i<2*N_MEA;i++) {
				if (i<t2i) srj = surv_rate[i][1] - 1.;
				else srj = surv_rate[i][1];
				srj /= 1.*abs[0]-1.;

				err_sr[i] += (srj - surv_rate[i][1]/(1.*abs[0])) * (srj - surv_rate[i][1]/(1.*abs[0]));
			}
		}
		for (i=0;i<2*N_MEA;i++) err_sr[i] = sqrt(err_sr[i] * (1.*abs[0] - 1.) / abs[0]);
		fclose(inFile);

		// Normalization
		for (i=0;i<2*N_MEA;i++) {
			surv_rate[i][0] /= 1.*norm;
			if (abs[0]>0) surv_rate[i][1] /= 1.*abs[0];
			surv_rate[i][2] /= 1.*norm;
		}

		// To help eliminate some initial termalization effects on the fit
		/*double tini[2] = {0.,0.};
		for (i=0;i<2*N_MEA;i++) {
			if (surv_rate[i][0] == 1.) tini[0] = 1.*i*N_SKIP*DT;
			if (surv_rate[i][1] == 1.) tini[1] = 1.*i*N_SKIP*DT;
		}*/

		// Correcting the "prediction" of MFPT for the fit
		// A rough estimation of MFPT to help the exponential fit converge faster
		/*meanb[0] = 1.*(meanb[0] + 1.*T_MAX*DT*(norm-abs[0]))/norm;
		meanb[1] = (abs[0]>0 ? 1.*(meanb[1] + 1.*T_MAX*DT*(abs[0]-abs[1]))/abs[0] : 0.);*/

		// Do we have enough data?
		double min_sr = 0.9;
		int enough_data;
		enough_data = (surv_rate[N_MEA-1][1] < min_sr ? 1 : 0);
		
		// Saving
		FILE *outFile;
		sprintf(name,FILE_NAME_SUR,R,rho,PE,ang,mu,KBT);
		outFile = fopen(name,"w");
		for (i=0;i<2*N_MEA;i++) {
			/*fprintf(outFile, "%f %f %f %f ", 1.*i*N_SKIP*DT, surv_rate[i][0], surv_rate[i][1], surv_rate[i][2]);
			fprintf(outFile, "%d %d %d ", norm, abs[0], abs[1]);
			fprintf(outFile, "%f %f %f %f\n", tini[0], meanb[0], tini[1], meanb[1]);*/
			/*fprintf(outFile, "%f %f ", 1.*i*N_SKIP*DT, surv_rate[i][1]);
			fprintf(outFile, "%d ", abs[1]);
			fprintf(outFile, "%f %f %f\n", tini[1], meanb[1], err_sr[i]);*/
			fprintf(outFile, "%f %f ", 1.*i*N_SKIP*DT, surv_rate[i][1]);
			fprintf(outFile, "%f %d %d\n", err_sr[i], abs[1], enough_data);

		}
		fclose(outFile);

	}

}

void main() {

	int i, j, k, l, m;

	double ba = M_PI/6.;

    double rho[15] = {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75};  // Packing fractions
    double ang[2] = {1.5*ba,2.*ba}; // Maximum angles
    double R[1] = {2.0}; // Rod's size
    double F[4] = {0.,30.,60.,90.}; // Magnitude of active force
    double mu[1] = {0.05}; // Diffusion coefficient

    /** SURVIVAL RATE **/
    for (i=0;i<15;i++)
    for (j=0;j<2;j++)
    for (k=0;k<1;k++)
    for (l=0;l<4;l++)
    for (m=0;m<1;m++)
    	surv_rate(R[k],rho[i],F[l],ang[j],mu[m]);

}