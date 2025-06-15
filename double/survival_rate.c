/*
survival_rate.c

Function to calculate the survival probability of the rod,
which will be used to estimate the MFPT.

Input: measurement files (obtained from running double_hinge.c)
Output: survival rate files
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define T_MAX 	50000000
#define N_SKIP 	1000
#define N_MEA 	T_MAX/N_SKIP
#define DT 		0.00005

#define FILE_NAME_MEA	"data/mea_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"
#define FILE_NAME_SUR 	"mfpt/surv_R%.1f-pf%.2f-f%d-ang%.3f-D%d.dat"

/** Calculating the survival rate **/
void surv_rate(double R, double rho, double F, double ang, int D) {
	int i, j, n;

	double t1, t2, tb;
	int t1i, t2i, tbi, nabs;

	char name[100];

	double surv_rate[2*N_MEA][3]; // It has to be 2*N_MEA to account for both rod's absorptions dynamics, the 3 comes from the 3 FPT (first rod, second rod, both rods)
	for (i=0;i<2*N_MEA;i++)
    for (j=0;j<3;j++)
    	surv_rate[i][j] = 0.;

    FILE *inFile;
    sprintf(name,FILE_NAME_MEA,10,10,R,rho,F,ang,D);
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

	// Normalization
	for (i=0;i<2*N_MEA;i++) {
    	surv_rate[i][0] /= 1.*norm;
    	if (abs[0]>0) surv_rate[i][1] /= 1.*abs[0];
    	surv_rate[i][2] /= 1.*norm;
	}

	// To help eliminate some initial termalization effects on the fit
	double tini[2] = {0.,0.};
    for (i=0;i<2*N_MEA;i++) {
    	if (surv_rate[i][0] == 1.) tini[0] = 1.*i*N_SKIP*DT;
    	if (surv_rate[i][1] == 1.) tini[1] = 1.*i*N_SKIP*DT;
    }

    // Correcting the "prediction" of MFPT for the fit
    // A rough estimation of MFPT to help the exponential fit converge faster
    meanb[0] = 1.*(meanb[0] + 1.*T_MAX*DT*(norm-abs[0]))/norm;
    meanb[1] = (abs[0]>0 ? 1.*(meanb[1] + 1.*T_MAX*DT*(abs[0]-abs[1]))/abs[0] : 0.);
	
	// Saving
	FILE *outFile;
	sprintf(name,FILE_NAME_SUR,R,rho,(int)F,ang,D);
    outFile = fopen(name,"w");
	for (i=0;i<2*N_MEA;i++) {
	    fprintf(outFile, "%f %f %f %f ", 1.*i*N_SKIP*DT, surv_rate[i][0], surv_rate[i][1], surv_rate[i][2]);
	    fprintf(outFile, "%d %d %d ", norm, abs[0], abs[1]);
	    fprintf(outFile, "%f %f %f %f\n", tini[0], meanb[0], tini[1], meanb[1]);
	}
	fclose(outFile);

}

void main() {

	int i, j, k, l, m;

	double ba = M_PI/6.;

    double rho[6] = {0.5,0.55,0.6,0.65,0.7,0.75}; // Packing fractions
    double ang[2] = {1.5*ba,2.*ba}; // Maximum angles
    double R[1] = {2.0}; // Rod's size
    double F[4] = {0.,30.,60.,90.}; // Magnitude of active force
    int D[1] = {200}; // Diffusion coefficient

    /** SURVIVAL RATE **/
    for (i=0;i<6;i++)
    for (j=0;j<2;j++)
    for (k=0;k<1;k++)
    for (l=0;l<4;l++)
    for (m=0;m<1;m++)
    	surv_rate(R[k],rho[i],F[l],ang[j],D[m]);

}