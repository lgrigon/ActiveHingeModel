/*
terms_surv_rate.c

Function to calculate the survival probability of each coarse-grained state N until absorption of the rod,
which will be used to estimate the T_N terms of the MFPT.

Input: measurement files (obtained from running single_hinge.c)
Output: survival rate file
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define T_MAX 	50000000
#define N_SKIP 	1000
#define N_MEA 	T_MAX/N_SKIP
#define DT 		0.00005
#define N_SIM	500

#define FILE_NAME_DETAILS	"data/ang_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"
#define FILE_NAME_TERMS	 	"mfpt/mfpt_terms_sr.dat"

void get_MFPT(char *fname, double rho, double ang, int F, double R, int D) {

	// Load data
	int sim, nR, t, n;

	double tn[N_SIM][4];
	for (n=0;n<4;n++)
		for (sim=0;sim<N_SIM;sim++)
			tn[sim][n] = 0.;

	FILE *inFile = fopen(fname, "r");
	while (fscanf(inFile, "%d %*d %*d %d %*lf %*d %*lf %d %*lf %*lf",&sim,&t,&nR) != EOF) {
		tn[sim][nR] += 1.;
	}
	fclose(inFile);

	double sr_tn[5][N_MEA];
	for (n=0;n<5;n++)
		for (t=0;t<N_MEA;t++)
			sr_tn[n][t] = 0.;

	for (sim=0;sim<N_SIM;sim++) {
		for (n=0;n<4;n++) {
			for (t=0;t<tn[sim][n];t++)
				sr_tn[n][t] += 1.;
		}
		for (t=0;t<(tn[sim][0]+tn[sim][1]+tn[sim][2]+tn[sim][3]);t++)
			sr_tn[4][t] += 1.;
	}

	for (n=0;n<5;n++) 
		for (t=0;t<N_MEA;t++)
			sr_tn[n][t] /= 1.*N_SIM;

	FILE *outFile = fopen(FILE_NAME_TERMS, "a");
	for (n=0;n<5;n++) 
		for (t=0;t<N_MEA;t++)
			fprintf(outFile,"%.2f %.3f %d %.1f %d %d %f %f\n",rho,ang,F,R,D,n,1.*t*N_SKIP*DT,sr_tn[n][t]);
	fclose(outFile);
	
}

void main() {

	double ba = M_PI/6.;

    double rho[1] = {0.75};
    double ang[1] = {2.0*ba};
    double R[1] = {2.0};
    double F[7] = {0.,15.,30.,45.,60.,75.,90.};
    double D[1] = {200.};

    double mfpt;
   
    char name[100];

    /** TRANSITION RATES **/
    int i,j,k,l,m;
    for (i=0;i<1;i++) {
        for (j=0;j<1;j++) {
            for (k=0;k<1;k++) {
                for (l=0;l<7;l++) {
                	for (m=0;m<1;m++) {
                		sprintf(name,FILE_NAME_DETAILS,10,10,R[k],rho[i],F[l],ang[j],(int)D[m]);
                		get_MFPT(name,rho[i],ang[j],(int)F[l],R[k],(int)D[m]);
                	}
                }
            }
        }
    }

}