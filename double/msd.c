/*
msd.c

Calculated the MSD using the files generated from the double_hinge.c
Outputs a file 'measures.dat' in the 'data/' folder.
Each line of 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define fileName 		"data/ang_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"
#define fileNameSave 	"data/measures.dat"

/* Calculating the MSD */
void calculate_MSD(double *pos, double *MSD, int *cnt, int tabs, int mode) {
	int t;

	if (mode==0) { // calculating ETAMSD

		int t0;
		double auxMSD[tabs];
		int auxCnt[tabs];

		for (t=0;t<tabs;t++) {
			auxMSD[t] = 0.;
			auxCnt[t] = 0;
		}
		for (t0=0;t0<tabs;t0++) {
			for (t=t0;t<tabs;t++) {
				auxMSD[t-t0] += pow(pos[t]-pos[t0],2);
				auxCnt[t-t0]++;
			}
		}
		for (t=0;t<tabs;t++) {
			MSD[t] += auxMSD[t]/(1.*auxCnt[t]);
			cnt[t]++;
		}

	} else { // calculating EAMSD

		for (t=0;t<tabs;t++) {
			MSD[t] += pow(pos[t]-pos[0],2);
			cnt[t]++;
		}

	}
	
}

void main() {

	int dt = 1000;
	int tmax = 50000000/dt;

	double dtsim = 0.00005;

	double ba = M_PI/6.;

    double rho[6] = {0.5,0.55,0.6,0.65,0.7,0.75};
    double ang[2] = {1.5*ba,2.*ba};
    double R[1] = {2.0};
    double F[4] = {0.,30.,60.,90.};
    double D[1] = {200.};

    double pos[2*tmax],pos_1[tmax],pos_2[tmax];
    pos[0] = 0.; pos_1[0] = 0.; pos_2[0] = 0.;

    /*
    MSD: for the whole dynamics;
    MSD_1: up until first rod is abdorbed;
	MSD_2: from absorption of first rod until absorption of the second;

	MSD: Ensemble Average MSD;
	MSD0 : Ensemble Time Average MSD;
    */
    double MSD[2*tmax],MSD_1[tmax],MSD_2[tmax];
    double MSD0[2*tmax],MSD0_1[tmax],MSD0_2[tmax];
    int cntMSD[2*tmax],cntMSD_1[tmax],cntMSD_2[tmax];
    int cntMSD0[2*tmax],cntMSD0_1[tmax],cntMSD0_2[tmax];

    char name[100];
    
    FILE *inFile; // In file
    FILE *outFile = fopen(fileNameSave, "w"); // out file

    int i, j, k, l, m;
    int sim, t, t_tot, tabs, t1abs, t2abs;
    double dx, dx_tot;

    for (i=0;i<6;i++) {
        for (j=0;j<2;j++) {
            for (k=0;k<1;k++) {
                for (l=0;l<4;l++) {
                	for (m=0;m<1;m++) {

                		for (t=0;t<2*tmax;t++) {
	                		MSD[t] = 0.;
						    MSD0[t] = 0.;
						   	cntMSD[t] = 0;
						    cntMSD0[t] = 0;
						    if (t<tmax) {
							    MSD_1[t] = 0.;
							    MSD0_1[t] = 0.;
							    cntMSD_1[t] = 0;
							    cntMSD0_1[t] = 0;
							    MSD_2[t] = 0.;
							    MSD0_2[t] = 0.;
							    cntMSD_2[t] = 0;
							    cntMSD0_2[t] = 0;
							}
						}

						t1abs = 0; t2abs = 0;

	                	sprintf(name,fileName,10,10,R[k],rho[i],F[l],ang[j],(int)D[m]);
	                	inFile = fopen(name, "r");
	                    while (fscanf(inFile, "%d %*d %d %d %*lf %*d %*lf %*d %lf %lf",&sim,&t,&t_tot,&dx,&dx_tot) != EOF) {

	                    	if (t_tot==dt&&sim!=0) {
	                    		calculate_MSD(pos,  MSD,   cntMSD,   tabs, 0);
	                    		calculate_MSD(pos,  MSD0,  cntMSD0,  tabs, 1);
	                    		calculate_MSD(pos_1,MSD_1, cntMSD_1, t1abs,0);
	                    		calculate_MSD(pos_1,MSD0_1,cntMSD0_1,t1abs,1);
	                    		calculate_MSD(pos_2,MSD_2, cntMSD_2, t2abs,0);
	                    		calculate_MSD(pos_2,MSD0_2,cntMSD0_2,t2abs,1);
	                    		t1abs = 0; t2abs = 0;
	                    	}

	                    	if (t==t_tot) {
	                    		pos_1[(int)(t/dt)] = dx;
	                    		t1abs = (int)(t/dt);
	                    	} else {
	                    		pos_2[(int)(t/dt)] = dx;
	                    		t2abs = (int)(t/dt);
	                    	}
	                    	pos[t1abs+t2abs] = dx_tot;
	                    	tabs = t1abs + t2abs;

					    }
					    fclose(inFile);
					    
					    calculate_MSD(pos,MSD,cntMSD,tabs,0);
	                    calculate_MSD(pos,MSD0,cntMSD0,tabs,1);
	                    calculate_MSD(pos_1,MSD_1,cntMSD_1,t1abs,0);
	                    calculate_MSD(pos_1,MSD0_1,cntMSD0_1,t1abs,1);
	                    calculate_MSD(pos_2,MSD_2,cntMSD_2,t2abs,0);
	                    calculate_MSD(pos_2,MSD0_2,cntMSD0_2,t2abs,1);

					    // Correctly normalize everything
					    for (t=0;t<2*tmax;t++) {
					    	if (cntMSD[t])		MSD[t] /= 1.*cntMSD[t];
					    	if (cntMSD0[t])		MSD0[t] /= 1.*cntMSD0[t];
					    	if (t<tmax) {
						    	if (cntMSD_1[t])	MSD_1[t] /= 1.*cntMSD_1[t];
						    	if (cntMSD0_1[t])	MSD0_1[t] /= 1.*cntMSD0_1[t];
						    	if (cntMSD_2[t])	MSD_2[t] /= 1.*cntMSD_2[t];
						    	if (cntMSD0_2[t])	MSD0_2[t] /= 1.*cntMSD0_2[t];
						    }
					    }

					    for (t=0;t<2*tmax;t++)
	                		fprintf(outFile, "%.2f %.3f %d %.1f %d %f %f %f %f %f %f %f\n",
	                			rho[i],ang[j],(int)F[l],R[k],(int)D[m],dt*t*dtsim,
	                			MSD[t],(t<tmax?MSD_1[t]:0.),(t<tmax?MSD_2[t]:0.),
	                			MSD0[t],(t<tmax?MSD0_1[t]:0.),(t<tmax?MSD0_2[t]:0.));

                	}
                }
            }
        }
    }  	

    fclose(outFile);

}