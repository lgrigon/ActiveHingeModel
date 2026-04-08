#include <stdio.h>
#include "params.h"
#include "particles.h"
#include "random.h"
#include "utils.h"
#include <math.h>

int main() {
    
    // Initialize random number generator
    srand(time(0));
	myrand = SEED;
    Init_Random();

    // Particle arrays
    passive_particle p[N_PASSIVE];	// passive particles
	passive_particle ll[N_BEADS];	// left active rod
	passive_particle lr[N_BEADS];	// right active rod

    // Cell lists
    int head[LA_X*LA_Y]; // head of the linked list for passive particles
	int lscl[N_PASSIVE]; // linked list for passive particles

	// Output files
	char name[100];
	FILE *outFileDet;
	if (CHECK_SAVE_DET) {
	  sprintf(name,FNAME_DET,L_ROD,PF,F_ACT,ANG_MAX,RATIO_MU*MU_BEADS,MU_PASSIVE,KBT);
	  outFileDet = fopen(name, "a");
	}
	FILE *outFileMea;
	sprintf(name,FNAME_MEA,L_ROD,PF,F_ACT,ANG_MAX,RATIO_MU*MU_BEADS,MU_PASSIVE,KBT);
	outFileMea = fopen(name, "a");

    int t, t_tot, sim;
	double dx, dx_tot, xant;
	double tabs[2], pos[2], spd[2];
	int nabs, Ns[2];

    sim = 0;
	do {

        // Initialize particles positions and orientations
		init_particles(p,ll,lr,head,lscl);

		dx = 0.; dx_tot = 0.; xant = ll[0].x;
		for (t=0;t<2;t++) {
			tabs[t] = 0.;
			pos[t]  = 0.;
			spd[t]  = 0.;
		}
		nabs = (CHECK_FROM_SECOND ? 1 : 0);
		t = 0;
		
		do {

			fill_list(p,head,lscl); // fill the cell lists
			update(p,ll,lr,1,head,lscl); // update particles positions

			if (t%N_SKIP==0) { 
				if (CHECK_LINUX_PLOT) plot_config(ll,lr,p);
				if (CHECK_SAVE_DET) {
					get_N(Ns,p,ll,lr);
					fprintf(outFileDet,"%d %d %d %d %f %d %f %d %f %f\n",sim,nabs,t,t_tot,ang[0],Ns[0],ang[1],Ns[1],dx,dx_tot);
				}
			}

			// Total displacement of the hinge
			dx += periodic(ll[0].x - xant);
			dx_tot += periodic(ll[0].x - xant);
			xant = ll[0].x;

			// If one of the rods is absorbed
			if ((ang[0] == 0. || ang[1] == 0.) && nabs == 0) {
				tabs[0] = 1.*t*DT; // Save FPT
				pos[0] = dx; // Save total hinge displacement up to absorption
				spd[0] = 1.*pos[0]/tabs[0]; // Save hinge's average speed
				nabs++;
				t = 0;
				dx = 0.;
			}

			// If the second rod is absorbed
			if (ang[0] == 0. && ang[1] == 0.) {
				tabs[1] = 1.*t*DT; // Save FPT
				pos[1] = dx; // Save total hinge displacement
				spd[1] = 1.*pos[1]/tabs[1]; // Save hinge's average speed
				nabs++;
			}

			t++;
			t_tot++;

		} while (t<T_MAX && nabs<2);

		/* Saving basic measurements into file */
		/* Saving: 	simulation id, FPT of FIRST rod, displacement of the hinge until FIRST rod was absorbed, average speed of hinge until FIRST rod was absorbed,
								FPT of SECOND rod, displacement of the hinge until SECOND rod was absorbed, average speed of hinge until SECOND rod was absorbed, 
								total simulation time, total displacement of the hinge, average hinge speed, number of absorbed rods
		*/
		fprintf(outFileMea,"%d %f %f %f %f %f %f ",sim,tabs[0],pos[0],spd[0],tabs[1],pos[1],spd[1]);
		fprintf(outFileMea,"%f %f %f %d\n",1.*t_tot*DT,dx_tot,1.*dx_tot/(1.*t_tot*DT),nabs);

		sim++;

	} while (sim<N_SIM);

	fclose(outFileMea);
	fclose(outFileDet);


    return 0;
}
