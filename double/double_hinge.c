/*
double_hinge.c

This algorithm simulates the active hinge model.
Outputs data from the simulations into files at 'data/' folder.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

/** SIMULATION PARAMETERS **/

// System size
#define LX				10.
#define LY				10.

// Active particles (hinge)
#define R 		 		2.0 // Hinge size
#define F_ACT		 	90. // Magnitude of active force
#define ANG_MAX		2.*M_PI/6. // Maximum angle

// Passive particles
#define PF			 0.65 // Packing fraction of passive particles

// Run parameters
#define T_EQ			500000 // Some steps for termalization of passive particles
#define T_MAX			20000000 // Total simulation steps
#define N_SKIP		1000 // Number of steps between measurements
#define N_SIM			1000 // Number of simulations

#define USE_ACT_FRAME 	0 // Use the hinge's frame of reference 

/** FILES **/
#define SAVE_CONFIGS		1 // Do we want to also save the configurations?
#define SAVE_DETAILS		1 // Do we want to save detailed information?
#define FILE_NAME_CONF	"data/conf_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d-t%d-N%d_%d.dat"
#define FILE_NAME_ANG		"data/ang_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"
#define FILE_NAME_MEA		"data/mea_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"

/** FIXED PARAMETERS **/

// System (General)
#define DT				0.00005

#define MU				0.05 // Motility coefficient
#define D					200. // Diffusion coefficient

#define K_WCA				50 // WCA potential strenght

// Active particles
#define R_LEVER		0.1 // Radius of active rod's particles
#define N_SUB			(int) round(2*R/R_LEVER) // Number of those particles
#define MU_R			0.01 // Rotational mobility coefficient 
#define D_R				0. // Rotational diffusion coefficient

// Passive particles
#define R_DISK   	0.5 // Passive particles radius
#define N        	(int) round(PF*LX*LY/(M_PI*R_DISK*R_DISK)) // Number of active particles
#define PADDING		0.9

// For linked cells
#define LA_X          (int) floor(LX/(2.*R_DISK))
#define RA_X          1.*LY/LA_X
#define LA_Y          (int) floor(LY/(2.*R_DISK))
#define RA_Y          1.*LY/LA_Y

/** STRUCT FOR THE PASSIVE PARTICLES **/
typedef struct {
	double x, y;
} passive_particle;

double ang[2]; // Global variable for the angle of each rod

/*---------------------------------------------------------------------*/
/*----------P. Random Number Generator by Parisi & Rapuano-------------*/
/*---------------------------------------------------------------------*/
#define SEED time(NULL)
#define FNORM   (2.3283064365e-10)
#define RANDOM  ((ira[ip++] = ira[ip1++] + ira[ip2++]) ^ ira[ip3++])
#define FRANDOM (FNORM * RANDOM)
#define pm1 ((FRANDOM > 0.5) ? 1 : -1)

unsigned myrand, ira[256];
unsigned char ip, ip1, ip2, ip3;

unsigned rand4init(void) {
  unsigned long long y;

  y = (myrand*16807LL);
  myrand = (y&0x7fffffff) + (y>>31);
  if (myrand&0x80000000)
    myrand = (myrand&0x7fffffff) + 1;
  return myrand;
}

void Init_Random(void) {
  unsigned i;

  ip=128;
  ip1=ip-24;
  ip2=ip-55;
  ip3=ip-61;

  for (i=ip3; i<ip; i++)
    ira[i] = rand4init();
}

// Box-Muller method
double norm_Parisi(double std_dev, double mean) {
  double   u, r, theta;
  double   x;
  double   norm_rv;

  u = 0.0;
  while (u == 0.0)
    u = FRANDOM;
  r = sqrt(-2.0 * log(u));
  theta = 0.0;
  while (theta == 0.0)
    theta = 2.0 * M_PI * FRANDOM;
  x = r * cos(theta);
  norm_rv = (x * std_dev) + mean;
  return(norm_rv);
}
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/

/* Simple function to account for the periodic boundary condition */
double periodic(double dis) {
	if (dis<-LX/2.)	return dis+LX;
	if (dis>LX/2.) return dis-LX;
}

/* Simple function to account for angle periodicity */
double mod_PI(double anga) {
	if (anga > M_PI) {
		while (anga>M_PI) anga -= 2.*M_PI;
	} else if (anga < -M_PI) {
		while (anga<-M_PI) anga += 2.*M_PI;
	}
	return anga;
}

/* Linked cells */
void fill_list(passive_particle *p, int *head, int *lscl) {
	int mc[2];
	int i,c;

	for (i=0;i<LA_X*LA_Y;i++) head[i] = -1;
	for (i=0;i<N;i++) {
		mc[0] = (int) floor((p[i].x + LX/2.)/(1.*RA_X));
		mc[1] = (int) floor((p[i].y + LY/2.)/(1.*RA_Y));
		c = mc[0]*LA_Y + mc[1];
		lscl[i] = head[c];
		head[c] = i;
	}
}

/* Function to update the positions of particles */
void update(passive_particle *p, passive_particle *ll, passive_particle *lr, int move, int *head, int *lscl) {

	lr[0].x = periodic(lr[1].x - periodic(lr[2].x-lr[1].x)); 
	lr[0].y = - 1.*LY/2.;
	ll[0].x = periodic(ll[1].x - periodic(ll[2].x-ll[1].x)); 
	ll[0].y = - 1.*LX/2.;

	int i,j,k,c,mc[2],c1,mc1[2];

	double r_0 = 2.*R_DISK*pow(2.,-1./6.);
	double r, r_a, fa, len;

	double force_p[N][2];
	double force_a[2];
	double torque[2];

	// passive: noise
	for (i=0;i<N;i++) {
		for (j=0;j<2;j++)
			force_p[i][j] = norm_Parisi(sqrt(2.*D/DT),0.);
	}

	// lever: noise
	force_a[0] = 0.;
	force_a[1] = 0.;
	torque[0] = norm_Parisi(sqrt(2.*D_R/DT),0.);
	torque[1] = norm_Parisi(sqrt(2.*D_R/DT),0.);

	// passive-passive interactions

	// Scan inner cells
	for (mc[0]=0; mc[0]<LA_X; (mc[0])++)
	for (mc[1]=0; mc[1]<LA_Y; (mc[1])++) {
		// Calculate a scalar cell index
		c = mc[0]*LA_Y + mc[1];
		// Scan the neighbor cells (including itself) of cell c
		for (mc1[0]=mc[0]-1; mc1[0]<=mc[0]+1; (mc1[0])++)
		for (mc1[1]=mc[1]-1; mc1[1]<=mc[1]+1; (mc1[1])++) {

			if (mc1[1]>=0&&mc1[1]<LA_Y) {
				// Calculate the scalar cell index of the neighbor cell
				c1 = ((mc1[0]+LA_X)%LA_X)*LA_Y + mc1[1];
				// Scan atom i in cell c
				i = head[c];
				while (i != -1) {
					// Scan atom j in cell c1
					j = head[c1];
					while (j != -1) {
						if (i < j) { // Avoid double counting of pair (i, j)
							r = sqrt(pow(periodic(p[i].x-p[j].x),2.)+pow(p[i].y-p[j].y,2.));
							if (r<2.*R_DISK) {
								fa = 4.*K_WCA*(-12.*pow(r,-13.)/pow(r_0,-12.)+ 6.*pow(r,-7.)/pow(r_0,-6.));
								force_p[i][0] -= fa*periodic(p[i].x-p[j].x)/r;
								force_p[i][1] -= fa*(p[i].y-p[j].y)/r;
								force_p[j][0] += fa*periodic(p[i].x-p[j].x)/r;
								force_p[j][1] += fa*(p[i].y-p[j].y)/r;
							}
						}
						j = lscl[j];
					}
					i = lscl[i];
				}
			}

		}
	}

	// passive-lever interaction
	r_0 = (R_DISK+R_LEVER)*pow(2.,-1./6.);
	double dx, dy;

	for (i=0;i<N;i++) {
		for (j=0;j<N_SUB;j++) {

			// Left lever
			dx = periodic(p[i].x-ll[j].x);
			dy = p[i].y-ll[j].y;
			r = sqrt(dx*dx+dy*dy);

			if (r<(R_DISK+R_LEVER)) {
				len = sqrt(pow(ll[j].x-ll[0].x,2.) + pow(ll[j].y-ll[0].y,2.));

				fa = 4.*(K_WCA/10.)*(-12.*pow(r,-13.)/pow(r_0,-12.)+ 6.*pow(r,-7.)/pow(r_0,-6.));

				force_p[i][0] -= fa*dx/r;
				force_p[i][1] -= fa*dy/r;
				force_a[0] += fa*dx/r;
				torque[0] += len*fa*(dy*cos(M_PI-ang[0])/r-dx*sin(M_PI-ang[0])/r);
			}

			// Right lever
			dx = periodic(p[i].x-lr[j].x);
			dy = p[i].y-lr[j].y;
			r = sqrt(dx*dx+dy*dy);

			if (r<(R_DISK+R_LEVER)) {
				len = sqrt(pow(lr[j].x-lr[0].x,2.) + pow(lr[j].y-lr[0].y,2.));

				fa = 4.*(K_WCA/10.)*(-12.*pow(r,-13.)/pow(r_0,-12.)+ 6.*pow(r,-7.)/pow(r_0,-6.));

				force_p[i][0] -= fa*dx/r;
				force_p[i][1] -= fa*dy/r;
				force_a[1] += fa*dx/r;
				torque[1] += len*fa*(dy*cos(ang[1])/r-dx*sin(ang[1])/r);
			}
		}
	}

	// passive-wall interaction
	r_0 = 2.*R_DISK*pow(2.,-1./6.);
	for (c=0;c<LA_X;c++) {
		i = head[c*LA_Y];

		while (i!=-1) {
			r = p[i].y + 1.*LY/2 + R_DISK;
			if (r < 2.*R_DISK) {
				fa = 4.*K_WCA*(-12.*pow(fabs(r),-13.)/pow(r_0,-12.)+ 6.*pow(fabs(r),-7.)/pow(r_0,-6.));
				force_p[i][1] -= fa*r/fabs(r);
			}
			i = lscl[i];
		}

		i = head[c*LA_Y + LA_Y-1];

		while (i!=-1) {
			r = 1.*LY/2 + R_DISK - p[i].y;
			if (r < 2.*R_DISK) {
				fa = 4.*K_WCA*(-12.*pow(fabs(r),-13.)/pow(r_0,-12.)+ 6.*pow(fabs(r),-7.)/pow(r_0,-6.));
				force_p[i][1] += fa*r/fabs(r);
			}
			i = lscl[i];
		}

	}

	// Updating the hinge
	if (move) {

		dx = MU*DT*F_ACT*(cos(ang[0]) - cos(ang[1]));
		dx += MU*DT*(force_a[0]+force_a[1]);
		ll[0].x = periodic(ll[0].x + dx);
		lr[0].x = periodic(lr[0].x + dx);

		if (ang[0] > ANG_MAX || ang[0] < 0.) torque[0] = 0.;
		if (ang[1] > ANG_MAX || ang[1] < 0.) torque[1] = 0.;
		if (ang[0] + ang[1] >= M_PI - 2.*R_DISK/R) {
			torque[0] = 0.;
			torque[1] = 0.;
		}
		ang[0] -= MU_R*torque[0]*DT;
		ang[1] += MU_R*torque[1]*DT;
		if (ang[0] > ANG_MAX) ang[0] = ANG_MAX;
		if (ang[0] < 0.) ang[0] = 0.;
		if (ang[1] > ANG_MAX) ang[1] = ANG_MAX;
		if (ang[1] < 0.) ang[1] = 0.;

		for (i=1;i<N_SUB;i++) {
			r = 1.*i*R_LEVER/2.;
	    ll[i].x = -1.*r*cos(ang[0]) + ll[0].x;
	    ll[i].x = periodic(ll[i].x);
	    ll[i].y = r*sin(ang[0]) + ll[0].y;
	    lr[i].x = r*cos(ang[1]) + lr[0].x;
	    lr[i].x = periodic(lr[i].x);
	    lr[i].y = r*sin(ang[1]) + lr[0].y;
		}
	}

	// position passive
	for (i=0;i<N;i++) {
		p[i].x += MU*force_p[i][0]*DT;
		p[i].x = periodic(p[i].x);
		p[i].y += MU*force_p[i][1]*DT;
	}

	if (USE_ACT_FRAME) {
		double disp = ll[0].x;

		// passive
		for (i=0;i<N;i++)
			p[i].x = periodic(p[i].x - disp);

		// hinge
		for (i=0;i<N_SUB;i++) {
			ll[i].x = periodic(ll[i].x - disp);
			lr[i].x = periodic(lr[i].x - disp);
		}
    		
	}

}

/* [Initial setup] Finding possible positions for the passive particles */
void find_possible_positions(passive_particle *pos, int *avaiPos, passive_particle *ll, passive_particle *lr) {
	int ix, iy, id;
	int nLines = (int) floor((LY-2.*R_DISK)/(2.*R_DISK*PADDING)) + 1;
	double dy = 1.*(LY-2.*R_DISK)/(nLines-1.);
	int nColumns = (int) floor(LX/(2.*R_DISK*PADDING));
	double dx = 1.*LX/nColumns;

	int nPos = nLines*nColumns;

	// All possible positions on the grid
	for (iy=0;iy<nLines;iy++) {
		for (ix=0;ix<nColumns;ix++) {
			pos[ix+iy*nColumns].x = periodic(ix*dx + (iy%2)*dx/2.);
			pos[ix+iy*nColumns].y = iy*dy - LY/2. + R_DISK;
			avaiPos[ix+iy*nColumns] = 1;
		}	
	}

	// Removing crosses with the initial position of the hinge
	int i;
	double dist;
	for (i=0;i<N_SUB;i++) {
		for (id=0;id<nPos;id++) {
			if (avaiPos[id]) {
				dist = sqrt( pow(pos[id].x - ll[i].x,2) + pow(pos[id].y - ll[i].y,2) );
				if (dist < R_DISK+R_LEVER)
					avaiPos[id] = 0;
				dist = sqrt( pow(pos[id].x - lr[i].x,2) + pow(pos[id].y - lr[i].y,2) );
				if (dist < R_DISK+R_LEVER)
					avaiPos[id] = 0;
			}
		}
	}

}

/* Find number of passive particles below the hinges */
void get_N(int *Ns, passive_particle *p, passive_particle *ll, passive_particle *lr) {
	int i;
  Ns[0] = 0; Ns[1] = 0;
  double r, th;
  for (i=0;i<N;i++) {
  	r = sqrt( pow(periodic(ll[0].x - p[i].x),2) + pow(ll[0].y - p[i].y,2) );
  	th = M_PI - atan2(p[i].y - ll[0].y, periodic(p[i].x - ll[0].x));
  	if (r < R && th < ang[0]) Ns[0]++;
  	r = sqrt( pow(periodic(lr[0].x - p[i].x),2) + pow(lr[0].y - p[i].y,2) );
  	th = atan2(p[i].y - lr[0].y, periodic(p[i].x - lr[0].x));
  	if (r < R && th < ang[1]) Ns[1]++;
  }
}

/* Function to initialize the system */
void init(passive_particle *p, passive_particle *ll, passive_particle *lr, int *head, int *lscl){

  int i,j;

	ang[0] = ((ANG_MAX > (M_PI/2.-1.*R_DISK/R)) ? (M_PI/2.-1.*R_DISK/R) : ANG_MAX);
	ang[1] = ((ANG_MAX > (M_PI/2.-1.*R_DISK/R)) ? (M_PI/2.-1.*R_DISK/R) : ANG_MAX);

  // Hinge
  double r;
  for (i=0;i<N_SUB;i++) {
  	r = 1.*i*R_LEVER/2.;
  	ll[i].x = -1.*r*cos(ang[0]);
  	ll[i].y = r*sin(ang[0]) - LY/2.;
  	lr[i].x = r*cos(ang[1]);
  	lr[i].y = r*sin(ang[1]) - LY/2.;
  }

	// Passive particles
	int nLines = (int) floor((LY-2.*R_DISK)/(2.*R_DISK*PADDING)) + 1;
	int nColumns = (int) floor(LX/(2.*R_DISK*PADDING));
	int Npos = nLines*nColumns;

	passive_particle initPos[Npos];
	int avaiPos[Npos];
	find_possible_positions(initPos,avaiPos,ll,lr);

	int found, id;
	for (i=0;i<N;i++) {
		found = 0;
		while (!found) {
			id = (int) floor(FRANDOM*Npos);
			if (avaiPos[id])
				found = 1;
		}
		p[i].x = initPos[id].x;
    p[i].y = initPos[id].y;
    avaiPos[id] = 0;
	}

	for (i=0;i<T_EQ;i++) {
		fill_list(p, head, lscl);
		update(p,ll,lr,0,head,lscl);
	}

}

/* Save particle's positions to a file */
void save_configuration(passive_particle *p, passive_particle *ll, passive_particle *lr) {
    int n;

    FILE *outFile;
    char name[100];
    sprintf(name,FILE_NAME_CONF,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,T_MAX/N_SKIP-1,N,N_SUB);
    outFile = fopen(name, "a");

    for (n=0;n<N_SUB;n++)
        fprintf(outFile,"%f %f\n",ll[n].x,ll[n].y); // Left rod

    for (n=0;n<N_SUB;n++)
      fprintf(outFile,"%f %f\n",lr[n].x,lr[n].y); // Right rod

    for (n=0;n<N;n++)
        fprintf(outFile,"%f %f\n",p[n].x,p[n].y); // Passive particles

    fclose(outFile);
}

/* Main function */
void main() {

	// Random numbers initialization
	srand(time(0));
	myrand = SEED;
	Init_Random();

	// Particles initialization
	passive_particle p[N]; // N passive particles
	passive_particle ll[N_SUB]; // 1 hinge structure (left)
	passive_particle lr[N_SUB]; // 1 hinge structure (right)
	passive_particle aux[N_SUB]; // Auxiliary

	// Cell-List
	int head[LA_X*LA_Y];
	int lscl[N];

	int t, t_tot, sim, i;

	// Clear save files
	FILE *outFile;
  char name[100];

	if (SAVE_DETAILS) {
    sprintf(name,FILE_NAME_ANG,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
    outFile = fopen(name, "w");
    fclose(outFile);
	}
	if (SAVE_CONFIGS) {
    sprintf(name,FILE_NAME_CONF,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D,T_MAX/N_SKIP-1,N,N_SUB);
    outFile = fopen(name, "w");
    fclose(outFile);
	}

  // Save ang
  sprintf(name,FILE_NAME_ANG,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
  outFile = fopen(name, "w");

  // Save mea
  FILE *outFileMea;
  sprintf(name,FILE_NAME_MEA,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
  outFileMea = fopen(name, "w");

  // Simulations
	double dist, dx, dx_tot, x0, xant;
	double tabs[2], pos[2], spd[2];
	int nabs, Ns[2];

	sim = 0;
	do {
		init(p,ll,lr,head,lscl);

    dx = 0.;
    dx_tot = 0.;
		x0 = ll[0].x;
		xant = x0;
		for (t=0;t<2;t++) {
			tabs[t] = 0.;
			pos[t]  = 0.;
			spd[t]  = 0.;
		}
		nabs = 0;
		t = 1;
		t_tot = 1;

		do {

			for (i=0;i<N_SUB;i++) aux[i] = ll[i];
			for (i=0;i<N_SUB;i++) aux[i] = lr[i];

			fill_list(p,head,lscl);
    	update(p,ll,lr,1,head,lscl);

    	// Save configurations if needed
			if (t%N_SKIP==0 && SAVE_CONFIGS) save_configuration(p,ll,lr);
			// Save extra details if needed (angles of active particles and hinge displacement after every N_SKIP step)
			if (t%N_SKIP==0 && SAVE_DETAILS) { //
				get_N(Ns,p,ll,lr);
				fprintf(outFile,"%d %d %d %d %f %d %f %d %f %f\n",sim,nabs,t,t_tot,ang[0],Ns[0],ang[1],Ns[1],dx,dx_tot);
			}

			dx += periodic(ll[0].x - xant); // Displacement since last rod was absorbed
			dx_tot += periodic(ll[0].x - xant); // Total hinge's displacement
			xant = ll[0].x;

			// If one of the rods is absorbed
			if ((ang[0] == 0. || ang[1] == 0.) && nabs == 0) {
				tabs[0] = 1.*t*DT; // Save FPT
				pos[0] = dx; // Save total hinge displacement
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

		} while (t<T_MAX && nabs<2); // Run either until T_MAX or until both rods get absorbed

		/* Saving basic measurements into file */
		/* Saving: 	simulation id, FPT of FIRST rod, displacement of the hinge until FIRST rod was absorbed, average speed of hinge until FIRST rod was absorbed,
								FPT of SECOND rod, displacement of the hinge until SECOND rod was absorbed, average speed of hinge until SECOND rod was absorbed, 
								total simulation time, total displacement of the hinge, average hinge speed, number of absorbed rods
		*/
		if (nabs==0) {
			// If neither rod gets absorbed
			fprintf(outFileMea,"%d %f %f %f %f %f %f ",sim,0.,0.,0.,0.,0.,0.);
		} else if (nabs==1) {
			// If one rod gets absorbed
			fprintf(outFileMea,"%d %f %f %f %f %f %f ",sim,tabs[0],pos[0],spd[0],0.,0.,0.);
		} else {
			// If both rods get absorbed
			fprintf(outFileMea,"%d %f %f %f %f %f %f ",sim,tabs[0],pos[0],spd[0],tabs[1],pos[1],spd[1]);
		}
		fprintf(outFileMea,"%f %f %f %d\n",1.*t_tot*DT,dx_tot,1.*dx_tot/(1.*t_tot*DT),nabs);

		sim++;

	} while (sim < N_SIM);

	fclose(outFileMea);
	fclose(outFile);

}
