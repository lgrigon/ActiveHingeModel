/*
double_hinge.c

This algorithm simulates the active hinge model.
This version is similar to the double_hinge.c, but simulates only one rod.
(assuming the first one was alreadyu absorbed).
Outputs data from the simulations into files at 'data/' folder.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

/** SIMULATION PARAMETERS **/

// System size
#define LX          10.
#define LY          10.

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

/** FILE SAVE **/
#define SAVE_CONFIGS		1 // Do we want to also save the configurations?
#define SAVE_DETAILS		1 // Do we want to save detailed information?
#define SAVE_MSD				0 // Do we want to save the MSD?
#define SAVE_FINALS 		0 // Do we want to save special information about the last
													// 50 steps of the simulation? (such as torque, config., etc)
#define SAVE_TORQUE  		0 // Do we want to save torque information?
#define SAVE_N0					1 // Do we want to save informations about the system when
													// the number of particles below the hinge is equal to 0 or 1?

#define FILE_NAME_CONF 			"data/conf_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"
#define FILE_NAME_ANG  			"data/ang_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"
#define FILE_NAME_MEA 			"data/mea_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"
#define FILE_NAME_FOR 			"data/forces_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"
#define FILE_NAME_MSD 			"data/msd_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"
#define FILE_NAME_FIN 			"data/finals_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"
#define FILE_NAME_FIN_CONF 	"data/finalsConf_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"
#define FILE_NAME_CONF_N0 	"data/confN0_L%dx%d-R%.1f-pf%.2f-f%.1f-ang%.3f-D%d.dat"

/** FIXED PARAMETERS **/

// System
#define DT				0.00005

#define MU				0.05 // Motility coefficient
#define D					200. // Diffusion coefficient

#define K_WCA				50 // WCA potential strenght

// Active particles
#define R_LEVER		0.1 // Radius of active rod's particles
#define N_SUB			(int) round(2*R/R_LEVER) // Number of those particles
#define MU_R			0.01 // Rotational mobility coefficient 
#define D_R				0. // Rotational diffusion coefficient

// Passive
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

// Some global variables
double ang[2]; // Global variable for the angle of each rod
double torque_rod[3]; // Global variable for the torque of each rod (+total torque)
int coll[2]; // Global variable for the number of collisions with each rod

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
	ll[0].y = - 1.*LY/2.;

	int i,j,k,c,mc[2],c1,mc1[2];

	double r_0 = 2.*R_DISK*pow(2.,-1./6.);
	double r, r_a, fa, len;

	double force_p[N][2];
	double force_a[2];
	double torque[2], trqa;

	torque_rod[0] = 0.;
	torque_rod[1] = 0.;
	coll[0] = 0;
	coll[1] = 0;

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

	// passive-lever interactions
	r_0 = (R_DISK+R_LEVER)*pow(2.,-1./6.);
	double dx, dy;

	for (i=0;i<N;i++) {
		for (j=0;j<N_SUB;j++) {

			// Left rod
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

			// Right rod
			dx = periodic(p[i].x-lr[j].x);
			dy = p[i].y-lr[j].y;
			r = sqrt(dx*dx+dy*dy);

			if (r<(R_DISK+R_LEVER)) {
				len = sqrt(pow(lr[j].x-lr[0].x,2.) + pow(lr[j].y-lr[0].y,2.));

				fa = 4.*(K_WCA/10.)*(-12.*pow(r,-13.)/pow(r_0,-12.)+ 6.*pow(r,-7.)/pow(r_0,-6.));

				force_p[i][0] -= fa*dx/r;
				force_p[i][1] -= fa*dy/r;
				force_a[1] += fa*dx/r;
				trqa = len*fa*(dy*cos(ang[1])/r-dx*sin(ang[1])/r);
				torque[1] += trqa;
				torque_rod[(trqa>0.?0:1)] += trqa;
				coll[(trqa>0.?0:1)]++;
			}
		}
	}

	// passive-wall interactions
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
		torque_rod[2] = torque[1];
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

/* Function to get local packing fraction underneath the rod */
void get_PF(passive_particle *p, passive_particle *lr, double *pfrod) {
  int i;
  double r, th;
  pfrod[0] = 0.;
  pfrod[1] = 0.;
  for (i=0;i<N;i++) {
  	r = sqrt( pow(periodic(lr[0].x - p[i].x),2) + pow(lr[0].y - p[i].y,2) );
  	th = atan2(p[i].y - lr[0].y, periodic(p[i].x - lr[0].x));

  	if (r < R + R_DISK) {

  		if (r <= R - R_DISK)
  			pfrod[th<ang[1]?0:1] += 1.*M_PI*R_DISK*R_DISK;
  		else {
  			pfrod[th<ang[1]?0:1] += 1.*R_DISK*R_DISK*acos(1.*(r*r+R_DISK*R_DISK-R*R)/(2.*r*R_DISK));
  			pfrod[th<ang[1]?0:1] += 1.*R*R*acos(1.*(r*r+R*R-R_DISK*R_DISK)/(2.*r*R));
  			pfrod[th<ang[1]?0:1] -= 0.5*sqrt((-r+R_DISK+R)*(r+R_DISK-R)*(r-R_DISK+R)*(r+R_DISK+R));
  		}

  	}
  	
  }
  if (ang[1]>0.) pfrod[0] /= 0.5*R*R*ang[1];
  if ((M_PI-ang[1])>0.) pfrod[1] /= 0.5*R*R*(M_PI-ang[1]);

}

/* Function to get the position of one particle when it is below the hinge */
int get_N1_pos(passive_particle *p, passive_particle *ll, passive_particle *lr) {
	int i;
  	double r, th;
  	for (i=0;i<N;i++) {
			r = sqrt( pow(periodic(lr[0].x - p[i].x),2) + pow(lr[0].y - p[i].y,2) );
			th = atan2(p[i].y - lr[0].y, periodic(p[i].x - lr[0].x));
			if (r < R && th < ang[1]) return i;
  	}
  	return 0;
}

/* Function to initialize the system */
void init(passive_particle *p, passive_particle *ll, passive_particle *lr, int *head, int *lscl){

  int i,j;

	ang[0] = 0.0;
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

/* store positions for msd */
void add_pos_to_list(int t, passive_particle *p, double *posX, double *posY) {
    int i;

    int nP = (int) N;
    int nMea = (int) T_MAX/N_SKIP;

    for (i=0;i<nP;i++) {
        posX[t*nP + i] = p[i].x;
        posY[t*nP + i] = p[i].y;
    }
}

/* Calculates the MSD */
void calc_MSD(double *posX, double *posY, double *MSD, int *count, int tm) {

    int t0, t, i;

    tm = (int)floor(tm/N_SKIP);

    double msdA[tm];
    int countA[tm];

    for (t=0;t<tm;t++) {
        msdA[t] = 0.;
        countA[t] = 0;
    }

    int nP = (int) N;

    double dx, dy;

    for (t0=0;t0<tm;t0++) {
        for (t=t0;t<tm;t++) {
            for (i=0;i<nP;i++) {
                dx = periodic(posX[t*nP+i] - posX[t0*nP+i]);
                dy = posY[t*nP+i] - posY[t0*nP+i];
                msdA[t-t0] += dx*dx+dy*dy;
                countA[t-t0]++;
            }
        }
    }

    for (t=0;t<tm;t++) {
    	if (countA[t]>0) msdA[t] /= 1.*countA[t];
    	MSD[t] += msdA[t];
    	count[t]++;
    }
}

/* To keep track of last 50 measurements to study final absorption */
void update_last(double *trq, double *angles, int *nbelow, int n) {
	int t;
	for (t=0;t<49;t++) {
		trq[2*t]   = trq[2*(t+1)];
		trq[2*t+1] = trq[2*(t+1)+1];
		angles[t] = angles[t+1];
		nbelow[t] = nbelow[t+1];
	}
	trq[2*49] 	= torque_rod[0];
	trq[2*49+1] = -1.*torque_rod[1];
	angles[49]	= ang[1];
	nbelow[49]	= n;
}

/* Keeping track of the particles closer to the hinge to study the dynamics */
int get_N7(int *ids, passive_particle *p, passive_particle *lr) {
	int i;
	int Ns = 0;
	double r, th;
	double rBelow[30];
	for (i=0;i<N;i++) {
		r = sqrt( pow(periodic(lr[0].x - p[i].x),2) + pow(lr[0].y - p[i].y,2) );
		th = atan2(p[i].y - lr[0].y, periodic(p[i].x - lr[0].x));
		if (r < R*2.5) {
			if (th < ang[1] || (r>0.9*R && th < ang[1]+0.25)) {
				ids[Ns] = i;
				rBelow[Ns] = r;
				Ns++;
			}
		} 
	}
	int j, id;
	double rmin;
	for (i=0;i<Ns;i++) {
		rmin = rBelow[i];
		for (j=i+1;j<Ns;j++) {
			if (rmin>rBelow[j]){
				id = j;
				rmin = rBelow[j];
			}
		}
		if (rBelow[i]>rmin) {
			rBelow[id] = rBelow[i];
			rBelow[i] = rmin;
			j = ids[id];
			ids[id] = ids[i];
			ids[i] = j;
		}
	}

	return Ns;
}

/* To keep track of last 50 particle configurations to study final absorption */
void update_pos_last(passive_particle *p, passive_particle *lr, double *posX, double *posY) {
	//Configurations
	int i, t;
	for (t=0;t<49;t++) {
		for (i=0;i<7;i++) {
			posX[t*7 + i] = posX[(t+1)*7 + i];
			posY[t*7 + i] = posY[(t+1)*7 + i];
		}
	}
	int idBelow[30];
	int Ns = get_N7(idBelow, p, lr);
	if (Ns>7) Ns = 7;
	for (i=0;i<Ns;i++) {
    	posX[49*7 + i] = p[idBelow[i]].x;
    	posY[49*7 + i] = p[idBelow[i]].y;
    }
}

/* Save final 50 measurements to a file */
void save_final(int sim, double *trq, double *angles, int *nbelow, double *posX, double *posY) {
	FILE *outFile;
	char name[100];

	sprintf(name,FILE_NAME_FIN,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
	outFile = fopen(name, "a");
	int t;
	for (t=0;t<50;t++)
		fprintf(outFile,"%d %d %f %f %f %d\n",sim,-49+t,trq[2*t],trq[2*t+1],angles[t],nbelow[t]);
	fclose(outFile);

	sprintf(name,FILE_NAME_FIN_CONF,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
	outFile = fopen(name, "a");
	int i;
	for (t=0;t<50;t++)
		for (i=0;i<7;i++)
		fprintf(outFile,"%d %d %f %f\n",sim,-49+t,posX[t*7+i],posY[t*7+i]);
	fclose(outFile);
}

/* Save particle's positions to a file */
void save_configuration(passive_particle *p, passive_particle *ll, passive_particle *lr) {
    int n;

    FILE *outFile;
    char name[100];
    sprintf(name,FILE_NAME_CONF,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
    outFile = fopen(name, "a");

    for (n=0;n<N_SUB;n++)
        fprintf(outFile,"%f %f %f\n",ll[n].x,ll[n].y,R_LEVER);

    for (n=0;n<N_SUB;n++)
      fprintf(outFile,"%f %f %f\n",lr[n].x,lr[n].y,R_LEVER);

    for (n=0;n<N;n++)
      	fprintf(outFile,"%f %f %f\n",p[n].x,p[n].y,R_DISK);

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

	// Clear save file
	FILE *outFile;
  char name[100];

	if (SAVE_CONFIGS) {
	    sprintf(name,FILE_NAME_CONF,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
	    outFile = fopen(name, "w");
	    fclose(outFile);
	}

	if (SAVE_FINALS) {
	    sprintf(name,FILE_NAME_FIN,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
	    outFile = fopen(name, "w");
	    fclose(outFile);
	    sprintf(name,FILE_NAME_FIN_CONF,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
	    outFile = fopen(name, "w");
	    fclose(outFile);
	}

	if (SAVE_DETAILS) {
	  sprintf(name,FILE_NAME_ANG,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
	  outFile = fopen(name, "w");
	}
  	
	FILE *outFileMea;
	sprintf(name,FILE_NAME_MEA,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
	outFileMea = fopen(name, "w");

	FILE *outFileTor;
	if (SAVE_TORQUE) {
		sprintf(name,FILE_NAME_FOR,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
		outFileTor = fopen(name, "w");
	}

	FILE *outFileMSD;
	if (SAVE_MSD) {
		sprintf(name,FILE_NAME_MSD,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
		outFileMSD = fopen(name, "w");
	}

	FILE *outFileN0;
	if (SAVE_N0) {
		sprintf(name,FILE_NAME_CONF_N0,(int)LX,(int)LY,R,PF,F_ACT,ANG_MAX,(int)D);
		outFileN0 = fopen(name, "w");
	}

	int nMea = (int) T_MAX/N_SKIP; // number of measurements
	// variables for the MSD
	double *posX = (double *) malloc(nMea*((int)N)*sizeof(double));
	double *posY = (double *) malloc(nMea*((int)N)*sizeof(double));
	double *MSD = (double *) malloc(nMea*sizeof(double));
	// to keep track of some (7) particles closest to the rod
	double posX_fin[50*7];
	double posY_fin[50*7];

	int cnt_msd[nMea];
	for (t=0;t<nMea;t++) {
		MSD[t] = 0.;
		cnt_msd[t] = 0;
	}

  // Simulations
	double dist, dx, dx_tot, x0, xant;
	double tabs[2], pos[2], spd[2], pfrod[2];
	int nabs, Ns[2], id_N1;
	double tqr_final[50*2], ang_final[50]; 
	int nbelow[50];
	int Nt;

	sim = 0;
	do {
		init(p,ll,lr,head,lscl);
		add_pos_to_list((int)t/N_SKIP,p,posX,posY);

		for (t=0;t<50;t++) {
			tqr_final[2*t] = 0.;
			tqr_final[2*t+1] = 0.;
			ang_final[t] = 0.;
			nbelow[t] = 0;
			for (i=0;i<7;i++) {
				posX_fin[7*t+i] = 0.;
				posY_fin[7*t+i] = 0.;
			}
		}

    dx = 0.;
    dx_tot = 0.;
		x0 = ll[0].x;
		xant = x0;
		for (t=0;t<2;t++) {
			tabs[t] = 0.;
			pos[t]  = 0.;
			spd[t]  = 0.;
		}
		nabs = 1;
		t = 1;
		t_tot = 1;

		do {

			for (i=0;i<N_SUB;i++) aux[i] = ll[i];
			for (i=0;i<N_SUB;i++) aux[i] = lr[i];

			fill_list(p,head,lscl);
    	update(p,ll,lr,1,head,lscl);

    		if (t%N_SKIP==0) {
    			// Save configurations if needed
    			if (SAVE_CONFIGS) save_configuration(p,ll,lr);
    			// Save extra details if needed (angles of active particles and hinge displacement after every N_SKIP step)
					if (SAVE_DETAILS) {
						get_N(Ns,p,ll,lr);
						fprintf(outFile,"%d %d %d %d %f %d %f %d %f %f\n",sim,nabs,t,t_tot,ang[0],Ns[0],ang[1],Ns[1],dx,dx_tot);
					}
					// Save torque information if needed
					if (SAVE_TORQUE) fprintf(outFileTor,"%d %f %f %f\n",sim,torque_rod[0],torque_rod[1],torque_rod[2]);
					// Save MSD if needed
					if (SAVE_MSD) add_pos_to_list((int)t/N_SKIP,p,posX,posY);\
					// Save extra information about the last 50 measurements of simulation if needed
					if (SAVE_FINALS) {
						get_N(Ns,p,ll,lr);
						update_last(tqr_final,ang_final,nbelow,Ns[1]);
						update_pos_last(p,lr,posX_fin,posY_fin);
					}
					// Save additional information about configurations with 0/1 particles below the hinge
					if (SAVE_N0) {
						get_N(Ns,p,ll,lr);
						if (Ns[1]==0) {
							//saveVisual_Configurations(p,ll,lr);
							get_PF(p,lr,pfrod);
							fprintf(outFileN0,"%d %d %d %f %f %f %f %f %f %d %d %f %f %f\n",
								sim,t,0,ang[1],p[id_N1].x,p[id_N1].y,torque_rod[0],torque_rod[1],torque_rod[2],coll[0],coll[1],
								pfrod[0],pfrod[1],periodic(ll[0].x - xant));
						} else if (Ns[1]==1) {
							id_N1 = get_N1_pos(p,ll,lr);
							get_PF(p,lr,pfrod);
							fprintf(outFileN0,"%d %d %d %f %f %f %f %f %f %d %d %f %f %f\n",
								sim,t,1,ang[1],p[id_N1].x,p[id_N1].y,torque_rod[0],torque_rod[1],torque_rod[2],coll[0],coll[1],
								pfrod[0],pfrod[1],periodic(ll[0].x - xant));
						}
					}
    		}

			dx += periodic(ll[0].x - xant);
			dx_tot += periodic(ll[0].x - xant);
			xant = ll[0].x;

			// If the rod is absorbed
			if (ang[0] == 0. && ang[1] == 0.) {
				tabs[1] = 1.*t*DT; // Save FPT
				pos[1] = dx; // Save total hinge displacement
				spd[1] = 1.*pos[1]/tabs[1]; // Save hinge's average speed
				nabs++;
			}

			t++;
			t_tot++;

		} while (t<T_MAX && nabs<2);

		if (SAVE_N0) {
			get_N(Ns,p,ll,lr);
			if (Ns[1]==0) save_configuration(p,ll,lr);
		}

		if (nabs==2 && SAVE_FINALS) {
			get_N(Ns,p,ll,lr);
			update_last(tqr_final,ang_final,nbelow,Ns[1]);
			update_pos_last(p,lr,posX_fin,posY_fin);
			save_final(sim,tqr_final,ang_final,nbelow,posX_fin,posY_fin);
		}

		if (nabs==0) {
			fprintf(outFileMea,"%d %f %f %f %f %f %f ",sim,0.,0.,0.,0.,0.,0.);
		} else if (nabs==1) {
			fprintf(outFileMea,"%d %f %f %f %f %f %f ",sim,tabs[0],pos[0],spd[0],0.,0.,0.);
		} else {
			fprintf(outFileMea,"%d %f %f %f %f %f %f ",sim,tabs[0],pos[0],spd[0],tabs[1],pos[1],spd[1]);
		}
		fprintf(outFileMea,"%f %f %f %d\n",1.*t_tot*DT,dx_tot,1.*dx_tot/(1.*t_tot*DT),nabs);

		if (SAVE_MSD) calc_MSD(posX,posY,MSD,cnt_msd,t);

		sim++;

	} while (sim < N_SIM);

	if (SAVE_TORQUE) fclose(outFileTor);

	if (SAVE_CONFIGS) save_configuration(p,ll,lr);
	if (SAVE_CONFIGS) printf("%f (%d)\n",tabs[1],(int)(tabs[1]/(DT*N_SKIP)));

	if (SAVE_N0) fclose(outFileN0);

	if (SAVE_MSD) {
		for (t=0;t<nMea;t++) {
			if (cnt_msd[t]>0) MSD[t] /= 1.*cnt_msd[t];
			fprintf(outFileMSD,"%f %f\n",1.*t*N_SKIP*DT,MSD[t]);
		}
		fclose(outFileMSD);
	}

	fclose(outFileMea);
	if (SAVE_DETAILS) fclose(outFile);

	free(posX); free(posY); free(MSD);

}
