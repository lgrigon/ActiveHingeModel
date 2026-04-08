#include "particles.h"
#include "params.h"
#include "random.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>

double ang[2];

// Function to update the positions of particles
void update(passive_particle *p, passive_particle *ll, passive_particle *lr, int move, int *head, int *lscl) {

	int i,j,k,c,mc[2],c1,mc1[2];

	double r_0 = 2.*R_DISK*pow(2.,-1./6.);
	double r, r_a, fa, len;

	double force_p[N_PASSIVE][2];
	double force_a[2];
	double torque[2];

	// passive: noise
	for (i=0;i<N_PASSIVE;i++) {
		for (j=0;j<2;j++)
			force_p[i][j] = rand_norm(sqrt(2.*KBT/(MU_PASSIVE*DT)),0.);
	}

	// lever: noise
	force_a[0] = 0.;
	force_a[1] = 0.;
	torque[0] = 0.;
	torque[1] = 0.;

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
								if (move)
									fa = 4.*K_WCA*(-12.*pow(r,-13.)/pow(r_0,-12.)+ 6.*pow(r,-7.)/pow(r_0,-6.));
								else
									fa = 0.1*K_WCA*(-12.*pow(r,-13.)/pow(r_0,-12.)+ 6.*pow(r,-7.)/pow(r_0,-6.));
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
	double r_min = 1e-4 * (R_DISK + R_LEVER);
	for (i=0;i<N_PASSIVE;i++) {
		for (j=0;j<N_BEADS;j++) {

			// Left rod
			dx = periodic(p[i].x-ll[j].x);
			dy = p[i].y-ll[j].y;
			r = sqrt(dx*dx+dy*dy);

			if (r < r_min) r = r_min;  // clamp to avoid INF

			if (r<(R_DISK+R_LEVER)) {
				len = sqrt(pow(periodic(ll[j].x-ll[0].x),2.) + pow(ll[j].y-ll[0].y,2.));

				fa = 4.*(K_WCA/10.)*(-12.*pow(r,-13.)/pow(r_0,-12.)+ 6.*pow(r,-7.)/pow(r_0,-6.));

				force_p[i][0] -= fa*dx/r;
				force_p[i][1] -= fa*dy/r;
				force_a[0] += fa*dx/r;
				//torque[0] += len*fa*(dy*cos(PI-ang[0])/r-dx*sin(PI-ang[0])/r);
				torque[0] += len*fa*(dy*cos(ang[0])/r+dx*sin(ang[0])/r);
			}

			// Right rod
			dx = periodic(p[i].x-lr[j].x);
			dy = p[i].y-lr[j].y;
			r = sqrt(dx*dx+dy*dy);

			if (r<r_min) r = r_min; // clamp to avoid INF

			if (r<(R_DISK+R_LEVER)) {
				len = sqrt(pow(periodic(lr[j].x-lr[0].x),2.) + pow(lr[j].y-lr[0].y,2.));

				fa = 4.*(K_WCA/10.)*(-12.*pow(r,-13.)/pow(r_0,-12.)+ 6.*pow(r,-7.)/pow(r_0,-6.));

				force_p[i][0] -= fa*dx/r;
				force_p[i][1] -= fa*dy/r;
				force_a[1] += fa*dx/r;
				torque[1] += len*fa*(dy*cos(ang[1])/r-dx*sin(ang[1])/r);
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
	double mu_rod, mur_rod;
	double xdotfree, tdotfree[2];
	double a_cnst, b_cnst;
	int constraint[2],constraint_0[2],constraint_max[2];
	double xdot, tdot[2];
	double Fx_total;
	if (move) {

		// Actually, this is the rod's drag coefficient
		mu_rod = (2.*N_BEADS+1. -1.5*N_BEADS*(N_BEADS+1.)*(pow(sin(ang[0]),2)+pow(sin(ang[1]),2))/(2.*N_BEADS+1.))/MU_BEADS;
		// Now getting mu_rod back to mobility
		mu_rod = 1./mu_rod;
		
		//force_a[0] += rand_norm(sqrt(2.*KBT/(mu_rod*DT)),0.);
		//force_a[1] += rand_norm(sqrt(2.*KBT/(mu_rod*DT)),0.);
		//torque[0] += rand_norm(sqrt(2.*KBT/(3.*mur_rod*DT)),0.);
		//torque[1] += rand_norm(sqrt(2.*KBT/(3.*mur_rod*DT)),0.);

		// Actually, this is the rod's drag coefficient
		mur_rod = 1.*L_ROD*L_ROD*(1.*N_BEADS+1.)*(2.*N_BEADS+1.)/(6.*N_BEADS*MU_BEADS);
		mur_rod = mur_rod / RATIO_MU;
		mur_rod = 1./mur_rod; // Now getting mur_rod
		//mur_rod = 0.01;

		// Constants for the equations of motion
		a_cnst = 3.*N_BEADS/(L_ROD*(2.*N_BEADS+1.));
		b_cnst = 1.*L_ROD*(1.*N_BEADS+1.)/(2.*MU_BEADS);

		Fx_total = force_a[0] + force_a[1] + F_ACT*(cos(ang[0]) - cos(ang[1]));

		// Free velocities (without constraints)
		xdotfree = Fx_total; // All external forces (x-component)
		xdotfree -= a_cnst*(torque[0]*sin(ang[0]) - torque[1]*sin(ang[1])); // correction due to torques
		xdotfree *= mu_rod;
		tdotfree[0] = mur_rod * (torque[0] - b_cnst*sin(ang[0])*xdotfree);
		tdotfree[1] = mur_rod * (torque[1] + b_cnst*sin(ang[1])*xdotfree);

		// Check for constraints
		constraint_0[0] = (ang[0]+tdotfree[0]*DT <= 0.);
		constraint_max[0] = (ang[0]+tdotfree[0]*DT >= ANG_MAX);
		constraint_0[1] = (ang[1]+tdotfree[1]*DT <= 0.);
		constraint_max[1] = (ang[1]+tdotfree[1]*DT >= ANG_MAX);
		constraint[0] = (constraint_0[0] || constraint_max[0]);
		constraint[1] = (constraint_0[1] || constraint_max[1]);

		// Apply constraints
		if (constraint[0] && constraint[1]) { // both constrained
			xdot = Fx_total;
			xdot /= (1./mu_rod + a_cnst*b_cnst*(pow(sin(ang[0]),2)+pow(sin(ang[1]),2)));
			//xdot *= mu_rod; // extra
			tdot[0] = 0.;
			tdot[1] = 0.;
		} else if (constraint[0] && !constraint[1]) { // only left rod is constrained
			xdot = Fx_total + a_cnst*torque[1]*sin(ang[1]);
			xdot /= (1./mu_rod + a_cnst*b_cnst*pow(sin(ang[0]),2));
			//xdot *= mu_rod; // extra
			tdot[0] = 0.;
			tdot[1] = mur_rod*(torque[1] + b_cnst*sin(ang[1])*xdot);	
		} else if (!constraint[0] && constraint[1]) { // only right rod constrained
			xdot = Fx_total - a_cnst*torque[0]*sin(ang[0]);
			xdot /= (1./mu_rod + a_cnst*b_cnst*pow(sin(ang[1]),2));
			//xdot *= mu_rod; // extra
			tdot[0] = mur_rod*(torque[0] - b_cnst*sin(ang[0])*xdot);
			tdot[1] = 0.;
		} else { // else no constraints
			xdot = xdotfree;
			tdot[0] = tdotfree[0];
			tdot[1] = tdotfree[1];
		}

		//xdot = xdotfree;
		//tdot[0] = tdotfree[0];
		//tdot[1] = tdotfree[1];

		dx = xdot*DT;
		ll[0].x = periodic(ll[0].x + dx);
		lr[0].x = periodic(lr[0].x + dx);
		ang[0] += tdot[0]*DT;
		ang[1] += tdot[1]*DT;

		if (constraint[0]) ang[0] = (constraint_0[0] ? 0.0 : ANG_MAX);
		if (constraint[1]) ang[1] = (constraint_0[1] ? 0.0 : ANG_MAX);

		/*if (ang[0] > ANG_MAX) ang[0] = ANG_MAX;
		if (ang[0] < 0.) ang[0] = 0.;
		if (ang[1] > ANG_MAX) ang[1] = ANG_MAX;
		if (ang[1] < 0.) ang[1] = 0.;*/
		
		for (i=1;i<N_BEADS;i++) {
			r = 1.*i*R_LEVER;
			ll[i].x = -1.*r*cos(ang[0]) + ll[0].x;
			ll[i].x = periodic(ll[i].x);
			ll[i].y = r*sin(ang[0]) + ll[0].y;
			lr[i].x = r*cos(ang[1]) + lr[0].x;
			lr[i].x = periodic(lr[i].x);
			lr[i].y = r*sin(ang[1]) + lr[0].y;
		}
	}

	// position passive
	for (i=0;i<N_PASSIVE;i++) {
		p[i].x += MU_PASSIVE*force_p[i][0]*DT;
		p[i].x = periodic(p[i].x);
		p[i].y += MU_PASSIVE*force_p[i][1]*DT;
	}

	if (USE_ACT_FRAME) {
		double disp = ll[0].x;

		// passive
		for (i=0;i<N_PASSIVE;i++)
			p[i].x = periodic(p[i].x - disp);

		// hinge
		for (i=0;i<N_BEADS;i++) {
			ll[i].x = periodic(ll[i].x - disp);
			lr[i].x = periodic(lr[i].x - disp);
		}
    		
	}

}

// [Initial setup] Finding possible positions for the passive particles
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
	for (i=0;i<N_BEADS;i++) {
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

// Function to initialize the system
void init_particles(passive_particle *p, passive_particle *ll, passive_particle *lr, int *head, int *lscl){

	int i,j;

	ang[0] = (CHECK_FROM_SECOND ? 0.0 : 1.*ANG_MAX);
	ang[1] = 1.*ANG_MAX;

	// Hinge's rods
	double r;
	for (i=0;i<N_BEADS;i++) {
		r = 1.*i*R_LEVER;
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
	for (i=0;i<N_PASSIVE;i++) {
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
