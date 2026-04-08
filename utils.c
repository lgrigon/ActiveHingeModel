#include "utils.h"
#include "params.h"
#include "particles.h"
#include "random.h"
#include <math.h>
#include <stdio.h>

// Fill lists used for the cell list method
void fill_list(passive_particle *p, int *head, int *lscl) {
	int mc[2];
	int i,c;

	for (i=0;i<LA_X*LA_Y;i++) head[i] = -1;
	for (i=0;i<N_PASSIVE;i++) {
		mc[0] = (int) floor((p[i].x + LX/2.)/(1.*RA_X));
		if (mc[0] < 0) mc[0] = 0; else if (mc[0] >= LA_X) mc[0] = LA_X - 1;
		mc[1] = (int) floor((p[i].y + LY/2.)/(1.*RA_Y));
		if (mc[1] < 0) mc[1] = 0; else if (mc[1] >= LA_Y) mc[1] = LA_Y - 1;
		c = mc[0]*LA_Y + mc[1];
		lscl[i] = head[c];
		head[c] = i;
	}
}

// Function to account for the periodic boundary condition
double periodic(double dis) {
	if (dis <= -LX/2.)	return dis+LX;
	if (dis >= LX/2.)	return dis-LX;
	return dis;
}

// Function to account for angle periodicit
double mod2pi(double ang) {
    ang = fmod(ang, 2.*PI);
    if (ang < 0) ang += 2.*PI;
    return ang;
}

// Convert coordinates to plot using DynamicSimulator
double pc(double x) {
	return (((x+1.*LX/2.)/LX)-1.)*(-1.);
}

// Plot the current configuration of particles
void plot_config(passive_particle *ll, passive_particle *lr, passive_particle *p) {
    int n, i;

    for (n=0;n<N_BEADS;n++)
    	printf("%f %f %f %f %f\n",pc(ll[n].x),pc(ll[n].y),2.*R_LEVER*PLOT_SCALE/LX,2.1*R_LEVER*PLOT_SCALE/LX,0.);

	for (n=0;n<N_BEADS;n++)
    	printf("%f %f %f %f %f\n",pc(lr[n].x),pc(lr[n].y),2.*R_LEVER*PLOT_SCALE/LX,2.1*R_LEVER*PLOT_SCALE/LX,0.);

    for (n=0;n<N_PASSIVE;n++)
    	printf("%f %f %f %f %f\n",pc(p[n].x),pc(p[n].y),2.*R_DISK*PLOT_SCALE/LX,2.*R_DISK*PLOT_SCALE/LX,0.);
}

// Find number of passive particles below the hinges
void get_N(int *Ns, passive_particle *p, passive_particle *ll, passive_particle *lr) {
  int i;
  Ns[0] = 0; Ns[1] = 0;
  double r, th;
  for (i=0;i<N_PASSIVE;i++) {
  	r = sqrt( pow(periodic(ll[0].x - p[i].x),2) + pow(ll[0].y - p[i].y,2) );
  	th = PI - atan2(p[i].y - ll[0].y, periodic(p[i].x - ll[0].x));
  	if (r < L_ROD && th < ang[0]) Ns[0]++;
  	r = sqrt( pow(periodic(lr[0].x - p[i].x),2) + pow(lr[0].y - p[i].y,2) );
  	th = atan2(p[i].y - lr[0].y, periodic(p[i].x - lr[0].x));
  	if (r < L_ROD && th < ang[1]) Ns[1]++;
  }
}