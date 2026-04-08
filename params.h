#ifndef PARAMS_H
#define PARAMS_H

#define PI 		3.14159265359

#define LX      10.*SIGMA // Size of the simulation box
#define LY      10.*SIGMA

/** SIMULATION PARAMETERS - SCALING UNITS **/
#define KBT		10.0 // Temperature in units of Boltzmann constant
#define SIGMA 	1.0 // Length unit (radius of the particles)

/** SIMULATION PARAMETERS - ACTIVE RODS **/
#define L_ROD       2.0*SIGMA// Length of the active rods
#define PE 			9. // Péclet number
#define ANG_MAX     2.*PI/6. // Maximum angle
#define MU_BEADS 	8.0 // Mobility of the beads
#define MUR_BEADS 	0.01 // Mobility of the beads
#define RATIO_MU 1.0

#define F_ACT       1.*PE*KBT // Active force magnitude

/** SIMULATION PARAMETERS - PASSIVE PARTICLES **/
#define R_DISK   	0.5*SIGMA // Radius of the passive particles
#define PF          0.65 // Packing fraction of the passive particles
#define MU_PASSIVE 	0.05 // Mobility of the passive particles

#define N_PASSIVE   (int) round(PF*LX*LY/(PI*R_DISK*R_DISK)) // Number of active particles

/** SIMULATION PARAMETERS - SYSTEM**/
#define DT      0.00005 //0.00001

#define N_SIM   200 // Number of simulations
#define T_MAX   100000000 // Simulation time in steps
#define T_EQ    500000 // Equilibration time in steps
#define N_SKIP  2000 // Number of steps to skip for measures

#define USE_ACT_FRAME 	0 // Use the hinge's frame of reference

#define K_WCA   50.//2.5*KBT  // WCA potential parameter

/** SIMULATION FLAGS **/
#define CHECK_FROM_SECOND   1 // Start with the first rod already absorbed?
#define CHECK_SAVE_DET      0 // Save detailed information?
#define CHECK_LINUX_PLOT    0 // Plot configuration using DynamicSimulator?

/** FILE NAMES **/
#define FNAME_MEA	"data/mea_L%.1f-pf%.2f-F%.1f-ang%.3f-mub%.2f-mup%.2f-T%.1f.dat"
#define FNAME_DET	"data/ang_L%.1f-pf%.2f-F%.1f-ang%.3f-mub%.2f-mup%.2f-T%.1f.dat"

/** FIXED PARAMETERS **/

#define R_LEVER	    0.1*SIGMA // Radius of active rod's beads
#define N_BEADS	    (int) round(1.*L_ROD/R_LEVER-1) // Number of those particles

#define PLOT_SCALE	1000. // Scale factor for the plot

#define PADDING		0.9

/** Linked-cells **/
#define LA_X          (int) floor(1.*LX/(2.*R_DISK))
#define RA_X          1.*LX/LA_X
#define LA_Y          (int) floor(1.*LY/(2.*R_DISK))
#define RA_Y          1.*LY/LA_Y

#endif
