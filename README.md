# Active Hinge Simulations: Symmetry Breaking in Crowded Channels

This repository contains simulation codes used in our study of **motile active particle clusters**. The system involves a pair of self-propelled rods connected by a hinge and embedded in a crowded medium.

These codes were used to generate the results for our paper titled:

> Symmetry-breaking motility of an active hinge in a crowded channel  
> Leonardo Garibaldi Rigon and Yongjoo Baek  
> _[Soft Matter (2026)](https://doi.org/10.1039/D5SM00622H)_  

---

## Repository Structure

* **main.c**: The primary entry point for the simulation loop.
* **params.h**: Global simulation parameters and physical constants.
* **particles.c**: Core logic for position updates, hinge constraints, and particle interactions.
* **random.c / utils.c**: Helper functions for the random number generator and periodic boundary conditions.
* **survival_rate.c**: Post-processing script to calculate survival curves from raw data.
* **data/**: (Required) Target folder for raw measurement files (.dat).
* **mfpt/**: (Required) Target folder for processed survival rate outputs.
* **temp/**: (Required) Temporary storage for MFPT fits.

---

## Getting Started

### Prerequisites
* **GCC compiler**
* **[Gnuplot](http://www.gnuplot.info/)** for fitting survival curves

### Installations & Setup
1. Clone the repository.
2. Ensure the output directories exist:
```
mkdir -p data mfpt temp
```
3. Configure physical parameters in params.h. Key flags include:
- CHECK_FROM_SECOND: Set to 1 to start simulations with the first rod already "absorbed".
- CHECK_SAVE_DET: Set to 1 to output angle data and extra simulation details.

### Running the Simulation
To compile and run the main simulation:
```
gcc -O3 main.c particles.c random.c utils.c -lm -o hinge_sim
./hinge_sim
```

### Estimating MFPT
After generating data, use survival_rate.c to process the results:
```
gcc -O3 survival_rate.c -lm -o survival_tool
./survival_tool
```

### Fit MFPT in Gnuplot
```
load 'fit_mfpt.gnu'
```
Note: Ensure the parameters inside fit_mfpt.gnu match the ranges defined in your C code.
