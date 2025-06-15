# Active Hinge Simulations

This repository contains simulation codes used in our study of **motile active particle clusters**. The system involves a pair of self-propelled rods connected by a hinge and embedded in a crowded medium. We provide two simulation setups:

- `double/`: Simulations where **both rods are alive**.
- `single/`: Simulations where **one rod is already absorbed** (asymmetric case).

📌 Note  
Output folders: (data/, mfpt/, temp/)

These codes were used to generate the results for our paper titled:

> Symmetry breaking in crowded channels  
> Leonardo Garibaldi Rigon and Yongjoo Baek  
> _Submitted to: Soft Matter Journal_  
> _[arXiv link]_

---

## Requirements

- GCC compiler
- [Gnuplot](http://www.gnuplot.info/) for fitting survival curves

---

## 📂 double/: Two Active Rods

▶️ Main Simulation (file: double_hinge.c)
```
gcc double_hinge.c -lm -O3
./a.out
```
Set parameters at the top of the file:
1. SIMULATION PARAMETERS section: typical values to change
2. FIXED PARAMETERS section: for advanced control

Output files are written to ../data/

Optional flags (set in code):
1. SAVE_CONFIGS = 1: saves particle positions
2. SAVE_DETAILS = 1: saves angles and extra details

📉 Estimating MFPT (files: survival_rate.c and fit_mfpt.gnu)
```
gcc survival_rate.c -lm -O3
./a.out
```
Reads measurement files from ../data/  
Outputs survival probabilities to ../mfpt/

Fitting with Gnuplot:
```
load 'fit_mfpt.gnu'
```
Adjust parameters in the .gnu file to match the simulation

📉 Calculating Mean Squared Displacement (MSD) (file: msd.c)

Requires: SAVE_DETAILS = 1 in double_hinge.c  
Outputs saved to ../data/
```
gcc msd.c -lm -O3
./a.out
```

## 📂 single/: One Rod Already Absorbed
▶️ Main Simulation (file: single_hinge.c)
```
gcc single_hinge.c -lm -O3
./a.out
```
Setup similar to double_hinge.c

Extra optional flags (set in code):
1. SAVE_MSD: computes MSD internally
2. SAVE_FINALS: saves data from the final 50 steps (forces, torque, etc.)
3. SAVE_TORQUE: saves torque data
4. SAVE_N0: logs information when particles below the hinge are 0 or 1

📉 Estimating MFPT

Same process as in the double/ folder

📉 Estimating Time in State T<sub>N</sub> (files: terms_surv_rate.c and fit_mfpt_terms.gnu)
```
gcc terms_surv_rate.c -lm -O3
./a.out
```
Outputs survival rates for each coarse-grained state to mfpt/ folder

Fitting with Gnuplot:
```
load 'fit_mfpt_terms.gnu'
```
