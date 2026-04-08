#ifndef UTILS_H
#define UTILS_H

#include "particles.h"

void fill_list(passive_particle *p, int *head, int *lscl);
double periodic(double dis);
double mod2pi(double ang);
double pc(double x);
void plot_config(passive_particle *ll, passive_particle *lr, passive_particle *p);
void get_N(int *Ns, passive_particle *p, passive_particle *ll, passive_particle *lr);

#endif