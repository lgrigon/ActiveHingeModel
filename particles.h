#ifndef PARTICLES_H
#define PARTICLES_H

typedef struct {
    double x, y;
} passive_particle;

extern double ang[2];

void update(passive_particle *p, passive_particle *ll, passive_particle *lr, int move, int *head, int *lscl);
void find_possible_positions(passive_particle *pos, int *avaiPos, passive_particle *ll, passive_particle *lr);
void init_particles(passive_particle *p, passive_particle *ll, passive_particle *lr, int *head, int *lscl);

#endif
