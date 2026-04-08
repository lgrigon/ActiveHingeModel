#ifndef RANDOM_H
#define RANDOM_H

#include <stdlib.h>
#include <time.h>
#include <math.h>

extern unsigned int myrand;
extern unsigned myrand, ira[256];
extern unsigned char ip, ip1, ip2, ip3;

#define SEED    time(NULL)
#define FNORM   (2.3283064365e-10)
#define RANDOM  ((ira[ip++] = ira[ip1++] + ira[ip2++]) ^ ira[ip3++])
#define FRANDOM (FNORM * RANDOM)
#define pm1     ((FRANDOM > 0.5) ? 1 : -1)

unsigned rand4init(void);
void Init_Random(void);
double rand_norm(double std_dev, double mean);

#endif
