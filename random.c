#include "random.h"

unsigned myrand, ira[256];
unsigned char ip, ip1, ip2, ip3;

/*---------------------------------------------------------------------*/
/*----------P. Random Number Generator by Parisi & Rapuano-------------*/
/*---------------------------------------------------------------------*/

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

//Box-Muller method
double rand_norm(double std_dev, double mean) {
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