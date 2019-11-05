#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <xwc.h>

struct sprite {
  PIC *P;
  PIC *Aux;
  MASK *Msk;
};

struct bullet{
    double m;
    double *c,*v,*a;
    int o;
    struct sprite* spt;
};

struct nave{
    char* nome;
    double m;
    double *c,*v,*a;
    int o;
    struct sprite* spt;
};

struct planet{
    double r,m,t;
};