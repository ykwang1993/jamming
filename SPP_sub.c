
#include "SPP.h"

void insertList(ParticleList **root_list, ParticleList *i) {
    i->next = *root_list;
    *root_list = i;
}

void deleteList(ParticleList **q) {
    *q = (*q)->next; // (*q)->next points to element to be removed
}

void freeLists_LC(Cell* grid, int *nc){

    ParticleList *current, *prev;
    int ic[DIM];
    for (ic[0]=0; ic[0]<nc[0]; ic[0]++)
        for (ic[1]=0; ic[1]<nc[1]; ic[1]++)
#if 3==DIM
            for (ic[2]=0; ic[2]<nc[2]; ic[2]++)


#endif
            {
                current = grid[index(ic,nc)];
                 while(current != NULL){
                    prev=current;
                    current=current->next;
                    free(prev);

                }
            }

}



#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


/* "Minimal" random number generator of Park and Miller with Bays-Durham
   shuffle and added safeguards. Returns a uniform random deviate between 0.0
   and 1.0 (exclusive of the endpoint values). Call with idum a negative
   integer to initialize; thereafter, do not alter idum between successive
   deviates in a sequence. RNMX should approximate the largest floating
   value that is less than 1. */
/* (C) Copr. 1986-92 Numerical Recipes Software 7MZ9%"W5:!+). */
real ran1(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    real temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum=1;
        else *idum=-(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum+=IM;
            if (j < NTAB) iv[j]=*idum;
        } iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum+=IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j]=*idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

/* Returns a normally distributed deviate with zero mean and unit variance,
   using ran1(idum) as the source of uniform deviates. */
/* (C) Copr. 1986-92 Numerical Recipes Software 7MZ9%"W5:!+). */
real gasdev(long *idum)
{
    static int iset=0;
    static real gset;
    real fac,rsq,v1,v2;

    if (*idum < 0) iset=0;
    if (iset == 0) {
        do {
            v1=2.*ran1(idum)-1.;
            v2=2.*ran1(idum)-1.;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1. || rsq == 0.);
        fac=sqrt(-2.*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    } else {
        iset=0;
        return gset;
    }
}




