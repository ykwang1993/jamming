
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>


#define DIM 2
#define sqr(x) ((x)*(x))
#define D 1.00 //bead diameter
#define ETA 0.89e-3
#define PI 3.141592653590
#define NS 4 //Length of long chain
#define NR 4 //Length of short chain
#define KT 4.12e-3
#define DT 5.0e-5
#define CP 0.05


#if 1==DIM
#define index(ic,nc) ((ic)[0])
#elif 2==DIM
#define index(ic,nc) ((ic)[0] + (nc)[0]*(ic)[1])
#elif 3==DIM
#define index(ic,nc) ((ic)[0] + (nc)[0]*((ic)[1] + (nc)[1]*(ic)[2]))
#endif



typedef double real;

typedef struct {
	real m; // mass
	real x[DIM]; // position
	real v[DIM]; // velocity
	real F[DIM]; // force
	
	real x_old[DIM];

	real F_fene[DIM];
	real F_bend[DIM];
	real F_wca[DIM];
} Particle;

typedef struct ParticleList {
    Particle p;
    struct ParticleList *next;
} ParticleList;

typedef ParticleList* Cell;

extern real rzta; //rzta=119.21718
extern char dir[80];
extern real kl_long,kl_short,k_fene;


void outputResults_LC(int N, Particle **order, FILE *file);
real ran1(long *idum);
real gasdev(long *idum);

void insertList(ParticleList **root_list, ParticleList *i);
void deleteList(ParticleList **q);
void freeLists_LC(Cell* grid, int *nc);

void inputParameters_LC(real *delta_t, real *t_end, int *N, int *nc, real *l, real *r_cut);
void RanInitial(Particle **order,int N);
void initData_LC(int N, Cell *grid, int *nc, real *l, Particle **order);
void PBC_fail(Cell *grid, int *nc);

real makePBC(real x1, real x2, real xsize);
void moveParticles_LC(Cell *grid, int *nc, int *nc_moving, real *l);

void force_WCA_wall(Particle *i, Particle *j, real r, int d, real *l);
void force_WCA(Particle *i, Particle *j, real r, real *l);
void force_FENE(Particle *i, Particle *j, real r, real* l);
real calBend(int n, int j, Particle **order, real *l, int a);
real vecAngle(real x, real y);
real getAngle(int n, Particle **order, real *l, long *idum);

void compF_Intr(Particle **order, int N, real *l, long *seed, real thrust);
void compF_LC(Cell *grid, int *nc, int *nc_moving, real *l, real r_cut, int PBC_state);
void updateX(Particle *p, real delta_t);
void updateV(Particle *p, real delta_t);
void compX_LC(Cell *grid, int *nc, int *nc_moving, real *l, real delta_t);
void compV_LC(Cell *grid, int *nc, int *nc_moving, real *l, real delta_t);
void moving_boundary(Cell* grid, int *nc, int *nc_moving, real* l, real dx, real r_cut);
void timeIntegration_LC(real t, real delta_t, real t_end, Cell* grid, int *nc, real* l, real r_cut, int N, Particle **order);















































