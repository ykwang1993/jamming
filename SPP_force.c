
#include "SPP.h"

void force_WCA_wall(Particle *i, Particle *j, real r, int d, real *l){
	
	real sigma=1.;
    real epsilon=4.12e-3;

    real s = sqr(sigma) / r;
    s = sqr(s) * s;
    real f = 24 * epsilon * s / r * (1 - 2 * s);
	
	i->F[d] += f * rzta * makePBC(j->x[d],i->x[d],l[d]);
	
	
}


void force_WCA(Particle *i, Particle *j, real r, real *l) {

    real sigma=1.;
    real epsilon=4.12e-3;

    real s = sqr(sigma) / r;
    s = sqr(s) * s;
    real f = 24 * epsilon * s / r * (1 - 2 * s);
    for (int d=0; d<DIM; d++){
        
		i->F_wca[d] = f * rzta * makePBC(j->x[d],i->x[d],l[d]);
		i->F[d] += i->F_wca[d];
        
    }

}


void force_FENE(Particle *i, Particle *j, real r, real* l) {

    real R=1.25*D;
    real k=k_fene*KT/(D*D);

    real f = k/(1.-(r/sqr(R)));
    if ((1.-(r/sqr(R)))<= 1.0e-6){
        f = 0.;
        printf("FENE: out of r_cut\n");
    }
    for (int d=0; d<DIM; d++){

		i->F_fene[d] = f * rzta * makePBC(j->x[d],i->x[d],l[d]);
        i->F[d] += i->F_fene[d];
        

    }
}

real vecAngle(real x, real y)
/* Return the angle of a vector with respect to the horizontal axis in
   a value within the close interval [-PI,PI]. The input values are the
   x and y component of the vector, respectively. */
{
    real ang;

    if (fabs(x) <= 1.0e-6) { /* avoid NaN when calculating y/x with x=0 */
        if (y > 0.) ang=0.5*PI;
        else if (y < 0.) ang=-0.5*PI;
    } else if (x > 1.0e-6) ang=atan(y/x);
    else if (x < -1.0e-6) {
        if (y >= 0.) ang=atan(y/x)+PI;
        else if (y < 0.) ang=atan(y/x)-PI;
    } return ang;
}

real getAngle(int n, Particle **order, real *l, long *idum)
{
    int rodNum;
    real x1,x2,vSubr,tmp,coefDev,ang,coefVel=CP;

    if (n > NS) rodNum=n%NR;
    else rodNum=n%NS;
    if (rodNum > 1 && rodNum != 0) {
        x1=makePBC(order[n]->x[0], order[n-2]->x[0], l[0]);
        x2=makePBC(order[n]->x[1], order[n-2]->x[1], l[1]);
    } else if (rodNum == 1) {
        x1=makePBC(order[n]->x[0], order[n-1]->x[0], l[0]);
        x2=makePBC(order[n]->x[1], order[n-1]->x[1], l[1]);
    } else if (rodNum == 0) {
        x1=makePBC(order[n-1]->x[0], order[n-2]->x[0], l[0]);
        x2=makePBC(order[n-1]->x[1], order[n-2]->x[1], l[1]);
    }

    ang=vecAngle(x1,x2);
    vSubr=sqrt(sqr(order[n-1]->v[0])+sqr(order[n-1]->v[1]));
    tmp=1./(1.+coefVel*vSubr);
    coefDev=2.*PI*tmp;
    //ang+=coefDev*(ran1(idum)-0.5);
    return ang;
}


real calBend(int n, int j, Particle **order, real *l, int a)
{
    real dx = 1.0e-4*D;
    int j2=(j+1)%2,rodNum;
    real dev=a*dx,x1,x2,y1,y2,rx,ry,dotVec,ang;
    real kl=0.;
	
	if (n > NS) 
	{
		rodNum=n%NR;
		kl=kl_short*D*KT;
	}
    else
	{
		rodNum=n%NS;
		kl=kl_long*D*KT;	
	}


    if (j2 == 0) j2=2;
    if (rodNum > 1 && rodNum != 0) {
        x1=makePBC(order[n-1]->x[j-1],order[n-2]->x[j-1],l[j-1]);
        x2=makePBC(order[n-1]->x[j2-1],order[n-2]->x[j2-1],l[j-1]);
        x1+=dev;
        rx=sqrt(x1*x1+x2*x2);
        y1=makePBC(order[n]->x[j-1],order[n-1]->x[j-1],l[j-1]);
        y2=makePBC(order[n]->x[j2-1],order[n-1]->x[j2-1],l[j-1]);
        y1-=dev;
        ry=sqrt(y1*y1+y2*y2);
        dotVec=(x1*y1+x2*y2)/(rx*ry);
        if (dotVec > 1.) dotVec=1.;
        else if (dotVec < -1.) dotVec=-1.;
        ang=acos(dotVec);
        return kl*ang*ang;
    } else return 0.;
}