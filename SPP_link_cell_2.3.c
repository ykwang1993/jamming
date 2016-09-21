
/*
Author: Kuan
Date: 2016.09.03
Version: 2.3
Brief: moving boundary

os: windows 7/10
Compiler: MinGw
Compile command: gcc SPP_link_cell_2.3.c SPP_force.c SPP_sub.c -o jam -std=c99

os: Ubuntu
Compiler: gcc
Compile command: gcc SPP_link_cell_2.3.c SPP_force.c SPP_sub.c -o jam -std=gnu99 -lm




*/

#include "SPP.h"

real rzta=1./(3.*PI*ETA*D); //rzta=119.21718
char dir[80];
real kl_long,kl_short,k_fene;


int main(int argc, char *argv[])
{

	char ch; 
	while ((ch = getopt(argc, argv, "p:f")) != EOF)
		switch (ch) {
			case 'p': 
				printf("%c %s\n", ch, optarg);
				strcpy(dir,optarg);
				printf("%s\n", dir);
				break;
			case 'f':
				printf("%c\n", ch);
				break;
		}

    int nc[DIM];
    int N, pnc;

    real l[DIM], r_cut;
    real delta_t, t_end;
    inputParameters_LC(&delta_t, &t_end, &N, nc, l, &r_cut);
    pnc=1;
    for (int d=0; d<DIM; d++)
        pnc *= nc[d];
    Cell *grid = (Cell*)malloc(pnc*sizeof(Cell));

    Particle **order = malloc(N*sizeof(Particle*));

    initData_LC(N, grid, nc, l, order);
    timeIntegration_LC(0, delta_t, t_end, grid, nc, l, r_cut, N, order);
    freeLists_LC(grid, nc);

    free(grid);
    free(order);

    

    printf("Finish");
    return 0;

}

void inputParameters_LC(real *delta_t, real *t_end, int *N, int *nc, real *l, real *r_cut){


	char path[80];
	
	scanf("%d %lf %lf %lf",N,&kl_long,&kl_short,&k_fene);

    l[0] = 120.;
    l[1] = 120.;
    *delta_t = 5e-5;
    *r_cut = 2.5;
    *t_end = 120.;

    nc[0] = (int)floor(l[0]/(*r_cut));
    nc[1] = (int)floor(l[1]/(*r_cut));
	
	int pnc=1;
    for (int d=0; d<DIM; d++)
        pnc *= nc[d];
	
	sprintf(path,"%s/para.txt",dir);
	
	FILE *file = fopen(path,"w");
	
	fprintf(file,"total N: %d\n",*N);
	fprintf(file,"t_end: %e\n",*t_end);
	fprintf(file,"delta_t: %e\n",*delta_t);
	fprintf(file,"kl_long: %e\n",kl_long);
	fprintf(file,"kl_short: %e\n",kl_short);
	fprintf(file,"k_fene: %e\n",k_fene);
	fprintf(file,"nc=[%d,%d] pnc=%d",nc[0],nc[1],pnc);
	
	fclose(file);

}

void RanInitial(Particle **order,int N)
{
	int num_chain = N/NR;
	int index_chain = 0;
	int index_bead = 0;
	Particle *temp_arr[NR+1];
	
	srand(time(NULL));
	
	for (index_chain=0;index_chain<=num_chain-1;index_chain++)
	{
		if(((rand()%2)==1) && (index_chain >= NS/NR))
		{			
			for (index_bead=1;index_bead<=NR;index_bead++)
			{
				temp_arr[index_bead] = order[NR*index_chain+index_bead-1];
				
			}
			for (index_bead=1;index_bead<=NR;index_bead++)
			{
				order[NR*index_chain+index_bead-1] = temp_arr[NR-index_bead+1];
				
			}
		}
	}
	
}

void initData_LC(int N, Cell *grid, int *nc, real *l, Particle **order){

	char path[80];
    int pnc=1,kc[DIM],iniCol=96,yCount=1,ySign=1;
    real x0,y0,db=1.12*D;
    ParticleList *current;


    for (int d=0; d<DIM; d++)
        pnc *= nc[d];

    for (int i=0; i<pnc; i++)
        grid[i]=NULL;

    x0=0.5*db;
    y0=60.;

    for (int i=1; i<=N; i++){

        current = (ParticleList*)malloc(sizeof(ParticleList));
        order[i-1]=&(current->p);

        if (i==1){
            current->p.x[0]=x0;
            current->p.x[1]=y0;
            current->p.v[0]=0.;
            current->p.v[1]=0.;
        }
        else{
            current->p.x[0]=order[i-2]->x[0]+db;
            current->p.x[1]=0.;
            current->p.v[0]=0.;
            current->p.v[1]=0.;


            if (i%iniCol == 1) {
                current->p.x[0]=x0;
                current->p.x[1]=y0-(ySign*yCount*db);
                y0=current->p.x[1];
				ySign*=-1;
				yCount+=1;
            } else current->p.x[1]=y0;
        }

        current->next=NULL;

        for (int d=0; d<DIM; d++)
            kc[d] = (int)floor(current->p.x[d] * nc[d] / l[d]);

        if(NULL==grid[index(kc,nc)])
            grid[index(kc,nc)]=current;
        else
            insertList(&grid[index(kc,nc)],current);

        
    }

    
	
	sprintf(path,"%s/ini.txt",dir);
	RanInitial(order, N);
	
	FILE *file = fopen(path,"w");
	outputResults_LC(N,order,file);
	fclose(file);
	
    printf("initData_LC\n");

}

real makePBC(real x1, real x2, real xsize){
    real x=x1-x2;

    if (fabs(x) > 0.5*xsize) {
        if (x > 0.) x=x-xsize;
        else if (x < 0.) x=x+xsize;
    } return x;
}



void compF_Intr(Particle **order, int N, real *l, long *seed, real thrust){

    int i, rodNum;
    real r_cut = 1.25*D, r;

    Particle *current, *neighbor;



    for (i=1;i<=N;i++){

        //FENE


        if (i > NS) rodNum=i%NR;
        else rodNum=i%NS;

        current = order[i-1];

        if (rodNum > 1 && rodNum != 0) {

                r=0.;
                neighbor = order[i-2];
                for (int d=0; d<DIM; d++)
                    r += sqr(makePBC(current->x[d],neighbor->x[d],l[d]));
                if (r<sqr(r_cut))
                    force_FENE(current,neighbor,r,l);

                r=0.;
                neighbor = order[i];
                for (int d=0; d<DIM; d++)
                    r += sqr(makePBC(current->x[d],neighbor->x[d],l[d]));
                if (r<sqr(r_cut))
                    force_FENE(current,neighbor,r,l);


        } else {
            if (rodNum == 1) {
                neighbor = order[i]; //left side

                r=0.;
                for (int d=0; d<DIM; d++)
                    r += sqr(makePBC(current->x[d],neighbor->x[d],l[d]));
                if (r<sqr(r_cut))
                    force_FENE(current,neighbor,r,l);

            } else if (rodNum == 0) { //right side
                neighbor = order[i-2];

                r=0.;
                for (int d=0; d<DIM; d++)
                    r += sqr(makePBC(current->x[d],neighbor->x[d],l[d]));
                if (r<sqr(r_cut))
                    force_FENE(current,neighbor,r,l);

            }

        }



        //Gaussian
        real ranForce;
        real dc=KT*rzta;
        real invDT=1./DT;


        for (int d=0; d<DIM; d++){
            ranForce=sqrt(2.*dc*invDT)*gasdev(seed);  //about 140
            order[i-1]->F[d] += ranForce;
        }
        //if (i == N) printf("%e\n",ranForce);

        //Thrusting force
        
        real thrustAng;

        thrustAng=getAngle(i,order,l,seed);
        order[i-1]->F[0] += rzta*thrust*cos(thrustAng);
        order[i-1]->F[1] += rzta*thrust*sin(thrustAng);


        //bending force
        //force_Bend(i,order,l);


        for (int j=1;j<=2;j++){
            real dx=1.0e-4*D; // infinitesimal length for spatial derivatives
            real dxdd=1./(2.*dx);
            real bending = calBend(i,j,order,l,1)-calBend(i,j,order,l,-1);
            bending *= dxdd;
			
			order[i-1]->F_bend[j-1] = -rzta*bending;
            order[i-1]->F[j-1] += order[i-1]->F_bend[j-1];
			
        }




    }
}

void compF_LC(Cell *grid, int *nc, int *nc_moving, real *l, real r_cut, int PBC_state) {
    int ic[DIM], kc[DIM], kc_temp[DIM],flag;
    for (ic[0]=0; ic[0]<nc_moving[0]; ic[0]++)
        for (ic[1]=0; ic[1]<nc_moving[1]; ic[1]++)
#if 3==DIM
            for (ic[2]=0; ic[2]<nc_moving[2]; ic[2]++)
#endif
            for (ParticleList *i=grid[index(ic,nc)]; NULL!=i; i=i->next) {
                for (int d=0; d<DIM; d++)  //set force as 0
                    i->p.F[d] = 0;
				
				//Reflecting Boundaries
				Particle wall_particle;
				real r_wall=0.;
				int flag_wall=0;
				if(PBC_state)
					for (int d=0; d<DIM; d++){
						flag_wall=0;
						if(ic[d]==(nc_moving[d]-1)){
							wall_particle.x[d] = 2*l[d]-(i->p.x[d]);
							r_wall = l[d]-(i->p.x[d]);
							flag_wall=1;
						}
						if(ic[d]==0){
							wall_particle.x[d] = -(i->p.x[d]);
							r_wall = i->p.x[d];
							flag_wall=1;
						}
						if(flag_wall)	
							if(2*r_wall<=r_cut)
								force_WCA_wall(&i->p, &wall_particle, sqr(2*r_wall), d, l);
					}	 
				
				
				//Deal with the neighbors
                for (kc[0]=ic[0]-1; kc[0]<=ic[0]+1; kc[0]++)
                    for (kc[1]=ic[1]-1; kc[1]<=ic[1]+1; kc[1]++)
#if 3==DIM
                        for (kc[2]=ic[2]-1; kc[2]<=ic[2]+1; kc[2]++)
#endif
                        {
                            //treat kc[d]<0 and kc[d]>=nc[d] according to boundary conditions;

                            //PBC
                            flag=0;
							for (int d=0; d<DIM; d++){
								kc_temp[d]=kc[d];
								if (kc[d]<0){
                                    kc_temp[d]=nc_moving[d]-1;
                                    flag=1;
								}

								if (kc[d]>=nc_moving[d]){
                                    kc_temp[d]=0;
                                    flag=1;
								}
							}

							//if (distance of i->p to cell kc <= r_cut)
								for (ParticleList *j=grid[index(kc_temp,nc)];NULL!=j; j=j->next)
									if (i!=j) {
										real r = 0;
										for (int d=0; d<DIM; d++)
                                            if (flag==0)
                                                r += sqr(j->p.x[d] - i->p.x[d]);
                                            else //PBC
                                                //r += sqr((j->p.x[d] - l[d]) - i->p.x[d]);
                                                r += sqr(makePBC(j->p.x[d], i->p.x[d], l[d]));
										if (r<=sqr(r_cut))
											force_WCA(&i->p, &j->p, r, l);

                                }
                }
            }
}

void PBC_fail(Cell *grid, int *nc)
{
	int ic[DIM];
	char path[80];
	sprintf(path,"%s/PBC_fail.txt",dir);
	
	FILE *file = fopen(path,"w");
	
	for (ic[0]=0; ic[0]<nc[0]; ic[0]++)
        for (ic[1]=0; ic[1]<nc[1]; ic[1]++)
#if 3==DIM
            for (ic[2]=0; ic[2]<nc[2]; ic[2]++)
#endif
                for (ParticleList *i=grid[index(ic,nc)]; NULL!=i; i=i->next)
				{
                    fprintf(file,"%e    %e    %e    %e    %e    %e    %e    %e\n"\
					,i->p.x_old[0],i->p.x_old[1]\
					,i->p.F_wca[0],i->p.F_wca[1]\
					,i->p.F_bend[0],i->p.F_bend[1]\
					,i->p.F_fene[0],i->p.F_fene[1]);
				
				}
	
	fclose(file);
	
}

void moveParticles_LC(Cell *grid, int *nc, int *nc_moving, real *l) {
    int ic[DIM], kc[DIM];
    for (ic[0]=0; ic[0]<nc_moving[0]; ic[0]++)
        for (ic[1]=0; ic[1]<nc_moving[1]; ic[1]++)
#if 3==DIM
            for (ic[2]=0; ic[2]<nc_moving[2]; ic[2]++)
#endif
        {   ParticleList **q = &grid[index(ic,nc)]; // pointer to predecessor
            ParticleList *i = *q;
            while (NULL != i) {

                //PBC
                for (int d=0; d<DIM; d++){
                    if (i->p.x[d] >= l[d]) i->p.x[d]-=l[d];
                    else if (i->p.x[d] < 0.) i->p.x[d]+=l[d];

                    if (i->p.x[d]<0. || i->p.x[d]>=l[d]){
                        printf("PBC fix fail\n");
                        printf("i->p.x[d] = %e\n",i->p.x[d]);
						printf("l[0]=%e\n",l[0]);
						PBC_fail(grid,nc);
                        exit(-2);
                    }
                }
				
				for (int d=0; d<DIM; d++){
                    kc[d] = (int)floor(i->p.x[d]*nc_moving[d]/l[d]);
                }


                
                if ((ic[0]!=kc[0])||(ic[1]!=kc[1])
                #if 3==DIM
                || (ic[2]!=kc[2])
                #endif
                ) {
                    deleteList(q);
                    insertList(&grid[index(kc,nc)], i);
                } else q = &i->next;
                i = *q;
            }
    }
}

void updateX(Particle *p, real delta_t) {

	for (int d=0; d<DIM; d++) {
		
		p->x_old[d] = p->x[d];
		p->x[d] += DT*p->F[d];
	}
}

void updateV(Particle *p, real delta_t) {
	real a = delta_t * .5 / p->m;
	for (int d=0; d<DIM; d++)
		p->v[d] += p->F[d];
}

void compX_LC(Cell *grid, int *nc, int *nc_moving, real *l, real delta_t) {
    int ic[DIM];
    for (ic[0]=0; ic[0]<nc_moving[0]; ic[0]++)
        for (ic[1]=0; ic[1]<nc_moving[1]; ic[1]++)
#if 3==DIM
            for (ic[2]=0; ic[2]<nc_moving[2]; ic[2]++)
#endif
                for (ParticleList *i=grid[index(ic,nc)]; NULL!=i; i=i->next)
                    updateX(&i->p, delta_t);
    moveParticles_LC(grid, nc, nc_moving, l);
}

void compV_LC(Cell *grid, int *nc, int *nc_moving, real *l, real delta_t) {
    int ic[DIM];
    for (ic[0]=0; ic[0]<nc_moving[0]; ic[0]++)
        for (ic[1]=0; ic[1]<nc_moving[1]; ic[1]++)
#if 3==DIM
            for (ic[2]=0; ic[2]<nc_moving[2]; ic[2]++)
#endif
                for (ParticleList *i=grid[index(ic,nc)]; NULL!=i; i=i->next)
                    updateV(&i->p, delta_t);
}


void outputResults_LC(int N, Particle **order, FILE *file){

    for (int i=0;i<N;i++){
        fprintf(file,"%e    %e    %e    %e\n",order[i]->x[0],order[i]->x[1],order[i]->v[0],order[i]->v[1]);
    }

}

void moving_boundary(Cell* grid, int *nc, int *nc_moving, real* l, real dx, real r_cut){
	
	int nc_old[2];
	for(int d=0;d<DIM;d++){
		nc_old[d]=nc_moving[d];
		l[d]-=dx;
		nc_moving[d] = (int)floor(l[d]/r_cut);	
	}
	
	int ic[DIM];
	int ic_target[DIM];
	if((nc_old[0]!=nc_moving[0])||(nc_old[1]!=nc_moving[1])){
		
				
		ParticleList *temp;
		ic[1]=nc_old[1]-1;
		ic_target[1]=nc_moving[1]-1;
		for(ic[0]=0;ic[0]<nc_old[0];ic[0]++){	
			ic_target[0]=ic[0];
			for (ParticleList *i=grid[index(ic,nc)]; NULL!=i;){
				temp=i->next;
				insertList(&grid[index(ic_target,nc)], i);
				i=temp;
			}
			grid[index(ic,nc)]=NULL;
		}
		
		ic[0]=nc_old[0]-1;
		ic_target[0]=nc_old[0]-2;
		for(ic[1]=0;ic[1]<nc_old[1]-1;ic[1]++){	
			ic_target[1]=ic[1];
			for (ParticleList *i=grid[index(ic,nc)]; NULL!=i;){
				temp=i->next;
				insertList(&grid[index(ic_target,nc)], i);
				i=temp;
			}
			grid[index(ic,nc)]=NULL;	
		}
			
		
			
				
	}
	
	

	
	
	
}

void timeIntegration_LC(real t, real delta_t, real t_end, Cell* grid, int *nc, real* l, real r_cut, int N, Particle **order) {

    int count=0,file_name=0;
    char path[12];
	int PBC_state=1; 
	real dx=0.0001;
	real thrust=0.01;
	int nc_moving[DIM];
	for(int d=0;d<DIM;d++)
		nc_moving[d]=nc[d];
	

    long storeTime,*seed=NULL;
    storeTime=-time(NULL);
    seed=&storeTime;

	//compF_LC(grid,nc,r_cut);
	while (t < t_end) {
		t += delta_t;
		
		
		//if(((float)count)/200==250.) PBC_state=0;
		if(l[0]<90. && PBC_state){
			PBC_state=0;
			thrust=0.05;
		}
			
		
		compF_LC(grid,nc,nc_moving,l,r_cut,PBC_state);
		compF_Intr(order,N,l,seed,thrust);
		compV_LC(grid,nc,nc_moving,l,delta_t);
		compX_LC(grid,nc,nc_moving,l,delta_t);
		
		if(l[0]>90.){
			if(l[0]<100.)
				dx=0.00002;
			moving_boundary(grid,nc,nc_moving,l,dx,r_cut);
		}
		
		
			



		count++;
		if(count%200 == 0){
            printf("%d\n",file_name);
            sprintf(path,"%s/%03d.txt",dir,file_name);
            FILE *file = fopen(path, "w");
            outputResults_LC(N, order, file);
            fclose(file);
            file_name++;
		}
	}
	

}


