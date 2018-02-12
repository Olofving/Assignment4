#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include "graphics.h"

//Quadrangle struct
typedef struct qtree{
	struct qtree* quad1; //northwest
	struct qtree* quad2; //southwest
	struct qtree* quad3; //southeast
	struct qtree* quad4; //northeast
	struct qinfo* qdata;
}qtree;

typedef struct qinfo{
	double mass;
	double posx;
	double posy;
	double posxm;
	double posym;
	double width;
}qinfo;

//Particle struct
typedef struct particle{
  double velx;
  double vely;
  double posx;
  double posy;
  double forcex;
  double forcey;
	double mass;
	double brightness;
	
}particle;


// function for insertion of particle into qtree
void insert(qtree ** head, qinfo * p, int quadrant){

	//Checks if quadrant is empty and adds first node if so	
	if((*head) == NULL){
	  (*head) = malloc(sizeof(qtree));
	  (*head)->qdata = malloc(sizeof(qinfo));

	if(quadrant ==1 ){	
 			(*head)->qdata->posx= p->posx - p->width/2.0;
	  	(*head)->qdata->posy= p->posy - p->width/2.0;	  	
	  }else if(quadrant == 2){
	  	(*head)->qdata->posx= p->posx - p->width/2.0;
	  	(*head)->qdata->posy= p->posy + p->width/2.0;
	  }else if(quadrant == 4){
	  	(*head)->qdata->posx= p->posx + p->width/2.0;
	  	(*head)->qdata->posy= p->posy - p->width/2.0;
	  }else{
	  	(*head)->qdata->posx= p->posx + p->width/2.0;
	  	(*head)->qdata->posy= p->posy + p->width/2.0;
	  }
	  
	  
	  (*head)->qdata->mass = p->mass;
	  (*head)->qdata->posxm = p->posxm;
	  (*head)->qdata->posym = p->posym;
	  (*head)->qdata->width = p->width;
	  (*head)->quad1 = NULL;
	  (*head)->quad2 = NULL;
	  (*head)->quad3 = NULL;
	  (*head)->quad4 = NULL;
	  return;
	}
	
	// recursive function call to place particle 
	 else {	  
	  
	  // checks if there is only one particle in square and "sends it down" if so
	  if((*head)->quad1 == NULL && (*head)->quad2 == NULL && (*head)->quad3 == NULL && (*head)->quad4 == NULL){
	  
		  qinfo * newp = NULL;
		  newp = malloc(sizeof(qinfo)); 
		  newp->posx = (*head)->qdata->posx; 
		  newp->posy = (*head)->qdata->posy; 
		  newp->posxm = (*head)->qdata->posxm; 
		  newp->posym = (*head)->qdata->posym;
		  newp->mass = (*head)->qdata->mass;
		  newp->width = (*head)->qdata->width/2.0;   
		  
		  
		  if((*head)->qdata-> posxm <= (*head)->qdata->posx){
		 		if((*head)->qdata-> posym <= (*head)->qdata->posy){
		 			// quadrant 1  
		 			printf("quad 1");
					insert(&((*head)->quad1),newp,1);
					}
		 	else{
		 		// quadrant 2
		 		printf("quad 2");
				insert(&((*head)->quad2),newp,2);
		 	}
		} 
		else{
			if((*head)->qdata-> posym >= (*head)->qdata->posym){
				//quadrant 3
				printf("quad 3");
				insert(&((*head)->quad3),newp,3);
			}
		else{
			// quadrant 4
			printf("quad 4");
			insert(&((*head)->quad4),newp,4);
				}
			}
		}	  
		
		// uppdating mass and center of mass for square
		(*head)->qdata-> posxm= (p->posxm * p->mass + (*head)->qdata-> posxm*(*head)->qdata->mass)/ ( (*head)->qdata-> mass + p->mass);
		(*head)->qdata-> posym= (p->posym *p->mass + (*head)->qdata-> posym*(*head)->qdata->mass)/ ( (*head)->qdata-> mass + p->mass);
		(*head)->qdata-> mass+= p->mass;	  
		  
		 
		 //storing current squares center in p
		 p->posx =(*head)->qdata->posx;
		 p->posy =(*head)->qdata->posy;
		 p->width = p->width/2.0; 
		 
		 
		  // Checks which quadrangle the particle belongs to and inserts it there
		 if(p->posxm <= (*head)->qdata->posx){
		 	if(p->posym <= (*head)->qdata->posy){
		 		// quadrant 1  
				insert(&((*head)->quad1),p,1);
			 }
		 	else{
		 		// quadrant 2
				insert(&((*head)->quad2),p,2);
		 	}
		} 
		else{
			if(p->posym >= (*head)->qdata->posy){
				//quadrant 3
				insert(&((*head)->quad3),p,3);
			}
			else{
				// quadrant 4
				insert(&((*head)->quad4),p,4);
			}
		}
	}
}


void force_qtree(qtree * node, qinfo particle, double theta, double * forceqt){
	double epsilon = 0.001;
	double dist;
	double r;
	
	
	if(node != NULL){
	dist = sqrt( (node->qdata->posx - particle.posxm)*(node->qdata->posx - particle.posxm) + (node->qdata->posy - particle.posym)*(node->qdata->posy - particle.posym));
		if( node->qdata->width/dist  <= theta){
		
			r = sqrt((node->qdata->posxm - particle.posxm)*(node->qdata->posxm - particle.posxm) + (node->qdata->posym - particle.posym)*(node->qdata->posym - particle.posym)) +epsilon;
			if(r > epsilon){
				forceqt[0] += (node->qdata->mass)/(r*r*r)*(particle.posxm-node->qdata->posxm);				
				forceqt[1] += (node->qdata->mass)/(r*r*r)*(particle.posym-node->qdata->posym);
			}
		}else
		
		//OPTIMERA!!
		if(node->quad1 != NULL || node->quad2 != NULL || node->quad3 != NULL || node->quad4 != NULL){
			force_qtree(node->quad1,particle,theta,forceqt);
			force_qtree(node->quad2,particle,theta,forceqt);
			force_qtree(node->quad3,particle,theta,forceqt);
			force_qtree(node->quad4,particle,theta,forceqt);
		}else{
			r = sqrt((node->qdata->posxm - particle.posxm)*(node->qdata->posxm - particle.posxm) + (node->qdata->posym - particle.posym)*(node->qdata->posym - particle.posym)) +epsilon;
				if(r > epsilon){
					forceqt[0] += (node->qdata->mass)/(r*r*r)*(particle.posxm-node->qdata->posxm);				
					forceqt[1] += (node->qdata->mass)/(r*r*r)*(particle.posym-node->qdata->posym);
				}
			}		
		}
	}


// print quadtree
void print_quadtree(qtree *node)
{
   if(node == NULL) {printf("Tree is empty!\n"); return;}

   if (node != NULL) printf("mass: %lf. posx: %lf. posy: %lf. posxm: %lf, posym %lf.", node->qdata->mass,node->qdata->posx,node->qdata->posy,node->qdata->posxm,node->qdata->posym);

   printf("\n");

   if (node->quad1 != NULL)
      print_quadtree(node->quad1);
   if (node->quad2 != NULL)
      print_quadtree(node->quad2);
   if (node->quad3 != NULL)
      print_quadtree(node->quad3);
   if (node->quad4 != NULL)
      print_quadtree(node->quad4);
}

// removal of qtree
void remove_qtree(qtree ** node){
	if((*node) != NULL){
	
		remove_qtree(&(*node)->quad1);
		remove_qtree(&(*node)->quad2);
		remove_qtree(&(*node)->quad3);
		remove_qtree(&(*node)->quad4);
	
		free((*node)->qdata);
		(*node)->qdata = NULL;
		free(*node);
		*node = NULL;
	}
}


// function for timetaking
static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}





// uppdates all the velocities and positions, with timestep delta_t
void timestep_quadtree(qinfo *particle_list, double *velx, double *vely,double *forcex, double * forcey, int N, double delta_t){
	double accx;
	double accy;
	
	for(int i=0;i<N;i++){
		accx = forcex[i]/particle_list[i].mass;
		accy = forcey[i]/particle_list[i].mass;		
		velx[i] += delta_t*accx;
		vely[i] += delta_t*accy;		
		particle_list[i].posxm += delta_t*velx[i];
		particle_list[i].posym += delta_t*vely[i];
	}
}


 
/*
*
*		Main Mentods
*
*/


//test of quadtree
int main(){
qinfo a1;
qinfo a2;
qinfo a3;
qinfo a4;
qinfo b;
qinfo c;
qinfo d;

a1.mass=1;
a1.posx=1;
a1.posy=1;
a1.posxm=0.2;
a1.posym=0.2;
a1.width=1;

a2.mass=1;
a2.posx=1;
a2.posy=1;
a2.posxm=0.2;
a2.posym=0.4;
a2.width=1;

a3.mass=1;
a3.posx=1;
a3.posy=1;
a3.posxm=0.4;
a3.posym=0.4;
a3.width=1;

a4.mass=1;
a4.posx=1;
a4.posy=1;
a4.posxm=0.4;
a4.posym=0.2;
a4.width=1;

b.mass=1;
b.posx=1;
b.posy=1;
b.posxm=0.25;
b.posym=0.8;
b.width=1;

c.mass=1;
c.posx=1;
c.posy=1;
c.posxm=0.75;
c.posym=0.75;
c.width=1;

d.mass=1;
d.posx=1;
d.posy=1;
d.posxm=0.75;
d.posym=0.25;
d.width=1;

qtree * head = NULL;

head = NULL;
printf("inserting d: ");
insert(&head,&d,1);
printf("insert d ok!\n");
insert(&head,&b,1);
printf("insert b ok!\n");
insert(&head,&c,1);
printf("insert c ok!\n");
//insert(&head,&a1,1);
//insert(&head,&a2,1);
//insert(&head,&a3,1);
//insert(&head,&a4,1);
print_quadtree(head);
printf("\n");
remove_qtree(&head);
print_quadtree(head);
//}

return 0;
}


/*
int main(int argc, const char** argv){
	
	// Timer
	double time1 = get_wall_seconds();

	//Checking for correct ammount of arguments
	if(argc != 7) {
		printf("ERROR! Expected five input arguments: N filename nsteps delta_t thetamax graphics\n");
  	return -1;
	}

	//declaring variables for input args
	const int N = atoi(argv[1]);
	const char* filename = argv[2];
	const int nsteps = atoi(argv[3]);
	const double delta_t = atof(argv[4]);
	const double theta = atof(argv[5]);
	const int graphics = atoi(argv[6]);
	const double G = 100.0/( (double) N);
	double brightness[N];
	double velx[N];
	double vely[N];
	double forcex[N];
	double forcey[N];

	// array of all the particles
	qinfo *particle_list = (qinfo*)malloc(N*sizeof(qinfo));

	//opening and reading of file
	FILE* fptr;
	fptr = fopen(filename,"rb");

	int size_double=sizeof(double);
	
	size_t returnval = 0;
	for(int i =0;i<N;i++){  
		//reading data
		returnval = fread(&particle_list[i].posxm,size_double,1,fptr);
		returnval = fread(&particle_list[i].posym,size_double,1,fptr);
		returnval = fread(&particle_list[i].mass,size_double,1,fptr);
		returnval = fread(&velx[i],size_double,1,fptr);
		returnval = fread(&vely[i],size_double,1,fptr);
		
		//setting initial values for quadtree asscociated variables
		particle_list[i].posx=1;
		particle_list[i].posy=1;
		particle_list[i].width=1;
		
		//brightness is not used in simulation so is saved in a separate vector
		returnval = fread(&brightness[i],size_double,1,fptr);
	}
returnval++;
fclose(fptr);


// graphics related stuff

	const float circleRadius=0.005, circleColor=0;
	float L=1, W=1;
	if(graphics == 1){
		InitializeGraphics("FATAL ERROR!",1400,800);
		SetCAxes(0,1);
	}
	
	//Simulation
  int i;
	for(i=0; i<nsteps; i++){
	
		qtree * head = NULL;
	
		// creation of tree
		for(int j=0;j<N;j++){  
	
			qinfo * copy = NULL;
			copy = malloc(sizeof(qinfo));
			copy->posxm = particle_list[j].posxm;
			copy->posym = particle_list[j].posym;
			copy->posx = particle_list[j].posx;
			copy->posy = particle_list[j].posy;
			copy->mass = particle_list[j].mass;
			copy->width = particle_list[j].width;
		
			insert(&head, copy,1);
			free(copy);
			copy = NULL;
			}
			// print_quadtree(head);
			
			// calculation of the total force on particle k
			double forcetemp[2];

			for(int k =0; k<N; k++){
				forcetemp[0]=0;
				forcetemp[1]=0;
				force_qtree(head, particle_list[k],theta,forcetemp);
				forcex[k] = -G*particle_list[k].mass*forcetemp[0];
				forcey[k] = -G*particle_list[k].mass*forcetemp[1];
				//printf("forcex: %f. forcey: %f \n",forcex[k], forcey[k]);
			}
			
			// free quadtree
			remove_qtree(&head);
			
			// updating the/ position of the particles
			timestep_quadtree(particle_list,velx, vely, forcex,forcey,N,delta_t);
	

    	// Call graphics routines. 
	  	if(graphics == 1){
			ClearScreen();
	  	  for(int j=0;j<N;j++){
	  	  	DrawCircle(particle_list[j].posxm, particle_list[j].posym, L, W, circleRadius, circleColor);
	  	  }   
	  	     Refresh();
	  	     
	  	  // Sleep a short while to avoid screen flickering. 
	  	  usleep(30000);   
			}     
}


// Writing data to .gal-file	
	FILE *fp;
	fp = fopen("result.gal","wb");

	for(int i =0;i<N;i++){
		fwrite(&particle_list[i].posxm,size_double,1,fp);
		fwrite(&particle_list[i].posym,size_double,1,fp);
		fwrite(&particle_list[i].mass,size_double,1,fp);
		fwrite(&velx[i],size_double,1,fp);
		fwrite(&vely[i],size_double,1,fp);
		fwrite(&brightness[i],size_double,1,fp);
	}
	
	free(particle_list);
	fclose(fp);
	// printing of computation time in wall seconds
	printf("%7.3f\n",get_wall_seconds()-time1);	
		
	return 0;
}
*/


