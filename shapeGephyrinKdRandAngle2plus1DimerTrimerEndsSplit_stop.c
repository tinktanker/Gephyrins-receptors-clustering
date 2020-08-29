//Author: Wenjun Xia
//Created: Sep., 2018

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "kdtree.h"

#define N 500 // number of gephyrins
#define L 320 // size of lattice
#define Nr 100 // number of receptors
#define geLength 11
#define cir 3.0 //the diameter of circle end
#define squ 4.9 //the diameter of square end
#define Pd 0.1 // probability of terminal dimers off-binding of gephyrin clusters
#define Pt 0.1 // probability of terminal timers off-binding of gephyrin clusters
#define Pr 0.1 // probablity of receptor off-binding
#define Pc 0.1 // probablity of gephyrins' 2 circle heads combine
#define	tmax 0// duration of random walk
#define	diff 10.0  // diffusion coefficient of gephyrins, use even numbers
#define	diffr 10.0 //diffusion coefficient of receptors
//#define	dist 10.0  // critical distance for clustering of gephyrins
#define	distgr 3.0*squ // critical distance for receptor-gephyrin binding
#define divider 10//every 'divider' steps to produce results files
#define pi acos(-1.)

void bubbleSort(int *arr, int n) {
	int i,j;
    for (i = 0; i<n - 1; i++)
        for (j = 0; j < n - i - 1; j++)
        {
            ////if the latter smaller than the former
            if (arr[j] < arr[j + 1]) {
				int temp;
                temp = arr[j];
				arr[j] = arr[j + 1];
				arr[j + 1] = temp;
            }
        }
}

static double (*kdtree_knn_dump(struct kdtree *tree))[2]
{
    int i;
		static double dis[3][2];
		int count;
		count=0;
    struct knn_list *p = tree->knn_list_head.next;
    while (p != &tree->knn_list_head) {
        for (i = 0; i < tree->dim; i++) {
            if (i == tree->dim - 1) {
			         dis[count][0]=sqrt(p->distance);
				       dis[count][1]=p->node->coord_index;
			         count++;
            }
        }
        p = p->next;
    }
		return dis;
}

static void knn_clear(struct kdtree *tree)
{
	  int i;
		tree->knn_list_head.next = &tree->knn_list_head;
		tree->knn_list_head.prev = &tree->knn_list_head;
        tree->knn_list_head.node = NULL;
		tree->knn_list_head.distance = 0;
		for(i=0;i<tree->capacity;i++){
			tree->coord_passed[i]=0;
		}
		tree->knn_num=0;
}


int main(int argc, char* argv[]){


	//printf("pi=%17.16f\n",pi);
	if(argc < 3){
	 fprintf(stderr, "Give 1/0/2! 1: print details to screen; 0: not print; 2: only print gephyrin clusters' splitting info.\nGive seed/0! seed: for random number; 0: produce seeds by computer\n");
	 exit(-1);
	}
  int t;
	int par2;
	par2=atoi(argv[2]);
	if(par2==0){
		t=time(NULL);
	}else{
	    t=par2;
	}
	printf("t=%d\n",t);
	srand(t);



  int sil;
	sil=atoi(argv[1]);


	double *x,*x1,*y,*y1;//(x,y),circle head; (x1,y1),square head
	x=(double*)malloc(N*sizeof(double));
	x1=(double*)malloc(N*sizeof(double));
	y=(double*)malloc(N*sizeof(double));
	y1=(double*)malloc(N*sizeof(double));
	int i,j,tau,l,m,n,p,q;
	double xi,yi; //x-intercept of the point of intersection
	int *pnRand;
	pnRand=(int*)malloc(N*sizeof(int));
	double u1,u2,z0,z1;
	int *node,*node1;//node:circle end taken=the index of the adjacent circle, ie. i/j/l; node1:square end taken=the index of gephyrin who took it. or else=-1.
	node=(int*)malloc(N*sizeof(int));
	node1=(int*)malloc(N*sizeof(int));
	int **cls;//mark which gephyrins each cluster contains
	cls=(int**)malloc(N*sizeof(int*));//1st dimension
	for(i=0;i<(N+1); i++){
      cls[i]=(int*)malloc((N+1)* sizeof(int));//2nd dimension
  }
	double Gx,Gy;//core of a triangle
	int *beCls;//label which cluster each gephyrin belongs to(eg. cluster 1, cluster 2,...)
	beCls=(int*)malloc(N*sizeof(int));
	int tolCls,totCls;//two ways to count total number of clusters(to check if the code is right)
	int newCls;//the number of newly created clusters. newCls>tolCls because some clusters combine to one.
	double tempd;
	double delta;
	int upperbound;
	int temp,temp1;
	int index,index0;
	double arta;//arctan


	int freeg,freeG;//two ways to count the number of free gephyrins(to check if the code is right)
	int totalG;//total number of gephyrins=N(to check if the code is right)
	int indexb,tempp;

	double *xr,*yr;//receptors coordinates
	xr=(double*)malloc(Nr*sizeof(double));
	yr=(double*)malloc(Nr*sizeof(double));
	int *freer;//mark if a receptor joined a cluster:joined,the index of gephyrin who took it;free,0.
	freer=(int*)malloc(Nr*sizeof(int));
	int *rec1;//mark if a square head taking a receptor:taken,the index of the taken receptor;not,-1.
	rec1=(int*)malloc(N*sizeof(int));
	int freeR,freeRe;//two ways to count the total number of free receptors



	int tolRcls;//total number of receptors clusters

  int dim=2;//dimension of data
  struct kdtree *tree;
  int kn;//kn_th nearest neighbour
  int mark[N];
  int count;
  int nj;//the index of gephyrin whose square head is the nearest to the i_th gephyrin
  int nj1,nj2;//the index of gephyrin whose circle head is the nearest to the i_th gephyrin
  double (*d)[2];

	double prob;


/**
 * initialization
 */


 	//tolCls=0;
	//newCls=0;
	//freeg=N;
	//freeR=Nr;

	FILE *fr1;
	fr1 = fopen("gephyrin0", "r");
	double g[N][10];

	while(!feof(fr1)){
		for(i=0;i<N;i++){
			for(j=0;j<10;j++){
				fscanf(fr1, "%lf", &g[i][j]);
			}
		}
	}

	FILE *fr2;
	fr2 = fopen("receptor0", "r");
	double r[Nr][3];

	while(!feof(fr2)){
		for(i=0;i<Nr;i++){
			for(j=0;j<3;j++){
				fscanf(fr2, "%lf", &r[i][j]);
			}
		}
	}

	FILE *fr3;
	fr3 = fopen("cluster0", "r");
	int c[N][N+1];

	while(!feof(fr3)){
		for(i=0;i<N;i++){
			for(j=0;j<N+1;j++){
				fscanf(fr3, "%d", &c[i][j]);
			}
		}
	}

	FILE *fr4;
	fr4 = fopen("stat0", "r");
	int s[4];

	while(!feof(fr4)){
		for(j=0;j<4;j++){
			fscanf(fr4, "%d", &s[j]);
		}
	}

	////gephyrins
	double *k,*alpha;
	k=(double*)malloc(N*sizeof(double));
	alpha=(double*)malloc(N*sizeof(double));
	for(i=0;i<N;i++){
        rec1[i]=(int)g[i][4];
		alpha[i]=g[i][9];;
		k[i]=tan(alpha[i]);
		x[i]=g[i][0];
		y[i]=g[i][1];
		x1[i]=g[i][2];
		y1[i]=g[i][3];
		pnRand[i]=(int)g[i][8];
		/*
		if(((alpha[i]>=0)&&(alpha[i]<pi/2.0))||((alpha[i]>=3.0*pi/2.0)&&(alpha[i]<2.0*pi))){
			pnRand[i]=0;
		}else{
			pnRand[i]=1;
		}
		if(pnRand[i]==0){
			x1[i]=x[i]+sqrt(geLength*geLength/(k[i]*k[i]+1.0));
			y1[i]=y[i]+k[i]*sqrt(geLength*geLength/(k[i]*k[i]+1.0));
		}else{
			x1[i]=x[i]-sqrt(geLength*geLength/(k[i]*k[i]+1.0));
			y1[i]=y[i]-k[i]*sqrt(geLength*geLength/(k[i]*k[i]+1.0));
		}
    */
		node[i]=(int)g[i][5];
		node1[i]=(int)g[i][6];
		for(j=0;j<(N+1);j++){
			cls[i][j]=c[i][j];
		}
		beCls[i]=(int)g[i][7];
	}

	////receptors
	for(i=0;i<Nr;i++){
		xr[i]=r[i][0];
		yr[i]=r[i][1];
		freer[i]=(int)r[i][2];
		//clsRsize[i]=0;
		//clsRfreq[i]=0;
	}

	////Statistics
	freeg=s[0];
	tolCls=s[1];
	freeR=s[2];
	newCls=s[3];

	fclose(fr1);
	fclose(fr2);
	fclose(fr3);
	fclose(fr4);


/*
  ////dealing with collision
	xi=0.0;
	delta=0.01;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			if(j!=i){
				xi=(k[i]*x[i]-k[j]*x[j]-y[i]+y[j])/(k[i]-k[j]);
				if((xi>x[i])&&(xi<x1[i])&&(xi>x[j])&&(xi<x1[j])){
					if(node[j]==-1&&node1[j]==-1){ //case: single goes across a cluster/single goes across single
						x[j]=x[i]+1.5*(cir+squ)/2.0;
						y[j]=y[i];
						k[j]=k[i]+delta;
						////ALPHA
						if(pnRand[j]==0){
							x1[j]=x[j]+sqrt(geLength*geLength/(k[j]*k[j]+1.0));
							y1[j]=y[j]+k[j]*sqrt(geLength*geLength/(k[j]*k[j]+1.0));
						}else{
							x1[j]=x[j]-sqrt(geLength*geLength/(k[j]*k[j]+1.0));
							y1[j]=y[j]-k[j]*sqrt(geLength*geLength/(k[j]*k[j]+1.0));
						}
						arta=atan2((y1[j]-y[j]),(x1[j]-x[j]));
						if(arta<0.0){
							alpha[j]=arta+2*pi;
						}else{
							alpha[j]=arta;
						}

					}
					//if((node[i]>0||node1[i]>0)&&(node[j]>0||node1[j]>0))//case: a cluster goes across a cluster-----for now, to be ignored.
				}
			}
			xi=0.0;
			yi=0.0;
		}
	}
*/






/**
 * self-organization
 */

	FILE *f1;
	FILE *f3;
	FILE *f6;
	FILE *f8;
	FILE *f9;

	FILE *f2;
	FILE *f4;
	FILE *f5;
	FILE *f7;



	char gefile[30],refile[30],clsfile[30],statfile[30];

	sprintf(gefile,"gephyrin%d",tmax);
	sprintf(refile,"receptor%d",tmax);
	sprintf(clsfile,"cluster%d",tmax);
	sprintf(statfile,"stat%d",tmax);

	f2 = fopen(gefile,"w");
	f4 = fopen(refile,"w");
	f5 = fopen(clsfile,"w");
	f7 = fopen(statfile,"w");

	for(tau=0;tau<tmax;tau++){
    /*
		for(i=0;i<N;i++){
			for(j=0;j<(N+1);j++){
				printf("%d\t",cls[i][j]);
			}
			printf("\n");
		}
    */
    /*
		for(i=0;i<N;i++){
			if(beCls[i]>0){
				printf("before------i=%d,x=%f,y=%f,x1=%f,y1=%f,node=%d,node1=%d,beCls=%d\n",i,x[i],y[i],x1[i],y1[i],node[i],node1[i],beCls[i]);
			}
		}
    */


		char csfile[50],rcfile[50],nnfile[50],fgfile[50],frfile[50];

		if(((tau+1)%divider)==0){

			sprintf(csfile,"results/clustersStatistics%d.dat",tau);
			sprintf(rcfile,"results/receptorsClusters%d.dat",tau);
			sprintf(nnfile,"results/nearestNeighbour%d.dat",tau);
	    sprintf(fgfile,"results/gephyrinCoor%d.dat",tau);
	    sprintf(frfile,"results/receptorCoor%d.dat",tau);

			f1 = fopen(csfile,"w");
			f3 = fopen(rcfile,"w");
			f6 = fopen(nnfile,"w");
			f8 = fopen(fgfile,"w");
			f9 = fopen(frfile,"w");
    }
		  /**
		   * random walk
		   */
		  ////free receptors' random walk
		  for(i=0;i<Nr;i++){
        if(freer[i]==0){
                   u1 = (rand()+0.1)/(RAND_MAX + 1.0);//to avoid log(0)
				   u2 = rand()/(RAND_MAX + 1.0);
				   z0 = 2*sqrt(-2.0 * log(u1)*diffr) * cos(2*pi * u2);
				   z1 = 2*sqrt(-2.0 * log(u1)*diffr) * sin(2*pi * u2);
				   xr[i]+= z0;
				   yr[i]+= z1;
				   //boundary condition: sticky
				   if(xr[i]<0.0)
					 xr[i]=0.0;
				   if(xr[i]>(float)L)
					 xr[i]=(float)L;
				   if(yr[i]<0.0)
					 yr[i]=0.0;
				   if(yr[i]>(float)L)
					 yr[i]=(float)L;
        }
		  }

		  ////free gephyrins' random walk
		   for(i=0;i<N;i++){
			   if((node[i]==-1)&&(node1[i]==-1)){

				   ////update step sizes
				   u1 = (rand()+0.1)/(RAND_MAX + 1.0);
				   u2 = rand()/(RAND_MAX + 1.0);
				   z0 = 2*sqrt(-2.0 * log(u1)*diff) * cos(2*pi * u2);
				   z1 = 2*sqrt(-2.0 * log(u1)*diff) * sin(2*pi * u2);
				   x[i]+= z0;
				   y[i]+= z1;
				   alpha[i]=2.0*pi*((rand()+0.1)/(RAND_MAX+1.0));
		           k[i]=tan(alpha[i]);
				   if(((alpha[i]>=0)&&(alpha[i]<pi/2.0))||((alpha[i]>=3.0*pi/2.0)&&(alpha[i]<2.0*pi))){
						pnRand[i]=0;
				   }else{
						pnRand[i]=1;
				   }
				   if(pnRand[i]==0){
					   x1[i]=x[i]+sqrt(geLength*geLength/(k[i]*k[i]+1.0));
					   y1[i]=y[i]+k[i]*sqrt(geLength*geLength/(k[i]*k[i]+1.0));
				   }else{
					   x1[i]=x[i]-sqrt(geLength*geLength/(k[i]*k[i]+1.0));
					   y1[i]=y[i]-k[i]*sqrt(geLength*geLength/(k[i]*k[i]+1.0));
				   }


				   ////boundary condition: sticky
				   if(x[i]<0||x1[i]<0){
					   if(x[i]<x1[i]){
						   x[i]=0;
						   x1[i]=x[i]+sqrt(geLength*geLength/(k[i]*k[i]+1.0));
					   }else{
						   x1[i]=0;
						   x[i]=x1[i]+sqrt(geLength*geLength/(k[i]*k[i]+1.0));
					   }
				   }
				   if(x[i]>(float)L||x1[i]>(float)L){
					   if(x1[i]>x[i]){
						   x1[i]=(float)L;
						   x[i]=x1[i]-sqrt(geLength*geLength/(k[i]*k[i]+1.0));
					   }else{
						   x[i]=(float)L;
						   x1[i]=x[i]-sqrt(geLength*geLength/(k[i]*k[i]+1.0));
					   }
				   }
				   if(y[i]<0||y1[i]<0){
						if(y[i]<y1[i]){
						   y[i]=0;
						   if(pnRand[i]==0){
							   y1[i]=y[i]+k[i]*sqrt(geLength*geLength/(k[i]*k[i]+1.0));
						   }else{
							   y1[i]=y[i]-k[i]*sqrt(geLength*geLength/(k[i]*k[i]+1.0));
						   }
						}else{
						   y1[i]=0;
						   if(pnRand[i]==0){
							   y[i]=y1[i]-k[i]*sqrt(geLength*geLength/(k[i]*k[i]+1.0));
						   }else{
							   y[i]=y1[i]+k[i]*sqrt(geLength*geLength/(k[i]*k[i]+1.0));
						   }
						}
				   }
				   if(y[i]>(float)L||y1[i]>(float)L){
						if(y1[i]>y[i]){
						   y1[i]=(float)L;
						   if(pnRand[i]==0){
							   y[i]=y1[i]-k[i]*sqrt(geLength*geLength/(k[i]*k[i]+1.0));
						   }else{
							   y[i]=y1[i]+k[i]*sqrt(geLength*geLength/(k[i]*k[i]+1.0));
						   }
						}else{
						   y[i]=(float)L;
						   if(pnRand[i]==0){
							   y1[i]=y[i]+k[i]*sqrt(geLength*geLength/(k[i]*k[i]+1.0));
						   }else{
							   y1[i]=y[i]-k[i]*sqrt(geLength*geLength/(k[i]*k[i]+1.0));
						   }
						}
				   }
		       }
		   }

       /*
			 for(i=0;i<N;i++){
	 			if(beCls[i]>0){
	 				printf("free------i=%d,x=%f,y=%f,x1=%f,y1=%f,node=%d,node1=%d,beCls=%d\n",i,x[i],y[i],x1[i],y1[i],node[i],node1[i],beCls[i]);
	 			}
	 		}
      */
			  ////clusters' random walk
				int temprec;
				temprec=0;
			  for(i=0;i<N;i++){
				  if(cls[i][0]>0){
					  u1 = (rand()+0.1)/(RAND_MAX + 1.0);
					  u2 = rand()/(RAND_MAX + 1.0);
					  z0 = 2*sqrt(-2.0 * log(u1)*diff) * cos(2*pi * u2);
					  z1 = 2*sqrt(-2.0 * log(u1)*diff) * sin(2*pi * u2);
					  for(j=1;j<=cls[i][0];j++){
						   x[cls[i][j]]+= z0;
						   y[cls[i][j]]+= z1;
						   if(pnRand[cls[i][j]]==0){
							   x1[cls[i][j]]=x[cls[i][j]]+sqrt(geLength*geLength/(k[cls[i][j]]*k[cls[i][j]]+1.0));
							   y1[cls[i][j]]=y[cls[i][j]]+k[cls[i][j]]*sqrt(geLength*geLength/(k[cls[i][j]]*k[cls[i][j]]+1.0));
						   }else{
							   x1[cls[i][j]]=x[cls[i][j]]-sqrt(geLength*geLength/(k[cls[i][j]]*k[cls[i][j]]+1.0));
							   y1[cls[i][j]]=y[cls[i][j]]-k[cls[i][j]]*sqrt(geLength*geLength/(k[cls[i][j]]*k[cls[i][j]]+1.0));
						   }
							 tempd=0.0;

						 	 if(x[cls[i][j]]<0){
						 		 tempd=-x[cls[i][j]];
						 		 for(l=1;l<=cls[i][0];l++){
						 			 x[cls[i][l]]+=tempd;
						 			 x1[cls[i][l]]+=tempd;
						 		 }
						 		 tempd=0.0;
						 	 }
						 	 if(x1[cls[i][j]]<0){
						 		 tempd=-x1[cls[i][j]];
						 		 for(l=1;l<=cls[i][0];l++){
						 			 x[cls[i][l]]+=tempd;
						 			 x1[cls[i][l]]+=tempd;
						 		 }
						 		 tempd=0.0;
						 	 }
						 	 if(x[cls[i][j]]>(float)L){
						 		 tempd=x[cls[i][j]]-(float)L;
						 		 for(l=1;l<=cls[i][0];l++){
						 			 x[cls[i][l]]-=tempd;
						 			 x1[cls[i][l]]-=tempd;
						 		 }
						 		 tempd=0.0;
						 	 }
						 	 if(x1[cls[i][j]]>(float)L){
						 		 tempd=x1[cls[i][j]]-(float)L;
						 		 for(l=1;l<=cls[i][0];l++){
						 			 x[cls[i][l]]-=tempd;
						 			 x1[cls[i][l]]-=tempd;
						 		 }
						 		 tempd=0.0;
						 	 }
						 	 if(y[cls[i][j]]<0){
						 		 tempd=-y[cls[i][j]];
						 		 for(l=1;l<=cls[i][0];l++){
						 			 y[cls[i][l]]+=tempd;
						 			 y1[cls[i][l]]+=tempd;
						 		 }
						 		 tempd=0.0;
						 	 }
						 	 if(y1[cls[i][j]]<0){
						 		 tempd=-y1[cls[i][j]];
						 		 for(l=1;l<=cls[i][0];l++){
						 			 y[cls[i][l]]+=tempd;
						 			 y1[cls[i][l]]+=tempd;
						 		 }
						 		 tempd=0.0;
						 	 }
						 	 if(y[cls[i][j]]>(float)L){
						 		 tempd=y[cls[i][j]]-(float)L;
						 		 for(l=1;l<=cls[i][0];l++){
						 			 y[cls[i][l]]-=tempd;
						 			 y1[cls[i][l]]-=tempd;
						 		 }
						 		 tempd=0.0;
						 	 }
						 	 if(y1[cls[i][j]]>(float)L){
						 		 tempd=y1[cls[i][j]]-(float)L;
						 		 for(l=1;l<=cls[i][0];l++){
						 			 y[cls[i][l]]-=tempd;
						 			 y1[cls[i][l]]-=tempd;
						 		 }
						 		 tempd=0.0;
						 	 }
						 	 tempd=0.0;


							 if(rec1[cls[i][j]]>=0){
								 temprec=rec1[cls[i][j]];
								 xr[temprec]=x1[cls[i][j]];
								 yr[temprec]=y1[cls[i][j]];
							 }
					  }

				  }

			  }
        /*
				for(i=0;i<N;i++){
					if(beCls[i]>0){
						printf("cluster------i=%d,x=%f,y=%f,x1=%f,y1=%f,node=%d,node1=%d,beCls=%d\n",i,x[i],y[i],x1[i],y1[i],node[i],node1[i],beCls[i]);
					}
				}
				*/
			    /**
				 * dealing with collision: if go across, one move a little bit, almost parallel to the other.
				 * single goes cross single
				 * single goes across a cluster
				 * a cluster goes across a cluster
				 */
				xi=0.0;
				delta=0.001;
				for(i=0;i<N;i++){
					for(j=0;j<N;j++){
						if(j!=i){
							xi=(k[i]*x[i]-k[j]*x[j]-y[i]+y[j])/(k[i]-k[j]);
							if((xi>x[i])&&(xi<x1[i])&&(xi>x[j])&&(xi<x1[j])){
								if(node[j]==-1&&node1[j]==-1){ //case: single goes across a cluster/single goes across single
									x[j]=x[i]+1.5*(cir+squ)/2.0;
									y[j]=y[i];
									k[j]=k[i]+delta;
									if(pnRand[j]==0){
										x1[j]=x[j]+sqrt(geLength*geLength/(k[j]*k[j]+1.0));
										y1[j]=y[j]+k[j]*sqrt(geLength*geLength/(k[j]*k[j]+1.0));
									}else{
										x1[j]=x[j]-sqrt(geLength*geLength/(k[j]*k[j]+1.0));
										y1[j]=y[j]-k[j]*sqrt(geLength*geLength/(k[j]*k[j]+1.0));
									}
									arta=atan2((y1[j]-y[j]),(x1[j]-x[j]));
									if(arta<0.0){
										alpha[j]=arta+2*pi;
									}else{
										alpha[j]=arta;
									}
								}
								//if((node[i]>0||node1[i]>0)&&(node[j]>0||node1[j]>0))//case: a cluster goes across a cluster-----for now, to be ignored.
							}
						}
						xi=0.0;
						yi=0.0;
					}
				}

        /*
				for(i=0;i<N;i++){
					if(beCls[i]>0){
						printf("collision------i=%d,x=%f,y=%f,x1=%f,y1=%f,node=%d,node1=%d,beCls=%d\n",i,x[i],y[i],x1[i],y1[i],node[i],node1[i],beCls[i]);
					}
				}
        */

		   /**
		    * Clustering
		    */
		      /**
		       * 3 gephyrins combine
		       */
		   double d1,d2,d3,r,rrr;
		   double d11,d12,d13;
		   d1=(float)L;
		   d2=(float)L;
		   d3=(float)L;
		   d11=(float)L;
		   d12=(float)L;
		   d13=(float)L;

		   r=cir;
		   rrr=1.5*r;
		   double epsilon;
		   epsilon=cir*0.5;
		   double alpha0;
		   alpha0=0.3*pi;

		   double alpha1;
		   alpha1=pi*0.9;

		   double distance1,distance2;
		   distance1=(float)L;
		   distance2=(float)L;
		   kn=3;
		   count=0;
		   tree = kdtree_init(dim);
		   if(tree==NULL){
			  exit(-1);
		   }
		   for(i=0;i<N;i++){
			   double circle[]={x[i],y[i]};
			   kdtree_insert(tree, circle);
			   count++;
		   }
		   kdtree_rebuild(tree);
		   kdtree_dump(tree);

		   for(i=0;i<N;i++){

				if(node[i]==-1){

						   double targetG[2] = {x[i],y[i]};
						   kdtree_knn_search(tree, targetG, kn);
						   d=kdtree_knn_dump(tree);
						   distance1=d[1][0];
						   distance2=d[2][0];
						   j=(int)d[1][1];
						   l=(int)d[2][1];
						   knn_clear(tree);
						   /*
						    * 3 circle heads combine
						   */
						   if((distance1<=rrr)&&(distance2<=rrr)){
							        ///1+1+1=3gephyrins
									if((node[j]==-1)&&(node[l]==-1)){
										   Gx=(x[i]+x[j]+x[l])/3.0;
										   Gy=(y[i]+y[j]+y[l])/3.0;
										   d1=(x[i]-Gx)*(x[i]-Gx)+(y[i]-Gy)*(y[i]-Gy);
										   d1=sqrt(d1);
										   d2=(x[j]-Gx)*(x[j]-Gx)+(y[j]-Gy)*(y[j]-Gy);
										   d2=sqrt(d2);
										   d3=(x[l]-Gx)*(x[l]-Gx)+(y[l]-Gy)*(y[l]-Gy);
										   d3=sqrt(d3);
										   d11=(x1[i]-Gx)*(x1[i]-Gx)+(y1[i]-Gy)*(y1[i]-Gy);
										   d11=sqrt(d11);
										   d12=(x1[j]-Gx)*(x1[j]-Gx)+(y1[j]-Gy)*(y1[j]-Gy);
										   d12=sqrt(d12);
										   d13=(x1[l]-Gx)*(x1[l]-Gx)+(y1[l]-Gy)*(y1[l]-Gy);
										   d13=sqrt(d13);
										   if((fabs(d1-d2)<epsilon)&&(fabs(d2-d3)<epsilon)&&(fabs(d3-d1)<epsilon)&&(fabs(d11-d12)<epsilon)&&(fabs(d12-d13)<epsilon)&&(fabs(d13-d11)<epsilon)&&(fabs(alpha[i]-alpha[j])>alpha0)&&(fabs(alpha[j]-alpha[l])>alpha0)&&(fabs(alpha[l]-alpha[i])>alpha0)&&(fabs(alpha[i]-alpha[j])<(2*pi-alpha0))&&(fabs(alpha[j]-alpha[l])<(2*pi-alpha0))&&(fabs(alpha[l]-alpha[i])<(2*pi-alpha0))){

												 if((beCls[i]==0&&beCls[j]==0)||(beCls[i]==0&&beCls[l]==0)||(beCls[j]==0&&beCls[l]==0)||(beCls[i]!=0&&beCls[j]!=0&&beCls[i]!=beCls[j]&&beCls[i]!=beCls[l]&&beCls[j]!=beCls[l])||(beCls[i]!=0&&beCls[l]!=0&&beCls[i]!=beCls[l]&&beCls[i]!=beCls[j]&&beCls[l]!=beCls[j])||(beCls[j]!=0&&beCls[l]!=0&&beCls[j]!=beCls[l]&&beCls[j]!=beCls[i]&&beCls[l]!=beCls[i])){


											 if(sil==1){
											 printf("333333gephyrin%f %f %f %f %f %f %d %d %d i%d j%d l%d,%d %d %d\n",d1,d2,d3,fabs(d1-d2),fabs(d2-d3),fabs(d3-d1),node1[i],node1[j],node1[l],beCls[i],beCls[j],beCls[l],i,j,l);
											 }

											 //fprintf(f4,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",x[i],y[i],x1[i],y1[i],x[j],y[j],x1[j],y1[j],x[l],y[l],x1[l],y1[l]);


											 ///update clusters information
												 if((node1[i]==-1)&&(node1[j]==-1)&&(node1[l]==-1)){
													 tolCls++;
													 newCls++;
													 freeg=freeg-3;

													 for(index0=1;index0<N;index0++){
														 if(cls[index0][0]==-1){
															 break;
														 }
													 }
													 beCls[i]=index0;
													 beCls[j]=index0;
													 beCls[l]=index0;
													 cls[index0][0]=1;//cls[][0] count how many gephyrins added in this cluster
													 cls[index0][cls[index0][0]]=i;//cls[][>0] record which gephyrin just added in this cluster
													 cls[index0][0]++;
													 cls[index0][cls[index0][0]]=j;
													 cls[index0][0]++;
													 cls[index0][cls[index0][0]]=l;

												 }else if((node1[i]>=0)&&(node1[j]==-1)&&(node1[l]==-1)){
													 freeg=freeg-2;
													 beCls[j]=beCls[i];
													 beCls[l]=beCls[i];
													 cls[beCls[i]][0]++;
													 cls[beCls[i]][cls[beCls[i]][0]]=j;
													 cls[beCls[i]][0]++;
													 cls[beCls[i]][cls[beCls[i]][0]]=l;
												 }else if((node1[i]==-1)&&(node1[j]>=0)&&(node1[l]==-1)){
													 freeg=freeg-2;
													 beCls[i]=beCls[j];
													 beCls[l]=beCls[j];
													 cls[beCls[j]][0]++;
													 cls[beCls[j]][cls[beCls[j]][0]]=i;
													 cls[beCls[j]][0]++;
													 cls[beCls[j]][cls[beCls[j]][0]]=l;
												 }else if((node1[i]==-1)&&(node1[j]==-1)&&(node1[l]>=0)){
													 freeg=freeg-2;
													 beCls[i]=beCls[l];
													 beCls[j]=beCls[l];
													 cls[beCls[l]][0]++;
													 cls[beCls[l]][cls[beCls[l]][0]]=i;
													 cls[beCls[l]][0]++;
													 cls[beCls[l]][cls[beCls[l]][0]]=j;
												 }else if((node1[i]==-1)&&(node1[j]>=0)&&(node1[l]>=0)){
													 freeg=freeg-1;
													 temp=beCls[j];
													 tolCls--;
													 cls[beCls[j]][0]++;
													 cls[beCls[j]][cls[beCls[j]][0]]=i;
													 upperbound=cls[beCls[l]][0];
													 for(p=1;p<=upperbound;p++){
														 cls[beCls[j]][0]++;
														 cls[beCls[j]][cls[beCls[j]][0]]=cls[beCls[l]][p];
													 }

													 temp=cls[beCls[j]][0];
													 temp=beCls[j];
													 beCls[i]=temp;
													 temp=beCls[l];
													 for(m=1;m<=upperbound;m++){
														beCls[cls[temp][m]]=beCls[j];
													 }
													 for(q=0;q<=upperbound;q++){
														 cls[temp][q]=-1;
													 }

													 tempp=beCls[j];
													 if(sil==1){
													 printf("%d contains:\t",tempp);
													 for(index=0;index<=N;index++){
														 printf("%d\t",cls[tempp][index]);
													 }
													 printf("\n");
													 }
													 tempp=beCls[l];
													 if(sil==1){
													 printf("%d contains:\t",tempp);
													 for(index=0;index<=N;index++){
														 printf("%d\t",cls[tempp][index]);
													 }
													 printf("\n");
													 }
												 }else if((node1[i]>=0)&&(node1[j]>=0)&&(node1[l]==-1)){
													 tolCls--;
													 freeg=freeg-1;
													 cls[beCls[i]][0]++;
													 cls[beCls[i]][cls[beCls[i]][0]]=l;
													 upperbound=cls[beCls[j]][0];
													 for(p=1;p<=upperbound;p++){
														 cls[beCls[i]][0]++;
														 cls[beCls[i]][cls[beCls[i]][0]]=cls[beCls[j]][p];
													 }
													 beCls[l]=beCls[i];
													 temp=beCls[j];
													 for(m=1;m<=upperbound;m++){
														beCls[cls[temp][m]]=beCls[i];
													 }
													 for(q=0;q<=upperbound;q++){
														 cls[temp][q]=-1;
													 }
													 tempp=beCls[i];
													 if(sil==1){
													 printf("%d contains:\t",tempp);
													 for(index=0;index<=N;index++){
														 printf("%d\t",cls[tempp][index]);
													 }
													 printf("\n");
													 }
													 tempp=beCls[j];
													 if(sil==1){
													 printf("%d contains:\t",tempp);
													 for(index=0;index<=N;index++){
														 printf("%d\t",cls[tempp][index]);
													 }
													 printf("\n");
													 }

												 }else if((node1[i]>=0)&&(node1[j]==-1)&&(node1[l]>=0)){
													 tolCls--;
													 freeg=freeg-1;
													 cls[beCls[i]][0]++;
													 cls[beCls[i]][cls[beCls[i]][0]]=j;
													 upperbound=cls[beCls[l]][0];
													 for(p=1;p<=upperbound;p++){
														 cls[beCls[i]][0]++;
														 cls[beCls[i]][cls[beCls[i]][0]]=cls[beCls[l]][p];
													 }
													 beCls[j]=beCls[i];
													 temp=beCls[l];
													 for(m=1;m<=upperbound;m++){
														beCls[cls[temp][m]]=beCls[i];
													 }
													 for(q=0;q<=upperbound;q++){
														 cls[temp][q]=-1;
													 }
													 tempp=beCls[i];
													 if(sil==1){
													 printf("%d contains:\t",tempp);
													 for(index=0;index<=N;index++){
														 printf("%d\t",cls[tempp][index]);
													 }
													 printf("\n");
													 }
													 tempp=beCls[l];
													 if(sil==1){
													 printf("%d contains:\t",tempp);
													 for(index=0;index<=N;index++){
														 printf("%d\t",cls[tempp][index]);
													 }
													 printf("\n");
													 }


												 }else if((node1[i]>=0)&&(node1[j]>=0)&&(node1[l]>=0)){
													 tolCls-=2;
													 temp=beCls[j];
													 temp1=beCls[l];
													 upperbound=cls[temp][0];

													 for(p=1;p<=upperbound;p++){
														 cls[beCls[i]][0]++;
														 cls[beCls[i]][cls[beCls[i]][0]]=cls[temp][p];
													 }
													 upperbound=cls[temp1][0];
													 for(p=1;p<=upperbound;p++){
														 cls[beCls[i]][0]++;
														 cls[beCls[i]][cls[beCls[i]][0]]=cls[temp1][p];
													 }
													 upperbound=cls[temp][0];
													 for(m=1;m<=upperbound;m++){
														beCls[cls[temp][m]]=beCls[i];
													 }
													 upperbound=cls[temp1][0];
													 for(m=1;m<=upperbound;m++){
														beCls[cls[temp1][m]]=beCls[i];
													 }
													 upperbound=cls[temp][0];
													 for(q=0;q<=upperbound;q++){
														 cls[temp][q]=-1;
														 if(sil==1){
														 printf("cls%d=%d\n",temp,cls[temp][q]);
														 }
													 }
													 upperbound=cls[temp1][0];
													 for(q=0;q<=upperbound;q++){
														 cls[temp1][q]=-1;
														 if(sil==1){
														 printf("cls%d=%d\n",temp1,cls[temp1][q]);
														 }
													 }
													 tempp=beCls[i];
													 if(sil==1){
													 printf("%d contains:\t",tempp);
													 for(index=0;index<=N;index++){
														 printf("%d\t",cls[tempp][index]);
													 }
													 printf("\n");
													 }

													 tempp=beCls[j];
													 if(sil==1){
													 printf("%d contains:\t",tempp);
													 for(index=0;index<=N;index++){
														 printf("%d\t",cls[tempp][index]);
													 }
													 printf("\n");
													 }
													 tempp=beCls[l];
													 if(sil==1){
													 printf("%d contains:\t",tempp);
													 for(index=0;index<=N;index++){
														 printf("%d\t",cls[tempp][index]);
													 }
													 printf("\n");
													 }
												 }

												 node[i]=j;
												 node[j]=l;
												 node[l]=i;
											  }

										   }
										   Gx=0.0;
										   Gy=0.0;
										   d1=(float)L;
										   d2=(float)L;
										   d3=(float)L;
										   d11=(float)L;
										   d12=(float)L;
										   d13=(float)L;

									}
									///2+1=3gephyrins
									if((node[j]==l)&&(node[l]==j)){
										   Gx=(x[i]+x[j]+x[l])/3.0;
										   Gy=(y[i]+y[j]+y[l])/3.0;
										   d1=(x[i]-Gx)*(x[i]-Gx)+(y[i]-Gy)*(y[i]-Gy);
										   d1=sqrt(d1);
										   d2=(x[j]-Gx)*(x[j]-Gx)+(y[j]-Gy)*(y[j]-Gy);
										   d2=sqrt(d2);
										   d3=(x[l]-Gx)*(x[l]-Gx)+(y[l]-Gy)*(y[l]-Gy);
										   d3=sqrt(d3);
										   d11=(x1[i]-Gx)*(x1[i]-Gx)+(y1[i]-Gy)*(y1[i]-Gy);
										   d11=sqrt(d11);
										   d12=(x1[j]-Gx)*(x1[j]-Gx)+(y1[j]-Gy)*(y1[j]-Gy);
										   d12=sqrt(d12);
										   d13=(x1[l]-Gx)*(x1[l]-Gx)+(y1[l]-Gy)*(y1[l]-Gy);
										   d13=sqrt(d13);
										   if((fabs(d1-d2)<epsilon)&&(fabs(d2-d3)<epsilon)&&(fabs(d3-d1)<epsilon)&&(fabs(d11-d12)<epsilon)&&(fabs(d12-d13)<epsilon)&&(fabs(d13-d11)<epsilon)&&(fabs(alpha[i]-alpha[j])>alpha0)&&(fabs(alpha[j]-alpha[l])>alpha0)&&(fabs(alpha[l]-alpha[i])>alpha0)&&(fabs(alpha[i]-alpha[j])<(2*pi-alpha0))&&(fabs(alpha[j]-alpha[l])<(2*pi-alpha0))&&(fabs(alpha[l]-alpha[i])<(2*pi-alpha0))){

												 if(beCls[i]!=beCls[j]){


											 if(sil==1){
											 printf("2+1gephyrin%f %f %f %f %f %f %d %d %d i%d j%d l%d,%d %d %d\n",d1,d2,d3,fabs(d1-d2),fabs(d2-d3),fabs(d3-d1),node1[i],node1[j],node1[l],beCls[i],beCls[j],beCls[l],i,j,l);
											 }

											 //fprintf(f4,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",x[i],y[i],x1[i],y1[i],x[j],y[j],x1[j],y1[j],x[l],y[l],x1[l],y1[l]);


											 ///update clusters information
                                             if(node1[i]==-1){
												 freeg=freeg-1;
												 beCls[i]=beCls[j];
												 cls[beCls[j]][0]++;
												 cls[beCls[j]][cls[beCls[j]][0]]=i;
											 }
											 if(node1[i]>=0){
												 tolCls--;
												 upperbound=cls[beCls[i]][0];
												 for(p=1;p<=upperbound;p++){
													 cls[beCls[j]][0]++;
													 cls[beCls[j]][cls[beCls[j]][0]]=cls[beCls[i]][p];
												 }
												 temp=beCls[i];
												 for(m=1;m<=upperbound;m++){
													beCls[cls[temp][m]]=beCls[j];
												 }
												 for(q=0;q<=upperbound;q++){
													 cls[temp][q]=-1;
												 }
												 tempp=beCls[j];
												 if(sil==1){
												 printf("%d contains:\t",tempp);
												 for(index=0;index<=N;index++){
													 printf("%d\t",cls[tempp][index]);
												 }
												 printf("\n");
												 }
												 tempp=beCls[i];
												 if(sil==1){
												 printf("%d contains:\t",tempp);
												 for(index=0;index<=N;index++){
													 printf("%d\t",cls[tempp][index]);
												 }
												 printf("\n");
												 }
											 }
												 node[i]=j;
												 node[j]=l;
												 node[l]=i;
											  }

										   }
										   Gx=0.0;
										   Gy=0.0;
										   d1=(float)L;
										   d2=(float)L;
										   d3=(float)L;
										   d11=(float)L;
										   d12=(float)L;
										   d13=(float)L;

									}
						   }
						   /*
						    * 2 circle heads combine
						   */
						   if((distance1<=rrr)&&(distance2>rrr)){
							 prob=(float)(rand()/(RAND_MAX + 1.0));
							 if(prob>(1-Pc)){
							  if(node[j]==-1){
								if((fabs(alpha[i]-alpha[j])>alpha1)&&(fabs(alpha[i]-alpha[j])<=pi)){
									if((beCls[i]==0||beCls[j]==0)||((beCls[i]!=0)&&(beCls[i]!=beCls[j]))){
										   if(sil==1){
										   printf("!!!22222222222gephy Circle Heads   %f %f ai=%f aj=%f i%d j%d, %d %d\n",distance1,(fabs(alpha[i]-alpha[j])),alpha[i]*180.0/pi,alpha[j]*180.0/pi,beCls[i],beCls[j],i,j);
										   }

								   //fprintf(f5,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",x[i],y[i],x1[i],y1[i],x[j],y[j],x1[j],y1[j]);


								   if((node1[i]==-1)&&(node1[j]==-1)){
									   tolCls++;
									   freeg=freeg-2;
									   newCls++;

									   for(index0=1;index0<N;index0++){
										   if(cls[index0][0]==-1){
											   break;
										   }
									   }
									   beCls[i]=index0;
									   beCls[j]=index0;
									   cls[index0][0]=1;
									   cls[index0][cls[index0][0]]=i;
									   cls[index0][0]++;
									   cls[index0][cls[index0][0]]=j;

								   }else if((node1[i]>=0)&&(node1[j]==-1)){
									   freeg=freeg-1;
									   beCls[j]=beCls[i];
									   cls[beCls[i]][0]++;
									   cls[beCls[i]][cls[beCls[i]][0]]=j;

								   }else if((node1[i]==-1)&&(node1[j]>=0)){
									   freeg=freeg-1;
									   beCls[i]=beCls[j];
									   cls[beCls[j]][0]++;
									   cls[beCls[j]][cls[beCls[j]][0]]=i;
								   }else if((node1[i]>=0)&&(node1[j]>=0)){
									   tolCls--;
									   temp=beCls[j];
									   upperbound=cls[temp][0];
									   for(p=1;p<=upperbound;p++){
										   cls[beCls[i]][0]++;
										   cls[beCls[i]][cls[beCls[i]][0]]=cls[temp][p];
									   }

									   for(m=1;m<=upperbound;m++){
										   beCls[cls[temp][m]]=beCls[i];
									   }
									   for(q=0;q<=upperbound;q++){
										   cls[temp][q]=-1;

									   }
									   tempp=beCls[i];
									   if(sil==1){
										 printf("%d contains:\t",tempp);
										 for(index=0;index<=N;index++){
											 printf("%d\t",cls[tempp][index]);
										 }
										 printf("\n");
									   }
									   tempp=beCls[j];
									   if(sil==1){
										 printf("%d contains:\t",tempp);
										 for(index=0;index<=N;index++){
											 printf("%d\t",cls[tempp][index]);
										 }
										 printf("\n");
									   }
								   }
								   node[i]=j;
								   node[j]=i;
									 }

								   }
						       }
						     }
						   }
				}
		   }

		   kdtree_destroy(tree);


		        /**
		         * 2 gephyrins combine
                */

		   double r1;
		   r1=squ;
		   double dgr;
		   dgr=(float)L;

		   count=0;
		   kn=2;
		   double distance;
		   distance=(float)L;

		   tree = kdtree_init(dim);
		   if(tree==NULL){
			  exit(-1);
		   }
		   for(i=0;i<N;i++){
				   double square[]={x1[i],y1[i]};
				   kdtree_insert(tree, square);
				   count++;
		   }
		   kdtree_rebuild(tree);
		   kdtree_dump(tree);


		   for(i=0;i<N;i++){
			   if(node1[i]==-1){

				   double target[2] = {x1[i],y1[i]};
		           kdtree_knn_search(tree, target, kn);

                   d=kdtree_knn_dump(tree);
				   distance=d[1][0];
				   j=(int)d[1][1];
		           knn_clear(tree);

				   if(node1[j]==-1){
					   if((distance<=r1)&&(fabs(alpha[i]-alpha[j])>alpha1)&&(fabs(alpha[i]-alpha[j])<=pi)){
							if((beCls[i]==0||beCls[j]==0)||((beCls[i]!=0)&&(beCls[i]!=beCls[j]))){
								   if(sil==1){
						           printf("22222222222gephy   %f %f ai=%f aj=%f i%d j%d, %d %d\n",distance,(fabs(alpha[i]-alpha[j])),alpha[i]*180.0/pi,alpha[j]*180.0/pi,beCls[i],beCls[j],i,j);
						           }

						   //fprintf(f5,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",x[i],y[i],x1[i],y1[i],x[j],y[j],x1[j],y1[j]);


						   if((node[i]==-1)&&(node[j]==-1)){
							   tolCls++;
							   freeg=freeg-2;
							   newCls++;

							   for(index0=1;index0<N;index0++){
								   if(cls[index0][0]==-1){
									   break;
								   }
							   }
							   beCls[i]=index0;
							   beCls[j]=index0;
							   cls[index0][0]=1;
							   cls[index0][cls[index0][0]]=i;
							   cls[index0][0]++;
							   cls[index0][cls[index0][0]]=j;

						   }else if((node[i]>=0)&&(node[j]==-1)){
							   freeg=freeg-1;
							   beCls[j]=beCls[i];
							   cls[beCls[i]][0]++;
							   cls[beCls[i]][cls[beCls[i]][0]]=j;

						   }else if((node[i]==-1)&&(node[j]>=0)){
							   freeg=freeg-1;
							   beCls[i]=beCls[j];
							   cls[beCls[j]][0]++;
							   cls[beCls[j]][cls[beCls[j]][0]]=i;
						   }else if((node[i]>=0)&&(node[j]>=0)){
							   tolCls--;
							   temp=beCls[j];
							   upperbound=cls[temp][0];
							   for(p=1;p<=upperbound;p++){
								   cls[beCls[i]][0]++;
								   cls[beCls[i]][cls[beCls[i]][0]]=cls[temp][p];
							   }

							   for(m=1;m<=upperbound;m++){
								   beCls[cls[temp][m]]=beCls[i];
							   }
							   for(q=0;q<=upperbound;q++){
							   	   cls[temp][q]=-1;
							   }
							   tempp=beCls[i];
							   if(sil==1){
								 printf("%d contains:\t",tempp);
								 for(index=0;index<=N;index++){
									 printf("%d\t",cls[tempp][index]);
								 }
								 printf("\n");
							   }
							   tempp=beCls[j];
							   if(sil==1){
								 printf("%d contains:\t",tempp);
								 for(index=0;index<=N;index++){
									 printf("%d\t",cls[tempp][index]);
								 }
								 printf("\n");
							   }
						   }
					       node1[i]=j;
					       node1[j]=i;
							 }

						   }

			       distance=(float)L;

				   }

			   }
		   }
		   ////receptors' binding
		   for(l=0;l<Nr;l++){

				double targetR[2] = {xr[l],yr[l]};
				kdtree_knn_search(tree, targetR, kn);

                d=kdtree_knn_dump(tree);
				dgr=d[0][0];
				j=(int)d[0][1];
				knn_clear(tree);
				i=node1[j];
				if((dgr<=distgr)&&(freer[l]==0)&&(i>=0)&&(rec1[j]==-1)&&(rec1[i]==-1)){
    					rec1[i]=l;
						rec1[j]=l;
						if(i==0){
							freer[l]=j;
						}else{
							freer[l]=i;
						}
						freeR--;
						if(sil==1){
						   printf("receptor%d binds %d and %d\n",l,i,j);
						}
				}

		   }
		   ////

		   kdtree_destroy(tree);




		    /**
		     *receptors' dropping off
		     */
		  int tempgr,tempgrd;//indices of dual squares which took a receptor
            for(i=0;i<Nr;i++){
                if(freer[i]>0){
                    prob=(float)(rand()/(RAND_MAX + 1.0));
                    if(prob>(1.0-Pr)){
                        tempgr=freer[i];
                        tempgrd=node1[tempgr];
						if(sil==1){
                        printf("receptor%d dropped from gephyrin%d and gephyrin%d!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",i,tempgr,tempgrd);
						}
                        freer[i]=0;
                        rec1[tempgr]=-1;
                        rec1[tempgrd]=-1;
                        freeR++;
                    }
                }
            }
			/**
			 * Clusters splitting
			 */

			int tempSize;
			tempSize=0;
			int ge;
            int tempG;
			int terB,terBi;//the gephyrin which a terminal gephyrin binded, its index in its cluster
			int terBB;//the gephyrin that the gephyrin which a terminal gephyrin binded binded
			for(i=0;i<N;i++){
				tempSize=cls[i][0];
				if(tempSize>3){
				   for(j=1;j<=tempSize;j++){
					   ge=cls[i][j];

					   ////Terminal dimers off-binding
					   if((node1[ge]>=0)&&(node[ge]==-1)&&(rec1[ge]==-1)){
						   prob=(float)(rand()/(RAND_MAX + 1.0));
						   if(prob>(1.0-Pd)){
							   terB=node1[ge];
							   node1[terB]=-1;
							   for(l=j;l<tempSize;l++){
								   cls[i][l]=cls[i][l+1];
							   }
							   cls[i][tempSize]=-1;
							   tempSize--;
							   cls[i][0]--;
							   node1[ge]=-1;
							   beCls[ge]=0;
							   freeg++;
							   if(sil==1||sil==2){
								printf("--------------------------------------------------------------------------------------------------------cls%d's %d---dimer dropped!!!!! tempSize=%d\n",i,ge,tempSize);
							   }
							break;
						   }
					   }
				       
					   
					   ////Terminal trimers off-binding
					   if((node[ge]>=0)&&(node1[ge]==-1)){
					       prob=(float)(rand()/(RAND_MAX + 1.0));
						   if(prob>(1.0-Pt)){
						       terB=node[ge];

							   ///2-1
							   if(node[terB]==ge){
							       node[terB]=-1;
								   for(l=j;l<tempSize;l++){
		                               cls[i][l]=cls[i][l+1];
							       }
								   cls[i][tempSize]=-1;
							       tempSize--;
							       cls[i][0]--;
								   node[ge]=-1;
								   beCls[ge]=0;
								   freeg++;
								   if(sil==1||sil==2){
								      printf("--------------------------------------------------------------------------------------------------------cls%d's %d---trimer 2-1 dropped!!!!! tempSize=%d\n",i,ge,tempSize);
							       }
								   break;
							   }

							   ////3-1
							   if(node[terB]!=ge){

								   if(node1[terB]==-1){
									   node[node[terB]]=terB;
									   for(l=j;l<tempSize;l++){
										  cls[i][l]=cls[i][l+1];
									   }
									   cls[i][tempSize]=-1;
									   tempSize--;
									   cls[i][0]--;
									   node[ge]=-1;
									   beCls[ge]=0;
									   freeg++;
									   if(sil==1||sil==2){
								          printf("--------------------------------------------------------------------------------------------------------cls%d's %d---trimer 3-2 dropped!!!!! tempSize=%d\n",i,ge,tempSize);
							           }
									   break;
								   }

								   if(node1[terB]!=-1){

									   if(node1[node[terB]]!=-1){
										   node[node[terB]]=terB;
										   for(l=j;l<tempSize;l++){
		                                      cls[i][l]=cls[i][l+1];
							               }
										   cls[i][tempSize]=-1;
										   tempSize--;
										   cls[i][0]--;
										   node[ge]=-1;
										   beCls[ge]=0;
										   freeg++;
										   if(sil==1||sil==2){
								               printf("--------------------------------------------------------------------------------------------------------cls%d's %d---trimer 3-1 dropped!!!!! tempSize=%d\n",i,ge,tempSize);
							               }
										   break;
									   }

									   if(node1[node[terB]]==-1){
										   node[node[terB]]=terB;
										   for(l=j;l<tempSize;l++){
											  cls[i][l]=cls[i][l+1];
										   }
										   cls[i][tempSize]=-1;
										   tempSize--;
										   cls[i][0]--;
										   node[ge]=-1;
										   beCls[ge]=0;
										   freeg++;
										   if(sil==1||sil==2){
											  printf("--------------------------------------------------------------------------------------------------------cls%d's %d---trimer 3-2b dropped!!!!! tempSize=%d\n",i,ge,tempSize);
										   }
                                           break;
									   }

								   }

							   }

						   }
					   }

				   }

				}

			}

	   if(((tau+1)%divider)==0){
	   /**
		* Statistics
		*/
		////gephyrins clusters
		int *clsSize,*clsFreq;//cluster size and corresponding frequency
	    clsSize=(int*)malloc(N*sizeof(int));
	    clsFreq=(int*)malloc(N*sizeof(int));
		int *clsOrder;//cluster size in increasing order
	    clsOrder=(int*)malloc(N*sizeof(int));
		for(i=0;i<N;i++){
		   clsSize[i]=0;
           clsFreq[i]=0;
		   clsOrder[i]=cls[i][0];
		   upperbound=cls[i][0]+1;
		   if(sil==1){
			 printf("cluster%d\t",i);
			 for(j=0;j<upperbound;j++){
				   printf("%d\t",cls[i][j]);
			   }
			 printf("\n");
		   }
		}




		freeG=0;
		for(i=0;i<N;i++){
		   if((node[i]==-1)&&(node1[i]==-1)){
			   freeG++;
			   if(sil==1){
			   printf("freeG=%d\t",i);
			   }
		   }
		}
		if(sil==1){
		printf("\n free gephyrins:%d,%d\n",freeg,freeG);
		printf("clsOrder before: ");
		for(i=0;i<N;i++){
		   printf("%d\t",clsOrder[i]);
		}
		printf("\n");
		}
		bubbleSort(clsOrder, N);
		if(sil==1){
		printf("clsOrder after: ");
		for(i=0;i<N;i++){
		   printf("%d\t",clsOrder[i]);
		}
		printf("\n");
		}



		int loc;
		loc=0;
		totCls=0;
	    totalG=0;
		for(i=0;i<N;i++){
		   j=i+1;

		   if(clsOrder[j]!=clsOrder[i]){
			   clsSize[i]=clsOrder[i];
			   clsFreq[i]=j-loc;
			   loc=j;
		   }
		   if(clsOrder[i]>0){
			   totCls++;
		   }
		   if(clsOrder[i]<0){
			   break;
		   }
		}
		for(i=0;i<N;i++){
		   if(clsSize[i]==0){
			   clsSize[i]=1;
			   clsFreq[i]=freeg;
			   break;
		   }
		}


		for(i=0;i<N;i++){
		   if(clsSize[i]>0){
			   if(((tau+1)%divider)==0){
			   fprintf(f1,"%d\t%f\n",clsSize[i],100.0*(float)clsFreq[i]/((float)tolCls+(float)freeg));
			   }
		   }
		   if(clsOrder[i]>0){
			   totalG=totalG+clsOrder[i];
		   }

		}
		//free(clsSize);
	    //free(clsFreq);
		//free(clsOrder);
		if(sil==1){
		printf("total gephyrins=%d,%d\n",totalG+freeG,totalG+freeg);
		printf("total clusters=%d,%d\n",tolCls,totCls);
		}
		////

		////receptors clusters
		int rnum;
		int tolR;
		tolR=0;
		int gi;
		gi=-1;
		int *clsR;//receptors clusters
	    clsR=(int*)malloc(N*sizeof(int));
	    int *clsRsize,*clsRfreq;
	    clsRsize=(int*)malloc(N*sizeof(int));
	    clsRfreq=(int*)malloc(N*sizeof(int));

		for(i=0;i<N;i++){
			clsRsize[i]=0;
			clsRfreq[i]=0;
			rnum=0;
			for(j=1;j<=N;j++){
				gi=cls[i][j];
				if((rec1[gi]>=0)&&(gi>=0)){
					if(sil==1){
					printf("gi=%d,i=%d,j=%d\n",gi,i,j);
					}
					rnum++;
					tolR++;
				}
				gi=-1;
			}
			clsR[i]=rnum/2;
		}
		tolR=tolR/2;
		if(sil==1){
		printf("clsR before: ");
		for(i=0;i<N;i++){
		   printf("%d\t",clsR[i]);
		}
		printf("\n");
		}
		bubbleSort(clsR, N);
		if(sil==1){
		printf("clsR after: ");
		for(i=0;i<N;i++){
		   printf("%d\t",clsR[i]);
		}
		printf("\n");
		}

		loc=0;
		tolRcls=0;
		for(i=0;i<N;i++){
		   j=i+1;

		   if(clsR[j]!=clsR[i]){
			   clsRsize[i]=clsR[i];
			   clsRfreq[i]=j-loc;
			   loc=j;
		   }
		   if(clsR[i]==0){
			   break;
		   }
		   tolRcls++;
		}
        //free(clsR);

		for(i=0;i<N;i++){
		   if(clsRsize[i]!=0){
			   if(((tau+1)%divider)==0){
			   fprintf(f3,"%d\t%f\n",clsRsize[i],100.0*(float)clsRfreq[i]/((float)tolRcls+(float)freeR));
			   }
		   }

		}
	    //free(clsRsize);
	    //free(clsRfreq);
		freeRe=0;
		for(i=0;i<Nr;i++){
			if(freer[i]==0){
				freeRe++;
			}
		}
		tolR=tolR+freeR;
		if(sil==1){
		printf("total free receptors=%d,%d\n",freeR,freeRe);
		printf("total receptors=%d,%d\n",tolR,Nr);
		}

		///Neareat Neighbour of receptors
		kn=2;
	    tree = kdtree_init(dim);
	    if(tree==NULL){
		   exit(-1);
	    }
	    for(i=0;i<Nr;i++){
			double receptor[]={xr[i],yr[i]};
			kdtree_insert(tree, receptor);
	    }
	    kdtree_rebuild(tree);
	    kdtree_dump(tree);
	    double *nn;//distance of nearest neighbor
	    nn=(double*)malloc(Nr*sizeof(double));
		for(l=0;l<Nr;l++){
		   double targetr[2] = {xr[l],yr[l]};
		   kdtree_knn_search(tree, targetr, kn);
		   d=kdtree_knn_dump(tree);
		   knn_clear(tree);
		   nn[i]=d[1][0];
		   if(((tau+1)%divider)==0){
		   fprintf(f6,"%f\n",nn[i]);
		   }
		}
		//free(nn);
		kdtree_destroy(tree);

    //if(((tau+1)%divider)==0){
		    for(i=0;i<N;i++){
			     fprintf(f8,"%f\t%f\t%f\t%f\n",x[i],y[i],x1[i],y1[i]);
			  }
			  for(i=0;i<Nr;i++){
				   fprintf(f9,"%f\t%f\t%d\n",xr[i],yr[i],freer[i]);
			  }
				fclose(f1);
				fclose(f3);
				fclose(f6);
				fclose(f8);
				fclose(f9);
		//}
	 }
   /*
	 for(i=0;i<N;i++){
		 if(beCls[i]>0){
			 printf("after------i=%d,x=%f,y=%f,x1=%f,y1=%f,node=%d,node1=%d,beCls=%d\n",i,x[i],y[i],x1[i],y1[i],node[i],node1[i],beCls[i]);
		 }
	 }

	 for(i=0;i<N;i++){
		 for(j=0;j<(N+1);j++){
			 printf("%d\t",cls[i][j]);
		 }
		 printf("\n");
	 }
   */
	}

	for(i=0;i<N;i++){
		fprintf(f2,"%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%f\n",x[i],y[i],x1[i],y1[i],rec1[i],node[i],node1[i],beCls[i],pnRand[i],alpha[i]);
	}

	for(i=0;i<Nr;i++){
		fprintf(f4,"%f\t%f\t%d\n",xr[i],yr[i],freer[i]);
	}

	for(i=0;i<N;i++){
		for(j=0;j<N+1;j++){
			fprintf(f5,"%d\t",cls[i][j]);
		}
		fprintf(f5,"\n");
	}

  fprintf(f7,"%d\t",freeg);
	fprintf(f7,"%d\t",tolCls);
	fprintf(f7,"%d\t",freeR);
	fprintf(f7,"%d\t",newCls);

  fclose(f2);
  fclose(f4);
  fclose(f5);
	fclose(f7);
	
    /*
	for(i=0;i<(N+1);i++){
       free(cls[i]);//释放第二维指针
    }
    free(cls);//释放第一维指针


	free(k);
	free(alpha);
	free(x);
	free(x1);
	free(y);
	free(y1);
	free(pnRand);
	free(node);
	free(node1);
	free(beCls);



	free(xr);
	free(yr);
	free(freer);
	free(rec1);
    */
	return 0;

} // close main
