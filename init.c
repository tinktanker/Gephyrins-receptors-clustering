#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#define N 500 // number of gephyrins
#define L 320 // size of lattice
#define Nr 100 // number of receptors
#define pi acos(-1.)
#define geLength 11.0

int main(int argc, char* argv[]){

  if(argc < 2){
   fprintf(stderr, "Give seed/0! seed: for random number; 0: produce seeds by computer\n");
   exit(-1);
  }
  int t;
  int par2;
  par2=atoi(argv[1]);
  if(par2==0){
    t=time(NULL);
  }else{
      t=par2;
  }
  printf("t=%d\n",t);
  srand(t);

  FILE *f1;
  FILE *f2;
  FILE *f3;
  FILE *f4;
  char gefile[30],refile[30],clsfile[30],statfile[30];

  sprintf(gefile,"gephyrin0");
  sprintf(refile,"receptor0");
  sprintf(clsfile,"cluster0");
  sprintf(statfile,"stat0");

  f1 = fopen(gefile,"w");
  f2 = fopen(refile,"w");
  f3 = fopen(clsfile,"w");
  f4 = fopen(statfile,"w");

  double *x,*x1,*y,*y1;//(x,y),circle head; (x1,y1),square head
	x=(double*)malloc(N*sizeof(double));
	x1=(double*)malloc(N*sizeof(double));
	y=(double*)malloc(N*sizeof(double));
	y1=(double*)malloc(N*sizeof(double));
  int *pnRand;
  pnRand=(int*)malloc(N*sizeof(int));
  int *rec1;//mark if a square head taking a receptor:taken,the index of the taken receptor;not,-1.
  rec1=(int*)malloc(N*sizeof(int));
  int i,j;
  int *node,*node1;//node:circle end taken=the index of the adjacent circle, ie. i/j/l; node1:square end taken=the index of gephyrin who took it. or else=-1.
  node=(int*)malloc(N*sizeof(int));
  node1=(int*)malloc(N*sizeof(int));
  int **cls;//mark which gephyrins each cluster contains
  cls=(int**)malloc(N*sizeof(int*));//1st dimension
  for(i=0;i<(N+1); i++){
      cls[i]=(int*)malloc((N+1)* sizeof(int));//2nd dimension
  }
  int *beCls;//label which cluster each gephyrin belongs to(eg. cluster 1, cluster 2,...)
  beCls=(int*)malloc(N*sizeof(int));
  double *xr,*yr;//receptors coordinates
  xr=(double*)malloc(Nr*sizeof(double));
  yr=(double*)malloc(Nr*sizeof(double));
  int *freer;//mark if a receptor joined a cluster:joined,the index of gephyrin who took it;free,0.
  freer=(int*)malloc(Nr*sizeof(int));


  ////gephyrins
  double *k,*alpha;
  k=(double*)malloc(N*sizeof(double));
  alpha=(double*)malloc(N*sizeof(double));
	for(i=0;i<N;i++){
    rec1[i]=-1;
    alpha[i]=2.0*pi*((rand()+0.1)/(RAND_MAX+1.0));
    k[i]=tan(alpha[i]);
		x[i]=rand()/(RAND_MAX + 1.0)*(float)L;
		y[i]=rand()/(RAND_MAX + 1.0)*(float)L;
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
		node[i]=-1;
		node1[i]=-1;
		beCls[i]=0;
    fprintf(f1,"%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%f\n",x[i],y[i],x1[i],y1[i],rec1[i],node[i],node1[i],beCls[i],pnRand[i],alpha[i]);
		for(j=0;j<(N+1);j++){
			cls[i][j]=-1;
      fprintf(f3,"%d\t",cls[i][j]);
		}
    fprintf(f3,"\n");
	}

	////receptors
	for(i=0;i<Nr;i++){
		xr[i]=rand()/(RAND_MAX + 1.0)*(float)L;
		yr[i]=rand()/(RAND_MAX + 1.0)*(float)L;
		freer[i]=0;
    fprintf(f2,"%f\t%f\t%d\n",xr[i],yr[i],freer[i]);
	}

  ////statistics
  int freeg,tolCls,freeR,newCls;
  freeg=N;
  tolCls=0;
  freeR=Nr;
  newCls=0;
  fprintf(f4,"%d\t%d\t%d\t%d\n",freeg,tolCls,freeR,newCls);



  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);

  return 0;
}
