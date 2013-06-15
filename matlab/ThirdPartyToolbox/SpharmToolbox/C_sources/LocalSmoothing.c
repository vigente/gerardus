/* 
** Spherical Harmonic Modeling and Analysis Toolkit (SPHARM-MAT) is a 3D 
** shape modeling and analysis toolkit. 
** It is a software package developed at Shenlab in Center for Neuroimaging, 
** Indiana University (email: SpharmMat@gmail.com)
** It is available to the scientific community as copyright freeware 
** under the terms of the GNU General Public Licence.
** 
** Copyright 2009, 2010, ShenLab, Center for Neuroimaging, Indiana University
** 
** This file is part of SPHARM-MAT.
** 
** SPHARM-MAT is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** SPHARM-MAT is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with SPHARM-MAT. If not, see <http://www.gnu.org/licenses/>.
*/

/* This program used some functions written by John Burkardt and publicly 
available at http://orion.math.iastate.edu/burkardt/c_src/geometryc/geometryc.html
*/

/*
** file:    locsmo.c
**
** created: 11/02/2003 by Li Shen
**
** purpose: Do local smoothing for parameterization
*/

#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>

# ifndef _GENCLASS_
# define _GENCLASS_

#define FALSE	0
#define TRUE	1
#define PI		3.14159265358979323846

#if defined(__OS2__) || defined(__WINDOWS__) || defined(WIN32) || defined(_MSC_VER)
#define BLASCALL(f) f
#else
#define BLASCALL(f) f ## _
#endif

/* for calculate metric   */
#define METRICNUM	11
#define BAREAN		0
#define SQEOBJ		1
#define SQEPRM		2
#define ASROBJ		3
#define ASRPRM		4
#define L2OBJ		5
#define L2PRM		6
#define L2BOBJ		7
#define L2BPRM		8
#define LINF		9
#define LINFB		10

#if !defined(MAX)
#define max(A, B) ((A) > (B) ? (A):(B))
#endif

#if !defined(MIN)
#define min(A, B) ((A) < (B) ? (A):(B))
#endif

/* define a mapping structure */
typedef struct mapping_struct{
	int num[3]; /* pnum vnum fnum */
	double *param_vs;
	double *obj_vs;
	int *faces;
	double ttobja; /*total object area */
	double *obj_area; /* scaled object area */
	double *obj_vtarea; /*scaled submesh area around a vertex */
	double *prm_area; /* scaled parameter area */
	double *prm_vtarea; /* parameter submesh area around a vertex  */
	int **nbs; /* #nbfs+nbfs+nbvs */
	double *P;
	double *Q;
	double *fs2d; /* 3dto2d faces */
	int *fcix; /* face indices */
	double *lmd1; /* lambda 1 */
	double *lmd2; /* lambda 2 */
} MAPPING;

/* define a 2d regular mesh structure */
typedef struct regmesh_struct{
	int    n;    /* n even number; */
	double **xs; /* size 0-n */
	double **ys; /* size 0-n  */
	double **zs; /* size 0-n  */
	double **gs; /* grid values size 0-n */
	double **X;  /* theta -pi/2--pi/2, size 0-n/2     */
	double **Y;  /* phi -pi--pi, size 0-n/2   */
} REGMESH;

/* define a edgelist structure */
typedef struct edgelist_struct{
	int v1;
	int v2;
	int ord; /* center's ordering  */
	int fc;
	struct edgelist_struct *next;
} EDGELIST;

/* define a linear equation structure */
typedef struct leqs_struct{
	int n;
	int m;
	double *A;
	double *b;
	double *x;
} LEQS;

/* define a submesh smoothing parameter structure */
typedef struct smoparam_struct{
	int usenew; /* use new position or old position  */
	int ow; int nw; /* (ow*old + nw*new)/(ow+nw) => new   */
} SMOPARAM;

/********************************************************************************************************/
/*void smooth(char *InName, char *OutName);  */

void glbsmo(MAPPING *map, int iter);

REGMESH *create_regmesh(int n);

MAPPING *create_remeshmap(MAPPING *map,REGMESH *regmesh);  

void adjust_param(MAPPING *remeshmap,REGMESH *regmesh);

void create_fs2d(MAPPING *map);

void rotate_3dto2d(double *rv,double x1,double x2,double x3,
				   double y1,double y2,double y3,double z1,double z2,double z3);

void create_nbs(MAPPING *map);                                 

void insert_elist(EDGELIST **elists, int v0, int v1, int v2, int ord, int fc);   

void move_elists_to_mapnbs(MAPPING *map, EDGELIST **elists);

void delete_mapnbs(MAPPING *map);

double triangle_area_3d ( double x1, double y1, double z1, double x2, 
                          double y2, double z2, double x3, double y3, double z3 );

double sph_angle ( double x1, double y1, double z1, double x2, 
	               double y2, double z2, double x3, double y3, double z3 );

void cross_product ( double Ux, double Uy, double Uz, double Vx, 
	                 double Vy, double Vz, double *N );

double enorm0_3d ( double x0, double y0, double z0, double x1, double y1, double z1 );


void solve_dgels (LEQS *leqs); /* solve LLS using QR or LQ factorization*/

void local_smooth(char *InName, char *OutName);

void locsmo(MAPPING *map, int iter);

unsigned relocate(double *newvs, double *verts, int ix, MAPPING *map);

unsigned cal_new_posi(double *nc, double *lvs,
				  MAPPING *map, int *nbfs, int *ords, double *stda);

double *lsq_area_posi(double *nc, double *lvs,
					  MAPPING *map, double *stda);

double max_stretch(double cx, double cy, double cz, double *bx, double *by, 
				   double *bz, int n, MAPPING *map, int *nbfs, int *ords, double *stda);

void locate_nc(double *nc, double cx, double cy, 
			   double *bx, double *by, int n, double *stda);

double *submesh_area_2d(double cx, double cy, 
						double *bx, double *by, int n);

double signed_2d_area(double x1,double y1,
					  double x2,double y2,double x3,double y3);

double *extract_locmesh_verts(double *verts, int vnum, int ix, 
							  int *nbvs, int nbnum);

void proj2plane(double nx,double ny,double nz,
				double lx,double ly,double lz,
				double *fx,double *fy,double *fz);

void set_smoparam(int usenew, int nw, int ow);    

void create_obj_area(MAPPING *map);

void locate_points(MAPPING *map);

void map_points(MAPPING *map);

MAPPING *new_mapping(void);

void delete_mapping(MAPPING *map);

double *new_double(int size);

LEQS *new_leqs(void);

void delete_leqs(LEQS *leqs);

int *new_int(int size);

EDGELIST *new_edgelist(int v1, int v2, int ord, int fc);

void delete_regmesh(REGMESH *regmesh);

REGMESH *new_regmesh(int n);

void delete_dbarray(double **A);

double **new_dbarray(int n);

MAPPING *get_locsmo_values(char *fname, int *reso, double *extents);

void put_vector_values(char *fname, double *vector, int len); 

void put_verts_values(char *fname, double *verts, int vnum);

void get_whole_metric(double *metric, MAPPING *map);

void get_local_metric(double *stch,double cx,double cy,double cz,
					  double *bx,double *by,double *bz,int n,
					  MAPPING *map,int *nbfs,int *ords); 

extern int LWORK; /* workspace dimension  */

extern double WORK[]; /* workspace  */

extern double SV[]; /* store single values m*n  */

int RCOND; /* if RCOND<0, machine precision is used instead    */           

extern SMOPARAM SPRM; /* parameter for smoothing   */

/*extern void dgels_();*/  /* solve LLS using QR or LQ factorization */
extern void BLASCALL(dgels)();  /* solve LLS using QR or LQ factorization */

/*******************************************************************************************************/

#define LWSIZE 10000
#define SVSIZE 10000
int LWORK = LWSIZE; /* workspace dimension  */
double WORK[LWSIZE]; /* workspace  */
double SV[SVSIZE]; /* store single values m*n */
int RCOND = -1; /* if RCOND<0, machine precision is used instead   */

SMOPARAM SPRM; /* parameter for smoothing */

/* writes the usage if no functions given */
void usage(char *MyName)
{
	printf("\nusage: %s infile outfile taskname\n\n",MyName);
}

/****************************************************************************************/
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
	time_t   tStart, tFinish;
	clock_t  cStart, cFinish;
	double   tDura, cDura;
	char     *InName, *OutName;
	

	/* check for proper number of arguments */
    if(nrhs!=2) 
      mexErrMsgTxt("Two input parameters required.");
    
    /* input must be a string */
    if ( mxIsChar(prhs[0]) != 1 || mxIsChar(prhs[1]) != 1)
      mexErrMsgTxt("Input parameters must be a string.");

    /* get InName and OutName from input parameters    */
    InName = mxArrayToString(prhs[0]);
	OutName = mxArrayToString(prhs[1]);

	printf("spa: <Li Shen 2003>\n");
	printf("InName: %s; OutName: %s. ", InName, OutName);

	cStart = clock();
	tStart = time(NULL);

	printf("Task: Local smoothing ...\n");
	local_smooth(InName, OutName);
	
	cFinish = clock();
	tFinish = time(NULL);
	
	cDura = (double) (cFinish - cStart) / CLOCKS_PER_SEC;
	tDura = difftime(tFinish, tStart);
	printf("The cpu time: %f \n", cDura);
	printf("The calender time: %f \n", tDura);
}


/*************************************************************************************************/

void local_smooth(char *InName, char *OutName)
{
	double metric[METRICNUM];
	MAPPING *map;
	int		 reso;
	double	 extents;

	set_smoparam(TRUE,1,10); /* usenew, nw, ow  */
	map = get_locsmo_values(InName,&reso,&extents);
	create_nbs(map);
	create_fs2d(map);
	create_obj_area(map); /* also create space for the areas in the parameter space */
	locsmo(map,reso);
	get_whole_metric(metric,map); /* calculate the metric of the mapping */
	put_vector_values("metric",metric,METRICNUM);  /* put metric value into file */
	put_verts_values(OutName, map->param_vs, map->num[1]);
	delete_mapping(map); 
}


void locsmo(MAPPING *map, int iter)
{
	int vnum, fnum, i, k;
	double *newpvs, *pvs;
	unsigned changed;
	double prevcost=0;

	printf("Local smoothing, %d iterations ...\n",iter);
	vnum = map->num[1]; fnum = map->num[2];

	if (SPRM.usenew){
		newpvs = map->param_vs;
	}
	else{
		newpvs = new_double(fnum*3);
	}

	for(k=0;k<iter;k++){ /* do iterations */
		changed=FALSE;
		for(i=0;i<vnum;i++){ /* traverse the surface */
			if (relocate(newpvs,map->param_vs,i,map)) changed=TRUE;
		}
		if (changed){
			pvs = newpvs;
			newpvs = map->param_vs;
			map->param_vs = pvs;
		}
		else{
			break;
		}
	}

	if (!SPRM.usenew) free(newpvs);
}

/*
** extract the local mesh and then calculate the new center vertex
**
** verts: in parameter space
** ix: the index of the current vertex
** stdarea: standard areas
*/

unsigned relocate(double *newvs, double *verts, int ix, MAPPING *map)
{
	int vnum=map->num[1], *nbs=map->nbs[ix];
	int i,nbnum=nbs[0]; /* lvs[0]=nbnum+1   */
	int *nbfs = nbs+1, *nbvs = nbfs+nbnum, *ords = nbvs+nbnum+1;
	double *lvs, *flatvs, newposi[3];
	double *stda,*stdarea=map->obj_area,ttstda=map->obj_vtarea[ix];
	unsigned changed=FALSE;

	/* get standard areas  */
	stda = new_double(nbnum);
	for(i=0;i<nbnum;i++){
		stda[i] = stdarea[nbfs[i]]/ttstda;
	}

	lvs = extract_locmesh_verts(verts,vnum,ix,nbvs,nbnum);
	flatvs = lsq_area_posi(newposi,lvs,map,stda); /* solve using least square  */
	if (cal_new_posi(newposi,lvs,map,nbfs,ords,stda)){
		changed=TRUE;
		/* update the vertex   */
		newvs[ix] = newposi[0]; ix+=vnum;
		newvs[ix] = newposi[1]; ix+=vnum;
		newvs[ix] = newposi[2];
	}

	free(lvs);
	free(flatvs);
	free(stda);

	return changed;
}


unsigned cal_new_posi(double *nc, double *lvs,
				  MAPPING *map, int *nbfs, int *ords, double *stda)
{
	int n=(int)lvs[0],nbnum=n-1;
	double *lx=lvs+1,*ly=lx+n,*lz=ly+n;
	double tdb,tdbx;
	double cc[3],bcc[3],bcc_cost; /* candidate center*/
	int w1,w2,w;
	unsigned changed=FALSE;

	w = 3;
	for(w1=0;w1<=3;w1++){
		/* weighted new position  */
		w2 = w-w1;
		cc[0] = (w1*nc[0]+w2*lx[0])/w;
		cc[1] = (w1*nc[1]+w2*ly[0])/w;
		cc[2] = (w1*nc[2]+w2*lz[0])/w;
		/* map to the sphere to get the new position  */
		tdb = sqrt(cc[0]*cc[0]+cc[1]*cc[1]+cc[2]*cc[2]);
		cc[0] = cc[0]/tdb; cc[1] = cc[1]/tdb; cc[2] = cc[2]/tdb;

		if (w1==0){
			bcc_cost = max_stretch(cc[0],cc[1],cc[2],lx+1,ly+1,lz+1,
				                   nbnum,map,nbfs,ords,stda);
			memcpy(bcc,cc,sizeof(double)*3);
		}
		else{
			tdb = max_stretch(cc[0],cc[1],cc[2],lx+1,ly+1,lz+1,
							  nbnum,map,nbfs,ords,stda);
			if (tdb>0 && (bcc_cost>tdb || bcc_cost<0)){
				memcpy(bcc,cc,sizeof(double)*3);
				bcc_cost = tdb;
				changed = TRUE;
			}
		}
	}
	nc[0]=bcc[0]; nc[1]=bcc[1]; nc[2]=bcc[2];
	*lx=bcc[0];   *ly=bcc[1];   *lz=bcc[2];

	tdbx = max_stretch(bcc[0],bcc[1],bcc[2],lx+1,ly+1,lz+1,
		               nbnum,map,nbfs,ords,stda);

	return changed;
}


double max_stretch(double cx, double cy, double cz, double *bx, double *by, 
				   double *bz, int n, MAPPING *map, int *nbfs, int *ords, double *stda)
{
	double cost,metric[METRICNUM];

	get_local_metric(metric,cx,cy,cz,bx,by,bz,n,map,nbfs,ords);

	if (metric[BAREAN]>0){
		cost = -1;
	}
	else{
		cost = metric[SQEOBJ]*0 + metric[SQEPRM]*0 +
			   metric[ASROBJ]*0 + metric[ASRPRM]*0 +
			   metric[L2OBJ] *0 + metric[L2PRM] *0 +
			   metric[L2BOBJ]*0 + metric[L2BPRM]*0 +
			   metric[LINF]  *0 + metric[LINFB] *1;
		}

	return cost;
}



/* note that the calculated nc is on the tangent plane instead of the sphere  */
double *lsq_area_posi(double *nc, double *lvs,
					  MAPPING *map, double *stda)
{
	int i,ic,n=(int)lvs[0];
	double *lx=lvs+1,*ly=lx+n,*lz=ly+n;
	double nx=*lx,ny=*ly,nz=*lz; /* normal direction, the center node  */
	double *flatvs,*fx,*fy,*fz;
	double ax=fabs(nx),ay=fabs(ny),az=fabs(nz);

    ic = 3;                     /* ignore z-coord */
    if (ax > ay) {
        if (ax > az) ic = 1;    /* ignore x-coord  */
    }
    else if (ay > az) ic = 2;   /* ignore y-coord */

	flatvs = new_double(n*3+1); 
	flatvs[0] = lvs[0];
	fx=flatvs+1; fy=fx+n; fz=fy+n;

	fx[0]=lx[0]; fy[0]=ly[0]; fz[0]=lz[0];
	for(i=1;i<n;i++){
		proj2plane(nx,ny,nz,lx[i],ly[i],lz[i],fx+i,fy+i,fz+i);
	}

	/* ignore coordinates  */
    switch (ic) {
    case 1: /* y,z */
		locate_nc(nc,fy[0],fz[0],fy+1,fz+1,n-1,stda);
		nc[2] = nc[1]; /* z  */
		nc[1] = nc[0]; /* y */
		nc[0]=nx-(ny*(nc[1]-ny)+nz*(nc[2]-nz))/nx; /* x */
		break;
    case 2: /* z,x */
		locate_nc(nc,fz[0],fx[0],fz+1,fx+1,n-1,stda);
		nc[2] = nc[0]; /* z */
		nc[0] = nc[1]; /* x */
		nc[1] = ny-(nx*(nc[0]-nx)+nz*(nc[2]-nz))/ny; /* y */
		break;
    case 3: /* x,y */
		locate_nc(nc,fx[0],fy[0],fx+1,fy+1,n-1,stda);
		nc[2] = nz-(nx*(nc[0]-nx)+ny*(nc[1]-ny))/nz; /* z */
		break;
    }
	
	return flatvs;
}


void locate_nc(double *nc, double cx, double cy, 
			   double *bx, double *by, int m, double *stda)
{
	int i,j,n=2;
	double *Ax,*Ay,*B,*area,ttarea=0,avarea,xd,yd,calarea;
	LEQS *leqs;

	area = submesh_area_2d(cx,cy,bx,by,m);

	leqs = new_leqs();
	leqs->m = m; leqs->n = n;

	leqs->A = new_double(m*n); Ax=leqs->A; Ay=Ax+m;
	leqs->b = new_double(m); B = leqs->b;
	leqs->x = new_double(m);

	for(i=0;i<m;i++){
		ttarea+=area[i];
	}
	avarea = ttarea/m;

	for(i=0;i<m;i++){
		j = (i+1)%m;
/*		area = (x1-x3)*(y2-y3)-(x2-x3)*(y1-y3);  */
		xd = (bx[i]-bx[j]);
		yd = (by[i]-by[j]);
		Ax[i] = yd; 
		Ay[i] = -xd;
		calarea = ttarea*stda[i];

		B [i] = calarea+bx[j]*yd-xd*by[j];
	}

	free(area);

	solve_dgels (leqs); /* solve LLS using QR or LQ factorization */

	nc[0] = leqs->x[0]; nc[1] = leqs->x[1];

	delete_leqs(leqs);

}


/* need absolute area   */
double *submesh_area_2d(double cx, double cy, double *bx, double *by, int n)
{
	int i,j;
	double *area=new_double(n);

	for(i=0;i<n;i++){
		j=(i+1)%n;
		area[i] = ((bx[i]-cx)*(by[j]-cy)-(bx[j]-cx)*(by[i]-cy))/2;
	}

	return area;
}

double signed_2d_area(double x1,double y1,
					  double x2,double y2,double x3,double y3)
{
	double area = (x1-x3)*(y2-y3)-(x2-x3)*(y1-y3);
	return area/2;
}


/* project l onto the plane of normal n to get f
   n has a length of 0   */
void proj2plane(double nx,double ny,double nz,
				double lx,double ly,double lz,
				double *fx,double *fy,double *fz)
{
	double k=nx*(lx-nx)+ny*(ly-ny)+nz*(lz-nz);

	*fx = lx-nx*k; *fy = ly-ny*k; *fz = lz-nz*k;
}


double *extract_locmesh_verts(double *verts, int vnum, int ix, int *nbvs, int nbnum)
{
	double *lvs;
	int i;
	double *lvsx,*lvsy,*lvsz;
	double *vsx=verts,*vsy=vsx+vnum,*vsz=vsy+vnum;

	lvs = new_double((vnum+1)*3+1);
	lvs[0] = nbnum+1; /* number of vertices  */
	lvsx=lvs+1; lvsy=lvsx+nbnum+1; lvsz=lvsy+nbnum+1;

	lvsx[0] = vsx[ix]; lvsy[0] = vsy[ix]; lvsz[0] = vsz[ix];
	lvsx++; lvsy++; lvsz++; /* skip the first element, which is the center  */
	for(i=0;i<nbnum;i++){
		lvsx[i] = vsx[nbvs[i]]; 
		lvsy[i] = vsy[nbvs[i]]; 
		lvsz[i] = vsz[nbvs[i]];
	}

	return lvs;
}


void set_smoparam(int usenew, int nw, int ow)
{
	SPRM.usenew = usenew; /* update vertices after a whole scan  */
	SPRM.nw = nw; /* weight on the new position  */
	SPRM.ow = ow; /* weight on the old position  */
	printf("Smooth parameters: usenew = %d, nw = %d, ow = %d\n",
		    usenew, nw, ow);
}

void create_obj_area(MAPPING *map)
{
	double *verts=map->obj_vs;
	int *faces=map->faces;
	int vnum=map->num[1], fnum=map->num[2];
	double ttarea,*area,*vtarea;
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,tdb;
	int i,a,b,c;

	/* also create space for areas in the parameter space  */
	area = new_double(fnum);   map->prm_area = area;
	vtarea = new_double(vnum); map->prm_vtarea = vtarea;

	/* create areas in the object space  */
	area   = new_double(fnum);
	vtarea = new_double(vnum);
	ttarea = 0;

	for(i=0;i<vnum;i++) vtarea[i] = 0;

	for(i=0;i<fnum;i++){
		a = faces[i];        x1 = verts[a]; y1 = verts[a+vnum]; z1 = verts[a+vnum*2];
		b = faces[i+fnum];   x2 = verts[b]; y2 = verts[b+vnum]; z2 = verts[b+vnum*2];
		c = faces[i+fnum*2]; x3 = verts[c]; y3 = verts[c+vnum]; z3 = verts[c+vnum*2];
		tdb = triangle_area_3d(x1,y1,z1,x2,y2,z2,x3,y3,z3);
		area[i]=tdb; vtarea[a]+=tdb; vtarea[b]+=tdb; vtarea[c]+=tdb;
		ttarea += tdb;
	}

	printf("Total area: %f\n",ttarea);

	for(i=0;i<fnum;i++) area[i] = area[i]/ttarea;

	for(i=0;i<vnum;i++) vtarea[i] = vtarea[i]/ttarea;

	map->obj_area   = area;
	map->obj_vtarea = vtarea;
	map->ttobja     = ttarea;
}


MAPPING *get_locsmo_values(char *fname, int *reso, double *extents)
{
	FILE *fid;
	MAPPING *map;
	int vsize, fsize;

	if(!(fid=fopen(fname,"rb"))) {
		printf("can't open file %s ...\n", fname);
		exit(1);
	}
	else{
		printf("succeed in opening file %s ... \n", fname);
	}

	map = new_mapping();

	fread(map->num,sizeof(int),3,fid); 
	*reso = map->num[0];

	printf("reso: %d; vnum: %d; fnum %d\n",map->num[0],map->num[1],map->num[2]);

	vsize = 3*map->num[1]; fsize = 3*map->num[2];

	map->param_vs = new_double(vsize);
	map->obj_vs = new_double(vsize);
	map->faces = new_int(fsize);

	fread(extents,sizeof(double),1,fid);
	fread(map->param_vs,sizeof(double),vsize,fid);
	fread(map->obj_vs,sizeof(double),vsize,fid);
	fread(map->faces,sizeof(int),fsize,fid);

	fclose(fid);

	return map;
}


void create_nbs(MAPPING *map)
{
	int i, vnum=map->num[1], fnum=map->num[2], *faces=map->faces, v0, v1, v2;
	EDGELIST **elists;

	if ((elists=(EDGELIST**)malloc(sizeof(EDGELIST*)*vnum))==NULL){
		fprintf(stderr, "not enough memory\n");
		exit(1);
	}
	for(i=0;i<vnum;i++) elists[i]=NULL;

	/* insert to elists */
	for(i=0;i<fnum;i++){
		v0 = faces[i]; v1 = faces[i+fnum]; v2 = faces[i+fnum*2];
		insert_elist(elists,v0,v1,v2,0,i);
		insert_elist(elists,v1,v2,v0,1,i);
		insert_elist(elists,v2,v0,v1,2,i);
	}
    
	/* move to mapnbs  */
	move_elists_to_mapnbs(map,elists);
	/* show_nbs(map->nbs, map->num[1]);  */

	free(elists);
}


void insert_elist(EDGELIST **elists, int v0, int v1, int v2, int ord, int fc)
{
	EDGELIST *elist;
	elist = new_edgelist(v1,v2,ord,fc);
	elist->next = elists[v0];
	elists[v0] = elist;
}


void move_elists_to_mapnbs(MAPPING *map, EDGELIST **elists)
{
	int i,vnum=map->num[1];
	int nbvs[200],nbfs[200],ords[200],n,*curnbs;
	EDGELIST *elist, *cedge, *prev;

	if ((map->nbs=(int**)malloc(sizeof(int*)*vnum))==NULL){
	    printf("matlab crashes at move_elists_to_mapnbs\n");
		fprintf(stderr, "not enough memory\n");
		exit(1);
	}
	/* move elists to nbs   */
	for(i=0;i<vnum;i++){
/*        printf("loop count: %d\n", i); */
		elist=elists[i];
		nbvs[0]=elist->v1; n=0;
		while (elist!=NULL){
			prev = NULL; cedge = elist;
			while (cedge->v1!=nbvs[n]){
				prev = cedge; cedge=cedge->next;
			}
			if (prev==NULL){
				elist = cedge->next;
			}
			else{
				prev->next = cedge->next;
			}
			cedge->next = NULL;
			nbfs[n] = cedge->fc;
			ords[n] = cedge->ord;
			nbvs[++n] = cedge->v2;
			free(cedge);
        /*    printf("end of nested while loop\n"); */
		}
   /*     printf("Outside the nested while loop\n"); */
		curnbs = new_int((n+1)*2+n);
		curnbs[0] = n;
		memcpy(curnbs+1,    nbfs,sizeof(int)*n);
		memcpy(curnbs+1+n,  nbvs,sizeof(int)*(n+1));
		memcpy(curnbs+2+n*2,ords,sizeof(int)*n);
		map->nbs[i]=curnbs;		
		elists[i]=NULL;
   /*     printf("end of for loop: %d\n", i);  */
	}
/*    printf("while for loop is running\n");  */
}


void create_fs2d(MAPPING *map)
{
	int vnum, fnum, i;
	double *fs2d,*u1,*u2,*u3,*v1,*v2,*v3,*x,*y,*z,rv[9];
	double x1,x2,x3,y1,y2,y3,z1,z2,z3;
	int *f1,*f2,*f3;

	vnum = map->num[1]; fnum = map->num[2];
	x =map->obj_vs; y =x+vnum;  z =y+vnum;
	f1=map->faces;  f2=f1+fnum; f3=f2+fnum;

	fs2d = new_double(fnum*6);
	u1=fs2d;    u2=u1+fnum; u3=u2+fnum;
	v1=u3+fnum; v2=v1+fnum; v3=v2+fnum;

	for(i=0;i<fnum;i++){

		x1 = x[f1[i]]; y1 = y[f1[i]]; z1 = z[f1[i]];
		x2 = x[f2[i]]; y2 = y[f2[i]]; z2 = z[f2[i]];
		x3 = x[f3[i]]; y3 = y[f3[i]]; z3 = z[f3[i]];

		rotate_3dto2d(rv,x1,x2,x3,y1,y2,y3,z1,z2,z3);

		u1[i]=rv[0]; u2[i]=rv[1]; u3[i]=rv[2];
		v1[i]=rv[3]; v2[i]=rv[4]; v3[i]=rv[5];
	}
	map->fs2d = fs2d;
}

void get_whole_metric(double *metric, MAPPING *map)
{
	int vnum = map->num[1], fnum = map->num[2];
	double *verts = map->param_vs;
	int *faces = map->faces;
	double *area=map->prm_area,ttarea=0,x1,y1,z1,x2,y2,z2,x3,y3,z3; /* in parameter space  */
	double *stda=map->obj_area; /* in object space  */
	double dist=0,d;
	int i,j;
	int barean=0;
	double sqeobj=0,sqeprm=0,asrobj=0,asrprm=0,stchobj=0,stchprm=0;
	double l2obj=0,l2prm=0,l2bobj=0,l2bprm=0,linf=0,linfb=0;

	double *fs2d=map->fs2d,ttobja=map->ttobja;
	double ttprma = 4*PI, sf=ttobja/ttprma; /* note that this is an approximation  */
	double *s1=fs2d,*s2=s1+fnum,*s3=s2+fnum;
	double *t1=s3+fnum,*t2=t1+fnum,*t3=t2+fnum;
	double A,Ss[3],St[3],a,b,c,maxsv,minsv;

	/* need to fix the parameter polyhedron area? use the approximation  */

	for(i=0;i<fnum;i++){
		/* prepare indices, calculate areas  */
		j = faces[i];        x1 = verts[j]; y1 = verts[j+vnum]; z1 = verts[j+vnum*2];
		j = faces[i+fnum];   x2 = verts[j]; y2 = verts[j+vnum]; z2 = verts[j+vnum*2];
		j = faces[i+fnum*2]; x3 = verts[j]; y3 = verts[j+vnum]; z3 = verts[j+vnum*2];
		area[i] = -PI;
		area[i] += sph_angle(x1,y1,z1,x2,y2,z2,x3,y3,z3);
		area[i] += sph_angle(x2,y2,z2,x3,y3,z3,x1,y1,z1);
		area[i] += sph_angle(x3,y3,z3,x1,y1,z1,x2,y2,z2);
		if (area[i]<=0 || area[i]>PI){
			barean++;
		}
		ttarea += area[i];

		A  = stda[i]*ttobja;
		Ss[0] = (x1*(t2[i]-t3[i])+x2*(t3[i]-t1[i])+x3*(t1[i]-t2[i]))/(2*A);
		St[0] = (x1*(s3[i]-s2[i])+x2*(s1[i]-s3[i])+x3*(s2[i]-s1[i]))/(2*A);
		Ss[1] = (y1*(t2[i]-t3[i])+y2*(t3[i]-t1[i])+y3*(t1[i]-t2[i]))/(2*A);
		St[1] = (y1*(s3[i]-s2[i])+y2*(s1[i]-s3[i])+y3*(s2[i]-s1[i]))/(2*A);
		Ss[2] = (z1*(t2[i]-t3[i])+z2*(t3[i]-t1[i])+z3*(t1[i]-t2[i]))/(2*A);
		St[2] = (z1*(s3[i]-s2[i])+z2*(s1[i]-s3[i])+z3*(s2[i]-s1[i]))/(2*A);
		a = sf*(Ss[0]*Ss[0]+Ss[1]*Ss[1]+Ss[2]*Ss[2]);
		b = sf*(Ss[0]*St[0]+Ss[1]*St[1]+Ss[2]*St[2]);
		c = sf*(St[0]*St[0]+St[1]*St[1]+St[2]*St[2]);
		maxsv = sqrt(((a+c)+sqrt((a-c)*(a-c)+4*b*b))/2);
		minsv = sqrt(((a+c)-sqrt((a-c)*(a-c)+4*b*b))/2);

		l2obj  += (a+c)*stda[i]/2;
		l2prm  += (a+c)*area[i]/2;
		linf    = max(linf,maxsv);
		if (maxsv<1) maxsv=1/maxsv;
		if (minsv<1) minsv=1/minsv;
		l2bobj += (maxsv*maxsv+minsv*minsv)*stda[i]/2;
		l2bprm += (maxsv*maxsv+minsv*minsv)*area[i]/2;
		linfb   = max(linfb,max(maxsv,minsv));
	}
	/* stda and area are normalized  */
 	l2obj=sqrt(l2obj);
 	l2prm=sqrt(l2prm/ttarea);
 	l2bobj=sqrt(l2bobj);
 	l2bprm=sqrt(l2bprm/ttarea);

	printf("total_area/pi = %f, %d bad areas\n",ttarea/PI,barean);

	for(i=0;i<fnum;i++){
	}

	/* calculate others   */
	for(i=0;i<fnum;i++){
		area[i] = area[i]/ttarea; /* normalize  */

		d = area[i]/stda[i];
		if (d<1) d=1/d;
		if (area[i]<=0) barean--;
		asrobj+=d*stda[i];
		asrprm+=d*area[i];
		d = area[i]-stda[i]; d=d*d;
		sqeobj+=d*stda[i];
		sqeprm+=d*area[i];
	}

	metric[BAREAN]	= barean;
	metric[SQEOBJ]	= sqeobj;
	metric[SQEPRM]	= sqeprm;
	metric[ASROBJ]  = asrobj;
	metric[ASRPRM]  = asrprm;
	metric[L2OBJ]   = l2obj;
	metric[L2PRM]   = l2prm;
	metric[L2BOBJ]  = l2bobj;
	metric[L2BPRM]  = l2bprm;
	metric[LINF]    = linf;
	metric[LINFB]   = linfb;

	printf("bad areas             : %d\n",barean);
	printf("square error (obj/prm): %g, %g\n",sqeobj,sqeprm);
	printf("area scaling (obj/prm): %g, %g\n",asrobj,asrprm);
	printf("L2  (obj/prm)         : %g, %g\n",l2obj, l2prm);
	printf("L2b (obj/prm)         : %g, %g\n",l2bobj,l2bprm);
	printf("Linf, Linfb           : %g, %g\n",linf,  linfb);
}

/* area has been normalized  */
void get_local_metric(double *metric,double cx,double cy,double cz,
					  double *bx,double *by,double *bz,
					  int n,MAPPING *map,int *nbfs,int *ords)
{
	double *area2d=new_double(n),ttarea2d=0;
	double *area3d=new_double(n),ttarea3d=0, N[3];
	double *obj_area=map->obj_area, ttobja=map->ttobja;
	double *vert2d=new_double(n*6),*fs2d=map->fs2d;
	double *s1=vert2d,*s2=s1+n,*s3=s2+n;
	double *t1=s3+n,*t2=t1+n,*t3=t2+n;
	int fnum = map->num[2];
	double *u1=fs2d,*u2=u1+fnum,*u3=u2+fnum;
	double *v1=u3+fnum,*v2=v1+fnum,*v3=v2+fnum;
	double A,Ss[3],St[3],a,b,c,maxsv,minsv,sum=0,sumarea=0;
	int barean=0;
	double sqeobj=0,sqeprm=0,asrobj=0,asrprm=0,stchobj=0,stchprm=0;
	double l2obj=0,l2prm=0,l2bobj=0,l2bprm=0,linf=0,linfb=0;
	int i,j;
	double d;
	double ttprma = 4*PI, sf=ttobja/ttprma; /* note that this is an approximation  */


	for(i=0;i<n;i++){
		/* determine parameter space submesh area */
		j=(i+1)%n;
		cross_product(bx[i]-cx,by[i]-cy,bz[i]-cz,bx[j]-cx,by[j]-cy,bz[j]-cz,N);
		/* signed only work for vertices close to the unit sphere */
		if (cx*N[0]+cy*N[1]+cz*N[2]<0){
			metric[BAREAN]=1;
			return;
		}
		area3d[i] = triangle_area_3d(cx,cy,cz,
			bx[i],by[i],bz[i],bx[j],by[j],bz[j]);
		ttarea3d+=area3d[i];

		/* determine object space submesh area   */
		area2d[i] = obj_area[nbfs[i]];
		ttarea2d += area2d[i];

		/* prepare for calculating stretch   */
		s1[i] = u1[nbfs[i]]; s2[i] = u2[nbfs[i]]; s3[i] = u3[nbfs[i]];
		t1[i] = v1[nbfs[i]]; t2[i] = v2[nbfs[i]]; t3[i] = v3[nbfs[i]];
	}

	/* calculating stretch */
	for (i=0;i<n;i++){
		j = (i+1)%n;

		A  = area2d[i]*ttobja;
		switch (ords[i]){
		case 0:
			Ss[0] = (cx*(t2[i]-t3[i])+bx[i]*(t3[i]-t1[i])+bx[j]*(t1[i]-t2[i]))/(2*A);
			St[0] = (cx*(s3[i]-s2[i])+bx[i]*(s1[i]-s3[i])+bx[j]*(s2[i]-s1[i]))/(2*A);
			Ss[1] = (cy*(t2[i]-t3[i])+by[i]*(t3[i]-t1[i])+by[j]*(t1[i]-t2[i]))/(2*A);
			St[1] = (cy*(s3[i]-s2[i])+by[i]*(s1[i]-s3[i])+by[j]*(s2[i]-s1[i]))/(2*A);
			Ss[2] = (cz*(t2[i]-t3[i])+bz[i]*(t3[i]-t1[i])+bz[j]*(t1[i]-t2[i]))/(2*A);
			St[2] = (cz*(s3[i]-s2[i])+bz[i]*(s1[i]-s3[i])+bz[j]*(s2[i]-s1[i]))/(2*A); 
			break;
		case 1:
			Ss[0] = (bx[j]*(t2[i]-t3[i])+cx*(t3[i]-t1[i])+bx[i]*(t1[i]-t2[i]))/(2*A);
			St[0] = (bx[j]*(s3[i]-s2[i])+cx*(s1[i]-s3[i])+bx[i]*(s2[i]-s1[i]))/(2*A);
			Ss[1] = (by[j]*(t2[i]-t3[i])+cy*(t3[i]-t1[i])+by[i]*(t1[i]-t2[i]))/(2*A);
			St[1] = (by[j]*(s3[i]-s2[i])+cy*(s1[i]-s3[i])+by[i]*(s2[i]-s1[i]))/(2*A);
			Ss[2] = (bz[j]*(t2[i]-t3[i])+cz*(t3[i]-t1[i])+bz[i]*(t1[i]-t2[i]))/(2*A);
			St[2] = (bz[j]*(s3[i]-s2[i])+cz*(s1[i]-s3[i])+bz[i]*(s2[i]-s1[i]))/(2*A);
			break;
		case 2:
			Ss[0] = (bx[i]*(t2[i]-t3[i])+bx[j]*(t3[i]-t1[i])+cx*(t1[i]-t2[i]))/(2*A);
			St[0] = (bx[i]*(s3[i]-s2[i])+bx[j]*(s1[i]-s3[i])+cx*(s2[i]-s1[i]))/(2*A);
			Ss[1] = (by[i]*(t2[i]-t3[i])+by[j]*(t3[i]-t1[i])+cy*(t1[i]-t2[i]))/(2*A);
			St[1] = (by[i]*(s3[i]-s2[i])+by[j]*(s1[i]-s3[i])+cy*(s2[i]-s1[i]))/(2*A);
			Ss[2] = (bz[i]*(t2[i]-t3[i])+bz[j]*(t3[i]-t1[i])+cz*(t1[i]-t2[i]))/(2*A);
			St[2] = (bz[i]*(s3[i]-s2[i])+bz[j]*(s1[i]-s3[i])+cz*(s2[i]-s1[i]))/(2*A);
			break;
		}
		a = sf*(Ss[0]*Ss[0]+Ss[1]*Ss[1]+Ss[2]*Ss[2]);
		b = sf*(Ss[0]*St[0]+Ss[1]*St[1]+Ss[2]*St[2]);
		c = sf*(St[0]*St[0]+St[1]*St[1]+St[2]*St[2]);
		maxsv = sqrt(((a+c)+sqrt((a-c)*(a-c)+4*b*b))/2);
		minsv = sqrt(((a+c)-sqrt((a-c)*(a-c)+4*b*b))/2);

		l2obj  += (a+c)*area2d[i]/2;
		l2prm  += (a+c)*area3d[i]/2;
		linf    = max(linf,maxsv);
		if (maxsv<1) maxsv=1/maxsv;
		if (minsv<1) minsv=1/minsv;
		l2bobj += (maxsv*maxsv+minsv*minsv)*area2d[i]/2;
		l2bprm += (maxsv*maxsv+minsv*minsv)*area3d[i]/2;
		linfb   = max(linfb,max(maxsv,minsv));
	}
	/* area3d is normalized normalized  */
 	l2obj=sqrt(l2obj/ttarea2d);
 	l2prm=sqrt(l2prm/ttarea3d);
 	l2bobj=sqrt(l2bobj/ttarea2d);
 	l2bprm=sqrt(l2bprm/ttarea3d);


	for(i=0;i<n;i++){
		area2d[i] = area2d[i]/ttarea2d;
		area3d[i] = area3d[i]/ttarea3d;

		d = area3d[i]/area2d[i];
		if (d<1) d=1/d;
		asrobj+=d*area2d[i];
		asrprm+=d*area3d[i];
		d = area3d[i]-area2d[i]; d=d*d;
		sqeobj+=d*area2d[i];
		sqeprm+=d*area3d[i];
	}

	metric[BAREAN]	= barean;
	metric[SQEOBJ]	= sqeobj;
	metric[SQEPRM]	= sqeprm;
	metric[ASROBJ]  = asrobj;
	metric[ASRPRM]  = asrprm;
	metric[L2OBJ]   = l2obj;
	metric[L2PRM]   = l2prm;
	metric[L2BOBJ]  = l2bobj;
	metric[L2BPRM]  = l2bprm;
	metric[LINF]    = linf;
	metric[LINFB]   = linfb;

	free(area2d);
	free(vert2d);
	free(area3d);

}


void put_vector_values(char *fname, double *vector, int len)
{
	FILE *fid;

	if(!(fid=fopen(fname,"wb"))) {
		printf("can't open file %s ...\n", fname);
		exit(1);
	}
	else{
		printf("succeed in opening file %s ... \n", fname);
	}

	fwrite(vector,sizeof(double),len,fid);

	fclose(fid);
} 


void delete_mapping(MAPPING *map)
{
	free(map->param_vs);
	free(map->obj_vs);
	free(map->faces);
	free(map->obj_area);
	free(map->obj_vtarea);
	free(map->prm_area);
	free(map->prm_vtarea);
	delete_mapnbs(map);
	free(map->P);
	free(map->Q);
	free(map->fcix);
	free(map->fs2d);
	free(map->lmd1);
	free(map->lmd2);
	free(map);
}


void put_verts_values(char *fname, double *verts, int vnum)
{
	FILE *fid;

	if(!(fid=fopen(fname,"wb"))) {
		printf("can't open file %s ...\n", fname);
		exit(1);
	}
	else{
		printf("succeed in opening file %s ... \n", fname);
	}

	fwrite(verts,sizeof(double),vnum*3,fid);

	fclose(fid);
}


double *new_double(int size)
{
	double *newdb;
	
	if ((newdb=(double*)malloc(sizeof(double)*size))==NULL){
		fprintf(stderr, "not enough memory\n");
		exit(1);
	}
	return newdb;
}


LEQS *new_leqs(void)
{
	LEQS *leqs;
	
	if ((leqs=(LEQS*)malloc(sizeof(LEQS)))==NULL){
		fprintf(stderr, "not enough memory\n");
		exit(1);
	}

	leqs->A = NULL;
	leqs->b = NULL;
	leqs->x = NULL;

	return leqs;
}


void solve_dgels(LEQS *leqs)
{
	int m=leqs->m, n=leqs->n;
	int NRHS=1, ok;
	double *A=leqs->A, *b=leqs->b, *x=leqs->x;
	int lwork;
	double *work;

	lwork = max(max(m,n),NRHS)*8+min(m,n);

	if (lwork<=LWORK){
		BLASCALL(dgels)("N", &m, &n, &NRHS, A, &m, b, &m, WORK, &LWORK, &ok);      
	}
	else{
		printf(">>> New lwork=%d\n",lwork);
		work = new_double(lwork);
		/* solve LLS using QR or LQ factorization
		dgels_(char *, integer *, integer *, integer *, 
		doublereal *, integer *, doublereal *, integer *, doublereal *, 
		integer *, integer *); */
		BLASCALL(dgels)("N", &m, &n, &NRHS, A, &m, b, &m, work, &lwork, &ok);      
		free(work);
	}
     							
	memcpy(x,b,sizeof(double)*m);
}


void delete_leqs(LEQS *leqs)
{
	free(leqs->A);
	free(leqs->b);
	free(leqs->x);
	free(leqs);
}


/**********************************************************************/

/* Purpose:

    TRIANGLE_AREA_3D computes the area of a triangle in 3D.

  Modified:

    22 April 1999

  Author:

    John Burkardt

  Parameters:

    Input, double X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the (X,Y,Z) 
    coordinates of the corners of the triangle.

    Output, double TRIANGLE_AREA_3D, the area of the triangle.
*/

double triangle_area_3d ( double x1, double y1, double z1, double x2, 
  double y2, double z2, double x3, double y3, double z3 )
{
  double a;
  double alpha;
  double area;
  double b;
  double base;
  double c;
  double dot;
  double height;
/*
  Find the projection of (P3-P1) onto (P2-P1).
*/
  dot = 
    ( x2 - x1 ) * ( x3 - x1 ) +
    ( y2 - y1 ) * ( y3 - y1 ) +
    ( z2 - z1 ) * ( z3 - z1 );

  base = enorm0_3d ( x1, y1, z1, x2, y2, z2 );
/*
  The height of the triangle is the length of (P3-P1) after its
  projection onto (P2-P1) has been subtracted.
*/
  if ( base == 0.0 ) {

    height = 0.0;

  }
  else {

    alpha = dot / ( base * base );
      
    a = x3 - x1 - alpha * ( x2 - x1 );
    b = y3 - y1 - alpha * ( y2 - y1 );
    c = z3 - z1 - alpha * ( z2 - z1 );

    height = sqrt ( a * a + b * b + c * c );

  }

  area = 0.5 * base * height;
 
  return area;
}


void glbsmo(MAPPING *map, int iter)
{
	double metric[METRICNUM];
	REGMESH *regmesh;
	MAPPING *remeshmap;
	double *tdb;
	int i;

	get_whole_metric(metric,map); /* calculate the metric of the mapping  */
	regmesh = create_regmesh(10);

	for(i=0;i<10;i++){
		/* do remeshing  */
		remeshmap = create_remeshmap(map,regmesh);
		locate_points(remeshmap);
		map_points(remeshmap);
		adjust_param(remeshmap,regmesh);

		/* switch the param_vs  */
		tdb = remeshmap->param_vs;
		remeshmap->param_vs = map->param_vs;
		map->param_vs = tdb;
		delete_mapping(remeshmap); 

		locsmo(map,10);

		get_whole_metric(metric,map); /* calculate the metric of the mapping  */
	}
	delete_regmesh(regmesh);
}


MAPPING *new_mapping(void)
{
	MAPPING *map;
	
	if ((map=(MAPPING*)malloc(sizeof(MAPPING)))==NULL){
		fprintf(stderr, "not enough memory\n");
		exit(1);
	}
	map->param_vs = NULL;
	map->obj_vs = NULL;
	map->faces = NULL;
	map->obj_area = NULL;
	map->obj_vtarea = NULL;
	map->prm_area = NULL;
	map->prm_vtarea = NULL;
	map->nbs = NULL;
	map->P = NULL;
	map->Q = NULL;
	map->fs2d = NULL;
	map->fcix = NULL; /* face indices */
	map->lmd1 = NULL; /* lambda 1 */
	map->lmd2 = NULL; /* lambda 2 */

	return map;
}


int *new_int(int size)
{
	int *newint;
	if ((newint=(int*)malloc(sizeof(int)*size))==NULL){
		fprintf(stderr, "not enough memory\n");
		exit(1);
	}
	return newint;
}


EDGELIST *new_edgelist(int v1, int v2, int ord, int fc)
{
	EDGELIST *elist;
	
	if ((elist=(EDGELIST*)malloc(sizeof(EDGELIST)))==NULL){
		fprintf(stderr, "not enough memory\n");
		exit(1);
	}

	elist->v1   = v1;
	elist->v2   = v2;
	elist->ord  = ord;
	elist->fc   = fc;
	elist->next = NULL;

	return elist;
}

void rotate_3dto2d(double *rv,double x1,double x2,double x3,
				   double y1,double y2,double y3,double z1,double z2,double z3)
{
	double N[3],theta,phi,y,z;
	double *nx=rv,*ny=nx+3,*nz=ny+3,a2,a3;

	cross_product(x2-x1,y2-y1,z2-z1,x3-x1,y3-y1,z3-z1,N);
	
	theta = atan2(N[1],N[0]);
	phi = atan2(N[2],sqrt(N[0]*N[0]+N[1]*N[1]));
	z = theta; y = PI/2-phi; /* first z then y  */

	nx[0] = x1*cos(z)*cos(y)+y1*sin(z)*cos(y)-z1*sin(y);
	ny[0] =-x1*sin(z)       +y1*cos(z);
	nz[0] = x1*cos(z)*sin(y)+y1*sin(z)*sin(y)+z1*cos(y); /* should comment  */
	nx[1] = x2*cos(z)*cos(y)+y2*sin(z)*cos(y)-z2*sin(y);
	ny[1] =-x2*sin(z)       +y2*cos(z);
	nz[1] = x2*cos(z)*sin(y)+y2*sin(z)*sin(y)+z2*cos(y); /* should comment  */
	nx[2] = x3*cos(z)*cos(y)+y3*sin(z)*cos(y)-z3*sin(y);
	ny[2] =-x3*sin(z)       +y3*cos(z);
	nz[2] = x3*cos(z)*sin(y)+y3*sin(z)*sin(y)+z3*cos(y); /* should comment  */

	a2=signed_2d_area(nx[0],ny[0],nx[1],ny[1],nx[2],ny[2]);
	a3=triangle_area_3d(x1,y1,z1,x2,y2,z2,x3,y3,z3);
	if (fabs(a2-a3)>10e-10){
		printf("Error: a2 %g != a3 %g",a2,a3);
		exit(1);
	}
}


/**********************************************************************/

/*
  Purpose:

    SPH_ANGLE computes the value of a spherical angle on unit sphere.

  Modified:

    05 October 2003

  Author:

    Li Shen

  Parameters:

    Input, double X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the (X,Y,Z) 
    coordinates of the corners of the triangle.

    Output, double SPH_ANGLE, the value of the spherical angle.   */
/**********************************************************************/

double sph_angle ( double x1, double y1, double z1, double x2, 
	double y2, double z2, double x3, double y3, double z3 )


{
	double angle,x,y;

	y = x1*y2*z3 - x1*z2*y3 + y1*z2*x3 - y1*x2*z3 + z1*x2*y3 - z1*y2*x3;
	x = (x2*x3+y2*y3+z2*z3) - (x1*x3+y1*y3+z1*z3) * (x1*x2+y1*y2+z1*z2);
	angle = atan2(y,x); 

	return angle;
}

void cross_product ( double Ux, double Uy, double Uz, double Vx, 
	double Vy, double Vz, double *N )
{
	N[0] = Uy*Vz-Uz*Vy;
	N[1] = Uz*Vx-Ux*Vz;
	N[2] = Ux*Vy-Uy*Vx;
}

void delete_mapnbs(MAPPING *map)
{
	int i,vnum=map->num[1];

	if (map->nbs!=NULL){
		for(i=0;i<vnum;i++){
			free(map->nbs[i]);
		}
		free(map->nbs);
	}	
}

/**********************************************************************/
/*
  Purpose:

    ENORM0_3D computes the Euclidean norm of (P1-P0) in 3D.

  Modified:

    18 April 1999

  Author:

    John Burkardt

  Parameters:

    Input, double X0, Y0, Z0, X1, Y1, Z1, the coordinates of the points 
    P0 and P1.

    Output, double ENORM0_3D, the Euclidean norm of (P1-P0).
*/
/**********************************************************************/

double enorm0_3d ( double x0, double y0, double z0, double x1, double y1, double z1 )


{
  double value;

  value = sqrt (
    ( x1 - x0 ) * ( x1 - x0 ) + 
    ( y1 - y0 ) * ( y1 - y0 ) + 
    ( z1 - z0 ) * ( z1 - z0 ) );
 
  return value;
}


REGMESH *create_regmesh(int n)
{
	REGMESH *regmesh=new_regmesh(n);
	int i,j,k=n/2;
	double **xs=regmesh->xs,**ys=regmesh->ys,**zs=regmesh->zs;
	double **X=regmesh->X,  **Y=regmesh->Y;
	double theta,phi; /* theta -pi/2--pi/2, size 0-n/2, phi -pi--pi, size 0-n/2   */

	for(i=-n;i<=n;i++){
		for(j=-n;j<=n;j++){
			theta=asin((double)j/n); phi=i*PI/n;
			if (i%2==0 && j%2==0){
				X[i/2+k][j/2+k] = theta;
				Y[i/2+k][j/2+k] = phi;
			}
			if (i%2!=0 && j%2!=0){
				xs[(i-1)/2+k][(j-1)/2+k] = cos(theta)*cos(phi);
				ys[(i-1)/2+k][(j-1)/2+k] = cos(theta)*sin(phi);
				zs[(i-1)/2+k][(j-1)/2+k] = sin(theta);
			}
		}
	}

	return regmesh;
}


MAPPING *create_remeshmap(MAPPING *map,REGMESH *regmesh)
{
	MAPPING *remeshmap = new_mapping();
	int i,j,ix;
	int n=regmesh->n,pnum=n*n,vnum=map->num[1],fnum=map->num[2];
	double *px,*py,*pz;
	double **xs=regmesh->xs,**ys=regmesh->ys,**zs=regmesh->zs;
	double w;
	int *v1=map->faces, *v2=v1+fnum, *v3=v2+fnum;
	int **nbs=map->nbs,*nbfs,nbfsn; /* #nbfs+nbfs+nbvs */
	double *prm_area=map->prm_area, *obj_vtarea=map->obj_vtarea;

	remeshmap->num[0] = pnum; /* pnum  */
	remeshmap->num[1] = vnum; /* vnum  */
	remeshmap->num[2] = fnum; /* fnum  */

	/* create regular mesh vertices  */
	remeshmap->P = new_double(pnum*3);
	px = remeshmap->P; py = px+pnum; pz = py+pnum;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			ix=j*n+i;
			px[ix]=xs[i][j]; py[ix]=ys[i][j]; pz[ix]=zs[i][j];
		}
	}

	/*======= the following needs to be accelerated!!! =======   */

	/* copy param_vs */ 
	remeshmap->param_vs = new_double(vnum*3);
	memcpy(remeshmap->param_vs,map->param_vs,sizeof(double)*vnum*3);

	/* area scaling ratio values */
	remeshmap->obj_vs = new_double(vnum*3);
	memcpy(remeshmap->obj_vs,map->param_vs,sizeof(double)*vnum*3);
	/* do scaling, assume that map->prm_area is updated and map->nbs is created  */
	px = remeshmap->obj_vs; py = px+vnum; pz = py+vnum;
	for(i=0;i<vnum;i++){
		nbfsn = nbs[i][0]; nbfs=nbs[i]+1; w = 0;
		for(j=0;j<nbfsn;j++) w+=prm_area[nbfs[j]];
		w = w/obj_vtarea[i];
		px[i] = px[i]*w; py[i] = py[i]*w; pz[i] = pz[i]*w;
	}

	/* copy faces */
	remeshmap->faces = new_int(fnum*3);
	memcpy(remeshmap->faces,map->faces,sizeof(int)*fnum*3);

	/* reserve space for the result  */
	remeshmap->Q = new_double(pnum*3);

	return remeshmap;
}


void locate_points(MAPPING *map)
{
	int i, j;
	int *fcix; /* face indices */
	double *lmd1; /* lambda 1 */
	double *lmd2; /* lambda 2 */
	int    pnum = map->num[0], vnum = map->num[1], fnum = map->num[2];
	double *Px  = map->P,        *Py  = Px+pnum,  *Pz  = Py+pnum;
	double *pvx = map->param_vs, *pvy = pvx+vnum, *pvz = pvy+vnum;
	int    *A   = map->faces,    *B   = A+fnum,   *C   = B+fnum;
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,x,y,z,xx,yy,zz,k,d1,d2;
	int count=0;

	fcix = new_int(pnum);
	lmd1 = new_double(pnum);
	lmd2 = new_double(pnum);

	for(i=0;i<pnum;i++){
		for(j=0;j<fnum;j++){
			/* points on surface */
			x1 = pvx[A[j]]; y1 = pvy[A[j]]; z1 = pvz[A[j]];
			x2 = pvx[B[j]]; y2 = pvy[B[j]]; z2 = pvz[B[j]];
			x3 = pvx[C[j]]; y3 = pvy[C[j]]; z3 = pvz[C[j]];
			x  = Px[i];     y  = Py[i];     z  = Pz[i];
			xx = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1); 
			yy = (z2-z1)*(x3-x1)-(z3-z1)*(x2-x1); 
			zz = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
			k = (x1*xx + y1*yy + z1*zz)/(x*xx + y*yy + z*zz);
			if (k<=1 && k>0.5){
				/* linear mapping */
				yy = (x1-x3)*(y2-y3)-(y1-y3)*(x2-x3);
				zz = (x1-x3)*(z2-z3)-(z1-z3)*(x2-x3);
				if (yy*yy>zz*zz) {
					d1 = ((x*k-x3)*(y2-y3)-(x2-x3)*(y*k-y3))/yy;
					d2 = ((x1-x3)*(y*k-y3)-(x*k-x3)*(y1-y3))/yy;
				}
				else{
					d1 = ((x*k-x3)*(z2-z3)-(x2-x3)*(z*k-z3))/zz;
					d2 = ((x1-x3)*(z*k-z3)-(x*k-x3)*(z1-z3))/zz;
				}

				/* note that these parameters need to be adjusted */
				if (d1>=-0.000001 && d2>=-0.000001 && (d1+d2)<=1.000001){
					fcix[i] = j; lmd1[i] = d1; lmd2[i] = d2; count++;			
				}
			}
		}
	}

	printf("after projection and linear mapping: %d\n",count);

	map->fcix = fcix;
	map->lmd1 = lmd1;
	map->lmd2 = lmd2;
}


void map_points(MAPPING *map)
{
	int i, j;
	int *fcix; /* face indices */
	double *lmd1; /* lambda 1 */
	double *lmd2; /* lambda 2 */
	int    pnum = map->num[0], vnum = map->num[1], fnum = map->num[2];
	double *Qx  = map->Q,        *Qy  = Qx+pnum,  *Qz  = Qy+pnum;
	double *ovx = map->obj_vs,   *ovy = ovx+vnum, *ovz = ovy+vnum;
	int    *A   = map->faces,    *B   = A+fnum,   *C   = B+fnum;
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,d1,d2;

	fcix = map->fcix;
	lmd1 = map->lmd1;
	lmd2 = map->lmd2;

	for(i=0;i<pnum;i++){
		j = fcix[i]; d1 = lmd1[i]; d2 = lmd2[i];

		x1 = ovx[A[j]]; y1 = ovy[A[j]]; z1 = ovz[A[j]];
		x2 = ovx[B[j]]; y2 = ovy[B[j]]; z2 = ovz[B[j]];
		x3 = ovx[C[j]]; y3 = ovy[C[j]]; z3 = ovz[C[j]];
	
		Qx[i] = d1*(x1-x3) + d2*(x2-x3) + x3;
		Qy[i] = d1*(y1-y3) + d2*(y2-y3) + y3;
		Qz[i] = d1*(z1-z3) + d2*(z2-z3) + z3;
	}
}


/* adjust parameterization  */
void adjust_param(MAPPING *remeshmap,REGMESH *regmesh)
{
	int i,j,n=regmesh->n,pnum=remeshmap->num[0],vnum=remeshmap->num[1],ix;
	double *qx=remeshmap->Q, *qy=qx+pnum, *qz=qy+pnum;
	double **gs=regmesh->gs,min_gs=1,max_gs=1,tdb,tdb2;
	double **ts=regmesh->X,**ps=regmesh->Y,*t=ts[0];
	double *xs=remeshmap->param_vs,*ys=xs+vnum,*zs=ys+vnum;
	double theta,phi;
	
	/* get grid value, meaning the height  */
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			ix = j*n+i;
			tdb = 1/sqrt(qx[ix]*qx[ix]+qy[ix]*qy[ix]+qz[ix]*qz[ix]);
			gs[i][j] = tdb;
			min_gs = min(min_gs,tdb);
			max_gs = max(max_gs,tdb);
		}
	}

	/* scale the value, can be removed later or modified  */
	max_gs = min_gs*10;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			gs[i][j] = min(gs[i][j],max_gs);
		}
	}

	/* assign t values */
	for(j=0;j<n+1;j++) t[j]=0;
	for(j=0;j<n;j++){
		for(i=0;i<n;i++){
			t[j+1] += gs[i][j];
		}
		t[j+1]+=t[j];
	}

	for(j=0;j<n+1;j++){
		t[j] = t[j]*2/t[n]-1;

	}

	/* do interpolation  */
	for(i=0;i<vnum;i++){
		phi = atan2(ys[i],xs[i]);
		tdb = (zs[i]+1)*n/2;
		j = (int) floor(tdb);
		tdb2 = 2.0*j/n-1;
		zs[i] = t[j]+(t[j+1]-t[j])*(zs[i]-tdb2)*n/2;
		theta = asin(zs[i]);
		xs[i] = cos(theta)*cos(phi);
		ys[i] = cos(theta)*sin(phi);
	}
}
	
void delete_regmesh(REGMESH *regmesh)
{
	delete_dbarray(regmesh->xs);
	delete_dbarray(regmesh->ys);
	delete_dbarray(regmesh->zs);
	delete_dbarray(regmesh->gs);
	delete_dbarray(regmesh->X);
	delete_dbarray(regmesh->Y);
	free(regmesh);
}


REGMESH *new_regmesh(int n)
{
	REGMESH *regmesh;
	
	if ((regmesh=(REGMESH*)malloc(sizeof(REGMESH)))==NULL){
		fprintf(stderr, "not enough memory\n");
		exit(1);
	}

	regmesh->n = n;
	regmesh->xs = new_dbarray(n);   /* center */
	regmesh->ys = new_dbarray(n);	/* center */
	regmesh->zs = new_dbarray(n);	/* center */
	regmesh->gs = new_dbarray(n);	/* grid value */
	regmesh->X  = new_dbarray(n+1); /* boundary */
	regmesh->Y  = new_dbarray(n+1); /* boundary */

	return regmesh;
}


void delete_dbarray(double **A)
{
	free(A[0]);
	free(A);
}


double **new_dbarray(int n)
{
	double **A,*data;
	int i;

	if ((A=(double**)malloc(sizeof(double*)*n))==NULL){
		fprintf(stderr, "not enough memory\n");
		exit(1);
	}
	if ((data=(double*)malloc(sizeof(double)*n*n))==NULL){
		fprintf(stderr, "not enough memory\n");
		exit(1);
	}

	A[0] = data;
	for(i=1;i<n;i++) A[i]=A[i-1]+n;

	return A;
}
#endif
