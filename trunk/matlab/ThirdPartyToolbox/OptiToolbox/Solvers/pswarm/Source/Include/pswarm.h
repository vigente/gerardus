#ifndef PSwarm_included

#define PSwarm_included


enum { /* Exit codes */ 
  EXIT_OK = 0,
  EXIT_KO = 1,
  EXIT_MEM = 2,
  EXIT_INITIAL =3
};

#define Inf 1.0E20



struct swarm {
  double *x;      /* particles */
  double *v;      /* velocities */
  double *y;      /* best particle so far */
  int *active; /* particle is active */
  double *fx;
  double *fy;     /* best particle so far objective function */
  double delta;   /* poll step size */
  double scale;   /* scale for active constraints */
};



struct Options {
  int s;              /* swarm size */
  double mu, nu;      /* cognitial and social parameters */
  double maxvfactor;  /* maximum velocity factor */
  int maxiter;        /* maximum of iterations */
  int maxf;           /* maximum of function evaluations */
  double iweight;     /* initial weight */
  double fweight;     /* final weight */
  double n2grd;       /* tolerance for gradient norm */
  double blim;        /* bound limit */
  double tol;         /* Tolerance for stopping */
  double delta;       /* initial delta */
  double fdelta;      /* factor for initial delta */
  double idelta;      /* increase delta */
  double ddelta;      /* decrease delta */
  int pollbasis;      /* Poll order */
  double EpsilonActive; /* Epsilon for active constraints */
  int IPrint;         /* print information each IPrint iterations */
  double (*outfcn)();      /* function pointer for printing */
  int vectorized;     /* objective function called with a vector of points */
};


struct Stats {
	int objfunctions;
	int pollsteps;	
	int sucpollsteps;
	int solveriters;
};

enum {N2};

#define SYS_RANDOM 1 

#endif
