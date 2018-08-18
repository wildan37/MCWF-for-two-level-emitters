#include "stdio.h"
#include "stdlib.h"
#include "complex.h"
#include "math.h"
#include "time.h"
#include "mt19937ar.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_odeiv2.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_sort_double.h"
#include "gsl/gsl_statistics_double.h"
#include "sys/time.h"

struct timeval tv;
time_t nowtime;
struct tm * nowtm;
char tmbuf[64], buf[64];

struct ODE_params {
    
    double Delta0;
};

//atoms, states and realisations
unsigned state;         //state
int spin_i;        //i-th spin
int spin_j;        //j-th spin
unsigned Nspins;        //spin-number
unsigned atomiclevel;   //Nlevel
unsigned Nstates;       //the number of state within the spin ensemble
unsigned Nbasis;        //the number of basis
unsigned ilattice;      //lattice_i
unsigned jlattice;      //lattice_j
unsigned realisation;   //realisation
unsigned i;             //additional indices
unsigned Nreal;    //number realisation
unsigned iter;
unsigned row;
unsigned col;
unsigned  Nrow;
unsigned Ncol;
unsigned base;
unsigned Niter;
int deph;
int d;
int Lx, Ly, Lz, L;

///time propagation
long int tstep;         //tstep for array indices
unsigned notimesteps;   //number timesteps
double dt;              //short time discrepancy during euler method
double pulselength;     //the length of pulse
double t;               //time
double t1;              //discrete time
double tmax;            //maximum time propagation
double p0;

///interaction range
double cut_off;         //interaction cutoff
double alpha;           //power law
double epsilon;         //power law range
double rsquare;         //vdW interaction

///Free parameters for normal timescale
double Omega;            //probe laser
double Delta;           //laser detuning
double V0;              //interaction strength
double friction;        //meant effect
double Gamma;           //Damping
double rad;
double specwidth;       //spectral width
int boundcond;          //boundary condition
int typeint;            //type of interaction
double Delta_E;         //distance between symmetric and non-symmetric state
double ree;
double mu;

///Free parameters for igor timescale
double xi;
double delta;
double R6;
double R3;
double frictionbar;

///Global noise
double noise;           //global noise
double initnoise;       //initialize noise
double dW;              //white noise processes

///observables
double Nryd;            //<n>
double norm;            //normalisation
double Nryd2;           //<n2>
double numerator;       //numerator for g2
double denom1;          //denominator for g2 atom_i
double denom2;          //denominator for g2 atom_j
double prob;            //normalisation


///experiments
double sumNryd;
double sumprob;
double sumNryd2;
double sumnumerator;
double sumdenom1;
double sumdenom2;
double sumree;

///averaging
double aveNryd;
double aveprob;
double aveNryd2;
double avenumerator;
double avedenom1;
double avedenom2;
double averee;

///storing elements
unsigned *excited_num;  //storing the number Rydberg atoms
unsigned **k;           //storing flipping atom
double *re;             //storing real part of state
double *im;             //storing imaginary part of state
double *localnoise;     //storing localnoise of state
//double *basis;
double *y;
double *xold;
double *dp;
double **collect;

//preparing the file
FILE *outfile;
//FILE *outfile1;
char filename[100];
char filenamestate[100];
char filestat[80];

//function
void init_params(void);
void free_alloc(void);
void init_configuration(void);
void counting_excited_state(void);
void flip_state(void);
void flip_interaction(void);
void euller(void);
void euller_igor(void);
void euller_three(void);
double rand_gauss(void);
double max (double data[], long unsigned Ndata);
void init_inverted_noise(void);
void init_ground_noise(void);
void trajectory_noise (void);
int funcnoise (double t, const double y[], double f[], void *params);
int func_wave (double t, const double y[], double f[], void *params);
void globalnoise(void);
void local_noise(void);
void eulermeruyama_global(int nruns);
void eulermeruyama_local(int nruns);
double fi (double x);
double gi (double x);
void steady_state_run (void);
void run_normal (void);
void spectrum_symmetric (int Nspins);
void steady_noise(void);

