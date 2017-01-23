/*               globalvar.cpp          */
//////////////////////////////////////////
// define global variables		//
//////////////////////////////////////////


#include "parameter.h"

//////////////////////////////////////////
// global variables

//////////////////////////////////////////
// cell listing
long *_list;
long _head = 0;
const int _pperc = 2;			// particles per biological cell
					// DO NOT CHANGE UNLESS YOU KNOW WHAT YOU DO!!!

//////////////////////////////////////////
// cell properties
double *_rx;
double *_ry;
double *_rz;
double *_vx;
double *_vy;
double *_vz;
double *_fx;
double *_fy;
double *_fz;
double *_fdx;
double *_fdy;
double *_fdz;
long *_id;
int *_spec;				// cell interaction parameters can be different for different species

//////////////////////////////////////////
// global properties
double LX   = BOXLENGTH_X;              // box dimensions
double LY   = BOXLENGTH_Y;
double LZ   = BOXLENGTH_Z;
double curtime;				// current time index
double dt   = TIMESTEP;                 // time step
double tend = SIMTIME;                  // simulation time
double kbt = KBT;

long next_free_id = 0;			// give every cell a unique id.
					// in case of division each newly created cell particle gets new id
					// so that lineages can be restored from ids
long lastcn;				// size of _list and therefore (size of _r,...)/pperc
long lastn;				// size of arrays
long maxn;				// max size of arrays
MTRand *drand;				// random number generator handle



//////////////////////////////////////////
// species properties
double *rpp;				// see parameter.h for explanation
double *r0;
double *rt;
double *rc;
double *rct;
double *rz;
double *b;
double *gammac;
double *gammat;
double *gammab;
double *gammas;
double *ka;
double *f0;
double *f1;
double *fa;
double *m;				// mass of particles


//////////////////////////////////////////
// other calculated properties of tissue
double com_rx = 0.;
double com_ry = 0.;
double com_rz = 0.;
long   com_n  = 0;


//////////////////////////////////////////
// other stuff
double max_ia_range;			// maximum interaction range
int trajdatastep = TRAJDATASTEP;
int savedatastep = SAVEDATASTEP;
unsigned long rngseed = RNGSEED;
double epsilon = EPSILON;		// epsilon added to distance calculation to circumvent divergence for small dr

////////////////////////////////////////
// Verlet list
const double skin = 0.3;		// max(skin+rpp, skin+rt) is verlet list cut off radius
double currc = RPP+skin;			// store currently used min skin radius
long *vlist;				// verlet list
long *vlist_offset;			// offset vector
long *vlist_nc;				// vlist used entry length vector (nc=neighbour count)
double *vlist_dx;			// displacement in x,y,z direction since last verlet list update
double *vlist_dy;
double *vlist_dz;

////////////////////////////////////////
// special stuff
double *rna;				// random number array
unsigned char *rnai;      // save which fields are used in rna (much faster to zero out)
double *randf_prefc;			// random force prefactor
double *randf_preft;			// random force prefactor
int vlist_rebuild = 0;

double boxlength[3] = {BOXLENGTH_X, BOXLENGTH_Y, BOXLENGTH_Z};


long t_begin;				// for actual run time measurement

//////////////////////////////////////////////////
// error value descriptions			//
//////////////////////////////////////////////////
long _errno = 0;			// set error value (for now only used by load_dump())
short _int_flag = 0;			// interrupt occured
int _int_no = 0;			// interrupt number

const char *_errno_desc[] = {
  "Sucess",
  "Error during fdump.open()",
  "INITSIZE smaller than number of particles in dump file",
  "MAXSPECNUM smaller than number of species in dump file",
  "Boundary conditions differ from those defined in dump file",
  "WARNING: No pressure measurements conducted during dump file creation",
  "WARNING: Pressure measurements were conducted during dump file creation but are deactivated now",
  "Memory allocation error",
  "WARNING: sdiv has reached SDIV_SIZE",
  "Unknown error",
  "Error during fopen(/dev/urandom)"
};
const char *_errno_add_info = "";


#ifdef PMEASUREMENT
double *pm_px0;
double *pm_pxL;
double *pm_py0;
double *pm_pyL;
double *pm_pz0;
double *pm_pzL;
#endif


#ifdef DENSITYPROFILE
 #if (defined (RECTANGULAR) || defined (SPHERICAL))
int *dens_rhoV;				// density profile
int *dens_rhoVsq;
int *dens_rhoVsqtmp;
int *dens_nmeas;				// number of measurements at s
int *dens_p_nmeas;
int *dens_knmeas;		// number of measurements for rates kd/ka (meas makes only sense when cells present)
int *dens_nkd;				// number of divisions at s
int *dens_nka;				// number of deaths at s
double *dens_ka;
double *dens_kd;
double dens_binsize = 0.1;
int dens_length = (int)(LZ/dens_binsize+1);
double dens_maxz = -1.;
double dens_mean_maxz = 0.;
double dens_lastmean_maxz = 0.;
double dens_p_maxz = -1.;
double dens_mean_p_maxz = 0.;
double dens_lastmean_p_maxz = 0.;
double dens_minz = -1.;
double dens_mean_minz = 0.;
double dens_lastmean_minz = 0.;
double dens_p_minz = -1.;
double dens_mean_p_minz = 0.;
double dens_lastmean_p_minz = 0.;
double *dens_kperp;			// k_d perpendicular to surface
double *dens_kpara;			// k_d parallel to surface
int *dens_curn;				// current number of cells (used for kd/a)
int *dens_curnka;			// current number of cell deaths (used for ka)
int *dens_curnkd;			// current number of cell divisions (used for kd)
double *dens_sq;				// nematic order parameter of division (zz component)
double *dens_cell_sq;			// nematic order parameter of cell alignment (zz component)
double dens_com_rx;			// Center of mass r_x
double dens_com_ry;			// -"- r_y
double dens_com_rz;			// -"- r_z
long dens_com_n;
double *dens_mean_cell_dr;		// Mean cell "size", i.e. distance between its two particles
double *dens_mean_NN_dr;			// Mean nearest neighbour distance
 #elif defined (ZYLINDRICAL)
int *dens_rhoV;
int *dens_rhoVsq;
int *dens_rhoVsqtmp;
int *dens_nmeas;				// number of measurements at s
int *dens_p_nmeas;
int *dens_knmeas;		// number of measurements for rates kd/ka (meas makes only sense when cells present)
int *dens_nkd;			// number of divisions at s
int *dens_nka;			// number of deaths at s 
double *dens_ka;
double *dens_kd;	
double dens_binsize_z = 0.2;
double dens_binsize_r = 0.5;
int dens_length_z = (int)(LZ/dens_binsize_z+1);
int dens_length_r = (int)(LX/(2*dens_binsize_r)+1);
double dens_maxz = -1.;
double dens_maxr = -1.;
double dens_mean_maxz = 0.;
double dens_mean_maxr = 0.;
double dens_lastmean_maxz = 0.;
double dens_lastmean_maxr = 0.;
double dens_p_maxz = -1.;
double dens_p_maxr = -1;
double dens_mean_p_maxz = 0.;
double dens_mean_p_maxr = 0.;
double dens_lastmean_p_maxz = 0.;
double dens_lastmean_p_maxr = 0.;
double dens_minz = -1.;
double dens_mean_minz = 0.;
double dens_lastmean_minz = 0.;
double dens_p_minz = -1.;
double dens_mean_p_minz = 0.;
double dens_lastmean_p_minz = 0.;
double *dens_kperp;			// k_d perpendicular to surface
double *dens_kpara;			// k_d parallel to surface
int *dens_curn;				// current number of cells (used for kd/a)
int *dens_curnka;			// current number of cell deaths (used for ka)
int *dens_curnkd;			// current number of cell divisions (used for kd)
double *dens_sq;				// nematic order parameter of division (zz component)
double *dens_cell_sq;			// nematic order parameter of cell alignment (zz component)
double dens_com_rx;			// Center of mass r_x
double dens_com_ry;			// -"- r_y
double dens_com_rz;			// -"- r_z
long dens_com_n;
double *dens_mean_cell_dr;		// Mean cell "size", i.e. distance between its two particles
double *dens_mean_NN_dr;			// Mean nearest neighbour distance
 #endif
               
 #ifdef RECTANGULAR
double *flux_vx;				// x component of mean velocity in layer i
double *flux_vy;				// y component of mean velocity in layer i
double *flux_vz;				// z component of mean velocity in layer i
double *flux_p_vx;			// x component of mean particle velocity in layer i
double *flux_p_vy;			// y component of mean particle velocity in layer i
double *flux_p_vz;			// z component of mean particle velocity in layer i
 #elif defined (SPHERICAL)
double *flux_vr; 			// radial component of mean velocity in layer i
double *flux_p_vr;			// radial component of mean particle velocity in layer i
 #elif defined (ZYLINDRICAL)
double *flux_vr; 			// radial component of mean velocity in layer i
double *flux_p_vr;			// radial component of mean particle velocity in layer i
double *flux_vz; 			// radial component of mean velocity in layer i
double *flux_p_vz;			// radial component of mean particle velocity in layer i
 #endif
 #ifdef REALFLUX
long *rflux_jz;				// real flux in z direction
double *rflux_crz;			// old cell rz position
 #endif
#endif


#ifdef LOCALSTRESS
double stress_binsize[3] = {1., 1., 1.};
double stress_ibinsize[3] = {1./stress_binsize[0], 1./stress_binsize[1], 1./stress_binsize[2]};
int stress_lengthx = (int)(LX/stress_binsize[0]);
int stress_lengthy = (int)(LY/stress_binsize[1]);
int stress_lengthz = (int)(LZ/stress_binsize[2]);
int stress_length  = stress_lengthx*stress_lengthy*stress_lengthz*9;	// 9 because full tensor
int stress_lengthi[3] = {stress_lengthx, stress_lengthy, stress_lengthz};
double *stress_fcur;
double *stress_fcurd;
double **stress_fnext;
double **stress_fnextd;
double *stress_pcur;
double **stress_pnext;
double *stress_rx;
double *stress_ry;
double *stress_rz;
int stress_meas_enabled = 1;
double stress_meas_start = 0.;
double stress_meas_end = std::numeric_limits<double>::infinity();
/////////////////////////////////
// put things specific to either of these here
 #ifdef LOCALSTRESSMOP
 #endif
 #ifdef LOCALSTRESSVOL
 #endif
#endif


#ifdef FIXEDCELLS
double fixed_deltaz = 1.5;
const double fixed_init_displacement = 1.;
long fixed_wlist_max = 0;
long *fixed_wlist;
double *fixed_wlist_rx0;
double *fixed_wlist_ry0;
double *fixed_wlist_rz0;
double fixed_kspringx = 1000.;		// spring constant
double fixed_kspringy = 1000.;		// spring constant
double fixed_kspringz = 1000.;		// spring constant
double *fixed_pkz0;			// forces due to spring for fixed particles at z=0 
double *fixed_pkzL;			// forces due to spring for fixed particles at z=L
const int fixed_spec = 1;		// species of fixed cells
#ifdef DOUBLEFIXEDCELLS
double dfixed_lzpos = 10.;		// position of upper fixed cells layer
#endif
#endif


#ifdef PGAS
const int gas_spec = MAXSPECNUM-1;			// species number for gas particles
double gas_kt = PGAS_KT;
 #ifdef SPHERICAL
double gas_sph_rx0 = LX/2.;
double gas_sph_ry0 = LY/2.;
double gas_sph_rz0 = LZ/2.;
double gas_sph_k   = 1.0;
 #endif
#endif


#ifdef CONSTP_PBC
double constp_pbc_chi3 = 1.;			// rescaling factor chi^3 (for volume)
double constp_pbc_chi  = 1.;			// rescaling factor chi (for length)
double constp_pbc_p    = CONSTP_PBC_TARGET_P;	// Target pressure
double constp_pbc_curp = 0.;
double constp_pbc_fac  = 1.*TIMESTEP;
double constp_pbc_fxcur = 0.;
double *constp_pbc_fxnext;
double constp_pbc_fycur = 0.;
double *constp_pbc_fynext;
double constp_pbc_fzcur = 0.;
double *constp_pbc_fznext;
double constp_pbc_fxcurd = 0.;
double *constp_pbc_fxnextd;
double constp_pbc_fycurd = 0.;
double *constp_pbc_fynextd;
double constp_pbc_fzcurd = 0.;
double *constp_pbc_fznextd;
double constp_pbc_pxcur = 0.;
double constp_pbc_pycur = 0.;
double constp_pbc_pzcur = 0.;
#endif


#ifdef __DEBUG_TIMER
struct timeval tv[NUM_TIMERS];
unsigned long tv_usec[NUM_TIMERS];
unsigned long tv_count[NUM_TIMERS];

unsigned long __dtmp_s = 0;
unsigned long __dtmp_us = 0;
unsigned short tv_run[NUM_TIMERS];
unsigned int tv_used = 0;
char tv_desc[NUM_TIMERS][TIMER_DESC_SIZE];
#endif


