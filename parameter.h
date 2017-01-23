/*		parameter.h		*/
// define all paramter in here via #define

#ifndef __PARAMETER_H_
#define __PARAMETER_H_

#define TIMESTEP 0.001		// time step size (dt)

#define SIMTIME 40		// simulation time

#define INITSIZE 1048576	// initial number of possible cell particles
#define VLIST_MES 64		// verlet list mean entry size
#define VLISTINITSIZE VLIST_MES*INITSIZE	// initial vlist size

#define MAXSPECNUM 2		// maximum number of species
#define TARGETNUMTHREADS 16	// try to get this many threads


#define RPP 1.	 		// range of pair potential
#define R0  1.			// cellular growth expansion pressure constant II
#define RT  1.			// range of dissipative forces
#define RC  0.00001		// distance at which new cell particles are placed after division
#define RCT 0.8			// distance threshold for cell division

#define B 50    		// growth coefficient

#define MASS 1.			// mass of each particle

#define GAMMAC 100.
#define GAMMAT 50.
#define GAMMAB 0.1

#define KA 0.01			// rate of cell death

#define KBT 0.1			// noise intensity in the tissue

#define F0 2.39566029905288	// repulsive cell-cell potential coefficient
#define F1 7.5			// attractive cell-cell potential coefficient

#define BOXLENGTH_X 6		// box length in x direction
#define BOXLENGTH_Y 6
#define BOXLENGTH_Z 12		// for now only cubic boxes are supported


#define OUTDIR        "./"  		// output directory
#define FTRAJNAME     "traj.dat"        // file for trajectory output
#define FDIVTNAME     "divt.dat"	// file for status output
#define FDUMPNAME     "_dump.dat"	// save all data on shutdown in this file

#define FNUMCELLSNAME "numcells.dat"	// file for number of cells output

#define TRAJDATASTEP 1000		// output trajectory data every ... timesteps
#define SAVEDATASTEP 1000		// save state every ... timesteps
#define NUMCELLDATASTEP 100		// output number of cells every ... timesteps


#define RNGSEED get_seed()              // get_seed uses /dev/random (better than time())

#define EPSILON 1e-10

#define MAXTHREADNUM 16


////////////////////////////////////
// Check run time periodically and stop
// after WTRUNTIME seconds
//#define WALLTIMESTOP

#ifdef WALLTIMESTOP

 #define WTSTOPCHECKSTEP 100
 #define WTRUNTIME 12*3600       // run time in seconds

#endif


////////////////////////////////////
// to enable pressure measurement uncomment #define PMEASUREMENT
#define PMEASUREMENT

#ifdef PMEASUREMENT

#define PMEANSTEP 1000			// average pressure over PMEANSTEP steps
#define FPRESSNAME  "pressure.dat"	// file for pressure output

#endif

////////////////////////////////////
// enable debugging code
//#define __DEBUG
//#define __DEBUG_VLIST
//#define __DEBUG_VALIDATE_VLIST
//#define __DEBUG_TIMER



/////////////////////////////////////
// in which geometry would you like to
// measure (important for DENSITYPROFILE
// and PGAS)
//#define RECTANGULAR           // measure in rectangular geometry
//#define SPHERICAL               // measure in radial geometry
#define ZYLINDRICAL		// measure in zylindrical geometry 

////////////////////////////////////
// enable density profile, division rate and death rate measurement (be careful with PBC)
#define DENSITYPROFILE

#ifdef DENSITYPROFILE
 #define DENS_PAD 16		// number of additional bins
 #define DENSDATASTEP 1000
 #define FDENSITYNAME "density.dat"

 /////////////////////////////////////
 // Which coordinate system to use?
 #define MEASUREINR
 //#define MEASUREINS		// distance to (some) front


 ////////////////////////////////////
 // enable also real particle flux to be measured
 //#define REALFLUX

#endif


////////////////////////////////////
// enable local stress measurement (slow)
// for now makes only sense for RECTANGULAR
// IMPORTANT: Pressure is measured, not stress!!!!
// for stress multiply results with -1.
#define LOCALSTRESS

#ifdef LOCALSTRESS
 ///////////////////////////////////
 // either define VOL or MOP
 #define LOCALSTRESSVOL
 //#define LOCALSTRESSMOP

 #define LSTRESSDATASTEP 1000

 #define FLOCALSTRESSMOPNAME "stressmop.dat"
 #define FLOCALSTRESSVOLNAME "stressvol.dat"
#endif



////////////////////////////////////
// enable fixed cells
//#define FIXEDCELLS
#ifdef FIXEDCELLS			// per default only fixed cells are kept, others are discarded
 #define WLISTSIZE	4096		// maximum number of fixed cells (change according to needs)

 ////////////////////////////////////
 // enable loading of steady state fixed cells (fixed cells + normal cells are kept)
 //#define STEADYSTATEFIXEDCELLS

 ////////////////////////////////////
 // enable loading of double walls (also discards normal cells)
 //#define DOUBLEFIXEDCELLS
#endif


////////////////////////////////////
// enable gas particles
//#define PGAS

#ifdef PGAS
 #define PGAS_F0 0.1
 #define PGAS_KT 0.1
 #define PGAS_DX 0.6
 #define PGAS_DR 16.
#endif


////////////////////////////////////
// Constant pressure ensemble but only
// for PBC
//#define CONSTP_PBC

#ifdef CONSTP_PBC
 #define CONSTP_PBC_TARGET_P 0.
 #define CONSTP_PBC_STEP 1000
 #define FCONSTPPBCNAME "__curp.dat"
 //#define CONSTP_PBC_NORESCALE		// uncomment if you only want to measure mean pressure
#endif


////////////////////////////////////
// uncomment one of these to implement the according boundary condition

//#define PBC		// Periodic boundary condition
//#define BBC		// Bounce back boundary condition
//#define RBC		// Reflective boundary condition
#define HBC		// Hybrid boundary condition, adjust all HBC blocks to your specific needs
			// (normaly a combination of the 3 above)
			// default: periodic in x/y and bounce back in z direction




#include "mtrand.h"


//////////////////////////////////////////
// cell listing
extern long *_list;
extern long _head;
extern const int _pperc;			// particles per biological cell
						// DO NOT CHANGE UNLESS YOU KNOW WHAT YOU DO!!!

//////////////////////////////////////////
// cell properties
extern double *_rx;
extern double *_ry;
extern double *_rz;
extern double *_vx;
extern double *_vy;
extern double *_vz;
extern double *_fx;
extern double *_fy;
extern double *_fz;
extern double *_fdx;
extern double *_fdy;
extern double *_fdz;
extern long *_id;
extern int *_spec;				// cell interaction parameters can be different for different species

//////////////////////////////////////////
// global properties
extern double LX;				// box dimensions
extern double LY;
extern double LZ;
extern double curtime;				// current time index
extern double dt;               		// time step
extern double tend;				// simulation time
extern double kbt;

extern long next_free_id;			// give every cell a unique id.
						// in case of division each newly created cell particle gets new id
						// so that lineages can be restored from ids
extern long lastcn;				// size of _list and therefore (size of _r,...)/pperc
extern long lastn;				// size of arrays
extern long maxn;				// max size of arrays
extern MTRand *drand;				// random number generator handle (gaussuan


//////////////////////////////////////////
// species properties
extern double *rpp;
extern double *r0;
extern double *rt;
extern double *rc;
extern double *rct;
extern double *b;
extern double *gammac;
extern double *gammat;
extern double *gammab;
extern double *ka;
extern double *f0;
extern double *f1;
extern double *m;			// mass of particles


//////////////////////////////////////////
// other calculated properties of tissue
extern double com_rx;
extern double com_ry;
extern double com_rz;
extern long   com_n;


//////////////////////////////////////////
// other stuff
extern double max_ia_range;
extern int trajdatastep;
extern int savedatastep;
extern unsigned long rngseed;
extern double epsilon;

////////////////////////////////////////
// Verlet list
extern const double skin;			// max(skin*rpp, skin*rt) is verlet list cut off radius
extern double currc;				// store currently used min skin radius
extern long *vlist;				// verlet list
extern long *vlist_offset;			// offset vector
extern long *vlist_nc;				// vlist used entry length vector (nc=neighbour count)
extern double *vlist_dx;			// displacement in x,y,z direction since last verlet list update
extern double *vlist_dy;
extern double *vlist_dz;

////////////////////////////////////////
// special stuff
extern double *rna;				// random number array
extern unsigned char *rnai;     // save which fields are used in rna (much faster to zero out)
extern double *randf_prefc;			// random force prefactor
extern double *randf_preft;			// random force prefactor
extern int vlist_rebuild;			// rebuild vlist?

extern double boxlength[];


extern long t_begin;				// for actual run time measurement


//////////////////////////////////////////////////
// error value descriptions			//
//////////////////////////////////////////////////
extern long _errno;				// set error value (for now only used by load_dump())
extern short _int_flag;				// interrupt occured
extern int _int_no;				// interrupt number

extern const char *_errno_desc[];
extern const char *_errno_add_info;


#ifdef PMEASUREMENT
extern double *pm_px0;
extern double *pm_pxL;
extern double *pm_py0;
extern double *pm_pyL;
extern double *pm_pz0;
extern double *pm_pzL;
#endif


#ifdef DENSITYPROFILE
 #if (defined (RECTANGULAR) || defined (SPHERICAL))
extern int *dens_rhoV;				// density profile
extern int *dens_rhoVsq;
extern int *dens_rhoVsqtmp;
extern int *dens_nmeas;				// number of measurements at s
extern int *dens_p_nmeas;
extern int *dens_knmeas;		// number of measurements for rates kd/ka (meas makes only sense when cells present)
extern int *dens_nkd;				// number of divisions at s
extern int *dens_nka;				// number of deaths at s
extern double *dens_ka;
extern double *dens_kd;
extern double dens_binsize;
extern int dens_length;
extern double dens_maxz;
extern double dens_mean_maxz;
extern double dens_lastmean_maxz;
extern double dens_p_maxz;
extern double dens_mean_p_maxz;
extern double dens_lastmean_p_maxz;
extern double dens_minz;
extern double dens_mean_minz;
extern double dens_lastmean_minz;
extern double dens_p_minz;
extern double dens_mean_p_minz;
extern double dens_lastmean_p_minz;
extern double *dens_kperp;			// k_d perpendicular to surface
extern double *dens_kpara;			// k_d parallel to surface
extern int *dens_curn;				// current number of cells (used for kd/a)
extern int *dens_curnka;			// current number of cell deaths (used for ka)
extern int *dens_curnkd;			// current number of cell divisions (used for kd)
extern double *dens_sq;				// nematic order parameter of division (zz component)
extern double *dens_cell_sq;			// nematic order parameter of cell alignment (zz component)
extern double dens_com_rx;			// Center of mass r_x
extern double dens_com_ry;			// -"- r_y
extern double dens_com_rz;			// -"- r_z
extern long dens_com_n;
extern double *dens_mean_cell_dr;		// Mean cell "size", i.e. distance between its two particles
extern double *dens_mean_NN_dr;			// Mean nearest neighbour distance

 #elif defined (ZYLINDRICAL)
extern int *dens_rhoV;
extern int *dens_rhoVsq;
extern int *dens_rhoVsqtmp;
extern int *dens_nmeas;				// number of measurements at s
extern int *dens_p_nmeas;
extern int *dens_knmeas;		// number of measurements for rates kd/ka (meas makes only sense when cells present)
extern int *dens_nkd;			// number of divisions at s
extern int *dens_nka;			// number of deaths at s 
extern double *dens_ka;
extern double *dens_kd;	
extern double dens_binsize_z;
extern double dens_binsize_r;
extern int dens_length_z;
extern int dens_length_r;
extern double dens_maxz;
extern double dens_maxr;
extern double dens_mean_maxz;
extern double dens_mean_maxr;
extern double dens_lastmean_maxz;
extern double dens_lastmean_maxr;
extern double dens_p_maxz;
extern double dens_p_maxr;
extern double dens_mean_p_maxz;
extern double dens_mean_p_maxr;
extern double dens_lastmean_p_maxz;
extern double dens_lastmean_p_maxr;
extern double dens_minz;
extern double dens_mean_minz;
extern double dens_lastmean_minz;
extern double dens_p_minz;
extern double dens_mean_p_minz;
extern double dens_lastmean_p_minz;
extern double *dens_kperp;			// k_d perpendicular to surface
extern double *dens_kpara;			// k_d parallel to surface
extern int *dens_curn;		// current number of cells (used for kd/a)
extern int *dens_curnka;			// current number of cell deaths (used for ka)
extern int *dens_curnkd;			// current number of cell divisions (used for kd)
extern double *dens_sq;				// nematic order parameter of division (zz component)
extern double *dens_cell_sq;			// nematic order parameter of cell alignment (zz component)
extern double dens_com_rx;			// Center of mass r_x
extern double dens_com_ry;			// -"- r_y
extern double dens_com_rz;			// -"- r_z
extern long dens_com_n;
extern double *dens_mean_cell_dr;		// Mean cell "size", i.e. distance between its two particles
extern double *dens_mean_NN_dr;			// Mean nearest neighbour distance
 #endif
               
#ifdef RECTANGULAR
extern double *flux_vx;				// x component of mean velocity in layer i
extern double *flux_vy;				// y component of mean velocity in layer i
extern double *flux_vz;				// z component of mean velocity in layer i
extern double *flux_p_vx;			// x component of mean particle velocity in layer i
extern double *flux_p_vy;			// y component of mean particle velocity in layer i
extern double *flux_p_vz;			// z component of mean particle velocity in layer i
#elif defined (SPHERICAL)
extern double *flux_vr; 			// radial component of mean velocity in layer i
extern double *flux_p_vr;			// radial component of mean particle velocity in layer i
#elif defined ZYLINDRICAL
extern double *flux_vr; 			// radial component of mean velocity in layer i
extern double *flux_p_vr;			// radial component of mean particle velocity in layer i
extern double *flux_vz; 			// radial component of mean velocity in layer i
extern double *flux_p_vz;			// radial component of mean particle velocity in layer i
#endif
#ifdef REALFLUX
extern long *rflux_jz;				// real flux in z direction
extern double *rflux_crz;			// old cell rz position
#endif
#endif


#ifdef LOCALSTRESS
extern double stress_binsize[];
extern double stress_ibinsize[];
extern int stress_lengthx;
extern int stress_lengthy;
extern int stress_lengthz;
extern int stress_length;
extern int stress_lengthi[];
extern double *stress_fcur;
extern double *stress_fcurd;
extern double **stress_fnext;
extern double **stress_fnextd;
extern double *stress_pcur;
extern double **stress_pnext;
extern double *stress_rx;
extern double *stress_ry;
extern double *stress_rz;
extern int stress_meas_enabled;
extern double stress_meas_start;
extern double stress_meas_end;
 #ifdef LOCALSTRESSMOP
 #endif
 #ifdef LOCALSTRESSVOL
 #endif
#endif


#ifdef FIXEDCELLS
extern double fixed_deltaz;
extern const double fixed_init_displacement;
extern long fixed_wlist_max;
extern long *fixed_wlist;
extern long *fixed_wlist_ids;
extern double *fixed_wlist_rx0;
extern double *fixed_wlist_ry0;
extern double *fixed_wlist_rz0;
extern double fixed_kspringx;			// spring constant
extern double fixed_kspringy;			// spring constant
extern double fixed_kspringz;			// spring constant
extern double *fixed_pkz0;			// forces due to spring for fixed particles at z=0 
extern double *fixed_pkzL;			// forces due to spring for fixed particles at z=L
extern const int fixed_spec;			// species of fixed cells
#ifdef DOUBLEFIXEDCELLS
extern double dfixed_lzpos;			// position of upper fixed cells layer
#endif
#endif


#ifdef PGAS
extern const int gas_spec;
extern double gas_kt;
 #ifdef SPHERICAL
extern double gas_sph_rx0;
extern double gas_sph_ry0;
extern double gas_sph_rz0;
extern double gas_sph_k;
 #endif
#endif


#ifdef CONSTP_PBC
extern double constp_pbc_chi3;			// rescaling factor chi^3 (for volume)
extern double constp_pbc_chi;			// rescaling factor chi (for length)
extern double constp_pbc_p;			// Target pressure
extern double constp_pbc_curp;
extern double constp_pbc_fxcur;
extern double *constp_pbc_fxnext;
extern double constp_pbc_fycur;
extern double *constp_pbc_fynext;
extern double constp_pbc_fzcur;
extern double *constp_pbc_fznext;
extern double constp_pbc_fxcurd;
extern double *constp_pbc_fxnextd;
extern double constp_pbc_fycurd;
extern double *constp_pbc_fynextd;
extern double constp_pbc_fzcurd;
extern double *constp_pbc_fznextd;
extern double constp_pbc_pxcur;
extern double constp_pbc_pycur;
extern double constp_pbc_pzcur;
#endif




#ifdef __DEBUG_TIMER
#define NUM_TIMERS 24
#define TIMER_DESC_SIZE 64
extern struct timeval tv[NUM_TIMERS];
extern unsigned long tv_usec[NUM_TIMERS];
extern unsigned long tv_count[NUM_TIMERS];

extern unsigned long __dtmp_s;
extern unsigned long __dtmp_us;
extern unsigned short tv_run[NUM_TIMERS];
extern unsigned int tv_used;
extern char tv_desc[NUM_TIMERS][TIMER_DESC_SIZE];
#endif



///////////////////////////////////
// definition save guards
//   -> some defines are not
//      allowed at the same time
#if (defined (PBC) && (defined (BBC) || defined (RBC) || defined (HBC))) || \
    (defined (BBC) && (defined (PBC) || defined (RBC) || defined (HBC))) || \
    (defined (RBC) && (defined (PBC) || defined (BBC) || defined (HBC))) || \
    (defined (HBC) && (defined (PBC) || defined (BBC) || defined (RBC)))
  #error More than one boundary condition is enabled
#endif



#if (defined (RECTANGULAR) && defined (SPHERICAL))
  #error RECTANGULAR and SPHERICAL are mutually exclusive
#endif

#if (defined (RECTANGULAR) && defined (ZYLINDRICAL))
  #error RECTANGULAR and ZYLINDRICAL are mutually exclusive
#endif

#if (defined (ZYLINDRICAL) && defined (SPHERICAL))
  #error ZYLINDRICAL and SPHERICAL are mutually exclusive
#endif

#if (defined (MEASUREINS) && defined (MEASUREINR))
  #error MEASUREINS and MEASUREINR are mutually exclusive
#endif

#if (!defined (RECTANGULAR) && !defined (SPHERICAL) && !defined(ZYLINDRICAL))
  #error You must either specify RECTANGULAR or SPHERICAL or ZYLINDRICAL
#endif

#if (defined (LOCALSTRESS) && defined (PGAS))
  #error The results of LOCALSTRESS with PGAS are not known, disable this only if you really know what you are doing
#endif

#if (defined (LOCALSTRESSMOP) && defined (LOCALSTRESSVOL))
  #error The results of LOCALSTRESSMOP in combination with LOCALSTRESSVOL are unknown!!!
#endif

#if (defined (CONSTP_PBC) && (!defined (PBC)))
  #error We need PBC for CONSTP_PBC to work!
#endif

#if (defined (CONSTP_PBC) && (!defined (CONSTP_PBC_NORESCALE)) && defined (DENSITYPROFILE))
  #error DENSITYPROFILE does not work with rescaling CONSTP_PBC, because of fixed array size!
#endif

#if (defined (CONSTP_PBC) && (!defined (CONSTP_PBC_NORESCALE)) && defined (LOCALSTRESS))
  #error LOCALSTRESS does not work with rescaling CONSTP_PBC, because of fixed array size!
#endif

#if (defined (CONSTP_PBC) && (!defined (CONSTP_PBC_NORESCALE)) && defined (PGAS))
  #warning "I don't see a reason why rescaling CONSTP_PBC should not work with PGAS but I never tested it!"
#endif



// Actually this is not necessary, works with non cubic boxes as well
// WRONG: rescaling method requires cubic box!!! otherwise different chi
// for each spatial coordinate...
#if (defined (CONSTP_PBC) && !(BOXLENGTH_X==BOXLENGTH_Y && BOXLENGTH_X==BOXLENGTH_Z))
  #error For CONSTP_PBC we need cubic box!
#endif


#endif



