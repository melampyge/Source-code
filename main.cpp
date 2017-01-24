/*		main.cpp		*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <ctime>
#include <limits>

#ifdef __DEBUG_TIMER
#include <sys/time.h>
#endif

#include "parameter.h"      // parameter file
#include "bc_func.h"        // boundary conditions functions
#include "verlet.h"         // velocity verlet algorithm functions
#include "io_func.h"        // input/output functions
#include "basics.h"	    // particle properties macros
#include "divdeath.h"	    // cell division and death
#include "shutdown.h"	    // signal handling and shutdown procedures
#include "mtrand.h"	    // random number generation
#include "debug.h"	    // debugging functions


using namespace std;



///////////////////////////////////////
// function declarations
int initialize();



//////////////////////////////////////
// global variables

#include "parameter.cpp"



///////////////////////////////////////
// main()
int main(int argc, char *argv[])
{
  long i, j, k, l;
  ifstream iflock;		// to check whether lock file exists
  ofstream oflock;		// to create the lock file
  long count = 0;		// simulation step counter


  ///////////////////////////////
  // first check whether lock file exists -> simulation is already running
  iflock.open((string(OUTDIR)+string("lock.dat")).c_str());
  if(iflock) {
    cerr << "lock file exists -> simulation is possibly running! Exiting..." << endl;
    return (-1);
  }
  iflock.close();


  ///////////////////////////////////
  // set up signal handler
  init_signal_handler();


  /////////////////////////////////////
  // set number of threads
  omp_set_num_threads(TARGETNUMTHREADS);


#ifdef __DEBUG_TIMER
  INIT_TIMERS
  SET_TIMER(0, "main()")
  SET_TIMER(1, "verlet_step_dpd(): check_verlet_list()")
  SET_TIMER(2, "verlet_step_dpd(): calcv, calcr, check_bc")
  SET_TIMER(3, "verlet_step_dpd(): cell_death()")
  SET_TIMER(4, "verlet_step_dpd(): cell_division()")
  SET_TIMER(5, "verlet_step_dpd(): set rna to NaN")
  SET_TIMER(6, "verlet_step_dpd(): force loop 1")
  SET_TIMER(7, "verlet_step_dpd(): populate other half of rna")
  SET_TIMER(8, "verlet_step_dpd(): calcrandf, calcv new and set fd=0")
  SET_TIMER(9, "verlet_step_dpd(): update_verlet_list()")
  SET_TIMER(10, "verlet_step_dpd(): force loop 2")
  SET_TIMER(11, "main(): verlet_step_dpd()")
  SET_TIMER(12, "main(): output and measurements")

  SET_TIMER(13, "verlet_step_dpd(): update_verlet_list()")
  SET_TIMER(14, "verlet_step_dpd(): update_verlet_list()")
  SET_TIMER(15, "verlet_step_dpd(): update_verlet_list()")
  START_TIMER(0)
#endif


  ///////////////////////////////////////
  // allocate space for particle arrays and species properties
  if (!allocate_init_arrays()) {
    cerr << "allocate_init_arrays() error: " << _errno << endl;
    cerr << _errno_desc[_errno] << endl;
    cerr << _errno_add_info << endl;
  }


  /////////////////////////////////
  // if filename is given as argument -> load parameters and last data set from that file
  // else -> initialize from parameter.h
  if (argc > 1) {
    if (!load_dump(argv[1])) {
      cout << "Error while load_dump(): " << _errno
           << "\nError description: " << _errno_desc[_errno] << endl;
      return (-1);
    }
    cout << "load_dump(): " << _errno_desc[_errno] << endl;
    _errno = 0;
  }
  else {
    ///////////////////////////////////
    // initialize output
    init_output();

    t_begin = 0;

    ///////////////////////////////////
    // initialize r, v
    initialize();
  }


  /////////////////////////////////////
  // create lock file to avoid data corruption due to overlapping runs
  oflock.open((string(OUTDIR)+string("lock.dat")).c_str(), ios::trunc);
  if (!oflock) {
    cerr << "Error while creating lock file!" << endl;
    return (-1);
  }
  oflock.close();


  ////////////////////////////////////////
  // init verlet list
  update_verlet_list();


  if (argc < 2) {
    ///////////////////////////////////
    // initialize verlet, i.e. calculate forces
    init_verlet();

    //////////////////////////////////////
    // output data
    output_trajectory();
  }


  /////////////////////////////////////
  // get time for time measurement
  t_begin = time(0) - t_begin;


  /////////////////////////////////////
  // start velocity verlet algorithm //
  /////////////////////////////////////
  for (count = (long)(curtime/dt+1); count < _ceil(tend/dt); count++) {  // loop time steps
    curtime = count*dt;

//#ifdef __DEBUG
    if (count%1000==0) {
      cout << "t=" << curtime << ": lastcn=" << lastcn << ": lastn=" << lastn
           << ": next_free_id=" << next_free_id;
#ifdef WALLTIMESTOP
      cout << ": difft= " << WTRUNTIME-(time(0)-t_begin);
#endif
      cout << endl;
    }
//#endif

    /////////////////////////////////////
    // do one verlet step, i.e. calc forces, then new positions and velocities
    START_TIMER(11)
    verlet_step_dpd();
    STOP_TIMER(11)

    ////////////////////////////////////
    // if interrupt has occured during parallel execution, resolve it here
    resolve_int();

    START_TIMER(12)

    ////////////////////////////////////
    // output data every TRAJDATASTEPth step
    if ((count%TRAJDATASTEP) == 0)
      output_trajectory();

    ////////////////////////////////////
    // save state every SAVEDATASTEPth step
    if ((count%SAVEDATASTEP) == 0)
      dump_all();


    //////////////////////////////////////
    // conduct measurements of interest	//
    //////////////////////////////////////

    //////////////////////////////////////
    // calculate center of mass
    com_rx = com_ry = com_rz = 0.;
    com_n  = 0;
    for (i = 0; i < lastcn; ++i) {
      if (_list[i] != -1
#ifdef FIXEDCELLS
                         || _spec[i*_pperc] == fixed_spec
#endif
#ifdef PGAS
                         || _spec[i*_pperc] == gas_spec
#endif
                        ) continue;
      for (k=i*_pperc; k < (i+1)*_pperc; ++k) {
        com_rx += foldx(_rx[k]);
        com_ry += foldy(_ry[k]);
        com_rz += foldz(_rz[k]);
        com_n++;
      }
    }
    com_rx /= com_n;
    com_ry /= com_n;
    com_rz /= com_n;

#ifdef SPHERICAL
 #ifdef PGAS
    ////////////////////////////////////
    // fix COM by weak harmonic potential
#pragma omp parallel
{
#pragma omp for //schedule(auto)
    for (i = 0; i < lastcn; ++i)
      if (_list[i] == -1
  #ifdef FIXEDCELLS
          && _spec[i*_pperc] != fixed_spec
  #endif
  #ifdef PGAS
          && _spec[i*_pperc] != gas_spec
  #endif
          ) {
        double fxtmp = gas_sph_k*(gas_sph_rx0 - com_rx);
        double fytmp = gas_sph_k*(gas_sph_ry0 - com_ry);
        double fztmp = gas_sph_k*(gas_sph_rz0 - com_rz);
        _fx[i*_pperc] += fxtmp;
        _fy[i*_pperc] += fytmp;
        _fz[i*_pperc] += fztmp;
        _fx[i*_pperc+1] += fxtmp;
        _fy[i*_pperc+1] += fytmp;
        _fz[i*_pperc+1] += fztmp;
      }
}	// omp parallel end
 #endif
#endif


#ifdef PMEASUREMENT
    ///////////////////////////////////
    // average pressure over PMEANSTEP steps
    if ((count%PMEANSTEP) == 0) {
      output_pressure();

      fill(pm_px0, pm_px0+TARGETNUMTHREADS, 0.);
      fill(pm_pxL, pm_pxL+TARGETNUMTHREADS, 0.);
      fill(pm_py0, pm_py0+TARGETNUMTHREADS, 0.);
      fill(pm_pyL, pm_pyL+TARGETNUMTHREADS, 0.);
      fill(pm_pz0, pm_pz0+TARGETNUMTHREADS, 0.);
      fill(pm_pzL, pm_pzL+TARGETNUMTHREADS, 0.);
 #ifdef FIXEDCELLS
      fill(fixed_pkz0, fixed_pkz0+TARGETNUMTHREADS, 0.);
      fill(fixed_pkzL, fixed_pkzL+TARGETNUMTHREADS, 0.);
 #endif
    }
#endif


    ///////////////////////////////////
    // output number of cells
    if ((count%NUMCELLDATASTEP) == 0)
      output_numcells();

#ifdef DENSITYPROFILE
    ///////////////////////////////////
    // cell_death() ensures correct value for dens_maxz
    // !!!! no it does not, as particles move after cell_death, so calc new maxz
    dens_maxz = dens_p_maxz = -1.;
    dens_minz = dens_p_minz = LZ+1.;
    dens_com_rx = com_rx;
    dens_com_ry = com_ry;
    dens_com_rz = com_rz;
 #ifdef ZYLINDRICAL
    dens_maxr = dens_p_maxr = -1.;
 #endif
 #ifdef RECTANGULAR
    for (i = 0; i < lastcn; ++i) {
      if (_list[i] != -1
  #ifdef FIXEDCELLS
                          || _spec[i*_pperc] == fixed_spec
  #endif
  #ifdef PGAS
                          || _spec[i*_pperc] == gas_spec
  #endif
                        ) continue;
      if ((_rz[i*_pperc]+_rz[i*_pperc+1])/2. > dens_maxz)
        dens_maxz = (_rz[i*_pperc]+_rz[i*_pperc+1])/2.;
      if ((_rz[i*_pperc]+_rz[i*_pperc+1])/2. < dens_minz)
        dens_minz = (_rz[i*_pperc]+_rz[i*_pperc+1])/2.;
      for (k=i*_pperc; k < (i+1)*_pperc; ++k) {
        if (_rz[k] > dens_p_maxz)
          dens_p_maxz = _rz[k];
        if (_rz[k] < dens_p_minz)
          dens_p_minz = _rz[k];
      }
    }
 #elif defined (SPHERICAL)
    for (i = 0; i < lastcn; ++i) {
      if (_list[i] != -1
  #ifdef FIXEDCELLS
                         || _spec[i*_pperc] == fixed_spec
  #endif
  #ifdef PGAS
                         || _spec[i*_pperc] == gas_spec
  #endif
                        ) continue;
      double tmpcr[3] = {0., 0., 0.};
      for (k=i*_pperc; k < (i+1)*_pperc; ++k) {
        double tmppr[3] = { _rx[k]-dens_com_rx, _ry[k]-dens_com_ry, _rz[k]-dens_com_rz };
        double tmpprsq = NORMSQ(tmppr[0], tmppr[1], tmppr[2]);
        if (tmpprsq > dens_p_maxz*dens_p_maxz)
          dens_p_maxz = sqrt(tmpprsq);
        tmpcr[0] += tmppr[0];
        tmpcr[1] += tmppr[1];
        tmpcr[2] += tmppr[2];
      }
      double tmpcrsq = NORMSQ(tmpcr[0]/_pperc, tmpcr[1]/_pperc, tmpcr[2]/_pperc);
      if (tmpcrsq > dens_maxz*dens_maxz)
        dens_maxz = sqrt(tmpcrsq);
    }
 #elif defined (ZYLINDRICAL)
    for (i = 0; i < lastcn; ++i) {
      if (_list[i] != -1
  #ifdef FIXEDCELLS
                         || _spec[i*_pperc] == fixed_spec
  #endif
  #ifdef PGAS
                         || _spec[i*_pperc] == gas_spec
  #endif
                        ) continue;
      if ((_rz[i*_pperc]+_rz[i*_pperc+1])/2. > dens_maxz)
        dens_maxz = (_rz[i*_pperc]+_rz[i*_pperc+1])/2.;
      if ((NORM(_rx[i*_pperc]-dens_com_rx,_ry[i*_pperc]-dens_com_ry, 0.) + NORM(_rx[i*_pperc+1]-dens_com_rx,_ry[i*_pperc+1]-dens_com_ry, 0.))/2. > dens_maxr) {
	  if((NORM(_rx[i*_pperc]-dens_com_rx,_ry[i*_pperc]-dens_com_ry, 0.) + NORM(_rx[i*_pperc+1]-dens_com_rx,_ry[i*_pperc+1]-dens_com_ry, 0.))/2. < BOXLENGTH_X/2)
	    dens_maxr = (NORM(_rx[i*_pperc]-dens_com_rx,_ry[i*_pperc]-dens_com_ry, 0.) + NORM(_rx[i*_pperc+1]-dens_com_rx,_ry[i*_pperc+1]-dens_com_ry, 0.))/2.;
      }
      if ((_rz[i*_pperc]+_rz[i*_pperc+1])/2. < dens_minz)
        dens_minz = (_rz[i*_pperc]+_rz[i*_pperc+1])/2.;
      for (k=i*_pperc; k < (i+1)*_pperc; ++k) {
        if (_rz[k] > dens_p_maxz)
          dens_p_maxz = _rz[k];
        if (_rz[k] < dens_p_minz)
          dens_p_minz = _rz[k];
	if (NORM(_rx[k]-dens_com_rx, _ry[k]-dens_com_ry,0.) > dens_p_maxr)
	  dens_p_maxr = NORM(_rx[k]-dens_com_rx, _ry[k]-dens_com_ry,0.);
      }
    }
 #endif		// ifdef REC elifdef SPH elifdef ZYL

 #if ( defined (RECTANGULAR) || defined (SPHERICAL) )
     fill(dens_rhoVsqtmp, dens_rhoVsqtmp+dens_length, 0);
 #elif defined (ZYLINDRICAL)
     fill(dens_rhoVsqtmp, dens_rhoVsqtmp+dens_length_r*dens_length_z, 0);
 #endif
    for (i = 0; i < lastcn; ++i) {
      if (_list[i] != -1
 #ifdef FIXEDCELLS
                         || _spec[i*_pperc] == fixed_spec
 #endif
 #ifdef PGAS
                         || _spec[i*_pperc] == gas_spec
 #endif
                        ) continue;
 #if ( defined (RECTANGULAR) || defined (SPHERICAL) )
      int __index = GET__INDEX(i, dens_binsize);		// see basics.h
      __ASSERT(__index, 0, dens_length, "__index in main()");
 #elif defined (ZYLINDRICAL)
       int __index = GET__INDEX_R(i, dens_binsize_r) + GET__INDEX_Z(i, dens_binsize_z)*dens_length_r;	// see basics.h
      __ASSERT(__index, 0, dens_length_z*dens_length_r, "__index in main()");
 #endif
      dens_rhoV[__index]++;
      dens_rhoVsqtmp[__index]++;
      {
        double tmprx  = (_rx[i*_pperc]-_rx[i*_pperc+1])/2.;
        double tmpry  = (_ry[i*_pperc]-_ry[i*_pperc+1])/2.;
        double tmprz  = (_rz[i*_pperc]-_rz[i*_pperc+1])/2.;
        double tmprpx = (_rx[i*_pperc]+_rx[i*_pperc+1])/2.;
        double tmprpy = (_ry[i*_pperc]+_ry[i*_pperc+1])/2.;
        double tmprpz = (_rz[i*_pperc]+_rz[i*_pperc+1])/2.;
        double tmp    = CALC_THETA(tmprx, tmpry, tmprz, dens_com_rx-tmprpx, dens_com_ry-tmprpy, dens_com_rz-tmprpz);
        ///////////////////////////////////////////
        // produces nan for only 1 existing cell //
        ///////////////////////////////////////////
        dens_cell_sq[__index] += 3./2.*tmp*tmp-0.5;
        ///////////////////////////////////////////
        // Calc mean cell size
        dens_mean_cell_dr[__index] += NORM(2.*tmprx, 2.*tmpry, 2.*tmprz);

        //////////////////////////////////////
        // Calc nearest neighbour distance
        double mindrsq = LZ*LZ;
        for (j = vlist_offset[i*_pperc]; j < vlist_nc[i*_pperc]; ++j) {
          int j1 = vlist[j];
          int j2 = __GETOTHERCELLPART(j1);
          if (j1 > j2)
            continue;
  #ifdef PGAS
          if (_spec[j1] == gas_spec)
            continue;
  #endif
          double tmprpx2 = (_rx[j1]+_rx[j1+1])/2.;
          double tmprpy2 = (_ry[j1]+_ry[j1+1])/2.;
          double tmprpz2 = (_rz[j1]+_rz[j1+1])/2.;
          double drsq    = NORMSQ(tmprpx-tmprpx2, tmprpy-tmprpy2, tmprpz-tmprpz2);
          if (drsq < mindrsq)
            mindrsq = drsq;
        }
        dens_mean_NN_dr[__index] += sqrt(mindrsq);
      }
 
 #ifdef RECTANGULAR
      flux_vx[__index]   += (_vx[i*_pperc]+_vx[i*_pperc+1])/2.;
      flux_vy[__index]   += (_vy[i*_pperc]+_vy[i*_pperc+1])/2.;
      flux_vz[__index]   += (_vz[i*_pperc]+_vz[i*_pperc+1])/2.;
 #elif defined SPHERICAL
      {
        double tmpvx = (_vx[i*_pperc]+_vx[i*_pperc+1])/2.;
        double tmpvy = (_vy[i*_pperc]+_vy[i*_pperc+1])/2.;
        double tmpvz = (_vz[i*_pperc]+_vz[i*_pperc+1])/2.;
        double tmprx = (_rx[i*_pperc]+_rx[i*_pperc+1])/2.;
        double tmpry = (_ry[i*_pperc]+_ry[i*_pperc+1])/2.;
        double tmprz = (_rz[i*_pperc]+_rz[i*_pperc+1])/2.;
        double tmpernorm = NORM(dens_com_rx-tmprx, dens_com_ry-tmpry, dens_com_rz-tmprz);
        flux_vr[__index] += (tmpvx*(tmprx-dens_com_rx) + tmpvy*(tmpry-dens_com_ry)
                             + tmpvz*(tmprz-dens_com_rz))/tmpernorm;
      }
 #elif defined ZYLINDRICAL
      {
	double tmpvx = (_vx[i*_pperc]+_vx[i*_pperc+1])/2.;
	double tmpvy = (_vy[i*_pperc]+_vy[i*_pperc+1])/2.;
	double tmprx = (_rx[i*_pperc]+_rx[i*_pperc+1])/2.;
	double tmpry = (_ry[i*_pperc]+_ry[i*_pperc+1])/2.;
	double tmpernorm = NORM(dens_com_rx-tmprx, dens_com_ry-tmpry, 0.);
	flux_vr[__index] +=  (tmpvx*(tmprx-dens_com_rx) + tmpvy*(tmpry-dens_com_ry));
	flux_vz[__index]   += (_vz[i*_pperc]+_vz[i*_pperc+1])/2.;
      }
 #endif  // ifdef RECTANGULAR
 #if ( defined(RECTANGULAR) || defined (SPHERICAL) )
      for (k = i*_pperc; k < (i+1)*_pperc; ++k) {
        int __index_p = GET__INDEX_P(k, dens_binsize);		// see basics.h
        __ASSERT(__index_p, 0, dens_length, "__index_p in main()");
  #ifdef RECTANGULAR
        flux_p_vx[__index_p] += _vx[k];
        flux_p_vy[__index_p] += _vy[k];
        flux_p_vz[__index_p] += _vz[k];
  #elif defined SPHERICAL
        {
          double tmpernorm = NORM(dens_com_rx-_rx[k], dens_com_ry-_ry[k], dens_com_rz-_rz[k]);
          flux_vr[__index] += (_vx[k]*(_rx[k]-dens_com_rx) + _vy[k]*(_ry[k]-dens_com_ry)
                             + _vz[k]*(_rz[k]-dens_com_rz))/tmpernorm;
        }
  #endif
      }
 #elif defined ZYLINDRICAL
       for (k = i*_pperc; k < (i+1)*_pperc; ++k) {
	    int __index_p = GET__INDEX_P_R(k, dens_binsize_r) + GET__INDEX_P_Z(k, dens_binsize_z)*dens_length_r;	// see basics.h
	    __ASSERT(__index_p, 0, dens_length_r*dens_length_z, "__index_p in main()");
	    
	    double tmpernorm = NORM(dens_com_rx-_rx[k], dens_com_ry-_ry[k], 0.);
	    flux_p_vr[__index_p] += (_vx[k]*(_rx[k]-dens_com_rx) + _vy[k]*(_ry[k]-dens_com_ry))/tmpernorm;
	    flux_p_vz[__index_p] += _vz[k];
	}
 #endif  // ifdef RECTANGULAR

 #ifdef REALFLUX
      double crz = GET_REALFLUX_CRZ(i);
      int fbinold = (int)(rflux_crz[i]/dens_binsize);
      int fbinnew = (int)(crz/dens_binsize);

      if (fbinnew-fbinold > 0)
        rflux_jz[fbinnew]++;
      if (fbinnew-fbinold < 0)
        rflux_jz[fbinold]--;

      rflux_crz[i] = crz;
 #endif
    }
    ///////////////////////////////////
    // calculate rhoVsq
 #if ( defined(RECTANGULAR) || defined (SPHERICAL) )
    for (i = 0; i < dens_maxz/dens_binsize+1; ++i)
      dens_rhoVsq[i] += dens_rhoVsqtmp[i]*dens_rhoVsqtmp[i];

    ///////////////////////////////////
    // count which z values have been reached and measured
    for (i = 0; i <= dens_maxz/dens_binsize+1; ++i)
      dens_nmeas[i]++;
    for (i = 0; i <= dens_p_maxz/dens_binsize+1; ++i)
      dens_p_nmeas[i] += 2;


    dens_mean_maxz += dens_maxz;
    dens_mean_p_maxz += dens_p_maxz;
    dens_mean_minz += dens_minz;
    dens_mean_p_minz += dens_p_minz;
    if ((count%DENSDATASTEP) == 0) {
      output_density();
      memset(dens_rhoV,  0x00, sizeof(int)*dens_length);
      memset(dens_rhoVsq,  0x00, sizeof(int)*dens_length);
      memset(dens_rhoVsqtmp,  0x00, sizeof(int)*dens_length);
      memset(dens_nmeas, 0x00, sizeof(int)*dens_length);
      memset(dens_knmeas, 0x00, sizeof(int)*dens_length);
      memset(dens_nkd,   0x00, sizeof(int)*dens_length);
      memset(dens_nka,   0x00, sizeof(int)*dens_length);
      fill(dens_kpara, dens_kpara+dens_length, 0.);
      fill(dens_kperp, dens_kperp+dens_length, 0.);
      fill(dens_ka, dens_ka+dens_length, 0.);
      fill(dens_kd, dens_kd+dens_length, 0.);
      fill(dens_p_nmeas, dens_p_nmeas+dens_length, 0);
      fill(dens_sq, dens_sq+dens_length, 0.);
      fill(dens_cell_sq, dens_cell_sq+dens_length, 0.);
      fill(dens_mean_cell_dr, dens_mean_cell_dr+dens_length, 0.);
      fill(dens_mean_NN_dr, dens_mean_NN_dr+dens_length, 0.);
      dens_lastmean_maxz = dens_mean_maxz;
      dens_mean_maxz = 0.;
      dens_lastmean_p_maxz = dens_mean_p_maxz;
      dens_mean_p_maxz = 0.;
      dens_lastmean_minz = dens_mean_minz;
      dens_mean_minz = 0.;
      dens_lastmean_p_minz = dens_mean_p_minz;
      dens_mean_p_minz = 0.;
 #elif defined (ZYLINDRICAL)
    for (i = 0; i <= dens_maxz/dens_binsize_z; ++i) {
	for(j = 0; j <= dens_maxr/dens_binsize_r; ++j) {
	    dens_rhoVsq[j + i*dens_length_r] += dens_rhoVsqtmp[j + i*dens_length_r]*dens_rhoVsqtmp[j + i*dens_length_r];
	}
    }
    ///////////////////////////////////
    // count which z and r values have been reached and measured
    for (i = 0; i <= dens_maxz/dens_binsize_z; ++i) {
	for(j = 0; j <= dens_maxr/dens_binsize_r; ++j) {
	    dens_nmeas[j + i*dens_length_r]++;
	}
    }
    for (i = 0; i <= dens_p_maxz/dens_binsize_z; ++i) {
	for(j = 0; j <= dens_p_maxr/dens_binsize_r; ++j) {
	    dens_p_nmeas[j + i*dens_length_r] += 2;
	}
    }

    dens_mean_maxz += dens_maxz;
    dens_mean_maxr += dens_maxr;
    dens_mean_p_maxz += dens_p_maxz;
    dens_mean_p_maxr += dens_p_maxr;
    dens_mean_minz += dens_minz;
    dens_mean_p_minz += dens_p_minz;
    if ((count%DENSDATASTEP) == 0) {
    output_density();
    memset(dens_rhoV,  0x00, sizeof(int)*dens_length_r*dens_length_z);
    memset(dens_rhoVsq,  0x00, sizeof(int)*dens_length_r*dens_length_z);
    memset(dens_rhoVsqtmp,  0x00, sizeof(int)*dens_length_r*dens_length_z);
    memset(dens_nmeas, 0x00, sizeof(int)*dens_length_r*dens_length_z);
    memset(dens_knmeas, 0x00, sizeof(int)*dens_length_r*dens_length_z);
    memset(dens_nkd,   0x00, sizeof(int)*dens_length_r*dens_length_z);
    memset(dens_nka,   0x00, sizeof(int)*dens_length_r*dens_length_z);
    fill(dens_kpara, dens_kpara+dens_length_r*dens_length_z, 0.);
    fill(dens_kperp, dens_kperp+dens_length_r*dens_length_z, 0.);
    fill(dens_ka, dens_ka+dens_length_r*dens_length_z, 0.);
    fill(dens_kd, dens_kd+dens_length_r*dens_length_z, 0.);
    fill(dens_p_nmeas, dens_p_nmeas+dens_length_r*dens_length_z, 0);
    fill(dens_sq, dens_sq+dens_length_r*dens_length_z, 0.);
    fill(dens_cell_sq, dens_cell_sq+dens_length_r*dens_length_z, 0.);
    fill(dens_mean_cell_dr, dens_mean_cell_dr+dens_length_r*dens_length_z, 0.);
    fill(dens_mean_NN_dr, dens_mean_NN_dr+dens_length_r*dens_length_z, 0.);
    dens_lastmean_maxz = dens_mean_maxz;
    dens_lastmean_maxr = dens_mean_maxr;
    dens_mean_maxz = 0.;
    dens_mean_maxr = 0.;
    dens_lastmean_p_maxz = dens_mean_p_maxz;
    dens_lastmean_p_maxr = dens_mean_p_maxr;
    dens_mean_p_maxz = 0.;
    dens_mean_p_maxr = 0.;
    dens_lastmean_minz = dens_mean_minz;
    dens_mean_minz = 0.;
    dens_lastmean_p_minz = dens_mean_p_minz;
    dens_mean_p_minz = 0.;
 #endif
 #ifdef RECTANGULAR
      fill(flux_vx, flux_vx+dens_length, 0.);
      fill(flux_vy, flux_vy+dens_length, 0.);
      fill(flux_vz, flux_vz+dens_length, 0.);
      fill(flux_p_vx, flux_p_vx+dens_length, 0.);
      fill(flux_p_vy, flux_p_vy+dens_length, 0.);
      fill(flux_p_vz, flux_p_vz+dens_length, 0.);
 #elif defined (SPHERICAL)
      fill(flux_vr, flux_vr+dens_length, 0.);
      fill(flux_p_vr, flux_p_vr+dens_length, 0.);
 #elif defined (ZYLINDRICAL)
      fill(flux_vz, flux_vz+dens_length_r*dens_length_z, 0.);
      fill(flux_vr, flux_vr+dens_length_r*dens_length_z, 0.);
      fill(flux_p_vz, flux_p_vz+dens_length_r*dens_length_z, 0.);
      fill(flux_p_vr, flux_p_vr+dens_length_r*dens_length_z, 0.);
 #endif
 #ifdef REALFLUX
      fill(rflux_jz, rflux_jz+dens_length, 0);
 #endif
 #ifdef ORDERINGTENSOR
      output_odtens();
      fill(odtens_Q, odtens_Q+odtens_length*6, 0.);
      fill(odtens_nmeas, odtens_nmeas+odtens_length, 0.);
 #endif
    }
#endif		// ifdef DENSITYPROFILE


    ///////////////////////////////////////////////
    // local stress measurement
#ifdef LOCALSTRESS
    if (!stress_meas_enabled && curtime >= stress_meas_start && curtime < stress_meas_end)
      stress_meas_enabled = 1;
    if (stress_meas_enabled) {
      if (curtime > stress_meas_end)
        stress_meas_enabled = 0;
      for (k = 0; k < TARGETNUMTHREADS; ++k)
        for (l = 0; l < stress_length; ++l)
          stress_pcur[l] += stress_pnext[k][l];

      for (i = 0; i < lastcn; ++i) {
        if (_list[i] != -1
 #ifdef PGAS
                           || _spec[i*_pperc] == gas_spec
 #endif
                          ) continue;
        for (k = i*_pperc; k < (i+1)*_pperc; ++k) {
          stress_rx[k] = _rx[k];
          stress_ry[k] = _ry[k];
          stress_rz[k] = _rz[k];
        }
      }

      if ((count%LSTRESSDATASTEP) == 0) {
        output_localstress();

        fill (stress_fcur, stress_fcur+stress_length, 0.);
        fill (stress_fcurd, stress_fcurd+stress_length, 0.);
        fill (stress_pcur, stress_pcur+stress_length, 0.);
      }

      for (k = 0; k < TARGETNUMTHREADS; ++k)
        for (l = 0; l < stress_length; ++l) {
          stress_fcur[l]  += stress_fnext[k][l];
          stress_fcurd[l] += stress_fnextd[k][l];
        }

#pragma omp parallel
{
      int tid = omp_get_thread_num();
      fill (stress_fnext[tid], stress_fnext[tid]+stress_length, 0.);
      fill (stress_fnextd[tid], stress_fnextd[tid]+stress_length, 0.);
      fill (stress_pnext[tid], stress_pnext[tid]+stress_length, 0.);
}   // parallel end
    }
#endif



    ////////////////////////////////////////////////
    // Calc p=-sigma for CONSTP_PBC
    // WARNING: Right now not correct, f is from t+1 and p is from t
#ifdef CONSTP_PBC
    for (i = 0; i < lastcn; ++i) {
      if (_list[i] != -1
 #ifdef PGAS
                         || _spec[i*_pperc] == gas_spec
 #endif
                        ) continue;
      for (k = i*_pperc; k < (i+1)*_pperc; ++k) {
        constp_pbc_pxcur += m[_spec[k]]*_vx[k]*_vx[k]*dt;	// + because p=-sigma and we want p
        constp_pbc_pycur += m[_spec[k]]*_vy[k]*_vy[k]*dt;	// not sigma
        constp_pbc_pzcur += m[_spec[k]]*_vz[k]*_vz[k]*dt;
      }
    }

    if ((count%CONSTP_PBC_STEP)==0) {
      constp_pbc_curp  = constp_pbc_fxcur+constp_pbc_fycur+constp_pbc_fzcur;
      constp_pbc_curp += constp_pbc_fxcurd+constp_pbc_fycurd+constp_pbc_fzcurd;
      constp_pbc_curp += (constp_pbc_pxcur+constp_pbc_pycur+constp_pbc_pzcur)/dt;
      constp_pbc_curp /= 3.*CONSTP_PBC_STEP*LX*LY*LZ;		// /3. because P=-1/3.*(s_xx+s_yy+s_zz)
      constp_pbc_chi3  = 1. - constp_pbc_fac*(constp_pbc_p-constp_pbc_curp);
      constp_pbc_chi   = pow(constp_pbc_chi3, 1./3.);

      output_curp();

      ////////////////////////////////////
      // Rescale positions
      // Actually everything would have to
      // be rescaled here (_f[d]{xyz}, _v{xyz}, constp_pbc_f{xyz},...) ???
 #ifndef CONSTP_PBC_NORESCALE
      for (i = 0; i < lastcn; ++i) {
        if (_list[i] != -1
                          ) continue;
        k = i*_pperc;
  #ifdef PGAS
        if (_spec[i*_pperc] == gas_spec) {
          _rx[k] *= constp_pbc_chi;
          _ry[k] *= constp_pbc_chi;
          _rz[k] *= constp_pbc_chi;
        }
        else {
  #endif
          double rxc = (_rx[k]+_rx[k+1])/2.;
          double ryc = (_ry[k]+_ry[k+1])/2.;
          double rzc = (_rz[k]+_rz[k+1])/2.;
          _rx[k]   += rxc*(constp_pbc_chi-1.);
          _rx[k+1] += rxc*(constp_pbc_chi-1.);
          _ry[k]   += ryc*(constp_pbc_chi-1.);
          _ry[k+1] += ryc*(constp_pbc_chi-1.);
          _rz[k]   += rzc*(constp_pbc_chi-1.);
          _rz[k+1] += rzc*(constp_pbc_chi-1.);
  #ifdef PGAS
        }
  #endif
      }

      /////////////////////////////////////
      // Rescale box size!
      LX *= constp_pbc_chi;
      LY *= constp_pbc_chi;
      LZ *= constp_pbc_chi;
 #endif		// CONSTP_PBC_NORESCALE


      constp_pbc_curp   = 0.;
      constp_pbc_pxcur  = 0.;
      constp_pbc_pycur  = 0.;
      constp_pbc_pzcur  = 0.;
      constp_pbc_fxcur  = 0.;
      constp_pbc_fycur  = 0.;
      constp_pbc_fzcur  = 0.;
      constp_pbc_fxcurd = 0.;
      constp_pbc_fycurd = 0.;
      constp_pbc_fzcurd = 0.;
    }

    for (i = 0; i < TARGETNUMTHREADS; ++i) {
      constp_pbc_fxcur  += constp_pbc_fxnext[i];
      constp_pbc_fxcurd += constp_pbc_fxnextd[i];
      constp_pbc_fycur  += constp_pbc_fynext[i];
      constp_pbc_fycurd += constp_pbc_fynextd[i];
      constp_pbc_fzcur  += constp_pbc_fznext[i];
      constp_pbc_fzcurd += constp_pbc_fznextd[i];
      constp_pbc_fxnext[i] = 0.; constp_pbc_fxnextd[i] = 0.;
      constp_pbc_fynext[i] = 0.; constp_pbc_fynextd[i] = 0.;
      constp_pbc_fznext[i] = 0.; constp_pbc_fznextd[i] = 0.;
    }
#endif



    STOP_TIMER(12)

#ifdef __DEBUG_TIMER
    if ((count%500) == 0)
      OUTPUT_TIMERDATA((string(OUTDIR)+string("__timing.dat")).c_str());
#endif

#ifdef WALLTIMESTOP
    if ((count%WTSTOPCHECKSTEP) == 0)
      if ((time(0)-t_begin) > WTRUNTIME) {
        process_signal_handler(0xF000);
        resolve_int();
      }
#endif

  }  // end loop time steps


#ifdef __DEBUG_TIMER
  STOP_TIMER(0)

  ofstream fout((string(OUTDIR)+string("__timing.dat")).c_str(), ios::trunc);
  if (!fout) {
    cerr << "Error while fopen('__timing.dat')" << endl;
    return (-1);
  }
  for (unsigned int __i = 0; __i < tv_used; ++__i) {
    STOP_TIMER(__i)
    fout << __i << "\t" << setw(16) << tv_usec[__i] << setw(12) << tv_count[__i]
         << setw(12) << ((double)tv_usec[__i])/tv_count[__i] << setw(12) << (double)tv_usec[__i]/tv_usec[0]*100.
         << "\t" << tv_desc[__i] << endl;
  }
  fout.close();
#endif


  long t_diff = time(0) - t_begin;
  long t_s = t_diff % 60;
  long t_m = (t_diff/60) % 60;
  long t_h = (t_diff/3600) % 24;
  long t_d = t_diff/(3600*24);
  cout << "Process ran for >> " << t_d << "d " << setw(2) << setfill('0') << t_h << ":"
       << setw(2) << setfill('0') << t_m << ":"
       << setw(2) << setfill('0') << t_s << " << " << t_diff << "s" << endl;

  ///////////////////////////////////
  // delete lock file
  unlink((string(OUTDIR)+string("lock.dat")).c_str());


  return 0;
}


//////////////////////////////////////
// function definitions             //
//////////////////////////////////////


// initialize(...)
// ---------------
// set initial positions and velocities
// ===
// Parameters:
int initialize()
{
  /*int i, j, k;
  int count;
  int numPperD;*/
  double tmpx, tmpy, tmpz;



  {
    _rx[0] = LX/2.;
    _ry[0] = LY/2.;
    _rz[0] = LZ/2.;
    double initdist = 0.4;
    tmpx = (drand[0]()*initdist*2.-initdist)/sqrt(3);
    tmpy = (drand[0]()*initdist*2.-initdist)/sqrt(3);
    tmpz = sqrt(initdist*initdist-tmpx*tmpx-tmpy*tmpy);
    _rx[1] = LX/2. + tmpx;
    _ry[1] = LY/2. + tmpy;
    _rz[1] = LZ/2. + tmpz;

    _vx[0] = _vy[0] = _vz[0] = 0.;
    _vx[1] = _vy[1] = _vz[1] = 0.;

    _id[0] = next_free_id++;
    _id[1] = next_free_id++;

    _spec[0] = _spec[1] = 0;

    _list[0] = -1;
    _head    = 1;
    lastn    = 2;
    lastcn   = 1;
  }//*/



  /*double tx, ty, tz;
  int i = 0;
  int j;
  double initdist = 0.2;		// ensure that inidist < dxyz, otherwise particles could be placed outside!!!
  lastn = lastcn = 0;
  const double dxyz = 0.6;
  for (tx = dxyz; tx < LX-dxyz; tx += dxyz) {
    for (ty = dxyz; ty < LY-dxyz; ty += dxyz) {
      for (tz = dxyz; tz < LZ-dxyz; tz += dxyz) {
        j = i*_pperc;
        _rx[j] = tx;
        _ry[j] = ty;
        _rz[j] = tz;
        _id[j] = next_free_id++;
        _spec[j] = 0;
        _vx[j] = get_ran_gaussian(0, 0., sqrt(kbt));
        _vy[j] = get_ran_gaussian(0, 0., sqrt(kbt));
        _vz[j] = get_ran_gaussian(0, 0., sqrt(kbt));

        j++;
        double tdr = initdist*drand[0]();
        double tdtheta = 3.1415926536*drand[0]();
        double tdphi = 2.*3.1415926536*drand[0]();
        double tdx = tdr*sin(tdtheta)*cos(tdphi);
        double tdy = tdr*sin(tdtheta)*sin(tdphi);
        double tdz = tdr*cos(tdtheta);
        _rx[j] = tx+tdx;
        _ry[j] = ty+tdy;
        _rz[j] = tz+tdz;
        _id[j] = next_free_id++;
        _spec[j] = 0;
        _vx[j] = get_ran_gaussian(0, 0., sqrt(kbt));
        _vy[j] = get_ran_gaussian(0, 0., sqrt(kbt));
        _vz[j] = get_ran_gaussian(0, 0., sqrt(kbt));

        _head = _list[i];
        _list[i] = -1;
        lastcn++;
        i++;
      }
    }
  }
  lastn = lastcn*_pperc;//*/


#ifdef PGAS
  {
    double tx, ty, tz;
    int i = 0;
    int j;
    lastn = lastcn = 0;
    const double dxyz = 1.0;
    for (tx = dxyz; tx < LX-dxyz; tx += dxyz) {
      for (ty = dxyz; ty < LY-dxyz; ty += dxyz) {
        for (tz = dxyz; tz < LZ-dxyz; tz += dxyz) {
          if (NORM(tx-LX/2.,ty-LY/2.,tz-LZ/2.) < 8.)
            continue;
          j = i*_pperc;
          _rx[j] = tx;
          _ry[j] = ty;
          _rz[j] = tz;
          _id[j] = next_free_id++;
          _spec[j] = gas_spec;
          _spec[j+1] = -1;

          _vx[j] = get_ran_gaussian(0, 0., sqrt(gas_kt));
          _vy[j] = get_ran_gaussian(0, 0., sqrt(gas_kt));
          _vz[j] = get_ran_gaussian(0, 0., sqrt(gas_kt));

          _head = _list[i];
          _list[i] = -1;
          lastcn++;
          i++;
        }
      }
    }

    j = i*_pperc;
    _rx[j] = LX/2.;
    _ry[j] = LY/2.;
    _rz[j] = LZ/2.;
    double initdist = 0.4;
    tmpx = (drand[0]()*initdist*2.-initdist)/sqrt(3);
    tmpy = (drand[0]()*initdist*2.-initdist)/sqrt(3);
    tmpz = sqrt(initdist*initdist-tmpx*tmpx-tmpy*tmpy);
    _rx[j+1] = LX/2. + tmpx;
    _ry[j+1] = LY/2. + tmpy;
    _rz[j+1] = LZ/2. + tmpz;

    _vx[j] = _vy[j] = _vz[j] = 0.;
    _vx[j+1] = _vy[j+1] = _vz[j+1] = 0.;

    _id[j] = next_free_id++;
    _id[j+1] = next_free_id++;

    _spec[j] = _spec[j+1] = 0;

    _head    = _list[i];
    _list[i] = -1;
    lastcn++;
  }


  /*i=0;
  j = i*_pperc;
  _rx[j] = 3.;
  _ry[j] = 3.;
  _rz[j] = 3.;
  _id[j] = next_free_id++;
  next_free_id++;
  _spec[j] = gas_spec;
  _spec[j+1] = -1;

  _vx[j] = 0.;
  _vy[j] = 0.;
  _vz[j] = 1.;

  _head = _list[i];
  _list[i] = -1;
  lastcn++;

   i = 1;
   j = i*_pperc;
  _rx[j] = 3.;
  _ry[j] = 3.;
  _rz[j] = 4.1;
  _rx[j+1] = 3.;
  _ry[j+1] = 3.;
  _rz[j+1] = 4.4;

  _vx[j] = _vy[j] = _vz[j] = 0.;
  _vx[j+1] = _vy[j+1] = _vz[j+1] = 0.;

  _id[j] = next_free_id++;
  _id[j+1] = next_free_id++;

  _spec[j] = _spec[j+1] = 0;

  _head = _list[i];
  _list[i] = -1;
  lastcn++;*/

  lastn = lastcn*_pperc;
#endif


  return 0;
}

