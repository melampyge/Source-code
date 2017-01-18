/*		verlet.cpp		*/
//////////////////////////////////////////
// function definitions for velocity	//
// verlet algorithm			//
//////////////////////////////////////////


#include "verlet.h"


////////////////////////////////////
// void init_verlet(...)
// ---------------
// calculate forces from inital particle positions and velocities
// ===
// Parameters:
void init_verlet()
{
  long i, j, k;


  for (k = 0; k < lastcn; ++k) {  // loop through particles 1
    if (_list[k] == -1) {
      /////////////////////////////////////
      // loop over all particles of one cell
      for (i = _pperc*k; i < _pperc*(k+1); ++i) {	// loop over k
#ifdef PGAS
        if (_spec[i] == -1)
          continue;
#endif
        for (j = vlist_offset[i]; j < vlist_nc[i]; ++j)
          calcforces(i, vlist[j], j, 0);
        //calcbackgroundforces(i);
        calcbackgroundfd(i);		// as long as there are no other backgroundforces
      }  // end loop over k
      calcintracellforces(_pperc*k, 0);
    } // end if
  }  // end loop through particles 1


  for (k = 0; k < lastcn; ++k) {  // loop through particles 3
    if (_list[k] == -1) {
      /////////////////////////////////////
      // loop over all particles of one cell
      for (i = _pperc*k; i < _pperc*(k+1); ++i) {	// loop over k
#ifdef PGAS
        if (_spec[i] == -1)
          continue;
#endif
        for (j = vlist_offset[i]; j < vlist_nc[i]; ++j)
          calcrandf(i, vlist[j], j);
      }  // end loop over k
    } // end if
  }  // end loop through particles 3


#ifdef LOCALSTRESS
  for (k = 0; k < lastcn; ++k) {
    if (_list[k] != -1)
      continue;
    for (i = _pperc*k; i < _pperc*(k+1); ++i) {
      stress_rx[i] = _rx[i];
      stress_ry[i] = _ry[i];
      stress_rz[i] = _rz[i];
    }
  }
#endif

}


////////////////////////////////////
// void verlet_step_dpd(...)
// ---------------
// calculate one verlet step using dpd:
//   - calculate intermediate velocities
//   - calculate new positions
//   - calculate forces from new positions/velocities
//   - calculate actual new velocities
//   - calculate dissipative forces once again
// ===
// Parameters:
void verlet_step_dpd()
{
  long i, j, k, l;
  long size;
  int tid;

  /////////////////////////////////////
  // check if verlet list is still valid and update accordingly
  START_TIMER(1)
  check_verlet_list();
  STOP_TIMER(1)

#pragma omp parallel shared(rna,_rx,_ry,_rz,_vx,_vy,_vz,_fx,_fy,_fz,_fdx,_fdy,_fdz) private(i,j,k,l,tid,size)
{
#pragma omp flush(_int_flag,_int_no)
  if (!_int_flag) {
  tid = omp_get_thread_num();

#pragma omp barrier		// there is an implied barrier at end of single block
#pragma omp single
{ START_TIMER(2) }
#pragma omp for //schedule(auto)
  for (k = 0; k < lastcn; ++k) {  // loop through particles 1
    if (_list[k] != -1)
      continue;
    for (i = _pperc*k; i < _pperc*(k+1); ++i) {
#ifdef PGAS
      if (_spec[i] == -1)
        continue;
#endif
      /////////////////////////////////////
      // calculate v'(t+1)=v(t)+1/2*f(t)
      calcv(i);
      /////////////////////////////////////
      // calculate r(t+1)=r(t)+v'(t+1)*dt
      calcr(i);
#if defined (BBC) || defined (RBC) || defined (HBC)
      ////////////////////////////////////
      // for BBC or RBC, we need to check for boundary collisions
      check_boundcollision(i, tid);
#endif
      ////////////////////////////////////
      // set forces to 0, to calculate new ones in next loop
      _fx[i] = _fy[i] = _fz[i] = _fdx[i] = _fdy[i] = _fdz[i] = 0.;
    } // end loop over k
  }  // end loop through particles 1

// implied barrier from end of omp for
#pragma omp single
{ STOP_TIMER(2) }


  /////////////////////////////////////
  // check for cell death and division
// implied barrier from end of omp single
#pragma omp single
{
  START_TIMER(3)
  cell_death();
  STOP_TIMER(3)

  START_TIMER(4)
  cell_division();
  STOP_TIMER(4)
}

#ifdef __DEBUG
/*#pragma omp master
{
    if (!check_rvffd()) {
      process_signal_handler(0xFFFFFFFF);
      //return ;
    }
} // end omp master*/
#endif

  /////////////////////////////////////
  // calculate new forces and store in f{c, r,d}
// implied barrier from end of omp single
#pragma omp single
{ START_TIMER(5) }

  /////////////////////////////////////
  // set rna to NaN
  size = (sizeof(unsigned char))*(vlist_nc[lastn-1]/omp_get_num_threads()+(((vlist_nc[lastn-1]%omp_get_num_threads()) == 0)?0:1));
  memset(rnai+size*tid, 0x00, size);

#pragma omp barrier		// array filling has to be finished first
#pragma omp single
{ STOP_TIMER(5)
  START_TIMER(6) }

#pragma omp for //schedule(auto)
  for (k = 0; k < lastcn; ++k) {  // loop through particles 2
    if (_list[k] == -1) {
      /////////////////////////////////////
      // loop over all particles of one cell
      for (i = _pperc*k; i < _pperc*(k+1); ++i) {	// loop over k
#ifdef PGAS
        if (_spec[i] == -1)
          continue;
#endif
        for (j = vlist_offset[i]; j < vlist_nc[i]; ++j)
          calcforces(i, vlist[j], j, tid);
        //calcbackgroundforces(i);
        calcbackgroundfd(i);		// as long as there are no other backgroundforces
      }  // end loop over k
      calcintracellforces(_pperc*k, tid);
    } // end if
  }  // end loop through particles 2

// implied barrier from end of omp for
#pragma omp single
{ STOP_TIMER(6)
  START_TIMER(8) }

#ifdef FIXEDCELLS		// has to be calculated before calcv()!!!!!!!!
#pragma omp for //schedule(auto)
  for (i = 0; i < fixed_wlist_max; ++i) {
    double fxtmp = fixed_kspringx*(fixed_wlist_rx0[i]-_rx[fixed_wlist[i]]);
    double fytmp = fixed_kspringy*(fixed_wlist_ry0[i]-_ry[fixed_wlist[i]]);
    double fztmp = fixed_kspringz*(fixed_wlist_rz0[i]-_rz[fixed_wlist[i]]);

    if (_rz[fixed_wlist[i]] < fixed_deltaz+fixed_init_displacement+1.)
      fixed_pkz0[tid] -= fztmp;
    else
      fixed_pkzL[tid] -= fztmp;
    _fx[fixed_wlist[i]] += fxtmp;
    _fy[fixed_wlist[i]] += fytmp;
    _fz[fixed_wlist[i]] += fztmp;
  }
#endif // ifdef FIXEDCELLS

  /////////////////////////////////////
  // now calculate inter cell random forces
#pragma omp for //schedule(auto)
  for (k = 0; k < lastcn; ++k) {  // loop through particles 2
    if (_list[k] == -1) {
      /////////////////////////////////////
      // loop over all particles of one cell
      for (i = _pperc*k; i < _pperc*(k+1); ++i) {	// loop over k
#ifdef PGAS
        if (_spec[i] == -1)
          continue;
#endif
        for (j = vlist_offset[i]; j < vlist_nc[i]; ++j)
          calcrandf(i, vlist[j], j);
        /////////////////////////////////////
        // at this point new velocities for i can be calculated (calcrandf only needs positions)
        calcv(i);
        ///////////////////////////////////
        // reset dissipative forces fd to zero
        _fdx[i] = _fdy[i] = _fdz[i] = 0.;
      }  // end loop over k
    } // end if
  }  // end loop through particles 2


#ifdef LOCALSTRESS
#pragma omp for
  for (k = 0; k < lastcn; ++k) {
    if (_list[k] != -1)
      continue;
    for (i = _pperc*k; i < _pperc*(k+1); ++i) {
      localstress_momentum_update(i, tid);
    }
  }
#endif	// ifdef LOCALSTRESS


// implied barrier from end of omp for
#pragma omp single
{ STOP_TIMER(8)
  START_TIMER(10) }

#ifdef LOCALSTRESS
  //////////////////////////////////////
  // set stress_fnextd[tid] to zero (will be newly calculated below)
  std::fill(stress_fnextd[tid], stress_fnextd[tid]+stress_length, 0.);
#endif

#ifdef CONSTP_PBC
  constp_pbc_fxnextd[tid] = 0.;
  constp_pbc_fynextd[tid] = 0.;
  constp_pbc_fznextd[tid] = 0.;
#endif

  //////////////////////////////////////
  // once again calculate dissipative forces with new velocities
#pragma omp for //schedule(auto)
  for (k = 0; k < lastcn; ++k) {  // loop through particles 4
    if (_list[k] == -1) {
      /////////////////////////////////////
      // loop over all particles of one cell
      for (i = _pperc*k; i < _pperc*(k+1); ++i) {	// loop over k
#ifdef PGAS
        if (_spec[i] == -1)
          continue;
#endif
        for (j = vlist_offset[i]; j < vlist_nc[i]; ++j) {
          calcfd(i, vlist[j], tid);
        }
        calcbackgroundfd(i);
      }  // end loop over k
      calcintracellfd(_pperc*k, tid);
    }
  } // end loop through particles 4

#pragma omp single
{ STOP_TIMER(10) }
  } // end if (!_int_flag)
} // end parallel
}



////////////////////////////////////
// void calcv(...)
// ---------------
// calculate intermediate velocity v' or new velocity v after force update
// ===
// Parameters:
//   long i: index of particle
inline void calcv(long i)
{
  int speci = _spec[i];
  _vx[i] += 0.5*(_fx[i] + _fdx[i])/m[speci]*dt;
  _vy[i] += 0.5*(_fy[i] + _fdy[i])/m[speci]*dt;
  _vz[i] += 0.5*(_fz[i] + _fdz[i])/m[speci]*dt;
}


////////////////////////////////////
// void calcr(...)
// ---------------
// calculate new positions r
// ===
// Parameters:
//   long i: index of particle
inline void calcr(long i)
{
  _rx[i] += _vx[i]*dt;
  vlist_dx[i] += _vx[i]*dt;
  _ry[i] += _vy[i]*dt;
  vlist_dy[i] += _vy[i]*dt;
  _rz[i] += _vz[i]*dt;
  vlist_dz[i] += _vz[i]*dt;
}


