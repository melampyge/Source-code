/*		bc_func.cpp		*/
//////////////////////////////////////////
// functions related to boundaries are  //
// defined here				//
//////////////////////////////////////////


#include "bc_func.h"



#if defined (BBC) || defined (RBC) || defined (HBC)


/////////////////////////////////////
// outofBound(...)
// -----------------
// check whether particle is within boundaries or not (BBC/RBC/HBC)
// ===
// Parameters:
//   long i: index of cell particle
// ===
// Returns:
//   TRUE  (!=0): if particle is out of boundaries
//   FALSE (==0): if particle is wihtin boundaries
int outofBound(long i)
{
#if defined (BBC) || defined (RBC)
  if (_rx[i] < 0. || _rx[i] > LX || _ry[i] < 0. || _ry[i] > LY || _rz[i] < 0. || _rz[i] > LZ)
    return 1;   // TRUE
  else
    return 0;   // FALSE
#endif
#ifdef HBC
  if (_rz[i] < 0. || _rz[i] > LZ)
    return 1;   // TRUE
  else
    return 0;   // FALSE
#endif
}



//////////////////////////////////
// check_boundcollision(...)
// -------------------------
// check if particle interacts with wall and calculate interaction (BBC/RBC)
// ===
// Parameters:
//   long i	: index of cell particle
//   int tid    : thread id
int check_boundcollision(long i, int tid)
{
  double st, st2;       // smallest t
  double dts = dt;      // available time during this step
  int wallnum = 0;      // number of wall at which first collision takes place
                        // 1: x=0 (perpendicular to x axis)
                        // 2: x=L
                        // 3: y=0
                        // 4: y=L
                        // 5: z=0
                        // 6: z=L
#ifdef __DEBUG
  int lc = 0;           // loop counter to detect/avoid endless loops
#endif

  ////////////////////////////
  // calculate boundary collision and make changes to r and v according to boundary conditions

  while (outofBound(i)) {       // while outofBound(i) loop
    ////////////////////////////
    // calculate boundary collision times and find smallest (stored in st)
    // -----
    // IMPORTANT: during the backstep calculation, it could happen that the particle stays
    // outside of the box due to floating point rounding errors.
    st = 0.;
#if defined (BBC) || defined (RBC)
    st2 = _rx[i]/_vx[i];
    ///////////////////////////////////////
    // st2 is time we need to get to x=0 with current velocity
    // check - whether particle is really outside (_rx[i]<0.,etc)
    //       - whether st2 > st
    if (_rx[i] < 0. && st2 > st) {          // x=0 wall
      st = st2;
      wallnum = 1;
    }
    st2 = (_rx[i]-LX)/_vx[i];
    if (_rx[i] > LX && st2 > st) {        // x=L wall
      st = st2;
      wallnum = 2;
    }
    st2 = _ry[i]/_vy[i];
    if (_ry[i] < 0. && st2 > st) {        // y=0 wall
      st = st2;
      wallnum = 3;
    }
    st2 = (_ry[i]-LY)/_vy[i];
    if (_ry[i] > LY && st2 > st) {        // y=L wall
      st = st2;
      wallnum = 4;
    }
    st2 = _rz[i]/_vz[i];
    if (_rz[i] < 0. && st2 > st) {        // z=0 wall
      st = st2;
      wallnum = 5;
    }
    st2 = (_rz[i]-LZ)/_vz[i];
    if (_rz[i] > LZ && st2 > st) {        // z=L wall
      st = st2;
      wallnum = 6;
    }
#elif defined (HBC)
    st2 = _rz[i]/_vz[i];
    if (_rz[i] < 0. && st2 > st) {        // z=0 wall
      st = st2;
      wallnum = 5;
    }
    st2 = (_rz[i]-LZ)/_vz[i];
    if (_rz[i] > LZ && st2 > st) {        // z=L wall
      st = st2;
      wallnum = 6;
    }
#endif

    /////////////////////////////
    // move back to collision position (EPSILON is for compensation of fp inaccuracy)
    // Could lead to particle placed outside of box (but VERY unlikely, never happened so far)
    // perhaps add another outofBound(i) test after position update
    _rx[i] -= (st+EPSILON)*_vx[i];
    _ry[i] -= (st+EPSILON)*_vy[i];
    _rz[i] -= (st+EPSILON)*_vz[i];

    dts = st;

#ifdef PMEASUREMENT
    /////////////////////////////
    // calc momentum transfer to wall
    update_pressure(i, wallnum, _spec[i], tid);
#endif

    /////////////////////////////
    // set new velocity according to BC and subtract st from dts
    bc_vchange(i, wallnum);

    /////////////////////////////
    // now move with altered v to new position v*st and check again for collision
    _rx[i] += _vx[i]*st;
    _ry[i] += _vy[i]*st;
    _rz[i] += _vz[i]*st;

#ifdef __DEBUG
    lc++;
    if (lc > 1000) {
#pragma omp critical
{
      std::cerr << "__DEBUG: Possible endless loop () at t=" << curtime << " and i=" << i << "\n"
                << "_rx[i]=" << _rx[i] << ": _ry[i]=" << _ry[i] << ": _rz[i]=" << _rz[i] << "\n"
                <<": _fx="<<_fx[i]<<": _fy="<<_fy[i]<<": _fz="<<_fz[i]
                <<": _fdx="<<_fdx[i]<<": _fdy="<<_fdy[i]<<": _fdz="<<_fdz[i] << std::endl;
      process_signal_handler(0xFFFFFFFF);
      resolve_int();
      _Exit(1);
} // end omp critical
    }
#endif
  }     // end while outofBound(i) loop

  return 0;
}



////////////////////////////////////
// bc_vchange (...)
// ----------------------
// calculate velocities due to boundary collision (BBC)
// ===
// Parameters:
//   long i	: index of cell particle
//   int wallnum: index of wall where first collision takes place
//                -> see check_boundcollision() for numbering
void bc_vchange(long i, int wallnum)
{
#ifdef BBC
  _vx[i] *= -1.;
  _vy[i] *= -1.;
  _vz[i] *= -1.;
#elif defined (RBC)
  switch (wallnum) {
    case 1:
    case 2:
      _vx[i] *= -1.;
      break;
    case 3:
    case 4:
      _vy[i] *= -1.;
      break;
    case 5:
    case 6:
      _vz[i] *= -1.;
      break;
  }
#elif defined (HBC)
  _vx[i] *= -1.;
  _vy[i] *= -1.;
  _vz[i] *= -1.;
#endif
}



#endif


