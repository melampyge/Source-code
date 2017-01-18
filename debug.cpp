/*		debug.cpp		*/
//////////////////////////////////////////
// Debugging functions etc		//
//////////////////////////////////////////


#include "debug.h"


using namespace std;


////////////////////////////////////
// void __assert_(...)
// ---------------
// check for i in [l,u)
// ===
// Parameters:
//	long i : index variable to check
//	long l : lower boundary
//	long u : upper boundary
//	char *s: descriptive string
void __assert_(long i, long l, long u, const char *s)
{
  if (i < l || i >= u) {
    cerr << "__DEBUG at t=" << curtime << ": '" << s << "' is out of bound: l=" << l << ": u=" << u << ": i=" << i << endl;
    process_signal_handler(0xFFFFFFFF);
  }
}
void __assert_(int i, int l, int u, const char *s)
{
  if (i < l || i >= u) {
    cerr << "__DEBUG at t=" << curtime << ": '" << s << "' is out of bound: l=" << l << ": u=" << u << ": i=" << i << endl;
    process_signal_handler(0xFFFFFFFF);
  }
}


void __assert_(double i, double l, double u, const char *s)
{
  if (i < l || i >= u) {
    cerr << "__DEBUG at t=" << curtime << ": '" << s << "' is out of bound: l=" << l << ": u=" << u << ": i=" << i << endl;
    process_signal_handler(0xFFFFFFFF);
  }
}


////////////////////////////////////
// int check_rvffd(...)
// ---------------
// check r, v, f, fd for nan and inf
// for non periodic boundaries, also check for outofBound
// ===
// Parameters:
//
// Returns:
//	0: on sucess
//	1: on failure
int check_rvffd()
{
  long __i, __j;
  for (__j = 0; __j < lastcn; ++__j) {
    if (_list[__j] != -1)
      continue;
    for (__i = __j*_pperc; __i < (__j+1)*_pperc; ++__i) {
      if (isnan(_rx[__i]) || !isfinite(_rx[__i])) {
        cout << "__DEBUG: at t=" << curtime << " _rx[" << __i << "] is " << _rx[__i] << endl;
        return 1;
      }
      if (isnan(_ry[__i]) || !isfinite(_ry[__i])) {
        cout << "__DEBUG: at t=" << curtime << " _ry[" << __i << "] is " << _ry[__i] << endl;
        return 1;
      }
      if (isnan(_rz[__i]) || !isfinite(_rz[__i])) {
        cout << "__DEBUG: at t=" << curtime << " _rz[" << __i << "] is " << _rz[__i] << endl;
        return 1;
      }

      if (isnan(_vx[__i]) || !isfinite(_vx[__i])) {
        cout << "__DEBUG: at t=" << curtime << " _vx[" << __i << "] is " << _vx[__i] << endl;
        return 1;
      }
      if (isnan(_vy[__i]) || !isfinite(_vy[__i])) {
        cout << "__DEBUG: at t=" << curtime << " _vy[" << __i << "] is " << _vy[__i] << endl;
        return 1;
      }
      if (isnan(_vz[__i]) || !isfinite(_vz[__i])) {
        cout << "__DEBUG: at t=" << curtime << " _vz[" << __i << "] is " << _vz[__i] << endl;
        return 1;
      }

      if (isnan(_fx[__i]) || !isfinite(_fx[__i])) {
        cout << "__DEBUG: at t=" << curtime << " _fx[" << __i << "] is " << _fx[__i] << endl;
        return 1;
      }
      if (isnan(_fy[__i]) || !isfinite(_fy[__i])) {
        cout << "__DEBUG: at t=" << curtime << " _fy[" << __i << "] is " << _fy[__i] << endl;
        return 1;
      }
      if (isnan(_fz[__i]) || !isfinite(_fz[__i])) {
        cout << "__DEBUG: at t=" << curtime << " _fz[" << __i << "] is " << _fz[__i] << endl;
        return 1;
      }

      if (isnan(_fdx[__i]) || !isfinite(_fdx[__i])) {
        cout << "__DEBUG: at t=" << curtime << " _fdx[" << __i << "] is " << _fdx[__i] << endl;
        return 1;
      }
      if (isnan(_fdy[__i]) || !isfinite(_fdy[__i])) {
        cout << "__DEBUG: at t=" << curtime << " _fdy[" << __i << "] is " << _fdy[__i] << endl;
        return 1;
      }
      if (isnan(_fdz[__i]) || !isfinite(_fdz[__i])) {
        cout << "__DEBUG: at t=" << curtime << " _fdz[" << __i << "] is " << _fdz[__i] << endl;
        return 1;
      }
#if defined (BBC) || defined (RBC) || defined (HBC)
      if (outofBound(__i)) {
        cout << "__DEBUG: at t=" << curtime << " particle " << __i << " is out of bounds!" << endl;
        cout << "_rx=" << _rx[__i] << ": _ry=" << _ry[__i] << ": _rz=" << _rz[__i] << endl;
        return 1;
      }
#endif
    }
  }
  return 0;
}



////////////////////////////////////
// int check_list_lastcn(...)
// ---------------
// check whether there are no cells beyond lastcn
// ===
// Parameters:
//
// Returns:
//	0: on sucess
//	1: on failure
int check_list_lastcn()
{
  long i;
  for (i = _ceil(INITSIZE/_pperc)-2; i >= lastcn; --i)
    if (_list[i] == -1)
      return 1;
  return 0;
}


