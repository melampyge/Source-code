/*		pressure.cpp		*/
//////////////////////////////////////////
// calculate pressure			//
//////////////////////////////////////////


#include "pressure.h"

#ifdef PMEASUREMENT



////////////////////////////////////
// void update_pressure(...)
// ---------------
// calculate pressure due to boundary collisions
// ===
// Parameter:
//   long i		: index of colliding particle
//   int wallnum	: wall index (important for RBC)
//   int speci		: species
//   int tid            : thread id
void update_pressure(long i, int wallnum, int speci, int tid)
{
  switch (wallnum) {
  case 1:
    pm_px0[tid] += 2*m[speci]*_vx[i];
    break;
  case 2:
    pm_pxL[tid] += 2*m[speci]*_vx[i];
    break;
  case 3:
    pm_py0[tid] += 2*m[speci]*_vy[i];
    break;
  case 4:
    pm_pyL[tid] += 2*m[speci]*_vy[i];
    break;
  case 5:
    pm_pz0[tid] += 2*m[speci]*_vz[i];
    break;
  case 6:
    pm_pzL[tid] += 2*m[speci]*_vz[i];
    break;
  default:
    break;
  }
}



#endif

