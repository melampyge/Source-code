/*		model.cpp		*/
//////////////////////////////////////////
// model specific function definitions	//
//////////////////////////////////////////


#include "model.h"




////////////////////////////////////
// void calcforces(...)
// ---------------
// calculate conservative, random and dissipative forces
// and store them in double *f and double *fd
// IMPORTANT: should only calculate forces between cell
// particles of different cells (verletlist should be constructed to ensure that)
// ===
// Parameters:
//   long i: index of cell particle 1
//   long j: index of cell particle 2
//   long vlist_j: index of j in vlist
//   int tid     : thread index
void calcforces(long i, long j, long vlist_j, int tid)
{
  double dx, dy, dz, dr, drsq, ftmp;
  int speci = _spec[i];
  int specj = _spec[j];
  double rpp_ = rpp[speci*MAXSPECNUM+specj];


  //////////////////////////////////
  // calculate distance between i and j
  dx = distx(i,j);
  dy = disty(i,j);
  dz = distz(i,j);
  drsq = dx*dx + dy*dy + dz*dz + EPSILON;

  /////////////////////////////////
  // caculate forces except dissipative ones for dr < RPP
  if (drsq < rpp_*rpp_) {		// if dr < RPP
    dr = sqrt(drsq);
    ////////////////////////////////////
    // volume exclusion force
    // F^CC=f0*{(RPP/dr)^5-1.}*e_rij
    ftmp = f0[speci*MAXSPECNUM+specj]*((rpp_*rpp_*rpp_*rpp_*rpp_)/(drsq*drsq*dr)-1.)/dr;


    ////////////////////////////////////
    // cell adhesion force
    // F^a=-f1*e_rij
    ftmp += -f1[speci*MAXSPECNUM+specj]/dr;
    _fx[i] -= dx*ftmp;
    _fy[i] -= dy*ftmp;
    _fz[i] -= dz*ftmp;

#ifdef PGAS
 #ifdef __DEBUG
    if (speci == gas_spec || specj == gas_spec) {
      if (fabs(f1[speci*MAXSPECNUM+specj]) > 0.) {
        std::cout << "__DEBUG at t=" << curtime << ": f1!=0 for gas_spec: i=" << i << ": j=" << j << std::endl;
      }
    }
 #endif
#endif


    //gammat_=MAX(gammat[speci],gammat[specj]);
    ////////////////////////////////////
    // random force
    // ------------
    // -> draw gaussian distributed random number
    // -> fr{xyz} = SIGMA*o^R(r_{ij})*xi/dr*d{xyz}*1/sqrt(dt)
    //    ONLY DRAW RANDOM NUMBER HERE
    // -> forces are added in calcrandf()
    if (i < j) {
      double fxi = randf_preft[speci*MAXSPECNUM+specj]*(1.-dr/rt[speci*MAXSPECNUM+specj])*get_ran_gaussian(tid)/dr;
      rna[vlist_j] = fxi;
      rnai[vlist_j] = 1;
      /////////////////////////////////////////
      // Now find vlist_i of j
      for (int k = vlist_offset[j]; k < vlist_nc[j]; ++k)
        if (vlist[k] == i) {
          rna[k] = fxi;
          rnai[k] = 1;
          break;
        }
      ftmp += fxi;
#ifdef LOCALSTRESS
      if (stress_meas_enabled)
        localstress_update(i, j, ftmp, tid);
#endif
#ifdef CONSTP_PBC
      constp_pbc_fxnext[tid] += dx*dx*(ftmp);
      constp_pbc_fynext[tid] += dy*dy*(ftmp);
      constp_pbc_fznext[tid] += dz*dz*(ftmp);
#endif
    }
    
    ////////////////////////////////////
    // calc dissipative forces
    calcfd(i,j,dx,dy,dz,drsq,dr,tid);
  }		// end if dr < RPP
  else
    calcfd(i,j,dx,dy,dz,drsq,sqrt(drsq),tid);
}



////////////////////////////////////
// void calcintracellforces(...)
// ---------------
// calculate conservative and random forces of intra cell particles
// calls calcfd as well
// ===
// Parameters:
//   long i  : index  of cell particle 1
//   long tid: thread index (for random number generation)
void calcintracellforces(long i, int tid)
{
  double dx, dy, dz, dr, drsq, ftmp, xi;
  int speci = _spec[i];
  long j = i+1;
  double r0_ = r0[speci*MAXSPECNUM+speci];
  double rppsq_ = rpp[speci*MAXSPECNUM+speci]*rpp[speci*MAXSPECNUM+speci];


#ifdef PGAS
  if (speci == gas_spec)
    return;
#endif

  //////////////////////////////////
  // calculate distance between i and j
  dx = distx(i,j);
  dy = disty(i,j);
  dz = distz(i,j);
  drsq = dx*dx + dy*dy + dz*dz + EPSILON;

  /////////////////////////////////
  // caculate forces except dissipative ones for dr < RPP
  if (drsq < rppsq_) {		// if dr < RPP
    dr = sqrt(drsq);
    ////////////////////////////////////
    // growth force (only for particles of same cell)
    // the if statement distinguishs between same cell (if) and other cell (else)
    ftmp = b[speci*MAXSPECNUM+speci]/((dr+r0_)*(dr+r0_)*(dr+EPSILON));


    ////////////////////////////////////
    // random force
    // ------------
    // -> draw gaussian distributed random number
    // -> fr{xyz} = SIGMA*o^R(r_{ij})*xi/dr*d{xyz}*1/sqrt(dt)
    // -> no weight function for intra cell random force
    xi = get_ran_gaussian(tid);
    ftmp += randf_prefc[speci*MAXSPECNUM+speci]*xi/dr;
    _fx[i] -= dx*ftmp;
    _fy[i] -= dy*ftmp;
    _fz[i] -= dz*ftmp;
    _fx[j] += dx*ftmp;
    _fy[j] += dy*ftmp;
    _fz[j] += dz*ftmp;
#ifdef LOCALSTRESS
    if (stress_meas_enabled)
      localstress_update(i, j, ftmp, tid);
#endif
#ifdef CONSTP_PBC
    constp_pbc_fxnext[tid] += dx*dx*ftmp;
    constp_pbc_fynext[tid] += dy*dy*ftmp;
    constp_pbc_fznext[tid] += dz*dz*ftmp;
#endif

    ////////////////////////////////////
    // calc dissipative forces
    calcintracellfd(i,dx,dy,dz,drsq,dr,tid);
  }
  else
    calcintracellfd(i,dx,dy,dz,drsq,sqrt(drsq),tid);
}



////////////////////////////////////
// void calcfd(...)
// ---------------
// calculate dissipative forces and store them in double *fd
// ===
// Parameters:
//   long i: index of cell particle 1
//   long j: index of cell particle 2
void calcfd(long i, long j, int tid)
{
  double dx, dy, dz, dr, drsq, ftmp;
  int speci = _spec[i];
  int specj = _spec[j];
  double rt_ = rt[speci*MAXSPECNUM+specj];

  //////////////////////////////////
  // calculate distance between i and j
  dx = distx(i,j);
  dy = disty(i,j);
  dz = distz(i,j);
  drsq = dx*dx + dy*dy + dz*dz + EPSILON;

  /////////////////////////////////
  // caculate dissipative forces for dr < RT
  if (drsq < rt_*rt_) {		// if dr < RT
    ////////////////////////////////////
    // dissipative force between cells
    // -----------------
    // -> fd{xyz} = -GAMMA*(o^R(r_{ij}))^2*(dvx*dx+dvy*dy+dvz*dz)/dr^2*d{xyz}
    dr = sqrt(drsq);
    ftmp = -gammat[speci*MAXSPECNUM+specj]*(1.-dr/rt_)*(1.-dr/rt_)*((_vx[j]-_vx[i])*dx+(_vy[j]-_vy[i])*dy
        + (_vz[j]-_vz[i])*dz)/drsq;
    _fdx[i] -= ftmp*dx;
    _fdy[i] -= ftmp*dy;
    _fdz[i] -= ftmp*dz;

#ifdef PGAS
 #ifdef __DEBUG
    if (speci == gas_spec || specj == gas_spec) {
      if (fabs(ftmp) > 0.) {
        std::cout << "__DEBUG at t=" << curtime << ": calcfd ftmp>0 for gas_spec: i=" << i << ": j=" << j << std::endl;
      }
    }
 #endif
#endif

#ifdef LOCALSTRESS
    if (i < j && stress_meas_enabled)
      localstressd_update(i, j, ftmp, tid);
#endif
#ifdef CONSTP_PBC
    if (i < j) {
      constp_pbc_fxnextd[tid] += dx*dx*ftmp;
      constp_pbc_fynextd[tid] += dy*dy*ftmp;
      constp_pbc_fznextd[tid] += dz*dz*ftmp;
    }
#endif
  }		// end if dr < RT 
}


void calcfd(long i, long j, double dx, double dy, double dz, double drsq, double dr, int tid)
{
  double ftmp;
  int speci = _spec[i];
  int specj = _spec[j];
  double rt_ = rt[speci*MAXSPECNUM+specj];

  /////////////////////////////////
  // caculate dissipative forces for dr < RT
  if (drsq < rt_*rt_) {		// if dr < RT
    ////////////////////////////////////
    // dissipative force between cells
    // -----------------
    // -> fd{xyz} = -GAMMA*(o^R(r_{ij}))^2*(dvx*dx+dvy*dy+dvz*dz)/dr^2*d{xyz}
    ftmp = -gammat[speci*MAXSPECNUM+specj]*(1.-dr/rt_)*(1.-dr/rt_)*((_vx[j]-_vx[i])*dx+(_vy[j]-_vy[i])*dy
        + (_vz[j]-_vz[i])*dz)/drsq;
    _fdx[i] -= ftmp*dx;
    _fdy[i] -= ftmp*dy;
    _fdz[i] -= ftmp*dz;

#ifdef PGAS
 #ifdef __DEBUG
    if (speci == gas_spec || specj == gas_spec) {
      if (fabs(ftmp) > 0.) {
        std::cout << "__DEBUG at t=" << curtime << ": calcfd2 ftmp>0 for gas_spec: i=" << i << ": j=" << j << std::endl;
      }
    }
 #endif
#endif

#ifdef LOCALSTRESS
    if (i < j && stress_meas_enabled)
      localstressd_update(i, j, ftmp, tid);
#endif
#ifdef CONSTP_PBC
    if (i < j) {
      constp_pbc_fxnextd[tid] += dx*dx*ftmp;
      constp_pbc_fynextd[tid] += dy*dy*ftmp;
      constp_pbc_fznextd[tid] += dz*dz*ftmp;
    }
#endif
  }		// end if dr < RT 
}



////////////////////////////////////
// void calcintracellfd(...)
// ---------------
// calculate dissipative forces for intra cell particles
// ===
// Parameters:
//   long i: index of cell particle
void calcintracellfd(long i, int tid)
{
  double dx, dy, dz, drsq, ftmp;
  int speci = _spec[i];
  long j = i+1;


#ifdef PGAS
  if (speci == gas_spec)
    return;
#endif

  //////////////////////////////////
  // calculate distance between i and j
  dx = distx(i,j);
  dy = disty(i,j);
  dz = distz(i,j);
  drsq = dx*dx + dy*dy + dz*dz + EPSILON;

  /////////////////////////////////
  // caculate dissipative forces for dr < RT
  if (drsq < rt[speci*MAXSPECNUM+speci]*rt[speci*MAXSPECNUM+speci]) {		// if dr < RT
    ////////////////////////////////////
    // -> no weight function for intra cell dissipation
    ftmp = -gammac[speci*MAXSPECNUM+speci]*((_vx[j]-_vx[i])*dx+(_vy[j]-_vy[i])*dy
      + (_vz[j]-_vz[i])*dz)/drsq;
    _fdx[i] -= dx*ftmp;
    _fdy[i] -= dy*ftmp;
    _fdz[i] -= dz*ftmp;
    _fdx[j] += dx*ftmp;
    _fdy[j] += dy*ftmp;
    _fdz[j] += dz*ftmp;

#ifdef LOCALSTRESS
    if (stress_meas_enabled)
      localstressd_update(i, j, ftmp, tid);
#endif
#ifdef CONSTP_PBC
    constp_pbc_fxnextd[tid] += dx*dx*ftmp;
    constp_pbc_fynextd[tid] += dy*dy*ftmp;
    constp_pbc_fznextd[tid] += dz*dz*ftmp;
#endif
  }		// end if dr < RT 
}


void calcintracellfd(long i, double dx, double dy, double dz, double drsq, double dr, int tid)
{
  double ftmp;
  int speci = _spec[i];
  long j = i+1;


#ifdef PGAS
  if (speci == gas_spec)
    return;
#endif

  /////////////////////////////////
  // caculate dissipative forces for dr < RT
  if (drsq < rt[speci*MAXSPECNUM+speci]*rt[speci*MAXSPECNUM+speci]) {		// if dr < RT
    ////////////////////////////////////
    // -> no weight function for intra cell dissipation
    ftmp = -gammac[speci*MAXSPECNUM+speci]*((_vx[j]-_vx[i])*dx+(_vy[j]-_vy[i])*dy
      + (_vz[j]-_vz[i])*dz)/drsq;
    _fdx[i] -= dx*ftmp;
    _fdy[i] -= dy*ftmp;
    _fdz[i] -= dz*ftmp;
    _fdx[j] += dx*ftmp;
    _fdy[j] += dy*ftmp;
    _fdz[j] += dz*ftmp;

#ifdef LOCALSTRESS
    if (stress_meas_enabled)
      localstressd_update(i, j, ftmp, tid);
#endif
#ifdef CONSTP_PBC
    constp_pbc_fxnextd[tid] += dx*dx*ftmp;
    constp_pbc_fynextd[tid] += dy*dy*ftmp;
    constp_pbc_fznextd[tid] += dz*dz*ftmp;
#endif
  }		// end if dr < RT 
}



////////////////////////////////////
// void calcrandf(...)
// ---------------
// add random forces for to particle forces
// ===
// Parameters:
//   long i : index of cell particle 1
//   long j : index of cell particle 2
//   long vlist_j: index of j in vlist
void calcrandf(long i, long j, long vlist_j)
{
  if (rnai[vlist_j]) {
    ////////////////////////////////////
    // random force
    // ------------
    // -> draw gaussian distributed random number
    // -> fr{xyz} = SIGMA*o^R(r_{ij})*xi/dr*d{xyz}*1/sqrt(dt)
    // xi is ftmp = sqrt(2*gammat_*kbt/dt)*(1.-dr/rt_)*xi/dr;
    double xi = rna[vlist_j];
    _fx[i] -= distx(i,j)*xi;
    _fy[i] -= disty(i,j)*xi;
    _fz[i] -= distz(i,j)*xi;

#ifdef PGAS
 #ifdef __DEBUG
    if (_spec[i] == gas_spec || _spec[j] == gas_spec) {
      if (fabs(xi) > 0.) {
        std::cout << "__DEBUG at t=" << curtime << ": calcrandf ftmp>0: i=" << i << ": j=" << j << std::endl;
      }
    }
 #endif
#endif
  }
}

////////////////////////////////////
// void calcbackgroundforces(...)
// ---------------
// calculate dissipative and random force
// on particle
// ===
// Parameters:
//   unsigned long i: index of cell particle
void calcbackgroundforces(long i)
{
  ////////////////////////////////////
  // background random force
  // ------------
  // -> draw gaussian distributed random number
  /*double sq_gbkbt = sqrt(2*GAMMAB*KBT/dt);
  _fx[i] += sq_gbkbt*gsl_ran_gaussian(rng, 1.0)/sqrt(dt);
  _fy[i] += sq_gbkbt*gsl_ran_gaussian(rng, 1.0)/sqrt(dt);
  _fz[i] += sq_gbkbt*gsl_ran_gaussian(rng, 1.0)/sqrt(dt);*/

  ////////////////////////////////////
  // background dissipation force
  calcbackgroundfd(i);

#ifdef HBC
  ////////////////////////////////////
  // in case of HBC add force field at z = L
  // f = 1/(z-z0)^8
  //_fz[i] += -1./( (_rz[i]-z0)*(_rz[i]-z0)*(_rz[i]-z0)*(_rz[i]-z0)*(_rz[i]-z0)*(_rz[i]-z0)*(_rz[i]-z0)*(_rz[i]-z0) );
#endif

}


////////////////////////////////////
// void calcbackgroundfd(...)
// ---------------
// calculate background dissipation once again with new velocity
// ===
// Parameters:
//   unsigned long i: index of cell particle
void calcbackgroundfd(long i)
{
  int index = _spec[i];
  ////////////////////////////////////
  // background dissipation force
  _fdx[i] += -gammab[index]*_vx[i];
  _fdy[i] += -gammab[index]*_vy[i];
  _fdz[i] += -gammab[index]*_vz[i];
}


