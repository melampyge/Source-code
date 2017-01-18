/*		localstress.cpp		*/
//////////////////////////////////////////
// functions related to determing the	//
// local stress				//
// IMPORTANT: Actually the pressure	//
// is determined and not the stress!!!	//
//////////////////////////////////////////



#include "localstress.h"


////////////////////////////////////
// Functions universal to VOL and MOP


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions specific for MOP
#ifdef LOCALSTRESSMOP

////////////////////////////////////
// void transform_coordinate_pair(...)
// ---------------
// fold coordinates into unit cell, adjust for periodic boundaries
// and put smaller into ri and larger into rj
// ---
// if PBC: the particle closer to the r[pli]=0 wall is fixed and
// called j, the other is moved around adjecent boxes and called i
// ===
// Parameters:
//   double *ri: coordinates of particle i
//   double *rj: coordinates of particle j
//   int pli   : index of plane (i.e. 0=x, 1=y, 2=z)
inline void transform_coordinate_pair(double *ri, double *rj, int pli)
{
  int i;
  double tmp;

  ri[0] = foldx(ri[0]);
  ri[1] = foldy(ri[1]);
  ri[2] = foldz(ri[2]);
  rj[0] = foldx(rj[0]);
  rj[1] = foldy(rj[1]);
  rj[2] = foldz(rj[2]);


  //////////////////////////////
  // check whether force is acting across periodic boundary
#if defined (PBC) || defined (HBC)
  if (fabs(ri[pli]-rj[pli]) > RPP) {
    if (ri[pli] > rj[pli])
      ri[pli] -= boxlength[pli];
    else
      rj[pli] -= boxlength[pli];
  }
#endif
  //////////////////////////////
  // now change i and j, if ri[pli] > rj[pli]
  if (ri[pli] > rj[pli]) {
    tmp = ri[0];
    ri[0] = rj[0];
    rj[0] = tmp;
    tmp = ri[1];
    ri[1] = rj[1];
    rj[1] = tmp;
    tmp = ri[2];
    ri[2] = rj[2];
    rj[2] = tmp;
  }

  //////////////////////////////
  // change other coordinates as well, if force acts across periodic boundary
#if defined (PBC) || defined (HBC)
  for (i = 0; i < 3; ++i) {
    if (i == pli)
      continue;
    if (fabs(ri[i]-rj[i]) > RPP) {
      if (rj[i] > ri[i])
        ri[i] += boxlength[i];
      else
        ri[i] -= boxlength[i];
    }
  }
#endif
}



////////////////////////////////////
// void localstress_update(...)
// ---------------
// find test area through which the force acts and
// update stress_fnext[] accordingly
// ===
// Parameters:
//   int i:           index of cell particle 1
//   int j:           index of cell particle 2 (i < j per "definition")
//   double ftmp:     force between i and j
//   int tid:         thread id
void localstress_update(int i, int j, double tdx, double tdy, double tdz, double ftmp, int tid)
{
  int l, k, o, __index;
  int ir[3];
  double ri[3];
  double rj[3];
  double dx[3];
  double m;
  double area[3];
  int ai = 0;

  if (_floor(_rx[i]/stress_binsize[0]) != _floor((_rx[i]+tdx)/stress_binsize[0]))
    area[ai++] = 0;
  if (_floor(_ry[i]/stress_binsize[1]) != _floor((_ry[i]+tdy)/stress_binsize[1]))
    area[ai++] = 1;
  if (_floor(_rz[i]/stress_binsize[2]) != _floor((_rz[i]+tdz)/stress_binsize[2]))
    area[ai++] = 2;


  for (o = 0; o < ai; ++o) {
    l = area[o];
    ri[0] = _rx[i];
    ri[1] = _ry[i];
    ri[2] = _rz[i];
    rj[0] = _rx[j];
    rj[1] = _ry[j];
    rj[2] = _rz[j];

    ir[0] = ir[1] = ir[2] = 0;

    transform_coordinate_pair(ri, rj, l);
    dx[0] = rj[0]-ri[0];
    dx[1] = rj[1]-ri[1];
    dx[2] = rj[2]-ri[2];

    ir[l] = (int)_floor(rj[l]/stress_binsize[l]);
    m = (ir[l]*stress_binsize[l] - ri[l])/dx[l];

    for (k = 0; k < 3; ++k) {
      if (k == l)
        continue;
      ir[k] = (int)_floor((ri[k]+m*dx[k])/stress_binsize[k]);
      //////////////////////////////////////
      // ensure that area lies within unitcell
      if (ir[k] < 0)
        ir[k] += stress_lengthi[k];
      else if (ir[k] >= stress_lengthi[k])
        ir[k] -= stress_lengthi[k];
    }

    __index = ir[0]*stress_lengthy*stress_lengthz*3 + ir[1]*stress_lengthz*3 + ir[2]*3 + l;
    __ASSERT(__index, 0, stress_length, "__index in localstress_update()");
    stress_fnext[tid][__index] += ftmp*dx[l];
  }
}



////////////////////////////////////
// void localstressd_update(...)
// ---------------
// find test area through which the force acts and
// update stress_fnextd[] accordingly
// ===
// Parameters:
//   int i:           index of cell particle 1
//   int j:           index of cell particle 2 (i < j per "definition")
//   double ftmp:     force between i and j
//   int tid:         thread id
void localstressd_update(int i, int j, double tdx, double tdy, double tdz, double ftmp, int tid)
{
  int l, k, o, __index;
  int ir[3];
  double ri[3];
  double rj[3];
  double dx[3];
  double m;
  double area[3];
  int ai = 0;

  if (_floor(_rx[i]/stress_binsize[0]) != _floor((_rx[i]+tdx)/stress_binsize[0]))
    area[ai++] = 0;
  if (_floor(_ry[i]/stress_binsize[1]) != _floor((_ry[i]+tdy)/stress_binsize[1]))
    area[ai++] = 1;
  if (_floor(_rz[i]/stress_binsize[2]) != _floor((_rz[i]+tdz)/stress_binsize[2]))
    area[ai++] = 2;


  for (o = 0; o < ai; ++o) {
    l = area[o];
    ri[0] = _rx[i];
    ri[1] = _ry[i];
    ri[2] = _rz[i];
    rj[0] = _rx[j];
    rj[1] = _ry[j];
    rj[2] = _rz[j];

    ir[0] = ir[1] = ir[2] = 0;

    transform_coordinate_pair(ri, rj, l);
    dx[0] = rj[0]-ri[0];
    dx[1] = rj[1]-ri[1];
    dx[2] = rj[2]-ri[2];

    ir[l] = (int)_floor(rj[l]/stress_binsize[l]);
    m = (ir[l]*stress_binsize[l] - ri[l])/dx[l];

    for (k = 0; k < 3; ++k) {
      if (k == l)
        continue;
      ir[k] = (int)_floor((ri[k]+m*dx[k])/stress_binsize[k]);
      //////////////////////////////////////
      // ensure that area lies within unitcell
      if (ir[k] < 0)
        ir[k] += stress_lengthi[k];
      else if (ir[k] >= stress_lengthi[k])
        ir[k] -= stress_lengthi[k];
    }

    __index = ir[0]*stress_lengthy*stress_lengthz*3 + ir[1]*stress_lengthz*3 + ir[2]*3 + l;
    __ASSERT(__index, 0, stress_length, "__index in localstressd_update()");
    stress_fnextd[tid][__index] += ftmp*dx[l];
  }
}



////////////////////////////////////
// void localstress_momentum_update(...)
// ---------------
// update momentum transfer through test areas
// ===
// Parameters:
//   int i   : index of cell particle
//   int tid : thread id
void localstress_momentum_update(int i, int tid)
{
  int k, l, __index;
  int binnew[3];
  int binold[3];
  int ir[3];
  double ri[3];
  double rj[3];
  double vel[3];
  double dx[3];
  double m;
  int bla;

  binnew[0] = (int)_floor(_rx[i]/stress_binsize[0]);
  binnew[1] = (int)_floor(_ry[i]/stress_binsize[1]);
  binnew[2] = (int)_floor(_rz[i]/stress_binsize[2]);
  binold[0] = (int)_floor(stress_rx[i]/stress_binsize[0]);
  binold[1] = (int)_floor(stress_ry[i]/stress_binsize[1]);
  binold[2] = (int)_floor(stress_rz[i]/stress_binsize[2]);
  vel[0] = _vx[i];
  vel[1] = _vy[i];
  vel[2] = _vz[i];

 

  for (l = 0; l < 3; ++l) {
    //////////////////////////////////////
    // ensure that force actualy penetrates this wall
    if (binnew[l] == binold[l])
      continue;

    ri[0] = stress_rx[i];
    ri[1] = stress_ry[i];
    ri[2] = stress_rz[i];
    rj[0] = _rx[i];
    rj[1] = _ry[i];
    rj[2] = _rz[i];

    ir[0] = ir[1] = ir[2] = 0;

    transform_coordinate_pair(ri, rj, l);
    dx[0] = rj[0]-ri[0];
    dx[1] = rj[1]-ri[1];
    dx[2] = rj[2]-ri[2];

    ir[l] = (int)_floor(rj[l]/stress_binsize[l]);
    m = (ir[l]*stress_binsize[l] - ri[l])/dx[l];

    for (k = 0; k < 3; ++k) {
      if (k == l)
        continue;
      ir[k] = (int)_floor((ri[k]+m*dx[k])/stress_binsize[k]);
      //////////////////////////////////////
      // ensure that area lies within unitcell
      if (ir[k] < 0)
        ir[k] += stress_lengthi[k];
      else if (ir[k] >= stress_lengthi[k])
        ir[k] -= stress_lengthi[k];
    }

    __index = ir[0]*stress_lengthy*stress_lengthz*3 + ir[1]*stress_lengthz*3 + ir[2]*3 + l;
    __ASSERT(__index, 0, stress_length, "__index in localstress_momentum_update()");
    stress_pnext[tid][__index] += (binnew[l]-binold[l])*(vel[l]);
  }
}




#endif	// ifdef LOCALSTRESSMOP



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Same functions as above but this time for
// volume measurement

#ifdef LOCALSTRESSVOL

////////////////////////////////////
// void transform_coordinate_pair(...)
// ---------------
// fold coordinates into unit cell, adjust for periodic boundaries
// and put smaller into ri and larger into rj
// ---
// if PBC: the particle closer to the r[pli]=0 wall is fixed and
// called j, the other is moved around adjecent boxes and called i
// ===
// Parameters:
//   double *ri: coordinates of particle i
//   double *rj: coordinates of particle j
inline void transform_coordinate_pair(double *ri, double *rj)
{
  ri[0] = foldx(ri[0]);
  ri[1] = foldy(ri[1]);
  ri[2] = foldz(ri[2]);
  rj[0] = foldx(rj[0]);
  rj[1] = foldy(rj[1]);
  rj[2] = foldz(rj[2]);

  //////////////////////////////
  // check whether force is acting across periodic boundary
#ifdef PBC
  for (int l = 0; l < 3; ++l)
    if (fabs(ri[l]-rj[l]) > RPP) {
      if (ri[l] > rj[l])
        ri[l] -= boxlength[l];
      else
        rj[l] -= boxlength[l];
    }
#endif
#ifdef HBC
  for (int l = 0; l < 2; ++l)
    if (fabs(ri[l]-rj[l]) > RPP) {
      if (ri[l] > rj[l])
        ri[l] -= boxlength[l];
      else
        rj[l] -= boxlength[l];
    }
#endif
}


////////////////////////////////////
// void localstress_update(...)
// ---------------
// find test area through which the force acts and
// update stress_fnext[] accordingly
// ===
// Parameters:
//   int i:           index of cell particle 1
//   int j:           index of cell particle 2 (i < j per "definition")
//   double ftmp:     force between i and j
//   int tid:         thread id
void localstress_update(int i, int j, double ftmp, int tid)
{
  int l, k, __index;
  double ri[3];
  double rj[3];
  double dx[3];
  int bini;
  int binj;
  int ir[3];
  double rbox, tmprev;
  std::vector<double> vtm;


  ///////////////////////////////////////////
  // Transform into unit cell or PBC/HBC
  ri[0] = _rx[i];
  ri[1] = _ry[i];
  ri[2] = _rz[i];
/*#ifdef MEASUREINS	// insert something like this for MEASUREINS
  ri[2] = stress_maxz-_rz[i];
#else
  ri[2] = _rz[i];
#endif*/
  rj[0] = _rx[j];
  rj[1] = _ry[j];
  rj[2] = _rz[j];

#if defined (PBC) || defined (HBC)
  transform_coordinate_pair(ri, rj);
#endif

  for (l = 0; l < 3; ++l)
    dx[l] = rj[l] - ri[l];

  ///////////////////////////////////////////
  // Check x/y/z walls
  for (l = 0; l < 3; ++l) {
    bini = (int)(ri[l]*stress_ibinsize[l]);
    binj = (int)(rj[l]*stress_ibinsize[l]);
    if (binj < bini) {
      int __tmp = binj;
      binj = bini;
      bini = __tmp;
    }
    for (k = bini; k < binj; ++k)
      vtm.push_back(((k+1)*stress_binsize[l]-ri[l])/dx[l]);
  }
  vtm.push_back(1.);		// necessary

  ///////////////////////////////////////////
  // sort vtm in order of ascending size
  std::sort(vtm.begin(), vtm.end());

  
  ///////////////////////////////////////////
  // Now walk along dx and add force contributions accordingly
  // vtm holds the segments of dx, which are in different boxes
  // Easiest case: vtm[0]==1., dx lies completely inside one box
  tmprev = 0.;
  for (std::vector<double>::const_iterator it = vtm.begin(); it != vtm.end(); ++it) {
    for (l = 0; l < 3; ++l) {
      rbox  = ri[l]+0.5*(*it+tmprev)*dx[l];		// calc box of current segment
      ir[l] = _floor(rbox*stress_ibinsize[l]);
      if (ir[l] < 0)   // in case of periodic boundaries...
        ir[l] += stress_lengthi[l];
      else if (ir[l] >= stress_lengthi[l])
        ir[l] -= stress_lengthi[l];
    }

    for (k = 0; k < 3; ++k) {
    __index = ir[0]*stress_lengthy*stress_lengthz*9 + ir[1]*stress_lengthz*9 + ir[2]*9 + 3*k;
      for (l = 0; l < 3; ++l) {
        __ASSERT(__index+l, 0, stress_length, "__index in localstress_update()");
        stress_fnext[tid][__index+l] += ftmp*dx[l]*(*it-tmprev)*dx[k];
      }
    }
    tmprev = *it;
  }
}



////////////////////////////////////
// void localstressd_update(...)
// ---------------
// find test area through which the force acts and
// update stress_fnextd[] accordingly
// ===
// Parameters:
//   int i:           index of cell particle 1
//   int j:           index of cell particle 2 (i < j per "definition")
//   double ftmp:     force between i and j
//   int tid:         thread id
void localstressd_update(int i, int j, double ftmp, int tid)
{
  int l, k, __index;
  double ri[3];
  double rj[3];
  double dx[3];
  int bini;
  int binj;
  int ir[3];
  double rbox, tmprev;
  std::vector<double> vtm;


  ///////////////////////////////////////////
  // Transform into unit cell or PBC/HBC
  ri[0] = _rx[i];
  ri[1] = _ry[i];
  ri[2] = _rz[i];
/*#ifdef MEASUREINS	// insert something like this for MEASUREINS
  ri[2] = stress_maxz-_rz[i];
#else
  ri[2] = _rz[i];
#endif*/
  rj[0] = _rx[j];
  rj[1] = _ry[j];
  rj[2] = _rz[j];

#if defined (PBC) || defined (HBC)
  transform_coordinate_pair(ri, rj);
#endif

  for (l = 0; l < 3; ++l)
    dx[l] = rj[l] - ri[l];

  ///////////////////////////////////////////
  // Check x/y/z walls
  for (l = 0; l < 3; ++l) {
    bini = (int)(ri[l]*stress_ibinsize[l]);
    binj = (int)(rj[l]*stress_ibinsize[l]);
    if (binj < bini) {
      int __tmp = binj;
      binj = bini;
      bini = __tmp;
    }
    for (k = bini; k < binj; ++k)
      vtm.push_back(((k+1)*stress_binsize[l]-ri[l])/dx[l]);
  }
  vtm.push_back(1.);		// necessary

  ///////////////////////////////////////////
  // sort tm[l] in order of ascending size
  std::sort(vtm.begin(), vtm.end());


  ///////////////////////////////////////////
  // Now walk along dx and add force contributions accordingly
  // vtm holds the pieces of dx, which are in different boxes
  tmprev = 0.;
  for (std::vector<double>::const_iterator it = vtm.begin(); it != vtm.end(); ++it) {
    for (l = 0; l < 3; ++l) {
      rbox  = ri[l]+0.5*(*it+tmprev)*dx[l];
      ir[l] = _floor(rbox*stress_ibinsize[l]);
      if (ir[l] < 0)   // in case of periodic boundaries...
        ir[l] += stress_lengthi[l];
      else if (ir[l] >= stress_lengthi[l])
        ir[l] -= stress_lengthi[l];
    }

    for (k = 0; k < 3; ++k) {
    __index = ir[0]*stress_lengthy*stress_lengthz*9 + ir[1]*stress_lengthz*9 + ir[2]*9 + 3*k;
      for (l = 0; l < 3; ++l) {
        __ASSERT(__index+l, 0, stress_length, "__index in localstressd_update()");
        stress_fnextd[tid][__index+l] += ftmp*dx[l]*(*it-tmprev)*dx[k];
      }
    }
    tmprev = *it;
  }
}



////////////////////////////////////
// void localstress_momentum_update(...)
// ---------------
// update momentum transfer through test areas
// ===
// Parameters:
//   int i   : index of cell particle
//   int tid : thread id
void localstress_momentum_update(int i, int tid)
{
  int l, k, __index;
  double ri[3];
  double rj[3];
  double vel[3];
  double dx[3];
  int bini;
  int binj;
  int ir[3];
  double rbox, tmprev;
  std::vector<double> vtm;

  vel[0] = _vx[i];
  vel[1] = _vy[i];
  vel[2] = _vz[i];

  ri[0] = stress_rx[i];
  ri[1] = stress_ry[i];
  ri[2] = stress_rz[i];
  rj[0] = _rx[i];
  rj[1] = _ry[i];
  rj[2] = _rz[i];


#if defined (PBC) || defined (HBC)
  transform_coordinate_pair(ri, rj);
#endif

  for (l = 0; l < 3; ++l)
    dx[l] = rj[l] - ri[l];

  ///////////////////////////////////////////
  // Check x/y/z walls
  for (l = 0; l < 3; ++l) {
    bini = (int)(ri[l]*stress_ibinsize[l]);
    binj = (int)(rj[l]*stress_ibinsize[l]);
    if (binj < bini) {
      int __tmp = binj;
      binj = bini;
      bini = __tmp;
    }
    for (k = bini; k < binj; ++k)
      vtm.push_back(((k+1)*stress_binsize[l]-ri[l])/dx[l]);
  }
  vtm.push_back(1.);		// necessary

  ///////////////////////////////////////////
  // sort tm[l] in order of ascending size
  std::sort(vtm.begin(), vtm.end());


  ///////////////////////////////////////////
  // Now walk along dx and add force contributions accordingly
  // vtm holds the pieces of dx, which are in different boxes
  tmprev = 0.;
  for (std::vector<double>::const_iterator it = vtm.begin(); it != vtm.end(); ++it) {
    for (l = 0; l < 3; ++l) {
      rbox  = ri[l]+0.5*(*it+tmprev)*dx[l];
      ir[l] = _floor(rbox*stress_ibinsize[l]);
      if (ir[l] < 0)   // in case of periodic boundaries...
        ir[l] += stress_lengthi[l];
      else if (ir[l] >= stress_lengthi[l])
        ir[l] -= stress_lengthi[l];
    }

    for (k = 0; k < 3; ++k) {
    __index = ir[0]*stress_lengthy*stress_lengthz*9 + ir[1]*stress_lengthz*9 + ir[2]*9 + 3*k;
      for (l = 0; l < 3; ++l) {
        __ASSERT(__index+l, 0, stress_length, "__index in localstress_momentum_update()");
        stress_pnext[tid][__index+l] += m[_spec[i]]*vel[l]*(*it-tmprev)*dx[k];
      }
    }
    tmprev = *it;
  }
}



#endif	// ifdef LOCALSTRESSVOL








