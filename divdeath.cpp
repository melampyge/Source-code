/*		divdeath.cpp		*/
//////////////////////////////////////////
// cell division and death		//
//////////////////////////////////////////


#include "divdeath.h"




////////////////////////////////////
// void cell_death(...)
// ---------------
// check for cell death
// ===
// Parameters:
void cell_death()
{
  long i, k;

#ifdef DENSITYPROFILE
  int __index;
  //////////////////////////////////
  // find maximum z and calc current density
  std::fill(dens_curn, dens_curn+dens_length, 0);
  std::fill(dens_curnka, dens_curnka+dens_length, 0);
  dens_maxz = dens_p_maxz = -1.;
  dens_minz = dens_p_minz = LZ+1.;
  dens_com_rx = dens_com_ry = dens_com_rz = 0.;
  dens_com_n = 0;
  for (i = 0; i < lastcn; ++i)
    if (_list[i] == -1
 #ifdef FIXEDCELLS
                       && _spec[i*_pperc] != fixed_spec
 #endif
 #ifdef PGAS
                       && _spec[i*_pperc] != gas_spec
 #endif
                      ) {
 #ifdef RECTANGULAR
      if ((_rz[i*_pperc]+_rz[i*_pperc+1])/2. > dens_maxz)
        dens_maxz = (_rz[i*_pperc]+_rz[i*_pperc+1])/2.;
      if ((_rz[i*_pperc]+_rz[i*_pperc+1])/2. < dens_minz)
        dens_minz = (_rz[i*_pperc]+_rz[i*_pperc+1])/2.;
 #endif
      for (k=i*_pperc; k < (i+1)*_pperc; ++k) {
        dens_com_rx += _rx[k];
        dens_com_ry += _ry[k];
        dens_com_rz += _rz[k];
        dens_com_n++;
 #ifdef RECTANGULAR
        if (_rz[k] > dens_p_maxz)
          dens_p_maxz = _rz[k];
        if (_rz[k] < dens_p_minz)
          dens_p_minz = _rz[k];
 #endif
      }
    }
  dens_com_rx /= dens_com_n;
  dens_com_ry /= dens_com_n;
  dens_com_rz /= dens_com_n;

 #ifdef SPHERICAL
  for (i = 0; i < lastcn; ++i)
    if (_list[i] == -1
  #ifdef FIXEDCELLS
                       && _spec[i*_pperc] != fixed_spec
  #endif
  #ifdef PGAS
                       && _spec[i*_pperc] != gas_spec
  #endif
                      ) {
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
 #endif


  for (i = 0; i < lastcn; ++i)
    if (_list[i] == -1
 #ifdef FIXEDCELLS
                       && _spec[i*_pperc] != fixed_spec
 #endif
 #ifdef PGAS
                       && _spec[i*_pperc] != gas_spec
 #endif
                      ) {
      __index = GET__INDEX(i, dens_binsize);
      __ASSERT(__index, 0, dens_length, "__index0 in cell_death()");
      dens_curn[__index]++;
    }
#endif		// ifdef DENSITYPROFILE

  for (i = 0; i < lastcn; ++i) {	// count over cells not particles!
    if (_list[i] != -1
#ifdef PGAS
                       || _spec[i*_pperc] == gas_spec
#endif
                      )		// fixed cells have ka==0, so no need to check here but checking
      continue;			// for pgas should improve performance (especially, with high pressure)
    if (drand[0]() < ka[_spec[i*_pperc]]*dt) {
#ifdef DENSITYPROFILE
      __index = GET__INDEX(i, dens_binsize);
      __ASSERT(__index, 0, dens_length, "__index1 in cell_death()");
      dens_nka[__index]++;
      dens_curnka[__index]++;
#endif
      ////////////////////////////////
      // delete cell
      delete_cell(_pperc*i);
    }
  }

#ifdef DENSITYPROFILE
  for (i = 0; i < dens_length; ++i)
    if (dens_curn[i] > 0) {
      dens_knmeas[i]++;
      if (dens_curnka[i] > 0)
        dens_ka[i] += (double)dens_curnka[i]/dens_curn[i];
    }
#endif
}



////////////////////////////////////
// void cell_division(...)
// ---------------
// check for cell division
// ===
// Parameters:
void cell_division()
{
  double drsq, tmp;
  long tmp_lastcn = lastcn;
  long i;

#ifdef DENSITYPROFILE
  int __index;
  //////////////////////////////////
  // find maximum z and calc current density (may have changed de to cell death)
  std::fill(dens_curnkd, dens_curnkd+dens_length, 0);
#endif

  for (i = 0; i < tmp_lastcn; ++i) {	// count over cells not particles!
    if (_list[i] != -1
#ifdef PGAS
                       || _spec[i*_pperc] == gas_spec
#endif
                      )
      continue;


    drsq = distsq(i*_pperc, i*_pperc+1);
    if (drsq > rct[_spec[i*_pperc]]*rct[_spec[i*_pperc]]) {
#ifdef DENSITYPROFILE
      __index = GET__INDEX(i, dens_binsize);
      __ASSERT(__index, 0, dens_length, "__index1 in cell_division()");
      dens_nkd[__index]++;
      {
        double tmprx  = (_rx[i*_pperc]-_rx[i*_pperc+1])/2.;
        double tmpry  = (_ry[i*_pperc]-_ry[i*_pperc+1])/2.;
        double tmprz  = (_rz[i*_pperc]-_rz[i*_pperc+1])/2.;
        double tmprpx = (_rx[i*_pperc]+_rx[i*_pperc+1])/2.;
        double tmprpy = (_ry[i*_pperc]+_ry[i*_pperc+1])/2.;
        double tmprpz = (_rz[i*_pperc]+_rz[i*_pperc+1])/2.;
        tmp = CALC_THETA(tmprx, tmpry, tmprz, dens_com_rx-tmprpx, dens_com_ry-tmprpy, dens_com_rz-tmprpz);
      }
      dens_kperp[__index] += tmp;
      dens_kpara[__index] += 1. - tmp;

      ////////////////////////////////////////
      // calculate nematic order parameter
      dens_sq[__index] += 3./2.*tmp*tmp-0.5;

      dens_curnkd[__index]++;

#endif
      //////////////////////////////////////////
      // perform division
      do_division(i*_pperc);
    }
  }

#ifdef DENSITYPROFILE
  for (i = 0; i < dens_length; ++i)
    if (dens_curn[i]>0 && dens_curnkd[i]>0.)
      dens_kd[i] += (double)dens_curnkd[i]/dens_curn[i];
#endif

  if (vlist_rebuild) {
    update_verlet_list();
    vlist_rebuild = 0;
  }
}



////////////////////////////////////
// int do_division(...)
// ---------------
// do cell division
// ===
// Parameters:
//   long i: index of first cell particle (i+1 is second one)
int do_division(long i)
{
  double tmpx, tmpy, tmpz;
  long j;
  long tmp, tmphead, tmplastn;


#ifdef __DEBUG_1
  std::cout << "__DEBUG at t=" << curtime << ": cell insertion started: i=" << i << ": j=" << _pperc*_head << std::endl;
#endif


  /////////////////////////////////
  // TODO
  // -> copy i+1 to _pperc*_head
  // -> copy i to i+1
  // -> create new particle at i+1
  // -> create new particle at _pperc*_head+1

  /////////////////////////////////
  // check whether there is space left for divison
  if (_head < 0) {
    // array full, for now abort
    process_signal_handler(0xFFFF);
    return -1;
  }

  /////////////////////////////////
  // -> copy i+1 to _pperc*_head and _pperc*_head+1
  j = _pperc*_head;
  _rx[j]  = _rx[i+1];
  _ry[j]  = _ry[i+1];
  _rz[j]  = _rz[i+1];
  _vx[j]  = _vx[i+1];
  _vy[j]  = _vy[i+1];
  _vz[j]  = _vz[i+1];
  _fx[j]  = _fx[i+1];
  _fy[j]  = _fy[i+1];
  _fz[j]  = _fz[i+1];
  _fdx[j] = _fdx[i+1];
  _fdy[j] = _fdy[i+1];
  _fdz[j] = _fdz[i+1];
  _spec[j] = _spec[i+1];
  _id[j]   = _id[i+1];

#ifdef REALFLUX
  rflux_crz[_head] = rflux_crz[i/_pperc];
#endif
#ifdef LOCALSTRESS
  stress_rx[j] = _rx[i+1];
  stress_ry[j] = _ry[i+1];
  stress_rz[j] = _rz[i+1];
#endif


  ++j;
  _rx[j]  = _rx[i+1];
  _ry[j]  = _ry[i+1];
  _rz[j]  = _rz[i+1];
  _vx[j]  = _vx[i+1];
  _vy[j]  = _vy[i+1];
  _vz[j]  = _vz[i+1];
  _fx[j]  = _fx[i+1];
  _fy[j]  = _fy[i+1];
  _fz[j]  = _fz[i+1];
  _fdx[j] = _fdx[i+1];
  _fdy[j] = _fdy[i+1];
  _fdz[j] = _fdz[i+1];
  _spec[j] = _spec[i+1];
  _id[j]   = next_free_id++;

#ifdef LOCALSTRESS
  stress_rx[j] = _rx[i+1];
  stress_ry[j] = _ry[i+1];
  stress_rz[j] = _rz[i+1];
#endif


  /////////////////////////////////
  // place new particles at random distance < RC
  // in case of BBC or RBC check if they are inside box!
#if defined (BBC) || defined (RBC) || defined (HBC)
  tmpx = RC*(drand[0]()-0.5);
  tmpy = RC*(drand[0]()-0.5);
  tmpz = RC*(drand[0]()-0.5);
  _rx[j] += tmpx;
  _ry[j] += tmpy;
  _rz[j] += tmpz;
 #ifdef __DEBUG
  int __lc = 0;
 #endif
  while (outofBound(j)) {
    _rx[j] -= tmpx;
    _ry[j] -= tmpy;
    _rz[j] -= tmpz;
    tmpx = RC*(drand[0]()-0.5);
    tmpy = RC*(drand[0]()-0.5);
    tmpz = RC*(drand[0]()-0.5);
    _rx[j] += tmpx;
    _ry[j] += tmpy;
    _rz[j] += tmpz;
 #ifdef __DEBUG
    __lc++;
    if (__lc>1000) {
      std::cout << "__DEBUG at t=" << curtime << ": possible endless loop1: j=" << j << std::endl;
      std::cout << "_rx[j]=" << _rx[j]-tmpx << ": _ry[j]=" << _ry[j]-tmpy << ": _rz[j]=" << _rz[j]-tmpz << std::endl;
      process_signal_handler(0xFFFFFFFF);
    }
 #endif
  }
#endif		// ifdef BBC || RBC || HBC
#ifdef PBC
  _rx[j] += RC*(drand[0]()-0.5);
  _ry[j] += RC*(drand[0]()-0.5);
  _rz[j] += RC*(drand[0]()-0.5);
#endif

  ///////////////////////////////////
  // i+1 moved to new position and //
  // linked particle was created   //
  ///////////////////////////////////

  /////////////////////////////////
  // -> copy i to i+1
  j = i+1;
  _rx[j]  = _rx[i];
  _ry[j]  = _ry[i];
  _rz[j]  = _rz[i];
  _vx[j]  = _vx[i];
  _vy[j]  = _vy[i];
  _vz[j]  = _vz[i];
  _fx[j]  = _fx[i];
  _fy[j]  = _fy[i];
  _fz[j]  = _fz[i];
  _fdx[j] = _fdx[i];
  _fdy[j] = _fdy[i];
  _fdz[j] = _fdz[i];
  _spec[j] = _spec[i];
  _id[j]   = next_free_id++;

#ifdef LOCALSTRESS
  stress_rx[j] = _rx[i];
  stress_ry[j] = _ry[i];
  stress_rz[j] = _rz[i];
#endif


  /////////////////////////////////
  // place new particles at random distance < RC
  // in case of BBC or RBC check if they are inside box!
#if defined (BBC) || defined (RBC) || defined (HBC)
  tmpx = RC*(drand[0]()-0.5);
  tmpy = RC*(drand[0]()-0.5);
  tmpz = RC*(drand[0]()-0.5);
  _rx[j] += tmpx;
  _ry[j] += tmpy;
  _rz[j] += tmpz;
#ifdef __DEBUG
  __lc = 0;
#endif
  while (outofBound(j)) {
    _rx[j] -= tmpx;
    _ry[j] -= tmpy;
    _rz[j] -= tmpz;
    tmpx = RC*(drand[0]()-0.5);
    tmpy = RC*(drand[0]()-0.5);
    tmpz = RC*(drand[0]()-0.5);
    _rx[j] += tmpx;
    _ry[j] += tmpy;
    _rz[j] += tmpz;
#ifdef __DEBUG
    __lc++;
    if (__lc>1000) {
      std::cout << "__DEBUG at t=" << curtime << ": possible endless loop2: j=" << j << std::endl;
      std::cout << "_rx[j]=" << _rx[j]-tmpx << ": _ry[j]=" << _ry[j]-tmpy << ": _rz[j]=" << _rz[j]-tmpz << std::endl;
      process_signal_handler(0xFFFFFFFF);
    }
#endif
  }
#endif		// ifdef BBC || RBC || HBC
#ifdef PBC
  _rx[j] += RC*(drand[0]()-0.5);
  _ry[j] += RC*(drand[0]()-0.5);
  _rz[j] += RC*(drand[0]()-0.5);
#endif

  /////////////////////////////////
  // adjust _list[]
  tmplastn = lastn;
  tmphead = _head;
  tmp = _list[_head];
  _list[_head] = -1;
  _head = tmp;
  if (lastcn <= tmphead)
    lastcn = tmphead+1;
  verlet_list_fill_offset(lastn, _pperc*lastcn);
  lastn = _pperc*lastcn;

  if (!vlist_rebuild) {
    verlet_list_move_entry(i+1, _pperc*tmphead, tmplastn);
    if (!vlist_rebuild) {
      verlet_list_insert_i(i+1, tmplastn);
      // The following does not really make a difference but for
      // correctness, we add RC to vlist_dx
      vlist_dx[i+1] += RC;
      vlist_dy[i+1] += RC;
      vlist_dz[i+1] += RC;
      if (!vlist_rebuild) {
        verlet_list_insert_i(_pperc*tmphead+1, tmplastn);
	// The following does not really make a difference but for
        // correctness, we add RC to vlist_dx
        vlist_dx[_pperc*tmphead+1] += RC;
        vlist_dy[_pperc*tmphead+1] += RC;
        vlist_dz[_pperc*tmphead+1] += RC;
      }
    }
  }

#ifdef __DEBUG_1
  std::cout << "__DEBUG at t=" << curtime << ": cell insertion ended: i=" << i << std::endl;
#endif

  return 0;
}



////////////////////////////////////
// int delete_cell(...)
// ---------------
// delete cell particles from linked list
// ===
// Parameters:
//   long i: index of first cell particle (i+1 is second one)
// ===
// Returns:
//   nothing important
int delete_cell(long i)
{
  long j = i/_pperc;
#ifdef __DEBUG_1
  std::cout << "__DEBUG at t=" << curtime << ": cell deletion started: i=" << i << std::endl;
#endif

  /////////////////////////////////
  // check whether array is full
  if (_head < 0) {
    // array full, for now abort
    process_signal_handler(0xFFFF);
    return -1;
  }

  /////////////////////////////////////
  // check whether cell is really existent
  if (_list[j] != -1) {
    process_signal_handler(0xFFFE);
    return -1;
  }

  /////////////////////////////////////
  // delete cell
  _list[j] = _head;
  _head = j;

  while (lastcn > 0 && _list[lastcn-1] != -1)
    --lastcn;
  lastn = lastcn*_pperc;

  verlet_list_delete_i(i);
  verlet_list_delete_i(i+1);

#ifdef __DEBUG_1
  std::cout << "__DEBUG at t=" << curtime << ": cell deletion ended: i=" << i << std::endl;
#endif

  return 0;
}



