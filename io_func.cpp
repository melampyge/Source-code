/*                  io_func.cpp             */
//////////////////////////////////////////////
// contains all input/output realted        //
// functions to easily change format specs  //
//////////////////////////////////////////////

#include "io_func.h"


using namespace std;


/////////////////////////////////////
// init_output(...)
// -----------------
// initialize output, i.e. empty existing
// output files
// ===
// Parameters:
void init_output()
{
  ofstream ftmp;

  /////////////////////////////////////
  // open trajectory file
  ftmp.open((string(OUTDIR)+string(FTRAJNAME)).c_str(), ios::trunc);
  if (!ftmp) {
    cerr << "Error while fopen('" << OUTDIR << FTRAJNAME << "') in init_output()" << endl;
    exit(-1);
  }
  ftmp.close();

  /////////////////////////////////////
  // open number of cells file
  ftmp.open((string(OUTDIR)+string(FNUMCELLSNAME)).c_str(), ios::trunc);
  if (!ftmp) {
    cerr << "Error while fopen('" << OUTDIR << FNUMCELLSNAME << "') in init_output()" << endl;
    exit(-1);
  }
  ftmp.close();

#ifdef PMEASUREMENT
  ///////////////////////////////////
  // delete pressure.dat on startup
  ftmp.open((string(OUTDIR)+string(FPRESSNAME)).c_str(), ios::trunc);
  if (!ftmp) {
    cerr << "Error while fopen('" << OUTDIR << FPRESSNAME << "') in init_output()" << endl;
    exit(-1);
  }
  ftmp.close();
#endif
#ifdef DENSITYPROFILE
  ///////////////////////////////////
  // delete density.dat on startup
  ftmp.open((string(OUTDIR)+string(FDENSITYNAME)).c_str(), ios::trunc);
  if (!ftmp) {
    cerr << "Error while fopen('" << OUTDIR << FDENSITYNAME << "') in init_output()" << endl;
    exit(-1);
  }
  ftmp.close();
#endif
#ifdef LOCALSTRESSMOP
  ftmp.open((string(OUTDIR)+string(FLOCALSTRESSMOPNAME)).c_str(), ios::trunc);
  if (!ftmp) {
    cerr << "Error while fopen('" << OUTDIR << FLOCALSTRESSMOPNAME << "') in init_output()" << endl;
    exit(-1);
  }
  ftmp.close();
#endif
#ifdef LOCALSTRESSVOL
  ftmp.open((string(OUTDIR)+string(FLOCALSTRESSVOLNAME)).c_str(), ios::trunc);
  if (!ftmp) {
    cerr << "Error while fopen('" << OUTDIR << FLOCALSTRESSVOLNAME << "') in init_output()" << endl;
    exit(-1);
  }
  ftmp.close();
#endif
#ifdef CONSTP_PBC
  ftmp.open((string(OUTDIR)+string(FCONSTPPBCNAME)).c_str(), ios::trunc);
  if (!ftmp) {
    cerr << "Error while fopen('" << OUTDIR << FCONSTPPBCNAME << "') in init_output()" << endl;
    exit(-1);
  }
  ftmp.close();
#endif
#ifdef __DEBUG_TIMER
  ftmp.open((string(OUTDIR)+string("__timing.dat")).c_str(), ios::trunc);
  if (!ftmp) {
    cerr << "Error while fopen(__timing.dat) in init_output()" << endl;
    exit(-1);
  }
  ftmp.close();
#endif
}



/////////////////////////////////////
// output_trajectory(...)
// -----------------
// output trajectory for all particles
// ===
// Parameters:
void output_trajectory()
{
  long i,j;
  ofstream ftraj((string(OUTDIR)+string(FTRAJNAME)).c_str(), ios::app);
  if (!ftraj) {
    cerr << "Error while fopen('" << OUTDIR << FTRAJNAME << "') in init_output()" << endl;
    exit(-1);
  }

  for (j = 0; j < lastcn; ++j) {
    if (_list[j] != -1)
      continue;
    for (i = _pperc*j; i < _pperc*(j+1); ++i) {
      if (_spec[i] == -1)
        continue;
      ftraj << setprecision(10) << curtime << "\t" << _id[i] << "\t" << _rx[i] << "\t" << _ry[i] << "\t" << _rz[i] << "\t"
            << _spec[i] << "\t" << _vx[i] << "\t" << _vy[i] << "\t" << _vz[i] << "\n";
    }
  }
  /////////////////////////////////////
  // this also ensures that data is actually written
  ftraj.close();
}



/////////////////////////////////////
// output_numcells(...)
// -----------------
// output number of living cells
// ===
// Parameters:
void output_numcells()
{
  long count[MAXSPECNUM];
  ofstream fnumcells((string(OUTDIR)+string(FNUMCELLSNAME)).c_str(), ios::app);

  fill(count, count+MAXSPECNUM, 0);
  if (!fnumcells) {
    cout << "Error while fopen('" << OUTDIR << FNUMCELLSNAME << "')" << endl;
    exit(-1);
  }
  for (long i = 0; i < lastcn; ++i)
    if (_list[i] == -1) ++count[_spec[i*_pperc]];

  fnumcells << curtime;
  for (int i = 0; i < MAXSPECNUM; ++i)
    fnumcells << "\t" << count[i];
  fnumcells << endl;
  fnumcells.close();
}



#ifdef PMEASUREMENT
/////////////////////////////////////
// output_pressure(...)
// -----------------
// output measured pressure
// ===
// Parameters:
void output_pressure()
{
  ofstream fpress((string(OUTDIR)+string(FPRESSNAME)).c_str(), ios::app);
  double mp[6] = {0., 0., 0., 0., 0., 0.};

  for (int i = 0; i < TARGETNUMTHREADS; ++i) {
    mp[0] += pm_px0[i];
    mp[1] += pm_pxL[i];
    mp[2] += pm_py0[i];
    mp[3] += pm_pyL[i];
    mp[4] += pm_pz0[i];
    mp[5] += pm_pzL[i];
  }

  if (!fpress) {
    cerr << "Error while fopen('" << OUTDIR << FPRESSNAME << "') at t=" << curtime << endl;
    exit(-1);
  }

  fpress << curtime << "\t" << mp[0]/(LY*LZ*dt*PMEANSTEP) << "\t" << mp[1]/(LY*LZ*dt*PMEANSTEP) << "\t"
         << mp[2]/(LX*LZ*dt*PMEANSTEP) << "\t" << mp[3]/(LX*LZ*dt*PMEANSTEP) << "\t"
         << mp[4]/(LX*LY*dt*PMEANSTEP) << "\t" << mp[5]/(LX*LY*dt*PMEANSTEP);
#ifdef FIXEDCELLS
  double mpf[2] = {0., 0.};
  for (int i = 0; i < TARGETNUMTHREADS; ++i) {
    mpf[0] += fixed_pkz0[i];
    mpf[1] += fixed_pkzL[i];
  }
  fpress << "\t" << mpf[0]/(LX*LY*PMEANSTEP) << "\t" << mpf[1]/(LX*LY*PMEANSTEP);
#endif
  fpress << endl;

  fpress.close();
}
#endif


#ifdef DENSITYPROFILE
/////////////////////////////////////
// output_density(...)
// -----------------
// output density
// ===
// Parameters:
void output_density()
{
  ofstream fdens((string(OUTDIR)+string(FDENSITYNAME)).c_str(), ios::app);
  double rho_i = 0.;
  double vol_i = 0.;
  double area_i = 0.;

  if (!fdens) {
    cout << "Error while fopen('" << OUTDIR << FDENSITYNAME << "')" << endl;
    exit(-1);
  }
  for (long i = 0; i < dens_length; ++i) {
 #ifdef RECTANGULAR
    vol_i  = LX*LY*dens_binsize;
    area_i = LX*LY;
 #elif defined (SPHERICAL)
  #ifdef MEASUREINS
    int ip = dens_mean_maxz/(DENSDATASTEP*dens_binsize) - i;
    if (ip < 0)
      vol_i = 0.;
    else
      vol_i  = 4./3.*3.1415926536*dens_binsize*dens_binsize*dens_binsize*((ip+1)*(ip+1)*(ip+1)-ip*ip*ip);
    area_i = 4.*3.1415926536*ip*ip*dens_binsize*dens_binsize;
  #elif defined (MEASUREINR)
    vol_i  = 4./3.*3.1415926536*dens_binsize*dens_binsize*dens_binsize*((i+1)*(i+1)*(i+1)-i*i*i);
    area_i = 4.*3.1415926536*i*i*dens_binsize*dens_binsize;
  #endif
 #endif
    fdens << curtime << "\t" << i*dens_binsize;		// 1, 2
    if (dens_nmeas[i] > 0 && vol_i > 0.) {
      rho_i = (double)dens_rhoV[i]/(vol_i*dens_nmeas[i]);
      fdens << "\t" << setprecision(12) << rho_i;		// 3
    }
    else {
      rho_i = 0.;
      fdens << "\t" << rho_i;				// 3
    }
    //////////////////////////////////////
    // kd/ka measurement does only make sense, when there are any cells at all
    if (dens_rhoV[i] > 0) 					// 4
      fdens << "\t" << setprecision(12) << (double)dens_nkd[i]/(dens_rhoV[i]*dt);	// nkd/(rho*LX*LY*dens_binsize*dt)
    else
      fdens << "\t" << -0.;					// 4
    if (dens_rhoV[i]+dens_nka[i] > 0)				// 5
      fdens << "\t" << setprecision(12) << (double)dens_nka[i]/((dens_rhoV[i]+dens_nka[i])*dt);
    else							// 5
      fdens << "\t" << -0.;
    fdens << "\t" << dens_nmeas[i];				// 6
    if (dens_nkd[i] > 0) {					// 7, 8
      fdens << "\t" << dens_kpara[i]/dens_nkd[i]
                << "\t" << dens_kperp[i]/dens_nkd[i];
    }
    else							// 7, 8
      fdens << "\t" << 0.0 << "\t" << 0.0;

    if (dens_nmeas[i] > 0 && vol_i > 0.)
      fdens << "\t" << dens_rhoVsq[i]/(vol_i*vol_i*dens_nmeas[i]);	// 9
    else
      fdens << "\t" << 0.;					// 9

    fdens << "\t" << dens_mean_maxz/DENSDATASTEP;		// 10
    if (dens_lastmean_maxz > 0.)
      fdens << "\t" << (dens_mean_maxz-dens_lastmean_maxz)/(dt*DENSDATASTEP*DENSDATASTEP);	// 11
    else
      fdens << "\t" << 0.;					// 11
    fdens << "\t" << dens_mean_p_maxz/DENSDATASTEP;		// 12
    if (dens_lastmean_p_maxz > 0.)
      fdens << "\t" << (dens_mean_p_maxz-dens_lastmean_p_maxz)/(dt*DENSDATASTEP*DENSDATASTEP);	// 13
    else
      fdens << "\t" << 0.;					// 13

    fdens << "\t" << dens_mean_minz/DENSDATASTEP;		// 14
    if (dens_lastmean_minz > 0.)
      fdens << "\t" << (dens_mean_minz-dens_lastmean_minz)/(dt*DENSDATASTEP*DENSDATASTEP);	// 15
    else
      fdens << "\t" << 0.;					// 15
    fdens << "\t" << dens_mean_p_minz/DENSDATASTEP;		// 16
    if (dens_lastmean_p_minz > 0.)
      fdens << "\t" << (dens_mean_p_minz-dens_lastmean_p_minz)/(dt*DENSDATASTEP*DENSDATASTEP);	// 17
    else
      fdens << "\t" << 0.;					// 17

    if (dens_knmeas[i] > 0)					// 18, 19
      fdens << "\t" << dens_kd[i]/(dens_knmeas[i]*dt) << "\t" << dens_ka[i]/(dens_knmeas[i]*dt);
    else
      fdens << "\t" << -0. << "\t" << -0.;

    if (dens_nkd[i] > 0) {					// 20
      fdens << "\t" << dens_sq[i]/dens_nkd[i];
    }
    else							// 20
      fdens << "\t" << -0.0;

    // We output the average flow not velocity here:
    // sum v_i/(#meas*V)=N*<v>/(#meas*V)=<v>*rho=J
 #ifdef RECTANGULAR
    if (dens_nmeas[i] > 0 && vol_i > 0.) {			// 21, 22, 23
      fdens << "\t" << flux_vx[i]/(dens_nmeas[i]*vol_i)
                << "\t" << flux_vy[i]/(dens_nmeas[i]*vol_i)
                << "\t" << flux_vz[i]/(dens_nmeas[i]*vol_i);
    }
    else							// 21, 22, 23
      fdens << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0;
    if (dens_p_nmeas[i] > 0 && vol_i > 0.) {			// 24, 25, 26
      fdens << "\t" << flux_p_vx[i]/(dens_p_nmeas[i]*vol_i)
                << "\t" << flux_p_vy[i]/(dens_p_nmeas[i]*vol_i)
                << "\t" << flux_p_vz[i]/(dens_p_nmeas[i]*vol_i);
    }
    else							// 24, 25, 26
      fdens << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0;
 #elif defined (SPHERICAL)
    if (dens_nmeas[i] > 0 && vol_i > 0.) {			// 21, 22, 23
      fdens << "\t" << 0.0
                << "\t" << 0.0
                << "\t" << flux_vr[i]/(dens_nmeas[i]*vol_i);
    }
    else							// 21, 22, 23
      fdens << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0;
    if (dens_p_nmeas[i] > 0 && vol_i > 0.) {			// 24, 25, 26
      fdens << "\t" << 0.0
                << "\t" << 0.0
                << "\t" << flux_p_vr[i]/(dens_p_nmeas[i]*vol_i);
    }
    else							// 24, 25, 26
      fdens << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0;
 #endif

    if (dens_rhoV[i] > 0) 					// 27
      fdens << "\t" << setprecision(12) << dens_cell_sq[i]/dens_rhoV[i];
    else
      fdens << "\t" << -0.;					// 27
    if (dens_rhoV[i] > 0) 					// 28
      fdens << "\t" << setprecision(12) << dens_mean_cell_dr[i]/dens_rhoV[i];
    else
      fdens << "\t" << -0.;					// 28
    if (dens_rhoV[i] > 0) 					// 29
      fdens << "\t" << setprecision(12) << dens_mean_NN_dr[i]/dens_rhoV[i];
    else
      fdens << "\t" << -0.;					// 29

 #ifdef REALFLUX
    if (area_i > 0.)
      fdens << "\t" << (double)rflux_jz[i]/(area_i*DENSDATASTEP*dt);	// 30
    else
      fdens << "\t" << 0.;	                                // 30
 #endif

    fdens << endl;
  }

  fdens << endl;

  fdens.close();
}
#endif



/////////////////////////////////////
// output_localstress(...)
// -----------------
// output local stress map
// ===
// Parameters:
#ifdef LOCALSTRESSMOP
void output_localstress()
{
  ofstream fstress((string(OUTDIR)+string(FLOCALSTRESSMOPNAME)).c_str(), ios::app);
  int index;
  double area = 0.;

  if (!fstress) {
    cout << "Error while fopen('" << OUTDIR << FLOCALSTRESSMOPNAME << "')" << endl;
    exit(-1);
  }

  for (int i = 0; i < stress_lengthx; ++i)
    for (int j = 0; j < stress_lengthy; ++j)
      for (int k = 0; k < stress_lengthz; ++k) {
        fstress << curtime << "\t" << i*stress_binsize[0] << "\t" << j*stress_binsize[1] << "\t" << k*stress_binsize[2];
        for (int l = 0; l < 3; ++l) {
          index = i*stress_lengthy*stress_lengthz*3+j*stress_lengthz*3+k*3;
          area = 1.;
          for (int q = 0; q < 3; ++q)
            if (q != l)
              area *= stress_binsize[q];
          fstress << "\t" << stress_pcur[index+l]/(area*dt*LSTRESSDATASTEP)
                  << "\t" << stress_fcur[index+l]/(area*LSTRESSDATASTEP)
                  << "\t" << stress_fcurd[index+l]/(area*LSTRESSDATASTEP);
        }
        fstress << endl;
      }
  fstress << endl;
  fstress.close();

  double sum1, sum2, sum3;
  fstress.open((string(OUTDIR)+string("sigmamzmop.dat")).c_str(), ios::app);
  for (int k = 0; k < stress_lengthz; ++k) {
    sum1 = sum2 = sum3 = 0.;
    for (int i = 0; i < stress_lengthx; ++i)
      for (int j = 0; j < stress_lengthy; ++j) {
        index = i*stress_lengthy*stress_lengthz*3+j*stress_lengthz*3+k*3+2;
        sum1 += stress_pcur[index]/dt;
        sum2 += stress_fcur[index];
        sum3 += stress_fcurd[index];
      }
    fstress << curtime << "\t" << k*stress_binsize[2]  << "\t" << sum1/(LX*LY*LSTRESSDATASTEP)
            << "\t" << sum2/(LX*LY*LSTRESSDATASTEP) << "\t" << sum3/(LX*LY*LSTRESSDATASTEP) << endl;
  }
  fstress << endl;
  fstress.close();


  fstress.open((string(OUTDIR)+string("sigmamxmop.dat")).c_str(), ios::app);
  for (int i = 0; i < stress_lengthx-1; ++i) {
    sum1 = sum2 = sum3 = 0.;
    for (int k = 0; k < stress_lengthz; ++k)
      for (int j = 0; j < stress_lengthy; ++j) {
        index = i*stress_lengthy*stress_lengthz*3+j*stress_lengthz*3+k*3+0;
        sum1 += stress_pcur[index]/dt;
        sum2 += stress_fcur[index];
        sum3 += stress_fcurd[index];
      }
    fstress << curtime << "\t" << i*stress_binsize[0]  << "\t" << sum1/(LY*LZ*LSTRESSDATASTEP)
            << "\t" << sum2/(LY*LZ*LSTRESSDATASTEP) << "\t" << sum3/(LY*LZ*LSTRESSDATASTEP) << endl;
  }
  fstress << endl;
  fstress.close();


  fstress.open((string(OUTDIR)+string("sigmamymop.dat")).c_str(), ios::app);
  for (int j = 0; j < stress_lengthy; ++j) {
    sum1 = sum2 = sum3 = 0.;
    for (int i = 0; i < stress_lengthx; ++i)
      for (int k = 0; k < stress_lengthz; ++k) {
        index = i*stress_lengthy*stress_lengthz*3+j*stress_lengthz*3+k*3+1;
        sum1 += stress_pcur[index]/dt;
        sum2 += stress_fcur[index];
        sum3 += stress_fcurd[index];
      }
    fstress << curtime << "\t" << j*stress_binsize[1]  << "\t" << sum1/(LX*LY*LSTRESSDATASTEP)
            << "\t" << sum2/(LX*LY*LSTRESSDATASTEP) << "\t" << sum3/(LX*LY*LSTRESSDATASTEP) << endl;
  }
  fstress << endl;
  fstress.close();
}
#endif


#ifdef LOCALSTRESSVOL
void output_localstress()
{
  ofstream fstress((string(OUTDIR)+string(FLOCALSTRESSVOLNAME)).c_str(), ios::app);
  int index;
  double vol = stress_binsize[0]*stress_binsize[1]*stress_binsize[2];

  if (!fstress) {
    cout << "Error while fopen('" << OUTDIR << FLOCALSTRESSVOLNAME << "')" << endl;
    exit(-1);
  }


  for (int i = 0; i < stress_lengthx; ++i)
    for (int j = 0; j < stress_lengthy; ++j)
      for (int k = 0; k < stress_lengthz; ++k) {
        fstress << curtime << "\t" << i*stress_binsize[0] << "\t" << j*stress_binsize[1] << "\t" << k*stress_binsize[2];
        index = i*stress_lengthy*stress_lengthz*9+j*stress_lengthz*9+k*9;
        for (int l = 0; l < 9; ++l) {
          fstress << "\t" << stress_pcur[index+l]/(vol*dt*LSTRESSDATASTEP)
                  << "\t" << stress_fcur[index+l]/(vol*LSTRESSDATASTEP)
                  << "\t" << stress_fcurd[index+l]/(vol*LSTRESSDATASTEP);
        }
        fstress << endl;
      }
  fstress.close();


  double sum1[9];
  double sum2[9];
  double sum3[9];
  fstress.open((string(OUTDIR)+string("sigmamz.dat")).c_str(), ios::app);
  for (int k = 0; k < stress_lengthz; ++k) {
    fill(sum1, sum1+9, 0.);
    fill(sum2, sum2+9, 0.);
    fill(sum3, sum3+9, 0.);
    for (int i = 0; i < stress_lengthx; ++i)
      for (int j = 0; j < stress_lengthy; ++j) {
        index = i*stress_lengthy*stress_lengthz*9+j*stress_lengthz*9+k*9;
        for (int l = 0; l < 9; ++l) {
          sum1[l] += stress_pcur[index+l];
          sum2[l] += stress_fcur[index+l];
          sum3[l] += stress_fcurd[index+l];
        }
      }
    fstress << curtime << "\t" << k*stress_binsize[2];
    for (int l = 0; l < 9; ++l) {
      fstress << "\t" << sum1[l]/(dt*LX*LY*stress_binsize[2]*LSTRESSDATASTEP)
              << "\t" << sum2[l]/(LX*LY*stress_binsize[2]*LSTRESSDATASTEP)
              << "\t" << sum3[l]/(LX*LY*stress_binsize[2]*LSTRESSDATASTEP);
    }
    fstress << endl;
  }
  fstress.close();

  fstress.open((string(OUTDIR)+string("sigmamx.dat")).c_str(), ios::app);
  for (int k = 0; k < stress_lengthx; ++k) {
    fill(sum1, sum1+9, 0.);
    fill(sum2, sum2+9, 0.);
    fill(sum3, sum3+9, 0.);
    for (int i = 0; i < stress_lengthz; ++i)
      for (int j = 0; j < stress_lengthy; ++j) {
        index = k*stress_lengthy*stress_lengthz*9+j*stress_lengthz*9+i*9;
        for (int l = 0; l < 9; ++l) {
          sum1[l] += stress_pcur[index+l];
          sum2[l] += stress_fcur[index+l];
          sum3[l] += stress_fcurd[index+l];
        }
      }
    fstress << curtime << "\t" << k*stress_binsize[0];
    for (int l = 0; l < 9; ++l) {
      fstress << "\t" << sum1[l]/(dt*LZ*LY*stress_binsize[0]*LSTRESSDATASTEP)
              << "\t" << sum2[l]/(LZ*LY*stress_binsize[0]*LSTRESSDATASTEP)
              << "\t" << sum3[l]/(LZ*LY*stress_binsize[0]*LSTRESSDATASTEP);
    }
    fstress << endl;
  }
  fstress.close();

  fstress.open((string(OUTDIR)+string("sigmamy.dat")).c_str(), ios::app);
  for (int k = 0; k < stress_lengthy; ++k) {
    fill(sum1, sum1+9, 0.);
    fill(sum2, sum2+9, 0.);
    fill(sum3, sum3+9, 0.);
    for (int i = 0; i < stress_lengthz; ++i)
      for (int j = 0; j < stress_lengthx; ++j) {
        index = j*stress_lengthy*stress_lengthz*9+k*stress_lengthz*9+i*9;
        for (int l = 0; l < 9; ++l) {
          sum1[l] += stress_pcur[index+l];
          sum2[l] += stress_fcur[index+l];
          sum3[l] += stress_fcurd[index+l];
        }
      }
    fstress << curtime << "\t" << k*stress_binsize[1];
    for (int l = 0; l < 9; ++l) {
      fstress << "\t" << sum1[l]/(dt*LZ*LX*stress_binsize[1]*LSTRESSDATASTEP)
              << "\t" << sum2[l]/(LZ*LX*stress_binsize[1]*LSTRESSDATASTEP)
              << "\t" << sum3[l]/(LZ*LX*stress_binsize[1]*LSTRESSDATASTEP);
    }
    fstress << endl;
  }
  fstress.close();
}
#endif



/////////////////////////////////////
// output_curp(...)
// -----------------
// output current pressure, box size, etc.
// ===
// Parameters:
#ifdef CONSTP_PBC
void output_curp()
{
  ofstream fcurp((string(OUTDIR)+string(FCONSTPPBCNAME)).c_str(), ios::app);

  fcurp << curtime << "\t" << constp_pbc_curp << "\t"
        << (constp_pbc_pxcur+constp_pbc_pycur+constp_pbc_pzcur)/(3.*CONSTP_PBC_STEP*dt*LX*LY*LZ) << "\t"
        << (constp_pbc_fxcur+constp_pbc_fycur+constp_pbc_fzcur)/(3.*CONSTP_PBC_STEP*LX*LY*LZ) << "\t"
        << (constp_pbc_fxcurd+constp_pbc_fycurd+constp_pbc_fzcurd)/(3.*CONSTP_PBC_STEP*LX*LY*LZ) << "\t"
        << constp_pbc_chi3 << "\t" << constp_pbc_chi << "\t"
        << LX << "\t" << LY << "\t" << LZ << "\t"
        << constp_pbc_pxcur/(CONSTP_PBC_STEP*dt*LX*LY*LZ) << "\t"	// 11
        << constp_pbc_pycur/(CONSTP_PBC_STEP*dt*LX*LY*LZ) << "\t"	// 12
        << constp_pbc_pzcur/(CONSTP_PBC_STEP*dt*LX*LY*LZ) << "\t"	// 13
        << constp_pbc_fxcur/(CONSTP_PBC_STEP*LX*LY*LZ) << "\t"	// 14
        << constp_pbc_fycur/(CONSTP_PBC_STEP*LX*LY*LZ) << "\t"	// 15
        << constp_pbc_fzcur/(CONSTP_PBC_STEP*LX*LY*LZ) << "\t"	// 16
        << constp_pbc_fxcurd/(CONSTP_PBC_STEP*LX*LY*LZ) << "\t"	// 17
        << constp_pbc_fycurd/(CONSTP_PBC_STEP*LX*LY*LZ) << "\t"	// 18
        << constp_pbc_fzcurd/(CONSTP_PBC_STEP*LX*LY*LZ) << "\t"	// 19
        << endl;
  fcurp.close();
}
#endif



/////////////////////////////////////
// dump_all(...)
// -----------------
// dump current state to FDUMPNAME
// ===
// Parameters:
void dump_all()
{
  long i;
  long newlastn;
  long ldummy;
  ofstream fdump((string(OUTDIR)+string(FDUMPNAME)).c_str(), ios::binary|ios::trunc);

  if (!fdump) {
    cerr << "Error while fopen('" << OUTDIR << FDUMPNAME << "')" << endl;
    exit(-1);
  }

  ///////////////////////////////////////
  // layout of dump file:
  //   -> every variable takes 8 bytes (sizeof(double))
  //     * lastn		(long==double)
  //   -> now output ... for every particles
  //     * _r{x,y,z}		(lastn*double)
  //     * _v{x,y,z}		(lastn*double)
  //     * _f{x,y,z}		(lastn*double)
  //     * _fd{x,y,z}		(lastn*double)
  //     * _id			(lastn*unsigned long)
  //     * _spec		(lastn*int)
  //     * all add. parameters like dt,...

  for (i = 0, newlastn = 0; i < lastcn; ++i)
    if (_list[i] == -1)
      newlastn += 2;
  fdump.write((const char *)&newlastn, sizeof(long));
  for (i = 0; i < lastn; ++i)
    if (_list[i/_pperc] == -1) {
      fdump.write((const char *)(_rx+i), sizeof(double));
        fdump.write((const char *)(_ry+i), sizeof(double)); fdump.write((const char *)(_rz+i), sizeof(double));
      fdump.write((const char *)(_vx+i), sizeof(double));
        fdump.write((const char *)(_vy+i), sizeof(double)); fdump.write((const char *)(_vz+i), sizeof(double));
      fdump.write((const char *)(_fx+i), sizeof(double));
        fdump.write((const char *)(_fy+i), sizeof(double)); fdump.write((const char *)(_fz+i), sizeof(double));
      fdump.write((const char *)(_fdx+i), sizeof(double));
        fdump.write((const char *)(_fdy+i), sizeof(double)); fdump.write((const char *)(_fdz+i), sizeof(double));
      fdump.write((const char *)(_id+i), sizeof(unsigned long));
      fdump.write((const char *)(_spec+i), sizeof(int));
    }

  fdump.write((const char *)&next_free_id, sizeof(unsigned long));
  fdump.write((const char *)&curtime, sizeof(double));

  //////////////////////////////////////////
  // everything above will be used when reading dump file, this part is optional
  fdump.write((const char *)&LX, sizeof(double));
  fdump.write((const char *)&LY, sizeof(double));
  fdump.write((const char *)&LZ, sizeof(double));
  fdump.write((const char *)&tend, sizeof(double));
  fdump.write((const char *)&dt, sizeof(double));

  ldummy = MAXSPECNUM;
  fdump.write((const char *)&ldummy, sizeof(long));
  fdump.write((const char *)rpp, ldummy*ldummy*sizeof(double));
  fdump.write((const char *)r0, ldummy*ldummy*sizeof(double));
  fdump.write((const char *)rt, ldummy*ldummy*sizeof(double));
  fdump.write((const char *)rc, ldummy*ldummy*sizeof(double));
  fdump.write((const char *)rct, ldummy*ldummy*sizeof(double));
  fdump.write((const char *)b, ldummy*ldummy*sizeof(double));
  fdump.write((const char *)gammac, ldummy*ldummy*sizeof(double));
  fdump.write((const char *)gammat, ldummy*ldummy*sizeof(double));
  fdump.write((const char *)gammab, ldummy*ldummy*sizeof(double));
  fdump.write((const char *)ka, ldummy*ldummy*sizeof(double));
  fdump.write((const char *)f0, ldummy*ldummy*sizeof(double));
  fdump.write((const char *)f1, ldummy*ldummy*sizeof(double));
  fdump.write((const char *)m, ldummy*ldummy*sizeof(double));

  fdump.write((const char *)&kbt, sizeof(double));
  fdump.write((const char *)&trajdatastep, sizeof(int));
  fdump.write((const char *)&savedatastep, sizeof(int));
  fdump.write((const char *)&epsilon, sizeof(double));

  long delta_trun = time(0)-t_begin;
  fdump.write((const char *)&delta_trun, sizeof(long));

  long bc = 0;
#ifdef PBC
  bc = 1;
#endif
#ifdef BBC
  bc = 2;
#endif
#ifdef RBC
  bc = 3;
#endif
#ifdef HBC
  bc = 4;
#endif
  fdump.write((const char *)&bc, sizeof(long));


  fdump.close();
}



/////////////////////////////////////
// load_dump(...)
// -----------------
// get parameters from data file
// ===
// Parameters:
//   char *fname       : dump file name
//   int load_add_param: switch whether to load additional parameters or not
// ===
// Returns:
//   1: on success
//   0: on failure
int load_dump(char *fname, int load_add_param)
{
  long tmplastn;
  long ldummy = 0;
  double ddummy = 0.;
  long i,k;
  ifstream fdump(fname, ios::binary);


  if (!fdump) {
    _errno = ERRNO_FDUMP_OPEN;
    return 0;
  }

  fdump.read((char *)&tmplastn, sizeof(long));
  if (tmplastn >= INITSIZE) {
    _errno = ERRNO_INITSIZE_TOO_SMALL;
    return 0;
  }

#ifdef STEADYSTATEFIXEDCELLS
  next_free_id = 0;
#endif

  for (k = 0; k < tmplastn/_pperc; ++k) {
    for (i = _pperc*k; i < _pperc*(k+1); ++i) {
      fdump.read((char *)(_rx+i), sizeof(double));
        fdump.read((char *)(_ry+i), sizeof(double)); fdump.read((char *)(_rz+i), sizeof(double));
      fdump.read((char *)(_vx+i), sizeof(double));
        fdump.read((char *)(_vy+i), sizeof(double)); fdump.read((char *)(_vz+i), sizeof(double));
      fdump.read((char *)(_fx+i), sizeof(double));
        fdump.read((char *)(_fy+i), sizeof(double)); fdump.read((char *)(_fz+i), sizeof(double));
      fdump.read((char *)(_fdx+i), sizeof(double));
        fdump.read((char *)(_fdy+i), sizeof(double)); fdump.read((char *)(_fdz+i), sizeof(double));
      fdump.read((char *)(_id+i), sizeof(unsigned long));
      fdump.read((char *)(_spec+i), sizeof(int));
    }
#ifdef FIXEDCELLS
 #ifdef STEADYSTATEFIXEDCELLS
    for (i = _pperc*k; i < _pperc*(k+1); ++i)
      _id[i] = next_free_id++;
    if (_spec[_pperc*k] == 1) {
      ///////////////////////////////
      // displace particle and change species
      for (i = _pperc*k; i < _pperc*(k+1); ++i, ++fixed_wlist_max) {
        fixed_wlist[fixed_wlist_max] = i;
        fixed_wlist_rx0[fixed_wlist_max] = _rx[i];
        fixed_wlist_ry0[fixed_wlist_max] = _ry[i];
        fixed_wlist_rz0[fixed_wlist_max] = _rz[i];
      }
    }
    _head = _list[k];
    _list[k] = -1;
 #else
    for (i = _pperc*k; i < _pperc*(k+1); ++i) {
      if (_rz[i] < fixed_deltaz
  #ifdef DOUBLEFIXEDCELLS
          || _rz[i] > dfixed_lzpos-fixed_deltaz
  #endif
          ) {
        _spec[i] = fixed_spec;
        _spec[__GETOTHERCELLPART(i)] = fixed_spec;
      }
    }


    if (_spec[k*_pperc] != fixed_spec) {
      --k;
      tmplastn -= _pperc;
      continue;
    }
    else {
      ///////////////////////////////
      // displace particle and change species
      for (i = _pperc*k; i < _pperc*(k+1); ++i) {
        if (_rz[i]-0.8 < fixed_deltaz) {		// only needed for DOUBLEFIXEDCELLS
          _rz[i] += fixed_init_displacement;
          fixed_wlist[fixed_wlist_max] = i;
          fixed_wlist_rx0[fixed_wlist_max] = _rx[i];
          fixed_wlist_ry0[fixed_wlist_max] = _ry[i];
          fixed_wlist_rz0[fixed_wlist_max] = _rz[i];
          ++fixed_wlist_max;
          _spec[i] = fixed_spec;
        }
  #ifdef DOUBLEFIXEDCELLS
        else {
          _rz[i] += BOXLENGTH_Z-dfixed_lzpos-2.*fixed_init_displacement;
          fixed_wlist[fixed_wlist_max] = i;
          fixed_wlist_rx0[fixed_wlist_max] = _rx[i];
          fixed_wlist_ry0[fixed_wlist_max] = _ry[i];
          fixed_wlist_rz0[fixed_wlist_max] = _rz[i];
          ++fixed_wlist_max;
          _spec[i] = fixed_spec;
        }
      }
  #else
      }
  #endif // end ifdef DOUBLEFIXEDCELLS
      _head = _list[k];
      _list[k] = -1;
    }

 #endif	// ifdef STEADYSTATEFIXEDCELLS
#else
    _head = _list[k];
    _list[k] = -1;
#endif	// ifdef FIXEDCELLS
  }
  lastcn = k;
  lastn = _pperc*k;

#ifdef STEADYSTATEFIXEDCELLS
  fdump.read((char *)&ldummy, sizeof(unsigned long));
#else
  fdump.read((char *)&next_free_id, sizeof(unsigned long));
#endif
  fdump.read((char *)&curtime, sizeof(double));

  if (load_add_param) {
    fdump.read((char *)&LX, sizeof(double));
    fdump.read((char *)&LY, sizeof(double));
    fdump.read((char *)&LZ, sizeof(double));
    fdump.read((char *)&tend, sizeof(double));
    fdump.read((char *)&dt, sizeof(double));

    boxlength[0] = LX;
    boxlength[1] = LY;
    boxlength[2] = LZ;

    fdump.read((char *)&ldummy, sizeof(long));
    if (ldummy > MAXSPECNUM) {
      _errno = ERRNO_MAXSPECNUM_TOO_SMALL;
      return 0;
    }
    ldummy *= ldummy;
    fdump.read((char *)rpp, ldummy*sizeof(double));
    fdump.read((char *)r0, ldummy*sizeof(double));
    fdump.read((char *)rt, ldummy*sizeof(double));
    fdump.read((char *)rc, ldummy*sizeof(double));
    fdump.read((char *)rct, ldummy*sizeof(double));
    fdump.read((char *)b, ldummy*sizeof(double));
    fdump.read((char *)gammac, ldummy*sizeof(double));
    fdump.read((char *)gammat, ldummy*sizeof(double));
    fdump.read((char *)gammab, ldummy*sizeof(double));
    fdump.read((char *)ka, ldummy*sizeof(double));
    fdump.read((char *)f0, ldummy*sizeof(double));
    fdump.read((char *)f1, ldummy*sizeof(double));
    fdump.read((char *)m, ldummy*sizeof(double));

    fdump.read((char *)&kbt, sizeof(double));
    fdump.read((char *)&trajdatastep, sizeof(int));
    fdump.read((char *)&savedatastep, sizeof(int));
    fdump.read((char *)&epsilon, sizeof(double));
  }
  else {
    fdump.read((char *)&ddummy, sizeof(double));
    fdump.read((char *)&ddummy, sizeof(double));
    fdump.read((char *)&ddummy, sizeof(double));
    fdump.read((char *)&ddummy, sizeof(double));
    fdump.read((char *)&ddummy, sizeof(double));

    fdump.read((char *)&ldummy, sizeof(long));

    ldummy *= ldummy;
    double *_ddummy = new double[ldummy];
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    fdump.read((char *)_ddummy, ldummy*sizeof(double));
    delete[] _ddummy;

    fdump.read((char *)&ddummy, sizeof(double));
    fdump.read((char *)&ldummy, sizeof(int));
    fdump.read((char *)&ldummy, sizeof(int));
    fdump.read((char *)&ddummy, sizeof(double));
  }

  fdump.read((char *)&t_begin, sizeof(long));

  long bc = 0;
  fdump.read((char *)&bc, sizeof(long));

#ifndef PGAS		// dirty work around!!!!
#ifdef PBC
  if (bc != 1) {
    _errno = ERRNO_DIFFERING_BC;
    return 0;
  }
#endif
#ifdef BBC
  if (bc != 2) {
    _errno = ERRNO_DIFFERING_BC;
    return 0;
  }
#endif
#ifdef RBC
  if (bc != 3) {
    _errno = ERRNO_DIFFERING_BC;
    return 0;
  }
#endif
#ifdef HBC
  if (bc != 4) {
    _errno = ERRNO_DIFFERING_BC;
    return 0;
  }
#endif
#endif



#ifdef FIXEDCELLS
    curtime = 0.;

 #ifndef STEADYSTATEFIXEDCELLS
    /////////////////////////////////
    // adjust cell ids
    _id[lastn] = next_free_id++;
    _id[lastn+1] = next_free_id++;

    /////////////////////////////////
    // adjust 'linked' list
    _vx[lastn] = _vy[lastn] = _vz[lastn] = _fx[lastn] = _fy[lastn] = _fz[lastn] = 0.;
    _fdx[lastn] = _fdy[lastn] = _fdz[lastn] = 0.;
    _vx[lastn+1] = _vy[lastn+1] = _vz[lastn+1] = _fx[lastn+1] = _fy[lastn+1] = _fz[lastn+1] = 0.;
    _fdx[lastn+1] = _fdy[lastn+1] = _fdz[lastn+1] = 0.;

    _rx[lastn] = LX/2.;
    _ry[lastn] = LY/2.;
    _rz[lastn] = fixed_init_displacement+fixed_deltaz+0.4;
    double initdist = 0.4;
    double tmpx = (drand[0]()*initdist*2.-initdist)/sqrt(3.);
    double tmpy = (drand[0]()*initdist*2.-initdist)/sqrt(3.);
    double tmpz = sqrt(initdist*initdist-tmpx*tmpx-tmpy*tmpy);
    _rx[lastn+1] = LX/2. + tmpx;
    _ry[lastn+1] = LY/2. + tmpy;
    _rz[lastn+1] = fixed_init_displacement+fixed_deltaz+0.4 + tmpz;

    _head = _list[lastcn];
    _list[lastcn] = -1;
    ++lastcn;
    lastn = _pperc*lastcn;
 #endif	// ifndef STEADYSTATEFIXEDCELLS
#endif	// ifdef FIXEDCELLS


#ifdef LOCALSTRESS
  for (i = 0; i < lastn; ++i) {
    stress_rx[i] = _rx[i];
    stress_ry[i] = _ry[i];
    stress_rz[i] = _rz[i];
  }
#endif


#ifdef PGAS
  curtime = 0.;
  double tx, ty, tz;
  i = lastcn;
  int j;
  const double dxyz = PGAS_DX;
  const double drx = 0.;

 #ifdef SPHERICAL
  for (int k = 0; k < lastcn; ++k)
    if (_list[k] == -1
  #ifdef FIXEDCELLS
        && _spec[k*_pperc] != fixed_spec
  #endif
  #ifdef PGAS
        && _spec[k*_pperc] != gas_spec
  #endif
        ) {
      _rx[k*_pperc] -= drx;
      _rz[k*_pperc] -= drx;
      _ry[k*_pperc] -= drx;
      _rx[k*_pperc+1] -= drx;
      _rz[k*_pperc+1] -= drx;
      _ry[k*_pperc+1] -= drx;
    }
  for (tx = 0.1; tx < LX; tx += dxyz) {
    for (ty = 0.1; ty < LY; ty += dxyz) {
      for (tz = 0.1; tz < LZ; tz += dxyz) {
        if (NORM(tx-LX/2.,ty-LY/2.,tz-LZ/2.) < PGAS_DR)
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
  lastn = _pperc*lastcn;
 #elif defined RECTANGULAR
  const double LZl = 17.;
  const double LZu = 63.;
  for (tx = dxyz; tx < LX-dxyz; tx += dxyz) {
    for (ty = dxyz; ty < LY-dxyz; ty += dxyz) {
      for (tz = dxyz; tz < LZl; tz += dxyz) {
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
      for (tz = LZu; tz < LZ-dxyz; tz += dxyz) {
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
  lastn = _pperc*lastcn;
 #endif // ifdef SPH elif REC
#endif	// ifdef PGAS




  fdump.close();

  return 1;
}







