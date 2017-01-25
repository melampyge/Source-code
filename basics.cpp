/*		basics.cpp		*/
//////////////////////////////////////////
// basic types and function 		//
// definitions				//
//////////////////////////////////////////

#include "basics.h"


////////////////////////////////////
// DWORD crc32buf(...)
// ---------------
// calculate crc32 check sum of buf
// ===
// Parameters:
//	char *buf		: buffer to compute crc32 from
//	unsigned long len	: length of buffer
DWORD crc32buf(char *buf, unsigned long len)
{
      DWORD crc32;

      crc32 = 0xFFFFFFFF;

      for ( ; len; --len, ++buf)
      {
            crc32 = crc_32_tab[((crc32) ^ ((BYTE)*buf)) & 0xFF] ^ ((crc32) >> 8);
      }

      return ~crc32;
}



////////////////////////////////////
// void allocate_init_arrays(...)
// ---------------
// allocate and initialize all global arrays
// ===
// Parameters:
int allocate_init_arrays()
{
  long i, j;

  curtime = 0;

  try {
    _rx   = new double[INITSIZE];
    _ry   = new double[INITSIZE];
    _rz   = new double[INITSIZE];
    _vx   = new double[INITSIZE];
    _vy   = new double[INITSIZE];
    _vz   = new double[INITSIZE];
    _fx   = new double[INITSIZE];
    _fy   = new double[INITSIZE];
    _fz   = new double[INITSIZE];
    _fdx  = new double[INITSIZE];
    _fdy  = new double[INITSIZE];
    _fdz  = new double[INITSIZE];
    _id   = new long[INITSIZE];
    _spec = new int[INITSIZE];
    vlist = new long[VLISTINITSIZE];	// make vlist larger by mean neighbour number
    vlist_offset = new long[INITSIZE+1];
    vlist_nc = new long[INITSIZE];
    vlist_dx = new double[INITSIZE];
    vlist_dz = new double[INITSIZE];
    vlist_dy = new double[INITSIZE];
    rna      = new double[VLISTINITSIZE];
    rnai     = new unsigned char[VLISTINITSIZE];

    std::fill(_rx, _rx+INITSIZE, 0.); std::fill(_ry, _ry+INITSIZE, 0.); std::fill(_rz, _rz+INITSIZE, 0.);
    std::fill(_vx, _vx+INITSIZE, 0.); std::fill(_vy, _vy+INITSIZE, 0.); std::fill(_vz, _vz+INITSIZE, 0.);
    std::fill(_fx, _fx+INITSIZE, 0.); std::fill(_fy, _fy+INITSIZE, 0.); std::fill(_fz, _fz+INITSIZE, 0.);
    std::fill(_fdx, _fdx+INITSIZE, 0.); std::fill(_fdy, _fdy+INITSIZE, 0.); std::fill(_fdz, _fdz+INITSIZE, 0.);
    std::fill(vlist_dx, vlist_dx+INITSIZE, 0.); std::fill(vlist_dy, vlist_dy+INITSIZE, 0.);
      std::fill(vlist_dz, vlist_dz+INITSIZE, 0.);
    std::fill(_id, _id+INITSIZE, 0);
    std::fill(_spec, _spec+INITSIZE, 0);
    std::fill(vlist, vlist+VLISTINITSIZE, 0);
    std::fill(vlist_nc, vlist_nc+INITSIZE, 0);
    std::fill(vlist_offset, vlist_offset+INITSIZE+1, 0);
    memset(rnai, 0x00, sizeof(unsigned char)*VLISTINITSIZE);

#ifdef PMEASUREMENT
    pm_px0 = new double[TARGETNUMTHREADS];
    pm_pxL = new double[TARGETNUMTHREADS];
    pm_py0 = new double[TARGETNUMTHREADS];
    pm_pyL = new double[TARGETNUMTHREADS];
    pm_pz0 = new double[TARGETNUMTHREADS];
    pm_pzL = new double[TARGETNUMTHREADS];

    std::fill(pm_px0, pm_px0+TARGETNUMTHREADS, 0.);
    std::fill(pm_pxL, pm_pxL+TARGETNUMTHREADS, 0.);
    std::fill(pm_py0, pm_py0+TARGETNUMTHREADS, 0.);
    std::fill(pm_pyL, pm_pyL+TARGETNUMTHREADS, 0.);
    std::fill(pm_pz0, pm_pz0+TARGETNUMTHREADS, 0.);
    std::fill(pm_pzL, pm_pzL+TARGETNUMTHREADS, 0.);
#endif

#ifdef DENSITYPROFILE
    dens_rhoV    = new int[dens_length];
    dens_rhoVsq  = new int[dens_length];
    dens_rhoVsqtmp = new int[dens_length];
    dens_nmeas   = new int[dens_length];
    dens_p_nmeas = new int[dens_length];
    dens_knmeas  = new int[dens_length];
    dens_nkd     = new int[dens_length];
    dens_nka     = new int[dens_length];
    dens_ka      = new double[dens_length];
    dens_kd      = new double[dens_length];
    dens_kperp   = new double[dens_length];
    dens_kpara   = new double[dens_length];
    dens_curn    = new int[dens_length];
    dens_curnka  = new int[dens_length];
    dens_curnkd  = new int[dens_length];
    dens_sq      = new double[dens_length];
    dens_cell_sq = new double[dens_length];
    dens_mean_cell_dr = new double[dens_length];
    dens_mean_NN_dr   = new double[dens_length];

    std::fill(dens_rhoV, dens_rhoV+dens_length, 0);
    std::fill(dens_rhoVsq, dens_rhoVsq+dens_length, 0);
    std::fill(dens_rhoVsqtmp, dens_rhoVsqtmp+dens_length, 0);
    std::fill(dens_nmeas, dens_nmeas+dens_length, 0);
    std::fill(dens_knmeas, dens_knmeas+dens_length, 0);
    std::fill(dens_nkd, dens_nkd+dens_length, 0);
    std::fill(dens_nka, dens_nka+dens_length, 0);
    std::fill(dens_p_nmeas, dens_p_nmeas+dens_length, 0);
    std::fill(dens_ka, dens_ka+dens_length, 0.);
    std::fill(dens_kd, dens_kd+dens_length, 0.);
    std::fill(dens_kpara, dens_kpara+dens_length, 0.);
    std::fill(dens_kperp, dens_kperp+dens_length, 0.);
    std::fill(dens_curn, dens_curn+dens_length, 0);
    std::fill(dens_curnka, dens_curnka+dens_length, 0);
    std::fill(dens_curnkd, dens_curnkd+dens_length, 0);
    std::fill(dens_sq, dens_sq+dens_length, 0.);
    std::fill(dens_cell_sq, dens_cell_sq+dens_length, 0.);
    std::fill(dens_mean_cell_dr, dens_mean_cell_dr+dens_length, 0.);
    std::fill(dens_mean_NN_dr, dens_mean_NN_dr+dens_length, 0.);
 #ifdef RECTANGULAR
    flux_vx      = new double[dens_length];
    flux_vy      = new double[dens_length];
    flux_vz      = new double[dens_length];
    flux_p_vx    = new double[dens_length];
    flux_p_vy    = new double[dens_length];
    flux_p_vz    = new double[dens_length];

    std::fill(flux_vx, flux_vx+dens_length, 0.);
    std::fill(flux_vy, flux_vy+dens_length, 0.);
    std::fill(flux_vz, flux_vz+dens_length, 0.);
    std::fill(flux_p_vx, flux_p_vx+dens_length, 0.);
    std::fill(flux_p_vy, flux_p_vy+dens_length, 0.);
    std::fill(flux_p_vz, flux_p_vz+dens_length, 0.);
 #elif defined (SPHERICAL)
    flux_vr      = new double[dens_length];
    flux_p_vr    = new double[dens_length];

    std::fill(flux_vr, flux_vr+dens_length, 0.);
    std::fill(flux_p_vr, flux_p_vr+dens_length, 0.);
 #endif
 #ifdef REALFLUX
    rflux_jz     = new long[dens_length];
    rflux_crz    = new double[(int)_ceil(INITSIZE/_pperc)];

    std::fill(rflux_jz, rflux_jz+dens_length, 0);
    std::fill(rflux_crz, rflux_crz+(int)_ceil(INITSIZE/_pperc), 0.);
 #endif
#endif  // ifdef DENSITYPROFILE

#ifdef LOCALSTRESS
    stress_fcur   = new double[stress_length];
    stress_fcurd  = new double[stress_length];
    stress_pcur   = new double[stress_length];
    stress_fnext  = new double *[TARGETNUMTHREADS];
    stress_fnextd = new double *[TARGETNUMTHREADS];
    stress_pnext  = new double *[TARGETNUMTHREADS];
    stress_rx     = new double[INITSIZE];
    stress_ry     = new double[INITSIZE];
    stress_rz     = new double[INITSIZE];
    for (i = 0; i < TARGETNUMTHREADS; ++i) {
      stress_fnext[i]  = new double[stress_length];
      stress_fnextd[i] = new double[stress_length];
      stress_pnext[i]  = new double[stress_length];

      std::fill(stress_fnext[i], stress_fnext[i]+stress_length, 0.);
      std::fill(stress_fnextd[i], stress_fnextd[i]+stress_length, 0.);
      std::fill(stress_pnext[i], stress_pnext[i]+stress_length, 0.);
    }

    std::fill(stress_fcur, stress_fcur+stress_length, 0.);
    std::fill(stress_fcurd, stress_fcurd+stress_length, 0.);
    std::fill(stress_pcur, stress_pcur+stress_length, 0.);
    std::fill(stress_rx, stress_rx+INITSIZE, 0.);
    std::fill(stress_ry, stress_ry+INITSIZE, 0.);
    std::fill(stress_rz, stress_rz+INITSIZE, 0.);
#endif

#ifdef FIXEDCELLS
    fixed_wlist     = new long[WLISTSIZE];
    fixed_wlist_rx0 = new double[WLISTSIZE];
    fixed_wlist_ry0 = new double[WLISTSIZE];
    fixed_wlist_rz0 = new double[WLISTSIZE];
    fixed_pkz0      = new double[TARGETNUMTHREADS];
    fixed_pkzL      = new double[TARGETNUMTHREADS];

    std::fill(fixed_wlist, fixed_wlist+WLISTSIZE, 0);
    std::fill(fixed_wlist_rx0, fixed_wlist_rx0+WLISTSIZE, 0.);
    std::fill(fixed_wlist_ry0, fixed_wlist_ry0+WLISTSIZE, 0.);
    std::fill(fixed_wlist_rz0, fixed_wlist_rz0+WLISTSIZE, 0.);
    std::fill(fixed_pkz0, fixed_pkz0+TARGETNUMTHREADS, 0.);
    std::fill(fixed_pkzL, fixed_pkzL+TARGETNUMTHREADS, 0.);
#endif // ifdef FIXEDCELLS

#ifdef CONSTP_PBC
    constp_pbc_fxnext  = new double[TARGETNUMTHREADS];
    constp_pbc_fxnextd = new double[TARGETNUMTHREADS];
    constp_pbc_fynext  = new double[TARGETNUMTHREADS];
    constp_pbc_fynextd = new double[TARGETNUMTHREADS];
    constp_pbc_fznext  = new double[TARGETNUMTHREADS];
    constp_pbc_fznextd = new double[TARGETNUMTHREADS];

    std::fill(constp_pbc_fxnext, constp_pbc_fxnext+TARGETNUMTHREADS, 0.);
    std::fill(constp_pbc_fxnextd, constp_pbc_fxnextd+TARGETNUMTHREADS, 0.);
    std::fill(constp_pbc_fynext, constp_pbc_fynext+TARGETNUMTHREADS, 0.);
    std::fill(constp_pbc_fynextd, constp_pbc_fynextd+TARGETNUMTHREADS, 0.);
    std::fill(constp_pbc_fznext, constp_pbc_fznext+TARGETNUMTHREADS, 0.);
    std::fill(constp_pbc_fznextd, constp_pbc_fznextd+TARGETNUMTHREADS, 0.);
#endif
  }
  catch (std::bad_alloc &memAllocExcp) {
    _errno = ERRNO_MEM_ALLOC;
    _errno_add_info = memAllocExcp.what();
    return 0;
  }
  catch(...) {
    _errno = ERRNO_UNKNOWN;
    return 0;
  }


  ///////////////////////////////////
  // initialize cell list
  _list = new long[(long)_ceil(INITSIZE/_pperc)];
  for (i = 0; i < (long)_ceil(INITSIZE/_pperc); ++i)
    _list[i] = i+1;
  _list[(long)_ceil(INITSIZE/_pperc)-1] = -1;
  _head = 0;
  lastn = 0;
  lastcn = 0;
  maxn = INITSIZE;


  try {
    rpp    = new double[MAXSPECNUM*MAXSPECNUM];
    r0     = new double[MAXSPECNUM*MAXSPECNUM];
    rt     = new double[MAXSPECNUM*MAXSPECNUM];
    rc     = new double[MAXSPECNUM*MAXSPECNUM];
    rct    = new double[MAXSPECNUM*MAXSPECNUM];
    b      = new double[MAXSPECNUM*MAXSPECNUM];
    gammac = new double[MAXSPECNUM*MAXSPECNUM];
    gammat = new double[MAXSPECNUM*MAXSPECNUM];
    gammab = new double[MAXSPECNUM*MAXSPECNUM];
    ka     = new double[MAXSPECNUM*MAXSPECNUM];
    f0     = new double[MAXSPECNUM*MAXSPECNUM];
    f1     = new double[MAXSPECNUM*MAXSPECNUM];
    m      = new double[MAXSPECNUM*MAXSPECNUM];
#ifdef SUBSTRATE
    gammas = new double[MAXSPECNUM*MAXSPECNUM];
    rz	   = new double[MAXSPECNUM*MAXSPECNUM];
    fa	   = new double[MAXSPECNUM*MAXSPECNUM];
#endif
    

    randf_prefc = new double[MAXSPECNUM*MAXSPECNUM];
    randf_preft = new double[MAXSPECNUM*MAXSPECNUM];
  }
  catch (std::bad_alloc &memAllocExcp) {
    _errno = ERRNO_MEM_ALLOC;
    _errno_add_info = memAllocExcp.what();
    return 0;
  }
  catch(...) {
    _errno = ERRNO_UNKNOWN;
    return 0;
  }

  ///////////////////////////////////
  // initialize species properties
  for (i = 0; i < MAXSPECNUM; ++i) {
    for (j = 0; j < MAXSPECNUM; ++j) {
      if (i==j) {
        rpp[i*MAXSPECNUM+j] = RPP;
        r0[i*MAXSPECNUM+j] = R0;
        rt[i*MAXSPECNUM+j] = RT;
        rc[i*MAXSPECNUM+j] = RC;
        rct[i*MAXSPECNUM+j] = RCT;
        b[i*MAXSPECNUM+j] = B;
        gammac[i*MAXSPECNUM+j] = GAMMAC;
        gammat[i*MAXSPECNUM+j] = GAMMAT;
        gammab[i*MAXSPECNUM+j] = GAMMAB;
        ka[i*MAXSPECNUM+j] = KA;
        f0[i*MAXSPECNUM+j] = F0;
        f1[i*MAXSPECNUM+j] = F1;
        m[i*MAXSPECNUM+j] = MASS;
#ifdef SUBSTRATE
	fa[i*MAXSPECNUM+j] = FA;
	gammas[i*MAXSPECNUM+j] = GAMMAS;
	rz[i*MAXSPECNUM+j] = RZ;
#endif
      }
      else {
        rpp[i*MAXSPECNUM+j] = RPP;
        r0[i*MAXSPECNUM+j] = R0;
        rt[i*MAXSPECNUM+j] = RT;
        rc[i*MAXSPECNUM+j] = RC;
        rct[i*MAXSPECNUM+j] = RCT;
        b[i*MAXSPECNUM+j] = B;
        gammac[i*MAXSPECNUM+j] = GAMMAC;
        gammat[i*MAXSPECNUM+j] = GAMMAT;
        gammab[i*MAXSPECNUM+j] = GAMMAB;;
        ka[i*MAXSPECNUM+j] = KA;
        f0[i*MAXSPECNUM+j] = F0;
        f1[i*MAXSPECNUM+j] = F1;
        m[i*MAXSPECNUM+j] = MASS;
#ifdef SUBSTRATE
	fa[i*MAXSPECNUM+j] = FA;
	gammas[i*MAXSPECNUM+j] = GAMMAS;
	rz[i*MAXSPECNUM+j] = RZ;
#endif
      }
      randf_prefc[i*MAXSPECNUM+j] = sqrt(2.*gammac[i*MAXSPECNUM+j]*kbt/dt);
      randf_preft[i*MAXSPECNUM+j] = sqrt(2.*gammat[i*MAXSPECNUM+j]*kbt/dt);
    }
  }

#ifdef FIXEDCELLS
  ka[fixed_spec]  = 0.;
  rct[fixed_spec] = std::numeric_limits<double>::infinity();
#endif

#ifdef PGAS
  ka[gas_spec]  = 0.;
  rct[gas_spec] = std::numeric_limits<double>::infinity();
  for (i = 0; i <= gas_spec; ++i) {
    f1[i*MAXSPECNUM+gas_spec] = f1[gas_spec*MAXSPECNUM+i] = 0.;
    gammac[i*MAXSPECNUM+gas_spec] = gammac[gas_spec*MAXSPECNUM+i] = 0.;
    gammat[i*MAXSPECNUM+gas_spec] = gammat[gas_spec*MAXSPECNUM+i] = 0.;
#ifdef SUBSTRATE
    gammas[i*MAXSPECNUM+gas_spec] = gammas[gas_spec*MAXSPECNUM+i] = 0.;
#endif
    randf_prefc[i*MAXSPECNUM+gas_spec] = randf_prefc[gas_spec*MAXSPECNUM+i] = 0.;
    randf_preft[i*MAXSPECNUM+gas_spec] = randf_preft[gas_spec*MAXSPECNUM+i] = 0.;
  }
  f0[gas_spec*MAXSPECNUM+gas_spec] = PGAS_F0;
  b[gas_spec] = 0.;
  gammab[gas_spec] = 0.;
#endif

  ///////////////////////////////////
  // get maximum interaction range
  max_ia_range = 0.;
  for (i = 0; i < MAXSPECNUM*MAXSPECNUM; ++i) {
    if (rpp[i] > max_ia_range)
      max_ia_range = rpp[i];
    if (rt[i] > max_ia_range)
      max_ia_range = rt[i];
  }
  currc = max_ia_range+skin;

  ///////////////////////////////////
  // seed rng with get_seed() or predefined value (see parameter.h)
  unsigned long trngseed = RNGSEED;
  if (_errno != 0)
    return 0;
  std::cout << "RNGSEED: " << trngseed << std::endl;
  trngseed = crc32buf((char *)&trngseed, sizeof(trngseed));	// crc32buf() calcs crc32 checksum of seed
  MTRand_int32 tmprand(trngseed);				// result for similar input is very different
  try { drand = new MTRand[MAXTHREADNUM]; }
  catch (std::bad_alloc &memAllocExcp) {
    _errno = ERRNO_MEM_ALLOC;
    _errno_add_info = memAllocExcp.what();
    return 0;
  }
  catch(...) {
    _errno = ERRNO_UNKNOWN;
    return 0;
  }

  for (i = 0; i < MAXTHREADNUM; ++i) {
    trngseed = tmprand();
    trngseed = crc32buf((char *)&trngseed, sizeof(trngseed));
    drand[i].seed(trngseed);
  }

  return 1;
}



////////////////////////////////////
// double get_ran_gaussian(...)
// ---------------
// Get gaussian distributed rn from uniform rn
// using the "Marsaglia polar method"
// We are using separate functions since for the
// DPD algorithm we have mu=0 and sigma=1
// ===
// Parameters:
//	int tid		: thread id
//	double mu	: mean
//	double sigma	: standard deviation
// ===
// Return:
//	double		: gaussian distributed rn
double get_ran_gaussian(int tid)
{
  double a1, a2;
  static double savednum[TARGETNUMTHREADS];
  static int saved[TARGETNUMTHREADS];		// luckily static variables are default initialzed to 0
  double p, q;

  if (saved[tid] == 0) {
    do {
      a1 = 2.*drand[tid]()-1.;
      a2 = 2.*drand[tid]()-1.;
      q = a1*a1 + a2*a2;
    } while (q >= 1. || q == 0.);
    p = sqrt(-2.*log(q)/q);
    savednum[tid] = a2*p;
    saved[tid] = 1;

    return a1*p;
  }
  else {
    saved[tid] = 0;
    return savednum[tid];
  }
}

double get_ran_gaussian(int tid, double mu, double sigma)
{
  double a1, a2;
  static double savednum[TARGETNUMTHREADS];
  static int saved[TARGETNUMTHREADS];		// luckily static variables are default initialzed to 0
  double p, q;

  if (saved[tid] == 0) {
    do {
      a1 = 2.*drand[tid]()-1.;
      a2 = 2.*drand[tid]()-1.;
      q = a1*a1 + a2*a2;
    } while (q >= 1. || q == 0.);
    p = sqrt(-2.*log(q)/q);
    savednum[tid] = a2*p;
    saved[tid] = 1;

    return mu - sigma*a1*p;
  }
  else {
    saved[tid] = 0;
    return mu - sigma*savednum[tid];
  }
}


////////////////////////////////////
// unsigned long get_seed()
// ---------------
// get good seed (from /dev/urandom), in case simulations are started at exact same time
// this still ensures different initial seeds
// ===
// Return:
//	unsigned long: seed
unsigned long get_seed()
{
  std::ifstream fin;
  unsigned long seed = 0;

  fin.open("/dev/urandom", std::ios::binary);
  if (!fin) {
    _errno = ERRNO_URANDOM;
    return seed;	// if error, return 0
  }

  fin.read((char *)&seed, sizeof(unsigned long));
  fin.close();

  return seed;
}



