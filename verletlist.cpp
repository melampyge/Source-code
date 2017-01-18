/*		verletlist.cpp		*/
//////////////////////////////////////////
// functions to create and maintain	//
// the verlet list			//
//////////////////////////////////////////


#include "verletlist.h"


////////////////////////////////////
// void check_verlet_list(...)
// ---------------
// check whether verlet list needs updating by calculating two maximum
// displacements since last update
// => actually the squares of all operands are calculated!
// ===
// Parameters:
void check_verlet_list()
{
  double vlist_maxdr1 = 0.;		// maximum displacement
  double vlist_maxdr2 = 0.;		// second maximum displacement
  double drsq;
  long i, k;

  for (k = 0; k < lastcn; ++k) {
    if (_list[k] != -1)
      continue;
    for (i = _pperc*k; i < _pperc*(k+1); ++i) {	// loop over k
      drsq = vlist_dx[i]*vlist_dx[i] + vlist_dy[i]*vlist_dy[i] + vlist_dz[i]*vlist_dz[i];
      if (drsq > vlist_maxdr1) {
        vlist_maxdr2 = vlist_maxdr1;
        vlist_maxdr1 = drsq;
      }
      else if (drsq > vlist_maxdr2)
        vlist_maxdr2 = drsq;
#ifdef PGAS
      if (_spec[k*_pperc] == gas_spec)
        break;
#endif
    } // end loop over k
  } // end for k < lastcn

#ifdef __DEBUG_VLIST
  //std::cout << "__DEBUG at t=" << curtime << ": maxdr1=" << vlist_maxdr1 << "; maxdr2=" << vlist_maxdr2 << std::endl;
#endif

  //////////////////////////////////////
  // if sum of two maximum displacements are bigger than skin
  //   -> update verlet list
  if ((sqrt(vlist_maxdr1) + sqrt(vlist_maxdr2)) > (currc-max_ia_range)) {
#ifdef __DEBUG_VLIST
    std::cout << "__DEBUG at t=" << curtime << ": update_verlet_list()" << std::endl;
#endif
    update_verlet_list();
  }
}



void update_verlet_list() {
  long i, k;
  long **l_vlist = new long*[TARGETNUMTHREADS];
  static int *list,first=1;
  static int *cx,*cy,*cz;
  int nthreads;

  START_TIMER(9)
  START_TIMER(13)

  if (first) {
    cx   = new int[INITSIZE];
    cy   = new int[INITSIZE];
    cz   = new int[INITSIZE];
    list = new int[INITSIZE];
    first=0;
  }
  double rc = max_ia_range+skin;
  double rc_[3] = {LX/_floor(LX/rc), LY/_floor(LY/rc), LZ/_floor(LZ/rc)};
  int cdx,cdy,cdz;
#ifdef PBC
  cdx=2+(int)_ceil(LX/rc_[0]);
  cdy=2+(int)_ceil(LY/rc_[1]);
  cdz=2+(int)_ceil(LZ/rc_[2]);
#elif defined (RBC) || defined (BBC)
  double minx=LX,miny=LY,minz=LZ;
  double maxx=0,maxy=0,maxz=0;
  for(k=0;k<lastcn;++k) {
    if (_list[k] != -1)
      continue;
    for (i = _pperc*k; i < _pperc*(k+1); ++i) {	// loop over k
      if (_rx[i]<minx)
        minx=_rx[i];
      if (_ry[i]<miny)
        miny=_ry[i];
      if (_rz[i]<minz)
        minz=_rz[i];
      if (_rx[i]>maxx)
        maxx=_rx[i];
      if (_ry[i]>maxy)
        maxy=_ry[i];
      if (_rz[i]>maxz)
        maxz=_rz[i];
 #ifdef PGAS
      if (_spec[i] == gas_spec)		// if gas particle -> don't check >k*_pperc!!!
        break;
 #endif
    } // end loop over k
  } // end loop over cells
  double _dx[3] = {maxx-minx, maxy-miny, maxz-minz};
  for (i = 0; i < 3; ++i) {
    rc_[i] = _dx[i]/_floor(_dx[i]/rc);
    if (!std::isfinite(rc_[i]))		// if we only have one box smaller than rc...
      rc_[i] = rc;
  }
  // In principle there is nothing wrong with rc_==infinity but the following
  // would be wrong (no box for particles...
  cdx=2+(int)_ceil(_dx[0]/rc_[0]);
  cdy=2+(int)_ceil(_dx[1]/rc_[1]);
  cdz=2+(int)_ceil(_dx[2]/rc_[2]);
#elif defined (HBC)
  double minz=LZ;
  double maxz=0;
  for(k=0;k<lastcn;++k) {
    if (_list[k] != -1)
      continue;
    for (i = _pperc*k; i < _pperc*(k+1); ++i) {	// loop over k
      if (_rz[i]<minz)
        minz=_rz[i];
      if (_rz[i]>maxz)
        maxz=_rz[i];
 #ifdef PGAS
      if (_spec[i] == gas_spec)		// if gas particle -> don't check >k*_pperc!!!
        break;
 #endif
    } // end loop over k
  } // end loop over cells
  double _dz = maxz-minz;
  rc_[2] = _dz/_floor(_dz/rc);
  if (!std::isfinite(rc_[2]))		// if we only have one box smaller than rc...
    rc_[2] = rc;
  cdx=2+(int)_ceil(LX/rc_[0]);
  cdy=2+(int)_ceil(LY/rc_[1]);
  cdz=2+(int)_ceil(_dz/rc_[2]);
#endif
  int *head,ncell,cell_;
  int ix,iy,iz,ic,nn;
  ncell=cdx*cdy*cdz;
  head = new int[ncell];
  for(i=0;i<ncell;++i)
    head[i]=-1;
  for(k=0;k<lastcn;++k) {
    if (_list[k] != -1)
      continue;
    for (i = _pperc*k; i < _pperc*(k+1); ++i) {	// loop over k
#ifdef PBC
      cx[i] = 1 + (int)((_rx[i] - LX*_floor(_rx[i]/LX))/rc_[0]);
      cy[i] = 1 + (int)((_ry[i] - LY*_floor(_ry[i]/LY))/rc_[1]);
      cz[i] = 1 + (int)((_rz[i] - LZ*_floor(_rz[i]/LZ))/rc_[2]);
#elif defined (RBC) || defined (BBC)
      cx[i] = 1 + (int)((_rx[i]-minx)/rc_[0]);
      if (_rx[i] == maxx && cx[i] > 1)	// special case
	cx[i]--;
      cy[i] = 1 + (int)((_ry[i]-miny)/rc_[1]);
      if (_ry[i] == maxy && cy[i] > 1)	// special case
	cy[i]--;
      cz[i] = 1 + (int)((_rz[i]-minz)/rc_[2]);
      if (_rz[i] == maxz && cz[i] > 1)	// special case
	cz[i]--;
#elif defined (HBC)
      cx[i] = 1 + (int)((_rx[i] - LX*_floor(_rx[i]/LX))/rc_[0]);
      cy[i] = 1 + (int)((_ry[i] - LY*_floor(_ry[i]/LY))/rc_[1]);
      cz[i] = 1 + (int)((_rz[i]-minz)/rc_[2]);
      if (_rz[i] == maxz && cz[i] > 1)	// special case
	cz[i]--;
#endif
      cell_=cz[i]*cdx*cdy+cy[i]*cdx+cx[i];
      list[i]=head[cell_];
      head[cell_]=i;
#ifdef PGAS
      if (_spec[i] == gas_spec)		// if gas particle -> don't check >k*_pperc!!!
        break;
#endif
    } // end loop over k
  } // end loop over cells
#ifdef PBC
  for(iz=1;iz<cdz-1;++iz)
    for(iy=1;iy<cdy-1;++iy) {
      cell_=iz*cdx*cdy+iy*cdx;
      head[cell_]=head[cell_+cdx-2];
      head[cell_+cdx-1]=head[cell_+1];
      }
  for(iy=1;iy<cdy-1;++iy)
    for(ix=0;ix<cdx;++ix) {
      cell_=iy*cdx+ix;
      head[cell_]=head[cell_+cdx*cdy*(cdz-2)];
      head[cell_+cdx*cdy*(cdz-1)]=head[cell_+cdx*cdy];
      }
  for(iz=0;iz<cdz;++iz)
    for(ix=0;ix<cdx;++ix) {
      cell_=iz*cdx*cdy+ix;
      head[cell_]=head[cell_+cdx*(cdy-2)];
      head[cell_+cdx*(cdy-1)]=head[cell_+cdx];
      }
#elif defined (HBC)
  for(iz=1;iz<cdz-1;++iz)
    for(iy=1;iy<cdy-1;++iy) {
      cell_=iz*cdx*cdy+iy*cdx;
      head[cell_]=head[cell_+cdx-2];
      head[cell_+cdx-1]=head[cell_+1];
      }
  for(iz=1;iz<cdz-1;++iz)
    for(ix=0;ix<cdx;++ix) {
      cell_=iz*cdx*cdy+ix;
      head[cell_]=head[cell_+cdx*(cdy-2)];
      head[cell_+cdx*(cdy-1)]=head[cell_+cdx];
      }
#endif


  /////////////////////////////////////////////
  // Now save min(rc_[]) into curskin so we can use it
  // in check_verlet_list()
  currc = std::numeric_limits<double>::infinity();
  for (i = 0; i < 3; ++i)
    if (rc_[i] < currc)
      currc = rc_[i];

  STOP_TIMER(13)
  START_TIMER(14)



//////////////////////////////////////////////////////
// The following construct makes the integration of a few
// particles very slow compared to optimized single thread
// but for larger system sizes this speeds up calculation considerably
#pragma omp parallel shared(cx,cy,cz,cdx,cdy,cdz,nthreads) private(k,i,ix,iy,iz,ic,cell_,nn)
{
#pragma omp master
{
  nthreads = omp_get_num_threads();
  if (lastcn < TARGETNUMTHREADS)	// only use 1 thread for few particles
    nthreads = 1;
}
#pragma omp barrier	// implies a flush so all threads see the changed value of nthreads

  long vl_i;
  int tid = omp_get_thread_num();
  if (lastcn < TARGETNUMTHREADS && tid != 0)
    tid = -1;


  if (tid != -1) {
    vl_i = nn = 0;
    int size = (int)( ((lastcn*(tid+1))/nthreads-(lastcn*tid)/nthreads)*_pperc*VLIST_MAXES );
    l_vlist[tid] = new long[size];
    for(k=(tid*lastcn)/nthreads;k<((tid+1)*lastcn)/nthreads;++k) {
      ///////////////////////////////
      // for non existent particles set empty vlist entry
      if (_list[k] != -1) {
        l_vlist[tid][nn++] = -1;
        l_vlist[tid][nn++] = -1;
        continue;
      }
      for (i = _pperc*k; i < _pperc*(k+1); ++i) {	// loop over k(1)
        for(ix=cx[i]-1;ix<=cx[i]+1;++ix) 
          for(iy=cy[i]-1;iy<=cy[i]+1;++iy) 
            for(iz=cz[i]-1;iz<=cz[i]+1;++iz) {
              cell_=iz*cdx*cdy+iy*cdx+ix;
              for(ic=head[cell_];ic!=-1;ic=list[ic])
                if (ic!=i && !__SAMECELL(i,ic)) {
                  if (distsq(i,ic) <= currc*currc) {
                    l_vlist[tid][nn++] = ic;
		    /////////////////////////////////////////////
		    // Small note: if nn overflows, this is most probably
		    // caused by some error so that cell density is large locally
		    // Otherwise increase VLIST_MAXES
		    if (nn >= size) {	// This is serious, so handle it NOW
#pragma omp critical
{
		      process_signal_handler(0xFFFB);
		      resolve_int();
}
		    }
                  }
                }
            }
//optional sorting loop
        long *to_sort=l_vlist[tid]+vl_i;
        long done=0;
        long nni=nn-vl_i;
        for(ix=0;!done;++ix) {
          done=1;
          for(iy=1;iy<nni-ix;++iy)
            if (to_sort[iy-1]>to_sort[iy]) {
              iz=to_sort[iy];
              to_sort[iy]=to_sort[iy-1];
              to_sort[iy-1]=iz;
              done=0;
              }
          }
        l_vlist[tid][nn++] = -1;
        vl_i = nn;
#ifdef PGAS
        if (_spec[i] == gas_spec) {		// if gas particle -> don't check k*_pperc+1!!!
          l_vlist[tid][nn++] = -1;
          vl_i = nn;
          break;
        }
#endif
      }		// end loop over k(1)
    }
  }
}	// end omp parallel (implies flush and barrier)


  resolve_int();	// check whether we got a SIGNAL and resolve it!

    ///////////////////////////////
    // Old single thread version //
    ///////////////////////////////

    ///////////////////////////////
    // for non existent particles set empty vlist entry
/*    if (_list[k] != -1) {
      vlist_offset[k*_pperc+1] = nn;
      vlist_offset[k*_pperc+2] = nn;
      vlist_nc[k*_pperc+1] = nn;
      vlist_nc[k*_pperc+2] = nn;
      continue;
      }
    for (i = _pperc*k; i < _pperc*(k+1); ++i) {	// loop over k
      for(ix=cx[i]-1;ix<=cx[i]+1;++ix) 
        for(iy=cy[i]-1;iy<=cy[i]+1;++iy) 
          for(iz=cz[i]-1;iz<=cz[i]+1;++iz) {
            cell_=iz*cdx*cdy+iy*cdx+ix;
            for(ic=head[cell_];ic!=-1;ic=list[ic])
              if (ic!=i && !__SAMECELL(i,ic)) {
                if (distsq(i,ic) <= rc*rc) {
                  vlist[nn++] = ic;
                  }
                }
          }
//optional sorting loop
      long *to_sort=vlist+vlist_offset[i];
      long done=0;
      long nni=nn-vlist_offset[i];
      for(ix=0;!done;++ix) {
        done=1;
        for(iy=1;iy<nni-ix;++iy)
          if (to_sort[iy-1]>to_sort[iy]) {
            iz=to_sort[iy];
            to_sort[iy]=to_sort[iy-1];
            to_sort[iy-1]=iz;
            done=0;
            }
        }
      vlist_nc[i] = nn;				// vlist_nc points to first free place in given verlet list entry
      if (nn-vlist_offset[i] <= VLIST_MES-VLIST_PADDING)
        nn = vlist_offset[i]+VLIST_MES;
      else
        nn += VLIST_PADDING;			// save some free places for new particles
      vlist_offset[i+1]=nn;
      if (nn > VLISTINITSIZE-2*VLIST_MES) {	// vlist is likely to grow beyond its size
        system("touch ./vlistinitsize_reached.dat");
        process_signal_handler(0xFFFD);
    	return ;
        }
      }
    }*/

  ///////////////////////////////////////////
  // now put together real vlist with
  // padding and mean entry size
  int cpn = 0;
  nn = 0;
  vlist_offset[cpn] = 0;

  for (i = 0; i < nthreads; ++i) {
    for (k = 0; k < (int)( (((lastcn*(i+1))/nthreads-(lastcn*i)/nthreads))*_pperc*VLIST_MAXES ) ; ++k) {
      if (l_vlist[i][k] == -1) {
        vlist_nc[cpn] = nn;
        /////////////////////////////////////////////////////////
        // The following ensures each entry is at least VLIST_MES in size
        // and hast at least VLIST_PADDING free entries. This should speed
        // up the introduction of new particles so that we do not need to
        // calculate a new vlist every cell division
        if (nn-vlist_offset[cpn] <= VLIST_MES-VLIST_PADDING)
          nn = vlist_offset[cpn]+VLIST_MES;	// ensure vlist entry is at least VLIST_MES long
        else
          nn += VLIST_PADDING;			// and has at least VLIST_PADDING free entries
        vlist_offset[cpn+1] = nn;
        cpn++;
        if (cpn/_pperc >= ((i+1)*lastcn)/nthreads)
          break;
      }
      else
        vlist[nn++] = l_vlist[i][k];
    }
    delete[] l_vlist[i];
  }
  delete[] l_vlist;

  for (i = 0; i < lastn; ++i)
    vlist_dx[i] = vlist_dy[i] = vlist_dz[i] = 0.;

  delete[] head;

  STOP_TIMER(14)
  STOP_TIMER(9)
}



////////////////////////////////////
// void verlet_list_delete_i(...)
// ---------------
// delete particle i from verlet list (cell death)
// ===
// Parameters:
//   long numi: index of cell particle, i.e offset into _r,...
void verlet_list_delete_i(long numi) {
  long i,j,k,l;
//printf("delete %d\n",numi);

  ///////////////////////////////
  // delete vlist entry for particle numi
  vlist_nc[numi] = vlist_offset[numi];

  for(k=0;k<lastcn;++k) {
    if (_list[k] != -1)
      continue;
    for (i = _pperc*k; i < _pperc*(k+1); ++i) {
      for(j=vlist_offset[i];j<vlist_nc[i];++j)
        if (vlist[j]==numi) {	// move all numbers in this entry one down and decrement vlist_nc[i]
          --vlist_nc[i];
          for (l = j; l < vlist_nc[i]; ++l)
            vlist[l] = vlist[l+1];
          break;
        } // end if
#ifdef PGAS
      if (_spec[i] == gas_spec)	// if gas particle -> don't check k*_pperc+1!!!
        break;
#endif
    } // end for i=k*_pperc
  } // end loop over cells

#ifdef __DEBUG_VLIST
  for (int __j = 0; __j <= lastn; ++__j) {
    if (vlist_offset[__j] < 0) {
      std::cout << "__DEBUG at t=" << curtime << " (del entry) with " << std::endl;
      std::cout << "lastn: " << lastn << std::endl;
      for (int __i = 0; __i <= lastn; ++__i)
        std::cout << __i << ": " << vlist_offset[__i] << std::endl;
      //abort_signal_handler(994);
    }
  }
#endif
}



////////////////////////////////////
// void verlet_list_fill_offset(...)
// ---------------
// fill vlist_offset so that all not defined indices
// are mapped to zero length vlist entries
// ===
// Parameters:
//   long src : index of cell particle to move
//   long dest: index of new cell particle position
//   long tmplastn: max used index of arrays _rc, ...
//   long tmplastn: max used index of arrays _rc, ...
void verlet_list_fill_offset(long oldlastn, long newlastn)
{
  long i;

  for (i = oldlastn; i < newlastn; ++i) {
    vlist_nc[i] = vlist_offset[i];
    vlist_offset[i+1] = vlist_offset[i] + VLIST_MES;
  }

#ifdef __DEBUG_VLIST
  for (int __j = 0; __j <= lastn; ++__j) {
    if (vlist_offset[__j] < 0) {
      std::cout << "__DEBUG at t=" << curtime << " (move entry) with" << std::endl;
      std::cout << "lastn: " << lastn << std::endl;
      for (int __i = 0; __i <= lastn; ++__i)
        std::cout << __i << ": " << vlist_offset[__i] << std::endl;
    }
  }
#endif
}



////////////////////////////////////
// void verlet_list_move_entry(...)
// ---------------
// move entry src to pos dest and adjust index occurrence
//   => here src and dest have to get each other into their vlist entries
// ===
// Parameters:
//   long src : index of cell particle to move
//   long dest: index of new cell particle position
// ===
// Returns:
//   1: if update_verlet_list() was called
//   0: otherwise
int verlet_list_move_entry(long src, long dest, long oldlastn)
{
  long i, j;
  long parent = src-1;
  long off = vlist_offset[dest] - vlist_offset[src];
  long shift;


  /////////////////////////////////
  // check whether there is enough space
  if (dest+1 >= maxn) {
    process_signal_handler(0xFFFD);
    return 1;
  }

  /////////////////////////////////
  // check whether there is enough space to copy entry, otherwise call update_verlet_list()
  // vlist_offset[dest+1]-1 because we need an additional entry for parent
  if (vlist_nc[src] - vlist_offset[src] > vlist_offset[dest+1]-1 - vlist_offset[dest] && dest < oldlastn) {
#ifdef __DEBUG_VLIST
    std::cout << "__DEBUG at t=" << curtime << ": update_vlist()1: dest=" << dest << ", oldlastn=" << oldlastn
              << ",vlist_nc[src]=" << vlist_nc[src] << ",vlist_offset[src]=" << vlist_offset[src]
              << ",vlist_offset[dest+1]-1=" << vlist_offset[dest+1]-1 << ",vlist_offset[dest]="
              << vlist_offset[dest] << std::endl;
#endif
    vlist_rebuild = 1;
    return 1;
  }

  /////////////////////////////////
  // check whether parent has 1 add entry space left
  if (vlist_nc[parent] >= vlist_offset[parent+1]) {
#ifdef __DEBUG_VLIST
    std::cout << "__DEBUG at t=" << curtime << ": update_vlist()2: vlist_nc[parent]=" << vlist_nc[parent]
              << ", vlist_offset[parent+1]=" << vlist_offset[parent+1] << std::endl;
#endif
    vlist_rebuild = 1;
    return 1;
  }

  /////////////////////////////////
  // move src to dest
  for (i = vlist_offset[src]; i < vlist_nc[src]; ++i) {
    vlist[i+off] = vlist[i];
  }
  vlist_nc[dest] = i+off;
  vlist_nc[src] = vlist_offset[src];
  if (dest >= oldlastn)
    for (j = dest+1; j <= lastn; ++j) {
      if (vlist_nc[j-1]-vlist_offset[j-1] <= VLIST_MES-VLIST_PADDING)
        shift = vlist_offset[j-1]+VLIST_MES;
      else
        shift = vlist_nc[j-1]+VLIST_PADDING;
      vlist_nc[j] = vlist_offset[j] = shift;
    }


  /////////////////////////////////
  // adjust indices (i.e. change src to dest)
  // one could also skip the not used entries but its probably more expensive than this
  for (i = 0; i < vlist_offset[lastn]; ++i)
    if (vlist[i] == src)
      vlist[i] = dest;

  /////////////////////////////////
  // insert dest to parent and vice versa
  for (i = vlist_nc[dest]; i > vlist_offset[dest]; --i) {
    if (vlist[i-1] < parent) {
      vlist[i] = parent;
      break;
    }
    else
      vlist[i] = vlist[i-1];
  }
  if (i == vlist_offset[dest])
    vlist[i] = parent;
  vlist_nc[dest]++;

  for (i = vlist_nc[parent]; i > vlist_offset[parent]; --i) {
    if (vlist[i-1] < dest) {
      vlist[i] = dest;
      break;
    }
    else
      vlist[i] = vlist[i-1];
  }
  if (i == vlist_offset[parent])
    vlist[i] = dest;
  vlist_nc[parent]++;

  /////////////////////////////////
  // adjust displacement vector
  vlist_dx[dest] = vlist_dx[src];
  vlist_dy[dest] = vlist_dy[src];
  vlist_dz[dest] = vlist_dz[src];

#ifdef __DEBUG_VLIST
  for (int __j = 0; __j <= lastn; ++__j) {
    if (vlist_offset[__j] < 0) {
      std::cout << "__DEBUG at t=" << curtime << " (move entry) with __j=" << j << std::endl;
      std::cout << "lastn: " << lastn << std::endl;
      for (int __i = 0; __i <= lastn; ++__i)
        std::cout << __i << ": " << vlist_offset[__i] << std::endl;
    }
  }
#endif

  return 0;
}



////////////////////////////////////
// void verlet_list_insert_i(...)
// ---------------
// insert particle i into vlist
// ===
// Parameters:
//   long dest: index of new cell particle position
//   long oldlastn: max used index of arrays _rc, ...
// ===
// Returns:
//   1: if update_verlet_list() was called
//   0: otherwise
int verlet_list_insert_i(long dest, long oldlastn) {
  int i,j,k;
  long src=dest-1;
  long off;
  long shift;


  ////////////////////////////////////
  // Check whether vlist has enough free space.
  // Increase VLISTINITSIZE accordingly
  if (vlist_offset[lastn]+VLIST_MAXES+vlist_offset[src+1]-vlist_offset[src] >= VLISTINITSIZE) {
    process_signal_handler(0xFFFD);
    return 1;
  }

  ////////////////////////////////////
  // check whether there is space left to insert i at dest
  if (vlist_offset[dest+1]-vlist_offset[dest] < vlist_nc[src] - vlist_offset[src] && dest < oldlastn) {
#ifdef __DEBUG_VLIST
    std::cout << "__DEBUG at t=" << curtime << ": update_vlist()3: dest=" << dest << ", oldlastn=" << oldlastn
              << ",vlist_offset[dest+1]=" << vlist_offset[dest+1] << ",vlist_offset[dest]=" << vlist_offset[dest]
              << ",vlist_nc[src]=" << vlist_nc[src] << ",vlist_offset[src]=" << vlist_offset[src] << std::endl;
#endif
    vlist_rebuild = 1;
    return 1;
  }


  ///////////////////////////////////
  // sweep through vlist and insert dest for every occurrence of src
  for (i = 0; i < lastn; ++i) {
    for (j = vlist_offset[i]; j < vlist_nc[i]; ++j) {
      if (vlist[j] == src) {
        if (vlist_nc[i]+1 > vlist_offset[i+1]) {
#ifdef __DEBUG_VLIST
          std::cout << "__DEBUG at t=" << curtime << ": update_vlist()4: vlist_nc[i]+1=" << vlist_nc[i]+1
                    << ", vlist_offset[i+1]=" << vlist_offset[i+1] << std::endl;
#endif
          vlist_rebuild = 1;
          return 1;
        }
        for (k = vlist_nc[i]-1; k > j; --k)
          vlist[k+1] = vlist[k];
        vlist[j+1] = dest;
        vlist_nc[i]++;
      }
    }
  }

  ///////////////////////////////////
  // copy new particle to its position inside vlist
  for (i = vlist_offset[src], off = vlist_offset[dest]-i; i < vlist_nc[src]; ++i)
    vlist[i+off] = vlist[i];
  vlist_nc[dest] = vlist_nc[src]+off;

  if (dest >= oldlastn)
    for (j = dest+1; j <= lastn; ++j) {
      if (vlist_nc[j-1]-vlist_offset[j-1] <= VLIST_MES-VLIST_PADDING)
        shift = vlist_offset[j-1]+VLIST_MES;
      else
        shift = vlist_nc[j-1]+VLIST_PADDING;
      vlist_nc[j] = vlist_offset[j] = shift;
    }

  ///////////////////////////////////
  // copy displacement vector
  vlist_dx[dest] = vlist_dx[src];
  vlist_dy[dest] = vlist_dy[src];
  vlist_dz[dest] = vlist_dz[src];

#ifdef __DEBUG_VLIST
  for (int __j = 0; __j <= lastn; ++__j) {
    if (vlist_offset[__j] < 0) {
      std::cout << "__DEBUG at t=" << curtime << " (insert i) with __j=" << j << std::endl;
      std::cout << "lastn: " << lastn << std::endl;
      for (int __i = 0; __i <= lastn; ++__i)
        std::cout << __i << ": " << vlist_offset[__i] << std::endl;
    }
  }
#endif

  return 0;
}


////////////////////////////////////
// void validate_vlist(...)
// ---------------
// creates new vlist and compares interaction partners to existing one
// I am not sure whether this function still gives sane results, so do
// a thorough check before you use it (and then remove this line)!!!
// ===
// Parameters:
void validate_vlist()
{
  long i,j,k;
  int found;
  double dx,dy,dz,drsq;

  ///////////////////////////////
  // calc interaction partners and check whether in vlist
  for (i = 0; i < lastn; ++i) {
    if (_list[i/_pperc] != -1)
      continue;
    for (j = 0; j < lastn; ++j) {
      if (_list[j/_pperc] != -1 || i==j || __SAMECELL(i,j))
        continue;
      dx = distx(i,j);
      dy = disty(i,j);
      dz = distz(i,j);
      drsq = dx*dx + dy*dy + dz*dz + EPSILON;
      found = 1;
      if (drsq < rpp[0]*rpp[0]) {
        found = 0;
        for (k = vlist_offset[i]; k < vlist_nc[i]; ++k)
          if (vlist[k] == j) {
            found = 1;
            break;
          }
      }
      if (!found) {
        std::cout << "__DEBUG at t=" << curtime << ": validate_vlist(): i=" << i << "; j="
                  << j << " not found in vlist: lastcn=" << lastcn << ": lastn=" << lastn << std::endl;

        for (k = 0; k < lastcn; ++k) {
	  if (_list[k] != -1)
	    continue;
	  for (int m = k*_pperc; m < _pperc*(k+1); ++m)
  	    for (int l = vlist_offset[m]; l < vlist_nc[m]; ++l)
	      std::cout << "vlist[]: i=" << m << ", j=" << vlist[l] << ": vo=" << vlist_offset[m]
	                << ": vlist_nc=" << vlist_nc[m] << ": vo[i+1]=" << vlist_offset[m+1] << std::endl;
	}
	process_signal_handler(0xFFF0);
	resolve_int();
	return ;
      }
    }

    //////////////////////////////////
    // check whether there are double entries in any verlet list entry
    found = 0;
    for (j = vlist_offset[i]; j < vlist_nc[i]; ++j)
      for (k = j+1; k < vlist_nc[i]; ++k)
        if (vlist[j] == vlist[k])
          found++;

    if (found) {
      std::cout << "__DEBUG at t=" << curtime << ": validate_vlist(): i=" << i << "; j="
                << j << " found double entry: lastcn=" << lastcn << ": lastn=" << lastn << std::endl;

      for (k = 0; k < lastcn; ++k) {
	if (_list[k] != -1)
	  continue;
	for (int m = k*_pperc; m < _pperc*(k+1); ++m)
  	  for (int l = vlist_offset[m]; l < vlist_nc[m]; ++l)
	    std::cout << "vlist[]: i=" << m << ", j=" << vlist[l] << ": vo=" << vlist_offset[m]
	              << ": vlist_nc=" << vlist_nc[m] << ": vo[i+1]=" << vlist_offset[m+1] << std::endl;
      }
      process_signal_handler(0xFFF0);
      resolve_int();
      return ;
    }

    //////////////////////////////////
    // check that not from same cell
    for (j = vlist_offset[i]; j < vlist_nc[i]; ++j)
      if (__SAMECELL(i,vlist[j])) {
        std::cout << "__DEBUG at t=" << curtime << ": validate_vlist(): i=" << i << "; j="
                  << j << " from same cell: lastcn=" << lastcn << ": lastn=" << lastn << std::endl;

        for (k = 0; k < lastcn; ++k) {
	  if (_list[k] != -1)
	    continue;
	  for (int m = k*_pperc; m < _pperc*(k+1); ++m)
  	    for (int l = vlist_offset[m]; l < vlist_nc[m]; ++l)
	      std::cout << "vlist[]: i=" << m << ", j=" << vlist[l] << ": vo=" << vlist_offset[m]
	                << ": vlist_nc=" << vlist_nc[m] << ": vo[i+1]=" << vlist_offset[m+1] << std::endl;
	}
	process_signal_handler(0xFFF0);
	resolve_int();
	return ;
      }
  }
}



