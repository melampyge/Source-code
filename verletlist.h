/*		verletlist.h		*/
//////////////////////////////////////////
// functions to create and maintain	//
// the verlet list			//
//////////////////////////////////////////


#ifndef __VERLETLIST_H_
#define __VERLETLIST_H_

#include <limits>

#include "parameter.h"
#include "model.h"
#include "shutdown.h"

#if defined (__DEBUG_VLIST) || defined (__DEBUG_1)
#include <iostream>
#endif


#define VLIST_PADDING	8	// how many additional free entries per particle
#define VLIST_MAXES	256	// defines max entry size for parallel verlet list calculation
// This is the AVERAGE maximum entry size, so individual entries may be bigger
// as long as average is not. Tune this if you change interaction range or
// more general the particle size (if its too small this causes a SIG(0xFFFB)).
// -----------------------
// For maximum performance one should adjust these values for a given parameter set,
// up until now they are more like educated guesses...
/////////////////////////////////////////


void check_verlet_list();
void check_list(long,char *);

void update_verlet_list();

void verlet_list_delete_i(long);

void verlet_list_fill_offset(long, long);

int verlet_list_move_entry(long, long, long);

int verlet_list_insert_i(long, long);

//////////////////////////////////////
// CUDA Version, which I never got running
//void update_verlet_list_cuda(double *, double *, double *, double, long);


void validate_vlist();

#endif
