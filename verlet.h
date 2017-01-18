/*		verlet.h		*/
//////////////////////////////////////////
// functions related to velocity	//
// verlet algorithm			//
//////////////////////////////////////////


#ifndef __VERLET_H_
#define __VERLET_H_

#include <cstring>
#include <cmath>
#include <omp.h>

#include "parameter.h"
#include "io_func.h"
#include "bc_func.h"
#include "model.h"
#include "divdeath.h"
#include "verletlist.h"
#include "debug.h"
#include "localstress.h"



void init_verlet();

void verlet_step_dpd();

void calcv(long);

void calcr(long);


#endif
