/*		localstress.h		*/
//////////////////////////////////////////
// functions related to determing the	//
// local stress				//
// IMPORTANT: Actually the pressure	//
// is determined and not the stress!!!	//
//////////////////////////////////////////


#ifndef __LOCALSTRESS_H_
#define __LOCALSTRESS_H_

#include <vector>
#include <algorithm>
#include "parameter.h"
#include "model.h"




void transform_coordinate_pair(double *, double *);

void localstress_update(int, int, double, int);
void localstressd_update(int, int, double, int);

//void localstress_meanvelocity(int, int);
void localstress_momentum_update(int, int);
























#endif
