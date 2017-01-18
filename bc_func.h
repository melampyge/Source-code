/*		bc_func.h		*/
//////////////////////////////////////////
// functions related to boundaries are  //
// declared here			//
//////////////////////////////////////////


#ifndef __BC_FUNC_H_
#define __BC_FUNC_H_

#include "parameter.h"
#include "pressure.h"

#ifdef __DEBUG
#include <iostream>
#include "shutdown.h"
#endif


int outofBound(long);

int check_boundcollision(long, int);

void bc_vchange(long, int);


#endif
