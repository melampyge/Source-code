/*		divdeath.cpp		*/
//////////////////////////////////////////
// cell division and death		//
//////////////////////////////////////////


#ifndef __DIVDEATH_H_
#define __DIVDEATH_H_

#include "parameter.h"
#include "bc_func.h"
#include "model.h"
#include "basics.h"
#include "shutdown.h"
#include "verletlist.h"
#include "debug.h"


void cell_death();
void cell_division();

int do_division(long);
int delete_cell(long);


#endif
