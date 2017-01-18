/*		debug.h		*/
//////////////////////////////////////////
// Debugging helper functions		//
//////////////////////////////////////////


#ifndef __DEBUG_H_
#define __DEBUG_H_

#include <iostream>
#include <cmath>

#include "parameter.h"
#include "shutdown.h"
#include "bc_func.h"


//////////////////////////////////////////////////
// if debugging, enable asserts
#ifdef __DEBUG
#define __ASSERT(i, j, k, l) __assert_((i),(j),(k),(l))
#else
#define __ASSERT(i, j, k, l) {}
#endif


void __assert_(long, long, long, const char *);
void __assert_(int, int, int, const char *);
void __assert_(double, double, double, const char *);

int check_rvffd();
int check_list_lastcn();



#endif

