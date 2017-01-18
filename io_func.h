/*                  io_func.h               */
//////////////////////////////////////////////
// contains all input/output realted        //
// functions to easily change format specs  //
//////////////////////////////////////////////


#ifndef __IO_FUNC_H_
#define __IO_FUNC_H_

#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>

#include "parameter.h"
#include "basics.h"




void init_output();

int load_dump(char *, int load_add_param = 0);

void output_trajectory();

void output_numcells();

void output_pressure();

void output_density();

void output_localstress();

void output_curp();

void dump_all();

#endif



