/*		shutdown.h		*/
//////////////////////////////////////////
// signal handling and sim. shutdown	//
//////////////////////////////////////////



#ifndef __SHUTDOWN_H_
#define __SHUTDOWN_H_

#include "parameter.h"
#include "io_func.h"


#include <iostream>
#include <string>
#include <csignal>
#include <cstdlib>
#include <unistd.h>
#include <omp.h>



void process_signal_handler(int);

void init_signal_handler();

int resolve_int();

void print_status();



#endif
