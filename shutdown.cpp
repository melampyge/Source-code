/*		shutdown.cpp		*/
//////////////////////////////////////////
// signal handling and sim shutdown	//
//////////////////////////////////////////


#include "shutdown.h"


using namespace std;


////////////////////////////////////
// process_signal_handler(...)
// ---------------
// function invoked on signal raise
// ===
// Parameters:
//   int s: signal number
void process_signal_handler(int s) {
  _int_flag = 1;
  _int_no   = s;

  if (!omp_in_parallel())
    resolve_int();
}



////////////////////////////////////
// void init_signal_handler(...)
// ---------------
// initialize signal handler(s)
// ===
// Parameters:
void init_signal_handler()
{
  struct sigaction sigIntHandler;

  sigIntHandler.sa_handler = process_signal_handler;
  sigemptyset(&sigIntHandler.sa_mask);
  sigIntHandler.sa_flags = 0;


  ////////////////////////////////////
  // install signal handling for SIGINT (CTRL-C) and SIGTERM
  // BUT: program can still be killed with kill -9
  sigaction(SIGINT, &sigIntHandler, NULL);
  sigaction(SIGTERM, &sigIntHandler, NULL);
  ////////////////////////////////////
  // install signal handling for SIGUSR1 (print program status to stdout)
  sigaction(SIGUSR1, &sigIntHandler, NULL);
}


////////////////////////////////////
// int resolve_int(...)
// ---------------
// check for interrupt occurrence and resolve
// ===
// Parameters:
int resolve_int() {
  string str;

  if (_int_flag == 0)
    return 0;
  _int_flag = 0;

  //////////////////////////////////
  // tell the world what happened
  switch (_int_no) {
  case 2:
    str = string("SIGINT");
    break;
  case 15:
    str = string("SIGTERM");
    break;
  case 0xFFFB:
    str = string("l_vlist too small");
    break;
  case 0xFFFC:
    str = string("Thread number larger or equal to MAXTHREADNUM");
    break;
  case 0xFFFD:
    str = string("max vlist size reached!");
    break;
  case 0xFFFE:
    str = string("Trying to delete non existent cell particle!");
    break;
  case 0xFFFF:
    str = string("Max array size reached!");
    break;
  case 0xFFFFFFFF:
    str = string("__DEBUG_TERM");
    break;
  default:
    stringstream stmp;
    stmp << _int_no;
    str = string("SIGNUM '") + stmp.str() + string("'");
    break;
  }
  str += " has been invoked, shutting down...";

  if (_int_no == SIGUSR1) {
    print_status();
    return 0;
  }

  cerr << str << endl;

  //////////////////////////////////
  // save all data in DUMPFILE
  dump_all();

#ifdef __DEBUG_TIMER
  OUTPUT_TIMERDATA((string(OUTDIR)+string("__timing.dat")).c_str());
#endif

  long t_diff = time(0) - t_begin;
  long t_s = t_diff % 60;
  long t_m = (t_diff/60) % 60;
  long t_h = (t_diff/3600) % 24;
  long t_d = t_diff/(3600*24);
  cout << "Process ran for >> " << t_d << "d " << setw(2) << setfill('0') << t_h << ":"
       << setw(2) << setfill('0') << t_m << ":"
       << setw(2) << setfill('0') << t_s << " << " << t_diff << "s" << endl;

  ///////////////////////////////////
  // delete lock file
  unlink((string(OUTDIR)+string("lock.dat")).c_str());

  //////////////////////////////////
  // exit
  //abort();
  _Exit(_int_no);
}



////////////////////////////////////
// int print_status(...)
// ---------------
// Print current status of simulation (add content as necessary)
// ===
// Parameters:
void print_status()
{
  long t_diff = time(0) - t_begin;
  long t_s = t_diff % 60;
  long t_m = (t_diff/60) % 60;
  long t_h = (t_diff/3600) % 24;
  long t_d = t_diff/(3600*24);

  cout << "Status: t=" << curtime << ": lastcn=" << lastcn << ": lastn=" << lastn
       << ": next_free_id=" << next_free_id << endl;

  cout << "Run time: " << t_d << "d " << setw(2) << setfill('0') << t_h << ":"
       << setw(2) << setfill('0') << t_m << ":"
       << setw(2) << setfill('0') << t_s << " << " << t_diff << "s" << endl;
#ifdef WALLTIMESTOP
  cout << "End Time: " << WTRUNTIME << ": diff= " << WTRUNTIME-t_diff << endl;
#endif

  cout << "Verlet list: _head=" << _head << ": vl_off[lastn]=" << vlist_offset[lastn] << endl;

}



