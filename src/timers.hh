#ifndef BEEPMBP__TIMERS_HH
#define BEEPMBP__TIMERS_HH

#include <string>
#include <iostream>
#include <vector>

using namespace std;

#include "utils.hh"

struct Timer { 
  void start();
  void stop();
  
  long val;
};

extern vector <Timer> timer;

void timers_init();
void output_timers(string file);
#endif
