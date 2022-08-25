#pragma once

#include <string>
#include <iostream>
#include <vector>

using namespace std;

#include "utils.hpp"

struct Timer {
	void start();
	void stop();

	long val;
};

extern vector <Timer> timer;

void timers_init();
void output_timers(string file);

