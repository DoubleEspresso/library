#ifndef DIAGNOSTICS_STOPWATCH_H_
#define DIAGNOSTICS_STOPWATCH_H_

#include "../system/types.h"

#ifdef _WIN32
#include <time.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip>
class Clock
{
	clock_t start_time;
	clock_t stop_time;
	uint64 start_ticks;
	uint64 stop_ticks;
public:
	Clock() {};
	~Clock() {};
	void start() { start_time = clock(); }
	void stop() { stop_time = clock(); }
	double seconds() { return (stop_time - start_time) / double(CLOCKS_PER_SEC); }
	double ms() { return seconds() * double(1000); }
	double elapsed_sec() { stop(); double secs = seconds(); start(); return secs; }
	double elapsed_ms() { stop(); double msecs = ms(); start(); return msecs; }

	void print(const char *desc)
	{
		stop();
		double secs = seconds();
		double ms = (secs > 0) ? (1000 * secs) : 0;
		printf("%s\t%7.3f(s)\t%7.0f(ms)\n", desc, secs, ms);
	}
};
#else
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>

class Clock
{
	struct timeval strt, ed;
public:
	Clock() {};
	~Clock() {};
	void start() { gettimeofday(&strt, NULL); }
	void stop() { gettimeofday(&ed, NULL); }
	double seconds() { return (double)(ed.tv_sec - strt.tv_sec); }
	double useconds() { return (double)(ed.tv_usec - strt.tv_usec); }
	double ms() { return ((seconds()) * 1000 + useconds() / 1000); }
	double elapsed_sec() { stop(); double secs = seconds(); start(); return secs; }
	double elapsed_ms() { stop(); double msecs = ms(); start(); return msecs; }
	void print(const char *desc)
	{
		stop();
		double secs = seconds();
		double ms = (secs > 0) ? (1000 * secs) : 0;
		printf("%s\t%7.3f(s)\t%7.0f(ms)\n", desc, secs, ms);
	}
};
#endif

#endif
