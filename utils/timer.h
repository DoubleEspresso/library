#ifndef UTILS_TIMER_H
#define UTILS_TIMER_H

#include <time.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip>

#include "../system/types.h"

#ifdef __linux
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
class lTimer
{
	struct timeval strt, ed;
public:
	lTimer() {};
	~lTimer() {};

	void start() { gettimeofday(&strt, NULL); }
	void stop() { gettimeofday(&ed, NULL); }
	double seconds() { return (double)(ed.tv_sec - strt.tv_sec); }
	double useconds() { return (double)(ed.tv_usec - strt.tv_usec); }
	double ms() { return ((seconds()) * 1000 + useconds() / 1000); }
	double elapsed_sec() { stop(); double secs = seconds(); start(); return secs; }
	double elapsed_ms() { stop(); double msecs = ms(); start(); return msecs; }
	void print(const char *desc)
	{
		double secs = seconds();
		double ms = (secs > 0) ? (1000 * secs) : 0;
		printf("%s\t%7.3f(s)\t%7.0f(ms)\n", desc, secs, ms);
	}
};
typedef lTimer Timer;

#else // windows

class wTimer
{
	clock_t start_time;
	clock_t stop_time;
	uint64 start_ticks;
	uint64 stop_ticks;
public:
	wTimer() {};
	~wTimer() {};

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
typedef wTimer Timer;
#endif

#endif
