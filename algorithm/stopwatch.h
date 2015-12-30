#ifndef STOPWATCH_H_
#define STOPWATHCH_H

#include <time.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip>

//#ifdef LINUX
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>

class Clock
{
public:
  Clock() {};
  ~Clock() {};
  void start();
  void stop();
  double seconds();
  double useconds();
  double ms();
  double elapsed_sec();
  double elapsed_ms();

 private:
  struct timeval strt, ed;
};
/*
#elif WIN32

class Clock
{
 public:
  Clock() {};
  ~Clock() {};

  void start() {start_time = clock();}
  void stop()  {stop_time = clock();}

  double seconds () { return (stop_time - start_time) / double(CLOCKS_PER_SEC) ;}
  double ms() { return seconds() * double ( 1000 ); }
  double elapsed_sec() { stop(); double secs = seconds(); start(); return secs;}
  double elapsed_ms() { stop(); double msecs = ms(); start(); return msecs; }

  void print(const char *desc)
    {
      stop();
      double secs=seconds();
      double ms = (secs > 0) ? (1000*secs) : 0;

      printf("%s\t%7.3f(s)\t%7.0f(ms)\n",
	     desc, secs, ms);
    }
 private:
  clock_t start_time;
  clock_t stop_time;
  U64 start_ticks;
  U64 stop_ticks;  
};
#endif
*/
#endif
