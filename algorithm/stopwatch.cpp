#include "stopwatch.h"


void Clock::start() { gettimeofday(&strt, NULL); }
void Clock::stop() { gettimeofday(&ed,NULL); }
double Clock::seconds()  { return (double) (ed.tv_sec - strt.tv_sec); }
double Clock::useconds()  { return (double) (ed.tv_usec - strt.tv_usec) ; }
double Clock::ms()  { return ( (seconds()) * 1000 + useconds()/1000) ;}
double Clock::elapsed_sec()  { stop(); double secs = seconds(); start(); return secs;}
double Clock::elapsed_ms()  { stop(); double msecs = ms(); start(); return msecs; }
