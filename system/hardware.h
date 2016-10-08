#pragma once

#ifndef SYSTEM_HARDWARE_H
#define SYSTEM_HARDWARE_H

#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif


namespace Hardware
{
  unsigned cpu_count(bool &hyperthreading);
};


inline unsigned Hardware::cpu_count(bool &hyperthreading)  
{
  uint32_t registers[4];
  unsigned logicalCPUS;
  unsigned physicalCPUS;
#ifdef _WIN32
  SYSTEM_INFO si;
  GetSystemInfo( &si );
  logicalCPUS = si.dwNumberOfProcessors;
#else
  logicalCPUS = sysconf ( _SC_NPROCESSORS_ONLN );

  __asm__ __volatile__("cpuid " :
		       "=a" (registers[0]),
		       "=b" (registers[1]),
		       "=c" (registers[2]),
		       "=d" (registers[3])
		       : "a" (1), "c" (0));
  
  unsigned CPUfeatureSet = registers[3];
  hyperthreading = CPUfeatureSet & ( 1 << 28 );
#endif  
  if (hyperthreading) physicalCPUS = logicalCPUS / 2;
  else physicalCPUS = logicalCPUS;  
  return physicalCPUS;
}

#endif
