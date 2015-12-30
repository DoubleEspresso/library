#ifndef CONCURRENT_THREADS_H
#define CONCURRENT_THREADS_H

#include <vector>
#include "../system/platform.h"


class NativeMutex
{
 public:
  NativeMutex() 
    {
#ifdef _WIN32
      InitializeCriticalSection(&mutex);
#else 
      mutex_init(mutex);
#endif
    };
  ~NativeMutex() {};
  
  void lock() { mutex_lock(mutex); }
  void unlock() { mutex_unlock(mutex); }
  MUTEX mutex;
 private:
  friend struct ConditionVariable;
};


struct ConditionVariable 
{
  ConditionVariable() { thread_cond_init(c); }
  ~ConditionVariable() { cond_destroy(c); }
  
  void wait(NativeMutex& m) { thread_wait(c, m.mutex); }
#ifdef __linux
  void wait_for(NativeMutex& m, timespec& d) { thread_timed_wait(c, m.mutex, d); }
#else
  void wait_for(NativeMutex& m, int ms) { thread_timed_wait(c, m.mutex, ms); }
#endif
  void notify_one() { thread_signal(c); }
  
private:
  CONDITION c;
};



class Thread
{
 public:
 Thread(int id) : idx(id), work_fnc(0) { };    
 Thread(int id, thread_fnc tf) : idx(id), work_fnc(tf) { };    
 Thread(int id, thread_fnc tf, void * dta) : idx(id), work_fnc(tf), data(dta) { };    
  ~Thread() 
    { 
      if (work_fnc) work_fnc = NULL; 
    };  

  thread_fnc work_fnc;  
  long start(void * dat) { return thread_create(handle, work_fnc, dat);}
  long start() { return thread_create(handle, work_fnc, data);}
  int join() { return thread_join(handle); }
  int id() { return idx; }

 private:
  size_t idx;
  THREAD handle;
  NativeMutex mutex;
  ConditionVariable cond;
  void * data;
};

class ThreadPool : public std::vector<Thread*>
{
public:
  void init();
  void exit();
  Thread* master() { return static_cast<Thread*>(at(0)); }
  Thread* available_slave(const Thread * master) const;
  
  NativeMutex mutex;
  ConditionVariable sleep_condition;  
};


namespace 
{
  extern "C"
  {
    long do_async(Thread *thread, void * data) { thread->work_fnc(data); return 0;}
  }
}

#endif
  
