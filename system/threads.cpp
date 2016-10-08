#pragma once

#include "../system/threads.h"


class Mutex
{
	friend struct Condition;
	MUTEX m;
public:
	Mutex() { mutex_init(m); };
	~Mutex() {};

	void lock() { mutex_lock(m); }
	void unlock() { mutex_unlock(m); }	
	MUTEX mutex() { return m; }
};


struct Condition
{
	CONDITION c;
public:
	Condition() { thread_cond_init(c); }
	~Condition() { cond_destroy(c); }

	void wait(Mutex& m) { thread_wait(c, m.mutex()); }
#ifdef __linux
	void wait_for(NativeMutex& m, timespec& d) { thread_timed_wait(c, m.mutex, d); }
#else
	void wait_for(Mutex& m, int ms) { thread_timed_wait(c, m.mutex(), ms); }
#endif
	void notify_one() { thread_signal(c); }
};


class Thread
{
	size_t idx;
	THREAD thread;
	Mutex mutex;
	Condition condition;
	void * data;
	thread_fnc work_fnc;
	THREAD_HANDLE h;

public:
	Thread(int id) : idx(id), work_fnc(0), h(0) { };
	Thread(int id, thread_fnc tf) : idx(id), work_fnc(tf), h(0) { };
	Thread(int id, thread_fnc tf, void * dta) : idx(id), work_fnc(tf), data(dta), h(0) { };
	~Thread() { if (work_fnc) { work_fnc = NULL; } };
	
	void start(void * dat) { h = thread_create(thread, work_fnc, dat); }
	void start() { h = thread_create(thread, work_fnc, data); }
	void join() { if (h) thread_join(&h); }
	int id() { return idx; }
};



