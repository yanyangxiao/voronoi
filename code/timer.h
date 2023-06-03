
#ifndef TIMER_H
#define TIMER_H

#include <time.h>

class Timer
{
protected:
	clock_t _start;
	clock_t _stop;

public:
	Timer()
	{
		_start = clock();
	}

	void start()
	{
		_start = clock();
	}

	void stop()
	{
		_stop = clock();
	}

	float elapsed_time() const
	{
		return float(_stop - _start) / 1000;
	}
};

#endif

