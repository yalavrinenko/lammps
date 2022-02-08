#ifndef GETTIMEOFDAY_H
#define GETTIMEOFDAY_H

//#include <config.h>
#include <time.h>
#include <sys/timeb.h>
//#include "../include/time.h"

//#include <sysinfoapi.h>

struct timeval 
{
    time_t tv_sec;
    time_t tv_usec;
};


inline int GetTickCount(){
  clock_t t = clock();
  return (int)((1000*(float)t)/CLOCKS_PER_SEC);
}





inline int gettimeofday(struct timeval *tp, void *tzp)
{
    
    struct _timeb timebuffer;
    
    _ftime(&timebuffer);
    tp->tv_sec = timebuffer.time;
    tp->tv_usec = timebuffer.millitm * 1000;
    
    return 0;
    

}

#endif /* GETTIMEOFDAY_H */