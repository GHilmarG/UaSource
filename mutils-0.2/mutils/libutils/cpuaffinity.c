/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#include "config.h"
#include "cpuaffinity.h"

#ifndef WINDOWS
#include <sys/syscall.h>
#include <unistd.h>
#include <sched.h>
#endif

#include <errno.h>
#include <string.h>

#ifdef USE_NUMA
#include <numa.h>
#endif

int affinity_bind(int thrid, int cpu)
{
  int retval = -1;

#ifdef WINDOWS
#pragma message ("CPU affinity not implemented on Windows")
  thrid = 0;
  cpu = 0;
#else
#ifdef USE_NUMA
    {
      nodemask_t nodemask;
      nodemask_t *mask=&nodemask;
      nodemask_zero(mask);
      nodemask_set(mask, cpu);
      numa_run_on_node_mask(mask);
      numa_set_membind(mask);
      numa_bind(mask);
      numa_set_bind_policy(1);
      numa_set_strict(1);
      MESSAGE("thread %li: numa_preferred: %i", thrid, numa_preferred());
      retval = 0;
    }
#else
    {
      pid_t tid = syscall(SYS_gettid);
      MESSAGE("thrid %d, tid %d, bind to CPU %d", thrid, tid, cpu);
      unsigned int cpusetsize = sizeof(cpu_set_t);
      cpu_set_t mask;
      CPU_ZERO(&mask);
      CPU_SET(cpu, &mask);
      retval = sched_setaffinity(tid, cpusetsize, &mask);
      MESSAGE("setaff %d: %s", retval, strerror(errno));
    }
#endif
#endif
    
  return retval;
}


