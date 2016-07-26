/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#include "tictoc.h"
#include "debug_defs.h"

static double flops = 0;
static double bytes = 0;
static long time_us = 0;

#ifdef WINDOWS
#else
static struct timeval tb, te;

# ifdef USE_OPENMP
#  pragma omp threadprivate(tb, te, flops, bytes)
# endif /* USE_OPENMP */

#endif



#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#endif

void stats_zero(void)
{
  flops = 0;
  time_us = 0;
}

void flops_add(double nflops)
{
  flops += nflops;
}

void bytes_add(double nbytes)
{
  bytes += nbytes;
}

double flops_get()
{
  return flops;
}

double bytes_get()
{
  return bytes;
}


#ifdef WINDOWS
#else
void _tic()
{
  gettimeofday(&tb, NULL);
  flops=0;
  bytes=0;
  fflush(stdout);
}

void _toc()
{
  long s,u;
  double tt;
  gettimeofday(&te, NULL);
  s=te.tv_sec-tb.tv_sec;
  u=te.tv_usec-tb.tv_usec;
  tt=((double)s)*1000000+u;
  MESSAGE("time:                  %li.%.6li", (s*1000000+u)/1000000, (s*1000000+u)%1000000);
  MESSAGE("MFLOPS:                %.3lf", flops/tt);
  MESSAGE("total fp operations:   %.0lf", flops);
  MESSAGE("MB/s:                  %.3lf", bytes/tt);
  MESSAGE("total memory traffic   %.0lf", bytes);
  fflush(stdout);
}

void _midtoc()
{
  long s,u;
  double tt;
  gettimeofday(&te, NULL);
  s=te.tv_sec-tb.tv_sec;
  u=te.tv_usec-tb.tv_usec;
  tt=((double)s)*1000000+u;
  MESSAGE("time:                  %li.%.6li", (s*1000000+u)/1000000, (s*1000000+u)%1000000);
  MESSAGE("MFLOPS:                %.3lf\n", flops/tt);
  fflush(stdout);
}

void _ntoc(const char *idtxt)
{
  long s,u;
  gettimeofday(&te, NULL);
  s=te.tv_sec-tb.tv_sec;
  u=te.tv_usec-tb.tv_usec;
  if(idtxt){
    MESSAGE("%s: time:                  %li.%.6li", idtxt, (s*1000000+u)/1000000, (s*1000000+u)%1000000);
  } else {
    MESSAGE("time:                  %li.%.6li", (s*1000000+u)/1000000, (s*1000000+u)%1000000);
  }
  fflush(stdout);
}

void _nntoc()
{
  long s,u;
  gettimeofday(&te, NULL);
  s=te.tv_sec-tb.tv_sec;
  u=te.tv_usec-tb.tv_usec;
  printf("%li.%.6li", (s*1000000+u)/1000000, (s*1000000+u)%1000000);
  fflush(stdout);
}

void _inctime()
{
  gettimeofday(&te, NULL);
  time_us += (te.tv_sec-tb.tv_sec)*1000000 + (te.tv_usec-tb.tv_usec);
}

void stats_print()
{
  MESSAGE("total time %li.%6li", time_us/1000000, time_us%1000000);
}

#endif /* WINDOWS */


