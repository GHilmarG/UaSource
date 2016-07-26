/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#ifndef DEBUG_DEFS_H
#define DEBUG_DEFS_H

#include "config.h"
#include <stdlib.h>
#include <stdio.h>

#define DEBUG_NONE 0
#define DEBUG_BASIC 1
#define DEBUG_DETAILED 5
#define DEBUG_MEMORY 6
#define DEBUG_ALL 10

#ifdef MATLAB_MEX_FILE
#include <mex.h>
#define WARNING_HDR ""
#define ERROR_HDR ""
#define MSG_NEWLINE ""
#define PRINTF_MSG(msg) mexPrintf("%s", msg)
#define PRINTF_ERR(msg, id) mexErrMsgIdAndTxt(id, msg);
#define PRINTF_WRN(msg, id) mexWarnMsgIdAndTxt(id, msg);
#else
#define WARNING_HDR " ** WARNING: "
#define ERROR_HDR " ** ERROR: "
#define MSG_NEWLINE "\n"
#define PRINTF_MSG(msg) printf("%s", msg)
#define PRINTF_ERR(msg, id) fprintf(stderr, "%s %s", id, msg);
#define PRINTF_WRN(msg, id) printf("%s %s", id, msg);
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  int get_debug_mode(void);
  void set_debug_mode(int debug);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#undef MESSAGE
#define MESSAGE(msg, ...) {						\
    char __buff[256];							\
    SNPRINTF(__buff, 255, " ** " msg "\n", ##__VA_ARGS__); \
    __buff[255] = 0;							\
    PRINTF_MSG(__buff);fflush(stdout);					\
  }

#undef EMESSAGE
#define EMESSAGE(msg, ...) {						\
    char __buff[256];							\
    SNPRINTF(__buff, 255, ERROR_HDR "%s: %s(): %d: " msg MSG_NEWLINE, __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    __buff[255] = 0;							\
    PRINTF_ERR(__buff, "");						\
  }

#undef WMESSAGE
#define WMESSAGE(msg, ...) {						\
    char __buff[256];							\
    SNPRINTF(__buff, 255, WARNING_HDR " %s: %s(): %d: " msg MSG_NEWLINE, __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    __buff[255] = 0;							\
    PRINTF_WRN(__buff, "");						\
  }

#undef ERROR
#define ERROR(msg,...) {			\
    EMESSAGE(msg, ##__VA_ARGS__);		\
    exit(1);					\
  }

#undef WARNING
#define WARNING(msg, ...) {			\
    WMESSAGE(msg, ##__VA_ARGS__);		\
  }

/* User-versions for mex files */
#undef USERERROR
#define USERERROR(msg, id, ...) {					\
    char __buff[256];							\
    SNPRINTF(__buff, 255, ERROR_HDR msg MSG_NEWLINE, ##__VA_ARGS__);	\
    __buff[255] = 0;							\
    PRINTF_ERR(__buff, id);						\
  }
#undef USERWARNING
#define USERWARNING(msg, id, ...) {					\
    char __buff[256];							\
    SNPRINTF(__buff, 255, WARNING_HDR msg MSG_NEWLINE, ##__VA_ARGS__); \
    __buff[255] = 0;							\
    PRINTF_WRN(__buff, id);						\
  }

#define HERE fprintf(stderr, "%s: %s(): %i: HERE\n", __FILE__, __FUNCTION__, __LINE__);

#ifndef DEBUG
#undef DMESSAGE
#define DMESSAGE(msg, level, ...)
#define TODO(msg, ...)
#define F_ENTER
#define F_EXIT
#else
#undef DMESSAGE
#define DMESSAGE(msg, level, ...) {					\
    if(level<=get_debug_mode()){					\
      char __buff[256];							\
      SNPRINTF(__buff, 255, " ** %s: %s(): %d: " msg "\n", __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__); \
      __buff[255] = 0;							\
      PRINTF_MSG(__buff);fflush(stdout);				\
    }									\
  }
#define TODO(msg, ...) {						\
    char __buff[256];							\
    SNPRINTF(__buff, 255, "TODO %s: %s(): %d: " msg "\n", __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    __buff[255] = 0;							\
    PRINTF_MSG(__buff);fflush(stdout);					\
  }


#define F_ENTER							\
  {								\
    printf("\n --------------------------\n");			\
    printf(" -- ENTER %s:%s\n", __FILE__, __FUNCTION__);	\
    printf(" --------------------------\n");			\
  }								\


#define F_EXIT						\
  {							\
    printf(" --------------------------\n");		\
    printf(" -- EXIT %s:%s\n", __FILE__, __FUNCTION__);	\
    printf(" --------------------------\n\n");		\
  }							\

#endif

#endif

