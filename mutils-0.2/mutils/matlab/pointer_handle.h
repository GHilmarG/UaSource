#ifndef _POINTER_HANDLE_H
#define _POINTER_HANDLE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* do not build mexFunction, only the below utility functions */
#ifndef LIB_POINTER_HANDLE
#define LIB_POINTER_HANDLE
#endif

  mxArray *set_pointer_handle(void *ptr, const char *name);
  void *get_pointer_handle(mxArray *ptr);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _POINTER_HANDLE_H */
