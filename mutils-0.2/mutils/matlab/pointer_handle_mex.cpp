#ifdef  _GNU_SOURCE
# define __USE_GNU      1
# define __USE_XOPEN2K8 1
#endif

#include <libutils/config.h>
#include <libutils/debug_defs.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mex.h>
#include <matrix.h>

#include <list>

static std::list<void*> ptr_list;

void pointer_handle_mex_cleanup(void) {
    
  /* find the pointer in the managed list */
  std::list<void*>::iterator i;
  for(i=ptr_list.begin(); i!=ptr_list.end(); i++){
    free(*i);
  }    
}


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

mxArray *set_pointer_handle(void *ptr, const char *name)
{
  mxArray *handle_ptr    = NULL;
  mxArray *handle_own    = NULL;
  mxArray *handle_name   = NULL;

  fprintf(stderr, "set_pointer_handle\n");

  /* set the pointer and the owner */
  handle_ptr = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
  ((long long*)mxGetData(handle_ptr))[0] = (long long)ptr;

  /* name of the pointer class */
  handle_name = mxCreateCharMatrixFromStrings(1, &name);

  /* owner of the pointer frees the memory */
  handle_own = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
  ((long*)mxGetData(handle_own))[0] = 1;

  /* create pointer handle object - call the class constuctor through MATLAB directly */
  mxArray *obj = NULL;
  mxArray *rhs[3] = {handle_ptr, handle_name, handle_own};
  mexCallMATLAB(1,&obj,3,rhs,"PointerHandle");

  if(!obj){
    /* not good - memory leak if pointer is not freed */
    free(ptr);
    mexErrMsgTxt("could not create a PointerHandle object.");
  }

  /* register pointer handle */
  const char *cmd = "register";
  rhs[0] = obj;
  rhs[1] = mxCreateCharMatrixFromStrings(1, &cmd);

  /* TODO not good - memory leak if the below call fails */
  mexCallMATLAB(0,NULL,2,rhs,"pointer_handle");

  return obj;
}


void *get_pointer_handle(mxArray *ptr)
{
  /* verify we deal with valid PointerHandle objects */
  if(strcmp(mxGetClassName(ptr), "PointerHandle")) {
    mexErrMsgTxt("object not of type PointerHandle");
  }

  /* validate the pointer handle */
  const char *cmd = "isvalid";
  mxArray *isvalid;
  mxArray *rhs[2];
  rhs[0] = ptr;
  rhs[1] = mxCreateCharMatrixFromStrings(1, &cmd);

  mexCallMATLAB(1,&isvalid,2,rhs,"pointer_handle");
  if(isvalid==NULL){
    mexErrMsgTxt("unable to validate the PointerHandle object. Verify that pointer_handle mex file is available.");
  }
  if(((long*)mxGetData(isvalid))[0] == 0){
    mexErrMsgTxt("the PointerHandle object points to an invalid memory area.");
  } 

  ptr = mxGetProperty(ptr, 0, "pointer");
  if(!ptr){
    mexErrMsgTxt("invalid PointerHandle object (possibly deallocated before?)");
  } 
  
  long *temp = (long*)mxGetData(ptr);
  if(!temp){
    mexErrMsgTxt("invalid PointerHandle object (possibly deallocated before?)");
  }

  return (void*)(temp[0]);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */


#ifndef LIB_POINTER_HANDLE
void mexFunction(int nargout, mxArray *pargout [ ], int nargin, const mxArray *pargin[])
{
  const mxArray *handle_ptr = NULL;
  long          *temp       = NULL;
  void          *pointer    = NULL;
  char           cmd[256];

  mexAtExit(pointer_handle_mex_cleanup);

  if(nargin<2) return;

  /* get the command */
  if(!mxIsChar(pargin[1])){
    ERROR("command parameter must be a string");
  }
  mxGetString(pargin[1], cmd, 255);
  

  /* verify we deal with valid PointerHandle objects */
  handle_ptr = pargin[0];
  if(strcmp(mxGetClassName(handle_ptr), "PointerHandle")) {
    if(strcmp(cmd, "isvalid")){
      ERROR("object not of type PointerHandle");
    }
  }

  handle_ptr = mxGetProperty(handle_ptr, 0, "pointer");
  if(!handle_ptr){
    if(strcmp(cmd, "isvalid")){
      ERROR("invalid PointerHandle object");
    }
  } else {
    temp = (long*)mxGetData(handle_ptr);
    if(!temp){
      if(strcmp(cmd, "isvalid")){
	ERROR("invalid PointerHandle object");
      }
    } else {
      pointer    = (void*)(temp[0]);
    }
  }


  /***************************************/
  /* register the pointer for future handling */
  /***************************************/
  if(!strcmp(cmd, "register")){

    std::list<void*>::iterator i;
    for(i=ptr_list.begin(); i!=ptr_list.end(); i++){
      if(pointer == *i) break;
    }

    if(i==ptr_list.end()){

      /* store the pointer in the list */
      ptr_list.push_back(pointer);
    }
    return;
  }


  /***************************************/
  /* delete the pointer, free the memory */
  /***************************************/
  if(!strcmp(cmd, "delete")){

    /* find the pointer in the list */
    std::list<void*>::iterator i;
    for(i=ptr_list.begin(); i!=ptr_list.end(); i++){
      if(pointer == *i){

	/* pointer found - delete and remove from the list */
	ptr_list.erase(i);
	/* free(pointer); */
	return;
      }
    }    
    
    /* invalid pointer. return without error */
    return;
  }
  

  /***************************************/
  /* check if pointer is valid */
  /***************************************/
  if(!strcmp(cmd, "isvalid")){
    
    mxArray *out = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
    ((long*)mxGetData(out))[0] = 0;

    std::list<void*>::iterator i;
    for(i=ptr_list.begin(); i!=ptr_list.end(); i++){
      if(pointer == *i){
	((long*)mxGetData(out))[0] = 1;
	break;
      }
    }
    pargout[0] = out;
    return;
  }
  
  ERROR("unknown command");
}

#endif /* POINTER_HANDLE_LIB */
