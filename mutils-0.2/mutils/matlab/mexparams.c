/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#include "mexparams.h"

#define get_matlab_int_class(type, class_out)	\
  {						\
    if(sizeof(type)==1){			\
      if(IS_TYPE_SIGNED(Int))			\
	class_out = mxINT8_CLASS;		\
      else					\
	class_out = mxUINT8_CLASS;		\
    } else if(sizeof(type)==2){			\
      if(IS_TYPE_SIGNED(Int))			\
	class_out = mxINT16_CLASS;		\
      else					\
	class_out = mxUINT16_CLASS;		\
    } else if(sizeof(type)==4){			\
      if(IS_TYPE_SIGNED(Int))			\
	class_out = mxINT32_CLASS;		\
      else					\
	class_out = mxUINT32_CLASS;		\
    } else { /* if(sizeof(Int)==8) */		\
      if(IS_TYPE_SIGNED(type))			\
	class_out = mxINT64_CLASS;		\
      else					\
	class_out = mxUINT64_CLASS;		\
    }						\
  }

#define get_matlab_int_class_name(type, class_out)	\
  {							\
    if(sizeof(type)==1){				\
      if(IS_TYPE_SIGNED(Int))				\
	class_out = "int8";				\
      else						\
	class_out = "uint8";				\
    } else if(sizeof(type)==2){				\
      if(IS_TYPE_SIGNED(Int))				\
	class_out = "int16";				\
      else						\
	class_out = "uint16";				\
    } else if(sizeof(type)==4){				\
      if(IS_TYPE_SIGNED(Int))				\
	class_out = "int32";				\
      else						\
	class_out = "uint32";				\
    } else { /* if(sizeof(Int)==8) */			\
      if(IS_TYPE_SIGNED(type))				\
	class_out = "int64";				\
      else						\
	class_out = "uint64";				\
    }							\
  }


Double *mex_get_matrix_Double(const mxArray *param, Int *m, Int *n, 
			      const char *varname, const char *sm, const char *sn, int can_be_empty)
{
  Int _m = 0, _n = 0;
  char buff[256] = {0};

  if(!param || mxIsEmpty(param)){
    if(can_be_empty) return NULL;
    USERERROR("'%s' can not be empty", MUTILS_INVALID_PARAMETER, varname);
  }

  if(sizeof(Double)==sizeof(double)){
    if(!mxIsClass(param, "double")) USERERROR("'%s' must be of type 'double'", MUTILS_INVALID_PARAMETER, varname);
  } else {
    if(!mxIsClass(param, "single")) USERERROR("'%s' must be of type 'single'", MUTILS_INVALID_PARAMETER, varname);
  }

  SNPRINTF(buff, 255, "No dimensions of '%s' can be larger than %lli", varname, (long long)MaxInt);
  managed_type_cast(Int, _m, mxGetM(param), buff);
  managed_type_cast(Int, _n, mxGetN(param), buff);

  SNPRINTF(buff, 255, "Dimensions of '%s' should be [%s X %s]", varname, sm, sn);

  if((*m) && (*m)!=_m) USERERROR("%s", MUTILS_INVALID_PARAMETER, buff);
  if((*n) && (*n)!=_n) USERERROR("%s", MUTILS_INVALID_PARAMETER, buff);

  *m = _m;
  *n = _n;

  return (Double*)mxGetData(param);
}


mxArray *mex_set_matrix_Double(Double *values, Int m, Int n)
{
  mxArray *outp;

  if(!values){
    if(sizeof(Double)==sizeof(double))
      outp = mxCreateNumericMatrix(m,n,mxDOUBLE_CLASS,mxREAL);
    else
      outp = mxCreateNumericMatrix(m,n,mxSINGLE_CLASS,mxREAL);
  } else {
    if(sizeof(Double)==sizeof(double))
      outp = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,mxREAL);
    else
      outp = mxCreateNumericMatrix(0,0,mxSINGLE_CLASS,mxREAL);
    mxSetM(outp, m);
    mxSetN(outp, n);
    mxSetData(outp, values);
  }
  
  return outp;
}


Int *mex_get_matrix_Int(const mxArray *param, Int *m, Int *n, 
			const char *varname, const char *sm, const char *sn, int can_be_empty)
{
  Int _m = 0, _n = 0;
  char buff[256] = {0};
  char *class_name;

  if(!param || mxIsEmpty(param)){
    if(can_be_empty) return NULL;
    USERERROR("'%s' can not be empty", MUTILS_INVALID_PARAMETER, varname);
  }

  get_matlab_int_class_name(Int, class_name);
  
  if(!mxIsClass(param, class_name))
    USERERROR("'%s' must be of type '%s'", MUTILS_INVALID_PARAMETER, varname, class_name);

  SNPRINTF(buff, 255, "No dimensions of '%s' can be larger than %lli", varname, (long long)MaxInt);
  managed_type_cast(Int, _m, mxGetM(param), buff);
  managed_type_cast(Int, _n, mxGetN(param), buff);

  SNPRINTF(buff, 255, "Dimensions of '%s' should be [%s X %s]", varname, sm, sn);

  if((*m) && (*m)!=_m) USERERROR("%s", MUTILS_INVALID_PARAMETER, buff);
  if((*n) && (*n)!=_n) USERERROR("%s", MUTILS_INVALID_PARAMETER, buff);

  *m = _m;
  *n = _n;

  return (Int*)mxGetData(param);
}


mxArray *mex_set_matrix_Int(Int *values, Int m, Int n)
{
  mxArray *outp;
  mxClassID class_id;

  get_matlab_int_class(Int, class_id);

  if(!values){
    outp = mxCreateNumericMatrix(m,n,class_id,mxREAL);
  } else {
    outp = mxCreateNumericMatrix(0,0,class_id,mxREAL);
    mxSetM(outp, m);
    mxSetN(outp, n);
    mxSetData(outp, values);
  }
  
  return outp;
}


Int mex_get_integer_scalar(const mxArray *param, const char *varname, int can_be_empty, Int def)
{
  Int out = 0;
  int sign;
  char buff[256] = {0};
  mxClassID class_id;
  long long stemp;
  unsigned long long utemp;
  double dtemp;

  if(IS_TYPE_SIGNED(Int)) sign = 1; else sign=0;

  if(!param || mxIsEmpty(param)){
    if(can_be_empty) return def;
    USERERROR("'%s' can not be empty", MUTILS_INVALID_PARAMETER, varname);
  }

  if(mxGetM(param)!=1 || mxGetN(param)!= 1)
    USERERROR("'%s' must be a scalar.", MUTILS_INVALID_PARAMETER, varname);

  SNPRINTF(buff, 255, 
	   "Value of '%s' does not match the internal Int type. Must be a maximum %lu bits %s integer.", 
	   varname, sizeof(Int)*8, sign ? "signed": "unsigned");
  
  class_id = mxGetClassID(param);
  switch(class_id){

  case mxINT8_CLASS:
  case mxINT16_CLASS:
  case mxINT32_CLASS:
  case mxINT64_CLASS:
    if(class_id == mxINT8_CLASS) stemp = ((int8_T*)mxGetData(param))[0];
    else if(class_id == mxINT16_CLASS) stemp = ((int16_T*)mxGetData(param))[0];
    else if(class_id == mxINT32_CLASS) stemp = ((int32_T*)mxGetData(param))[0];
    else stemp = ((int64_T*)mxGetData(param))[0];
    managed_type_cast(Int, out, stemp, buff);
    return out;

  case mxUINT8_CLASS:
  case mxUINT16_CLASS:
  case mxUINT32_CLASS:
  case mxUINT64_CLASS:
    if(class_id == mxUINT8_CLASS) utemp = ((uint8_T*)mxGetData(param))[0];
    else if(class_id == mxUINT16_CLASS) utemp = ((uint16_T*)mxGetData(param))[0];
    else if(class_id == mxUINT32_CLASS) utemp = ((uint32_T*)mxGetData(param))[0];
    else utemp = ((uint64_T*)mxGetData(param))[0];
    managed_type_cast(Int, out, utemp, buff);
    return out;
    
  case mxDOUBLE_CLASS:
  case mxSINGLE_CLASS:
    if(class_id == mxDOUBLE_CLASS) dtemp = ((double*)mxGetData(param))[0];
    else dtemp = ((float*)mxGetData(param))[0];
    out = (Int)dtemp;
    if(dtemp != (double)out) USERERROR("%s", MUTILS_INVALID_PARAMETER, buff);
    return out;

  default:
    USERERROR("'%s' must be either of an integer, or a real type.", MUTILS_INVALID_PARAMETER, varname);
    break;
  }

  return 0;
}
