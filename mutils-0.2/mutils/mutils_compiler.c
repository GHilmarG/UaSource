#include <matrix.h>
#include <stdlib.h>
#include <stdio.h>

void mexFunction(int nargout, mxArray *pargout [ ], int nargin, const mxArray *pargin[])
{
#ifdef __ICC
  pargout[0] = mxCreateString("icc");
#endif
#ifdef _MSC_VER
  pargout[0] = mxCreateString("cl");
#endif
#if defined(__GNUC__) & !defined(__ICC)
  pargout[0] = mxCreateString("gcc");
#endif
}
