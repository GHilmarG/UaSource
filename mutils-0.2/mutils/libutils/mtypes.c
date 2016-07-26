/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#include "mtypes.h"

tfloat macheps(void)
{
  tfloat machEps = 1.0;
  do {
    machEps /= (tfloat)2.0;
  } while (( (tfloat)1.0 + machEps/(tfloat)2.0) != 1.0 );
  return machEps;
}
