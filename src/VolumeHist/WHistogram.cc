/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, John G. Sled, 
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- 
$RCSfile: WHistogram.cc,v $
$Revision: 1.1 $
$Author: bert $
$Date: 2003-04-16 14:31:39 $
$State: Exp $
--------------------------------------------------------------------------*/
#include "WHistogram.h"
#include <assert.h>

#ifdef HAVE_MATLAB
extern "C" {    
  #include"mat.h" 
}
#endif

