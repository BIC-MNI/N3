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
$RCSfile: evaluateFieldArgs.h,v $
$Revision: 1.1 $
$Author: bert $
$Date: 2003-04-16 14:25:59 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : evaluateFieldArgs.h,v
@INPUT      : 
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : July 9, 1996 John G. Sled
@MODIFIED   : 
@COPYRIGHT  : 1996
---------------------------------------------------------------------------- */
#ifndef _ARGS_H
#define _ARGS_H

#include <ParseArgv.h>
#include <EBTKS/Path.h>		/* (bert) - Added EBTKS subdirectory */
#include "../SplineSmooth/splineSmooth.h"

class evaluateArgs {
public:
  MString command;
  Path    maskPath;
  Path    likePath;
  Path    inputPath;
  Path    outputPath;

  static int verbose;
  static int clobber;
  static int use_mask;
  static nc_type datatype;
  static int signed_flag;

  static ArgvInfo ArgTable[];

  evaluateArgs(int argc, char **argv);
  ~evaluateArgs() {};
};

#endif




