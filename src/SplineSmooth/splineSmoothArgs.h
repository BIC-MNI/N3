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
$RCSfile: splineSmoothArgs.h,v $
$Revision: 1.1 $
$Author: bert $
$Date: 2003-04-16 14:30:54 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : splineSmoothArgs.h,v
@INPUT      : 
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: argument parsing class for SplineSmooth
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 21, 1996 (John G. Sled)
@MODIFIED   : Log: splineSmoothArgs.h,v 
 * Revision 1.2  1996/04/23  13:37:00  jgsled
 * Working version.  Problems with thin plate splines have been fixed.
 *
 * Revision 1.1  1996/04/21  23:41:51  jgsled
 * Initial version of SplineSmooth tool
 * - B spline implementation appears to work
 *
@COPYRIGHT  : 1996
---------------------------------------------------------------------------- */
#ifndef _ARGS_H
#define _ARGS_H

#include <ParseArgv.h>
#include <EBTKS/Path.h>		/* (bert) - Added EBTKS subdirectory */
#include "splineSmooth.h"

class args {
public:
  MString command;
  Path    maskPath;
  Path    outputMaskPath;
  Path    inputPath;
  Path    outputPath;
  Path    compactPath;

  static int verbose;
  static int clobber;
  static double lambda;      // scale invariant smoothing parameter
  static double distance;    // distance between basis functions
  static enum spline_type spline;  
  static int use_mask;
  static int use_output_mask;
  static int extrapolate;
  static int full_support;
  static int subsample;       // skip every nth point
  static int produce_volume;
  static int produce_compact;

  static ArgvInfo argTable[];

  args(int argc, char **argv);
  ~args() {};
};

#endif

