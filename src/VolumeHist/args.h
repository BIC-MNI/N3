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
$RCSfile: args.h,v $
$Revision: 1.1 $
$Author: bert $
$Date: 2003-04-16 14:31:39 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _ARGS_H
#define _ARGS_H

#include <ParseArgv.h>
#include <EBTKS/Path.h>		/* (bert) - Added EBTKS subdirectory */

class args {
public:
  MString command;
  Path    maskPath;
  Path    inputPath;
  Path    outputPath;

  static int verbose;
  static int clobber;
  static int matlab_format;  // produce files in matlab format
  static int number_of_bins;
  static int number_of_classes;
  static int limit_classes;  // flag for using number of classes
  static int auto_range;     // flag for computing a range for each class
  static int selected_class; // only compute histogram for selected class
  static int select_class_flag; // only compute one histogram if true
  static int mask_flag;      // true -> use mask
  static int window_flag;       // true -> use Parzen window
  

  static ArgvInfo argTable[];

  args(int argc, char **argv);
  ~args() {};
};

#endif

