#ifndef _ARGS_H
#define _ARGS_H

#include <ParseArgv.h>
#include <EBTKS/Path.h>		/* (bert) - Added EBTKS subdirectory */

class args {
public:
  MString command;
  Path    inputPath;
  Path    outputPath;

  static int verbose;
  static int clobber;

  static double fwhm;   // filter parameters
  static double noise;  //
  static int given_min;  // true -> user supplied minimum 
  static int given_max;
  static int blur_flag;  // true -> skip deblurring step
  static double range[2];
  static int debug_flag;

  static ArgvInfo argTable[];

  args(int argc, char **argv);
  ~args() {};
};

#endif




