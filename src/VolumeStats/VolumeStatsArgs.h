/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, Alex P. Zijdenbos, 
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
$RCSfile: VolumeStatsArgs.h,v $
$Revision: 1.1 $
$Author: bert $
$Date: 2003-04-16 14:32:14 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef AVERAGE_VOLUMES_ARGS_H
#define AVERAGE_VOLUMES_ARGS_H

#include <ParseArgv.h>
#include <EBTKS/Path.h>		/* (bert) - Added EBTKS subdirectory */

const unsigned _MAX_ITEMS_IN_LIST = 256;

class VolumeStatsArgs {
public:
  MString     command;
  Array<Path> inputPath;
  Path        maskPath;

  static int verbose;
  static int quiet;
  static int debug;
  static int cached;
  static int useWorldCoord;
  static int ignoreNaN;
  static int all;
  static int calcVolume;
  static int calcMin;
  static int calcMax;
  static int calcSum;
  static int calcSum2;
  static int calcMean;
  static int calcVariance;
  static int calcStddev;
  static int calcMedian;
  static int calcMajority;
  static int calcBiModalT;
  static int calcPctT;
  static int calcEntropy;
  static int calcCoM;

  static char    *axisOrder[3];
  static int      sliceList[_MAX_ITEMS_IN_LIST];
  static double   pctT;
  static double   binValue[_MAX_ITEMS_IN_LIST], maskRange[2];
  static int      maskBinValue[_MAX_ITEMS_IN_LIST];
  static double   floor[_MAX_ITEMS_IN_LIST], ceiling[_MAX_ITEMS_IN_LIST];
  static double   maskFloor[_MAX_ITEMS_IN_LIST], maskCeil[_MAX_ITEMS_IN_LIST];
  static unsigned nExtrema, nMaskExtrema;
  static char    *maskString;

  static ArgvInfo argTable[];

  static int getList(char *dst, char *key, char *nextArg);
  static int getIntList(char *dst, char *key, char *nextArg);
  static int getAxisOrder(char *dst, char *key, char *nextArg);

  // Constructors/destructor
  VolumeStatsArgs(int argc, char **argv);
  ~VolumeStatsArgs() {}
};

#endif

