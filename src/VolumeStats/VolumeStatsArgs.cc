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
$RCSfile: VolumeStatsArgs.cc,v $
$Revision: 1.2 $
$Author: stever $
$Date: 2003-11-17 04:30:56 $
$State: Exp $
--------------------------------------------------------------------------*/
#include "VolumeStatsArgs.h"
#include <assert.h>
#include <ctype.h>
#include <EBTKS/Minc.h>		// (bert) - Added EBTKS subdirectory
#include <version.h>

using namespace std;


int VolumeStatsArgs::verbose = FALSE;
int VolumeStatsArgs::quiet = FALSE;
int VolumeStatsArgs::debug = FALSE;
int VolumeStatsArgs::cached = TRUE;
int VolumeStatsArgs::useWorldCoord = FALSE;
int VolumeStatsArgs::ignoreNaN = TRUE;
int VolumeStatsArgs::all = TRUE;
int VolumeStatsArgs::calcVolume = FALSE;
int VolumeStatsArgs::calcMin = FALSE;
int VolumeStatsArgs::calcMax = FALSE;
int VolumeStatsArgs::calcSum = FALSE;
int VolumeStatsArgs::calcSum2 = FALSE;
int VolumeStatsArgs::calcMean = FALSE;
int VolumeStatsArgs::calcVariance = FALSE;
int VolumeStatsArgs::calcStddev = FALSE;
int VolumeStatsArgs::calcMedian = FALSE;
int VolumeStatsArgs::calcMajority = FALSE;
int VolumeStatsArgs::calcBiModalT = FALSE;
int VolumeStatsArgs::calcPctT = FALSE;
int VolumeStatsArgs::calcEntropy = FALSE;
int VolumeStatsArgs::calcCoM = FALSE;

char    *VolumeStatsArgs::axisOrder[3] = { ANY_SPATIAL_DIMENSION, 
					   ANY_SPATIAL_DIMENSION, 
					   ANY_SPATIAL_DIMENSION };
int      VolumeStatsArgs::sliceList[];
double   VolumeStatsArgs::pctT = -MAXDOUBLE;
double   VolumeStatsArgs::binValue[];
double   VolumeStatsArgs::floor[];
double   VolumeStatsArgs::ceiling[];
int      VolumeStatsArgs::maskBinValue[];
double   VolumeStatsArgs::maskRange[2] = { -MAXDOUBLE, MAXDOUBLE };
double   VolumeStatsArgs::maskFloor[];
double   VolumeStatsArgs::maskCeil[];
unsigned VolumeStatsArgs::nMaskExtrema = 0;
unsigned VolumeStatsArgs::nExtrema     = 0;
char    *VolumeStatsArgs::maskString   = 0;

ArgvInfo VolumeStatsArgs::argTable[] = {
  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "General options:"},
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::verbose, 
   "Print out log messages as processing is being done."},
  {"-quiet", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::quiet, 
   "Print requested values only, without labels."},
  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::debug, 
   "Print debugging junk."},
  {"-version", ARGV_FUNC, (char *) print_version_info, 
   (char *) MNI_LONG_VERSION, "Print out version info and exit."},

  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "\nVolume input options:"},
  {"-nocache", ARGV_CONSTANT, (char *) FALSE, (char *) &VolumeStatsArgs::cached,
   "Do not use volume caching."},
  {"-transverse", ARGV_FUNC, (char *) getAxisOrder, 
   (char *) &VolumeStatsArgs::axisOrder,
   "Read in transverse slices"},
  {"-sagittal", ARGV_FUNC, (char *) getAxisOrder, 
   (char *) &VolumeStatsArgs::axisOrder,
   "Read in sagittal slices"},
  {"-coronal", ARGV_FUNC, (char *) getAxisOrder, 
   (char *) &VolumeStatsArgs::axisOrder,
   "Read in coronal slices"},
  {"-dimorder", ARGV_FUNC, (char *) getAxisOrder, 
   (char *) &VolumeStatsArgs::axisOrder,
   "Specify dimension order (<dim1>,<dim2>,<dim3>,...)."},
  
  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "\nVoxel selection options:"},
  {"-include_nan", ARGV_CONSTANT, (char *) FALSE, (char *) &VolumeStatsArgs::ignoreNaN,
   "Include NaN values in counts and calculations."},
  {"-slice", ARGV_FUNC, (char *) VolumeStatsArgs::getIntList,
   (char *) &VolumeStatsArgs::sliceList,
   "Only collect data from the supplied slice(s)"},
  {"-ceil", ARGV_FUNC, (char *) VolumeStatsArgs::getList, 
   (char *) &VolumeStatsArgs::ceiling,
   "Ignore voxels with a value above this argument.\n\t\t (a quoted or comma-separated list produces multiple outputs)"},
  {"-floor", ARGV_FUNC, (char *) VolumeStatsArgs::getList, 
   (char *) &VolumeStatsArgs::floor,
   "Ignore voxels with a value below this argument.\n\t\t (a quoted or comma-separated list produces multiple outputs)"},
  {"-binvalue", ARGV_FUNC, (char *) VolumeStatsArgs::getList,
   (char *) &VolumeStatsArgs::binValue,
   "Same as -floor <value>-0.5 -ceil <value>+0.5.\n\t\t (a quoted or comma-separated list produces multiple outputs)"},
  {"-mask", ARGV_STRING, (char *) 1, (char *) &VolumeStatsArgs::maskString,
   "Ignore voxels for which the argument is zero."},
  {"-useWorldCoord", ARGV_CONSTANT, (char *) TRUE, 
   (char *) &VolumeStatsArgs::useWorldCoord,
   "Use world coordinates when matching the mask."},
  {"-maskbinvalue", ARGV_FUNC, (char *) VolumeStatsArgs::getIntList, 
   (char *) &VolumeStatsArgs::maskBinValue,
   "Only use mask voxels with this integer value (+/- 0.5).\n\t\t (a quoted or comma-separated list produces multiple outputs)"},
  {"-maskfloor", ARGV_FUNC, (char *) VolumeStatsArgs::getList,
   (char *) VolumeStatsArgs::maskFloor,
   "Only use mask voxels greater than or equal to <value>.\n\t\t (a quoted or comma-separated list produces multiple outputs)"},
  {"-maskceil", ARGV_FUNC, (char *) VolumeStatsArgs::getList, 
   (char *) VolumeStatsArgs::maskCeil,
   "Only use mask voxels less than or equal to <value>.\n\t\t (a quoted or comma-separated list produces multiple outputs)"},
  {"-maskrange", ARGV_FLOAT, (char *) 2, (char *) VolumeStatsArgs::maskRange,
   "Only use mask voxels in this range."},

  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "\nOutput options (values are always printed in this order):"},
  {"-all", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::all,
   "Print all stats (default)."},
  {"-none", ARGV_CONSTANT, (char *) FALSE, (char *) &VolumeStatsArgs::all,
   "Print voxel counts only."},
  {"-volume", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcVolume,
   "Print volume (in mm3)."},
  {"-min", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcMin,
   "Print minimum value."},
  {"-max", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcMax,
   "Print maximum value."},
  {"-sum", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcSum,
   "Print sum."},
  {"-sum2", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcSum2,
   "Print sum of squares."},
  {"-mean", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcMean,
   "Print mean value."},
  {"-var", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcVariance,
   "Print variance."},
  {"-stddev", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcStddev,
   "Print standard deviation."},
  {"-median", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcMedian,
   "Print median value."},
  {"-majority", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcMajority,
   "Print most frequently occurring value."},
  {"-biModalT", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcBiModalT,
   "Print optimal threshold that separates the volume into two classes."},
  {"-pctT", ARGV_FLOAT, (char *) 1, (char *) &VolumeStatsArgs::pctT,
   "Print thresholds at the supplied % of data."},
  {"-entropy", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcEntropy,
   "Print the entropy of the volume."},
  {"-CoM", ARGV_CONSTANT, (char *) TRUE, (char *) &VolumeStatsArgs::calcCoM,
   "Print the center of mass of the volume."},

  {NULL, ARGV_END, NULL, NULL, NULL}
};
   
/* ----------------------------- MNI Header -----------------------------------
@NAME       : Constructor of VolumeStatsArgs
@INPUT      : argc - number of command-line arguments
              argv - command-line arguments
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Builds a VolumeStatsArgs object from commandline arguments
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 15, 1995 (Alex Zijdenbos)
@MODIFIED   : 
---------------------------------------------------------------------------- */

VolumeStatsArgs::VolumeStatsArgs(int argc, char **argv)
: command(argv[0]),
  inputPath(0)
{
  unsigned i;

  for (i = 1; i < argc; i++) {
    command += " ";
    command += argv[i];
  }

  // Initialize option lists
  for (i = 0; i < _MAX_ITEMS_IN_LIST; i++) {
    sliceList[i]    = -1;
    maskBinValue[i] = -MAXINT;
    maskFloor[i]    = floor[i]   = binValue[i] = -MAXDOUBLE;
    maskCeil[i]     = ceiling[i] = MAXDOUBLE;
  }

  set_program_name(argv[0]);  // for version info

  // Call ParseArgv
  if (ParseArgv(&argc, argv, argTable, 0) || (argc < 1)) {
    cerr << endl << "Usage: " << argv[0] 
	 << " [<options>] <file1> [<file2> ...]" << endl;
    cerr <<         "       " << argv[0] << " [-help]" << endl << endl;
    exit(EXIT_FAILURE);
  }

  if (cached) {
    set_n_bytes_cache_threshold(0);           // Always cache volume
    set_default_max_bytes_in_cache(524288);   // 0.5M cache size
    set_cache_block_sizes_hint(SLICE_ACCESS); // Cache volume slices
  }

  nExtrema = 0;
  while ((nExtrema < _MAX_ITEMS_IN_LIST) && 
	 ((floor[nExtrema] > -MAXDOUBLE) || (ceiling[nExtrema] < MAXDOUBLE)))
    nExtrema++;

  // Append binvalue list to floor/ceil lists
  i = 0;
  while (binValue[i] > -MAXDOUBLE) {
    floor[nExtrema]   = binValue[i] - 0.5;
    ceiling[nExtrema] = binValue[i] + 0.5;
    nExtrema++;
    i++;
  }

  if (maskString) {
    maskPath = maskString;

    nMaskExtrema = 0;
    while ((nMaskExtrema < _MAX_ITEMS_IN_LIST) && 
	   ((maskFloor[nMaskExtrema] > -MAXDOUBLE) || 
	    (maskCeil[nMaskExtrema] < MAXDOUBLE)))
      nMaskExtrema++;
    
    if ((maskRange[0] != -MAXDOUBLE) || (maskRange[1] != MAXDOUBLE)) {
      maskFloor[nMaskExtrema] = maskRange[0];
      maskCeil[nMaskExtrema]  = maskRange[1];
      nMaskExtrema++;
    }
    
    unsigned j = 0;
    while ((nMaskExtrema < _MAX_ITEMS_IN_LIST) && 
	   (maskBinValue[nMaskExtrema] > -MAXINT)) {
      maskFloor[nMaskExtrema] = maskBinValue[j] - 0.5;
      maskCeil[nMaskExtrema]  = maskBinValue[j] + 0.5;
      nMaskExtrema++;
      j++;
    }
  }

  for (i = 0; i < argc - 1; i++)
    inputPath.append(argv[i + 1]);

  if (pctT > -MAXDOUBLE)
    calcPctT = TRUE;

  if (all && 
      (calcVolume || 
       calcMin || calcMax || calcSum || calcSum2 || calcMean || calcVariance || 
       calcStddev || calcMedian || calcMajority || calcBiModalT || calcPctT ||
       calcEntropy || calcCoM))
    all = FALSE;
}

int
VolumeStatsArgs::getList(char *dst, char *key, char *nextArg)
{
  // Check for next argument
  if (nextArg == NULL) {
    cerr << key << " option requires an additional argument" << endl;
    exit(EXIT_FAILURE);
  }

  double *values = (double *) dst;

  // Set up pointers to end of string and first non-space character
  char *cur = nextArg;
  while (isspace(*cur)) cur++;
  unsigned n = 0;

  // Loop through string looking for space or comma-separated values
  while ((n < _MAX_ITEMS_IN_LIST) && (*cur != '\0')) {
    *values++ = atof(cur);
    n++;

    // Similar to getIntList: without the following line, this function 
    // will ignore arguments on SGI Origin200 when the binary is compiled 
    // for IRIX 5.3. Don't ask me why.
    cout << "";

    // Search for end of value
    while (!isspace(*cur) && (*cur != ',') && (*cur != '\0')) cur++;
    
    // Skip any spaces
    while (isspace(*cur) || (*cur == ',')) cur++;
  }

  return TRUE;
}

int
VolumeStatsArgs::getIntList(char *dst, char *key, char *nextArg)
{
  // Check for next argument
  if (nextArg == NULL) {
    cerr << key << " option requires an additional argument" << endl;
    exit(EXIT_FAILURE);
  }

  int *values = (int *) dst;

  // Set up pointers to end of string and first non-space character
  char *cur = nextArg;
  while (isspace(*cur)) cur++;
  unsigned n = 0;

  // Loop through string looking for space or comma-separated values
  while ((n < _MAX_ITEMS_IN_LIST) && (*cur != '\0')) {
    *values++ = atoi(cur);
    n++;

    // Without the following line, this function will choke when a 
    // value "0" is given on the command line. Don't ask me why.
    cout << "";

    // Search for end of value
    while (!isspace(*cur) && (*cur != ',') && (*cur != '\0')) cur++;
    
    // Skip any spaces
    while (isspace(*cur) || (*cur == ',')) cur++;
  }

  return TRUE;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_axis_order
@INPUT      : dst - Pointer to client data from argument table
              key - argument key
              nextArg - argument following key
@OUTPUT     : (nothing) 
@RETURNS    : TRUE or FALSE (so that ParseArgv will discard nextArg only
              when needed)
@DESCRIPTION: Routine called by ParseArgv to set the axis order
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : March 16, 1994 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int
VolumeStatsArgs::getAxisOrder(char *dst, char *key, char *nextArg)
     /* ARGSUSED */
{
   char **axis_order;
   char *cur;
   int ndims;

   /* Get pointer to client data */
   axis_order = (char **) dst;

   /* Check key */
   if (strcmp(key, "-transverse") == 0) {
      axis_order[0] = MIzspace;
      axis_order[1] = MIyspace;
      axis_order[2] = MIxspace;
      return FALSE;
   }
   if (strcmp(key, "-sagittal") == 0) {
      axis_order[0] = MIxspace;
      axis_order[1] = MIzspace;
      axis_order[2] = MIyspace;
      return FALSE;
   }
   if (strcmp(key, "-coronal") == 0) {
      axis_order[0] = MIyspace;
      axis_order[1] = MIzspace;
      axis_order[2] = MIxspace;
      return FALSE;
   }

   /* Make sure that we have a "-dimorder" argument */
   if (strcmp(key, "-dimorder") != 0) {
     cerr << "Unrecognized option \"" << key << "\": internal program error" << endl;
     exit(EXIT_FAILURE);
   }

   /* Check for next argument */
   if (nextArg == NULL) {
     cerr << "\"" << key << "\" option requires an additional argument" << endl;
     exit(EXIT_FAILURE);
   }

   /* Set up pointers to end of string and first non-space character */
   cur = nextArg;

   while (isspace(*cur)) cur++;
   ndims = 0;

   /* Loop through string looking for space or comma-separated names */
   while (*cur!='\0') {

     /* Get string */
     if (ndims < 3)
       axis_order[ndims] = cur;
     else {
       cerr << "You can specify up to 3 dimensions using -dimorder, not more" << endl;
       exit(EXIT_FAILURE);
     }

     /* Search for end of dimension name */
     while (!isspace(*cur) && (*cur != ',') && (*cur != '\0')) cur++;
     if (*cur != '\0') {
       *cur = '\0';
       cur++;
     }

     if (strcmp(axis_order[ndims], MIxspace) &&
	 strcmp(axis_order[ndims], MIyspace) &&
	 strcmp(axis_order[ndims], MIzspace)) {
       cerr << "Invalid spatial dimension " << axis_order[ndims] << " specified" << endl;
       exit(EXIT_FAILURE);
     }
     
     ndims++;
     
     /* Skip any spaces */
     while (isspace(*cur)) cur++;
     
   }

   return TRUE;
}
