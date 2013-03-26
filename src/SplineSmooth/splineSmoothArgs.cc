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
$RCSfile: splineSmoothArgs.cc,v $
$Revision: 1.2 $
$Author: stever $
$Date: 2003-11-17 04:30:56 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : splineSmoothArgs.c,v
@INPUT      : 
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: argument parsing class for SplineSmooth
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 21, 1996 (John G. Sled)
@MODIFIED   : Log: splineSmoothArgs.c,v 
 * Revision 1.2  1996/04/23  13:36:59  jgsled
 * Working version.  Problems with thin plate splines have been fixed.
 *
 * Revision 1.1  1996/04/21  23:41:50  jgsled
 * Initial version of SplineSmooth tool
 * - B spline implementation appears to work
 *
@COPYRIGHT  : 1996
---------------------------------------------------------------------------- */

#include <config.h>
#include <EBTKS/Matrix.h>	// (bert) - Added EBTKS subdirectory
#include <EBTKS/MString.h>	// (bert)
#undef VIO_ROUND
#undef SIGN
extern "C" {
#include <volume_io.h>
}
#include <version.h>

#undef VIO_ROUND // Added to avoid conflict between Volume_io's and
             // AZgen's definition of VIO_ROUND, Alex Zijdenbos 97/12/05
#include "splineSmoothArgs.h"


using namespace std;


int      args::clobber  = FALSE;
int      args::verbose  = TRUE;
double   args::lambda = 0.01;  // chosen to provide minimal smoothing except in 
                               //  regions of missing data 
double   args::distance = 60;  // distance in mm
enum spline_type args::spline = b_spline;
int      args::use_mask = FALSE;
int      args::use_output_mask = FALSE;
int      args::extrapolate = FALSE;
int      args::full_support = FALSE;
int      args::subsample = 1;  // no subsampling
int      args::produce_volume = 1;  // produce minc volume
int      args::produce_compact = 0; // don't produce coefficient representation

char *maskString = NULL;
char *outputMaskString = NULL;
char *compactString = NULL;

ArgvInfo args::argTable[] = {
  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "General options:"},
  {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &args::clobber, 
   "Overwrite existing file."},
  {"-noclobber", ARGV_CONSTANT, (char *) FALSE, (char *) &args::clobber, 
   "Do not overwrite existing file (default)."},
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &args::verbose, 
   "Print out log messages as processing is being done (default)."},
  {"-quiet", ARGV_CONSTANT, (char *) FALSE, (char *) &args::verbose, 
   "Do not print out any log messages."},
  {"-version", ARGV_FUNC, (char *) print_version_info, 
   (char *)MNI_LONG_VERSION, "Print out version info and exit."},

   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "\n Options:"},
  {"-lambda", ARGV_FLOAT, (char *) 1, (char *) &args::lambda, 
   "Scale invariant smoothing parameter.  The default is chosen\n to provide minimal"
   " smoothing except in regions of missing data."},
  {"-distance", ARGV_FLOAT, (char *) 1, (char *) &args::distance, 
   "Distance between basis functions in mm.  This parameter\n determines the overall"
   " smoothness."},
  {"-b_spline", ARGV_CONSTANT, (char *) b_spline, (char *) &args::spline,
   "Basis functions are tensor B splines (default)."},
  {"-tp_spline", ARGV_CONSTANT, (char *) thin_plate_spline,
   (char *) &args::spline,
   "Basis functions are thin plate splines."},
  {"-mask", ARGV_STRING, (char *) 1, (char *) &maskString, 
   "Specify data region to which splines are fit."},
  {"-output_mask", ARGV_STRING, (char *) 1, (char *) &outputMaskString, 
   "Specify region on which to evaluate smoothed function."},
  {"-extrapolate", ARGV_CONSTANT, (char *) TRUE, (char *) &args::extrapolate, 
   "Extrapolate to region outside mask (implies full_support)."},
  {"-noextrapolate", ARGV_CONSTANT, (char *)FALSE, (char *)&args::extrapolate, 
   "Do not extrapolate to region outside mask. (default)"},
  {"-full_support", ARGV_CONSTANT, (char *) TRUE, (char *) &args::full_support,   "Cover region outside mask with basis functions."},
  {"-subsample", ARGV_INT, (char *) 1, (char *) &args::subsample,
   "Only fit to every nth data point."},
  {"-compact", ARGV_STRING, (char *) 1, (char *) &compactString,
   "Specify file in which to store compact field representation."},
  {"-novolume", ARGV_CONSTANT, (char *) FALSE, (char *) &args::produce_volume,
   "Do not produce smoothed minc volume"},
  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   " "},

  {NULL, ARGV_END, NULL, NULL, NULL}
};
   
/* ----------------------------- MNI Header -----------------------------------
@NAME       : Constructor of args
@INPUT      : argc - number of command-line arguments
              argv - command-line arguments
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Builds a args object from commandline arguments
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 21, 1996
@MODIFIED   : 
---------------------------------------------------------------------------- */

args::args(int argc, char **argv)
: command(argv[0])
{
  for (unsigned i = 1; i < argc; i++) {
    command += " ";
    command += argv[i];
  }
  set_program_name(argv[0]);  // for version info

  // Call ParseArgv
  if (ParseArgv(&argc, argv, argTable, 0) || 
      (argc != 3 && produce_volume == TRUE) || 
      (argc != 2 && produce_volume == FALSE)) {
    cerr << endl << "Usage: " << argv[0] 
	 << " [<options>] <infile> <outfile>" << endl;
    cerr <<         "       " << argv[0] << " [-help]" << endl << endl;
    exit(EXIT_FAILURE);
  }

  // check arguments
  if(distance <=  0) 
    {
      cerr << "Distance parameter must be positive\n";
      exit(EXIT_FAILURE);
    }

  if(lambda < 0 || (spline == b_spline && lambda == 0) )
    {
      cerr << "Smoothing parameter lambda must be positive\n";
      exit(EXIT_FAILURE);
    }

  if(subsample <= 0)
    {
      cerr << "Subsampling factor must be positive.\n";
      exit(EXIT_FAILURE);
    }

  // Check and expand paths
  MString extension;
  if(maskString != NULL)
    {
      use_mask = TRUE;
      maskPath = maskString;
      if (!maskPath.expanded().existsCompressed(&extension)) {
	cerr << "Can't find " << maskPath << endl;
	exit(EXIT_FAILURE);
      }
      maskPath = maskPath.expanded() + extension;
    }

  if(compactString != NULL)
    {
      produce_compact = TRUE;
      compactPath = compactString;
      if (!clobber && compactPath.expanded().existsCompressed(&extension)) {
        cerr << compactPath << " exists! (use -clobber to overwrite)" << endl;
        exit(EXIT_FAILURE);
      }   
      compactPath = compactPath.expanded() + extension;
    }

  if(outputMaskString != NULL)
    {
      use_output_mask = TRUE;
      outputMaskPath = outputMaskString;
      if (!outputMaskPath.expanded().existsCompressed(&extension)) {
	cerr << "Can't find " << outputMaskPath << endl;
	exit(EXIT_FAILURE);
      }
      outputMaskPath = outputMaskPath.expanded() + extension;
    }

  inputPath = argv[1];
  if (!inputPath.expanded().existsCompressed(&extension)) {
    cerr << "Can't find " << inputPath << endl;
    exit(EXIT_FAILURE);
  }
  inputPath = inputPath.expanded() + extension;

  if(argc > 2) {
    outputPath = argv[2];
    if (!clobber && outputPath.expanded().existsCompressed(&extension)) {
      cerr << outputPath << " exists! (use -clobber to overwrite)" << endl;
      exit(EXIT_FAILURE);
    }
    outputPath = outputPath.expanded().removeCompressedExtension();
  }

  if(produce_volume == FALSE && produce_compact == FALSE)
    {
      cerr << "Nothing to do.  Specify an output format.\n";
      exit(EXIT_FAILURE);
    }
}





