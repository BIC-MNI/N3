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
$RCSfile: evaluateFieldArgs.cc,v $
$Revision: 1.1 $
$Author: bert $
$Date: 2003-04-16 14:25:59 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : evaluateFieldArgs.c
@INPUT      : 
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: argument parsing class for EvaluateField
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : July 9, 1996        John G. Sled
@MODIFIED   : 
@COPYRIGHT  : 1996
---------------------------------------------------------------------------- */

#include <config.h>
#include <EBTKS/Matrix.h>	// (bert) - Added EBTKS subdirectory
#include <EBTKS/MString.h>	// (bert)
#undef ROUND
#undef SIGN
extern "C" {
#include <volume_io.h>
}
#include <version.h>
#undef ROUND // Added to avoid conflict between Volume_io's and
             // AZgen's definition of ROUND, Alex Zijdenbos 97/12/05
#include "evaluateFieldArgs.h"

int evaluateArgs::clobber  = FALSE;
int evaluateArgs::verbose  = TRUE;
int evaluateArgs::use_mask = FALSE;
nc_type evaluateArgs::datatype = NC_UNSPECIFIED;
int evaluateArgs::signed_flag  = -MAXINT;

char *maskString = NULL;
char *likeString = NULL;

ArgvInfo evaluateArgs::ArgTable[] = {
  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "General options:"},
  {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &evaluateArgs::clobber, 
   "Overwrite existing file."},
  {"-noclobber", ARGV_CONSTANT, (char *) FALSE, (char *) &evaluateArgs::clobber, 
   "Do not overwrite existing file (default)."},
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &evaluateArgs::verbose, 
   "Print out log messages as processing is being done (default)."},
  {"-quiet", ARGV_CONSTANT, (char *) FALSE, (char *) &evaluateArgs::verbose, 
   "Do not print out any log messages."},
  {"-version", ARGV_FUNC, (char *) print_version_info, 
   (char *)MNI_LONG_VERSION, "Print out version info and exit."},
  
  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "\n Output options:"},
  {"-like", ARGV_STRING, (char *) 1, (char *) &likeString, 
   "Specify volume that field is to be like."},
  {"-mask", ARGV_STRING, (char *) 1, (char *) &maskString, 
   "Specify region on which field is evaluated."},
  {"-byte", ARGV_CONSTANT, (char *) NC_BYTE, (char *) &evaluateArgs::datatype,
   "Write out byte data."},
  {"-short", ARGV_CONSTANT, (char *) NC_SHORT, (char *) &evaluateArgs::datatype,
   "Write out short integer data."},
  {"-long", ARGV_CONSTANT, (char *) NC_LONG, (char *) &evaluateArgs::datatype,
   "Write out long integer data."},
  {"-float", ARGV_CONSTANT, (char *) NC_FLOAT, (char *) &evaluateArgs::datatype,
   "Write out single-precision floating-point data."},
  {"-double", ARGV_CONSTANT, (char *) NC_DOUBLE, (char *) &evaluateArgs::datatype,
   "Write out double-precision floating-point data."},
  {"-signed", ARGV_CONSTANT, (char *) TRUE, (char *) &evaluateArgs::signed_flag,
   "Write signed integer data."},
  {"-unsigned", ARGV_CONSTANT, (char *) FALSE, (char *) &evaluateArgs::signed_flag,
   "Write unsigned integer data."},
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

evaluateArgs::evaluateArgs(int argc, char **argv)
: command(argv[0])
{
  for (unsigned i = 1; i < argc; i++) {
    command += " ";
    command += argv[i];
  }

  set_program_name(argv[0]);  // for version info

  // Call ParseArgv
  if (ParseArgv(&argc, argv, ArgTable, 0) || argc != 3) {
    cerr << endl << "Usage: " << argv[0] 
         << " [<options>] -like <volume.mnc> <infile.fld> <outfile.mnc>\n";
    cerr <<         "       " << argv[0] << " [-help]" << endl << endl;
    exit(EXIT_FAILURE);
  }

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

  if(likeString == NULL)
    {
      cerr << "Must specify a -like volume\n";
      exit(EXIT_FAILURE);
    }
  else
    {
      likePath = likeString;
      if (!likePath.expanded().existsCompressed(&extension)) {
	cerr << "Can't find " << likePath << endl;
	exit(EXIT_FAILURE);
      }
      likePath = likePath.expanded() + extension;
    }

  inputPath = argv[1];
  if (!inputPath.expanded().existsCompressed(&extension)) {
    cerr << "Can't find " << inputPath << endl;
    exit(EXIT_FAILURE);
  }
  inputPath = inputPath.expanded() + extension;

  outputPath = argv[2];
  if (!clobber && outputPath.expanded().existsCompressed(&extension)) {
    cerr << outputPath << " exists! (use -clobber to overwrite)" << endl;
    exit(EXIT_FAILURE);
  }
  outputPath = outputPath.expanded().removeCompressedExtension();
}







