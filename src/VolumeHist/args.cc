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
$RCSfile: args.cc,v $
$Revision: 1.2 $
$Author: stever $
$Date: 2003-11-17 04:30:56 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "args.h"
#include <version.h>

using namespace std;


int   args::clobber  = FALSE;
int   args::verbose  = TRUE;
int   args::number_of_bins   = 200;
int   args::number_of_classes   = -1;
int   args::selected_class = -1;
int   args::select_class_flag = FALSE;
int   args::limit_classes = FALSE;
int   args::auto_range = FALSE;
#ifdef HAVE_MATLAB
int   args::matlab_format = TRUE;
#else
int   args::matlab_format = FALSE;
#endif
int   args::mask_flag = FALSE;
int   args::window_flag = FALSE;

char *mask_string = NULL;

ArgvInfo args::argTable[] = {
  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "General options:"},
  {"-clobber", ARGV_CONSTANT, (char *)(int)TRUE, (char *) &args::clobber, 
   "Overwrite existing file."},
  {"-noclobber", ARGV_CONSTANT, (char *)(int)FALSE, (char *) &args::clobber, 
   "Do not overwrite existing file (default)."},
  {"-verbose", ARGV_CONSTANT, (char *)(int)TRUE, (char *) &args::verbose, 
   "Print out log messages as processing is being done (default)."},
  {"-quiet", ARGV_CONSTANT, (char *)(int)FALSE, (char *) &args::verbose, 
   "Do not print out any log messages."},
  {"-version", ARGV_FUNC, (char *) print_version_info, 
   (char *)MNI_LONG_VERSION, "Print out version info and exit."},

  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "\n Options:"},
  {"-bins", ARGV_INT, (char *) 1, (char *) &args::number_of_bins, 
   "Number of bins in each histogram."},
  {"-mask", ARGV_STRING, (char *) 1, (char *) &mask_string, 
   "<mask.mnc> Use mask volume to make one or more histograms."},
  {"-classes", ARGV_INT, (char *) 1, (char *) &args::number_of_classes, 
   "<number> Limit number of histograms to make."},
  {"-select", ARGV_INT, (char *) 1, (char *) &args::selected_class, 
   "<class> Only compute histogram for selected class."},
  {"-auto_range", ARGV_CONSTANT, (char *)(int)TRUE, (char *) &args::auto_range, 
   "Compute histogram range for each class."},
  {"-window", ARGV_CONSTANT, (char *)(int)TRUE, (char *) &args::window_flag, 
   "Use triangular Parzen window."},
#ifdef HAVE_MATLAB
  {"-matlab", ARGV_CONSTANT, (char *)(int)TRUE, (char *) &args::matlab_format, 
   "Produce files in matlab format. (default)"},
  {"-text", ARGV_CONSTANT, (char *)(int)FALSE, (char *) &args::matlab_format, 
   "Produce files in text format."},
#else
  {"-text", ARGV_CONSTANT, (char *)(int)FALSE, (char *) &args::matlab_format, 
   "Produce files in text format. (default)"},
  {"-matlab", ARGV_CONSTANT, (char *)(int)TRUE, (char *) &args::matlab_format, 
   "Produce files in matlab format. (not available)"},
#endif
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
@CREATED    : January 19, 1996
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
  if (ParseArgv(&argc, argv, argTable, 0) || (argc != 3)) {
    cerr << endl << "Usage: " << argv[0] 
	 << " [<options>] <infile> <outfile>" << endl;
    cerr <<         "       " << argv[0] << " [-help]" << endl << endl;
    exit(EXIT_FAILURE);
  }

  // check arguments
  if(number_of_bins < 1) 
    {
      cerr << "Must have one or bins per histogram\n";
      exit(EXIT_FAILURE);
    }

  if(number_of_classes > 0)
    limit_classes = TRUE;

  if(selected_class >= 0)
    {
      select_class_flag = TRUE;
    }
  
  // Check and expand paths
  MString extension;
  if(mask_string != NULL)
    {
      mask_flag = TRUE;

      maskPath = mask_string;
      if (!maskPath.expanded().existsCompressed(&extension)) {
        cerr << "Can't find " << maskPath << endl;
        exit(EXIT_FAILURE);
      }
      maskPath = maskPath.expanded() + extension;
    }
  else if(select_class_flag == TRUE)
    cerr << "Warning select option has no effect without mask option\n";


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




