#include <config.h>
#include <limits.h>
#include "args.h"
#include <version.h>

using namespace std;


int    args::clobber  = FALSE;
int    args::verbose  = TRUE;
double args::fwhm = 1.0;
double args::noise = 0.1;
double args::range[2] = { -DBL_MAX, DBL_MAX };
int    args::given_min = FALSE;
int    args::given_max = FALSE;
int    args::debug_flag = FALSE;
int    args::blur_flag = FALSE;

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
   "\n options:"},
  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &args::fwhm, 
   "Width of deconvolution kernel."},
  {"-noise", ARGV_FLOAT, (char *) 1, (char *) &args::noise, 
   "Noise level for Weiner filter."},
  {"-blur", ARGV_CONSTANT, (char *)(int)TRUE, (char *) &args::blur_flag, 
   "Skip deblurring step."},
  {"-min", ARGV_FLOAT, (char *) NULL, (char *) &args::range[0], 
   "Bin value mapped to zero in lookup table."},
  {"-max", ARGV_FLOAT, (char *) NULL, (char *) &args::range[1], 
   "Bin value mapped to one in lookup table."},
  {"-range", ARGV_FLOAT, (char *) 2, (char *) args::range, 
   "Range of bin values mapped from zero to one in lookup table."},
  {"-debug", ARGV_CONSTANT, (char *)(int)TRUE, (char *) &args::debug_flag, 
   "Save intermediate steps as matlab files."},

  {NULL, ARGV_END, NULL, NULL, NULL}
};
   
/* ----------------------------- MNI Header -----------------------------------
@NAME       : Constructor of args
@INPUT      : argc - number of command-line arguments
              argv - command-line arguments
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Builds an args object from command line arguments
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Febuary 27, 1996
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
  if (ParseArgv(&argc, argv, argTable, 0 ) || (argc != 3)) {
    cerr << endl << "Usage: " << argv[0] 
	 << " [<options>] <infile> <outfile>" << endl;
    cerr <<         "       " << argv[0] << " -help" << endl << endl;
    exit(EXIT_FAILURE);
  }

  // check arguments
  if(fwhm  <= 0) 
    {
      cerr << "Error: FWHM must be postive.\n";
      exit(EXIT_FAILURE);
    }

  if(noise < 0)
    {
      cerr << "Error: Noise threshold must be postive.\n";
      exit(EXIT_FAILURE);
    }
  if(range[0] > -DBL_MAX) given_min = TRUE;
  if(range[1] < DBL_MAX) given_min = TRUE;
 
  // Check and expand paths
  MString extension;

  inputPath = argv[1];
  if (!inputPath.expanded().existsCompressed(&extension)) {
    cerr << "Can't find " << inputPath << endl;
    exit(EXIT_FAILURE);
  }
  inputPath = inputPath.expanded() + extension;

  outputPath = argv[2];
  if (!clobber && outputPath.expanded().existsCompressed(&extension)) {
    cerr << outputPath << " exists. (use -clobber to overwrite)" << endl;
    exit(EXIT_FAILURE);
  }
  outputPath = outputPath.expanded().removeCompressedExtension();
}






