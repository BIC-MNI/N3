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
$RCSfile: evaluateField.cc,v $
$Revision: 1.1 $
$Author: bert $
$Date: 2003-04-16 14:25:59 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : evaluateField.c,v
@INPUT      : 
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: Tool for evaluating compact spline representations as minc volumes
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : July 9, 1996      John G. Sled 
@MODIFIED   : Log: evaluateField.c,v 
@COPYRIGHT  : 1996
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/INSECT/N3/src/EvaluateField/evaluateField.cc,v 1.1 2003-04-16 14:25:59 bert Exp $";
#endif


#include "EBTKS/Spline.h"	/* (bert) - Added EBTKS subdirectory */
#include "EBTKS/TBSpline.h"	/* (bert) */
#include "EBTKS/MString.h"	/* (bert) */
#undef ROUND
#undef SIGN
extern "C" {
#include <volume_io.h>
}
#undef ROUND // Added to avoid conflict between Volume_io's and
             // AZgen's definition of ROUND, Alex Zijdenbos 97/12/05
#include "evaluateFieldArgs.h"
#include "../SplineSmooth/fieldIO.h"
#include <iostream>		/* (bert) */
using namespace std;		/* (bert) */

//#define DEBUG_SPLINESMOOTH

#define N_DIMENSIONS 3


//----------------------------------------------------------------------------
//  Implementation notes
/*
  The spline basis functions are defined in a world coordinate system aligned
  with the voxel coordinate system and sharing the same origin.

 */

//----------------------------------------------------------------------------
// main program
int main( int argc,  char *argv[] )
     
{
  Volume volume, mask_volume; 
  enum spline_type spline_type;   
  nc_type output_type; /* type for output volume */
  BOOLEAN signed_flag;
  evaluateArgs args(argc, argv);
  DblMat domain;   // region in world coordinates on which splines are defined

  /* use slice based volume caching to reduce memory usage */
  set_n_bytes_cache_threshold(0);           // Always cache volume
  set_cache_block_sizes_hint(SLICE_ACCESS); // Cache volume slices
  set_default_max_bytes_in_cache(0);        // keep only one slice at a time

  // load data volume as float type
  volume = loadEmptyFloatVolume(args.likePath, &output_type, &signed_flag);

  if(args.use_mask == TRUE)     
  { // open mask volume as well
    mask_volume = loadVolume(args.maskPath);
    if(compareVolumes(volume, mask_volume) == FALSE)
      {
	cerr << "Mask volume and input volume must be the same size.\n";
	return(2);
      }
  }

  // read compact field representation from file
  Spline *theSplines;
  Status status = inputCompactField(args.inputPath, &theSplines, 
                                    &spline_type, volume);
  if(status != OK) {
    cerr << "Failure reading field file: " << args.inputPath << endl;
    return(ERROR);
  }

  // write smooth function to volume
  Real real_min, real_max;
  
  //#ifdef DEBUG_FIELDIO
  double value = (*(TBSplineVolume *) theSplines)(0,0,0);
  float point[3] = {0.0, 0.0, 0.0};
  value = (*theSplines)(point);
  //#endif
  if(args.use_mask == TRUE)
    {
      if(spline_type == b_spline)
        smoothVolumeLookup((TBSplineVolume *) theSplines, volume, mask_volume,
                           &real_min, &real_max);
      else
        smoothVolume(theSplines, volume, mask_volume,
                     &real_min, &real_max);
    }
  else
    {
      if(spline_type == b_spline)
        smoothVolumeLookup((TBSplineVolume *) theSplines, volume,
                           &real_min, &real_max);
      else
        smoothVolume(theSplines, volume, &real_min, &real_max); 
    }

  if (args.datatype != NC_UNSPECIFIED)
    output_type = args.datatype;

  if (args.signed_flag != -MAXINT)
    signed_flag = args.signed_flag;
  
  // write smooth volume to disk in original data type 
  outputVolume(volume, args.outputPath, output_type, signed_flag, real_min, real_max, 
               args.command);

  return(0);
} 

