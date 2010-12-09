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
$RCSfile: splineSmooth.cc,v $
$Revision: 1.5 $
$Author: claude $
$Date: 2010-12-09 19:35:01 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : splineSmooth.c,v
@INPUT      : 
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: Tool for smoothing and extrapolating data in minc volumes
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 21, 1996 (John G. Sled)
@MODIFIED   : Log: splineSmooth.c,v 
 * Revision 1.2  1996/04/23  13:36:58  jgsled
 * Working version.  Problems with thin plate splines have been fixed.
 *
 * Revision 1.1  1996/04/21  23:41:50  jgsled
 * Initial version of SplineSmooth tool
 * - B spline implementation appears to work
 *
@COPYRIGHT  : 1996
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <iostream>		// (bert)
using namespace std;		// (bert)
#include <math.h>
#include "../Splines/Spline.h"
#include <EBTKS/Matrix.h>	// (bert) 
#include "../Splines/TBSpline.h"
#include <EBTKS/MString.h>	// (bert)
#undef ROUND
#undef SIGN
extern "C" {
#include <volume_io.h>
}
#undef ROUND // Added to avoid conflict between Volume_io's and
             // AZgen's definition of ROUND, Alex Zijdenbos 97/12/05

#include "splineSmoothArgs.h"
#include "fieldIO.h"

//#define DEBUG_SPLINESMOOTH

//---------------------------------------------------------------------------------
//  Implementation notes
/*
  The spline basis functions are defined in a world coordinate system aligned
  with the voxel coordinate system and sharing the same origin.

 */

//---------------------------------------------------------------------------------
// Declarations
DblMat volume_domain(Volume volume);
DblMat reduced_domain(Volume volume);
Spline *createThinPlateSpline(const DblMat &domain, double distance,
			      double lambda, int verbose);
void fitSplinesToVolume(Spline *spline, Volume volume, int subsample);
void fitSplinesToVolume(Spline *spline, Volume volume, Volume mask_volume,
			const DblMat &domain, int subsample);
void fitSplinesToVolumeLookup(TBSplineVolume *spline, Volume volume,
                              int subsample);
void fitSplinesToVolumeLookup(TBSplineVolume *spline, Volume volume, 
                              Volume mask_volume, const DblMat &domain,
                              int subsample);

//--------------------------------------------------------------------------------
// main program
int main( int argc,  char *argv[] )
     
{
  Volume volume = NULL, mask_volume = NULL, output_mask_volume = NULL;     
  nc_type output_type; /* type for output volume */
  args args(argc, argv);
  DblMat domain;   // region in world coordinates on which splines are defined

  /* use slice based volume caching to reduce memory usage */
  set_n_bytes_cache_threshold(0);           // Always cache volume
  set_cache_block_sizes_hint(SLICE_ACCESS); // Cache volume slices
  set_default_max_bytes_in_cache(0);        // keep only one slice at a time

  // load data volume as float type
  volume = loadFloatVolume(args.inputPath, &output_type);

  if(args.use_mask == TRUE)     
  { // open mask volume as well
    mask_volume = loadVolume(args.maskPath);
    if(compareVolumes(volume, mask_volume) == FALSE)
      {
	cerr << "Mask volume and input volume must be the same size.\n";
	return(2);
      }
  }

  if(args.use_output_mask == TRUE)     
  { // open output mask volume as well
    output_mask_volume = loadVolume(args.outputMaskPath);
    if(compareVolumes(volume, output_mask_volume) == FALSE)
      {
	cerr << "Output mask volume and input volume must be the same size.\n";
	return(2);
      }
  }

  // determine domain on whick splines are defined
  if((args.extrapolate == FALSE || args.spline == thin_plate_spline)
     && args.use_mask == TRUE && args.full_support == FALSE)
    // shrink bounding box around mask
    domain = reduced_domain(mask_volume);
  else
    // domain is whole volume
    domain = volume_domain(volume);


  // create spline basis
  Spline *theSplines;
  
  if(args::spline == b_spline) // create B spline basis
    {
      Real separations[N_DIMENSIONS];
      int sizes[N_DIMENSIONS];
      Real start[N_DIMENSIONS] = { 0.0, 0.0, 0.0 };
      get_volume_separations(volume, separations);
      get_volume_sizes(volume, sizes);
      theSplines = new TBSplineVolume(domain, start, separations, sizes,
                                      args.distance, args.lambda);
    }
  else
    theSplines = createThinPlateSpline(domain, args.distance, args.lambda,
                                       args.verbose);
  
  // do least squares fit to data 
  if(args::spline == b_spline) {
    if(args.use_mask)
      fitSplinesToVolumeLookup((TBSplineVolume *)theSplines, volume,
                               mask_volume, domain, args.subsample);
    else
      fitSplinesToVolumeLookup((TBSplineVolume *)theSplines,
                               volume, args.subsample);
  }
  else {
    if(args.use_mask)
      fitSplinesToVolume(theSplines, volume, mask_volume, 
                         domain, args.subsample);
    else
      fitSplinesToVolume(theSplines, volume, args.subsample);
  }

#ifdef DEBUG_SPLINESMOOTH
  theSplines->saveState("spline.mat");
#endif
    
  if(args.produce_volume == TRUE) 
    {
      // write smooth function to volume
      Real real_min, real_max;
      if(args::spline == b_spline)
      {
        if(args.use_output_mask == TRUE && args.extrapolate == FALSE) {
          if(!args.full_support && args.use_mask && args.verbose) 
            printf("Warning: if output mask is larger than input mask, the"
                   " -full_support\nis necessary to avoid clipping"
                   " of the result.\n");
          smoothVolumeLookup((TBSplineVolume *) theSplines, volume,
                             output_mask_volume, &real_min, &real_max);
        }
        else if(args.use_mask == TRUE && args.extrapolate == FALSE)
          smoothVolumeLookup((TBSplineVolume *) theSplines, volume,
                             mask_volume, &real_min, &real_max);
        else
          smoothVolumeLookup((TBSplineVolume *)theSplines, volume,
                             &real_min, &real_max);
      }
      else   // this plate spline
        {
          if(args.use_output_mask == TRUE && args.extrapolate == FALSE)
            smoothVolume(theSplines, volume, output_mask_volume, 
                         &real_min, &real_max);
          else if(args.use_mask == TRUE && args.extrapolate == FALSE)
            smoothVolume(theSplines, volume, mask_volume,
                        &real_min, &real_max);
          else
            smoothVolume(theSplines, volume, &real_min, &real_max);
        }  
      // write smooth volume to disk in original data type 
      outputVolume(volume, args.outputPath, output_type, TRUE, real_min, real_max, 
                   args.command);
    }
  if(args.produce_compact == TRUE)
    {
      outputCompactField(args.compactPath, domain, args.distance,
                         theSplines->getCoefficients(), args.spline,
                         args.command, volume);
    }

  delete theSplines;
  delete_volume( volume );
  if( mask_volume ) delete_volume( mask_volume );
  if( output_mask_volume ) delete_volume( output_mask_volume );

  return(0);
} 


//-----------------------------------------------------------------------------
// Supporting functions



// determine domain from size of volume in world coordinates
// Returns an N_DIMENSIONS by 2 matrix
DblMat
volume_domain(Volume volume)
{
  int sizes[N_DIMENSIONS];
  Real separations[N_DIMENSIONS];
  DblMat domain(N_DIMENSIONS,2);

  get_volume_separations(volume, separations);
  get_volume_sizes(volume, sizes);
  
  for(int i = 0; i < N_DIMENSIONS; i++)
    {
      if(separations[i] > 0) {
        domain(i,0) = -0.5*separations[i];
        domain(i,1) = (sizes[i]-0.5)*separations[i];
      }
      else {
        domain(i,1) = -0.5*separations[i];
        domain(i,0) = (sizes[i]-0.5)*separations[i];
      }
    }
  return domain;
}

// determine smallest domain that contains all non-zero values in volume
DblMat reduced_domain(Volume volume)
{
  int sizes[N_DIMENSIONS];
  Real separations[N_DIMENSIONS];
  DblMat domain(N_DIMENSIONS,2);
  int i,j,k;

  get_volume_separations(volume, separations);
  get_volume_sizes(volume, sizes);

  // limits
  int lower[N_DIMENSIONS], upper[N_DIMENSIONS];
  for(i = 0; i < N_DIMENSIONS; i++)
    {
      lower[i] = sizes[i]-1;
      upper[i] = 0;
    }

  // proceed through volume slicewise starting at 0,0,0
  for(i = 0; i < sizes[0]; i++)
    for(j = 0; j < sizes[1]; j++)
      for(k = 0; k < sizes[2]; k++)
	if(get_volume_real_value(volume, i, j, k, 0, 0) > 0)
	  {
	    if(i > upper[0]) upper[0] = i;
	    if(j > upper[1]) upper[1] = j;
	    if(k > upper[2]) upper[2] = k;
	    if(i < lower[0]) lower[0] = i;
	    if(j < lower[1]) lower[1] = j;
	    if(k < lower[2]) lower[2] = k;
	  }

  if(lower[0] > upper[0])
    {
      cerr << "Mask is empty, no processing done.\n";
      exit(0);
    }

  for(i = 0; i < N_DIMENSIONS; i++)
    {
      if(separations[i] > 0) {
        domain(i,0) = separations[i]*(lower[i]-0.5);
        domain(i,1) = separations[i]*(upper[i]+0.5);
      }
      else {
        domain(i,1) = separations[i]*(lower[i]-0.5);
        domain(i,0) = separations[i]*(upper[i]+0.5);
      }
    }
  return domain;
}  


// do fit on all data in volume
void
fitSplinesToVolume(Spline *spline, Volume volume, int subsample = 1)
{
  int sizes[N_DIMENSIONS];
  Real separations[N_DIMENSIONS];
  int i,j,k;
  float point[N_DIMENSIONS];
  Real value;

  get_volume_separations(volume, separations);
  get_volume_sizes(volume, sizes);

  progress_struct progress;
  initialize_progress_report(&progress, FALSE, sizes[0], 
			    "Fitting splines");

  for(i = 0; i < sizes[0]; i += subsample)
    {
      for(j = 0; j < sizes[1]; j += subsample)
	for(k = 0; k < sizes[2]; k += subsample)
	  {
	    point[0] = i*separations[0];
	    point[1] = j*separations[1];
	    point[2] = k*separations[2];

	    value = get_volume_real_value(volume, i, j, k, 0, 0);
	    spline->addDataPoint(point, value);
	  }
      update_progress_report(&progress, i+1);
    }
  terminate_progress_report(&progress);

  if(spline->fit() == FALSE) // fit splines to the data
    {
      cerr << "Fatal Error: Spline fit failed.\n"
        "A possible cause of this problem is an empty of nearly\n"
        "empty mask volume.\n";
      exit(3);
    }
}

// fit splines to data points that are non-zero in the mask
void
fitSplinesToVolume(Spline *spline, Volume volume, Volume mask_volume, 
		   const DblMat &domain, int subsample = 1)
{
  int sizes[N_DIMENSIONS];
  Real separations[N_DIMENSIONS];
  int i,j,k;
  float point[N_DIMENSIONS];
  Real value;

  get_volume_separations(volume, separations);
  get_volume_sizes(volume, sizes);

  // only look at values within domain
  int lower[N_DIMENSIONS], upper[N_DIMENSIONS];
  for(i = 0; i < N_DIMENSIONS; i++)
    {
      lower[i] = (int) ceil(domain(i,0)/separations[i]);
      upper[i] = (int) floor(domain(i,1)/separations[i]);
    }

  progress_struct progress;
  initialize_progress_report(&progress, FALSE, upper[0]-lower[0]+1, 
			    "Fitting splines");
  for(i = lower[0]; i <= upper[0]; i += subsample)
    {
      for(j = lower[1]; j <= upper[1]; j += subsample)
	for(k = lower[2]; k <= upper[2]; k += subsample)
	  {
	    if(get_volume_real_value(mask_volume, i, j, k, 0, 0) > 0.5)
	      {
		point[0] = i*separations[0];
		point[1] = j*separations[1];
		point[2] = k*separations[2];

		value = get_volume_real_value(volume, i, j, k, 0, 0);
		spline->addDataPoint(point, value);
	      }
	  }
      update_progress_report(&progress, i-lower[0]+1);
    }

  terminate_progress_report(&progress);
  if(spline->fit() == FALSE) // fit splines to the data
    {
      cerr << "Fatal Error: Spline fit failed.\n";
      exit(3);
    }
}

// do fit on all data in volume
void
fitSplinesToVolumeLookup(TBSplineVolume *spline, Volume volume,
                         int subsample = 1)
{
  int sizes[N_DIMENSIONS];
  Real separations[N_DIMENSIONS];
  int i,j,k;
  Real value;

  get_volume_separations(volume, separations);
  get_volume_sizes(volume, sizes);

  progress_struct progress;
  initialize_progress_report(&progress, FALSE, sizes[0], 
			    "Fitting splines");

  for(i = 0; i < sizes[0]; i += subsample)
    {
      for(j = 0; j < sizes[1]; j += subsample)
	for(k = 0; k < sizes[2]; k += subsample)
	  {
	    value = get_volume_real_value(volume, i, j, k, 0, 0);
	    spline->addDataPoint(i,j,k, value);
	  }
      update_progress_report(&progress, i+1);
    }
  terminate_progress_report(&progress);

  if(spline->fit() == FALSE) // fit splines to the data
    {
      cerr << "Fatal Error: Spline fit failed.\n";
      exit(3);
    }
}

// fit splines to data points that are non-zero in the mask
void
fitSplinesToVolumeLookup(TBSplineVolume *spline, Volume volume,
                         Volume mask_volume, 
                         const DblMat &domain, int subsample = 1)
{
  int sizes[N_DIMENSIONS];
  Real separations[N_DIMENSIONS];
  int i,j,k;
  Real value;

  get_volume_separations(volume, separations);
  get_volume_sizes(volume, sizes);

  // only look at values within domain
  int lower[N_DIMENSIONS], upper[N_DIMENSIONS];
  for(i = 0; i < N_DIMENSIONS; i++)
    {
      if(separations[i] > 0) {
        lower[i] = (int) ceil(domain(i,0)/separations[i]);
        upper[i] = (int) floor(domain(i,1)/separations[i]);
      }
      else {
        upper[i] = (int) floor(domain(i,0)/separations[i]);
        lower[i] = (int) ceil(domain(i,1)/separations[i]);
      }
    }

  progress_struct progress;
  initialize_progress_report(&progress, FALSE, upper[0]-lower[0]+1, 
			    "Fitting splines");
  for(i = lower[0]; i <= upper[0]; i += subsample)
    {
      for(j = lower[1]; j <= upper[1]; j += subsample)
	for(k = lower[2]; k <= upper[2]; k += subsample)
	  {
	    if(get_volume_real_value(mask_volume, i, j, k, 0, 0) > 0.5)
	      {
		value = get_volume_real_value(volume, i, j, k, 0, 0);
		spline->addDataPoint(i,j,k, value);
	      }
	  }
      update_progress_report(&progress, i-lower[0]+1);
    }

  terminate_progress_report(&progress);
  if(spline->fit() == FALSE) // fit splines to the data
    {
      cerr << "Fatal Error: Spline fit failed.\n";
      exit(3);
    }
}










