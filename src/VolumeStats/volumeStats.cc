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
$RCSfile: volumeStats.cc,v $
$Revision: 1.2 $
$Author: stever $
$Date: 2003-11-17 04:30:56 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : volumeStats
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : error status
@DESCRIPTION: Program that spits out a bunch of stats for a given volume
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 15, 1995 (Alex Zijdenbos)
@MODIFIED   : $Log: volumeStats.cc,v $
@MODIFIED   : Revision 1.2  2003-11-17 04:30:56  stever
@MODIFIED   : *** empty log message ***
@MODIFIED   :
@MODIFIED   : Revision 1.1  2003/04/16 14:32:14  bert
@MODIFIED   : Initial checkin
@MODIFIED   :
@MODIFIED   : Revision 1.11  1999/03/16 05:21:54  alex
@MODIFIED   : Fixed handling of volumes with more than 3 dimensions
@MODIFIED   :
@MODIFIED   : Revision 1.10  1997/10/06 15:15:52  alex
@MODIFIED   : Added -debug options and dimension checks
@MODIFIED   :
@MODIFIED   : Revision 1.9  1997/09/23 21:55:23  alex
@MODIFIED   : Fixed return value bug in getList functions and removed ROUND from CoM print
@MODIFIED   :
@MODIFIED   : Revision 1.8  1997/08/22 17:08:06  alex
@MODIFIED   : Added check for a complete lack of valid voxels
@MODIFIED   :
@MODIFIED   : Revision 1.7  1997/07/30 19:10:08  alex
@MODIFIED   : Added pctT and dimension ordering options
@MODIFIED   :
@MODIFIED   : Revision 1.6  1997/02/07 15:26:20  alex
@MODIFIED   : Added -slice option
@MODIFIED   :
@MODIFIED   : Revision 1.5  1997/01/13 20:38:29  alex
@MODIFIED   : Added list facility to -binvalue, added -volume switch
@MODIFIED   :
@MODIFIED   : Revision 1.4  1996/12/19 20:43:57  alex
@MODIFIED   : Various changes directed at g++ compatibility
@MODIFIED   :
@MODIFIED   : Revision 1.3  1996/11/28 21:55:21  alex
@MODIFIED   : Added list of extrema handling
@MODIFIED   :
@MODIFIED   : Revision 1.2  1996/10/12 02:29:22  alex
@MODIFIED   : extended -maskfloor and -maskceil to accept a list of values
@MODIFIED   :
 * Revision 1.1.1.1  1996/08/29  19:12:23  alex
 * Source for volume_stats
 *
@COPYRIGHT  :
---------------------------------------------------------------------------- */

#include <assert.h>
#include <time.h>
#include <EBTKS/Minc.h>		// (bert) - Added EBTKS subdirectory
#include <EBTKS/CachedArray.h>	// (bert)
#include <EBTKS/Histogram.h>	// (bert)
#include "VolumeStatsArgs.h"

using namespace std;


Volume loadVolume(const Path& path, char **axisOrder = 0, int verbose = 0);

//
// Main program
//

int
main(int argc, char *argv[])
{
  // Get argument information
  VolumeStatsArgs args(argc, argv);

  Volume mask = 0;
  unsigned maskD1, maskD2, maskD3;
  IntArray sizes(3);
  // Load mask, if specified
  if ((char *)args.maskPath) {
    mask = loadVolume(args.maskPath, args.axisOrder, args.verbose);
    if (mask) {
      get_volume_sizes(mask, sizes.contents());
      maskD1 = sizes[0];
      maskD2 = sizes[1];
      maskD3 = sizes[2];
    }

    if (args.useWorldCoord) {
      IntArray blockSizes(sizes);
      blockSizes[0] = blockSizes[1] = blockSizes[2] = 8;
      set_volume_cache_block_sizes(mask, blockSizes);
    }
  }

  unsigned nExtremaRuns     = MAX(args.nExtrema, 1);
  unsigned nMaskExtremaRuns = MAX(args.nMaskExtrema, 1);

  SimpleArray<float> *voxels = 0;
  if (args.all || args.calcMedian || args.calcMajority || args.calcBiModalT || 
      args.calcPctT || args.calcEntropy) {
    if (args.cached) 
      voxels = new CachedArray<float>(0);
    else
      voxels = new SimpleArray<float>(0);
    assert(voxels);
  }

  if (args.debug)
    Array<float>::debug(TRUE);

  unsigned nVolumes = args.inputPath.size(); 
  for (unsigned volumeCtr = 0; volumeCtr < nVolumes; volumeCtr++) {
    Path path(args.inputPath[volumeCtr]);

    Volume volume = loadVolume(path, args.axisOrder, args.verbose);
    
    if (volume) {
      Real voxelMin = get_volume_voxel_min(volume);
      Real voxelMax = get_volume_voxel_max(volume);

      get_volume_sizes(volume, sizes.contents());

      unsigned D1 = sizes[0];
      Boolean haveSliceList = Boolean(args.sliceList[0] >= 0);
      for (unsigned D1ctr = 0; D1ctr < D1; D1ctr++) {
	if (!haveSliceList)
	  args.sliceList[D1ctr] = D1ctr;
      }

      unsigned D2 = sizes[1];
      unsigned D3 = sizes[2];

      Real separations[3];
      get_volume_separations(volume, separations);
	
      for (unsigned extremaRunCtr = 0; extremaRunCtr < nExtremaRuns; extremaRunCtr++) {
	for (unsigned maskExtremaRunCtr = 0; 
	     maskExtremaRunCtr < nMaskExtremaRuns; 
	     maskExtremaRunCtr++) {
	  unsigned nVoxels = D1*D2*D3;
	  
	  if (voxels) {
	    voxels->newSize(nVoxels);
	    voxels->resetIterator();
	  }
	  
	  double minVal = MAXDOUBLE;
	  double maxVal = -MAXDOUBLE;
	  double sumVal = 0; 
	  double sum2val = 0;
	  double d1Center = 0, d2Center = 0, d3Center = 0;

	  double ceiling = args.ceiling[extremaRunCtr];
	  double floor   = args.floor[extremaRunCtr];
  
	  double maskFloor = args.maskFloor[maskExtremaRunCtr];
	  double maskCeil  = args.maskCeil[maskExtremaRunCtr];
	  
	  nVoxels = 0;
	  
	  Boolean convertVoxelToWorld = (mask && args.useWorldCoord);
	  if (mask && !args.useWorldCoord && 
	      ((maskD1 != D1) || (maskD2 != D2) || (maskD3 != D3))) {
	    cerr << "WARNING: Dimensions of " << args.maskPath 
		 << " do not match those of " << path 
		 << ". Using world coordinates." << endl;
	    convertVoxelToWorld = TRUE;
	  }
	  
	  if (args.verbose)
	    cout << "Scanning " << path << flush;
	  
	  for (unsigned d1i = 0; d1i < D1; d1i++) {
	    int d1 = args.sliceList[d1i];
	    if ((d1 >= 0) && (d1 < D1)) {
	      for (unsigned d2 = 0; d2 < D2; d2++)
		for (unsigned d3 = 0; d3 < D3; d3++) {
		  Boolean valid = TRUE;
		  if (mask) {
		    Real value = 0;
		    if (convertVoxelToWorld) {
		      Real xWorld, yWorld, zWorld;
		      convert_3D_voxel_to_world(volume, d1, d2, d3, &xWorld, &yWorld, &zWorld);
		      Real voxel1, voxel2, voxel3;
		      convert_3D_world_to_voxel(mask, xWorld, yWorld, zWorld,
						&voxel1, &voxel2, &voxel3);
		      int v1 = ROUND(voxel1);
		      int v2 = ROUND(voxel2);
		      int v3 = ROUND(voxel3);
		      if ((v1 >= 0) && (v2 >= 0) && (v3 >= 0) &&
			  (v1 < maskD1) && (v2 < maskD2) && (v3 < maskD3))
			GET_VALUE_3D(value, mask, v1, v2, v3);
		    }
		    else
		      GET_VALUE_3D(value, mask, d1, d2, d3);
		    
		    if (args.nMaskExtrema)
		      valid = (value >= maskFloor) && (value <= maskCeil);
		    else
		      valid = (value != 0);
		  }
		  
		  if (valid) {
		    Real value;
		    GET_VOXEL_3D(value, volume, d1, d2, d3);
		    // Eliminate NaN's
		    if (!args.ignoreNaN || ((value >= voxelMin) && (value <= voxelMax))) { 
		      value = convert_voxel_to_value(volume, value);
		      
		      if ((value >= floor) && (value <= ceiling)) {
			if (value < minVal)
			  minVal = value;
			if (value > maxVal)
			  maxVal = value;
			sumVal += value;
			sum2val += SQR(value);
			d1Center += value*d1;
			d2Center += value*d2;
			d3Center += value*d3;
			if (voxels)
			  (*voxels)++ = value;
			nVoxels++;
		      }
		    }
		  }
		}
	    }
	    
	    if (args.verbose)
	      cout << "." << flush;
	  }
	  if (args.verbose)
	    cout << "Done" << endl;

	  if (nVoxels) {
	    if (voxels)
	      voxels->newSize(nVoxels);
	    
	    double meanVal  = sumVal/nVoxels;
	    double variance = sum2val/nVoxels - SQR(meanVal);
	    
	    double majorityVal = 0;
	    double biModalT = 0;
	    double entropy = 0;
	    double pctT1 = 0;
	    double pctT2 = 0;
	    
	    if (args.all || args.calcMajority || args.calcBiModalT || args.calcPctT || 
		args.calcEntropy) {
	      unsigned nBins = unsigned(ceil(voxelMax - voxelMin + 1));
	      Histogram hist(minVal, maxVal, nBins);
	      if (args.verbose) 
		cout << "Creating histogram with " << nBins << " bins, bin width "
		     << hist.binWidth() << "..." << flush;
	      add(hist, *voxels);
	      if (args.verbose)
		cout << "Done" << endl;
	      
	      if (args.all || args.calcMajority)	  
		majorityVal = hist.majority();
	      
	      if (args.all || args.calcBiModalT)	  
		biModalT = hist.biModalThreshold();
	      
	      if (args.all || args.calcEntropy)	  
		entropy = hist.entropy();
	      
	      if (args.all || args.calcPctT) {
		pctT1 = hist.pctThreshold(args.pctT);
		pctT2 = hist.pctThreshold(100 - args.pctT);
	      }
	    }
	    
	    if (!args.quiet) {
	      cout << "File:         " << path << endl;
	      if (args.nExtrema)
		cout << "Extrema:      [" << floor << ", " << ceiling << "]" << endl;
	      if (mask) {
		cout << "Mask:         " << args.maskPath << endl;
		if (args.nMaskExtrema)
		  cout << "Mask extrema: [" << maskFloor << ", " << maskCeil << "]" <<endl;
	      }
	      cout << "# voxels:     " << nVoxels << endl;
	      cout << "% of total:   " << 100.0*nVoxels/(D1*D2*D3) << endl;
	    }
	    
	    if (args.all || args.calcVolume || !args.quiet) {
	      if (!args.quiet)
		cout << "Volume (mm3): ";
	      cout << nVoxels*fabs(separations[0]*separations[1]*separations[2]) 
		   << endl;
	    }
	    
	    if (args.all || args.calcMin) {
	      if (!args.quiet)
		cout << "Min:          ";
	      cout << minVal << endl;
	    }
	    
	    if (args.all || args.calcMax) {
	      if (!args.quiet)
		cout << "Max:          ";
	      cout << maxVal << endl;
	    }
	    
	    if (args.all || args.calcSum) {
	      if (!args.quiet)
		cout << "Sum:          ";
	      cout << sumVal << endl;
	    }
	    
	    if (args.all || args.calcSum2) {
	      if (!args.quiet)
		cout << "Sum2:         ";
	      cout << sum2val << endl;
	    }
	    
	    if (args.all || args.calcMean) {
	      if (!args.quiet)
		cout << "Mean:         ";
	      cout << meanVal << endl;
	    }
	    
	    if (args.all || args.calcVariance) {
	      if (!args.quiet)
		cout << "Variance:     ";
	      cout << variance << endl;
	    }
	    
	    if (args.all || args.calcStddev) {
	      if (!args.quiet)
		cout << "Stddev:       ";
	      cout << sqrt(variance) << endl;
	    }
	    
	    // NOTE: THIS DESTROYS THE VOXELS ARRAY!
	    if (args.all || args.calcMedian) {
	      if (!args.quiet)
		cout << "Median:       ";
	      cout << flush << voxels->medianVolatile() << endl;
	    }
	    
	    if (args.all || args.calcMajority) {
	      if (!args.quiet)
		cout << "Majority:     ";
	      cout << majorityVal << endl;
	    }
	    
	    if (args.all || args.calcBiModalT) {
	      if (!args.quiet)
		cout << "BiModalT:     ";
	      cout << biModalT << endl;
	    }
	    
	    if (args.all || args.calcPctT) {
	      if (!args.quiet) {
		cout << "PctT (";
		cout.width(2);
		cout << args.pctT << "):    ";
	      }
	      cout << pctT1 << ' ' << pctT2 << endl;
	    }
	    
	    if (args.all || args.calcEntropy) {
	      if (!args.quiet)
		cout << "Entropy:      ";
	      cout << entropy << endl;
	    }
	    
	    if (args.all || args.calcCoM) {
	      d1Center = d1Center/sumVal;
	      d2Center = d2Center/sumVal;
	      d3Center = d3Center/sumVal;
	      
	      STRING *dimNames = get_volume_dimension_names(volume);
	      
	      if (!args.quiet)
		cout << "CoM_voxel:    ";
	      cout << dimNames[0] << ":" << d1Center << " " 
		   << dimNames[1] << ":" << d2Center << " " 
		   << dimNames[2] << ":" << d3Center << endl;
	      
	      if (!args.quiet)
		cout << "CoM_world:    ";
	      Real xWorld, yWorld, zWorld;
	      convert_3D_voxel_to_world(volume, d1Center, d2Center, d3Center, 
					&xWorld, &yWorld, &zWorld);
	      cout << "zspace:" << zWorld << " " 
		   << "yspace:" << yWorld << " " 
		   << "xspace:" << xWorld << endl;
	      
	      delete_dimension_names(volume, dimNames);
	    }
	  }
	  else {
	    cerr << "Sorry, found no valid voxels in " << path << endl;
	  }
	}
      }

      delete_volume(volume);
    }
  }

  if (voxels)
    delete voxels;

  return EXIT_SUCCESS;
}

//
// Load a volume
//
Volume loadVolume(const Path& rawPath, char **axisOrder, int verbose)
{
  Real amountDone;
  volume_input_struct inputInfo;
  Volume volume = 0;

  Path path(rawPath.expanded());

  MString extension;
  if (!path.existsCompressed(&extension)) {
    cerr << "Couldn't find " << path << endl;
    exit(EXIT_FAILURE);
  }
  path += extension;

  int nDimensions = get_minc_file_n_dimensions(path);
  if (nDimensions < 3) {
    cerr << "ERROR: Volume " << path << " has only " << nDimensions  << " dimensions."
	 << endl << "  (volume_stats requires 3)" << endl;
    exit(EXIT_FAILURE);
  }    
  else if (nDimensions > 3) {
    cerr << "WARNING: Volume " << path << " has " << nDimensions << " dimensions."
	 << endl << " => attempting to read spatial dimensions only." << endl 
	 << "You may want to specify -dimorder, or mincreshape the volume." << endl;
  }
    
  if (verbose)
    cout << "Reading volume " << path << flush;
  if (start_volume_input(path, 3, axisOrder, NC_UNSPECIFIED, 
			 TRUE, 0.0, 0.0, TRUE, &volume, (minc_input_options *) NULL, 
			 &inputInfo) != OK)
    return 0;
  while (input_more_of_volume(volume, &inputInfo, &amountDone))
    if (verbose)
      cout << "." << flush;
  if (verbose)
    cout << "Done" << endl;
  delete_volume_input(&inputInfo);

  return volume;
}
