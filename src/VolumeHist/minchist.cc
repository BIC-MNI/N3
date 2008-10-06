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
$RCSfile: minchist.cc,v $
$Revision: 1.3 $
$Author: rotor $
$Date: 2008-10-06 02:10:05 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : minchist
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : error status
@DESCRIPTION: creates histograms of intensity in a minc volume
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : January 19, 1996 (J.G.Sled)
@MODIFIED   : Log: minchist.c,v 
 * Revision 1.1.1.1  1996/02/27  23:45:21  jgsled
 * volumeHist sources
 *
@COPYRIGHT  : 1996 John G. Sled   
---------------------------------------------------------------------------- */

#include <config.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <iostream>		/* (bert) */
using namespace std;		/* (bert) */
extern "C" {
#include <volume_io.h>
}
#undef ROUND
#include "DHistogram.h"
#include "WHistogram.h"
#include "args.h"

void write_histogram_to_text_file(FILE *fp, int select, DHistogram *histogram);
void bin_volume(Volume volume, Volume mask_volume, args &args, 
           DHistogram *histogram[], int n_histograms);

int main( int argc,  char *argv[] )
{
  Volume volume, mask_volume;     
  Real real_min, real_max;
  args args(argc, argv);
  int i, j, k;

   /* use slice based volume caching to reduce memory usage */
  set_n_bytes_cache_threshold(0);           // Always cache volume
  set_cache_block_sizes_hint(SLICE_ACCESS); // Cache volume slices
  set_default_max_bytes_in_cache(0);        // keep only one slice at a time

  // open mask volume
  if(args.mask_flag == TRUE)
  {
    if (input_volume((char *)args.maskPath.string(), N_DIMENSIONS, 
                     File_order_dimension_names,
                     NC_UNSPECIFIED, /* data type */ FALSE,
                     NC_UNSPECIFIED, NC_UNSPECIFIED, /* min, max */
                     TRUE, &mask_volume, (minc_input_options *) NULL) != OK)
      return(1);

    // check that mask is a label volume
    nc_type mask_type;
    BOOLEAN signed_flag;
    Real voxel_min, voxel_max;
    mask_type = get_volume_nc_data_type(mask_volume, &signed_flag);
    get_volume_voxel_range(mask_volume, &voxel_min, &voxel_max);
    get_volume_real_range(mask_volume, &real_min, &real_max);
    
    if((mask_type != NC_SHORT && mask_type != NC_BYTE) /* || 
        voxel_min != real_min || voxel_max != real_max */ )
      {
        cerr << "Error: " << args.maskPath.string() 
             << " is not a label volume." << endl;
        exit(EXIT_FAILURE);
      }
  }

  // open source volume
  if (input_volume((char *)args.inputPath.string(), N_DIMENSIONS, 
		   File_order_dimension_names,
		   NC_UNSPECIFIED, /* data type */ FALSE,
		   NC_UNSPECIFIED, NC_UNSPECIFIED, /* min, max */
		   TRUE, &volume, (minc_input_options *) NULL) != OK)
    return(1);

  // allocate histograms
  int n_histograms = (args.mask_flag) ? ROUND(real_max + 1) : 1;
  if(args.select_class_flag && (real_max < args.selected_class ||
                                  real_min > args.selected_class)) {
    cerr << "Error: no labels of selected class in mask volume." << endl;
    exit(EXIT_FAILURE);
  }
  if(args.limit_classes == TRUE)
    n_histograms = MIN(n_histograms, args.number_of_classes);
  if(n_histograms > 100) 
    cerr << "Warning: mask volume may not me a label volume" << endl;

  // check that volumes are same size
  int sizes[N_DIMENSIONS], mask_sizes[N_DIMENSIONS];
  get_volume_sizes(volume, sizes);
  if(args.mask_flag == TRUE) 
    {
      get_volume_sizes(mask_volume, mask_sizes);
      for(i = 0; i < 3; i++)
        if(sizes[i] != mask_sizes[i])
          {
            cerr << "Mask volume and source volume must be same size" << endl;
            exit(EXIT_FAILURE);
          }
    }
  
  // setup bins for histograms
  get_volume_real_range(volume, &real_min, &real_max);
  DHistogram **histogram = (DHistogram **) 
    malloc(sizeof(DHistogram *)*n_histograms);
  Real *class_min = new Real[n_histograms];
  Real *class_max = new Real[n_histograms];
  
  if(args.auto_range == FALSE)
    {
      for(i = 0; i < n_histograms; i++)
        {
          if(args.window_flag == FALSE) 
            histogram[i] = new DHistogram(real_min, real_max, 
                                          (unsigned)args.number_of_bins);
          else
            histogram[i] = new WHistogram(real_min, real_max, 
                                          (unsigned)args.number_of_bins);
          class_min[i] = real_min;
          class_max[i] = real_max;
        }
    }
  else // compute range for each class
    {
      progress_struct progress;
      int bin;
      Real value;

      initialize_progress_report(&progress, FALSE, sizes[0], 
			    "Computing ranges for each class.");
      for(i = 0; i < n_histograms; i++)
	{
	  class_min[i] = real_max;
	  class_max[i] = real_min;
	}
      for(i = 0; i < sizes[0]; i++)
	{
	  for(j = 0; j < sizes[1]; j++)
            if(args.mask_flag == TRUE) 
              {
                for(k = 0; k < sizes[2]; k++)
                  {
                    bin = ROUND(get_volume_real_value(mask_volume,
                                                      i, j, k, 0, 0));
                    if(bin >= 0 && bin < n_histograms)
                      {
                        value = get_volume_real_value(volume,i,j,k,0,0);
                        if(value < class_min[bin]) 
                          class_min[bin] = value;
                        else if(value > class_max[bin])
                          class_max[bin] = value;
                      }
                  }
              }
            else
              {
                for(k = 0; k < sizes[2]; k++)
                  {
                    value = get_volume_real_value(volume,i,j,k,0,0);
                    if(value < class_min[0]) 
                      class_min[0] = value;
                    else if(value > class_max[0])
                      class_max[0] = value;
                  }
              }
	  update_progress_report(&progress, i+1);
	}
      terminate_progress_report(&progress);

      for(i = 0; i < n_histograms; i++)
	{
	  if(class_max[i] <= class_min[i]) 
	    class_max[i] = class_min[i] + 1.0;
          if(args.window_flag == FALSE) 
            histogram[i] = new DHistogram(class_min[i], class_max[i], 
                                      (unsigned)args.number_of_bins);
          else
            histogram[i] = (DHistogram *) new WHistogram(class_min[i],
                             class_max[i], (unsigned)args.number_of_bins);
	}
    }

  // fill out histograms
  bin_volume(volume, mask_volume, args, histogram, n_histograms);

  // delete output file if it exists
  if(args.outputPath.exists())
    remove_file(args.outputPath.string());

#ifdef HAVE_MATLAB
  if(args.matlab_format == TRUE)
    {
      // save histograms as matlab files
      MString bin_name = "bins_";
      MString count_name = "counts_";
      char number[50];

      if(args.select_class_flag == TRUE)  // write out one histogram
 	{
	  sprintf(number,"%u", args.selected_class);
          histogram[args.selected_class]->saveMatlab(args.outputPath.string(),
                         bin_name + number, count_name + number, "update");
	}
      else  // write out all histograms
	{
	  for(i = 0; i < n_histograms; i++)
	    {
	      sprintf(number,"%u", i);
              histogram[i]->saveMatlab(args.outputPath.string(),
                                      bin_name + number,
                                      count_name + number, "update");
            }
	}
    }
  else  // if text format is desired
#endif
    {
      FILE *fp;
      
      fp = fopen(args.outputPath.string(),"w");
      if(fp == NULL) 
	{
	  cerr << "Can't open output file: " << args.outputPath.string()
	       << endl;
	  exit(EXIT_FAILURE);
	}
      if(args.select_class_flag == TRUE)  // write out one histogram
        write_histogram_to_text_file(fp,args.selected_class,
                                       histogram[args.selected_class]);
      else  // write out all histograms
	{
	  for(i = 0; i < n_histograms; i++)
            write_histogram_to_text_file(fp, i, histogram[i]);
	}
      fclose(fp);
    }
  
  return(0);
} 



// write bin centers and count to open text file
void 
write_histogram_to_text_file(FILE *fp, int select, DHistogram *histogram)
{
  fprintf(fp, "# histogram for class %d\n"
          "#  domain: %lf  %lf\n"
          "#  entropy: %lg\n",
          select, histogram->binCenter(0),
          histogram->binCenter(histogram->nBins()-1), 
          histogram->entropy());
  fprintf(fp, "#  bin centers     counts\n");
  for(int i = 0; i < histogram->nBins(); i++)
    {
      fprintf(fp, "  %lf       %lf\n", histogram->binCenter(i), 
	      (*histogram)[i]);
    }
}


// bin all values in source volume
void
bin_volume(Volume volume, Volume mask_volume, args &args, 
           DHistogram *histogram[], int n_histograms)
{
  progress_struct progress;
  int bin, i,j,k;
  int sizes[N_DIMENSIONS];
  get_volume_sizes(volume, sizes);
  
  initialize_progress_report(&progress, FALSE, sizes[0], 
			    "Putting values into bins");
  for(i = 0; i < sizes[0]; i++)
    {
      if(args.mask_flag == TRUE)
        {
          // case: with mask
          for(j = 0; j < sizes[1]; j++)
            for(k = 0; k < sizes[2]; k++)
              {
                bin = ROUND(get_volume_real_value
                            (mask_volume, i, j, k, 0, 0));
                if(bin >= 0 && bin < n_histograms)
                  {
                    histogram[bin]->
                      add(get_volume_real_value(volume,i,j,k,0,0));
                  }
              }
        }
      else
        { // case: without mask
          for(j = 0; j < sizes[1]; j++)
            for(k = 0; k < sizes[2]; k++)
              histogram[0]->add(get_volume_real_value(volume,i,j,k,0,0));
          update_progress_report(&progress, i+1);
        }
    }
  terminate_progress_report(&progress);
}





