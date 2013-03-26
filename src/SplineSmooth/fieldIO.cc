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
$RCSfile: fieldIO.cc,v $
$Revision: 1.5 $
$Author: claude $
$Date: 2011-01-13 20:15:17 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : fieldIO.c,v
@INPUT      : 
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: functions for reading and writing field file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : July 9, 1996 (John G. Sled)
@MODIFIED   : 
@COPYRIGHT  : 1996
---------------------------------------------------------------------------- */

#include "../Splines/Spline.h"
#include <EBTKS/MString.h>	// (bert)
#include "../Splines/TBSpline.h"
#undef VIO_ROUND
#undef SIGN
extern "C" {
#include <volume_io.h>
}
#undef VIO_ROUND // Added to avoid conflict between Volume_io's and
             // AZgen's definition of VIO_ROUND, Alex Zijdenbos 97/12/05
#include "splineSmooth.h"
#include "fieldIO.h"
#include <time.h>
#include <iostream>		// (bert)
using namespace std;		// (bert)
// #define DEBUG_FIELDIO

/*--------------------- file format keywords ------------------------------ */

static   const VIO_STR      FIELD_FILE_HEADER = "MNI Field File";
static   const VIO_STR      VERSION_STRING = "Version";
static   const VIO_STR      TYPE_STRING = "Field_Type";
static   const VIO_STR      B_SPLINE_STRING = "B_Spline";
static   const VIO_STR      THIN_PLATE_SPLINE_STRING = "Thin_plate_Spline";
static   const VIO_STR      DISTANCE_STRING = "Distance";
static   const VIO_STR      DOMAIN_STRING = "Domain";
static   const VIO_STR      COEFFICIENTS_STRING = "Coefficients";

static   const VIO_STR      CURRENT_VERSION_STRING = "0.9.0";
/*------------------------------------------------------------------------- */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_transform
@INPUT      : 
@OUTPUT     : 
@RETURNS    : VIO_OK or VIO_ERROR
@DESCRIPTION: Outputs a spline field to the file in MNI field format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : July 9, 1996            
@MODIFIED   : 
---------------------------------------------------------------------------- */
VIO_Status outputCompactField(const MString filename, 
                          const DblMat &domain,  // 3 rows, 2 columns 
                          double distance, const DblArray &coef, 
                          enum spline_type spline, const MString command,
                          VIO_Volume volume)
{
  FILE *file;
  char *date;
  int i;

  file = fopen(filename, "w");
    
  /* --- parameter checking */

  if( file == NULL )
    {
      cerr <<  "outputCompactField(): unable to open file " << filename 
           << endl;
      return( VIO_ERROR );
    }
  
  (void) fprintf( file, "%s\n", FIELD_FILE_HEADER );
  
  // create comment 
  time_t t;
  t = time(NULL);
  date = asctime(localtime(&t));
  date[24] = 0;

  output_comments(file, MString(date) + ">>> " + command);

  (void) fprintf( file, "%s = %s;\n", VERSION_STRING, CURRENT_VERSION_STRING);
  // write out field data
  (void) fprintf( file, "%s = %s;\n", TYPE_STRING, (spline == b_spline) ?
                  B_SPLINE_STRING : THIN_PLATE_SPLINE_STRING);

  (void) fprintf( file, "%s = %.15g;\n", DISTANCE_STRING, distance);

  // convert to world coordinates 
  VIO_Real voxel0[VIO_N_DIMENSIONS], voxel1[VIO_N_DIMENSIONS];
  VIO_Real world0[VIO_N_DIMENSIONS], world1[VIO_N_DIMENSIONS];
  VIO_Real separations[VIO_N_DIMENSIONS];
  get_volume_separations(volume, separations);
  for(i = 0; i < VIO_N_DIMENSIONS; i++) {
    voxel0[i] = domain(i,0)/separations[i];
    voxel1[i] = domain(i,1)/separations[i];
  }
  convert_voxel_to_world(volume,voxel0,&world0[0],&world0[1],&world0[2]);
  convert_voxel_to_world(volume,voxel1,&world1[0],&world1[1],&world1[2]);

  (void) fprintf( file, "%s =", DOMAIN_STRING);
  for(i = 0; i < VIO_N_DIMENSIONS; i++)
    (void) fprintf( file, "\n     %.15g     %.15g", world0[i], world1[i]);
  (void) fprintf( file, ";\n");

  (void) fprintf( file, "%s =", COEFFICIENTS_STRING);
  for(i = 0; i < coef.size(); i++)
    (void) fprintf( file, "\n%25.15g", coef[i]);
  (void) fprintf( file, ";\n");

  fclose(file);
  return VIO_OK;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_transform
@INPUT      : filename
@OUTPUT     : splines 
@RETURNS    : VIO_OK or VIO_ERROR
@DESCRIPTION: Inputs the field from the file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : July 9, 1996            
@MODIFIED   : 
---------------------------------------------------------------------------- */
VIO_Status inputCompactField(VIO_STR filename, Spline **splines, 
                         enum spline_type *type,  
                         VIO_Volume volume)
{
  FILE *file;
  VIO_STR line, type_name, version_name;
  double distance;
  DblMat domain(VIO_N_DIMENSIONS,2);
  VIO_Real separations[VIO_N_DIMENSIONS];
  DblArray coef;
  int i, n;
  VIO_Real *reals;
  VIO_Status status;

  file = fopen(filename, "r");
    
  /* parameter checking */
  
  if( file == NULL )
    {
      cerr <<  "inputCompactField(): unable to open file " << filename 
           << endl;
      return( VIO_ERROR );
    }

  /* okay read the header */

  if( mni_input_string( file, &line, (char) 0, (char) 0 ) != VIO_OK )
    {
      delete_string( line );
      print_error( "inputCompactField(): could not read header in file.\n");
      return( VIO_ERROR );
    }
  if( !equal_strings( line, FIELD_FILE_HEADER ) )
    {
      delete_string( line );
      print_error( "inputCompactField(): invalid header in file.\n");
      return( VIO_ERROR );
    }
  delete_string( line );

  /* --- read the version of the file if it there */
  status = mni_input_string( file, &line, (char) '=', (char) 0 );
  if( status != VIO_OK || mni_skip_expected_character( file, (char) '=' ) != VIO_OK )
    return( status );  
  
  version_name = NULL;  // indicates unversioned file
  if(equal_strings( line, VERSION_STRING )) 
    {
      if( mni_input_string( file, &version_name, (char) ';', (char) 0 ) != VIO_OK )
        {
          print_error( "inputCompactField(): missing version name.\n");
          return( VIO_ERROR );
        }
      if( mni_skip_expected_character( file, (char) ';' ) != VIO_OK )
        return( VIO_ERROR );

      /* read next field */
      status = mni_input_string( file, &line, (char) '=', (char) 0 );
      if( status != VIO_OK || 
          mni_skip_expected_character( file, (char) '=' ) != VIO_OK )
        return( status );
    }

  /* --- read the type of field */
  if(!equal_strings( line, TYPE_STRING )) 
    return(VIO_ERROR);
  
  if( mni_input_string( file, &type_name, (char) ';', (char) 0 ) != VIO_OK )
    {
      print_error( "inputCompactField(): missing field type.\n");
      return( VIO_ERROR );
    }
  if( mni_skip_expected_character( file, (char) ';' ) != VIO_OK )
    return( VIO_ERROR );

  if( equal_strings( type_name, B_SPLINE_STRING ) )
    *type = b_spline;
  else if( equal_strings( type_name, THIN_PLATE_SPLINE_STRING ) )
    *type = thin_plate_spline;
  else {
    cerr << "inputCompactField(): unknown field type " << type_name << endl;  
    return VIO_ERROR;
  }

  delete_string( type_name );
  
  /*  read distance parameter */
  status = mni_input_keyword_and_equal_sign( file, DISTANCE_STRING, TRUE );
  if( status != VIO_OK )
    return( status );
  
  status = mni_input_real(file, &distance);
  if( status != VIO_OK )  return( status );
  if( mni_skip_expected_character( file, (char) ';' ) != VIO_OK )
    return( VIO_ERROR );

  /* read domain parameters */
  status = mni_input_keyword_and_equal_sign( file, DOMAIN_STRING, TRUE );
  if( status != VIO_OK )  return( status );

  status = mni_input_reals(file, &n, &reals);
  if( status != VIO_OK )  return( status );

  if(n != 2*VIO_N_DIMENSIONS) {
    cerr << "inputCompactField(): Incorrect number of domain parameters\n";
    return(VIO_ERROR);
  }

  if(version_name == NULL) // if file is unversioned then don't use world
    {                      //   coordinates 
      for(i = 0; i < VIO_N_DIMENSIONS; i++) {
        domain(i,0) = reals[i*2];
        domain(i,1) = reals[i*2+1];
      }
    }
  else   // convert from world coordinates 
    {
      VIO_Real voxel0[VIO_N_DIMENSIONS], voxel1[VIO_N_DIMENSIONS];
      VIO_Real world0[VIO_N_DIMENSIONS], world1[VIO_N_DIMENSIONS];
      get_volume_separations(volume, separations);
      for(i = 0; i < VIO_N_DIMENSIONS; i++) {
        world0[i] = reals[i*2];
        world1[i] = reals[i*2+1];
      }
      convert_world_to_voxel(volume,world0[0],world0[1],world0[2],voxel0);
      convert_world_to_voxel(volume,world1[0],world1[1],world1[2],voxel1);
      for(i = 0; i < VIO_N_DIMENSIONS; i++) {
        domain(i,0) = voxel0[i]*separations[i];
        domain(i,1) = voxel1[i]*separations[i];
#ifdef DEBUG_FIELDIO
        printf("Domain: [ %lf   %lf ]\n", domain(i,0),domain(i,1));  
#endif
      }
    }
  FREE(reals);

  /* read coefficients */
  status = mni_input_keyword_and_equal_sign( file, COEFFICIENTS_STRING, TRUE );
  if( status != VIO_OK )  return( status );

  status = mni_input_reals(file, &n, &reals);
  if( status != VIO_OK )  return( status );

  coef = DblArray(reals, n);   // copy values to DblArray
  FREE(reals);
  

  // create Spline Basis
  if(*type == b_spline) {
    double voxel[VIO_N_DIMENSIONS] = { 0.0, 0.0, 0.0 };
    int sizes[VIO_N_DIMENSIONS];
    get_volume_sizes(volume, sizes);
    get_volume_separations(volume, separations);

    *splines = new TBSplineVolume(domain, voxel, separations, sizes,
                            distance, 1.0, FALSE);
  }
  else
    *splines = createThinPlateSpline(domain, distance, 1.0, TRUE);

  status = (VIO_Status) (*splines)->putCoefficients(coef);
  if(status == FALSE) {
    cerr << "inputCompactField(): Incorrect number of coefficients\n";
    return(VIO_ERROR);
  }
  return VIO_OK;
}

// layout basis functions for thin plate spline basis
Spline *createThinPlateSpline(const DblMat &domain, double distance,
			      double lambda, int verbose)
{
  DblMat knots;
  IntArray n;
  int i, j, k;
  FloatArray xKnots, yKnots, zKnots;
  TPSpline *tpspline;

  // determine number of basis functions in each dimension
  n.newSize(VIO_N_DIMENSIONS);
  for(i = 0; i < VIO_N_DIMENSIONS; i++)
    n[i] = (int) ceil((domain(i,1) - domain(i,0))/distance) + 1;

  // create array of knot locations for each dimension
  knots.resize(n.max(),VIO_N_DIMENSIONS);
  for(i = 0; i < VIO_N_DIMENSIONS; i++)
    {
      double start = 0.5*(domain(i,0) + domain(i,1) - distance*(n[i]-1));
      for(int j = 0; j < n[i]; j++)
        knots(j,i) = start + distance*j;
    }

  // fill in arrays of knot locations corresponding to a grid
  int nProduct = n[0]*n[1]*n[2];
  xKnots.newSize(nProduct);
  yKnots.newSize(nProduct);
  zKnots.newSize(nProduct);
  int index;
  
  for(i = 0; i < n[0]; i++)
    for(j = 0; j < n[1]; j++)
      for(k = 0; k < n[2]; k++)
	{
	  index = (i*n[1]+j)*n[2] + k;
	  xKnots[index] = knots(i,0);
	  yKnots[index] = knots(j,1);
	  zKnots[index] = knots(k,2);
	}

  Spline::verbose = verbose;
  tpspline = new TPSpline(xKnots, yKnots, zKnots, 1.0/distance);
  tpspline->lambda(lambda); // set lambda
  //  tpspline->radialFunction(rSquare);

  return tpspline;
}


// Reads volume as float volume and return actual data type
VIO_Volume
loadFloatVolume(const MString filename, nc_type *data_type)
{
  VIO_Volume volume;

  /**** READ MINC INPUT VOLUME ****/
  volume_input_struct input_info;
  if (start_volume_input((char *)(const char *)filename, VIO_N_DIMENSIONS, 
		   File_order_dimension_names,
		   NC_UNSPECIFIED, /* data type */ FALSE,
		   NC_UNSPECIFIED, NC_UNSPECIFIED, /* min, max */
		   TRUE, &volume, (minc_input_options *) NULL,
			 &input_info) != VIO_OK)
    {
      cerr << "Failed to read volume: " << filename << ".\n";
      exit(1);
    }
 
  VIO_BOOL signed_flag;
  *data_type = get_volume_nc_data_type(volume, &signed_flag);
  delete_volume_input( &input_info);
  delete_volume( volume );
  
  // open this time using float type
  if (input_volume((char *)(const char *)filename, VIO_N_DIMENSIONS, 
		   File_order_dimension_names,
		   NC_FLOAT, /* data type */ FALSE,
		   NC_UNSPECIFIED, NC_UNSPECIFIED, /* min, max */
		   TRUE, &volume, (minc_input_options *) NULL) != VIO_OK)
    {
      cerr << "Failed to read volume: " << filename << ".\n";
      exit(1);
    }
  return(volume);
}


// Create empty volume as float volume and return actual data type
VIO_Volume
loadEmptyFloatVolume(const MString filename, nc_type *data_type, VIO_BOOL *signed_flag)
{
  VIO_Volume volume;

  /**** READ MINC INPUT VOLUME ****/
  volume_input_struct input_info;
  if (start_volume_input((char *)(const char *)filename, VIO_N_DIMENSIONS, 
		   File_order_dimension_names,
		   NC_UNSPECIFIED, /* data type */ FALSE,
		   NC_UNSPECIFIED, NC_UNSPECIFIED, /* min, max */
		   TRUE, &volume, (minc_input_options *) NULL,
			 &input_info) != VIO_OK)
    {
      cerr << "Failed to read volume: " << filename << ".\n";
      exit(1);
    }
 
  *data_type = get_volume_nc_data_type(volume, signed_flag);
  delete_volume_input( &input_info);
  
  // open this time using float type
  if (start_volume_input((char *)(const char *)filename, VIO_N_DIMENSIONS, 
		   File_order_dimension_names,
		   NC_FLOAT, /* data type */ FALSE,
		   NC_UNSPECIFIED, NC_UNSPECIFIED, /* min, max */
		   TRUE, &volume, (minc_input_options *) NULL,
			 &input_info) != VIO_OK)
    {
      cerr << "Failed to read volume: " << filename << ".\n";
      exit(1);
    }
  delete_volume_input( &input_info);
  alloc_volume_data(volume);

  return(volume);
}

// read in minc volume 
VIO_Volume
loadVolume(const MString filename)
{
  VIO_Volume volume;

  if (input_volume((char *)(const char *)filename, VIO_N_DIMENSIONS, 
		   File_order_dimension_names,
		   NC_UNSPECIFIED, /* data type */ FALSE,
		   NC_UNSPECIFIED, NC_UNSPECIFIED, /* min, max */
		   TRUE, &volume, (minc_input_options *) NULL) != VIO_OK)
    {
      cerr << "Failed to read volume: " << filename << ".\n";
      exit(1);
    }
  return(volume);
}

// returns TRUE if volumes are same size and have same voxel spacing
Boolean compareVolumes(VIO_Volume v1, VIO_Volume v2)
{
  int sizes[2][VIO_N_DIMENSIONS];
  VIO_Real separations[2][VIO_N_DIMENSIONS];
  
  get_volume_separations(v1, separations[0]);
  get_volume_separations(v2, separations[1]);
  get_volume_sizes(v1, sizes[0]);
  get_volume_sizes(v2, sizes[1]);

  for(int i = 0; i < VIO_N_DIMENSIONS; i++)
    if(fabs(separations[0][i] - separations[1][i]) > 1e-7 ||
       sizes[0][i] != sizes[1][i]) {
      return(FALSE);
    }
  return(TRUE);
}



// compute spline function at every point within volume
//  returns minimum and maximum value
void 
smoothVolume(Spline *spline, VIO_Volume volume, double *real_min, double *real_max)
{
  int sizes[VIO_N_DIMENSIONS];
  VIO_Real separations[VIO_N_DIMENSIONS];
  int i,j,k;
  float point[VIO_N_DIMENSIONS];
  VIO_Real value, max, min;

  get_volume_separations(volume, separations);
  get_volume_sizes(volume, sizes);

  point[0] = 0; point[1] = 0; point[2] = 0;
  min = (*spline)(point);
  max = min;

  VIO_progress_struct progress;
  initialize_progress_report(&progress, FALSE, sizes[0], "Smoothing volume");

  for(i = 0; i < sizes[0]; i++)
    {
      for(j = 0; j < sizes[1]; j++)
	for(k = 0; k < sizes[2]; k++)
	  {
	    point[0] = i*separations[0];
	    point[1] = j*separations[1];
	    point[2] = k*separations[2];
	  
	    value = (*spline)(point); 
	    set_volume_real_value(volume, i, j, k, 0, 0, value);
	    // keep track of min and max
	    if(value > max) 
	      max = value;
	    else if(value < min)
	      min = value;
	  }
      update_progress_report(&progress, i+1);
    }
  terminate_progress_report(&progress);

  *real_max = max;
  *real_min = min;
} 

// compute spline function at every point within mask volume, 
// other values are set to zero
//  returns minimum and maximum value
void 
smoothVolume(Spline *spline, VIO_Volume volume, VIO_Volume mask_volume,
	double *real_min, double *real_max)
{
  int sizes[VIO_N_DIMENSIONS];
  VIO_Real separations[VIO_N_DIMENSIONS];
  int i,j,k;
  float point[VIO_N_DIMENSIONS];
  VIO_Real value, max, min;

  get_volume_separations(volume, separations);
  get_volume_sizes(volume, sizes);

  max = 0;
  min = 0;

  VIO_progress_struct progress;
  initialize_progress_report(&progress, FALSE, sizes[0],
			     "Smoothing volume");

  for(i = 0; i < sizes[0]; i++)
    {
      for(j = 0; j < sizes[1]; j++)
	for(k = 0; k < sizes[2]; k++)
	  {
	    if(get_volume_real_value(mask_volume, i, j, k, 0, 0) > 0.5)
	      {
		point[0] = i*separations[0];
		point[1] = j*separations[1];
		point[2] = k*separations[2];
		
		value = (*spline)(point); 
		set_volume_real_value(volume, i, j, k, 0, 0, value);

		// keep track of min and max
		if(value < min) 
		  min = value;
		else if(value > max) 
		  max = value;
	      }
	    else
	      set_volume_real_value(volume, i, j, k, 0, 0, 0.0);
	  }
      update_progress_report(&progress, i+1);
    }
  terminate_progress_report(&progress);

  *real_min = min;
  *real_max = max;
}


// compute spline function at every point within volume
//  returns minimum and maximum value
void 
smoothVolumeLookup(TBSplineVolume *spline, VIO_Volume volume, 
             double *real_min, double *real_max)
{
  int sizes[VIO_N_DIMENSIONS];
  int i,j,k;
  //  float point[VIO_N_DIMENSIONS];
  VIO_Real value, max, min;

  get_volume_sizes(volume, sizes);

  min = (*spline)(0,0,0);
  max = min;

  VIO_progress_struct progress;
  initialize_progress_report(&progress, FALSE, sizes[0], "Smoothing volume");

  for(i = 0; i < sizes[0]; i++)
    {
      for(j = 0; j < sizes[1]; j++)
	for(k = 0; k < sizes[2]; k++)
	  {
	    value = (*spline)(i,j,k); 
	    set_volume_real_value(volume, i, j, k, 0, 0, value);
	    // keep track of min and max
	    if(value > max) 
	      max = value;
	    else if(value < min)
	      min = value;
	  }
      update_progress_report(&progress, i+1);
    }
  terminate_progress_report(&progress);

  *real_max = max;
  *real_min = min;
} 

// compute spline function at every point within mask volume, 
// other values are set to zero
//  returns minimum and maximum value
void 
smoothVolumeLookup(TBSplineVolume *spline, VIO_Volume volume, VIO_Volume mask_volume,
	double *real_min, double *real_max)
{
  int sizes[VIO_N_DIMENSIONS];
  int i,j,k;
  //float point[VIO_N_DIMENSIONS];
  VIO_Real value, max, min;

  get_volume_sizes(volume, sizes);

  max = 0;
  min = 0;

  VIO_progress_struct progress;
  initialize_progress_report(&progress, FALSE, sizes[0],
			     "Smoothing volume");

  for(i = 0; i < sizes[0]; i++)
    {
      for(j = 0; j < sizes[1]; j++)
	for(k = 0; k < sizes[2]; k++)
	  {
	    if(get_volume_real_value(mask_volume, i, j, k, 0, 0) > 0.5)
	      {
		value = (*spline)(i,j,k); 
		set_volume_real_value(volume, i, j, k, 0, 0, value);

		// keep track of min and max
		if(value < min) 
		  min = value;
		else if(value > max) 
		  max = value;
	      }
	    else
	      set_volume_real_value(volume, i, j, k, 0, 0, 0.0);
	  }
      update_progress_report(&progress, i+1);
    }
  terminate_progress_report(&progress);

  *real_min = min;
  *real_max = max;
}


// write volume to disk with specified data type
void 
outputVolume(VIO_Volume volume, const MString filename, nc_type output_type, 
	     VIO_BOOL signed_flag, VIO_Real real_min,  VIO_Real real_max, const MString command)
{

  /* This seems like a hack. However, it fixes some round off problems. */ 
  if(output_type == NC_FLOAT || output_type == NC_DOUBLE) {
    //set_volume_voxel_range(volume, 0, 0); // set to maximum range
    //VF: no it isn't
    set_volume_real_range(volume, real_min, real_max);
  }
  else {  // for fixed point types
    set_volume_real_range(volume, real_min, real_max);
  }

  if(output_volume((char *)(const char *)filename, output_type,
		   signed_flag, 0, 0,
		   volume, (char *)(const char *)command,
		   (minc_output_options *) NULL) != VIO_OK )
    {
      cerr << "\nError: failed to write volume: " << filename << endl;
      exit( 1 );
    }
}




