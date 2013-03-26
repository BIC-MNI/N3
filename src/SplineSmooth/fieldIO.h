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
$RCSfile: fieldIO.h,v $
$Revision: 1.2 $
$Author: rotor $
$Date: 2008-10-06 02:10:05 $
$State: Exp $
--------------------------------------------------------------------------*/

#ifndef FIELDIO_H
#define FIELDIO_H

VIO_Status outputCompactField(const MString filename, 
                          const DblMat &domain,  // 3 rows, 2 columns 
                          double distance, const DblArray &coef, 
                          enum spline_type spline, const MString command,
                          VIO_Volume volume);
VIO_Status inputCompactField(VIO_STR filename, Spline **splines,
                         enum spline_type *type,  VIO_Volume volume);
Spline *createThinPlateSpline(const DblMat &domain, double distance,
			      double lambda, int verbose);
VIO_Volume loadFloatVolume(const MString filename, nc_type *data_type);
VIO_Volume loadEmptyFloatVolume(const MString filename, nc_type *data_type, VIO_BOOL *signed_flag);
VIO_Volume loadVolume(const MString filename);
Boolean compareVolumes(VIO_Volume v1, VIO_Volume v2);
void smoothVolume(Spline *spline, VIO_Volume volume, 
		  double *real_min, double *real_max);
void smoothVolume(Spline *spline, VIO_Volume volume, VIO_Volume mask_volume,
		   double *real_min, double *real_max);
void smoothVolumeLookup(TBSplineVolume *spline, VIO_Volume volume, 
		  double *real_min, double *real_max);
void smoothVolumeLookup(TBSplineVolume *spline, VIO_Volume volume,
                        VIO_Volume mask_volume,
                        double *real_min, double *real_max);
void outputVolume(VIO_Volume volume, const MString filename, nc_type output_type,
		  VIO_BOOL signed_flag, VIO_Real real_min, VIO_Real real_max,
		  const MString command);

#endif














