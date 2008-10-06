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

Status outputCompactField(const MString filename, 
                          const DblMat &domain,  // 3 rows, 2 columns 
                          double distance, const DblArray &coef, 
                          enum spline_type spline, const MString command,
                          Volume volume);
Status inputCompactField(STRING filename, Spline **splines,
                         enum spline_type *type,  Volume volume);
Spline *createThinPlateSpline(const DblMat &domain, double distance,
			      double lambda, int verbose);
Volume loadFloatVolume(const MString filename, nc_type *data_type);
Volume loadEmptyFloatVolume(const MString filename, nc_type *data_type, BOOLEAN *signed_flag);
Volume loadVolume(const MString filename);
Boolean compareVolumes(Volume v1, Volume v2);
void smoothVolume(Spline *spline, Volume volume, 
		  double *real_min, double *real_max);
void smoothVolume(Spline *spline, Volume volume, Volume mask_volume,
		   double *real_min, double *real_max);
void smoothVolumeLookup(TBSplineVolume *spline, Volume volume, 
		  double *real_min, double *real_max);
void smoothVolumeLookup(TBSplineVolume *spline, Volume volume,
                        Volume mask_volume,
                        double *real_min, double *real_max);
void outputVolume(Volume volume, const MString filename, nc_type output_type,
		  BOOLEAN signed_flag, Real real_min, Real real_max, 
		  const MString command);

#endif














