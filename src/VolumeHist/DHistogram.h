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
$RCSfile: DHistogram.h,v $
$Revision: 1.1 $
$Author: bert $
$Date: 2003-04-16 14:31:39 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ------------------------------ MNI Header ----------------------------------
#@NAME       : DHistogram.h
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: header for floating point histogram class
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : May 30, 1996      J.G.Sled  (based on Histogram class created
#            :   by Alex Zijdenbos)
#@MODIFIED   :
#  $Id: DHistogram.h,v 1.1 2003-04-16 14:31:39 bert Exp $
-----------------------------------------------------------------------------*/
#ifndef DHISTOGRAM_H
#define DHISTOGRAM_H

#include <iostream>		/* (bert) */
using namespace std;		/* (bert) */
#include <EBTKS/MTypes.h>	/* (bert) - Added EBTKS subdirectory */
#include <EBTKS/ValueMap.h>	/* (bert) */
#include <EBTKS/SimpleArray.h>	/* (bert) */
#include <math.h>

class DHistogram  : protected SimpleArray<double> 
{
protected:
  double _min, _max; // True extrema (i.e., not the bin centers)
  double _cmin, _cmax; // extrema based on bin centers
  double _binWidth;
  LinearMap _valueToBinMap;
public:
  // Constructors/destructor
  // The min and max values in these constructors are bin centers,
  // not absolute extrema.
  DHistogram();
  // min and max refer to bin centers
  DHistogram(double min, double max, unsigned nBins = 0);
  DHistogram(double min, double max, double binWidth);
  DHistogram(unsigned nBins, double min = 0.0, double binWidth = 1.0);
  DHistogram(const DHistogram&);
  virtual ~DHistogram() {}

  DHistogram& operator = (const DHistogram&);
  // Changes ranges; keeps binWidth
  DHistogram& newRange(double min, double max); 

  // Get functions
  double  operator [] (int i) const { 
    return SimpleArray<double>::operator [] (i); }
  double& operator [] (int i) { 
    return SimpleArray<double>::operator [] (i); }
  double    binWidth() const            { return _binWidth; }
  unsigned  nBins() const               { return _size; }
  double    binStart(unsigned i) const  { return _min + i*_binWidth; }
  DblArray  binStarts() const;
  double    binCenter(unsigned i) const { return binStart(i) + _binWidth/2; }
  DblArray  binCenters() const          { return binStarts() + _binWidth/2; }
  double    n() const                   { return sum(); }
  double    count(double value) const   { return _contents[bin(value)]; }
  unsigned  bin(double value) const;
  double    max(unsigned *bin = 0) const;
  double    mean() const;
  double    median(unsigned *bin = 0, double nBelow = 0, double nAbove = 0) const;
  double    majority(unsigned *bin = 0) const;
  double    biModalThreshold() const;
  double    varianceThreshold() const;
  double    kullbackThreshold() const;
  double    entropy() const;
  DblArray  pdf() const;
  DblArray  cdf() const;

  // Set functions
  // Set functions
  virtual Boolean add(double value) {
    if ((value < _min) || (value > _max))
      return FALSE;
    unsigned index = (unsigned) _valueToBinMap(value);
    if (index >= _size)
      index = _size - 1;
    _contents[index]++;
    return TRUE;
  }

// Other operators
  DHistogram&  operator += (const DHistogram& hist);
  LUT<double> equalize(const DHistogram& hist) const;

// I/O
  ostream& printHeadAndTail(ostream& os, unsigned n = 10) const;
#ifdef HAVE_MATLAB
  Boolean  saveMatlab(const char *fileName, const char *binVarName = "X", 
		      const char *countVarName = "N", const char *option = "u") const;
#endif

// Friends
  friend ostream&    operator << (ostream& os, const DHistogram& hist);
  friend DblArray    pdf(const DHistogram& hist) { return hist.pdf(); }
  friend DblArray    cdf(const DHistogram& hist) { return hist.cdf(); }
  friend LUT<double> equalize(const DHistogram& hist1, const DHistogram& hist2)   { return hist1.equalize(hist2); }

  friend SimpleArray<double>  asDblArray(const DHistogram& hist);
};

#endif










