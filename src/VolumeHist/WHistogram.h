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
$RCSfile: WHistogram.h,v $
$Revision: 1.1 $
$Author: bert $
$Date: 2003-04-16 14:31:39 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ------------------------------ MNI Header ----------------------------------
#@NAME       : WHistogram.h
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: header for triangular Parsen window histogram class
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : May 19, 1996      J.G.Sled  (based on Histogram class created
#            :   by Alex Zijdenbos)
#@MODIFIED   :
#  $Id: WHistogram.h,v 1.1 2003-04-16 14:31:39 bert Exp $
-----------------------------------------------------------------------------*/
#ifndef WHISTOGRAM_H
#define WHISTOGRAM_H

#include <iostream>		/* (bert) */
using namespace std;		/* (bert) */
#include <EBTKS/MTypes.h>	/* (bert) - Added EBTKS subdirectory */
#include <EBTKS/ValueMap.h>	/* (bert) */
#include <EBTKS/SimpleArray.h>	/* (bert) */
#include "DHistogram.h"
#include <math.h>

class WHistogram  : public DHistogram
{
protected:
  // Standard deviation of the Parzen window, in bin widths.  Zero selects the
  // original behaviour: each sample split linearly between the two bin centers
  // it falls between, which is a triangular kernel exactly one bin wide.
  double _sigma;

public:
  // Constructors/destructor
  // The min and max values in these constructors are bin centers,
  // not absolute extrema.
  WHistogram() : DHistogram(), _sigma(0.0) {;}
  // min and max refer to bin centers
  WHistogram(double min, double max, unsigned nBins = 0) :
    DHistogram(min, max, nBins), _sigma(0.0) {;}
  WHistogram(double min, double max, double binWidth) :
    DHistogram(min, max, binWidth), _sigma(0.0) {;}
  WHistogram(unsigned nBins, double min = 0.0, double binWidth = 1.0) :
    DHistogram(nBins, min, binWidth), _sigma(0.0) {;}
  WHistogram(const WHistogram& hist) : DHistogram(hist), _sigma(hist._sigma) {;}

  WHistogram& operator = (const WHistogram&);
  // Changes ranges; keeps binWidth
  WHistogram& newRange(double min, double max);

  // Get functions
  double windowSigma() const { return _sigma; }

  // Set functions
  // Width of the Gaussian window in bin widths; zero restores the triangular
  // one.  The width is in bins rather than in intensity units, so it follows
  // the range -auto_range picks rather than being fixed in the data's units.
  void setWindowSigma(double sigma) { _sigma = sigma; }

  virtual Boolean add(double value)
  {
    // throw away value beyond top half of last bin and bottom half of first
    if ((value < _cmin) || (value > _cmax))
      return FALSE;

    if (_sigma > 0.0)
      return addWindowed(value);

    double loc = (value - _min)/_binWidth;
    int index = int(::floor(loc));
    double offset = loc - index - 0.5;
    if(offset == 0)   // deal with this seperately to avoid 
      {                  //  problems at end points
        _contents[index]++;
        return TRUE;
      }
    else if(offset > 0 && index <= _size-2)
      {
        _contents[index] += 1.0-offset;
        _contents[index+1] += offset;
        return TRUE;
      }
    else if(index >= 1) 
      {
        _contents[index] += 1.0 + offset;
        _contents[index-1] -= offset;
        return TRUE;
      }
    else
      return FALSE;
  }

protected:
  // Distribute one sample over the bins with a Gaussian kernel of standard
  // deviation _sigma bin widths, in place of the linear split above.  The
  // kernel is evaluated at the bin centers, truncated at four standard
  // deviations, and normalized per sample, so that every sample the range
  // admits contributes exactly one count however near an end of the range it
  // falls.  Only the shape of the kernel differs from add() above; which
  // samples are kept does not.
  Boolean addWindowed(double value)
  {
    // Bin centers sit at the integers of this coordinate, and the caller has
    // already rejected anything outside the outermost two, so the nearest bin
    // is always one of this histogram's own.
    double loc = (value - _cmin)/_binWidth;
    int nearest = int(::floor(loc + 0.5));

    int radius = int(::ceil(4.0*_sigma));
    if(radius < 1)
      radius = 1;
    int last  = int(_size) - 1;
    int lower = (nearest - radius > 0) ? nearest - radius : 0;
    int upper = (nearest + radius < last) ? nearest + radius : last;

    // Exponents are taken relative to the nearest bin's, as a softmax is, so
    // that a window much narrower than one bin cannot underflow to no counts
    // at all.  That bin then carries a weight of exactly one and the total
    // cannot be zero.
    double reference = _exponent(loc, nearest);
    double total = 0.0;
    int i;

    for(i = lower; i <= upper; i++)
      total += ::exp(reference - _exponent(loc, i));
    for(i = lower; i <= upper; i++)
      _contents[i] += ::exp(reference - _exponent(loc, i))/total;

    return TRUE;
  }

  double _exponent(double loc, int bin) const
  {
    double distance = (loc - bin)/_sigma;
    return 0.5*distance*distance;
  }

public:
// Other operators
  /*WHistogram&  operator += (const WHistogram& hist);*/
  /*LUT<double> equalize(const WHistogram& hist) const*/

// Friends
  friend ostream&    operator << (ostream& os, const WHistogram& hist);
  friend DblArray    pdf(const WHistogram& hist) { return hist.pdf(); }
  friend DblArray    cdf(const WHistogram& hist) { return hist.cdf(); }
  friend LUT<double> equalize(const WHistogram& hist1, const WHistogram& hist2)   { return hist1.equalize(hist2); }

  friend SimpleArray<double>  asDblArray(const WHistogram& hist);
};

#endif










