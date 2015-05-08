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
public:
  // Constructors/destructor
  // The min and max values in these constructors are bin centers,
  // not absolute extrema.
  WHistogram() : DHistogram() {;}
  // min and max refer to bin centers
  WHistogram(double min, double max, unsigned nBins = 0) : 
    DHistogram(min, max, nBins) {;}
  WHistogram(double min, double max, double binWidth) : 
    DHistogram(min, max, binWidth) {;}
  WHistogram(unsigned nBins, double min = 0.0, double binWidth = 1.0) :
    DHistogram(nBins, min, binWidth) {;}
  WHistogram(const WHistogram& hist) : DHistogram(hist) {;}

  WHistogram& operator = (const WHistogram&);
  // Changes ranges; keeps binWidth
  WHistogram& newRange(double min, double max); 

  // Set functions
  virtual Boolean add(double value) 
  {
    // throw away value beyond top half of last bin and bottom half of first 
    if ((value < _cmin) || (value > _cmax))
      return FALSE;
    
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










