/* ------------------------------ MNI Header ----------------------------------
#@NAME       : GHistogram.cc
#@DESCRIPTION: Gaussian Parzen window histogram class
#@CREATED    : August 3, 2026
-----------------------------------------------------------------------------*/
#include <config.h>
#include "GHistogram.h"

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GHistogram::add
@INPUT      : value - the sample
@RETURNS    : TRUE if the sample fell inside the range
@DESCRIPTION: Spreads one sample over the bins with a Gaussian kernel.
@METHOD     : The rejection rule is WHistogram's: outside the outermost bin
              centers a sample is dropped, not clipped.  Only the shape of the
              kernel differs.
---------------------------------------------------------------------------- */
Boolean
GHistogram::add(double value)
{
  // throw away value beyond top half of last bin and bottom half of first
  if((value < _cmin) || (value > _cmax))
    return FALSE;

  // Bin centers sit at the integers of this coordinate.  Note that this is not
  // WHistogram's coordinate, which puts the bin edges there instead.
  double loc = (value - _cmin)/_binWidth;
  double center = ::nearbyint(loc);  // ties to even
  int n = 2*_radius + 1;
  int j;

  // Bins outside the array are given an infinite exponent, and so no weight.
  double minExponent = HUGE_VAL;
  for(j = 0; j < n; j++)
    {
      double index = center + double(j - _radius);
      if(index >= 0 && index < double(_size))
        {
          double offset = (loc - index)/_sigma;
          _weight[j] = 0.5*(offset*offset);
          if(_weight[j] < minExponent)
            minExponent = _weight[j];
        }
      else
        _weight[j] = HUGE_VAL;
    }

  // Shifted by the nearest bin's exponent before exponentiating, the way a
  // softmax is: the sample sits inside the range, so that bin is one of the
  // ones being kept and the shifted weights cannot all underflow.  Unshifted
  // they do, as soon as sigma falls well below a bin width.
  double sum = 0.0;
  for(j = 0; j < n; j++)
    {
      _weight[j] = (_weight[j] == HUGE_VAL) ? 0.0
                                            : ::exp(minExponent - _weight[j]);
      sum += _weight[j];
    }

  // Normalizing per sample rather than truncating leaves every retained sample
  // contributing exactly 1, as under the linear split.
  for(j = 0; j < n; j++)
    if(_weight[j] != 0.0)
      _contents[int(center) + j - _radius] += _weight[j]/sum;

  return TRUE;
}
