/* ------------------------------ MNI Header ----------------------------------
#@NAME       : GHistogram.h
#@DESCRIPTION: header for Gaussian Parzen window histogram class
#@CREATED    : August 3, 2026
-----------------------------------------------------------------------------*/
#ifndef GHISTOGRAM_H
#define GHISTOGRAM_H

#include <vector>
#include "DHistogram.h"
#include <math.h>

// The kernel is evaluated out to this many standard deviations and truncated
// there.  At 4 sigma the tails carry 6e-5 of a sample, well under the six
// decimals the counts are written with.
#define GHISTOGRAM_RADIUS 4.0

// A histogram whose samples are spread over the bins by a Gaussian kernel of
// stated width, rather than by the triangle one bin wide that WHistogram uses.
// The width is in bin widths rather than in intensity units, so it follows
// -auto_range as the range moves from iteration to iteration.
//
// Weights are normalized per sample over the bins present, so every retained
// sample contributes exactly 1 to the total, as under the linear split.  Near
// the ends of the range the kernel is renormalized rather than truncated.
class GHistogram : public DHistogram
{
public:
  // min and max refer to bin centers; sigma is in bin widths
  GHistogram(double min, double max, unsigned nBins, double sigma) :
    DHistogram(min, max, nBins) { setSigma(sigma); }
  GHistogram(const GHistogram& hist) :
    DHistogram(hist), _sigma(hist._sigma), _radius(hist._radius),
    _weight(hist._weight) {;}

  double sigma() const { return _sigma; }

  // Set functions
  virtual Boolean add(double value);

protected:
  double _sigma;               // standard deviation, in bin widths
  int    _radius;              // kernel truncation, in bins
  std::vector<double> _weight; // scratch, 2*_radius+1 long

  void setSigma(double sigma)
  {
    _sigma = sigma;
    _radius = int(::ceil(GHISTOGRAM_RADIUS*sigma));
    if(_radius < 1)   // a narrower kernel still has to reach both neighbours,
      _radius = 1;    //   so that a sample halfway between centers is halved
    _weight.resize(2*_radius + 1);
  }
};

#endif
