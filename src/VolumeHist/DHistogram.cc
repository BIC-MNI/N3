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
$RCSfile: DHistogram.cc,v $
$Revision: 1.1 $
$Author: bert $
$Date: 2003-04-16 14:31:39 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "DHistogram.h"
#include <assert.h>

#ifdef HAVE_MATLAB
extern "C" {    
  #include"mat.h" 
}
#endif

//
// Constructors/desctructor
//

DHistogram::DHistogram()
: SimpleArray<double>(0)
{
  _min = _max = 0;
  _cmin = _cmax = 0;
}

DHistogram::DHistogram(double min, double max, unsigned nBins)
: SimpleArray<double>(nBins)
{
  if (max < min)
    swap(min, max);

// Hm. These exception handlers may lead to unexpected results in some cases...

  if (!nBins) {
    nBins = unsigned(::ceil(max - min + 1));
    newSize(nBins);
    _binWidth = 1.0;
  }
  else {
    if (max == min)
      _binWidth = 1.0/nBins;
    else if (nBins <= 1) 
      _binWidth = max - min;
    else
      _binWidth = (max - min)/(nBins - 1);
  }
  _min = min - _binWidth/2;
  _max = max + _binWidth/2;
  _valueToBinMap(_min, _max, 0, nBins);
  clear(0);
  _cmin = min;
  _cmax = max;
}

DHistogram::DHistogram(double min, double max, double binWidth)
: SimpleArray<double>((unsigned) ::ceil((max - min)/binWidth + 1))
{
  _binWidth = binWidth;
  _min = min - _binWidth/2;
  _max = max + _binWidth/2;
  _valueToBinMap(_min, _max, 0, _size);
  clear(0);
  _cmin = min;
  _cmax = max;
}

DHistogram::DHistogram(unsigned nBins, double min, double binWidth)
: SimpleArray<double>(nBins)
{
  _binWidth = binWidth;
  _min = min - _binWidth/2;
  _max = _min + _binWidth*nBins;
  _valueToBinMap(_min, _max, 0, nBins);
  clear(0);
  _cmin = min;
  _cmax = _max - _binWidth/2;
}

DHistogram::DHistogram(const DHistogram& hist)
{
  *this = hist;
}

DHistogram&
DHistogram::operator = (const DHistogram& hist)
{
  if (this == &hist) return *this;

  newSize(hist._size);
  
  double *sourcePtr = hist._contents;
  double *destPtr   = _contents;
  for (unsigned i = _size; i; i--)
    *destPtr++ = *sourcePtr++;

  _min = hist._min;
  _max = hist._max;
  _cmin = hist._cmin;
  _cmax = hist._cmax;
  _binWidth = hist._binWidth;
  _valueToBinMap = hist._valueToBinMap;

  return *this;
}

DHistogram&
DHistogram::newRange(double min, double max)
{
  if (!_binWidth) {
    cerr << "Error: DHistogram::newRange() called but bin width not set" << endl;
    return *this;
  }

  _min = min - _binWidth/2;
  _max = max + _binWidth/2;
  _cmin = min;
  _cmax = max;

  unsigned nBins = unsigned(::ceil((max - min)/_binWidth + 1));
  int firstBin = (int) _valueToBinMap(min);
  int lastBin  = (int) _valueToBinMap(max);

  DblArray newHist(nBins);
  newHist.clear(0);
  
  double *sourcePtr = _contents;
  double *destPtr   = newHist.contents();

  if (firstBin < 0) {
    destPtr -= firstBin;
    nBins += firstBin;
  }
  else
    sourcePtr += firstBin;

  if (lastBin < 0)
    nBins = 0;
  else if (lastBin > _size - 1)
    nBins -= lastBin - _size + 1;
  
  for (unsigned i = nBins; i; i--)
    *destPtr++ = *sourcePtr++;

  this->absorb(newHist);

  return *this;
}

//
// Get functions
//

DblArray
DHistogram::binStarts() const
{
  DblArray starts(_size);
  double *startPtr = starts.contents();
  double binStart = _min;
  for (unsigned i = _size; i; i--) {
    *startPtr++ = binStart;
    binStart += _binWidth;
  }

  return starts;
}

unsigned
DHistogram::bin(double value) const
{
  if (!_size) {
    cerr << "DHistogram::bin() called on empty DHistogram" << endl;
    return 0;
  }

  if (value <= _min)
    return 0;

  if (value >= _max)
    return _size - 1;

  unsigned index = (unsigned) _valueToBinMap(value);
  if (index >= _size)
    index = _size - 1;

  return index;
}

double
DHistogram::max(unsigned *) const
{
  cerr << "DHistogram::max invalid; should use DHistogram::majority" << endl;
  return 0;
}

double
DHistogram::mean() const
{
  if (!_size) {
    cerr << "Warning! DHistogram::mean() called on empty DHistogram" << endl;
    return 0.0;
  }

  double meanValue = 0.0;
  double *binPtr = _contents;
  double offset = _binWidth/2;
  for (unsigned i = 0; i < _size; i++){
    meanValue += (*binPtr++)*(_valueToBinMap.reverse(i) + offset);
  }

  return (meanValue/n());
}

double
DHistogram::median(unsigned *bin, double nBelow, double nAbove) const
{
  if (!_size) {
    cerr << "Warning! DHistogram::median() called on empty DHistogram" << endl;
    return 0.0;
  }

  double N = double(n() + nBelow + nAbove)/2;
  unsigned int i;

  double *binPtr = _contents;
  double  sum    = nBelow + *binPtr++;
  for (i = 0; (i < _size) && (sum < N); i++, sum += *binPtr++);

  assert(sum >= N);

  if (bin)
    *bin = i;

  return binStart(i + 1) - double(sum - N)/_contents[i]*_binWidth;
}

double
DHistogram::majority(unsigned *bin) const
{
  if (!_size) {
    cerr << "Warning! DHistogram::majority() called on empty DHistogram" << endl;
    return 0.0;
  }

  unsigned i;
  SimpleArray<double>::max(&i);

  if (bin)
    *bin = i;

  return binCenter(i);
}

// function which is used to find the threshold of a histogram
// the function takes the histogram (greylevel), 
double
DHistogram::biModalThreshold() const
{
  if (!_size) {
    cerr << "Warning! DHistogram::biModalThreshold() called on empty Histogram" << endl;
    return 0.0;
  }

  double N = n();
  double   meanValue = mean();
  double   varMax = 0.0;
  unsigned i = 0;
  double   offset = _binWidth/2;

  double zeroMoment0  = double(_contents[0])/N;
  double firstMoment0 = (_valueToBinMap.reverse(0) + offset)*_contents[0]/N;

  for (unsigned k = 1; k < _size; k++) {
    // 0th and 1st cumulative moments of the histogram up to the kth level
    double zeroMoment1  = zeroMoment0 + double(_contents[k])/N;
    double firstMoment1 = 
      firstMoment0 + (_valueToBinMap.reverse(k) + offset)*_contents[k]/N;

    if ((zeroMoment1 > 0.0) && (zeroMoment1 < 1.0)) {
      double var = 
	SQR(meanValue*zeroMoment1 - firstMoment1)/(zeroMoment1*(1 - zeroMoment1));
      if (var > varMax) {
	i = k;
	varMax = var;
      }
    }

    zeroMoment0  = zeroMoment1;
    firstMoment0 = firstMoment1;
  }
  
  return _valueToBinMap.reverse(i) + offset;
}


// function which is used to find the threshold of a histogram
// based on minimum variance threshold obtained from Haralick and Shapiro book

double
DHistogram::varianceThreshold() const
{
  if (!_size) {
    cerr << "Warning! DHistogram::varianceThreshold() called on empty Histogram" << endl;
    return 0.0;
  }

  double var = 0;
  double varMin = MAXDOUBLE;
  double N = n();
  double   mean1;
  double   mean2;
  double   var1;
  double   var2;
  double   q1;
  double   q2;
  double   offset = _binWidth/2;
  
  int i = 0;
  int t = 0;
  
  for (unsigned k = 0 ; k < _size ; k++) {
     mean1 = mean2 = var1 = var2 = q1 = q2 = 0;
     for ( i=0 ;  i<=k ; i++)   q1 += double(_contents[i])/N;
     for ( i=k+1; i<_size ; i++) q2 += double(_contents[i])/N;

     if (q1 != 0) {
       for ( i=0; i<=k ; i++) 
	 mean1 += _valueToBinMap.reverse(i) * (double(_contents[i])/N) / q1; 
       for ( i=0; i<=k ; i++)
	 var1  += SQR(_valueToBinMap.reverse(i)-mean1)*(double(_contents[i])/N) / q1;
     }

     if (q2 != 0){
       for ( i=k+1; i<_size ; i++)
	 mean2 += _valueToBinMap.reverse(i) * (double(_contents[i])/N) / q2; 
       for ( i=k+1; i<_size ; i++)
	 var2  += SQR(_valueToBinMap.reverse(i)-mean2)*(double(_contents[i])/N) / q2;
     }
     
     var = q1 * var1 +  q2 * var2;
     
     if (var < varMin) {
	t = k;
	varMin = var;
     }
  }
  
  return _valueToBinMap.reverse(t) + offset;
}



// function which is used to find the threshold of a histogram
// based on minimum kullback function obtained from Haralick and Shapiro book

double
DHistogram::kullbackThreshold() const
{
  if (!_size) {
    cerr << "Warning! DHistogram::kullbackThreshold() called on empty Histogram" << endl;
    return 0.0;
  }

  double H = 0;
  double HMin = MAXDOUBLE;
  double N = n();
  double   mean1;
  double   mean2;
  double   var1;
  double   var2;
  double   q1;
  double   q2;
  double   offset = _binWidth/2;
  double    pi=4*atan(1.0);	/* constant, 3.14159.... */
  
  int i = 0;
  int t = 0;
  
  for (unsigned k = 0 ; k < _size ; k++) {
     mean1 = mean2 = var1 = var2 = q1 = q2 = 0;
     for ( i=0 ;  i<=k ; i++)   q1 += double(_contents[i])/N;
     for ( i=k+1; i<_size ; i++) q2 += double(_contents[i])/N;

     if (q1 != 0){
       for ( i=0; i<=k ; i++)
	 mean1 += _valueToBinMap.reverse(i) * (double(_contents[i])/N) / q1; 
       for ( i=0; i<=k ; i++) 
	 var1  += SQR(_valueToBinMap.reverse(i)-mean1)*(double(_contents[i])/N) / q1;
     }

     if (q2 != 0){
       for ( i=k+1; i<_size ; i++)
	 mean2 += _valueToBinMap.reverse(i) * (double(_contents[i])/N) / q2; 
       for ( i=k+1; i<_size ; i++) 
	 var2  += SQR(_valueToBinMap.reverse(i)-mean2)*(double(_contents[i])/N) / q2;
     }
     
     
     if (var1 > 0 && var2 > 0 && q1 > 0 && q2 > 0) {
       H  = (1+log10(2*pi))/2;
       H += ( - q1*log10(q1) - (q2*log10(q2)) );
       H += (q1*log10(var1) + q2*log10(var2) )/2;
     }        
     else if (var1 == 0 && var2 == 0)  H = -MAXDOUBLE;
     else  H =  MAXDOUBLE; 
     
     if (H < HMin) {
       t = k;
       HMin = H;
     }
  }
  
  return _valueToBinMap.reverse(t) + offset;
}

double
DHistogram::entropy() const
{
  double *binPtr = _contents;
  double  N = n();
  double    entropy = 0;
  double    log2 = ::log(2.0);
  for (unsigned i = _size; i; i--, binPtr++) {
    double p = double(*binPtr) / N;
    if (p)
      entropy += p * ::log(p)/log2;
  }

  return -entropy;
}

DblArray
DHistogram::pdf() const
{
  return asDblArray(*this)/n();
}

DblArray
DHistogram::cdf() const
{
  return pdf().cumSum();
}

//
// Other operators
//

DHistogram&
DHistogram::operator += (const DHistogram& hist)
{
  //printHeadAndTail(cout);
  //hist.printHeadAndTail(cout);

  if ((_min != hist._min) || (_max != hist._max) || (_size != hist._size))
    cerr << "DHistogram::operator +=() : cannot merge histograms with different mappings"
	 << endl;
  else
    this->SimpleArray<double>::operator += (hist);

  //printHeadAndTail(cout);

  return *this;
}

LUT<double>
DHistogram::equalize(const DHistogram& hist) const
{
  DblArray wx(pdf().cumSum());
  DblArray wy(hist.pdf().cumSum());

  LUT<double> lut(_size);

  int len = hist.nBins() - 1;

  int j = 0;
  for (unsigned i = 0; i < _size; i++) {
    double wxi = wx[i];
    while ((wy[j] < wxi) && (j < len))
      j++;
    
	lut.add(binCenter(i), hist.binCenter(j));
  }

  return lut;
}

//
// I/O
//

ostream&
DHistogram::printHeadAndTail(ostream& os, unsigned n) const
{
  os << "Hist size: " << _size << " min: " << _min << " max: " << _max << endl;
  n = MIN(n, _size/2);
  SimpleArray<double> temp((*this)(n));
  temp.print(os);
  if (_size > 2*n)
    os << " ... ";
  else
    os << " ";
  temp = (*this)(_size - n - 1, _size - 1);
  temp.print(os) << endl;

  return os;
}

#ifdef HAVE_MATLAB
Boolean
DHistogram::saveMatlab(const char *fileName, const char *binVarName, 
		      const char *countVarName, const char *option) const
{
  if ((option[0] != 'u') && (option[0] != 'w')) {
    cerr<<"Incorrect file write option, use u for update or w to";
    cerr<<" create a new file or overwrite the existing one."<<endl;
    return FALSE;
  }
   
  //open a Matfile for the object, the selected file option is update
  //This creates a new file if necessary or updates an existing one without
  //erasing its previous contents.
 
  MATFile *outFile = matOpen((char *) fileName, (char *) option); 
  if (!outFile) {
    cerr << "Couldn't open file" << fileName << endl;
    return FALSE;
  }
   
  if (matPutFull(outFile, (char *) binVarName, _size, 1, binCenters(), NULL)) {
    cerr << "Error in writing bins array to file"<<endl;
    matClose(outFile);
    return FALSE;
  }
  
  if (matPutFull(outFile, (char *) countVarName, _size, 1, asDblArray(*this), NULL)) {
    cerr<<"Error in writing counts array to file"<<endl;
    matClose(outFile);
    return FALSE;
  }
  
  if (matClose(outFile)==EOF) {
    cerr<<"Error in file write"<<endl;
    return FALSE;
  }

  return TRUE;
}
#endif

//
// Friends
//

ostream&
operator << (ostream& os, const DHistogram& hist)
{
  if (hist._size)
    hist.print(os);

  return os;
}

SimpleArray<double> 
asDblArray(const DHistogram& hist)
{
  DblArray bins(hist._size);

  const double *histPtr = hist._contents;
  double         *binPtr  = bins.contents();

  for (unsigned i = hist._size; i; i--)
    *binPtr++ = (double) *histPtr++;

  return bins;
}

//
// Private functions
//






