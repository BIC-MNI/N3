/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, Alex P. Zijdenbos, 
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
$RCSfile: Spline.cc,v $
$Revision: 1.1 $
$Author: claude $
$Date: 2010-12-09 19:35:01 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "Spline.h"
#include "EBTKS/Matrix.h"
#include "EBTKS/MPoint.h"
#include "EBTKS/trivials.h"
#include <assert.h>

using namespace std;


/*******************
 * Spline base class
 *******************/

Boolean Spline::verbose = FALSE;

//
// Constructors/destructor
//

Spline::Spline(unsigned nDimensions)
  : _coef(0) 
{
  _nKnots = 0;
  _nDimensions = nDimensions;
  _initialize();
}

Spline::Spline(const FloatArray& xKnots) 
  : _knots(size(xKnots), 1), _coef(0) 
{ 
  _nKnots = size(xKnots);
  xKnots.resetIterator();

  for (unsigned i = 0; i < _nKnots; i++)
    _knots(i, 0) = xKnots++;

  _nDimensions = 1; 

  _initialize();
}

Spline::Spline(const FloatArray& xKnots, const FloatArray& yKnots) 
  : _knots(size(xKnots), 2), _coef(0) 
{ 
  _nKnots = size(xKnots);
  assert(_nKnots == size(yKnots));
  xKnots.resetIterator();
  yKnots.resetIterator();

  for (unsigned i = 0; i < _nKnots; i++) {
    _knots(i, 0) = xKnots++;
    _knots(i, 1) = yKnots++;
  }

  _nDimensions = 2; 

  _initialize();
}

Spline::Spline(const FloatArray& xKnots, const FloatArray& yKnots, 
	       const FloatArray& zKnots) 
  : _knots(size(xKnots), 3), _coef(0) 
{ 
  _nKnots = size(xKnots);
  assert((_nKnots == size(yKnots)) && (_nKnots == size(zKnots)));
  xKnots.resetIterator();
  yKnots.resetIterator();
  zKnots.resetIterator();

  for (unsigned i = 0; i < _nKnots; i++) {
    _knots(i, 0) = xKnots++;
    _knots(i, 1) = yKnots++;
    _knots(i, 2) = zKnots++;
  }

  _nDimensions = 3; 

  _initialize();
}

Spline::Spline(const FlMat& knots)
  : _knots(knots), _coef(0) 
{
  _nKnots = knots.getrows();
  _nDimensions = knots.getcols(); 
  assert(_nDimensions);

  _initialize();
}

Boolean
Spline::J(const DblMat& J)
{
  if ((J.getrows() != _J.getrows()) || (J.getcols() != _J.getcols())) {
    cerr << "Matrix J should be " << _J.getrows() << "x" << _J.getcols() << endl;
    return FALSE;
  }

  _J = J;

  return TRUE;
}

//
// Fitting functions
//

Boolean
Spline::clearDataPoints()
{
  assert(_nCoef);

  _fitted = FALSE;
  _nsamples = 0;

  _AtA = Zeros<double>(_nCoef, _nCoef);
  _AtF = Zeros<double>(_nCoef, 1);
  _Aj  = DblMat(1, _nCoef);

  return TRUE;
}

Boolean 
Spline::fit(const FlMat& coord, const DblArray& values)
{
  if (_nDimensions != coord.getcols()) {
    cerr << "Cannot fit a " << _nDimensions << "-dimensional Spline with "
	 << coord.getcols() << "-dimensional coordinates" << endl;
    return FALSE;
  }
  
  unsigned nPoints = coord.getrows();

  if (values.size() != nPoints) {
    cerr << "Spline::fit: coord and values do not match!" << endl;
    return FALSE;
  }

  clearDataPoints();

  for (unsigned i = 0; i < nPoints; i++)
    addDataPoint(coord[i], values[i]);

  return fit();
}

Boolean
Spline::fit()
{
  DblMat W(inv(_AtA));
  if (!W)
    return FALSE;

  W *= _AtF;

  _coef = array(W);
  assert(_nCoef == _coef.size());

  _fitted = TRUE;

  return TRUE;
}

void
Spline::saveState(char *filename) const
{
#ifdef HAVE_MATLAB
  _knots.saveMatlab(filename, "knots");
  _J.saveMatlab(filename, "J");
  _AtA.saveMatlab(filename, "AtA");
  _AtF.saveMatlab(filename, "AtF");
  _coef.saveMatlab(filename, "coef");
#else
  cerr << "Warning: saveState function not available in class Spline\n";
#endif 
}

#ifdef HAVE_INTERFACE
FlMat 
Spline::coord2FlMat(const MRegion& coord)
{
  FlMat coordMat(coord.nPoints(), 2);

  MRegionIterator pointsIt(coord);
  MPoint         *point;
  unsigned        i = 0;
  while (point = pointsIt++) {
    coordMat(i, 0) = point->x;
    coordMat(i, 1) = point->y;
    i++;
  }

  return coordMat;
}
#endif

FlMat 
Spline::coord2FlMat(const OrderedCltn& coord)
{
  FlMat coordMat(coord.size(), 3);

  ocIterator pointsIt(coord);
  MPoint3D  *point;
  unsigned   i = 0;
  while (point = (MPoint3D *) pointsIt++) {
    coordMat(i, 0) = point->x;
    coordMat(i, 1) = point->y;
    coordMat(i, 2) = point->z;
    i++;
  }

  return coordMat;
}

void
Spline::_initialize()
{
  _nCoef = 0;
  _fitted = FALSE; 
  _lambda = 0.0;
  _nsamples = 0;
  assert(_nDimensions);
  _tempPoint = new float[_nDimensions];
  assert(_tempPoint);

  _prune();

  if (verbose)
    cout << "Created " << _nDimensions << "-dimensional Spline" << endl;
}

void
Spline::_prune()
{
  if (_nKnots <= 1)
    return;

  Boolean knotsRemoved = FALSE;

  for (unsigned i = 0; i < _nKnots; i++)
    for (unsigned j = i + 1; j < _nKnots; j++)
      if (_knots.row(i) == _knots.row(j)) {
	cerr << "Warning! Identical knots supplied: eliminating knot " 
	     << j << "." << endl;
	knotsRemoved = TRUE;
	_knots.swapRows(j, _nKnots - 1);
	_nKnots--;
	j--;
      }

  if (knotsRemoved)
    _knots = _knots.rows(0, _nKnots - 1);
}

/****************
 * B-spline class
 ****************/

/*********************
 * Radial spline class
 *********************/

RSpline::RSpline(const FloatArray& xKnots, double scale)
  : Spline(xKnots)
{
  _coordScale  = scale;
  _coordScale2 = SQR(scale);
  _initialize();
}

RSpline::RSpline(const FloatArray& xKnots, const FloatArray& yKnots, double scale)
  : Spline(xKnots, yKnots)
{
  _coordScale  = scale;
  _coordScale2 = SQR(scale);
  _initialize();
}

RSpline::RSpline(const FloatArray& xKnots, const FloatArray& yKnots, 
		 const FloatArray& zKnots, double scale)
  : Spline(xKnots, yKnots, zKnots)
{
  _coordScale  = scale;
  _coordScale2 = SQR(scale);
  _initialize();
}

RSpline::RSpline(const FlMat& knots, double scale)
  : Spline(knots)
{
  _coordScale  = scale;
  _coordScale2 = SQR(scale);
  _initialize();
}

#ifdef HAVE_INTERFACE
RSpline::RSpline(const MRegion& knots, double scale)
  : Spline(coord2FlMat(knots))
{
  _coordScale  = scale;
  _coordScale2 = SQR(scale);
  _initialize();
}
#endif

RSpline::RSpline(const OrderedCltn& knots, double scale)
  : Spline(coord2FlMat(knots))
{
  _coordScale  = scale;
  _coordScale2 = SQR(scale);
  _initialize();
}

Boolean 
RSpline::addDataPoint(const float *point, double value)
{
  if (!point)
    return FALSE;

  _nsamples++;

  double *AjPtr = (double *) _Aj[0];
  for (unsigned spline = 0; spline < _nKnots; spline++)
    *AjPtr++ = _radialFunction(_r2(point, _knots[spline]));

  _AtA += _Aj.transposeXself();
  _AtF += _Aj * value;

  return TRUE;
}

double
RSpline::operator () (const float *point) const
{
  double result = 0.0;

  if (!_fitted)
    return result;

  if (!point)
    point = _tempPoint;
  
  const double *coefPtr = _coef.contents();

  for (unsigned spline = 0; spline < _nKnots; spline++)
    result += (*coefPtr++) * _radialFunction(_r2(point, _knots[spline]));

  return result;
}

void
RSpline::_initialize()
{
  switch(_nDimensions) {
  case 1:
    _radialFunction = rToTheThird;
    break;
  case 2:
    _radialFunction = rSquareLogR;
    break;
  default:
    _radialFunction = r;
    break;
  }

  _nCoef = _nKnots;

  _J = Eye<double>(_nCoef, _nCoef); 

  clearDataPoints();

  if (verbose)
    cout << "Created RSpline with " << _nKnots << " basis functions" << endl;
}

/*************************
 * Thin-plate spline class
 *************************/

//
// Constructors/destructor
//

TPSpline::TPSpline(const FloatArray& xKnots, double scale)
  : RSpline(xKnots, scale) 
{ 
  _initialize();
}

TPSpline::TPSpline(const FloatArray& xKnots, const FloatArray& yKnots, double scale)
  : RSpline(xKnots, yKnots, scale) 
{ 
  _initialize();
}

TPSpline::TPSpline(const FloatArray& xKnots, const FloatArray& yKnots, 
		   const FloatArray& zKnots, double scale)
  : RSpline(xKnots, yKnots, zKnots, scale) 
{ 
  _initialize();
}

TPSpline::TPSpline(const FlMat& knots, double scale)
  : RSpline(knots, scale)
{ 
  _initialize();
}

#ifdef HAVE_INTERFACE
TPSpline::TPSpline(const MRegion& knots, double scale)
  : RSpline(knots, scale)
{ 
  _initialize();
}
#endif

TPSpline::TPSpline(const OrderedCltn& knots, double scale)
  : RSpline(knots, scale)
{ 
  _initialize();
}

TPSpline::~TPSpline() 
{ 
}

Boolean
TPSpline::clearDataPoints()
{
  if (!_nKnots || (_nKnots <= _nDimensions + 1)) {
    cerr << "TPSpline::clearDataPoints: #knots (" << _nKnots 
         << ") must be larger than " << _nDimensions + 1 << endl;
    return FALSE;
  }

  unsigned i, d;

  _fitted = FALSE;
  _nsamples = 0;

  _AtA = Zeros<double>(_nKnots, _nKnots);
  _AtF = Zeros<double>(_nKnots, 1);

  // Build Psi_inv matrix
  Ones<double> Psi_inv(_nDimensions + 1, _nDimensions + 1);
  for (d = 0; d < _nDimensions; d++) {
    for (i = 0; i <= _nDimensions; i++)
      Psi_inv(d + 1, i) = _knots(i, d);
  }

  // Replace knots in Psi_inv while matrix is singular
  unsigned knotToTest = 1;
  unsigned knot       = _nKnots - 1;
  Boolean  failure    = FALSE;

  while (!failure && (knotToTest <= _nDimensions)) {
    double D = det(Psi_inv(0, knotToTest, 0, knotToTest));

    if (verbose)
      cout << "D: " << D << endl;

    if (D)
      knotToTest++;
    else {
      if (verbose)
        cout << "Knot " << knotToTest << " failed. Trying knot " << knot << " ..." 
        << endl;

      for (d = 0; d < _nDimensions; d++) {
        swap(_knots(knotToTest, d), _knots(knot, d));
        Psi_inv(d + 1, knotToTest) = _knots(knotToTest, d);
      }

      failure = (--knot <= _nDimensions);
    }
  }

  if (failure) {
    cerr << "TPSpline::clearDataPoints: could not build matrix Psi_inv" << endl;
    return FALSE;
  }

  // Invert Psi_inv
  _Psi = inv(Psi_inv);

  // Build Psi_r matrix
  unsigned nKnotsRed = _nKnots - _nDimensions - 1;

  _Psi_r = Ones<double>(_nDimensions + 1, nKnotsRed);
  
  for (i = 0; i < nKnotsRed; i++)
    for (d = 0; d < _nDimensions; d++)
      _Psi_r(d + 1, i) = _knots(i + _nDimensions + 1, d);

  // Adding zeros to allow for _Aj -= _Ft*_Psi_r; in addDataPoints
  _Psi_r = (_Psi * _Psi_r).appendRight(Zeros<double>(_nDimensions+1, _nDimensions+1));

  _Ft = DblMat(1, _nDimensions + 1);
  _Aj = DblMat(1, _nKnots);

  return TRUE;
}

Boolean 
TPSpline::addDataPoint(const float *point, double value)
{
  if (!point)
    return FALSE;

  _nsamples++;

  unsigned spline;
  // Calculate F twiddle transposed
  double *FtPtr = (double *) _Ft[0];
  for (spline = 0; spline <= _nDimensions; spline++)
    *FtPtr++ = _radialFunction(_r2(point, _knots[spline]));

  // Calculate row j of A
  double *AjPtr = (double *) _Aj[0];
  for (; spline < _nKnots; spline++)
    *AjPtr++ = _radialFunction(_r2(point, _knots[spline]));

  *AjPtr++ = 1.0;
  for (unsigned d = 0; d < _nDimensions; d++)
    *AjPtr++ = _coordScale*(*point++);

  _Aj -= _Ft*_Psi_r;

  // Add Aj to AtA and AtF
  _AtA += _Aj.transposeXself();
  _AtF += _Aj * value;

  return TRUE;
}

Boolean 
TPSpline::fit()
{
  // Calculate C
  DblMat R(Ones<double>(_nDimensions + 1, _nKnots - _nDimensions - 1));
  for (unsigned spline = _nDimensions + 1, i = 0; spline < _nKnots; spline++, i++)
    for (unsigned d = 0; d < _nDimensions; d++)
      R(d + 1, i) = _knots(spline, d);

  DblMat J = _J.cols(0, _nDimensions) * _Psi * R + _J.cols(_nDimensions + 1, _nKnots-1);
  J = R.t() * _Psi.t() * J.rows(0, _nDimensions) + J.rows(_nDimensions + 1, _nKnots-1);
  J.pad(_nKnots, _nKnots, 0, 0);

  // Calculate W
  DblMat W(inv(_AtA + (_lambda*_nsamples)*J) * _AtF);

  _coef = array(-_Psi * R * W.rows(0, _nKnots - _nDimensions - 2));
  _coef.append(array(W));
  assert(_nCoef == _coef.size());

  _fitted = TRUE;

  return TRUE;
}

double
TPSpline::operator () (const float *point) const
{
  double result = 0.0;

  if (!_fitted)
    return result;

  if (!point)
    point = _tempPoint;

  const double *coefPtr = _coef.contents();

  for (unsigned spline = 0; spline < _nKnots; spline++)
    result += (*coefPtr++) * _radialFunction(_r2(point, _knots[spline]));

  result += *coefPtr++;
  for (unsigned d = 0; d < _nDimensions; d++)
    result += (*coefPtr++) * _coordScale*(*point++);

  return result;
}

void
TPSpline::_initialize()
{
  _nCoef = _nKnots + _nDimensions + 1;
  if(!clearDataPoints())
    cerr<<"ERROR: Can't initialize matrix!"<<std::endl;

  if (verbose)
    cout << "Created TPSpline with " << _nKnots << " basis functions" << endl;
}
