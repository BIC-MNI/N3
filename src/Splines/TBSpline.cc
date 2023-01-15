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
$RCSfile: TBSpline.cc,v $
$Revision: 1.1 $
$Author: claude $
$Date: 2010-12-09 19:35:01 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : TBSpline.c,v
@INPUT      : 
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: Tensor cubic B-splines in N dimensions with limited bending energy
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 21, 1996 (John G. Sled)
@MODIFIED   : Log: TBSpline.c,v 
 * Revision 1.4  1996/07/29  16:04:06  jgsled
 * last version before modification for inclusion in the AI source tree
 *
 * Revision 1.3  1996/07/09  15:13:59  jgsled
 * Working version
 *
 * Revision 1.2  1996/04/23  13:39:12  jgsled
 * Optimized addDataPoint function, modified to be compatible with changes to
 * Spline.h
 *
@COPYRIGHT  : 1996
---------------------------------------------------------------------------- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "TBSpline.h"

#undef min
#undef max

#ifdef DEBUG_TBSPLINE
#include <stdio.h>
#endif

using namespace std;

// For rounding errors.
#ifndef EPSILON
#define EPSILON ( 1.0e-14 )
#endif

//-----------------------------------------------------------------------
// LAPACK definitions
typedef long int _integer;
typedef float _real;
typedef double _doublereal;

extern "C" {
int dsysv_(char *uplo, _integer *n, _integer *nrhs, _doublereal
        *a, _integer *lda, _integer *ipiv, _doublereal *b, _integer *ldb,
        _doublereal *work, _integer *lwork, _integer *info);
}
//-----------------------------------------------------------------------
// static initializations
double TBSpline::_default_lambda = 0.1;

//-----------------------------------------------------------------------
// Methods for TBSpline Class
//-----------------------------------------------------------------------
// Constructor
TBSpline::TBSpline(const DblMat &domain, 
		   // an n by 2 matrix specifying a region in R^n
		   // upon which the spline is defined
		   // eg. [ x1 x2; y1 y2 ] 
		   double distance,  // the distance between adjacent knots 
		   double lambda,     // 0 < lambda < inf 
                                     // lambda is a limit on the
                                     // bending energy of the spline
                   int allocate_flag) :  // TRUE -> normal usage
                                         // FALSE -> cannot call add_data 
            Spline(domain.getrows())
{
  _thisTBSpline = this;
  _nDimensions = 0;

  if(domain.getcols() != 2 || distance <= 0 )
    return;  // fatal error

  _nDimensions = domain.getrows();
 
  _lambda = lambda;
  _distance = distance;
  _scale = 1.0/(distance*distance*distance);
  _nsamples = 0;
  _domain = domain;
  unsigned int i;
  // check that domain bounds are in correct order
  for(i = 0; i < _nDimensions; i++)
    {

      if(_domain(i,0) > _domain(i,1))
	{ // reverse order
	  double temp;
	  temp = _domain(i,0);
	  _domain(i,0) = _domain(i,1);
	  _domain(i,1) = temp;
	}
    }

  // determine number of basis functions in each dimension
  _n.newSize(_nDimensions);
  for(i = 0; i < _nDimensions; i++){
     _n[i] = (int) ceil( (_domain(i,1) - _domain(i,0))/(_distance * ( 1.0 + EPSILON ) )  ) + 3;
  }
   
  // create array of knot locations for each dimension
  _knots.resize(_n.max()+4,_nDimensions);
  for(i = 0; i < _nDimensions; i++)
    {
      double start = 0.5*(_domain(i,0) + _domain(i,1) - _distance*(_n[i]+3));
      for(int j = 0; j < _n[i]+4; j++)
	_knots(j,i) = start + _distance*j;
    }

  // location of knot zero
  _zero = _knots.row(3);

  // allocate (large) matrices
  _nProduct = 1;
  _four = 1;
  _smallN.newSize(4);
  for(i = 0; i < _nDimensions; i++)
  {
    _nProduct *= _n[i];
    _four *= 4;
    _smallN[i] = 4;
  }
  
  if(allocate_flag) {
    _AtA.resize(_nProduct,_nProduct);
    _AtF.resize(_nProduct,1);
  }
  _coef.newSize(_nProduct);
  _nCoef = _nProduct;
  _terms.resize(4,_nDimensions);
  _blockIndex.newSize(_nDimensions);
  _values.newSize(_four);
  _locations.newSize(_four);
  _dloc_i = new int[_four];
  _dloc_j = new int[_four];

  clearDataPoints();

}


// Public Methods of TBSpline

// remove all data points from fit
Boolean 
TBSpline::clearDataPoints(void)
{
  _fitted = FALSE;
  _haveData = FALSE;
  _nsamples = 0;
  if(_AtA.getrows() > 0) {
    _AtA.fill(0);
    _AtF.fill(0);
  }

  return TRUE;
}

#define for_each_dimension(i) for((i) = 0; (i) < _nDimensions; (i)++)

// incorporate new data point in fit
// returns FALSE if data is outside domain
Boolean
TBSpline::addDataPoint(const float *point, double value)
{
  int i,j,k,l;			// (bert) - changed from unsigned to signed.
  
  // check that point is within domain
  for_each_dimension(i)
    {
      if(point[i] > _domain(i,1) || point[i] < _domain(i,0))
	return(FALSE);
    }
    
  _fitted = FALSE;
  _haveData = TRUE;
  _nsamples++;

  // locate nearest knot location greater than or equal to point
  for_each_dimension(i)
    {
      _blockIndex[i] = (int) ceil((point[i] - _zero(i))/_distance) - 1;
      if(_blockIndex[i] < 0)
        _blockIndex[i] = 0;
      else if(_blockIndex[i] > _n[i] - 4)
	_blockIndex[i] = _n[i] - 4;
    }
	
  // compute terms used in tensor product
  for_each_dimension(i)
    {
      double temp;
      temp = _scale*cube(_knots(_blockIndex[i]+4,i)-point[i]);
      _terms(0,i) = temp;
      _terms(1,i) = 
        _scale*cube(_knots(_blockIndex[i]+5,i)-point[i]) - 4.0*temp;
      temp = _scale*cube(point[i]-_knots(_blockIndex[i]+3,i));
      _terms(2,i) = _scale*cube(point[i]-_knots(_blockIndex[i]+2,i)) - 
			   4.0*temp;
      _terms(3,i) = temp;
    }
    
  // compute tensor product
  TSubIndex subIndex(TIndex(_blockIndex,_n), _smallN);  // index to small tensor

  for(j = 0; j < _four; j++) // for each value in small tensor
    {
      _locations[j] = subIndex.flat(); // record index into large tensor
      _values[j] = 1;
      for_each_dimension(i)   // compute product
	_values[j] *= _terms(subIndex[i],i);

      ++subIndex;
    }

  // create list of differences between locations for fast indexing
  for(i = _four-2; i >= 0; i--)
    {
      _dloc_i[i] = _locations[i+1] - _locations[i];
      _dloc_j[i] = _dloc_i[i]*_nProduct;
    }
  

  // add new data to AtA and AtF
  double incr, value_k;
  double *AtAel, *AtFel;         // pointers to raw data
  double *U, *L;    // pointers to upper and lower triangles
  double *D;          // pointer to diagonal
  double *F;        // pointer to AtF

  AtAel = (double *) *(_AtA.getEl());
  AtFel = (double *) *(_AtF.getEl());
  // start in bottom right corner of block within AtA
  D = AtAel + _locations[(unsigned int)_four-1]*_nProduct + 
    _locations[(unsigned int)_four-1];
  F = AtFel + _locations[(unsigned int)_four-1];
  for(k = _four-1; k >= 0; )
    {
      L = U = D;
      value_k = _values[k];
      *F += value*value_k;
      for(l = k-1; l >= 0; l--)  // only do half since AtA is symmetric
	{
	  incr = value_k*_values[l];
	  U -= _dloc_i[l];
	  L -= _dloc_j[l];
	  *U += incr;
	  *L += incr;
	}
      *D += value_k*value_k; // do diagonal elements separately

      --k;
      if( k >= 0 ) {
        D -= _dloc_i[k] + _dloc_j[k];
        F -= _dloc_i[k];
      }
    }
  return TRUE;
}

// fit splines to data
//  returns true if sucessful
Boolean
TBSpline::fit(void)
{
  if(_haveData == FALSE)
    return FALSE;

  int info;
  DblMat A;
  bendingEnergyTensor(_n, A);  // create bending matrix on fly to save memory
  A *= ( _lambda * _nsamples );
  A += _AtA;

  _coef = solveSymmetricSystem(A, _AtF,&info).array();

  _fitted = (info == 0);
  return (_fitted);
}

// evaluate tensor B-Spline at specified point
// Note: evaluates to zero outside domain
double 
TBSpline::operator () (const float *point) const
{
  unsigned int i;
  
  if(_fitted == FALSE)
    return 0;

  // check that point is within domain
  for_each_dimension(i)
    {
      if(point[i] > _domain(i,1) || point[i] < _domain(i,0))
	return 0.0;  // outside domain
    }

  // locate nearest knot location greater than or equal to point
  for_each_dimension(i)
    {
      _thisTBSpline->_blockIndex[i] = (int) ceil((point[i] - _zero(i))/_distance)-1;
      if(_blockIndex[i] < 0)
	_thisTBSpline->_blockIndex[i] = 0;
      else if(_blockIndex[i] > _n[i] - 4)
	_thisTBSpline->_blockIndex[i] = _n[i] - 4;
    }
	
  // compute terms used in tensor product
  for_each_dimension(i)
    {
      double temp;
      temp = _scale*cube(_knots(_blockIndex[i]+4,i)-point[i]);
      _thisTBSpline->_terms(0,i) = temp;
      _thisTBSpline->_terms(1,i) = 
        _scale*cube(_knots(_blockIndex[i]+5,i)-point[i]) - 4.0*temp;
      temp = _scale*cube(point[i]-_knots(_blockIndex[i]+3,i));
      _thisTBSpline->_terms(2,i) = 
        _scale*cube(point[i]-_knots(_blockIndex[i]+2,i)) - 4.0*temp;
      _thisTBSpline->_terms(3,i) = temp;
    }

  // sum weighted terms of tensor product
  double sum = 0;
  TSubIndex subIndex(TIndex(_blockIndex,_n), _smallN);  // index to small tensor

  for(int j = 0; j < _four; j++) // for each value in small tensor
    {
     double product = 1;
      for_each_dimension(i)   // compute product
	product *= _terms(subIndex[i],i);

      sum += product*_coef[(unsigned int)subIndex.flat()];
      ++subIndex;
    }
  return sum;
}

// save internal state for analysis
void
TBSpline::saveState(char *filename) const
{
#ifdef HAVE_MATLAB
  _knots.saveMatlab(filename, "knots");
  //  _J.saveMatlab(filename, "J");
  _AtA.saveMatlab(filename, "AtA");
  _AtF.saveMatlab(filename, "AtF");
  DblMat coef(_coef);
  coef.saveMatlab(filename, "coef");
  _domain.saveMatlab(filename, "domain");
#else
  cerr << "Warning: saveState function not available in class Spline\n";
#endif 
}

#undef for_each_dimension


//-----------------------------------------------------------------------
// private and protected methods of TBSpline

// compute N dimensional bending energy matrix
//
// matrix is in block form with n[0] by n[0] blocks of 
//  n[1] by n[1] blocks of ...
//
void 
TBSpline::bendingEnergyTensor(const IntArray &n, DblMat &J)
{
  int nDimensions = n.size();
  int nProduct = 1;
  IntArray step(nDimensions);  
  unsigned int i,j,k,l,m;

  // deal with special case first
  if(nDimensions == 1)
    {
      J = bendingEnergy(n[(unsigned int)0], 2);
      return;
    }

  // calculate step sizes for indexing an array as a tensor
  for(i = 0; i < nDimensions; i++)
    {
      nProduct *= n[i];
      step[i] = 1;
      for(j = i+1; j < nDimensions; j++)
	step[i] *= n[j];
    }

  // compute 1-D bending energy matrices
  DblMat (*bendingMatrix)[3] = new DblMat[nDimensions][3];  
  for(i = 0; i < nDimensions; i++)
    for(j = 0; j < 3; j++)
      bendingMatrix[i][j] = bendingEnergy(n[i],j);
  
  // allocate large matrix 
  J.resize(nProduct,nProduct);
  J.fill(0.0);

  //  form kronecker products
  //  eg. in 3D   x"yz + xy"z + xyz" + 2x'y'z + 2xy'z' + 2x'yz'  
  //
  
  TIndex ti(n), tj(n);  // indices into tensor
  IntArray derivative(nDimensions);
  
  // do terms with second derivatives
  for(k = 0; k < nDimensions; k++) 
    {
      derivative.clear(0);
      derivative[k] = 2;
      for(i = 0; i < nProduct; i++, ++ti)
	for(j = 0; j < nProduct; j++, ++tj)
	  {
	    double product = 1.0;
	    for(m = 0; m < nDimensions; m++)
	      product *= bendingMatrix[m][derivative[m]](ti[m],tj[m]);
	    J(ti.flat(),tj.flat()) += product;
	  }
    } 

  ti.reset();
  tj.reset();
  // do terms with first derivatives
  for(k = 0; k < nDimensions; k++)
    for(l = k+1; l < nDimensions; l++)
    {
      derivative.clear(0);
      derivative[k] = 1;
      derivative[l] = 1;
      for(i = 0; i < nProduct; i++, ++ti)
	for(j = 0; j < nProduct; j++, ++tj)
	  {
	    double product = 1.0;
	    for(m = 0; m < nDimensions; m++)
	      product *= bendingMatrix[m][derivative[m]](ti[m],tj[m]);
	    J(ti.flat(),tj.flat()) += 2.0*product;
	  }
    } 

  delete [] bendingMatrix;
}


// compute 1-D bending energy matrices
// 
// returns a size by size matrix formed with order'th derivatives
//
DblMat 
TBSpline::bendingEnergy(int size, int order)
{
  int i,j,k,l;  // (vivek) changed all from unsigned to signed
  int interval, offset, region;
  const int spline = 4;  // ie  splines are cubic

  if(size < 4 || order < 0 || order > 2)  
    // bending energy not defined for size < 4
    return(*(new DblMat));


  // standardized cubic B-spline defined on [-4 0] with knots at
  //  -4 -3 -2 -1 0
  // each unit segement is written as a cubic defined on [0 1]
  //  eg. -x^3 + 3x^2 -3x + 1
  double B[spline][spline] = { { -1, 3, -3, 1}, { 3, -6, 0, 4}, 
			       { -3, 3,  3, 1}, { 1,  0, 0, 0}};
  // take order'th derivative of B
  for(i = 0; i < order; i++)
    for(j = 0; j < spline; j++)
      {
	for(k = 1; k < spline; k++) 
	  {
	    B[j][spline-k] = k*B[j][spline-k-1];
	  }
	B[j][0] = 0.0;
      }

  // compute product integral for each pair of segments
  // actually we only compute upper triangle since D will be symmetric
  DblArray C(2*spline-1);
  unsigned int sizeC = C.size();
  DblMat D(spline,spline);
  for(i = 0; i < spline; i++)
    for(j = i; j < spline; j++)
      {
	// convolve each pair of polynomials
	for(k = 0; k < sizeC; k++)
	  {
	    C[k] = 0;
	    for(l = 0; l < MIN(k+1,sizeC-k); l++)
	      C[k] += B[i][MAX(0,k-spline+1)+l]  //(vivek) - removed cast of index
		*B[j][spline-1-l-MAX(0,spline-1-k)];
	  }
	// evalute integral on [0,1]
	D(i,j) = 0;
	for(k = 0; k < sizeC; k++)
	  D(i,j) += C[k]/double((unsigned int)sizeC-k);
      }

  // define 6 regions:  [-1,0] [-2,0] [-3,0] [-4,0] [-2,-1] [-3,-1]
  // splines can be shifted with respect to each by 0 1 2 or 3
  // Compute product interval on each interval for each offset
  DblMat integral(6,4);  // regions by offset

  // start by computing integrals for first four regions
  for(offset = 0; offset < 4; offset++)
    {
      for(region = 0; region < 4 - offset; region++)
	{
	  integral(region,offset) = 0;
	    for(interval = 0; interval <= region; interval++)
	      integral(region,offset) +=  D(interval,interval+offset);
	}
      // last few regions are same since large offset is equivalent to a
      //   small region
      for(region = 4-offset; region < 4; region++)
	{
	  integral(region,offset) = integral(region-1,offset);
	}
    }

  // compute integrals for regions five and six separately
  integral(4,0) = D(1,1);
  integral(4,1) = D(1,2);
  integral(4,2) = D(1,3);
  integral(4,3) = 0;
  
  integral(5,0) = D(1,1) + D(2,2);
  integral(5,1) = D(1,2) + D(2,3);
  integral(5,2) = D(1,3);
  integral(5,3) = 0;

  // form integrals into bending energy matrix
  DblMat energy(size,size);

  // set elements common to matrices of all sizes first
  for(i = 0; i < 3; i++)
    {
      energy(0,i) = integral(0,i);
      energy(i,0) = integral(0,i);
      energy(size-1,size-i-1) = integral(0,i);
      energy(size-i-1,size-1) = integral(0,i);
    }
  
  if(size == 4) // deal with special cases
    {
      energy(1,1) = integral(4,0);
      energy(size-2,size-2) = integral(4,0);
      energy(1,2) = integral(4,1);
      energy(2,1) = integral(4,1);
      energy(3,0) = integral(3,3);
      energy(0,3) = integral(3,3);
    }
  else if(size == 5)
    {
      energy(1,1) = integral(1,0);
      energy(size-2,size-2) = integral(1,0);
      energy(2,2) = integral(5,0);
      energy(2,1) = integral(1,1);
      energy(1,2) = integral(1,1);
      energy(size-3,size-2) = integral(1,1);
      energy(size-2,size-3) = integral(1,1);
      energy(1,3) = integral(3,2);
      energy(3,1) = integral(3,2);

      for(i = 0; i < size-3; i++)
	{
	  energy(i,i+3) = integral(3,3);
	  energy(i+3,i) = integral(3,3);
	}
    }
  else  // n > 5  (general case)
    { 
      // fill in bulk region
      for(j = 0; j < 4; j++) 
	for(i = 3-j; i < size-3; i++)
	  {
	    energy(i,i+j) = integral(3,j);
	    energy(i+j,i) = integral(3,j);
	  }

      // fill in boundary region
      energy(2,2) = integral(2,0);
      energy(size-3,size-3) = integral(2,0);
      energy(1,2) = integral(1,1);
      energy(2,1) = integral(1,1);
      energy(size-2,size-3) = integral(1,1);
      energy(size-3,size-2) = integral(1,1);
      energy(1,1) = integral(1,0);
      energy(size-2,size-2) = integral(1,0);
    }
  return(energy);
}


// solve a system of equations A x = b for x where A is symmetric
// 
// returns x
// info == 0 indicates computation was sucessful
//
// Note: A is undefined after function call
DblMat
TBSpline::solveSymmetricSystem(DblMat &A, DblMat b, int *info)
{
  if(A.getrows() != A.getcols() || A.getrows() != b.getrows() ||
     b.getcols() != 1)
    { 
      DblMat x;
      *info = -1;
      return x;
    }
  
  _integer n = A.getcols();
  _integer nrhs = 1;
  _integer lda = n;
  _integer ldb = n;
  _integer *ipiv = new _integer[n];
  _doublereal work;
  _integer lwork = 1;
  _integer w_info;

  dsysv_("U", &n, &nrhs, (_doublereal *) *A.getEl(), &lda, ipiv, 
	 (_doublereal *) *b.getEl(), &ldb,
	 &work, &lwork, &w_info);

  delete [] ipiv;
  *info = (int) w_info;
  return(b);
}


//-------------------------------------------------------------------------
// Methods for TIndex Class

// Constructors
TIndex::TIndex(const IntArray &n)
{
  init(n);
}

TIndex::TIndex(const IntArray &index, const IntArray &n)
{
  unsigned int i;
  init(n);

  _index = index;
  _flat = _index[(unsigned int)0];
  for(i = 1; i < _nDimensions; i++)
    _flat = _flat*_n[i]+_index[i];
}

// initialization of variables common to both constructors
void 
TIndex::init(const IntArray &n)
{
  unsigned int i, j;
  _nDimensions = n.size();
  _n = n;
  _flat = 0;
  _index.newSize(_nDimensions);
  _index.clear(0);

  // calculate step sizes for indexing an array as a tensor
  _step.newSize(_nDimensions);  
  for(i = 0; i < _nDimensions; i++)
    {
      _step[i] = 1;
      for(j = i+1; j < _nDimensions; j++)
	_step[i] *= _n[j];
    }
}

void 
TIndex::operator ++ (void)
{
  for(int i = _nDimensions-1; i >= 0; i--) // (bert) changed 'i' to signed.
    {
      if(_index[i] < _n[i]-1)
	{
	  _index[i]++;
	  _flat += _step[i]; 
	  break;
	}
      else 
	{
	  _index[i] = 0;
	  _flat -= (_n[i]-1)*_step[i];
	}
    }
}


//-------------------------------------------------------------------------
// Methods for TSubIndex Class


// Note: small tensor is assumed to fit within large tensor
//    when placed at location index  
TSubIndex::TSubIndex(const TIndex &index, const IntArray &smallN) 
  : TIndex(index)
{
  init(smallN);
}

TSubIndex::TSubIndex(const TIndex &index, const IntArray &smallIndex,
		     const IntArray &smallN) :
   TIndex(index)
{
  init(smallN);
  _smallIndex = smallIndex;
}

// initialization of variables common to both constructors
void
TSubIndex::init(const IntArray &smallN)
{
  _smallN = smallN;
  _smallIndex.newSize(_nDimensions);
  _smallIndex.clear(0);
  _smallStep.newSize(_nDimensions);
  unsigned int i, j;

  // calculate step sizes for indexing an array as a tensor
  _smallStep.newSize(_nDimensions);  
  for(i = 0; i < _nDimensions; i++)
    {
      _smallStep[i] = 1;
      for(j = i+1; j < _nDimensions; j++)
	_smallStep *= _smallN[j];
    }
}

void 
TSubIndex::operator ++ ()
{
  int i;			// (bert) - changed from unsigned to signed.

  // increment indices
  for(i = _nDimensions-1; i >= 0 ; i--)
    {
      if(_smallIndex[i] < _smallN[i]-1)
	{
	  _smallIndex[i]++;   // increment index into small tensor
	  _index[i]++;        // increment index into large tensor
	  _flat += _step[i];  
	  break;
	}
      else 
	{
	  _smallIndex[(unsigned int)i] = 0;
	  _index[(unsigned int)i] -= _smallN[(unsigned int)i] - 1;
	  _flat -= (_smallN[(unsigned int)i] - 1)*_step[(unsigned int)i];
	}
    }
}

//------------------------------------------------------------------------
//  code for testing TBSpline 

#ifdef DEBUG_TBSPLINE
void 
test_index(void)
{
  int nvalues[2] = {2,3};
  IntArray n(nvalues,2); 
  TIndex index(n);
  int i;

  printf("Testing TIndex: size = {2 , 3 }\n");
  printf("%d %d: %d\n", index[0], index[1], index.flat());
  for(i = 0; i < 4; i++)
    ++index;
  printf("%d %d: %d\n", index[0], index[1], index.flat());

  int Nvalues[2] = {7,9};
  int atvalues[2] = {2,5};
  IntArray N(Nvalues,2);
  IntArray at(atvalues,2);
  TIndex bigIndex(at, N);
  TSubIndex subIndex(bigIndex, n);

  printf("Testing TSubIndex: size = {2 , 3 } with size { 7, 9 } at { 2, 5}\n");
  printf("%d %d: %d\n", subIndex[0], subIndex[1], subIndex.flat());
  ++subIndex;
  printf("%d %d: %d\n", subIndex[0], subIndex[1], subIndex.flat());
  for(i = 0; i < 3; i++)
    ++subIndex;
  printf("%d %d: %d\n", subIndex[0], subIndex[1], subIndex.flat());
}

void 
TBSpline::test_J(void)
{
  printf("Testing Bending Energy Matrix funtions\n");
  
  for(int j = 0; j < 3; j++)
    for(int i = 4; i < 8; i++)
      {
	DblMat J = bendingEnergy(i, j);
	printf("Size %d, derivative %d\n", i, j);
	cout << J << endl;
      }

  printf("Bending Energy Tensor: size = {4 , 5}\n");
  DblMat J;
  int nvalues[2] = { 4, 5};
  IntArray n(nvalues,2); 
  bendingEnergyTensor(n, J);
  cout << J << endl;

#ifdef HAVE_MATLAB
  J.saveMatlab("bending.mat","J");
#endif 
}

void 
TBSpline::test_solver(void)
{
  DblMat A(10,10);
  DblMat b(10,1);
  int i, j;
#ifdef HAVE_MATLAB
  int info;			// (bert) - made conditional to avoid warning
#endif // HAVE_MATLAB

  for(i = 0; i < 10; i++)
    {
      for(j = i; j < 10; j++)
	{
	  A(i,j) = i*i-j;
	  A(j,i) = A(i,j);
	}
      b(i,0) = i;
    }
#ifdef HAVE_MATLAB
  A.saveMatlab("solver.mat","A");
  b.saveMatlab("solver.mat","b");

  DblMat x = solveSymmetricSystem(A,b,&info);
  x.saveMatlab("solver.mat","x");
#endif
}

//-----------------------------------------------------------------------
// Methods for TBSplineVolume Class

//-----------------------------------------------------------------------
// Constructors
TBSplineVolume::TBSplineVolume(const DblMat &domain, 
	                  // an n by 2 matrix specifying a region in R^n
	                  // upon which the spline is defined
	                  // eg. [ x1 x2; y1 y2 ] 
                 const double start[VDIM],  // center of (0,0,0)th voxel 
                 const double step[VDIM],   // voxel separations
                 const int count[VDIM],     // number of voxels along each axis
                 double distance,  // the distance between adjacent knots 
                 double lambda, // 0 < lambda < inf 
                                             // lambda is a limit on the
                                             // bending energy of the spline
                 int allocate_flag) : // TRUE -> normal usage
                                             // FALSE -> cannot call add_data
  TBSpline(domain, distance, lambda, allocate_flag)
{
  createLookup(start, step, count);
}


TBSplineVolume::TBSplineVolume(const double start[VDIM],  
                               // center of (0,0,0)th voxel 
                 const double step[VDIM],   // voxel separations
                 const int count[VDIM],     // number of voxels along each axis
                 double distance,  // the distance between adjacent knots 
                 double lambda, // 0 < lambda < inf 
                                             // lambda is a limit on the
                                             // bending energy of the spline
                 int allocate_flag) : // TRUE -> normal usage
                                             // FALSE -> cannot call add_data
  TBSpline(computeDomain(start,step,count), distance, lambda, allocate_flag)
{
  createLookup(start, step, count);
}

// cache the individual spline evaluations for each voxel
void
TBSplineVolume::createLookup(const double start[VDIM],  
                             // center of (0,0,0)th voxel 
                             const double step[VDIM],   // voxel separations
                             const int count[VDIM])
{
  int i, j, k;
  int blockIndex;
  double x;
  double *pSpline;   // ptr into _1dSpline
  unsigned short *pOffset;      // ptr into offset
  unsigned short offsetStep, currentOffset, first, last;
  double temp;

  _sizes = IntArray(count, VDIM);

  // allocate space for cached values
  for(i = 0; i < VDIM; i++) {
    _1dSpline[i] = new double[count[i]*4];
    _offset[i] = new unsigned short[count[i]*4];
  }

  // evaluate 1D B splines at voxel centers
  for(i = 0; i < VDIM; i++)
    {
      pSpline = _1dSpline[i];
      pOffset = _offset[i];
      offsetStep = 1;
      for(j = 0; j < (VDIM - i - 1); j++)
        offsetStep *= _n[(unsigned int)VDIM-j-1];

      // skip any voxel outside the domain
      if(step[i] > 0) {
        first = (_domain(i,0) > start[i]) ?      
          // this may create a round off problem
          (unsigned short) ceil((_domain(i,0) - start[i])/step[i]) : 0;
        // compute index of last voxel within domain plus one
        last = (_domain(i,1) < (start[i]+(count[i]-1)*step[i])) ?
          (unsigned short) floor((_domain(i,1) - start[i])/step[i])+1 : count[i];      
      }
      else {
        first = (_domain(i,1) < start[i]) ?      
          // this may create a round off problem
          (unsigned short) ceil((_domain(i,1) - start[i])/step[i]) : 0;
        // compute index of last voxel within domain plus one
        last = (_domain(i,0) > (start[i]+(count[i]-1)*step[i])) ?
          (unsigned short) floor((_domain(i,0) - start[i])/step[i])+1 : count[i];
      }


      // fill outside domain with zeros
      for(k = 0; k < first*4; k++) {
        *pSpline++ = 0.0;
        *pOffset++ = 0;
      }

      // for each voxel along axis
      x = start[i]+first*step[i];
      for(j = first; j < last; j++, x += step[i])  
        {
          // locate nearest knot location greater than or equal to point
          blockIndex = (int) ceil((x - _zero(i))/_distance) - 1;
          if(blockIndex < 0)
            blockIndex = 0;
          else if(blockIndex > _n[(unsigned int)i] - 4)
            blockIndex = _n[(unsigned int)i] - 4;

          // set offset into _coef array
          currentOffset = blockIndex*offsetStep;
          *pOffset++ = currentOffset;
          for(k = 3; k > 0; k--) {
            currentOffset += offsetStep;
            *pOffset++ = currentOffset;
          }
	
          // compute terms used in tensor product
          temp = _scale*cube(_knots(blockIndex+4,i)-x);
          *pSpline++ = temp;
          *pSpline++ = 
            _scale*cube(_knots(blockIndex+5,i)-x) - 4.0*temp;
          temp = _scale*cube(x-_knots(blockIndex+3,i));
          *pSpline++ = _scale*cube(x-_knots(blockIndex+2,i)) - 
            4.0*temp;
          *pSpline++ = temp;

        }
      // fill outside domain with zeros
      for(k = last*4; k < count[i]*4; k++) {
        *pSpline++ = 0.0;
        *pOffset++ = 0;
      }
    }
}

//-----------------------------------------------------------------------
// Public Methods of TBSplineVolume


// evaluate B spline based on voxel indices
//  Warning: no range checking is done on x, y, z
double 
TBSplineVolume::operator () (int x, int y, int z) const
{
  double *xp  =  _1dSpline[0]+x*4;   // pointers into 1D spline arrays 
  double *yp  =  _1dSpline[1]+y*4;
  double *zp  =  _1dSpline[2]+z*4;
  unsigned short *xop = _offset[0]+x*4;  // pointers into offset arrays
  unsigned short *yop = _offset[1]+y*4;  // pointers into offset arrays
  unsigned short *zop = _offset[2]+z*4;  // pointers into offset arrays

  double value, xv, xyv;
  int xo, xyo, i, j, k;
  const double *coef = _coef.contents();  // direct access to array

  // compute 4 by 4 by 4 tensor and inner product with coefficients
  value = 0;
  for(i = 0; i < 4; i++) {
    xv = xp[i];
    xo = xop[i];
    for(j = 0; j < 4; j++) {
      xyv = xv*yp[j];
      xyo = xo + yop[j];
      for(k = 0; k < 4; k++)
        value += xyv*zp[k]*coef[xyo  + zop[k]];
    }
  }

  return value;
}


// incorporate new data point in fit
// returns FALSE if data is outside domain
// warning: there is no range checking on x, y, z
Boolean
TBSplineVolume::addDataPoint(int x, int y, int z, double value)
{
  int i,j,k,l;
  double *pValue = _values.contents();
  int *pLocation = _locations.contents();
  double *xp  =  _1dSpline[0]+x*4;   // pointers into 1D spline arrays 
  double *yp  =  _1dSpline[1]+y*4;
  double *zp  =  _1dSpline[2]+z*4;
  unsigned short *xop = _offset[0]+x*4;  // pointers into offset arrays
  unsigned short *yop = _offset[1]+y*4;  // pointers into offset arrays
  unsigned short *zop = _offset[2]+z*4;  // pointers into offset arrays

  double xv, xyv;
  int xo, xyo;

  _fitted = FALSE;
  _haveData = TRUE;
  _nsamples++;

  // compute 4 by 4 by 4 tensor
  for(i = 0; i < 4; i++) {
    xv = xp[i];
    xo = xop[i];
    for(j = 0; j < 4; j++) {
      xyv = xv*yp[j];
      xyo = xo + yop[j];
      for(k = 0; k < 4; k++) {
        *pValue++ = xyv*zp[k];
        *pLocation++ = xyo+zop[k];
      }
    }
  }

  // create list of differences between locations for fast indexing
  int *pDloc_i = _dloc_i;
  int *pDloc_j = _dloc_j;
  {
    pLocation = _locations.contents();
    int *pLocation2 = pLocation+1;
    for(i = _four-2; i >= 0; i--)
      {
        *pDloc_i = *pLocation2++ - *pLocation++;
        *pDloc_j++ = (*pDloc_i++)*_nProduct;
      }
  }

  // add new data to AtA and AtF
  double incr, value_k;
  double *AtAel, *AtFel;         // pointers to raw data
  double *U; //, *L;    // pointers to upper and lower triangles
  double *D;          // pointer to diagonal
  double *F;        // pointer to AtF

  AtAel = (double *) *(_AtA.getEl());
  AtFel = (double *) *(_AtF.getEl());
  // start in bottom right corner of block within AtA
  D = AtAel + _locations[(unsigned int)_four-1]*_nProduct + _locations[(unsigned int)_four-1];
  F = AtFel + _locations[(unsigned int)_four-1];
  for(k = _four-1; k >= 0; )
    {
      U = D;  // L = 
      value_k = _values[(unsigned int)k];
      *F += value*value_k;
      l = k-1;
      pValue = _values.contents() + l;
      pDloc_i = _dloc_i + l;
      //      pDloc_j = _dloc_j + l;
      for(; l >= 0; l--)  // only do half since AtA is symmetric
	{
          incr = value_k*(*pValue--);
	  U -= *pDloc_i--;
          //  L -= *pDloc_j--;
	  *U += incr;
	  // *L += incr;
	}
      *D += value_k*value_k; // do diagonal elements separately

      --k;
      if( k >= 0 ) {
        D -= _dloc_i[k] + _dloc_j[k];
        F -= _dloc_i[k];
      }
    }

  return TRUE;
}

// save internal state for analysis
void
TBSplineVolume::saveState(char *filename) const
{
#ifdef HAVE_MATLAB
  TBSpline::saveState(filename);
  DblArray xspline(_1dSpline[0], _sizes[0]*4);
  DblArray yspline(_1dSpline[1], _sizes[1]*4);
  DblArray zspline(_1dSpline[2], _sizes[2]*4);
//   SimpleArray<unsigned short> xoffset(_offset[0], _sizes[0]*4);
//   SimpleArray<unsigned short> yoffset(_offset[1], _sizes[1]*4);
//   SimpleArray<unsigned short> zoffset(_offset[2], _sizes[2]*4);
  xspline.saveMatlab(filename, "xspline");
  yspline.saveMatlab(filename, "yspline");
  zspline.saveMatlab(filename, "zspline");
//   xoffset.saveMatlab(filename, "xoffset");
//   yoffset.saveMatlab(filename, "yoffset");
//   zoffset.saveMatlab(filename, "zoffset");
  

  _knots.saveMatlab(filename, "knots");
  //  _J.saveMatlab(filename, "J");
  _AtA.saveMatlab(filename, "AtA");
  _AtF.saveMatlab(filename, "AtF");
  DblMat coef(_coef);
  coef.saveMatlab(filename, "coef");
#else
  cerr << "Warning: saveState function not available in class TBSplineVolume\n";
#endif 
}

//-----------------------------------------------------------------------
// Protected Methods of TBSplineVolume

// compute domain dimensions for TBSpline constructor
DblMat TBSplineVolume::computeDomain(const double start[VDIM],
                                     const double step[VDIM],
                                     const int count[VDIM])
{
  DblMat domain(VDIM,2);

  for(int i = 0; i < VDIM; i++)
    {
      if(step[i] > 0) { 
        domain(i,0) = start[i] - 0.5*step[i];
        domain(i,1) = start[i]+(count[i]-0.5)*step[i];
      }
      else {
        domain(i,1) = start[i] - 0.5*step[i];
        domain(i,0) = start[i]+(count[i]-0.5)*step[i];
      }        
    }
  return domain;
}



#endif

