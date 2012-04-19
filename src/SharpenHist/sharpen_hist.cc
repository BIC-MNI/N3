/* ----------------------------- MNI Header -----------------------------------
@NAME       : sharpen_hist.c
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : error status
@DESCRIPTION: Creates a mapping for minclookup to sharpen a volume's histogram
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 27, 1996   J.G.Sled
@MODIFIED   : $Log: sharpen_hist.cc,v $
@MODIFIED   : Revision 1.4  2005-03-10 22:09:47  bert
@MODIFIED   : Fix sharpen_hist.cc and finite() once more
@MODIFIED   :
@MODIFIED   : Revision 1.3  2005/03/08 15:55:11  bert
@MODIFIED   : Add checks for finite() vs. isfinite()
@MODIFIED   :
@MODIFIED   : Revision 1.2  2003/11/17 04:30:56  stever
@MODIFIED   : *** empty log message ***
@MODIFIED   :
@MODIFIED   : Revision 1.1  2003/04/16 14:30:17  bert
@MODIFIED   : Initial checkin
@MODIFIED   :
@MODIFIED   : Revision 1.7  2001/03/22 19:35:10  jgsled
@MODIFIED   : ifdef added for location of finite() on Solaris
@MODIFIED   :
@MODIFIED   : Revision 1.6  2001/03/22 16:15:40  jgsled
@MODIFIED   : modified for gcc compliance
@MODIFIED   :
@MODIFIED   : Revision 1.5  1999/08/10 18:51:13  jgsled
@MODIFIED   : removes non-portable reference to nan.h
@MODIFIED   :
 * Revision 1.4  1997/10/03  19:26:49  jgsled
 * Removed superfluous range warning.
 *
@MODIFIED   : Revision 1.3  1997/08/12 21:32:48  jgsled
@MODIFIED   : Modified for release of N3 package
@MODIFIED   :
@MODIFIED   : Revision 1.2  1997/02/06 15:57:58  jgsled
@MODIFIED   : Fixed bug re: NANs
@MODIFIED   :
 * Revision 1.1.1.1  1996/08/23  19:58:08  jgsled
 * Sources for SharpenHist
 *
 * Revision 1.3  1996/07/29  16:01:49  jgsled
 * last version before modification for inclusion in the AI source tree
 *
 * Revision 1.2  1996/02/28  17:35:58  jgsled
 * Working version.  Results verifiied using matlab implementation.
 *
 * Revision 1.1.1.1  1996/02/28  02:25:14  jgsled
 * SharpenHist sources
 *
@COPYRIGHT  :
---------------------------------------------------------------------------- */

#include <config.h>
#include <stdio.h>
#include <math.h>
#include <EBTKS/Matrix.h>	// (bert) - Added EBTKS subdirectory
//#include <nan.h>
#include "args.h"
#ifdef NEED_IEEEFP
#include <ieeefp.h>
#endif

using namespace std;


void load_histogram(char *filename, DblMat &X, 
		    double *min_bin, double *max_bin);
void  save_lookup(char *filename, DblMat &Y, double zero, double one,
		  double min_bin, double max_bin); 
DblMat gaussian(double fwhm, int size);
CompMat weiner(CompMat &blur, double noise);
void  non_negative(DblMat *X);

#if HAVE_FINITE
#ifndef finite
extern "C" int finite(double);
#endif /* finite() not defined (as macro) */
#define N3FINITE(x) finite(x)
#elif HAVE_ISFINITE
#ifndef isfinite
extern "C" int isfinite(double);
#endif /* isfinite() not defined (as macro) */
#define N3FINITE(x) isfinite(x)
#else
#error "Neither finite() nor isfinite() is defined on your system"
#endif /* HAVE_ISFINITE */

//
// Main program
//

int 
main(int argc, char *argv[])
{
  Mat<double> X, Y;    // X and Y are interpreted throughout as column vectors
  double min_bin, max_bin;  // center value for least and most bins
  args args(argc, argv);
  int i;

  load_histogram(args.inputPath.string(), X, &min_bin, &max_bin);
  if(min_bin == max_bin) {
    cerr << "Histogram " << args.inputPath.string() 
         << " is empty or degenerate" << endl;
    exit(-1);  // Failure
  }
  if(args.verbose == TRUE) cout << "Loaded histogram with " << X.getrows()
			    << " bins." << endl;

  // determine size of working matrices
  int padded_size = int(pow(2,ceil(log((double)X.getrows())/log(2.0))+1) + .5);
  int offset = (padded_size - X.getrows())/2;

  // create filters
  double slope = (max_bin - min_bin) / double(X.getrows() - 1);
  CompMat blur = fft(gaussian( args.fwhm/slope, padded_size),0,1); 
  if(args.debug_flag == TRUE)
    gaussian( args.fwhm/slope, padded_size).saveAscii("gaussian.txt");
  CompMat filter = weiner(blur, args.noise);

  if(args.debug_flag == TRUE)
    {
      DblMat kernel = real(ifft(filter,0,1));
      kernel.saveAscii("filter_ifftr.txt");
      real(filter).saveAscii("filter_r.txt");
      imag(filter).saveAscii("filter_i.txt");
    }

  // zero pad X 
  DblMat X_padded(padded_size,1, 0.0);
  X_padded.insert(X,offset,0);
 
  DblMat f;
  if(args.blur_flag == FALSE) 
    {
      // compute filtered distribution f
      f = real(ifft(pmultEquals(asCompMat(X_padded).fft(0,1), filter),0,1));
      if(args.debug_flag == TRUE) 
        f.saveAscii("distribution.txt");

      // make sure that f is non-negative
      non_negative(&f);
    }
  else
    {
      f = X_padded;
    }

  // create moment array
  DblMat moment(padded_size,1);
  for(i = 0; i < padded_size; i++)
    {
      moment(i,0) = (min_bin +  (i - offset)*slope)*f(i,0);
    }

  // compute mapping
  Mat<double> Y_padded = 
    pdiv(real(ifft(pmultEquals(asCompMat(moment).fft(0,1), blur),0,1)),
			real(ifft(pmultEquals(asCompMat(f).fft(0,1),blur),0,1)));
  Y = Y_padded(offset,offset+X.getrows()-1,0,0);

  // remove any NANs of INFs from Y 
  for(i = 0; i < Y.getrows(); i++)
    {
      if(!N3FINITE(Y(i,0)))
         Y(i,0) = 0.0;
    }
  
  // check that range has shrunk
  //  if(Y(0,0) < min_bin || Y(Y.getrows()-1,0) > max_bin)
  //  {
  //    cout << "Warning: range of mapping is not within domain." << endl; 
  //  }

  // save as a lookup table
  if(args.given_min == FALSE) args.range[0] = min_bin;
  if(args.given_max == FALSE) args.range[1] = max_bin;
  save_lookup(args.outputPath.string(), Y, args.range[0],
	      args.range[1], min_bin, max_bin);

  return 0;
}


// create gaussian function in column vector with center at (0,0)
// domain of mapping is 0 to size-1
// normalized to unit area
DblMat 
gaussian(double fwhm, int size)
{
  DblMat X(size,1,0.0);
  
  double factor = 4.0*log(2.0)/(fwhm*fwhm);
  double scale = 2.0*sqrt(log(2.0)/M_PI) / fwhm;
  X(0,0) = scale;
  
  for(int i = 1; i <= (size-1)/2; i++)
    {
      X(i,0) = X(size-i,0) = scale*exp(-i*i*factor);
    }
  
  if((size-1)/2 != size/2) // if size is even
    // fill in middle value
    X(size/2,0) = scale*exp(-size*size*factor/4.0);

  return(X);
}

// create Weiner restoration filter in frequency domain assuming
//  white noise model
CompMat 
weiner(CompMat &blur, double noise)
{
  CompMat filter(blur.getrows(),1);
  int size = blur.getrows();

  for(int i = 0; i < size; i++)
    {
      filter(i,0) = conj(blur(i,0))/(conj(blur(i,0))*blur(i,0) + noise);
    }

  return(filter);
}

// make negative elements zero
void 
non_negative(DblMat *X)
{
  int rows = X->getrows();
  int cols = X->getcols();
  for(int i = 0; i < rows; i++)
    for(int j = 0; j < cols; j++)
      if((*X)(i,j) < 0)
	(*X)(i,j) = 0;
}


// read in text file in two columns of bins and counts
void
load_histogram(char *filename, DblMat &X, double *min_bin, double *max_bin)
{
  FILE *fp;
  const int SIZE = 300;
  char str[SIZE];
  SimpleArray<double> counts;
  double bin, count;

  if((fp = fopen(filename,"r")) == NULL)
    {
      cerr << "Unable to open file: " << filename << endl;
      exit(EXIT_FAILURE);
    }

  while(!feof(fp)) 
    {
      char *ptr;

      *str = '\0';
      fgets(str, SIZE, fp);
      
      // remove everything following a '#' character from the string
      for(ptr = str; *ptr; ptr++)
	if(*ptr == '#')
	  {
	    *ptr = 0;
	    break;
	  }

      // read values
      if(sscanf(str, "%lf %lf", &bin, &count) == 2)
	{
	  if(counts.size() == 0)
	    *min_bin = bin;

	  counts.append(count);
	}
    }

  // set bin range
  if(counts.size() > 0)
    {
      *max_bin = bin;
    }
  else
    {
      *max_bin = *min_bin = 0;
    }

  // copy counts array to matrix X
  if(counts.size() > 0)  X.resize(counts.size(),1);
  for(int i = 0; i < counts.size(); i++)
    X(i,0) = counts[i];

  fclose(fp);
} 



// write out text file in two columns suitable for minclookup
void 
save_lookup(char *filename, DblMat &Y, 
	    double zero, double one,  // values mapped to zero and one 
	    double min_bin, double max_bin) // domain of mapping
{
  FILE *fp;

  if((fp = fopen(filename,"w")) == NULL)
    {
      cerr << "Unable to open file: " << filename << endl;
      exit(EXIT_FAILURE);
    }

  double new_min = (min_bin - zero)/(one - zero); // this is not a joke
  double new_max = (max_bin - zero)/(one - zero); 
  double slope = (new_max - new_min) / double(Y.getrows() - 1);

  if(zero < min_bin)
    fprintf(fp, "0.0     %lf\n", Y(0,0));

  for(int i = 0; i < Y.getrows(); i++)
    {
      fprintf(fp, "%lf     %lf\n",  new_min +  i*slope, Y(i,0));
    }

  if(one > max_bin)
    fprintf(fp, "1.0     %lf\n", Y(Y.getrows()-1,0));

  fclose(fp);
}









