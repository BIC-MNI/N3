/**
 * Smooth a N3 field outside a mask.
 * Author: Claude Lepage, 2009.
 **/

#include <iostream>
#include <iomanip>
extern "C" { 
#include <volume_io.h>
#include <time_stamp.h> 
}

using namespace std;

static short debug = 0;
#ifdef DEBUG
  debug = 1;
#endif

//
// Smooth outside of volume in the far field by solving a 
// Laplacian with non-reflective bc.
//

void smooth( int sizes[MAX_DIMENSIONS], Real seps[MAX_DIMENSIONS],
             Volume volume, Volume mask ) {

    int    i, j, k, ii;
    float  fval;

    int count;

    float SOR = 1.9;

    // Copy the volume data in temporary arrays, for speed.
    // This is faster than calling get_volume_real_value
    // many times.

    float * val = new float[sizes[0]*sizes[1]*sizes[2]];
    char * mm = new char[sizes[0]*sizes[1]*sizes[2]];

    count = 0;
    for( i = 0;  i < sizes[0];  ++i ) {
      for( j = 0;  j < sizes[1];  ++j ) {
        for( k = 0;  k < sizes[2];  ++k ) {
          if( get_volume_real_value( mask, i, j, k, 0, 0 ) > 0.5 ) {
            val[count] = get_volume_real_value( volume, i, j, k, 0, 0 );
            mm[count] = 1;
          } else {
            mm[count] = 0;
            val[count] = 0;
          }
          count++;
        }
      }
    }

    float fx = 1.0 / ( seps[0] *seps[0] );
    float fy = 1.0 / ( seps[1] *seps[1] );
    float fz = 1.0 / ( seps[2] *seps[2] );

    /* Multi-level resolution of Laplacian operator by Gauss-Seidel approach */

    /* We want the coarsest grid to be at least 4mm with 4*200 iterations */
    float min_sep = MIN( MIN( ABS(seps[0]), ABS(seps[1]) ), ABS(seps[2]) );
    int inc = 2;
    while( inc < (int)( 4.0 / min_sep ) ) inc *= 2;

    float thresh = 1.0e-10;
    int n_iters = 4 * 200;

    for( inc = inc; inc >=2; inc /= 2 ) {
      for( int iter = 0; iter < n_iters; iter++ ) {
        count = 0;
        float res = 0.0;
        for( i = 0;  i < sizes[0]; i += inc ) {
          for( j = 0;  j < sizes[1]; j += inc ) {
            for( k = 0;  k < sizes[2]; k += inc ) {
              if( mm[(i*sizes[1]+j)*sizes[2]+k] == 0 ) {
                float norm = 0.0;
                float uxx = 0.0, uyy = 0.0, uzz = 0.0;
                if( i > (inc-1) ) {
                  uxx += fx * val[((i-inc)*sizes[1]+j)*sizes[2]+k];
                  norm += fx;
                }
                if( i < sizes[0]-inc ) {
                  uxx += fx * val[((i+inc)*sizes[1]+j)*sizes[2]+k];
                  norm += fx;
                }
                if( j > (inc-1) ) {
                  uyy += fy * val[(i*sizes[1]+j-inc)*sizes[2]+k];
                  norm += fy;
                }
                if( j < sizes[1]-inc ) {
                  uyy += fy * val[(i*sizes[1]+j+inc)*sizes[2]+k];
                  norm += fy;
                }
                if( k > (inc+1) ) {
                  uzz += fz * val[(i*sizes[1]+j)*sizes[2]+k-inc];
                  norm += fz;
                }
                if( k < sizes[2]-inc ) {
                  uzz += fz * val[(i*sizes[1]+j)*sizes[2]+k+inc];
                  norm += fz;
                }

                float oldValue = val[(i*sizes[1]+j)*sizes[2]+k];
                float tmpValue = ( uxx + uyy + uzz ) / norm;
                float newValue = oldValue + SOR * ( tmpValue - oldValue );
                val[(i*sizes[1]+j)*sizes[2]+k] = newValue;
                count++;
                res += ABS( oldValue - newValue );
              }
            }
          }
        }
        res /= (float)count;
        if( debug ) cout << "Iter = " << iter << " res = " << res << endl;
        if( res < thresh ) break;
      }
      n_iters = n_iters / 2 + 1;
      // thresh *= 2.0;

      // Extension operator for next level.
      if( inc == 1 ) break;
      for( i = 0; i < sizes[0]; i += inc ) {
        for( j = 0; j < sizes[1]; j += inc ) {
          for( k = inc/2; k < sizes[2]-inc/2; k += inc ) {
            if( mm[(i*sizes[1]+j)*sizes[2]+k] == 0 ) {
              val[(i*sizes[1]+j)*sizes[2]+k] = 
                0.5 * ( val[(i*sizes[1]+j)*sizes[2]+k-inc/2] +
                        val[(i*sizes[1]+j)*sizes[2]+k+inc/2] );
            }
          }
          if( k < sizes[2] ) {
            if( mm[(i*sizes[1]+j)*sizes[2]+k] == 0 ) {
              val[(i*sizes[1]+j)*sizes[2]+k] = val[(i*sizes[1]+j)*sizes[2]+k-inc/2];
            }
          }
        }

        for( j = inc/2; j < sizes[1]-inc/2; j += inc ) {
          for( k = 0; k < sizes[2]; k += inc/2 ) {
            if( mm[(i*sizes[1]+j)*sizes[2]+k] == 0 ) {
              val[(i*sizes[1]+j)*sizes[2]+k] = 
                0.5 * ( val[(i*sizes[1]+j-inc/2)*sizes[2]+k] +
                        val[(i*sizes[1]+j+inc/2)*sizes[2]+k] );
            }
          }
        }
        if( j < sizes[1] ) {
          for( k = 0; k < sizes[2]; k += inc/2 ) {
            if( mm[(i*sizes[1]+j)*sizes[2]+k] == 0 ) {
              val[(i*sizes[1]+j)*sizes[2]+k] = val[(i*sizes[1]+j-inc/2)*sizes[2]+k];
            }
          }
        }
      }

      for( i = inc/2; i < sizes[0]-inc/2; i += inc ) {
        for( j = 0; j < sizes[1]; j += inc/2 ) {
          for( k = 0; k < sizes[2]; k += inc/2 ) {
            if( mm[(i*sizes[1]+j)*sizes[2]+k] == 0 ) {
              val[(i*sizes[1]+j)*sizes[2]+k] = 
                0.5 * ( val[((i+inc/2)*sizes[1]+j)*sizes[2]+k] +
                        val[((i-inc/2)*sizes[1]+j)*sizes[2]+k] );
            }
          }
        }
      }
      if( i < sizes[0] ) {
        for( j = 0; j < sizes[1]; j += inc/2 ) {
          for( k = 0; k < sizes[2]; k += inc/2 ) {
            if( mm[(i*sizes[1]+j)*sizes[2]+k] == 0 ) {
              val[(i*sizes[1]+j)*sizes[2]+k] = val[((i-inc/2)*sizes[1]+j)*sizes[2]+k];
            }
          }
        }
      }
    }

    // Obtain global min/max for field.

    float min_val = val[0];
    float max_val = val[0];

    for( i = 1;  i < sizes[0]*sizes[1]*sizes[2]; i++ ) {
      if( val[i] < min_val ) min_val = val[i];
      if( val[i] > max_val ) max_val = val[i];
    }
    set_volume_real_range( volume, (Real)min_val, (Real)max_val );

    // Save back to volume.

    count = 0;
    for( int i = 0; i < sizes[0]; i++ ) {
      for( int j = 0; j < sizes[1]; j++ ) {
        for( int k = 0; k < sizes[2]; k++ ) {
          set_volume_real_value( volume, i, j, k, 0, 0, (Real)val[count] );
          count++;
        }
      }
    }

    delete[] val;
    delete[] mm;
}


int  main( int ac, char* av[] ) {
  
    if( ac < 4 ) {
      cerr << "Usage: " << av[0] << " input.mnc mask.mnc output.mnc " 
           << endl;
      cerr << "       input.mnc = input intensity image" << endl;
      cerr << "       mask.mnc = mask of image" << endl;
      cerr << "       output.mnc = output intensity image" << endl
           << endl;
      return 1;
    }

    // Read the volume. 
    Volume in_volume;
    if ( input_volume( av[1], 3, NULL, 
                       MI_ORIGINAL_TYPE, 0, 0, 0,
                       TRUE, &in_volume, NULL ) != OK ) {
      cerr << "Error: cannot read volume file " << av[1] << endl;
      return 1;
    }

    if ( get_volume_n_dimensions( in_volume ) != 3 ) {
      cerr << "Error: volume in " << av[1] 
           << " does not have three dimensions." << endl;
      return 1;
    }

    // Read the mask. 
    Volume mask;
    if ( input_volume( av[2], 3, NULL, 
                       MI_ORIGINAL_TYPE, 0, 0, 0,
                       TRUE, &mask, NULL ) != OK ) {
      cerr << "Error: cannot read mask file " << av[2] << endl;
      return 1;
    }

    if ( get_volume_n_dimensions( mask ) != 3 ) {
      cerr << "Error: volume in " << av[2] 
           << " does not have three dimensions." << endl;
      return 1;
    }

    int sizes[MAX_DIMENSIONS];
    Real seps[MAX_DIMENSIONS];
    get_volume_sizes( in_volume, sizes );
    get_volume_separations( in_volume, seps );

    smooth( sizes, seps, in_volume, mask );

    char * history = time_stamp( ac, av );

    int rv = output_modified_volume( av[3], MI_ORIGINAL_TYPE,
                                     0, 0, 0, in_volume, av[1],
                                     history, NULL );

    if( history ) free( history );
    delete_volume( in_volume );
    delete_volume( mask );

    return ( rv != OK );

}

