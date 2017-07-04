/* ----------------------------- MNI Header -----------------------------------
@NAME       : extracttag
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 28, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : $Log: extracttag.c,v $
@MODIFIED   : Revision 1.2  2005-03-08 15:35:26  bert
@MODIFIED   : Two minor updates
@MODIFIED   :
@MODIFIED   : Revision 1.1  2003/04/16 14:27:22  bert
@MODIFIED   : Initial checkin
@MODIFIED   :
@MODIFIED   : Revision 1.5  1998/01/30 19:10:00  jgsled
@MODIFIED   : modified to reflect change to C++ compilation throughout
@MODIFIED   :
@MODIFIED   : Revision 1.4  1997/10/03 19:29:56  jgsled
@MODIFIED   : Changed version name.
@MODIFIED   :
@MODIFIED   : Revision 1.3  1997/10/01 23:40:34  jgsled
@MODIFIED   : Fixed bug related to compilation with gcc and volume_io macros.
@MODIFIED   :
@MODIFIED   : Revision 1.2  1997/09/29 19:10:55  jgsled
@MODIFIED   : Added -version option
@MODIFIED   :
@MODIFIED   : Revision 1.1  1997/08/12 21:39:11  jgsled
@MODIFIED   : Vasco's extracttag program.
@MODIFIED   : Included for release of N3 package
@MODIFIED   :
 * Revision 1.8  1996/09/09  02:21:41  vasco
 * fixed a minor bug that disallowed -random, and -maxtags 0, since random traversals are not allowed
 * when the entire qualifying set is required.
 *
 * Revision 1.7  1996/09/09  02:14:11  vasco
 * Added -random option, that traverses a volume in random (rather than fixed) stepsizes
 * thus allowing for multiple realization of the same number of tags traversed differently throught
 * the qualified tag point sets.
 *
 * Revision 1.6  1996/02/24  07:36:26  vasco
 * Few dramatic improvements like. A true -mask option, that loads in a mask file
 * rather that taking mask file limits, (previously -limit, now replaced by -mask).
 * A -world switch guarentees traversal of mask in world coordinate space.
 * The establishment of mintag threshold is done through a histogram, created from
 * a single traversal, rather than multiple.
 *
 * Revision 1.5  1995/11/04  17:29:26  vasco
 * Added -top_threshold to start minimum tag selection at a threshold other
 * than 1.0,  Also fixed a bug in stepsize for non-minimum tag requirements.
 *
 * Revision 1.4  1995/10/23  01:57:53  vasco
 * Dramatic improvements to this program.  It now taked only one
 * probability mask at a time. It has the following new switches:
 * -label to specify a class, -min to get a minimum number of
 * tagpoints while modifying thresholds the stepsize of which can
 * be controlled by -stepsize.  There is an append switch that
 * appends tag points to an existing tagfile or a tag volume.
 * There is a -limit to make sure subvolumes can be specified to
 * restrict probability volume traversals.
 *
---------------------------------------------------------------------------- */
#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */
#include <volume_io.h>
#if defined(__cplusplus)
}
#endif /* __cplusplus */
#include <ParseArgv.h>
#include <sys/time.h>
#include <time_stamp.h>
#include <version.h>

#define MAX_RAND 2147483647


/* function prototypes */

void parse_arguments(int argc, char* argv[]);
void load_volume(char *);
void load_mask_volume(char *);
void initialize_tag_volume(void);
long scan_input_volume(VIO_Real, VIO_Real);
void get_tag_points(VIO_Real, VIO_Real);
void write_tag_file(void);
void write_tag_volume(void);
void load_tag_file(char *);
void load_tag_volume(char *);
void build_histogram(void);
int mask_value_set(int vox1, int vox2, int vox3);
void find_proper_threshold(VIO_Real *lo, VIO_Real *hi);
int voxel_is_in_mask_volume( VIO_Real vox1, VIO_Real vox2, VIO_Real vox3);
int volume_size_is_ok( VIO_Volume loaded_volume);
long gen_seed(void);
float gen_random_step(void);


/* global variables */

VIO_Status     status;                
int        verbose = FALSE;
int        clobber = FALSE;
int        debug;
int        append = FALSE;                  /* true if to append tags */
VIO_Real       threshold[2] = { -1.0, -1.0 };   /* tag selection criterion */
VIO_Real       top_thresh = 1.001;              /* upper threshold for mintags */
VIO_Real       low = -1.0, high = -1.0;         /* starting thresh for mintag */ 

char       *tag_filename = NULL;           /* name of tag file */
char       *tag_volume_filename = NULL;    /* name of tag volume */
char       *prob_class;                    /* class label */
char       *prob_filename = NULL;          /* name of probability mask */
char       *mask_filename = NULL;          /* name of limiting size file */
char       *comment = NULL;                /* string to denote a comment */
char       *history = NULL;                /* string to denote a history */

VIO_Volume     prob_volume, tag_volume, mask_volume;

int        prob_volume_sizes[VIO_MAX_DIMENSIONS];  /* prob volume sizes */
int        mask_volume_sizes[VIO_MAX_DIMENSIONS];  /* mask volume sizes */
int        tag_volume_sizes[VIO_MAX_DIMENSIONS];   /* size of tag volume */

VIO_Real       **tags;              /* structure to hold tags */
char       **labels;            /* structure to hold labels*/

VIO_Real       **new_tags;              /* structure to hold tags */
char       **new_labels;            /* structure to hold labels*/

VIO_Real       wx, wy, wz;              /* world and voxel coordinates */
int        v1, v2, v3;           

int        prob_volume_num_dims; 
char       **prob_volume_dim_names;  

VIO_Real       value;

VIO_Real       output_tag_voxel = 0.0;
VIO_Real       user_mask_floor = 1.0;
VIO_Real       user_mask_ceil = -1;
VIO_Real       user_mask_binvalue = -1;
long       num_points = 0;
long       points_chosen = 0;
long       max_tags = -1;
long       min_tags = 0;

int        i = 0;      /* number of tag points */
int        new_i = 0;  /* number of new tag points */
long       histogram[101];  /* histogram of number of points / threshold */
int        world_coordinate = FALSE; /* specify world coordinate mask traversal */
int        random_step = FALSE;  /* choose fixed (rather than random) stepsize */
float      stepsize; /* the qualifying tag point stepsize */


ArgvInfo argTable[] = {

  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Show progress"},
  
  {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber,
     "Overwrite output file - if it exists."},
  
  {"-debug", ARGV_INT, (char *) NULL, (char *) &debug,
     "Show debugging information, depending on level specified"},
  
  {"-volume", ARGV_STRING, (char *) NULL, (char *) &tag_volume_filename,
     "Generate a tag volume out of the qualifying tag points"},      
  
  {"-tag", ARGV_STRING, (char *) NULL, (char *) &tag_filename,
     "Generate a tag file out of the qualifying tag points"},      
  
  {"-mintags", ARGV_LONG, (char *) NULL, (char *) &min_tags,
     "Set the minimum number of tag points to be chosen - threshold adjusted"},

  {"-maxtags", ARGV_LONG, (char *) NULL, (char *) &max_tags,
     "Set the max number of tag points to be chosen - zero for all"},
  
  {"-threshold", ARGV_FLOAT, (char *) 2, (char *) &threshold,
     "Set the threshold of selection criterion (probability level)"},

  {"-top_thres", ARGV_FLOAT, (char *) NULL, (char *) &top_thresh,
     "Set the upper threshold for minimum number of tagpoint selection"},

  {"-label", ARGV_STRING, (char *) NULL, (char *) &prob_class, 
     "Specify the voxel label of the extracted tag points"},

  {"-mask", ARGV_STRING, (char *) NULL, (char *) &mask_filename, 
     "Specify a mask volume, to limit search space."},

  {"-user_mask_value", ARGV_FLOAT, (char *) NULL, (char *) &user_mask_floor,
     "A synonym for -mask_floor"},

  {"-mask_floor", ARGV_FLOAT, (char *) NULL, (char *) &user_mask_floor,
     "Exclude mask voxels below this value"},

  {"-mask_ceil", ARGV_FLOAT, (char *) NULL, (char *) &user_mask_ceil,
     "Exclude mask voxels above this value"},

  {"-mask_binvalue", ARGV_FLOAT, (char *) NULL, (char *) &user_mask_binvalue,
     "Include mask voxels within 0.5 of this value"},

  {"-world", ARGV_CONSTANT, (char *) TRUE, (char *) &world_coordinate,
     "Specify that the mask be traversed in world coordinate space."},

  {"-random", ARGV_CONSTANT, (char *) TRUE, (char *) &random_step,
     "Specify that qualified tag points be traversed in random stepsizes."},

  {"-comment", ARGV_STRING, (char *) NULL, (char *) &comment, 
     "Specify a comment to be included in the tag file."},

  {"-append", ARGV_CONSTANT, (char *) TRUE, (char *) &append,
     "Append to an existing tagfile or tagvolume (thus clobbering)."},

  {"-version", ARGV_FUNC, (char *) print_version_info, 
   (char *)MNI_LONG_VERSION, "Print out version info and exit."},
 


  {NULL, ARGV_END, NULL, NULL, NULL}
};




int main(int argc, char *argv[])
{


  /* initialize the random seed */
  srandom( - (int) gen_seed() );

  set_program_name(argv[0]);  /* for version info */

  parse_arguments(argc, argv);

  /* load the probability map */
  load_volume(prob_filename);

  if ( mask_filename ) 
    load_mask_volume( mask_filename );

  /* if tags are being appended, load corresponding tag file or volume */
  if ( append ) {

    if ( tag_filename )
      load_tag_file( tag_filename );

    if ( tag_volume_filename )
      load_tag_volume( tag_volume_filename );
  
  }

  /* if tagvolume is specified without appending, get the size and create volume */ 
  if ( tag_volume_filename && !append ) { 
    
    initialize_tag_volume();

  }

  
  if ( min_tags > 0 ) {

    /* scan volume building a histogram in the process */
    build_histogram();

    find_proper_threshold(&low, &high);

    if ( debug > 4 ) 
      fprintf ( stdout, "Qualifying low = %f, high = %f\n", low, high);

    get_tag_points(low, high);    

  }
  else {

    num_points = scan_input_volume(threshold[0], threshold[1]);
    get_tag_points(threshold[0], threshold[1]);    

  } /* if ( mintag > 0 ) */


  if ( tag_volume_filename )
    write_tag_volume();

  if ( tag_filename )
    write_tag_file();
      
  exit(EXIT_SUCCESS);

} /* main */





/* ----------------------------- MNI Header -----------------------------------
@NAME       : parse_arguments
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 28, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : Oct. 19, 1995 (Vasco KOLLOKIAN)
---------------------------------------------------------------------------- */
void parse_arguments(int argc, char* argv[])
{

  /* form the history string - before parseargv is called */
  history = time_stamp(argc, argv);

  /* Call ParseArgv */
  if ( ParseArgv(&argc, argv, argTable, 0) || (argc < 2 )) {
    (void) fprintf(stderr, 
		   "\nUsage: %s <options> <prob_mask_filename>\n", 
		   argv[0]);
    (void) fprintf(stderr,   
		   "       %s [-help]\n\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if (tag_volume_filename && !clobber && 
      file_exists(tag_volume_filename)) {
    
    printf("File %s exists.\n", tag_volume_filename);
    printf("Use -clobber to overwrite.\n");
    exit(EXIT_FAILURE);
  }
  
  if (tag_filename && !clobber && file_exists(tag_filename)) {
    
    printf("File %s exists.\n", tag_filename);
    printf("Use -clobber to overwrite.\n");
    exit(EXIT_FAILURE);
  }
  
  /* make sure the specified tag file/volume exists, if -append is used*/
  if ( append ) {

    if ( tag_filename != NULL && !file_exists(tag_filename) ) {
      
      (void) fprintf(stderr,"File `%s' doesn't exist !\n ", tag_filename);
      exit(EXIT_FAILURE);
    }

    if ( tag_volume_filename != NULL && !file_exists(tag_volume_filename) ) {
      
      (void) fprintf(stderr,"File `%s' doesn't exist !\n ", tag_volume_filename);
      exit(EXIT_FAILURE);
    }

  }

  if ( mask_filename != NULL && !file_exists(mask_filename) ) {
    
    (void) fprintf(stderr,"File `%s' doesn't exist !\n ", mask_filename);
    exit(EXIT_FAILURE);
  }

  if ( user_mask_binvalue > -1 ){
     if ( user_mask_ceil > -1 ){
        printf("Set only one of floor/ceil or binvalue\n");
        exit(EXIT_FAILURE);
     }
     user_mask_floor = user_mask_binvalue - 0.5;
     user_mask_ceil = user_mask_binvalue + 0.5;
  }

  if ( !tag_filename  &&  !tag_volume_filename  ) {
    
    printf("Please specify one of `-volume', `-tag' options\n");
    exit(EXIT_FAILURE);
  }
     
  if ( !prob_class ) {

    printf("Please specify a class label with the -label option.\n");   
    exit(EXIT_FAILURE);
  }
     
  
  if (max_tags < 0) {
    
    print( "Please specify number of tags to choose with -maxtags  \n");
    exit(EXIT_FAILURE);
  }
  
  
  if ( min_tags > max_tags ) {
    
    print( "mintags is greater than maxtags, use -maxtag to increase default.\n");
    exit(EXIT_FAILURE);
  }
 


  if ( threshold[0] > threshold[1] ) {
    
    printf("Lower threshold limit specified is higher than upper limit\n");
    exit(EXIT_FAILURE);
  }
  
  
  if ( ( threshold[0] == -1.0) && (threshold[1] == -1.0) && ( min_tags == 0 ) ) {
    
    printf("You have to specify a -threshold, or -mintags \n");
    exit(EXIT_FAILURE);
  }

  if ( ( threshold[0] != -1.0) && (threshold[1] != -1.0) && ( min_tags != 0 ) ) {
    
    printf("-threshold, and -mintags are mutually exclusive\n");
    exit(EXIT_FAILURE);
  }
  
  
  if (( tag_filename ) && (tag_volume_filename ) ) {
    
    printf("-volume and -tag are mutually exclusive options\n");
    exit(EXIT_FAILURE);
  }


  /* set the minimum tag selection thresholds */
  if ( top_thresh == 1.001 ) {

    low  = 1.0;
    high = top_thresh;
  }
  else {

    low  = top_thresh;
    high = top_thresh;
  }

  /* if inproper top_thresh is given, complain */
  if ( low > high ) {
    
    printf("invalid -top_thresh value\n");
    exit(EXIT_FAILURE);
  }



  /* get the filename */
  prob_filename = argv[1];
  
  if ( !file_exists(prob_filename) ) {

    (void) fprintf(stderr, "filename `%s' not found. \n", prob_filename);
    exit(EXIT_FAILURE);
  }

  if ( (max_tags == 0) && random_step ) {

    (void) fprintf(stderr, "-random  and -maxtags 0, are mutually exclusive.\n");
    exit(EXIT_FAILURE);
  }




} /* parse_arguments() */



/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 28, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void load_volume(char *in_filename)
{


  if (verbose)
    printf ("Processing volume %s\n", prob_filename);


  /* load the volume */
  status = input_volume( in_filename, 3, NULL, NC_UNSPECIFIED, 
			FALSE, 0.0, 0.0,
			TRUE, &prob_volume, (minc_input_options *) NULL ) ;

  if( status != VIO_OK )
    exit(EXIT_FAILURE);

  /* get prob map volume sizes, num dims, and dim orders */
  get_volume_sizes(prob_volume, prob_volume_sizes);
  prob_volume_num_dims = get_volume_n_dimensions(prob_volume);
  prob_volume_dim_names = get_volume_dimension_names(prob_volume);


} /* load_volume(...) */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_tag_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 28, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void initialize_tag_volume(void)
{


  if (verbose)
    printf("Creating tag volume...\n");
  
  /* generate a blank volume here */

  tag_volume = copy_volume_definition(prob_volume, NC_BYTE, FALSE, 0.0, 0.0);
  set_volume_voxel_range(tag_volume, 0.0, 255.0);
  set_volume_real_range(tag_volume, 0.0, 255.0); 
  
  alloc_volume_data(tag_volume);
  
  /* initialize it to zeros */
  if (verbose)
    printf("Initializing tag volume\n");
  
  output_tag_voxel = 0.0;
  
  /* here prob_volume_sizes should be the same as tag_volume_sizes */
  for ( v1 = 0; v1 < prob_volume_sizes[0]; v1++) {
    for ( v2 = 0; v2 < prob_volume_sizes[1]; v2++) {
      for ( v3 = 0; v3 < prob_volume_sizes[2]; v3++) {

	set_volume_real_value(tag_volume, v1, v2, v3, 0, 0, output_tag_voxel);

      }
    }
  }
  
} /* initialize_tag_volume */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : scan_input_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mar. 1, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : Feb. 21, 1996 (Vasco KOLLOKIAN)
---------------------------------------------------------------------------- */
long scan_input_volume(VIO_Real low, VIO_Real high)
{

  long n_points = 0;
  VIO_Real probability;          /* prob variable */


  for ( v1 = 0; v1 < prob_volume_sizes[0]; v1++) {
    for ( v2 = 0; v2 < prob_volume_sizes[1]; v2++) {
      for ( v3 = 0; v3 < prob_volume_sizes[2]; v3++) {

	if ( mask_value_set(v1, v2, v3 ) ) {

	  probability = get_volume_real_value( prob_volume, v1, v2, v3, 0, 0);
	  
	  if ( debug > 15 ) 
	    fprintf( stdout, "prob = %f\n", probability);
	  
	  if (( probability >=  low) && ( probability <= high )) {
	    /* count number of points meeting criterion */
	    n_points++;
	  }
	  
	}
	
      } /* v3 */
    } /* v2 */
  } /* v1 */
  
  return n_points;
  
} /* scan_input_volume */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_tag_points
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mar. 1, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : Sep. 9, 1996 (Vasco KOLLOKIAN)
---------------------------------------------------------------------------- */
void get_tag_points(VIO_Real lo, VIO_Real hi)
{

  long point_count = 0;

  float goodpoint, good_random;


  int k;
 
  if ( max_tags != 0 && num_points > max_tags)

    stepsize  = ( (float) num_points / (float) max_tags); /*stepsize*/         

  else                                                                         

    stepsize = 1.0;


  goodpoint = stepsize;  /* the starting point */                              


  if ( debug > 1 ) {

    fprintf( stdout, "stepsize = %f \n", stepsize);
    fprintf( stdout, "num_points = %ld \n", num_points);

  }


  /* now append the tags from the previous tag file */
  if ( append && tag_filename ) {

    if ( verbose ) 
      fprintf( stdout, "Now appending old tags\n");

    for_less ( k, 0, i ) {

      SET_ARRAY_SIZE( new_tags, new_i, new_i+1, 1000);

      ALLOC( new_tags[new_i], 3);

      new_tags[new_i][VIO_X] = tags[k][VIO_X];
      new_tags[new_i][VIO_Y] = tags[k][VIO_Y];
      new_tags[new_i][VIO_Z] = tags[k][VIO_Z];
      
      SET_ARRAY_SIZE( new_labels, new_i, new_i+1, 1000);
      ALLOC( new_labels[new_i], strlen(labels[k])+1);

      strcpy(new_labels[new_i], labels[k]);

      if ( debug >= 15 ) 
	fprintf( stdout, "%d : %f, %f, %f, <-- %s\n", 
		k, 
		new_tags[new_i][VIO_X],
		new_tags[new_i][VIO_Y],
		new_tags[new_i][VIO_Z],
		new_labels[new_i]);

      new_i++;

    } /* for_less (k,... */
  } /* if (append) */

  if ( verbose && append && tag_filename ) 
    fprintf( stdout, "Now starting the new set\n");


  for ( v1 = 0; v1 < prob_volume_sizes[0]; v1++) {
    for ( v2 = 0; v2 < prob_volume_sizes[1]; v2++) {
      for ( v3 = 0; v3 < prob_volume_sizes[2]; v3++) {

        value = get_volume_real_value(prob_volume, v1, v2, v3, 0, 0);
	
	if ( mask_value_set(v1, v2, v3 ) ) {	
	  
	  if (( value >=  lo ) && ( value <= hi )) {
	    
	    /* do your thing here */
	    point_count++;
	    
	    if ( debug >= 10 ) 
	      fprintf( stdout, "goodpoint\n");
	    
	    
	    if ( ( point_count == (long) VIO_FLOOR( goodpoint)) && 
		( max_tags ? (points_chosen < max_tags): TRUE ) ) {
	      
	      
	      if ( debug >= 8 ) 
		fprintf( stdout, "tagpoint selected\n");
	      
	      points_chosen++;

	      if ( random_step ) {

		while ( (good_random = gen_random_step()) < 1 );
		goodpoint += good_random;
	      }

	      else
		goodpoint += stepsize;

	      if ( tag_filename != NULL) {
		convert_3D_voxel_to_world(prob_volume, (VIO_Real) v1, (VIO_Real) v2, 
					  (VIO_Real) v3, &wx, &wy, &wz); 
		
		SET_ARRAY_SIZE( new_tags, new_i, new_i+1, 1000);
		ALLOC( new_tags[new_i], 3);
		
		new_tags[new_i][VIO_X] = wx;
		new_tags[new_i][VIO_Y] = wy;
		new_tags[new_i][VIO_Z] = wz;
		
		SET_ARRAY_SIZE( new_labels, new_i, new_i+1, 1000);
		  ALLOC( new_labels[new_i], strlen(prob_class)+1);
		
		strcpy(new_labels[new_i], prob_class);
		
		if ( debug >= 15 ) {
		  fprintf( stdout, "%d : %f, %f, %f, --> %s\n", 
			  new_i, 
			  new_tags[new_i][VIO_X],
			  new_tags[new_i][VIO_Y],
			  new_tags[new_i][VIO_Z],
			  new_labels[new_i]);
		  
		}
		
		new_i++;
	      }
	      
	      
	      if ( tag_volume_filename != NULL) {
		
		/* put the voxel into the mask volume */
		output_tag_voxel = atof(prob_class);
		set_volume_real_value(tag_volume, v1, v2, v3, 0, 0, output_tag_voxel);
		
	      }
	    } /* if (point_count == (long) VIO_FLOOR ...) */
	  } /* if ( value >= threshold ...) */
	} /* if ( mask_value_set(..)) */
      } /* for v3 */
    } /* for v2 */
  } /* for v1 */
  
  if ( verbose ) {

    fprintf( stdout, "Number of qualifying tags = %ld\n", point_count);
    fprintf( stdout, "Number of tags chosen = %ld\n", points_chosen);

  }




} /* get_tag_points */
    
/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_tag_file
@INPUT      : name of tag file
@OUTPUT     : 
@RETURNS    : number of tag points read
@DESCRIPTION: opens and loads a tag file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 29, 1995 ( Vasco KOLLOKIAN)
@MODIFIED   : Oct 22, 1995 ( Vasco KOLLOKIAN)
---------------------------------------------------------------------------- */
void load_tag_file ( char *tag_filename )
{

  int num_vol, k;

  if (verbose)
    (void) fprintf(stdout, "Loading  tagfile %s\n", tag_filename);

  /* read the tag file */
  if ( input_tag_file(tag_filename, &num_vol, &i,
		      &tags, NULL, NULL, NULL, NULL, &labels ) != VIO_OK ) {

    fprintf(stderr, "Error reading the tag file.\n");
    exit(EXIT_FAILURE);
  }


  if ( debug >= 20 ) {

    fprintf( stdout, "Dumping loaded tag file...\n");

    for_less ( k, 0, i ) 
      fprintf( stdout, "%d : %f, %f, %f, <-- %s\n", 
	      k, 
	      tags[k][VIO_X],
	      tags[k][VIO_Y],
	      tags[k][VIO_Z],
	      labels[k]);
    
  }

}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : write_tag_file
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mar. 1, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void write_tag_file(void)
{

  int k;

  if ( debug >= 10 ) {

    for_less ( k, 0, new_i ) 
      fprintf( stdout, "%d : %f, %f, %f, = %s\n", 
	      k, 
	      new_tags[k][VIO_X],
	      new_tags[k][VIO_Y],
	      new_tags[k][VIO_Z],
	      new_labels[k]);
    

  }

  if ( points_chosen > 0 ) {

    if (verbose)
      printf("Writing tag file %s...\n", tag_filename);
    
    if (output_tag_file(tag_filename, comment, 1, new_i, new_tags, 
			NULL, NULL, NULL, NULL, new_labels) != VIO_OK)
      exit(EXIT_FAILURE);
  
  }
  else {

    fprintf( stdout, "No points were chosen by the specified criteria\n");
    exit(EXIT_FAILURE);
  }


}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_tag_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : ?
@MODIFIED   : Oct 22, 1995 ( Vasco KOLLOKIAN )
---------------------------------------------------------------------------- */
void load_tag_volume(char *tag_volume_filename)
{

    
  if (verbose) 
    fprintf (stdout, "Loading tag volume %s\n", tag_volume_filename);
    
  /* load the volume */
  status = input_volume(tag_volume_filename,
			3, 
			NULL,
			NC_BYTE, 
			FALSE, 
			0.0, 0.0,
			TRUE, 
			&tag_volume,
			(minc_input_options *) NULL ) ;
    

  if ( status != VIO_OK )
    exit(EXIT_FAILURE);
    
  /* get the volume sizes */
  get_volume_sizes(tag_volume, tag_volume_sizes);
    
  if ( !volume_size_is_ok( tag_volume) ) {

    fprintf( stderr, "in tag volume %s\n", tag_volume_filename);
    exit(EXIT_FAILURE);
  }
 
            
} /* load_tag_volume */
        

/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_mask_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : ?
@MODIFIED   : Oct 22, 1995 ( Vasco KOLLOKIAN )
---------------------------------------------------------------------------- */
void load_mask_volume(char *mask_filename)
{

  if (verbose) 
    fprintf (stdout, "Loading mask volume %s\n", mask_filename);
    
  /* load the volume */
  status = input_volume(mask_filename,
			3, 
			NULL,
			NC_BYTE, 
			FALSE, 
			0.0, 0.0,
			TRUE, 
			&mask_volume,
			(minc_input_options *) NULL ) ;
    

  if ( status != VIO_OK )
    exit(EXIT_FAILURE);

  get_volume_sizes(mask_volume, mask_volume_sizes);
    
  /* if world coordinates are not specified, and size is different, die */
  if ( !world_coordinate &&  !volume_size_is_ok( mask_volume) ) {
      
    fprintf( stderr, "in mask volume %s\n", mask_filename);
    exit(EXIT_FAILURE);
  }

} /* load_mask_volume */



/* ----------------------------- MNI Header -----------------------------------
@NAME       : write_tag_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mar. 1, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void write_tag_volume(void)
{
      

  if ( points_chosen > 0 ) {


    if (verbose)
      printf("Writing tag volume to %s\n", tag_volume_filename);
    
    /* write out the mask volume here */
    status = output_volume( tag_volume_filename, NC_BYTE, FALSE, 0.0, 0.0,
			   tag_volume, history,
			 (minc_output_options *) NULL ) ;
    
    if ( status != VIO_OK )
      exit(EXIT_FAILURE);
  }
  else {

    fprintf( stdout, "No points were chosen by the specified criteria\n");
    exit(EXIT_FAILURE);
  }
}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : build_histogram
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: builds a histogram of prob points at different thresholds
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 21, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void build_histogram(void)
{

  int k;                     /* counter */
  VIO_Real probability;          /* prob variable */
  VIO_Real mask_value;           /* mask variable */
  long total_test = 0;       /* check the total voxels */

  /* reset the histogram */
  for_less ( k, 0, 101 )
    histogram[k] = 0;

  for ( v1 = 0; v1 < prob_volume_sizes[0]; v1++) {
    for ( v2 = 0; v2 < prob_volume_sizes[1]; v2++) {
      for ( v3 = 0; v3 < prob_volume_sizes[2]; v3++) {
	
	/* see is the mask is set at this coordinate */
	if ( mask_value_set(v1, v2, v3 ) ) {
	  
	  probability = get_volume_real_value( prob_volume, v1, v2, v3, 0, 0);

	  if ( debug > 15 ) 
	    fprintf( stdout, "prob = %f\n", probability);
	  
	  /* increase the corresponding threshold bin */
	  histogram[ (int) (probability * 100) ]++;
	  
	} 
      } /* v3 */
    } /* v2 */
  } /* v1 */

  if ( debug > 5 ) {

    for_less ( k, 0, 101 ) {

      fprintf( stdout, "histogram[%d] = %ld\n", k, histogram[k]);;
      total_test += histogram[k];
    }
    
    fprintf( stdout, "total voxels  = %ld\n", total_test);
  }

} /* build_histogram */



/* ----------------------------- MNI Header -----------------------------------
@NAME       : mask_value_set
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: checks to see if vox1,2,3 have a mask set value
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 21, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int mask_value_set(int vox1, int vox2, int vox3)
{

  VIO_Real mask_wx, mask_wy, mask_wz;  /* world coordinates of mask */
  VIO_Real mask_v1, mask_v2, mask_v3;  /* voxel coordinates of mask */
  VIO_Real mask_value;                 /* var to hold the masking state */

  /* if no mask filename is specified, then it is always set */
  if ( !mask_filename ) 
    return TRUE;
	
  /* if a mask is specified, make sure the voxel you are considering,
     is turned on in the mask. */
 
  if ( world_coordinate ) {

  /* if world coordinate is specified, it is checked against the world
     coordinate of the mask at the very same location */
    convert_3D_voxel_to_world(prob_volume, 
			      (VIO_Real)vox1,
			      (VIO_Real)vox2,
			      (VIO_Real)vox3,
			      &mask_wx,
			      &mask_wy,
			      &mask_wz);
    

    convert_3D_world_to_voxel(mask_volume, 
			      mask_wx,
			      mask_wy,
			      mask_wz,
			      &mask_v1,
			      &mask_v2,
			      &mask_v3);
  }
  else {

    /* otherwise the voxel coordinates are considered */
    mask_v1 = vox1;
    mask_v2 = vox2;
    mask_v3 = vox3;
  }

  if ( voxel_is_in_mask_volume( mask_v1, mask_v2, mask_v3) ) {

    mask_value = get_volume_real_value(mask_volume, 
				       VIO_ROUND(mask_v1), 
				       VIO_ROUND(mask_v2), 
				       VIO_ROUND(mask_v3), 
				       0, 0);
	  
    /* if mask is on, return TRUE */
    if ( mask_value >= user_mask_floor ){
      if( user_mask_ceil > -1 && mask_value > user_mask_ceil ){
        return FALSE;
      } else{
        return TRUE;
      }
    } else{
      return FALSE;
    }
  }
  else
    return FALSE;

} /* mask_value_set */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_proper_threshold
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: given min_tags, top_thres, gives the threshold the meets criterion
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 21, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void find_proper_threshold(VIO_Real *lo, VIO_Real *hi)
{

  int k;

  int start;

  /* set the starting point as the top threshold */
  *hi = top_thresh;

 /* round the starting point to integer form, from 0 - 100 */
  for ( k = (int) top_thresh * 100; k >= 0; k-- ) {

    num_points += histogram[k];

    if ( debug > 3 ) 
      fprintf( stdout, "Scanning at %d%% theshold\n", k);

    if ( num_points >= min_tags ) {

      *lo = ((VIO_Real) k)/100.0;

      if ( verbose ) 
	fprintf( stdout, "Met mintags criterion at %d%% theshold\n", k);

      return;
    }

  }

  fprintf( stderr, "Could not find %ld tag points\n", min_tags);
  exit(EXIT_FAILURE);

} /* find_proper_threshold */
    


/* ----------------------------- MNI Header -----------------------------------
@NAME       : volume_size_is_ok
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: verifies that volume sizes are VIO_OK
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 10, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int volume_size_is_ok( VIO_Volume loaded_volume) 
{


  int    *loaded_volume_sizes;
  int    loaded_volume_num_dims;
  VIO_STR *loaded_volume_dim_names;
  
  /* allocate memory for first volume sizes */
  ALLOC(loaded_volume_sizes, VIO_MAX_DIMENSIONS); 

  /* get dim size, nums, order */
  get_volume_sizes(loaded_volume, loaded_volume_sizes);
  loaded_volume_num_dims = get_volume_n_dimensions(loaded_volume);
  loaded_volume_dim_names = get_volume_dimension_names(loaded_volume);
  
  if ( debug > 2 ) {

    int k; /* local counter */      
    
    fprintf(stdout, "Vol number of dims. = %d\n", loaded_volume_num_dims);
    
    fprintf(stdout, "Vol dimension names = ");
    for_less ( k, 0, loaded_volume_num_dims ) 
      fprintf(stdout, "%s ", loaded_volume_dim_names[k]);
    fprintf(stdout, "\n");
    
  }

  /* all the volume dimensions should be the same as the first volume */

  /* check for number of dimensions mismatch */
  if (loaded_volume_num_dims != prob_volume_num_dims ) {

    (void) fprintf(stderr,"Error - Number of dimensions mismatch ");
    return FALSE;    
  }

     
  /* check for volume size mismatches */
  if (loaded_volume_sizes[VIO_X] != prob_volume_sizes[VIO_X]) {

    (void) fprintf(stderr,"Error - VIO_Volume size mismatch in X dimension ");
    return FALSE;
  }

  if (loaded_volume_sizes[VIO_Y] != prob_volume_sizes[VIO_Y]) {

    (void) fprintf(stderr,"Error - VIO_Volume size mismatch in Y dimension ");
    return FALSE;
  }

  if (loaded_volume_sizes[VIO_Z] != prob_volume_sizes[VIO_Z]) {

    (void) fprintf(stderr,"Error - VIO_Volume size mismatch in Z dimension ");
    return FALSE;
  }


  /* check for dimensions order mismatch */
  if ( strcmp(loaded_volume_dim_names[0], prob_volume_dim_names[0]) ) {
      
    (void) fprintf(stderr,"Error - First dimension order mismatch ");
    return FALSE;
  }

  /* if there are more than 1 dimension - check dim order of second dim */
  if ( loaded_volume_num_dims > 1) 
    if ( strcmp(loaded_volume_dim_names[1], prob_volume_dim_names[1]) ) {
      
      (void) fprintf(stderr,"Error - Second dimension order mismatch ");
      return FALSE;
    }

  /* if there are more than 2 dimensions - check dim order of third dim*/
  if ( loaded_volume_num_dims > 2) 
    if ( strcmp(loaded_volume_dim_names[2], prob_volume_dim_names[2]) ) {
      
      (void) fprintf(stderr,"Error - Third dimension order mismatch ");
      return FALSE;
    }

  /* if there are more then 3 dimensions - warn and die */
  if ( loaded_volume_num_dims > 3) {

    (void) fprintf(stderr,"Support is limited to 3 spatial dimensions ");
    exit(EXIT_FAILURE);
  }

  /* free the reserved memory locations */
  delete_dimension_names( loaded_volume, loaded_volume_dim_names );
  FREE(loaded_volume_sizes);

  return TRUE;
            
} /* volume_size_is_ok */
        

/* ----------------------------- MNI Header -----------------------------------
@NAME       : voxel_is_in_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: check to see if a voxel is in the volume (vol 0 is same as all)
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep 22, 1995 ( Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int voxel_is_in_mask_volume( VIO_Real vox1, VIO_Real vox2, VIO_Real vox3)
{
 
  if ( vox1 < -0.5 || vox1 >= (VIO_Real) mask_volume_sizes[0] - 0.5) {

    return FALSE;
  }
  
  else if ( vox2 < -0.5 || vox2 >= (VIO_Real) mask_volume_sizes[1] - 0.5) {
    
    return FALSE;
  }
  
  else if ( vox3 < -0.5 || vox3 >= (VIO_Real) mask_volume_sizes[2] - 0.5) {
    
    return FALSE;
  }
  else
    return TRUE;
}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : gen_seed
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: generate time-varying seed for random
@METHOD     : 
@GLOBALS    :
@CALLS      : 
@CREATED    : February 6, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
long int gen_seed(void)
{
  struct timeval   tp;
  long   tempo, divisor, res, tempo_sec;
  int    x;

  tempo = gettimeofday(&tp, NULL);
  tempo_sec = tp.tv_sec;
   /*    printf("time of day = %d, time in sec = %ld\n",
        tempo, tempo_sec); */
  res = tempo_sec;
  divisor = 10.0;
  while (res >= 10) {
    divisor = divisor *10.0;
    res = tempo_sec / divisor;
  }
  res = tempo_sec;
  while (divisor >1000) {
    x = res/divisor;
    res = res - x * divisor;
    divisor /= 10;
  }
  return res;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : gen_random_step
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: generate a random number between 0 and stepsize
@METHOD     : 
@GLOBALS    :
@CALLS      : 
@CREATED    : February 19, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : Sep 9, 1996 (Vasco KOLLOKIAN)
---------------------------------------------------------------------------- */
float gen_random_step(void)
{

  float i, r, temp;

  i = (float) MAX_RAND /  (float) stepsize;
  temp = (float) random();
  r = temp / i;

  if ( debug >= 9 ) {

    /* fprintf ( stdout, "i = %f\n", i);
    fprintf ( stdout, "temp = %f\n", temp); */
    fprintf ( stdout, "random_step = %f\n", r);
  }
  return  r;
}
