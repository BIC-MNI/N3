#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  "version.h"

static   char *prog_name;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : print_version_info
@INPUT      : version info string
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: prints out program version information and exits
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1996 Louis Collins
@MODIFIED   : 

---------------------------------------------------------------------------- */

void  print_version_info( char *version_string )
{

    (void) printf( "The program <%s> was built from:\n", prog_name);
   
    (void) printf( "%s\n", version_string );

    exit(-1);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_program_name
@INPUT      : program name  (assumed persistent and constant)
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: set program name used in version info string
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1997 John Sled
@MODIFIED   : 

---------------------------------------------------------------------------- */

void  set_program_name( char *name )
{
  prog_name = name;
}



