/*
 * Program to write out the macros defined in version.h as a 'sed' script.
 */

#include "version.h"

int main (void)
{
   printf ("-e s/@VERSION@/%s/g\n", MNI_VERSION);
   printf ("-e s/@LONG_VERSION@/%s/g\n", MNI_LONG_VERSION);
   exit (0);
}
