#!/usr/bin/env perl


############################# MNI Header #####################################
#@NAME       :  compare_nu_result.pl
#@DESCRIPTION:  check if the N3 output is withing a range
#@COPYRIGHT  :
#              Vladimir S. Fonov  2009
#              Montreal Neurological Institute, McGill University.
#              Permission to use, copy, modify, and distribute this
#              software and its documentation for any purpose and without
#              fee is hereby granted, provided that the above copyright
#              notice appear in all copies.  The author and McGill University
#              make no representations about the suitability of this
#              software for any purpose.  It is provided "as is" without
#              express or implied warranty.
###############################################################################


use strict;
use Getopt::Long;
use File::Basename;
use File::Temp qw/ tempdir /;

my $verbose=0;
my $clobber=1;
my $fake=0;
my $mask;
my $me = basename ($0);

GetOptions(
      'verbose'           => \$verbose,
      'clobber'           => \$clobber,
      'mask=s'            => \$mask,
     );

my $Help = <<HELP;
  Usage: $me <input.mnc> <reference.mnc> [tolerance]
    --verbose be verbose
    --mask <mask.mnc>
  Problems or comments should be sent to: vladimir.fonov\@gmail.com
HELP

die $Help if $#ARGV < 2;

my ($in,$ref,$tol)=@ARGV;

$tol=1e-5 if !$tol;

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );
#my $tmpdir="/tmp";

my @args=("nu_estimate", $in, "$tmpdir/brain.imp");
push(@args,'-mask',$mask) if $mask;
push(@args,'-verbose') if $verbose;

do_cmd(@args);

do_cmd("nu_evaluate", $in, '-mapping', "$tmpdir/brain.imp", "$tmpdir/brain_nu.mnc",'-verbose');

do_cmd('mincmath','-sub','-float',"$tmpdir/brain_nu.mnc" , $ref, "$tmpdir/diff.mnc");

my $mean=`mincstats -q -mean $ref`;
my $count=`mincstats -q -count $tmpdir/diff.mnc`;
my $sum2=`mincstats -q -sum2 $tmpdir/diff.mnc`;

chomp($mean);chomp($count);chomp($sum2);

my $rms=sqrt($sum2/$count);

my $rms_pc=$rms/$mean;

if($rms_pc>$tol) {
  die "relative RMS difference: $rms_pc larger then $tol\n";
} 

sub do_cmd { 
    print STDOUT "@_\n" if $verbose;
    if(!$fake){
      system(@_) == 0 or die "DIED: @_\n";
    }
}
