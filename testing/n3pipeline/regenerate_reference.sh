#!/bin/sh
# Record the oracle answers the N3Pipeline tests compare against.
#
# This is the only thing here that runs an N3 program.  The tests read what it
# writes; re-running it and finding an empty `git diff` is itself a check.
#
#   sh regenerate_reference.sh
#
# Needs the installed N3 and MINC tools on PATH.  Where an option exists only
# in this tree (-gaussian_window, -parzen_sigma) the locally built and
# installed copy is used instead, and N3_LOCAL_BIN must point at it.

set -e

here=`cd \`dirname $0\` && pwd`
data=`cd $here/.. && pwd`
out=$here/reference
work=`mktemp -d`
trap "rm -rf $work" 0

N3_LOCAL_BIN=${N3_LOCAL_BIN:-/app/legacy/_install/bin}

mkdir -p $out
chunk=$data/chunk.mnc.gz
mask=$data/chunk_mask.mnc.gz

# Geometry of a volume, one number per line, in the order the tests read it:
# nx ny nz (file order), then start, step, dircos per file dimension.
geometry()
{
  mincinfo -vardims image $1 | tr ' ' '\n' | grep . > $work/dims
  for d in `cat $work/dims`; do mincinfo -dimlength $d $1; done
  for d in `cat $work/dims`; do mincinfo -attval $d:start $1; done
  for d in `cat $work/dims`; do mincinfo -attval $d:step $1; done
}

# ---------------------------------------------------------------- cycle 1
# ShrinkVolume, taken from the driver itself rather than reimplemented here:
# -tmpdir with -keeptmp leaves the shrunk volume behind.
#
# Two oracles per factor.  The Perl's own output carries the 12-bit
# quantisation of the file mincresample wrote, which is the very divergence
# the in-memory pipeline exists to remove; running the same resample to
# -double gives the sampling with that quantisation absent.  The tests hold
# the geometry and the double sampling exactly, and report the quantised one.

#
# The mask goes through the same driver call.  It is resampled onto the shrunk
# grid by resample_labels, which is trilinear thresholded at 0.5 and not
# nearest neighbour (resample_labels.in:176-177), so it is recorded here as
# cycle 2's oracle.  At an integer factor the two grids coincide and trilinear
# degenerates to selection; factor 2.5 is included so that the interpolation is
# actually exercised.

for factor in 2 3 4 2.5; do
  rm -rf $work/shrink && mkdir -p $work/shrink
  nu_estimate_np_and_em -tmpdir $work/shrink -keeptmp -shrink $factor \
      -iterations 1 -stop 0.001 -sharpen 0.15 0.01 -parzen -log \
      -distance 200 -mask $mask $chunk $work/shrink.imp -clobber > /dev/null 2>&1
  perl_out=$work/shrink/`basename $chunk .gz`
  perl_mask=$work/shrink/`basename $mask .gz`

  geometry $perl_out > $out/shrink${factor}.geom
  mincextract -double $perl_out > $out/shrink${factor}_quantised.f64
  geometry $perl_mask > $out/mask${factor}.geom
  mincextract -double $perl_mask > $out/mask${factor}.f64

  # The same resample without the file's quantisation.  -step and -nelements
  # are read back from the Perl's own output, so the rule under test is not
  # duplicated here.
  nel=`mincinfo -dimlength xspace $perl_out``echo -n ' '``mincinfo -dimlength yspace $perl_out``echo -n ' '``mincinfo -dimlength zspace $perl_out`
  stp=`mincinfo -attval xspace:step $perl_out``echo -n ' '``mincinfo -attval yspace:step $perl_out``echo -n ' '``mincinfo -attval zspace:step $perl_out`
  mincresample -clobber -quiet -nearest_neighbour -double \
      -nelements $nel -step $stp $chunk $work/shrink_double.mnc
  mincextract -double $work/shrink_double.mnc > $out/shrink${factor}_double.f64
done

geometry $chunk > $out/chunk.geom

# ---------------------------------------------------------------- cycle 3
# Masked statistics.  volume_stats prints through cout at its default six
# significant digits, which is the bound the test holds these to.  Its mask
# rule is value != 0 (volumeStats.cc:236) and its variance is the population
# one (:276), which is also what the stopping rule needs.

for stat in mean stddev min max; do
  volume_stats -quiet -$stat -mask $mask $chunk > $out/stats_masked_$stat.txt
  volume_stats -quiet -$stat $chunk > $out/stats_whole_$stat.txt
done

geometry $chunk > $out/chunk.geom

echo "reference written to $out"
ls -l $out
