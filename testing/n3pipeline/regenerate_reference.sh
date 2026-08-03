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

# Volume-sized oracles are recorded on every fourth voxel in file order.  A
# relative RMS over a systematic quarter of a volume says what one over all of
# it says, and the whole thing would put ten megabytes of float64 in a source
# tree.  The tests apply the same stride.
STRIDE=4

dump_strided()
{
  mincextract -double $1 | perl -e '
    binmode STDIN; binmode STDOUT;
    my $stride = shift;
    my $i = 0;
    while(read(STDIN, my $buf, 8) == 8) { print $buf if $i % $stride == 0; $i++; }
  ' $STRIDE > $2
}

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

# ---------------------------------------------------------------- cycle 4
# The masked histogram, in all three estimators.  -window is N3's linear split
# between the two nearest bin centres; -gaussian_window is the modification
# this tree carries, so its oracle is the locally built volume_hist and not
# the installed one, which has no such option.
#
# The text carries six decimals of everything (%lf, minchist.cc), which is
# what the comparisons are held to.

hist="-bins 200 -auto_range -mask $mask -clobber -text -select 1 -quiet"
volume_hist $hist $chunk $out/hist_plain.txt
volume_hist $hist -window $chunk $out/hist_window.txt
$N3_LOCAL_BIN/volume_hist $hist -gaussian_window 2 $chunk $out/hist_gauss2.txt

# ---------------------------------------------------------------- cycle 5
# The sharpening, run on the histogram recorded above.  Taking sharpen_hist's
# input from the same text the test reads keeps the comparison to the
# deconvolution itself: the six decimals that text carries are a separate
# divergence, measured where the pipeline is assembled.

range=`grep domain: $out/hist_window.txt | sed 's/.*domain: *//'`
sharpen_hist -clobber -quiet -fwhm 0.15 -noise 0.01 -range $range \
    $out/hist_window.txt $out/sharp_window.txt
sharpen_hist -clobber -quiet -blur -fwhm 0.15 -noise 0.01 -range $range \
    $out/hist_window.txt $out/sharp_window_blur.txt

# ---------------------------------------------------------------- cycle 6
# minclookup applied to the whole volume, written as double so that the
# comparison is against the interpolation and not against a storage type.
# Both sides read the same lookup table text, for the same reason as cycle 5.

minclookup -clobber -quiet -double -continuous -range $range \
    -lookup_table $out/sharp_window.txt $chunk $work/looked_up.mnc
dump_strided $work/looked_up.mnc $out/lookup_applied.f64

# ---------------------------------------------------------------- cycle 7
# Two different bimodal thresholds, both of which the pipeline needs.
# volume_stats builds an EBTKS histogram over ceil(voxelMax-voxelMin+1) bins of
# the file's *voxel* range (volumeStats.cc:286) and takes its biModalThreshold;
# mincstats runs 2000-bin Otsu and reports the winning bin centre.  They do not
# agree, and nothing in the drivers suggests they should.

volume_stats -quiet -biModalT $chunk > $out/bimodal_volume_stats.txt
volume_stats -quiet -biModalT -mask $mask $chunk > $out/bimodal_volume_stats_masked.txt
mincstats -quiet -biModalT $chunk > $out/bimodal_mincstats.txt
mincinfo -attval image:valid_range $chunk > $out/chunk_valid_range.txt

# ---------------------------------------------------------------- cycle 8
# The spline fit.  spline_smooth reads every volume as float
# (loadFloatVolume, splineSmooth.cc:101) and writes the result in the input's
# storage type, so the input here is a float copy of chunk.mnc: both sides then
# fit identical numbers and the only remaining difference is that the oracle's
# output is rounded to float on the way out.
#
# -full_support goes with -b_spline and not with -tp_spline, which is how the
# driver calls it and which gives the two different domains.

mincreshape -clobber -quiet -float $chunk $out/chunk_float.mnc
dump_strided $out/chunk_float.mnc $out/chunk_float.f64

for distance in 200 100 50; do
  spline_smooth -clobber -quiet -full_support -b_spline -lambda 1e-7 \
      -distance $distance -subsample 1 -mask $mask \
      $out/chunk_float.mnc $work/fit.mnc
  dump_strided $work/fit.mnc $out/fit_b${distance}.f64
done

spline_smooth -clobber -quiet -tp_spline -lambda 1e-7 -distance 200 \
    -subsample 1 -mask $mask $out/chunk_float.mnc $work/fit.mnc
dump_strided $work/fit.mnc $out/fit_tp200.f64

# The same fit evaluated on the full grid through the .imp file, which is the
# path nu_evaluate takes.
spline_smooth -clobber -quiet -full_support -b_spline -lambda 1e-7 \
    -distance 200 -subsample 1 -mask $mask -novolume \
    $out/chunk_float.mnc -compact $work/fit.imp
evaluate_field -clobber -quiet -like $out/chunk_float.mnc $work/fit.imp $work/field.mnc
dump_strided $work/field.mnc $out/field_b200_full.f64

# ---------------------------------------------------------------- cycle 9
# correct_field, which extends the field beyond the mask by relaxing Laplace's
# equation on the outside.  Its input is the masked b_spline fit at 200 mm,
# kept here as a float volume so that the transcription can be given exactly
# the numbers the oracle was given: the solve amplifies, so recomputing the
# input would not do.

spline_smooth -clobber -quiet -full_support -b_spline -lambda 1e-7 \
    -distance 200 -subsample 1 -mask $mask \
    $out/chunk_float.mnc $out/field_masked.mnc
correct_field $out/field_masked.mnc $mask $work/extended.mnc
dump_strided $work/extended.mnc $out/field_extended.f64

# ---------------------------------------------------------------- cycle 11
# One iteration of the estimation loop, with the driver's own intermediates
# kept by -save_fields: chunk_est0.mnc is exp(log - sharpened) before
# smoothing, chunk_field0.mnc is exp(residue) after it.  A per-iteration oracle
# rather than an end-to-end one.
#
# -shrink 1 so that the estimation runs on the full grid and the fewest MINC
# round trips stand between the two implementations.  -stop 0.0 never fires
# (0 alone is not a float to the driver's own regex), so the iteration count is
# fixed at one on both sides.

rm -rf $work/est && mkdir -p $work/est
nu_estimate_np_and_em -shrink 1 -iterations 1 -stop 0.0 -distance 200 \
    -b_spline 1.0e-7 -spline_subsample 1 -sharpen 0.15 0.01 -parzen -log \
    -save_fields -mask $mask $chunk $work/est/chunk.imp -clobber \
    > $work/est/log 2>&1
grep 'CV of field change' $work/est/log | sed 's/.*: *//' > $out/estimate_change.txt
dump_strided $work/est/chunk_est0.mnc $out/estimate_est0.f64
dump_strided $work/est/chunk_field0.mnc $out/estimate_field0.f64

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
