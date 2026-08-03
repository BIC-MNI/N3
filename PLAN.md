# `nu_correct` in C++, entirely in memory

## Context

N3's pipeline is a Perl layer over separate executables that communicate only through
files. `nu_correct` (`nu_estimate.in`) calls `nu_estimate_np_and_em.in`, which calls
`sharpen_volume.in`, `spline_smooth`, `volume_stats` and a dozen `mincmath` invocations,
then `nu_evaluate.in` calls `evaluate_field`, `correct_field` and `mincmath` again. Every
intermediate is written as a MINC volume and read back, quantised to the file's storage
type in transit, and the histogram and lookup table pass as text written with `%lf`.

This task replaces that layer with a single C++ program holding every intermediate in
`double` buffers in memory. The numerical blocks are not rewritten: the original
translation units are linked and called, which is the same approach
`torch_n3/_legacy/n3_shim.h` already takes for the PyTorch port's oracle. What is written
new is the driver — the logic that lives in Perl today — plus the handful of MINC
utilities the drivers lean on.

Scope decisions taken with the user:

- **A new binary alongside** the Perl drivers, which stay installed and unchanged. Both
  must be runnable on the same input; that comparison is the only way to validate this.
- **The N3 path plus thin-plate splines**: np method, `b_spline` and `tp_spline`,
  `-sharpen`, `-parzen`/`-parzen_sigma`, `-shrink`, `-mask`, `-distance`, `-lambda`,
  `-iterations`, `-stop`, `-normalize_field`, `-auto_mask`/`-bimodalT`, `-floor`,
  `-subsample`, `-mapping_dir`. Not the EM branch (it `die`s without `-sharpen`), not the
  WM branch (needs `class_statistics`, `estimate`, `lgmask`, none shipped), not `fir`
  smoothing, `-real`, `-differential`, `-initial` or `-islands`.
- **A double-precision copy of `correct_field`'s `smooth()`**, rather than linking the
  original, which does its solve on `float` arrays.
- **An estimate-only mode**: the binary writes only the `.imp` when invoked under a name
  matching `nu_estimate`, as the Perl driver decides it (`nu_estimate.in:420`).
- Install to `/app/legacy/_install`, as before. `/opt/minc/1.9.18.13` is the oracle and is
  not written to. **Nothing under `/app/torch_n3` is touched** — the vendored copies in
  `torch_n3/_legacy/n3/` must stay byte-identical to `legacy/N3/src`, which is automatic
  here because the reused sources are only linked, never edited.

## 1. What is reused, and what is new

Linked unmodified from `legacy/N3/src`:

| Source | Provides |
|---|---|
| `Splines/{Spline,TBSpline}.cc` | `TBSplineVolume`: the regularized fit, `getCoefficients`/`setCoefficients` (`Spline.h:131-134`), and evaluation through the per-grid lookup tables |
| `SplineSmooth/fieldIO.cc` | `createThinPlateSpline`, `outputCompactField` (the `.imp` writer) |
| `VolumeHist/{DHistogram,WHistogram,GHistogram}.cc` | the three histogram estimators, including the Gaussian window added in the previous task |
| `SharpenHist/sharpen_hist.cc` + `args.cc` | `gaussian()`, `weiner()`, `non_negative()` — all three are non-static free functions (`sharpen_hist.cc:71-74`) |
| EBTKS `Histogram` | `biModalThreshold()`, which is what `volume_stats -biModalT` reports |

`sharpen_hist.cc` carries a `main()`. Compile it into an `OBJECT` library with
`main` renamed, exactly as `torch_n3/_legacy/build_legacy.py:211` does:

```cmake
ADD_LIBRARY(sharpenhist_core OBJECT src/SharpenHist/sharpen_hist.cc src/SharpenHist/args.cc)
TARGET_COMPILE_DEFINITIONS(sharpenhist_core PRIVATE main=n3_sharpen_hist_unused_main)
```

`args.cc` is linked only to satisfy the references that unused `main()` makes. The
existing `sharpen_hist` executable keeps its own compilation of the file and is unaffected.

Written new, under `legacy/N3/src/N3Pipeline/`:

| File | Contents |
|---|---|
| `VolumeD.{h,cc}` | `struct VolumeD { std::vector<double> data; int count[3]; double start[3], step[3], dircos[3][3]; nc_type store_type; }`, `load()` via `input_volume(..., NC_DOUBLE, ...)`, `saveLike()` via `output_modified_volume` (the equivalent of `mincmath -copy_header`), `resampleNearest`, `resampleLabel`, `shrink`, and the masked `mean`/`stddev`/`min`/`max` that `volume_stats` prints |
| `MincTools.{h,cc}` | `applyLookup()` (`minclookup -continuous`) and `bimodalThresholdMincstats()` (2000-bin Otsu returning the winning bin centre) |
| `Sharpen.{h,cc}` | `autoRange()` and `sharpenLookup()`: `sharpen_hist.cc:100-190`'s sequence over buffers |
| `SmoothField.{h,cc}` | the double-precision transcription of `correctField.cc`'s `smooth()` |
| `FitField.{h,cc}` | build the spline over a masked buffer, fit, evaluate on any grid |
| `NuEstimate.{h,cc}` | the iteration of `nu_estimate_np_and_em.in:101-196` |
| `NuEvaluate.{h,cc}` | `nu_evaluate.in:49-80` |
| `nu_correct_cxx.cc` | `ParseArgv` argument table (the `src/VolumeHist/args.cc` pattern) and `main()` |

`torch_n3/blocks/` and `torch_n3/minc_tools.py` are a validated specification for the
three utilities with no C++ original — `apply_lut` at `minc_tools.py:12`,
`bimodal_threshold` at `:48`, `resample_like`/`shrink` at `volume.py:54-96`. Read them;
do not modify them.

## 2. Stage by stage

| Perl | C++ |
|---|---|
| `ShrinkVolume` (`:955`) | keep `start`, multiply `step` by the factor and take `ceil((n-1)/factor)+1` samples, on axes finer than `factor·min|step|` only; sample nearest neighbour |
| `CheckSampling(..., isLabel=1)` → `resample_labels` | **trilinear, then threshold at 0.5** (`resample_labels.in:176-177`), not nearest neighbour |
| `log_transform` (`:657`) | `log(max(v, 1.0))` |
| `CreateMask` (`:297`) | `input > background_threshold`, intersected with the user mask; `-bimodalT`, or `-auto_mask` on a non-Talairach volume, takes the threshold from EBTKS `Histogram::biModalThreshold` over `ceil(voxelMax-voxelMin+1)` bins (`volumeStats.cc:286-299`); `-auto_mask` on a Talairach volume uses the average brain mask, label-resampled |
| `mincmath -sub`/`-add`/`-mult` | buffer arithmetic |
| `sharpen_estimate` (`:500`) | `autoRange` → histogram class → `sharpenLookup` → `applyLookup`, then re-masked (`:527`) |
| `spline_smooth_volume` (`:567`) | `TBSplineVolume(domain, start={0,0,0}, step, sizes, distance, lambda)` with `domain` from `volume_domain` under `-full_support`; `addDataPoint(i,j,k,v)` over masked voxels stepping by `-subsample`; `fit()`; evaluate on the same grid, zero outside the mask. `-tp_spline` takes `createThinPlateSpline` and `addDataPoint(point, v)` instead |
| `field_CV` (`:702`) | population standard deviation of the field change inside the mask — the name is wrong, it is not a CV |
| `normalize_field_volume` (`:595`) | divide by the in-mask mean |
| `compact_spline_volume` (`:628`) | refit on the exponentiated field; keep the `Spline` in memory. `outputCompactField` runs only when a `.imp` was asked for |
| `nu_evaluate` auto mask (`:51`) | `mincstats -biModalT`: 2000-bin Otsu, returning the bin **centre** — a different rule from `volume_stats -biModalT` above, and both are needed |
| `evaluate_field` | a second `TBSplineVolume` on the **same domain** over the full-resolution grid, `setCoefficients`, evaluate. The two grids share voxel (0,0,0), so the spline's `index × step` coordinates carry across |
| `correct_field` | the double transcription of `smooth()` |
| field floor + `mincmath -div` | `output = input / max(field, floor)` |

The estimation grid and the output grid stay separate objects throughout: the field is
fitted on the shrunken grid and evaluated at full resolution.

## 3. Where this cannot match the Perl, and how to tell

Four deliberate divergences, all of them consequences of the task:

1. **No quantisation between stages.** Every intermediate is `double` rather than a MINC
   file. This is the point of the change and it moves the result.
2. **No `%lf` rounding** of the histogram, its domain, or the lookup table.
3. **`smooth()` in double** rather than float.
4. **`spline_smooth` reads volumes as float** (`loadFloatVolume`, `splineSmooth.cc:101`);
   the C++ fit sees double.

Add a `-legacy_rounding` flag that rounds the histogram counts, the range and the LUT to
six decimals, reproducing (2) alone. It is a verification instrument, not a feature: with
it on, a residual difference against the Perl is either the volume quantisation or a
mistake, and the two can be told apart.

**The stopping rule quantises everything downstream** (`change < 0.001`): two
implementations differing in the fifth decimal can run different numbers of iterations,
which moves the output by far more than any block difference. Every comparison below fixes
the iteration count first, and only then repeats at the default protocol.

## 4. Build

- New `ADD_EXECUTABLE(nu_correct_cxx ...)` in `legacy/N3/CMakeLists.txt` listing the
  `N3Pipeline/` sources, `Splines/{Spline,TBSpline}.cc`, `SplineSmooth/fieldIO.cc`,
  `VolumeHist/{DHistogram,WHistogram,GHistogram}.cc` and
  `$<TARGET_OBJECTS:sharpenhist_core>`.
- `INSTALL` it twice, as `nu_correct_cxx` and `nu_estimate_cxx`, so the argv[0] rule
  (`/nu_estimate/i` → estimate only) works as it does for the Perl. `-estimate_only` and
  `-correct` override it.
- The `LIBMINC_LIBRARIES` repair loop added in the previous task already lets this tree
  configure against the installed libminc.

## 5. Verification

Each new routine is checked against its own oracle before anything end to end. On
`tests/data/chunk.mnc` + `chunk_mask.mnc` (91×52×50, fast), driving the **installed**
programs in `/opt/minc/1.9.18.13/bin`:

| Routine | Oracle | Bound |
|---|---|---|
| `autoRange` + histogram | `volume_hist -bins 200 -auto_range -mask … -window` | exact to 1e-6, the `%lf` limit |
| `sharpenLookup` | `sharpen_hist -fwhm 0.15 -noise 0.01` on that histogram | 1e-6 |
| `applyLookup` | `minclookup -continuous -range` | 1e-6 |
| `shrink` | `mincresample -nearest_neighbour -nelements … -step …` | exact |
| `resampleLabel` | `resample_labels -resample -like` | exact |
| `biModalThreshold` (both) | `volume_stats -biModalT`, `mincstats -biModalT` | exact |
| spline fit | `spline_smooth -full_support -b_spline -lambda 1e-7 -mask` | fields to ~1e-6; **compare fields, never coefficients** — the normal equations are near-singular at 200 mm and the coefficients are determined to ~1e-4 by no solver |
| `smooth()` in double | `correct_field` | ~5e-6 relative RMS, the accuracy of the solve; the two sweep in different orders |

Then end to end:

1. `nu_correct` vs `nu_correct_cxx` on `chunk.mnc`, `-iterations 3 -stop 0` (fixed count),
   with and without `-legacy_rounding`. Report relative RMS and the per-iteration field
   change from both.
2. The same at the default protocol, reporting **both iteration counts**. A difference in
   count is the expected failure mode, not a defect.
3. `brain.mnc` at the default protocol, against `nu_correct` and against
   `legacy/N3/testing/brain_nu_ref.mnc.gz`. The reference will not be met at the legacy
   suite's `1e-4`: a double pipeline cannot reproduce a 16-bit-quantised one, and
   `torch_n3`'s own legacy backend sits 3.7e-3 from it.
4. `-tp_spline` on `chunk.mnc`, same protocol, against the Perl with `-tp_spline`.
5. `-estimate_only`: the `.imp` against the Perl's, compared by running `evaluate_field`
   on both and diffing the fields, not by diffing the text.
6. `-parzen_sigma 2` end to end, since that option was added to the C++ layer in the
   previous task and must reach the new driver too.

The harness is a script under `legacy/N3/testing/`, which is that project's own test
directory; it drives both pipelines and prints a table. It needs the installed MINC tools,
which is already true of everything else in there.

## 6. Order of work

1. `VolumeD` with load/save and the two resamplers; check §5's `shrink` and
   `resampleLabel` rows immediately — a wrong grid invalidates everything after it.
2. `MincTools`, `Sharpen`, and the histogram glue; check their §5 rows.
3. `FitField` (both spline types) and `SmoothField`; check their §5 rows.
4. `NuEstimate` — the iteration, the mask construction, the stopping rule.
5. `NuEvaluate`, the argument table, `main()`, and the two install names.
6. §5's end-to-end comparisons and the harness script.

Commit in `legacy/N3` at steps 1, 3, 5 and 6.

## 7. Files

New: `legacy/N3/src/N3Pipeline/{VolumeD,MincTools,Sharpen,SmoothField,FitField,NuEstimate,NuEvaluate}.{h,cc}`,
`legacy/N3/src/N3Pipeline/nu_correct_cxx.cc`, a comparison script under
`legacy/N3/testing/`.

Edited: `legacy/N3/CMakeLists.txt` only.

Untouched: every existing source in `legacy/N3/src` (they are linked, not modified), the
Perl drivers, `legacy/EBTKS`, everything under `/app/torch_n3`, `/opt/minc`.

## Not in scope, recorded

- Carrying back the port's other findings: the better-conditioned QR spline fit
  (`torch_n3/blocks/spline.py`, `solver="qr"`, cond 2.3e6 against the normal equations'
  5.3e12) and the `--denoise` prefilter. The C++ pipeline keeps `dsysv` on `AtA`.
- The EM and WM branches, `fir` smoothing, `-real`, `-differential`, `-initial`,
  `-islands`.
- Threading or GPU. The program stays single-threaded; the gain being sought here is the
  removal of file round trips, and a second change at the same time would make it
  unmeasurable.
- One divergence found while reading, not fixed and not this task's business:
  `torch_n3`'s `pipeline.py:122` resamples the user mask onto the estimation grid with
  **nearest neighbour**, where the Perl uses `resample_labels`, which is trilinear
  thresholded at 0.5. The C++ program follows the Perl.
