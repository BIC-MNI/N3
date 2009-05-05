# Makefile.am for N3
#
AUTOMAKE_OPTIONS = check-news subdir-objects
#

ACLOCAL_AMFLAGS = -I m4 -I .
#
AM_CPPFLAGS = -DUSE_COMPMAT -DUSE_DBLMAT
#
AM_LDFLAGS = 

SUBDIRS = testing

# Libraries which must be built (but not installed)
#
noinst_LIBRARIES = libmincprog.a

libmincprog_a_SOURCES = \
	src/MincProg/print_version.c

bin_PROGRAMS = \
	evaluate_field \
	extracttag \
	sharpen_hist \
	spline_smooth \
	volume_hist \
	volume_stats

LDADD = libmincprog.a
lapack_link_flags = $(LDADD) $(LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

m4_files = acx_blas.m4 acx_lapack.m4 \
	m4/mni_REQUIRE_LIB.m4 \
	m4/mni_REQUIRE_MNILIBS.m4 \
	m4/mni_REQUIRE_OPENINVENTOR.m4 \
	m4/mni_cxx_have_koenig_lookup.m4 \
	m4/smr_CGAL_MAKEFILE.m4 \
	m4/smr_OPTIONAL_LIB.m4 \
	m4/smr_REQUIRED_LIB.m4 \
	m4/smr_WITH_BUILD_PATH.m4

evaluate_field_SOURCES = \
	src/EvaluateField/evaluateField.cc \
	src/EvaluateField/evaluateFieldArgs.cc \
	src/SplineSmooth/fieldIO.cc
evaluate_field_LDADD = $(lapack_link_flags)

extracttag_SOURCES = \
	src/ExtractTag/extracttag.c

sharpen_hist_SOURCES = \
	src/SharpenHist/sharpen_hist.cc \
	src/SharpenHist/args.cc

spline_smooth_SOURCES = \
	src/SplineSmooth/splineSmooth.cc \
	src/SplineSmooth/splineSmoothArgs.cc \
	src/SplineSmooth/fieldIO.cc
spline_smooth_LDADD = $(lapack_link_flags)

volume_hist_SOURCES = \
	src/VolumeHist/minchist.cc \
	src/VolumeHist/DHistogram.cc \
	src/VolumeHist/WHistogram.cc \
	src/VolumeHist/args.cc

volume_stats_SOURCES = \
	src/VolumeStats/volumeStats.cc \
	src/VolumeStats/VolumeStatsArgs.cc

bin_SCRIPTS = \
	src/NUcorrect/field2imp \
	src/NUcorrect/imp2field \
	src/NUcorrect/make_template \
	src/NUcorrect/nu_correct \
	src/NUcorrect/nu_estimate \
	src/NUcorrect/nu_estimate_np_and_em \
	src/NUcorrect/nu_evaluate \
	src/NUcorrect/resample_labels \
	src/NUcorrect/sharpen_volume

# Yep, this is simply a copy (for backward compatibility)
src/NUcorrect/nu_correct: src/NUcorrect/nu_estimate
	cp -f -p src/NUcorrect/nu_estimate src/NUcorrect/nu_correct

noinst_HEADERS = \
	./src/EvaluateField/evaluateFieldArgs.h \
	./src/SharpenHist/args.h \
	./src/SplineSmooth/fieldIO.h \
	./src/SplineSmooth/splineSmooth.h \
	./src/SplineSmooth/splineSmoothArgs.h \
	./src/VolumeHist/DHistogram.h \
	./src/VolumeHist/WHistogram.h \
	./src/VolumeHist/args.h \
	./src/VolumeStats/VolumeStatsArgs.h

EXTRA_DIST = \
	INSTALL.GENERIC \
	README \
	NEWS \
	epm-header.in \
	include/version.h.in \
	src/NUcorrect/field2imp.in \
	src/NUcorrect/imp2field.in \
	src/NUcorrect/make_template.in \
	src/NUcorrect/nu_estimate.in \
	src/NUcorrect/nu_estimate_np_and_em.in \
	src/NUcorrect/nu_evaluate.in \
	src/NUcorrect/resample_labels.in \
	src/NUcorrect/sharpen_volume.in \
	$(m4_files)

# dist_model_DATA = model_data/N3/average_305_mask_1mm.mnc.gz

pkgdatadir = $(prefix)/share/N3

dist_pkgdata_DATA := model_data/N3/average_305_mask_1mm.mnc.gz

CLEANFILES = $(bin_SCRIPTS)

