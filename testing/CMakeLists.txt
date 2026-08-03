# TODO add tests
SET(PERL "${PERL_EXECUTABLE}")

MACRO(ADD_PERL_TEST name command )
 ADD_TEST( ${name}
    ${CMAKE_CURRENT_BINARY_DIR}/../${command} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${ARGV7} ${ARGV8} ${ARGV9} ${ARGV10} ${ARGV11}
      ${ARGV12} ${ARGV13} ${ARGV14} ${ARGV15} ${ARGV16} ${ARGV17} ${ARGV18} ${ARGV19} ${ARGV20} ${ARGV21} ${ARGV22}
      ${ARGV23} ${ARGV24} ${ARGV25} ${ARGV26}
 )

 IF(MINC_TEST_ENVIRONMENT)
  set_tests_properties( ${name} PROPERTIES ENVIRONMENT "${MINC_TEST_ENVIRONMENT}")
 ENDIF(MINC_TEST_ENVIRONMENT)
ENDMACRO(ADD_PERL_TEST)


MACRO(ADD_CUSTOM_PERL_TEST name command )

 ADD_TEST( ${name}
    ${CMAKE_CURRENT_SOURCE_DIR}/${command} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${ARGV7} ${ARGV8} ${ARGV9} ${ARGV10} ${ARGV11}
      ${ARGV12} ${ARGV13} ${ARGV14} ${ARGV15} ${ARGV16} ${ARGV17} ${ARGV18} ${ARGV19} ${ARGV20} ${ARGV21} ${ARGV22}
      ${ARGV23} ${ARGV24} ${ARGV25} ${ARGV26}
 )
 
 IF(MINC_TEST_ENVIRONMENT)
  set_tests_properties( ${name} PROPERTIES ENVIRONMENT "${MINC_TEST_ENVIRONMENT}")
 ENDIF(MINC_TEST_ENVIRONMENT)
ENDMACRO(ADD_CUSTOM_PERL_TEST)


ADD_PERL_TEST(nu_field2imp_1 field2imp ${CMAKE_CURRENT_SOURCE_DIR}/block.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/test.imp -clobber)
ADD_PERL_TEST(nu_field2imp_2 field2imp ${CMAKE_CURRENT_SOURCE_DIR}/block.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/test.imp -clobber)

ADD_PERL_TEST(nu_imp2field_1 imp2field ${CMAKE_CURRENT_BINARY_DIR}/test.imp -like ${CMAKE_CURRENT_SOURCE_DIR}/block.mnc.gz -shrink 2 ${CMAKE_CURRENT_BINARY_DIR}/test.mnc -clobber)
ADD_PERL_TEST(nu_imp2field_2 imp2field ${CMAKE_CURRENT_BINARY_DIR}/test.imp -like ${CMAKE_CURRENT_SOURCE_DIR}/block.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/test.mnc -clobber)

ADD_PERL_TEST(nu_sharpen_volume_1 sharpen_volume -fwhm 0.15 -noise 0.01 -bins 150 -parzen -verbose ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk_sharp.mnc -clobber -verbose)
ADD_PERL_TEST(nu_sharpen_volume_2 sharpen_volume -fwhm 0.15 -noise 0.01 -bins 150 -deblur -save_histogram ${CMAKE_CURRENT_BINARY_DIR}/sharp.hist -parzen -verbose ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk_sharp.mnc -clobber)

#ADD_PERL_TEST(nu_estimate_np_and_em_1 nu_estimate_np_and_em -help)
ADD_PERL_TEST(nu_estimate_np_and_em_2 nu_estimate_np_and_em -distance 100 -b_spline 1 -spline_subsample 2 -auto_mask -iterations 2 -stop 0.001 -shrink 2 -sharpen 0.15 0.01 -parzen  ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/brain.imp -verbose -clobber)
ADD_PERL_TEST(nu_estimate_np_and_em_3 nu_estimate_np_and_em -distance 100 -b_spline 1 -spline_subsample 2 -auto_mask -iterations 2 3 -stop 0.01 0.001 -shrink 2 -sharpen 0.15 0.01 -parzen ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -verbose -clobber)
ADD_PERL_TEST(nu_estimate_np_and_em_4 nu_estimate_np_and_em -tp_spline 1 -bimodalT -iterations 2 3 -stop 0.01 0.001 -shrink 2 -sharpen 0.15 0.01 -parzen ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -verbose -clobber -mask ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz -verbose )
ADD_PERL_TEST(nu_estimate_np_and_em_5 nu_estimate_np_and_em -fwhm 100 -fir -spline_subsample 2 -background 100000 -iterations 2 3 -stop 0.01 0.001 -shrink 2 -sharpen 0.15 0.01 -parzen ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -bins 100 -clobber -verbose )

#ADD_PERL_TEST(nu_evaluate_1 nu_evaluate -help)
ADD_PERL_TEST(nu_evaluate_2 nu_evaluate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz -mapping ${CMAKE_CURRENT_BINARY_DIR}/brain.imp ${CMAKE_CURRENT_BINARY_DIR}/brain_nu.mnc -clobber )

#ADD_PERL_TEST(nu_estimate_1 nu_estimate -help)
#ADD_PERL_TEST(nu_estimate_2 ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/brain.imp -clobber -verbose -V0.9)
#ADD_PERL_TEST(nu_estimate_3 ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/brain.imp -clobber -verbose -newest)
ADD_PERL_TEST(nu_estimate_4 nu_estimate -iterations 1 -stop 0.001 ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -clobber -verbose -mask ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz -clobber)
ADD_PERL_TEST(nu_estimate_5 nu_estimate -iterations 1 -stop 0.001 ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -clobber -verbose -options "-deblur -log" -clobber )
ADD_PERL_TEST(nu_estimate_6 nu_estimate ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -clobber -verbose  -stop 0.001 -iterations 1 -shrink 4  -distance 100 )

#ADD_PERL_TEST(nu_correct_7 nu_correct -help)
ADD_PERL_TEST(nu_correct_8  nu_correct -fwhm 0.2 -mapping_dir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk_nu.mnc -mask ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz -shrink 3 -verbose -stop 0.001 -iterations 1 -clobber  -verbose ) 
ADD_PERL_TEST(nu_estimate_7 nu_estimate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/brain.imp -clobber -verbose )
ADD_PERL_TEST(nu_evaluate_8 nu_evaluate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz -mapping ${CMAKE_CURRENT_BINARY_DIR}/brain.imp ${CMAKE_CURRENT_BINARY_DIR}/brain_nu.mnc -clobber  -verbose )

ADD_CUSTOM_PERL_TEST(nu_reference_1 compare_nu_result.pl ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz -mask ${CMAKE_CURRENT_SOURCE_DIR}/../model_data/N3/icbm_avg_152_t1_tal_nlin_symmetric_VI_mask.mnc.gz ${CMAKE_CURRENT_SOURCE_DIR}/brain_nu_ref.mnc.gz 0.0001  --verbose)

# ---- Test dependencies and resource locks for parallel ctest ----
# Several tests write to the same temp files in CMAKE_CURRENT_BINARY_DIR.
# Use DEPENDS for sequential producer/consumer relationships and
# RESOURCE_LOCK to prevent concurrent writes to the same files.

# test.imp / test.mnc: field2imp produces, imp2field consumes
set_tests_properties(nu_imp2field_1 PROPERTIES DEPENDS "nu_field2imp_1")
set_tests_properties(nu_imp2field_2 PROPERTIES DEPENDS "nu_field2imp_2")
set_tests_properties(nu_field2imp_1 nu_field2imp_2
                     nu_imp2field_1 nu_imp2field_2
                     PROPERTIES RESOURCE_LOCK "n3_test_imp")

# chunk_sharp.mnc: both sharpen_volume tests write to same file
set_tests_properties(nu_sharpen_volume_1 nu_sharpen_volume_2
                     PROPERTIES RESOURCE_LOCK "n3_chunk_sharp")

# brain.imp / brain_nu.mnc: estimate produces, evaluate consumes
set_tests_properties(nu_evaluate_2 PROPERTIES DEPENDS "nu_estimate_np_and_em_2")
set_tests_properties(nu_evaluate_8 PROPERTIES DEPENDS "nu_estimate_7")
set_tests_properties(nu_estimate_np_and_em_2 nu_evaluate_2
                     nu_estimate_7 nu_evaluate_8
                     PROPERTIES RESOURCE_LOCK "n3_brain_imp")

# chunk.imp: multiple tests write to same file
set_tests_properties(nu_estimate_np_and_em_3 nu_estimate_np_and_em_4
                     nu_estimate_np_and_em_5
                     nu_estimate_4 nu_estimate_5 nu_estimate_6
                     nu_correct_8
                     PROPERTIES RESOURCE_LOCK "n3_chunk_imp")

# ---- N3Pipeline: the in-memory C++ reimplementation of the Perl drivers ----
#
# These are plain executables, not Perl.  Each links the same object libraries
# as nu_correct_cxx and reads its oracle answers from n3pipeline/reference/,
# recorded once by n3pipeline/regenerate_reference.sh.

MACRO(ADD_N3PIPELINE_TEST name)
  ADD_EXECUTABLE(${name} n3pipeline/${name}.cc $<TARGET_OBJECTS:n3pipeline_core> ${ARGN})
  TARGET_INCLUDE_DIRECTORIES(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/n3pipeline)
  TARGET_COMPILE_DEFINITIONS(${name} PRIVATE
    N3_REFERENCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/n3pipeline/reference"
    N3_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  ADD_TEST(n3cxx_${name} ${CMAKE_CURRENT_BINARY_DIR}/${name})
  IF(MINC_TEST_ENVIRONMENT)
    set_tests_properties(n3cxx_${name} PROPERTIES ENVIRONMENT "${MINC_TEST_ENVIRONMENT}")
  ENDIF(MINC_TEST_ENVIRONMENT)
ENDMACRO(ADD_N3PIPELINE_TEST)

ADD_N3PIPELINE_TEST(test_harness)
ADD_N3PIPELINE_TEST(test_shrink)
ADD_N3PIPELINE_TEST(test_resample_label)
ADD_N3PIPELINE_TEST(test_buffers)

# A failing check must exit non-zero.  Asserted permanently rather than
# demonstrated once and deleted.
ADD_TEST(n3cxx_test_harness_red ${CMAKE_CURRENT_BINARY_DIR}/test_harness --fail)
set_tests_properties(n3cxx_test_harness_red PROPERTIES WILL_FAIL TRUE)
