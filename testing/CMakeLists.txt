# TODO add tests

SET(PERL "${PERL_EXECUTABLE}")

MACRO(ADD_PERL_TEST name command )
 ADD_TEST( ${name}
    ${CMAKE_CURRENT_BINARY_DIR}/../${command} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${ARGV7} ${ARGV8} ${ARGV9} ${ARGV10} ${ARGV11}
      ${ARGV12} ${ARGV13} ${ARGV14} ${ARGV15} ${ARGV16} ${ARGV17} ${ARGV18} ${ARGV19} ${ARGV20} ${ARGV21} ${ARGV22}
      ${ARGV23} ${ARGV24} ${ARGV25} ${ARGV26}
 )
 
 set_tests_properties( ${name} PROPERTIES ENVIRONMENT
  "PATH=${CMAKE_CURRENT_BINARY_DIR}/..:${CMAKE_CURRENT_BINARY_DIR}/../../minc/progs:$ENV{PATH};PERL5LIB=${CMAKE_SOURCE_DIR}/mni-perllib:${CMAKE_CURRENT_BINARY_DIR}/../../mni-perllib/perl:${CMAKE_CURRENT_BINARY_DIR}/../../perl-modules:$ENV{PERL5LIB}" 
 )
ENDMACRO(ADD_PERL_TEST)


MACRO(ADD_CUSTOM_PERL_TEST name command )

 ADD_TEST( ${name}
    ${CMAKE_CURRENT_SOURCE_DIR}/${command} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${ARGV7} ${ARGV8} ${ARGV9} ${ARGV10} ${ARGV11}
      ${ARGV12} ${ARGV13} ${ARGV14} ${ARGV15} ${ARGV16} ${ARGV17} ${ARGV18} ${ARGV19} ${ARGV20} ${ARGV21} ${ARGV22}
      ${ARGV23} ${ARGV24} ${ARGV25} ${ARGV26}
 )
 
 set_tests_properties( ${name} PROPERTIES ENVIRONMENT
  "PATH=${CMAKE_CURRENT_BINARY_DIR}/..:${CMAKE_CURRENT_BINARY_DIR}/../../minc/progs:$ENV{PATH};PERL5LIB=${CMAKE_SOURCE_DIR}/mni-perllib:${CMAKE_CURRENT_BINARY_DIR}/../../mni-perllib/perl:${CMAKE_CURRENT_BINARY_DIR}/../../perl-modules:$ENV{PERL5LIB}" 
 )
 
ENDMACRO(ADD_CUSTOM_PERL_TEST)


#ADD_PERL_TEST2(nu_env "/usr/bin/env")

ADD_PERL_TEST(nu_field2imp_1 field2imp ${CMAKE_CURRENT_SOURCE_DIR}/block.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/test.imp -clobber)
ADD_PERL_TEST(nu_field2imp_2 field2imp ${CMAKE_CURRENT_SOURCE_DIR}/block.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/test.imp -clobber)

ADD_PERL_TEST(nu_imp2field_1 imp2field ${CMAKE_CURRENT_BINARY_DIR}/test.imp -like ${CMAKE_CURRENT_SOURCE_DIR}/block.mnc.gz -shrink 2 ${CMAKE_CURRENT_BINARY_DIR}/test.mnc -clobber)
ADD_PERL_TEST(nu_imp2field_2 imp2field ${CMAKE_CURRENT_BINARY_DIR}/test.imp -like ${CMAKE_CURRENT_SOURCE_DIR}/block.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/test.mnc -clobber)

ADD_PERL_TEST(nu_sharpen_volume_1 sharpen_volume -fwhm 0.15 -noise 0.01 -bins 150 -parzen -verbose ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk_sharp.mnc -clobber)
ADD_PERL_TEST(nu_sharpen_volume_2 sharpen_volume -fwhm 0.15 -noise 0.01 -bins 150 -deblur -save_histogram ${CMAKE_CURRENT_BINARY_DIR}/sharp.hist -parzen -verbose ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk_sharp.mnc -clobber)


#ADD_PERL_TEST(nu_estimate_np_and_em_1 nu_estimate_np_and_em -help)
ADD_PERL_TEST(nu_estimate_np_and_em_2 nu_estimate_np_and_em -distance 100 -b_spline 1 -spline_subsample 2 -auto_mask -iterations 2 -stop 0.001 -shrink 2 -sharpen 0.15 0.01 -parzen  ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/brain.imp -verbose -clobber)
ADD_PERL_TEST(nu_estimate_np_and_em_3 nu_estimate_np_and_em -distance 100 -b_spline 1 -spline_subsample 2 -auto_mask -iterations 2 3 -stop 0.01 0.001 -shrink 2 -sharpen 0.15 0.01 -parzen ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -verbose -clobber)
ADD_PERL_TEST(nu_estimate_np_and_em_4 nu_estimate_np_and_em -tp_spline 1 -bimodalT -iterations 2 3 -stop 0.01 0.001 -shrink 2 -sharpen 0.15 0.01 -parzen ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -verbose -clobber -mask ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz )
#ADD_PERL_TEST(nu_estimate_np_and_em -fwhm 100 -fir -spline_subsample 2 -background 100000 -iterations 2 3 -stop 0.01 0.001 -shrink 2 -sharpen 0.15 0.01 -parzen ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -bins 100 -clobber)

#ADD_PERL_TEST(nu_evaluate_1 nu_evaluate -help)
ADD_PERL_TEST(nu_evaluate_2 nu_evaluate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz -mapping ${CMAKE_CURRENT_BINARY_DIR}/brain.imp ${CMAKE_CURRENT_BINARY_DIR}/brain_nu.mnc -clobber )

#ADD_PERL_TEST(nu_estimate_1 nu_estimate -help)
#ADD_PERL_TEST(nu_estimate_2 ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/brain.imp -clobber -verbose -V0.9)
#ADD_PERL_TEST(nu_estimate_3 ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/brain.imp -clobber -verbose -newest)
ADD_PERL_TEST(nu_estimate_4 nu_estimate -iterations 1 -stop 0.001 ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -clobber -verbose -mask ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz -clobber)
ADD_PERL_TEST(nu_estimate_5 nu_estimate -iterations 1 -stop 0.001 ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -clobber -verbose -options "-deblur -log" -clobber )
ADD_PERL_TEST(nu_estimate_6 nu_estimate ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -clobber -verbose  -stop 0.001 -iterations 1 -shrink 4  -distance 100 )

#ADD_PERL_TEST(nu_correct_7 nu_correct -help)
ADD_PERL_TEST(nu_correct_8  nu_correct -fwhm 0.2 -mapping_dir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk_nu.mnc -mask ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz -shrink 3 -verbose -stop 0.001 -iterations 1 -clobber) 

ADD_PERL_TEST(nu_estimate_7 nu_estimate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/brain.imp -clobber)
ADD_PERL_TEST(nu_evaluate_8 nu_evaluate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz -mapping ${CMAKE_CURRENT_BINARY_DIR}/brain.imp ${CMAKE_CURRENT_SOURCE_DIR}/brain_nu.mnc -clobber)

ADD_CUSTOM_PERL_TEST(nu_reference_1 compare_nu_result.pl ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz -mask ${CMAKE_CURRENT_SOURCE_DIR}/../model_data/N3/icbm_avg_152_t1_tal_nlin_symmetric_VI_mask.mnc.gz ${CMAKE_CURRENT_SOURCE_DIR}/brain_nu_ref.mnc.gz 0.0001)
