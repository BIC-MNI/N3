# TODO add tests

IF(MINC_TOOLKIT_BUILD)

SET(CTEST_ENVIRONMENT 
 "PATH=${CMAKE_CURRENT_BINARY_DIR}/..:$ENV{PATH}" 
 "PERL5LIB=${CMAKE_SOURCE_DIR}/mni-perllib:${CMAKE_CURRENT_BINARY_DIR}/../mni-perllib:$ENV{PERL5LIB}"
)
ENDIF(MINC_TOOLKIT_BUILD)


ADD_TEST(field2imp_1 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../field2imp ${CMAKE_CURRENT_SOURCE_DIR}/block.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/test.imp -clobber)
ADD_TEST(field2imp_2 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../field2imp ${CMAKE_CURRENT_SOURCE_DIR}/block.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/test.imp -clobber)

ADD_TEST(imp2field_1 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../imp2field ${CMAKE_CURRENT_BINARY_DIR}/test.imp -like ${CMAKE_CURRENT_SOURCE_DIR}/block.mnc.gz -shrink 2 ${CMAKE_CURRENT_BINARY_DIR}/test.mnc -clobber)
ADD_TEST(imp2field_2 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../imp2field ${CMAKE_CURRENT_BINARY_DIR}/test.imp -like ${CMAKE_CURRENT_SOURCE_DIR}/block.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/test.mnc -clobber)

ADD_TEST(sharpen_volume_1 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../sharpen_volume -fwhm 0.15 -noise 0.01 -bins 150 -parzen -verbose ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk_sharp.mnc -clobber)
ADD_TEST(sharpen_volume_2 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../sharpen_volume -fwhm 0.15 -noise 0.01 -bins 150 -deblur -save_histogram ${CMAKE_CURRENT_BINARY_DIR}/sharp.hist -parzen -verbose ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk_sharp.mnc -clobber)


#ADD_TEST(nu_estimate_np_and_em_1 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate_np_and_em -help)
ADD_TEST(nu_estimate_np_and_em_2 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate_np_and_em -distance 100 -b_spline 1 -spline_subsample 2 -auto_mask -iterations 2 -stop 0.001 -shrink 2 -sharpen 0.15 0.01 -parzen  ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/brain.imp -verbose -clobber)
ADD_TEST(nu_estimate_np_and_em_3 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate_np_and_em -distance 100 -b_spline 1 -spline_subsample 2 -auto_mask -iterations 2 3 -stop 0.01 0.001 -shrink 2 -sharpen 0.15 0.01 -parzen ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -verbose -clobber)
ADD_TEST(nu_estimate_np_and_em_4 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate_np_and_em -tp_spline 1 -bimodalT -iterations 2 3 -stop 0.01 0.001 -shrink 2 -sharpen 0.15 0.01 -parzen ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -verbose -clobber -mask ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz)
#ADD_TEST(nu_estimate_np_and_em -fwhm 100 -fir -spline_subsample 2 -background 100000 -iterations 2 3 -stop 0.01 0.001 -shrink 2 -sharpen 0.15 0.01 -parzen ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -bins 100 -clobber)

#ADD_TEST(nu_evaluate_1 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_evaluate -help)
ADD_TEST(nu_evaluate_2 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_evaluate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz -mapping ${CMAKE_CURRENT_BINARY_DIR}/brain.imp ${CMAKE_CURRENT_BINARY_DIR}/brain_nu.mnc -clobber )

#ADD_TEST(nu_estimate_1 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate -help)
#ADD_TEST(nu_estimate_2 ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/brain.imp -clobber -verbose -V0.9)
#ADD_TEST(nu_estimate_3 ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/brain.imp -clobber -verbose -newest)
ADD_TEST(nu_estimate_4 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate -iterations 1 -stop 0.001 ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -clobber -verbose -mask ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz -clobber)
ADD_TEST(nu_estimate_5 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate -iterations 1 -stop 0.001 ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -clobber -verbose -options "-deblur -log" -clobber )
ADD_TEST(nu_estimate_6 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk.imp -clobber -verbose  -stop 0.001 -iterations 1 -shrink 4  -distance 100 )

#ADD_TEST(nu_correct_7 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_correct -help)
ADD_TEST(nu_correct_8 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_correct -fwhm 0.2 -mapping_dir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/chunk.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/chunk_nu.mnc -mask ${CMAKE_CURRENT_SOURCE_DIR}/chunk_mask.mnc.gz -shrink 3 -verbose -stop 0.001 -iterations 1 -clobber) 

ADD_TEST(nu_estimate_7 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_estimate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz ${CMAKE_CURRENT_BINARY_DIR}/brain.imp -clobber)
ADD_TEST(nu_evaluate_8 ${PERL_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/../nu_evaluate ${CMAKE_CURRENT_SOURCE_DIR}/brain.mnc.gz -mapping ${CMAKE_CURRENT_BINARY_DIR}/brain.imp ${CMAKE_CURRENT_SOURCE_DIR}/brain_nu.mnc -clobber)
