all-test:
	cd $(srcdir) && chmod +x do_test

TESTS = do_test

check_SCRIPTS = do_test rms_diff

EXTRA_DIST = $(check_SCRIPTS) $(data_files) test_ref.log do_test.in rms_diff.in

CLEANFILES =

data_files = block.mnc.gz brain_mask.mnc.gz brain.mnc.gz brain_nu_ref.mnc.gz \
             chunk_mask.mnc.gz chunk.mnc.gz 

SEDFLAGS = -e 's,\@PERL\@,$(PERL),g'

do_test: do_test.in
	sed < $< $(SEDFLAGS) > $@
	chmod +x $@

rms_diff: rms_diff.in
	sed < $< $(SEDFLAGS) > $@
	chmod +x $@


