#! /bin/sh
# Originally from Steve Robbins, with small modications by Bert Vincent.
#

set -e

if [ ! -r m4/mni_REQUIRE_LIB.m4 ]; then
    cat <<EOF
The required m4 files were not found.
You need to check these out from their repository
using

    cvs -d /software/source checkout -d m4 libraries/mni-acmacros

(yes, two '-d' options)
Then re-run autogen.sh.
EOF
    exit 1
fi

cat <<EOF
Messages of the following type may be safely ignored.
Any other diagnostics may be a sign of trouble.  

    automake: configure.in: installing [...]
    warning: AC_TRY_RUN called without default to allow cross compiling

Let me (bert@bic.mni.mcgill.ca) know if something goes wrong!
EOF

test -d ac_config_aux || mkdir ac_config_aux

aclocal -I m4
autoheader
automake --add-missing --copy
autoconf

