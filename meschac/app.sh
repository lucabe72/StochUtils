#! /bin/sh

# Guess values for system-dependent variables and create Makefiles.
# Generated automatically using autoconf version 2.13 
# Copyright (C) 1992, 93, 94, 95, 96 Free Software Foundation, Inc.
#
# This configure script is free software; the Free Software Foundation
# gives unlimited permission to copy, distribute and modify it.

# Defaults:
ac_help="$ac_help
cat >> confdefs.h <<\EOF
#define REAL_DBL 1
EOF
"

# Initialize some variables set by options.
# The variables have the same names as the options, with
# dashes changed to underlines.
build=NONE
cache_file=./config.cache
exec_prefix=NONE
host=NONE
no_create=
nonopt=NONE
no_recursion=
prefix=NONE
program_prefix=NONE
program_suffix=NONE
program_transform_name=s,x,x,
silent=
site=
srcdir=
target=NONE
verbose=
x_includes=NONE
x_libraries=NONE
bindir='${exec_prefix}/bin'
sbindir='${exec_prefix}/sbin'
libexecdir='${exec_prefix}/libexec'
datadir='${prefix}/share'
sysconfdir='${prefix}/etc'
sharedstatedir='${prefix}/com'
localstatedir='${prefix}/var'
libdir='${exec_prefix}/lib'
includedir='${prefix}/include'
oldincludedir='/usr/include'
infodir='${prefix}/info'
mandir='${prefix}/man'

# Initialize some other variables.
subdirs=
MFLAGS= MAKEFLAGS=
SHELL=${CONFIG_SHELL-/bin/sh}
# Maximum number of lines to put in a shell here document.
ac_max_here_lines=12


# Extract the first word of "gcc", so it can be a program name with args.
set dummy gcc; ac_word=$2
echo $ac_n "checking for $ac_word""... $ac_c" 1>&6
echo "configure:635: checking for $ac_word" >&5
if eval "test \"`echo '$''{'ac_cv_prog_CC'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  if test -n "$CC"; then
  ac_cv_prog_CC="$CC" # Let the user override the test.
else
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS=":"
  ac_dummy="$PATH"
  for ac_dir in $ac_dummy; do
    test -z "$ac_dir" && ac_dir=.
    if test -f $ac_dir/$ac_word; then
      ac_cv_prog_CC="gcc"
      break
    fi
  done
  IFS="$ac_save_ifs"
fi
fi
CC="$ac_cv_prog_CC"
if test -n "$CC"; then
  echo "$ac_t""$CC" 1>&6
else
  echo "$ac_t""no" 1>&6
fi

if test -z "$CC"; then
  # Extract the first word of "cc", so it can be a program name with args.
set dummy cc; ac_word=$2
echo $ac_n "checking for $ac_word""... $ac_c" 1>&6
echo "configure:665: checking for $ac_word" >&5
if eval "test \"`echo '$''{'ac_cv_prog_CC'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  if test -n "$CC"; then
  ac_cv_prog_CC="$CC" # Let the user override the test.
else
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS=":"
  ac_prog_rejected=no
  ac_dummy="$PATH"
  for ac_dir in $ac_dummy; do
    test -z "$ac_dir" && ac_dir=.
    if test -f $ac_dir/$ac_word; then
      if test "$ac_dir/$ac_word" = "/usr/ucb/cc"; then
        ac_prog_rejected=yes
	continue
      fi
      ac_cv_prog_CC="cc"
      break
    fi
  done
  IFS="$ac_save_ifs"
if test $ac_prog_rejected = yes; then
  # We found a bogon in the path, so make sure we never use it.
  set dummy $ac_cv_prog_CC
  shift
  if test $# -gt 0; then
    # We chose a different compiler from the bogus one.
    # However, it has the same basename, so the bogon will be chosen
    # first if we set CC to just the basename; use the full file name.
    shift
    set dummy "$ac_dir/$ac_word" "$@"
    shift
    ac_cv_prog_CC="$@"
  fi
fi
fi
fi
CC="$ac_cv_prog_CC"
if test -n "$CC"; then
  echo "$ac_t""$CC" 1>&6
else
  echo "$ac_t""no" 1>&6
fi

  if test -z "$CC"; then
    case "`uname -s`" in
    *win32* | *WIN32*)
      # Extract the first word of "cl", so it can be a program name with args.
set dummy cl; ac_word=$2
echo $ac_n "checking for $ac_word""... $ac_c" 1>&6
echo "configure:716: checking for $ac_word" >&5
if eval "test \"`echo '$''{'ac_cv_prog_CC'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  if test -n "$CC"; then
  ac_cv_prog_CC="$CC" # Let the user override the test.
else
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS=":"
  ac_dummy="$PATH"
  for ac_dir in $ac_dummy; do
    test -z "$ac_dir" && ac_dir=.
    if test -f $ac_dir/$ac_word; then
      ac_cv_prog_CC="cl"
      break
    fi
  done
  IFS="$ac_save_ifs"
fi
fi
CC="$ac_cv_prog_CC"
if test -n "$CC"; then
  echo "$ac_t""$CC" 1>&6
else
  echo "$ac_t""no" 1>&6
fi
 ;;
    esac
  fi
  test -z "$CC" && { echo "configure: error: no acceptable cc found in \$PATH" 1>&2; exit 1; }
fi

echo $ac_n "checking whether the C compiler ($CC $CFLAGS $LDFLAGS) works""... $ac_c" 1>&6
echo "configure:748: checking whether the C compiler ($CC $CFLAGS $LDFLAGS) works" >&5

ac_ext=c
# CFLAGS is not in ac_cpp because -g, -O, etc. are not valid cpp options.
ac_cpp='$CPP $CPPFLAGS'
ac_compile='${CC-cc} -c $CFLAGS $CPPFLAGS conftest.$ac_ext 1>&5'
ac_link='${CC-cc} -o conftest${ac_exeext} $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
cross_compiling=$ac_cv_prog_cc_cross

cat > conftest.$ac_ext << EOF

#line 759 "configure"
#include "confdefs.h"

main(){return(0);}
EOF
if { (eval echo configure:764: \"$ac_link\") 1>&5; (eval $ac_link) 2>&5; } && test -s conftest${ac_exeext}; then
  ac_cv_prog_cc_works=yes
  # If we can't run a trivial program, we are probably using a cross compiler.
  if (./conftest; exit) 2>/dev/null; then
    ac_cv_prog_cc_cross=no
  else
    ac_cv_prog_cc_cross=yes
  fi
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  ac_cv_prog_cc_works=no
fi
rm -fr conftest*
ac_ext=c
# CFLAGS is not in ac_cpp because -g, -O, etc. are not valid cpp options.
ac_cpp='$CPP $CPPFLAGS'
ac_compile='${CC-cc} -c $CFLAGS $CPPFLAGS conftest.$ac_ext 1>&5'
ac_link='${CC-cc} -o conftest${ac_exeext} $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
cross_compiling=$ac_cv_prog_cc_cross

echo "$ac_t""$ac_cv_prog_cc_works" 1>&6
if test $ac_cv_prog_cc_works = no; then
  { echo "configure: error: installation or configuration problem: C compiler cannot create executables." 1>&2; exit 1; }
fi
echo $ac_n "checking whether the C compiler ($CC $CFLAGS $LDFLAGS) is a cross-compiler""... $ac_c" 1>&6
echo "configure:790: checking whether the C compiler ($CC $CFLAGS $LDFLAGS) is a cross-compiler" >&5
echo "$ac_t""$ac_cv_prog_cc_cross" 1>&6
cross_compiling=$ac_cv_prog_cc_cross

echo $ac_n "checking whether we are using GNU C""... $ac_c" 1>&6
echo "configure:795: checking whether we are using GNU C" >&5
if eval "test \"`echo '$''{'ac_cv_prog_gcc'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  cat > conftest.c <<EOF
#ifdef __GNUC__
  yes;
#endif
EOF
if { ac_try='${CC-cc} -E conftest.c'; { (eval echo configure:804: \"$ac_try\") 1>&5; (eval $ac_try) 2>&5; }; } | egrep yes >/dev/null 2>&1; then
  ac_cv_prog_gcc=yes
else
  ac_cv_prog_gcc=no
fi
fi

echo "$ac_t""$ac_cv_prog_gcc" 1>&6

if test $ac_cv_prog_gcc = yes; then
  GCC=yes
else
  GCC=
fi

ac_test_CFLAGS="${CFLAGS+set}"
ac_save_CFLAGS="$CFLAGS"
CFLAGS=
echo $ac_n "checking whether ${CC-cc} accepts -g""... $ac_c" 1>&6
echo "configure:823: checking whether ${CC-cc} accepts -g" >&5
if eval "test \"`echo '$''{'ac_cv_prog_cc_g'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  echo 'void f(){}' > conftest.c
if test -z "`${CC-cc} -g -c conftest.c 2>&1`"; then
  ac_cv_prog_cc_g=yes
else
  ac_cv_prog_cc_g=no
fi
rm -f conftest*

fi

echo "$ac_t""$ac_cv_prog_cc_g" 1>&6
if test "$ac_test_CFLAGS" = set; then
  CFLAGS="$ac_save_CFLAGS"
elif test $ac_cv_prog_cc_g = yes; then
  if test "$GCC" = yes; then
    CFLAGS="-g -O2"
  else
    CFLAGS="-g"
  fi
else
  if test "$GCC" = yes; then
    CFLAGS="-O2"
  else
    CFLAGS=
  fi
fi

echo $ac_n "checking how to run the C preprocessor""... $ac_c" 1>&6
echo "configure:855: checking how to run the C preprocessor" >&5
# On Suns, sometimes $CPP names a directory.
if test -n "$CPP" && test -d "$CPP"; then
  CPP=
fi
if test -z "$CPP"; then
if eval "test \"`echo '$''{'ac_cv_prog_CPP'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
    # This must be in double quotes, not single quotes, because CPP may get
  # substituted into the Makefile and "${CC-cc}" will confuse make.
  CPP="${CC-cc} -E"
  # On the NeXT, cc -E runs the code through the compiler's parser,
  # not just through cpp.
  cat > conftest.$ac_ext <<EOF
#line 870 "configure"
#include "confdefs.h"
#include <assert.h>
Syntax Error
EOF
ac_try="$ac_cpp conftest.$ac_ext >/dev/null 2>conftest.out"
{ (eval echo configure:876: \"$ac_try\") 1>&5; (eval $ac_try) 2>&5; }
ac_err=`grep -v '^ *+' conftest.out | grep -v "^conftest.${ac_ext}\$"`
if test -z "$ac_err"; then
  :
else
  echo "$ac_err" >&5
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  CPP="${CC-cc} -E -traditional-cpp"
  cat > conftest.$ac_ext <<EOF
#line 887 "configure"
#include "confdefs.h"
#include <assert.h>
Syntax Error
EOF
ac_try="$ac_cpp conftest.$ac_ext >/dev/null 2>conftest.out"
{ (eval echo configure:893: \"$ac_try\") 1>&5; (eval $ac_try) 2>&5; }
ac_err=`grep -v '^ *+' conftest.out | grep -v "^conftest.${ac_ext}\$"`
if test -z "$ac_err"; then
  :
else
  echo "$ac_err" >&5
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  CPP="${CC-cc} -nologo -E"
  cat > conftest.$ac_ext <<EOF
#line 904 "configure"
#include "confdefs.h"
#include <assert.h>
Syntax Error
EOF
ac_try="$ac_cpp conftest.$ac_ext >/dev/null 2>conftest.out"
{ (eval echo configure:910: \"$ac_try\") 1>&5; (eval $ac_try) 2>&5; }
ac_err=`grep -v '^ *+' conftest.out | grep -v "^conftest.${ac_ext}\$"`
if test -z "$ac_err"; then
  :
else
  echo "$ac_err" >&5
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  CPP=/lib/cpp
fi
rm -f conftest*
fi
rm -f conftest*
fi
rm -f conftest*
  ac_cv_prog_CPP="$CPP"
fi
  CPP="$ac_cv_prog_CPP"
else
  ac_cv_prog_CPP="$CPP"
fi
echo "$ac_t""$CPP" 1>&6

echo $ac_n "checking for AIX""... $ac_c" 1>&6
echo "configure:935: checking for AIX" >&5
cat > conftest.$ac_ext <<EOF
#line 937 "configure"
#include "confdefs.h"
#ifdef _AIX
  yes
#endif

EOF
if (eval "$ac_cpp conftest.$ac_ext") 2>&5 |
  egrep "yes" >/dev/null 2>&1; then
  rm -rf conftest*
  echo "$ac_t""yes" 1>&6; cat >> confdefs.h <<\EOF
#define _ALL_SOURCE 1
EOF

else
  rm -rf conftest*
  echo "$ac_t""no" 1>&6
fi
rm -f conftest*


AC_HPUX
ac_safe=`echo "minix/config.h" | sed 'y%./+-%__p_%'`
echo $ac_n "checking for minix/config.h""... $ac_c" 1>&6
echo "configure:961: checking for minix/config.h" >&5
if eval "test \"`echo '$''{'ac_cv_header_$ac_safe'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  cat > conftest.$ac_ext <<EOF
#line 966 "configure"
#include "confdefs.h"
#include <minix/config.h>
EOF
ac_try="$ac_cpp conftest.$ac_ext >/dev/null 2>conftest.out"
{ (eval echo configure:971: \"$ac_try\") 1>&5; (eval $ac_try) 2>&5; }
ac_err=`grep -v '^ *+' conftest.out | grep -v "^conftest.${ac_ext}\$"`
if test -z "$ac_err"; then
  rm -rf conftest*
  eval "ac_cv_header_$ac_safe=yes"
else
  echo "$ac_err" >&5
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  eval "ac_cv_header_$ac_safe=no"
fi
rm -f conftest*
fi
if eval "test \"`echo '$ac_cv_header_'$ac_safe`\" = yes"; then
  echo "$ac_t""yes" 1>&6
  MINIX=yes
else
  echo "$ac_t""no" 1>&6
MINIX=
fi

if test "$MINIX" = yes; then
  cat >> confdefs.h <<\EOF
#define _POSIX_SOURCE 1
EOF

  cat >> confdefs.h <<\EOF
#define _POSIX_1_SOURCE 2
EOF

  cat >> confdefs.h <<\EOF
#define _MINIX 1
EOF

fi

echo $ac_n "checking for POSIXized ISC""... $ac_c" 1>&6
echo "configure:1009: checking for POSIXized ISC" >&5
if test -d /etc/conf/kconfig.d &&
  grep _POSIX_VERSION /usr/include/sys/unistd.h >/dev/null 2>&1
then
  echo "$ac_t""yes" 1>&6
  ISC=yes # If later tests want to check for ISC.
  cat >> confdefs.h <<\EOF
#define _POSIX_SOURCE 1
EOF

  if test "$GCC" = yes; then
    CC="$CC -posix"
  else
    CC="$CC -Xp"
  fi
else
  echo "$ac_t""no" 1>&6
  ISC=
fi

# Extract the first word of "ranlib", so it can be a program name with args.
set dummy ranlib; ac_word=$2
echo $ac_n "checking for $ac_word""... $ac_c" 1>&6
echo "configure:1032: checking for $ac_word" >&5
if eval "test \"`echo '$''{'ac_cv_prog_RANLIB'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  if test -n "$RANLIB"; then
  ac_cv_prog_RANLIB="$RANLIB" # Let the user override the test.
else
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS=":"
  ac_dummy="$PATH"
  for ac_dir in $ac_dummy; do
    test -z "$ac_dir" && ac_dir=.
    if test -f $ac_dir/$ac_word; then
      ac_cv_prog_RANLIB="ranlib"
      break
    fi
  done
  IFS="$ac_save_ifs"
  test -z "$ac_cv_prog_RANLIB" && ac_cv_prog_RANLIB=":"
fi
fi
RANLIB="$ac_cv_prog_RANLIB"
if test -n "$RANLIB"; then
  echo "$ac_t""$RANLIB" 1>&6
else
  echo "$ac_t""no" 1>&6
fi

for ac_hdr in memory.h
do
ac_safe=`echo "$ac_hdr" | sed 'y%./+-%__p_%'`
echo $ac_n "checking for $ac_hdr""... $ac_c" 1>&6
echo "configure:1063: checking for $ac_hdr" >&5
if eval "test \"`echo '$''{'ac_cv_header_$ac_safe'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  cat > conftest.$ac_ext <<EOF
#line 1068 "configure"
#include "confdefs.h"
#include <$ac_hdr>
EOF
ac_try="$ac_cpp conftest.$ac_ext >/dev/null 2>conftest.out"
{ (eval echo configure:1073: \"$ac_try\") 1>&5; (eval $ac_try) 2>&5; }
ac_err=`grep -v '^ *+' conftest.out | grep -v "^conftest.${ac_ext}\$"`
if test -z "$ac_err"; then
  rm -rf conftest*
  eval "ac_cv_header_$ac_safe=yes"
else
  echo "$ac_err" >&5
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  eval "ac_cv_header_$ac_safe=no"
fi
rm -f conftest*
fi
if eval "test \"`echo '$ac_cv_header_'$ac_safe`\" = yes"; then
  echo "$ac_t""yes" 1>&6
    ac_tr_hdr=HAVE_`echo $ac_hdr | sed 'y%abcdefghijklmnopqrstuvwxyz./-%ABCDEFGHIJKLMNOPQRSTUVWXYZ___%'`
  cat >> confdefs.h <<EOF
#define $ac_tr_hdr 1
EOF
 
else
  echo "$ac_t""no" 1>&6
fi
done

echo $ac_n "checking for ANSI C header files""... $ac_c" 1>&6
echo "configure:1100: checking for ANSI C header files" >&5
if eval "test \"`echo '$''{'ac_cv_header_stdc'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  cat > conftest.$ac_ext <<EOF
#line 1105 "configure"
#include "confdefs.h"
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
EOF
ac_try="$ac_cpp conftest.$ac_ext >/dev/null 2>conftest.out"
{ (eval echo configure:1113: \"$ac_try\") 1>&5; (eval $ac_try) 2>&5; }
ac_err=`grep -v '^ *+' conftest.out | grep -v "^conftest.${ac_ext}\$"`
if test -z "$ac_err"; then
  rm -rf conftest*
  ac_cv_header_stdc=yes
else
  echo "$ac_err" >&5
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  ac_cv_header_stdc=no
fi
rm -f conftest*

if test $ac_cv_header_stdc = yes; then
  # SunOS 4.x string.h does not declare mem*, contrary to ANSI.
cat > conftest.$ac_ext <<EOF
#line 1130 "configure"
#include "confdefs.h"
#include <string.h>
EOF
if (eval "$ac_cpp conftest.$ac_ext") 2>&5 |
  egrep "memchr" >/dev/null 2>&1; then
  :
else
  rm -rf conftest*
  ac_cv_header_stdc=no
fi
rm -f conftest*

fi

if test $ac_cv_header_stdc = yes; then
  # ISC 2.0.2 stdlib.h does not declare free, contrary to ANSI.
cat > conftest.$ac_ext <<EOF
#line 1148 "configure"
#include "confdefs.h"
#include <stdlib.h>
EOF
if (eval "$ac_cpp conftest.$ac_ext") 2>&5 |
  egrep "free" >/dev/null 2>&1; then
  :
else
  rm -rf conftest*
  ac_cv_header_stdc=no
fi
rm -f conftest*

fi

if test $ac_cv_header_stdc = yes; then
  # /bin/cc in Irix-4.0.5 gets non-ANSI ctype macros unless using -ansi.
if test "$cross_compiling" = yes; then
  :
else
  cat > conftest.$ac_ext <<EOF
#line 1169 "configure"
#include "confdefs.h"
#include <ctype.h>
#define ISLOWER(c) ('a' <= (c) && (c) <= 'z')
#define TOUPPER(c) (ISLOWER(c) ? 'A' + ((c) - 'a') : (c))
#define XOR(e, f) (((e) && !(f)) || (!(e) && (f)))
int main () { int i; for (i = 0; i < 256; i++)
if (XOR (islower (i), ISLOWER (i)) || toupper (i) != TOUPPER (i)) exit(2);
exit (0); }

EOF
if { (eval echo configure:1180: \"$ac_link\") 1>&5; (eval $ac_link) 2>&5; } && test -s conftest${ac_exeext} && (./conftest; exit) 2>/dev/null
then
  :
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -fr conftest*
  ac_cv_header_stdc=no
fi
rm -fr conftest*
fi

fi
fi

echo "$ac_t""$ac_cv_header_stdc" 1>&6
if test $ac_cv_header_stdc = yes; then
  cat >> confdefs.h <<\EOF
#define STDC_HEADERS 1
EOF

fi

ac_safe=`echo "complex.h" | sed 'y%./+-%__p_%'`
echo $ac_n "checking for complex.h""... $ac_c" 1>&6
echo "configure:1205: checking for complex.h" >&5
if eval "test \"`echo '$''{'ac_cv_header_$ac_safe'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  cat > conftest.$ac_ext <<EOF
#line 1210 "configure"
#include "confdefs.h"
#include <complex.h>
EOF
ac_try="$ac_cpp conftest.$ac_ext >/dev/null 2>conftest.out"
{ (eval echo configure:1215: \"$ac_try\") 1>&5; (eval $ac_try) 2>&5; }
ac_err=`grep -v '^ *+' conftest.out | grep -v "^conftest.${ac_ext}\$"`
if test -z "$ac_err"; then
  rm -rf conftest*
  eval "ac_cv_header_$ac_safe=yes"
else
  echo "$ac_err" >&5
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  eval "ac_cv_header_$ac_safe=no"
fi
rm -f conftest*
fi
if eval "test \"`echo '$ac_cv_header_'$ac_safe`\" = yes"; then
  echo "$ac_t""yes" 1>&6
  cat >> confdefs.h <<\EOF
#define HAVE_COMPLEX_H 1
EOF

else
  echo "$ac_t""no" 1>&6
fi

ac_safe=`echo "malloc.h" | sed 'y%./+-%__p_%'`
echo $ac_n "checking for malloc.h""... $ac_c" 1>&6
echo "configure:1241: checking for malloc.h" >&5
if eval "test \"`echo '$''{'ac_cv_header_$ac_safe'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  cat > conftest.$ac_ext <<EOF
#line 1246 "configure"
#include "confdefs.h"
#include <malloc.h>
EOF
ac_try="$ac_cpp conftest.$ac_ext >/dev/null 2>conftest.out"
{ (eval echo configure:1251: \"$ac_try\") 1>&5; (eval $ac_try) 2>&5; }
ac_err=`grep -v '^ *+' conftest.out | grep -v "^conftest.${ac_ext}\$"`
if test -z "$ac_err"; then
  rm -rf conftest*
  eval "ac_cv_header_$ac_safe=yes"
else
  echo "$ac_err" >&5
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  eval "ac_cv_header_$ac_safe=no"
fi
rm -f conftest*
fi
if eval "test \"`echo '$ac_cv_header_'$ac_safe`\" = yes"; then
  echo "$ac_t""yes" 1>&6
  cat >> confdefs.h <<\EOF
#define HAVE_MALLOC_H 1
EOF

else
  echo "$ac_t""no" 1>&6
fi

ac_safe=`echo "varargs.h" | sed 'y%./+-%__p_%'`
echo $ac_n "checking for varargs.h""... $ac_c" 1>&6
echo "configure:1277: checking for varargs.h" >&5
if eval "test \"`echo '$''{'ac_cv_header_$ac_safe'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  cat > conftest.$ac_ext <<EOF
#line 1282 "configure"
#include "confdefs.h"
#include <varargs.h>
EOF
ac_try="$ac_cpp conftest.$ac_ext >/dev/null 2>conftest.out"
{ (eval echo configure:1287: \"$ac_try\") 1>&5; (eval $ac_try) 2>&5; }
ac_err=`grep -v '^ *+' conftest.out | grep -v "^conftest.${ac_ext}\$"`
if test -z "$ac_err"; then
  rm -rf conftest*
  eval "ac_cv_header_$ac_safe=yes"
else
  echo "$ac_err" >&5
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  eval "ac_cv_header_$ac_safe=no"
fi
rm -f conftest*
fi
if eval "test \"`echo '$ac_cv_header_'$ac_safe`\" = yes"; then
  echo "$ac_t""yes" 1>&6
  cat >> confdefs.h <<\EOF
#define VARARGS 1
EOF

else
  echo "$ac_t""no" 1>&6
fi

cat >> confdefs.h <<\EOF
#define NOT_SEGMENTED 1
EOF

echo $ac_n "checking for size_t""... $ac_c" 1>&6
echo "configure:1316: checking for size_t" >&5
if eval "test \"`echo '$''{'ac_cv_type_size_t'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  cat > conftest.$ac_ext <<EOF
#line 1321 "configure"
#include "confdefs.h"
#include <sys/types.h>
#if STDC_HEADERS
#include <stdlib.h>
#include <stddef.h>
#endif
EOF
if (eval "$ac_cpp conftest.$ac_ext") 2>&5 |
  egrep "(^|[^a-zA-Z_0-9])size_t[^a-zA-Z_0-9]" >/dev/null 2>&1; then
  rm -rf conftest*
  ac_cv_type_size_t=yes
else
  rm -rf conftest*
  ac_cv_type_size_t=no
fi
rm -f conftest*

fi
echo "$ac_t""$ac_cv_type_size_t" 1>&6
if test $ac_cv_type_size_t = no; then
  cat >> confdefs.h <<\EOF
#define size_t unsigned
EOF

fi

echo $ac_n "checking for working const""... $ac_c" 1>&6
echo "configure:1349: checking for working const" >&5
if eval "test \"`echo '$''{'ac_cv_c_const'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  cat > conftest.$ac_ext <<EOF
#line 1354 "configure"
#include "confdefs.h"

int main() {

/* Ultrix mips cc rejects this.  */
typedef int charset[2]; const charset x;
/* SunOS 4.1.1 cc rejects this.  */
char const *const *ccp;
char **p;
/* NEC SVR4.0.2 mips cc rejects this.  */
struct point {int x, y;};
static struct point const zero = {0,0};
/* AIX XL C 1.02.0.0 rejects this.
   It does not let you subtract one const X* pointer from another in an arm
   of an if-expression whose if-part is not a constant expression */
const char *g = "string";
ccp = &g + (g ? g-g : 0);
/* HPUX 7.0 cc rejects these. */
++ccp;
p = (char**) ccp;
ccp = (char const *const *) p;
{ /* SCO 3.2v4 cc rejects this.  */
  char *t;
  char const *s = 0 ? (char *) 0 : (char const *) 0;

  *t++ = 0;
}
{ /* Someone thinks the Sun supposedly-ANSI compiler will reject this.  */
  int x[] = {25, 17};
  const int *foo = &x[0];
  ++foo;
}
{ /* Sun SC1.0 ANSI compiler rejects this -- but not the above. */
  typedef const int *iptr;
  iptr p = 0;
  ++p;
}
{ /* AIX XL C 1.02.0.0 rejects this saying
     "k.c", line 2.27: 1506-025 (S) Operand must be a modifiable lvalue. */
  struct s { int j; const int *ap[3]; };
  struct s *b; b->j = 5;
}
{ /* ULTRIX-32 V3.1 (Rev 9) vcc rejects this */
  const int foo = 10;
}

; return 0; }
EOF
if { (eval echo configure:1403: \"$ac_compile\") 1>&5; (eval $ac_compile) 2>&5; }; then
  rm -rf conftest*
  ac_cv_c_const=yes
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  ac_cv_c_const=no
fi
rm -f conftest*
fi

echo "$ac_t""$ac_cv_c_const" 1>&6
if test $ac_cv_c_const = no; then
  cat >> confdefs.h <<\EOF
#define const 
EOF

fi

echo $ac_n "checking whether byte ordering is bigendian""... $ac_c" 1>&6
echo "configure:1424: checking whether byte ordering is bigendian" >&5
if eval "test \"`echo '$''{'ac_cv_c_bigendian'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  ac_cv_c_bigendian=unknown
# See if sys/param.h defines the BYTE_ORDER macro.
cat > conftest.$ac_ext <<EOF
#line 1431 "configure"
#include "confdefs.h"
#include <sys/types.h>
#include <sys/param.h>
int main() {

#if !BYTE_ORDER || !BIG_ENDIAN || !LITTLE_ENDIAN
 bogus endian macros
#endif
; return 0; }
EOF
if { (eval echo configure:1442: \"$ac_compile\") 1>&5; (eval $ac_compile) 2>&5; }; then
  rm -rf conftest*
  # It does; now see whether it defined to BIG_ENDIAN or not.
cat > conftest.$ac_ext <<EOF
#line 1446 "configure"
#include "confdefs.h"
#include <sys/types.h>
#include <sys/param.h>
int main() {

#if BYTE_ORDER != BIG_ENDIAN
 not big endian
#endif
; return 0; }
EOF
if { (eval echo configure:1457: \"$ac_compile\") 1>&5; (eval $ac_compile) 2>&5; }; then
  rm -rf conftest*
  ac_cv_c_bigendian=yes
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  ac_cv_c_bigendian=no
fi
rm -f conftest*
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
fi
rm -f conftest*
if test $ac_cv_c_bigendian = unknown; then
if test "$cross_compiling" = yes; then
    { echo "configure: error: can not run test program while cross compiling" 1>&2; exit 1; }
else
  cat > conftest.$ac_ext <<EOF
#line 1477 "configure"
#include "confdefs.h"
main () {
  /* Are we little or big endian?  From Harbison&Steele.  */
  union
  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  exit (u.c[sizeof (long) - 1] == 1);
}
EOF
if { (eval echo configure:1490: \"$ac_link\") 1>&5; (eval $ac_link) 2>&5; } && test -s conftest${ac_exeext} && (./conftest; exit) 2>/dev/null
then
  ac_cv_c_bigendian=no
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -fr conftest*
  ac_cv_c_bigendian=yes
fi
rm -fr conftest*
fi

fi
fi

echo "$ac_t""$ac_cv_c_bigendian" 1>&6
if test $ac_cv_c_bigendian = yes; then
  cat >> confdefs.h <<\EOF
#define WORDS_BIGENDIAN 1
EOF

fi

# Check whether --with-complex or --without-complex was given.
if test "${with_complex+set}" = set; then
  withval="$with_complex"
  :
fi

# Check whether --with-sparse or --without-sparse was given.
if test "${with_sparse+set}" = set; then
  withval="$with_sparse"
  :
fi

# Check whether --with-all or --without-all was given.
if test "${with_all+set}" = set; then
  withval="$with_all"
  :
fi

# Check whether --with-all or --without-all was given.
if test "${with_all+set}" = set; then
  withval="$with_all"
  :
fi

# Check whether --with-unroll or --without-unroll was given.
if test "${with_unroll+set}" = set; then
  withval="$with_unroll"
  :
fi

# Check whether --with-munroll or --without-munroll was given.
if test "${with_munroll+set}" = set; then
  withval="$with_munroll"
  :
fi

# Check whether --with-segmem or --without-segmem was given.
if test "${with_segmem+set}" = set; then
  withval="$with_segmem"
  :
fi

# Check whether --with-float or --without-float was given.
if test "${with_float+set}" = set; then
  withval="$with_float"
  :
fi

# Check whether --with-double or --without-double was given.
if test "${with_double+set}" = set; then
  withval="$with_double"
  :
fi

LIBS="$LIBS -lm"
cat > conftest.$ac_ext <<EOF
#line 1569 "configure"
#include "confdefs.h"
u_int
int main() {
#include <stdio.h>
#ifdef __STDC__
#include <stdlib.h>
#endif
; return 0; }
EOF
if { (eval echo configure:1579: \"$ac_compile\") 1>&5; (eval $ac_compile) 2>&5; }; then
  rm -rf conftest*
  u_int i; i = 1;
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  cat >> confdefs.h <<\EOF
#define U_INT_DEF 1
EOF

fi
rm -f conftest*
echo 'computing machine epsilon(s)'
echo $CC -o macheps dmacheps.c
$CC -o macheps dmacheps.c
cat >> confdefs.h <<EOF
#define D_MACHEPS `macheps`
EOF

echo $CC -o macheps fmacheps.c
$CC -o macheps fmacheps.c
cat >> confdefs.h <<EOF
#define F_MACHEPS `macheps`
EOF

echo computing M_MAX_INT
echo $CC -o maxint maxint.c
$CC -o maxint maxint.c
cat >> confdefs.h <<EOF
#define M_MAX_INT `maxint`
EOF

echo checking char '\\0' vs. float zeros
cat > conftest.$ac_ext <<EOF
#line 1614 "configure"
#include "confdefs.h"
main() {
    char	*cp = "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";
    double	*dp;
    dp = (double *)cp;
    if ( *dp == 0.0 ) printf("yes\n");  }

EOF
if (eval "$ac_cpp conftest.$ac_ext") 2>&5 |
  egrep "yes" >/dev/null 2>&1; then
  rm -rf conftest*
  cat >> confdefs.h <<\EOF
#define CHAR0ISDBL0 1
EOF

fi
rm -f conftest*

for ac_func in bcopy bzero
do
echo $ac_n "checking for $ac_func""... $ac_c" 1>&6
echo "configure:1636: checking for $ac_func" >&5
if eval "test \"`echo '$''{'ac_cv_func_$ac_func'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&6
else
  cat > conftest.$ac_ext <<EOF
#line 1641 "configure"
#include "confdefs.h"
/* System header to define __stub macros and hopefully few prototypes,
    which can conflict with char $ac_func(); below.  */
#include <assert.h>
/* Override any gcc2 internal prototype to avoid an error.  */
/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
char $ac_func();

int main() {

/* The GNU C library defines this for functions which it implements
    to always fail with ENOSYS.  Some functions are actually named
    something starting with __ and the normal name is an alias.  */
#if defined (__stub_$ac_func) || defined (__stub___$ac_func)
choke me
#else
$ac_func();
#endif

; return 0; }
EOF
if { (eval echo configure:1664: \"$ac_link\") 1>&5; (eval $ac_link) 2>&5; } && test -s conftest${ac_exeext}; then
  rm -rf conftest*
  eval "ac_cv_func_$ac_func=yes"
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  eval "ac_cv_func_$ac_func=no"
fi
rm -f conftest*
fi

if eval "test \"`echo '$ac_cv_func_'$ac_func`\" = yes"; then
  echo "$ac_t""yes" 1>&6
    ac_tr_func=HAVE_`echo $ac_func | tr 'abcdefghijklmnopqrstuvwxyz' 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'`
  cat >> confdefs.h <<EOF
#define $ac_tr_func 1
EOF
 
else
  echo "$ac_t""no" 1>&6
fi
done

cat > conftest.$ac_ext <<EOF
#line 1689 "configure"
#include "confdefs.h"
function prototypes
int main() {

; return 0; }
EOF
if { (eval echo configure:1696: \"$ac_compile\") 1>&5; (eval $ac_compile) 2>&5; }; then
  rm -rf conftest*
  extern int test (int i, double x);
else
  echo "configure: failed program was:" >&5
  cat conftest.$ac_ext >&5
  rm -rf conftest*
  cat >> confdefs.h <<\EOF
#define HAVE_PROTOTYPES 1
EOF

fi
rm -f conftest*
trap '' 1 2 15
cat > confcache <<\EOF
# This file is a shell script that caches the results of configure
# tests run on this system so they can be shared between configure
# scripts and configure runs.  It is not useful on other systems.
# If it contains results you don't want to keep, you may remove or edit it.
#
# By default, configure uses ./config.cache as the cache file,
# creating it if it does not exist already.  You can give configure
# the --cache-file=FILE option to use a different cache file; that is
# what configure does when it calls configure scripts in
# subdirectories, so they share the cache.
# Giving --cache-file=/dev/null disables caching, for debugging configure.
# config.status only pays attention to the cache file if you give it the
# --recheck option to rerun configure.
#
EOF
# The following way of writing the cache mishandles newlines in values,
# but we know of no workaround that is simple, portable, and efficient.
# So, don't put newlines in cache variables' values.
# Ultrix sh set writes to stderr and can't be redirected directly,
# and sets the high bit in the cache file unless we assign to the vars.
(set) 2>&1 |
  case `(ac_space=' '; set | grep ac_space) 2>&1` in
  *ac_space=\ *)
    # `set' does not quote correctly, so add quotes (double-quote substitution
    # turns \\\\ into \\, and sed turns \\ into \).
    sed -n \
      -e "s/'/'\\\\''/g" \
      -e "s/^\\([a-zA-Z0-9_]*_cv_[a-zA-Z0-9_]*\\)=\\(.*\\)/\\1=\${\\1='\\2'}/p"
    ;;
  *)
    # `set' quotes correctly as required by POSIX, so do not add quotes.
    sed -n -e 's/^\([a-zA-Z0-9_]*_cv_[a-zA-Z0-9_]*\)=\(.*\)/\1=${\1=\2}/p'
    ;;
  esac >> confcache
if cmp -s $cache_file confcache; then
  :
else
  if test -w $cache_file; then
    echo "updating cache $cache_file"
    cat confcache > $cache_file
  else
    echo "not updating unwritable cache $cache_file"
  fi
fi
rm -f confcache

trap 'rm -fr conftest* confdefs* core core.* *.core $ac_clean_files; exit 1' 1 2 15

test "x$prefix" = xNONE && prefix=$ac_default_prefix
# Let make expand exec_prefix.
test "x$exec_prefix" = xNONE && exec_prefix='${prefix}'

# Any assignment to VPATH causes Sun make to only execute
# the first set of double-colon rules, so remove it if not needed.
# If there is a colon in the path, we need to keep it.
if test "x$srcdir" = x.; then
  ac_vpsub='/^[ 	]*VPATH[ 	]*=[^:]*$/d'
fi

trap 'rm -f $CONFIG_STATUS conftest*; exit 1' 1 2 15

DEFS=-DHAVE_CONFIG_H

# Without the "./", some shells look in PATH for config.status.
: ${CONFIG_STATUS=./config.status}

echo creating $CONFIG_STATUS
rm -f $CONFIG_STATUS
cat > $CONFIG_STATUS <<EOF
#! /bin/sh
# Generated automatically by configure.
# Run this file to recreate the current configuration.
# This directory was configured as follows,
# on host `(hostname || uname -n) 2>/dev/null | sed 1q`:
#
# $0 $ac_configure_args
#
# Compiler output produced by configure, useful for debugging
# configure, is in ./config.log if it exists.

ac_cs_usage="Usage: $CONFIG_STATUS [--recheck] [--version] [--help]"
for ac_option
do
  case "\$ac_option" in
  -recheck | --recheck | --rechec | --reche | --rech | --rec | --re | --r)
    echo "running \${CONFIG_SHELL-/bin/sh} $0 $ac_configure_args --no-create --no-recursion"
    exec \${CONFIG_SHELL-/bin/sh} $0 $ac_configure_args --no-create --no-recursion ;;
  -version | --version | --versio | --versi | --vers | --ver | --ve | --v)
    echo "$CONFIG_STATUS generated by autoconf version 2.13"
    exit 0 ;;
  -help | --help | --hel | --he | --h)
    echo "\$ac_cs_usage"; exit 0 ;;
  *) echo "\$ac_cs_usage"; exit 1 ;;
  esac
done

ac_given_srcdir=$srcdir

trap 'rm -fr `echo "makefile machine.h" | sed "s/:[^ ]*//g"` conftest*; exit 1' 1 2 15
EOF
cat >> $CONFIG_STATUS <<EOF

# Protect against being on the right side of a sed subst in config.status.
sed 's/%@/@@/; s/@%/@@/; s/%g\$/@g/; /@g\$/s/[\\\\&%]/\\\\&/g;
 s/@@/%@/; s/@@/@%/; s/@g\$/%g/' > conftest.subs <<\\CEOF
$ac_vpsub
$extrasub
s%@SHELL@%$SHELL%g
s%@CFLAGS@%$CFLAGS%g
s%@CPPFLAGS@%$CPPFLAGS%g
s%@CXXFLAGS@%$CXXFLAGS%g
s%@FFLAGS@%$FFLAGS%g
s%@DEFS@%$DEFS%g
s%@LDFLAGS@%$LDFLAGS%g
s%@LIBS@%$LIBS%g
s%@exec_prefix@%$exec_prefix%g
s%@prefix@%$prefix%g
s%@program_transform_name@%$program_transform_name%g
s%@bindir@%$bindir%g
s%@sbindir@%$sbindir%g
s%@libexecdir@%$libexecdir%g
s%@datadir@%$datadir%g
s%@sysconfdir@%$sysconfdir%g
s%@sharedstatedir@%$sharedstatedir%g
s%@localstatedir@%$localstatedir%g
s%@libdir@%$libdir%g
s%@includedir@%$includedir%g
s%@oldincludedir@%$oldincludedir%g
s%@infodir@%$infodir%g
s%@mandir@%$mandir%g
s%@PROGS@%$PROGS%g
s%@CC@%$CC%g
s%@CPP@%$CPP%g
s%@RANLIB@%$RANLIB%g

CEOF
EOF

cat >> $CONFIG_STATUS <<\EOF

# Split the substitutions into bite-sized pieces for seds with
# small command number limits, like on Digital OSF/1 and HP-UX.
ac_max_sed_cmds=90 # Maximum number of lines to put in a sed script.
ac_file=1 # Number of current file.
ac_beg=1 # First line for current file.
ac_end=$ac_max_sed_cmds # Line after last line for current file.
ac_more_lines=:
ac_sed_cmds=""
while $ac_more_lines; do
  if test $ac_beg -gt 1; then
    sed "1,${ac_beg}d; ${ac_end}q" conftest.subs > conftest.s$ac_file
  else
    sed "${ac_end}q" conftest.subs > conftest.s$ac_file
  fi
  if test ! -s conftest.s$ac_file; then
    ac_more_lines=false
    rm -f conftest.s$ac_file
  else
    if test -z "$ac_sed_cmds"; then
      ac_sed_cmds="sed -f conftest.s$ac_file"
    else
      ac_sed_cmds="$ac_sed_cmds | sed -f conftest.s$ac_file"
    fi
    ac_file=`expr $ac_file + 1`
    ac_beg=$ac_end
    ac_end=`expr $ac_end + $ac_max_sed_cmds`
  fi
done
if test -z "$ac_sed_cmds"; then
  ac_sed_cmds=cat
fi
EOF

cat >> $CONFIG_STATUS <<EOF

CONFIG_FILES=\${CONFIG_FILES-"makefile"}
EOF
cat >> $CONFIG_STATUS <<\EOF
for ac_file in .. $CONFIG_FILES; do if test "x$ac_file" != x..; then
  # Support "outfile[:infile[:infile...]]", defaulting infile="outfile.in".
  case "$ac_file" in
  *:*) ac_file_in=`echo "$ac_file"|sed 's%[^:]*:%%'`
       ac_file=`echo "$ac_file"|sed 's%:.*%%'` ;;
  *) ac_file_in="${ac_file}.in" ;;
  esac

  # Adjust a relative srcdir, top_srcdir, and INSTALL for subdirectories.

  # Remove last slash and all that follows it.  Not all systems have dirname.
  ac_dir=`echo $ac_file|sed 's%/[^/][^/]*$%%'`
  if test "$ac_dir" != "$ac_file" && test "$ac_dir" != .; then
    # The file is in a subdirectory.
    test ! -d "$ac_dir" && mkdir "$ac_dir"
    ac_dir_suffix="/`echo $ac_dir|sed 's%^\./%%'`"
    # A "../" for each directory in $ac_dir_suffix.
    ac_dots=`echo $ac_dir_suffix|sed 's%/[^/]*%../%g'`
  else
    ac_dir_suffix= ac_dots=
  fi

  case "$ac_given_srcdir" in
  .)  srcdir=.
      if test -z "$ac_dots"; then top_srcdir=.
      else top_srcdir=`echo $ac_dots|sed 's%/$%%'`; fi ;;
  /*) srcdir="$ac_given_srcdir$ac_dir_suffix"; top_srcdir="$ac_given_srcdir" ;;
  *) # Relative path.
    srcdir="$ac_dots$ac_given_srcdir$ac_dir_suffix"
    top_srcdir="$ac_dots$ac_given_srcdir" ;;
  esac


  echo creating "$ac_file"
  rm -f "$ac_file"
  configure_input="Generated automatically from `echo $ac_file_in|sed 's%.*/%%'` by configure."
  case "$ac_file" in
  *Makefile*) ac_comsub="1i\\
# $configure_input" ;;
  *) ac_comsub= ;;
  esac

  ac_file_inputs=`echo $ac_file_in|sed -e "s%^%$ac_given_srcdir/%" -e "s%:% $ac_given_srcdir/%g"`
  sed -e "$ac_comsub
s%@configure_input@%$configure_input%g
s%@srcdir@%$srcdir%g
s%@top_srcdir@%$top_srcdir%g
" $ac_file_inputs | (eval "$ac_sed_cmds") > $ac_file
fi; done
rm -f conftest.s*

# These sed commands are passed to sed as "A NAME B NAME C VALUE D", where
# NAME is the cpp macro being defined and VALUE is the value it is being given.
#
# ac_d sets the value in "#define NAME VALUE" lines.
ac_dA='s%^\([ 	]*\)#\([ 	]*define[ 	][ 	]*\)'
ac_dB='\([ 	][ 	]*\)[^ 	]*%\1#\2'
ac_dC='\3'
ac_dD='%g'
# ac_u turns "#undef NAME" with trailing blanks into "#define NAME VALUE".
ac_uA='s%^\([ 	]*\)#\([ 	]*\)undef\([ 	][ 	]*\)'
ac_uB='\([ 	]\)%\1#\2define\3'
ac_uC=' '
ac_uD='\4%g'
# ac_e turns "#undef NAME" without trailing blanks into "#define NAME VALUE".
ac_eA='s%^\([ 	]*\)#\([ 	]*\)undef\([ 	][ 	]*\)'
ac_eB='$%\1#\2define\3'
ac_eC=' '
ac_eD='%g'

if test "${CONFIG_HEADERS+set}" != set; then
EOF
cat >> $CONFIG_STATUS <<EOF
  CONFIG_HEADERS="machine.h"
EOF
cat >> $CONFIG_STATUS <<\EOF
fi
for ac_file in .. $CONFIG_HEADERS; do if test "x$ac_file" != x..; then
  # Support "outfile[:infile[:infile...]]", defaulting infile="outfile.in".
  case "$ac_file" in
  *:*) ac_file_in=`echo "$ac_file"|sed 's%[^:]*:%%'`
       ac_file=`echo "$ac_file"|sed 's%:.*%%'` ;;
  *) ac_file_in="${ac_file}.in" ;;
  esac

  echo creating $ac_file

  rm -f conftest.frag conftest.in conftest.out
  ac_file_inputs=`echo $ac_file_in|sed -e "s%^%$ac_given_srcdir/%" -e "s%:% $ac_given_srcdir/%g"`
  cat $ac_file_inputs > conftest.in

EOF

# Transform confdefs.h into a sed script conftest.vals that substitutes
# the proper values into config.h.in to produce config.h.  And first:
# Protect against being on the right side of a sed subst in config.status.
# Protect against being in an unquoted here document in config.status.
rm -f conftest.vals
cat > conftest.hdr <<\EOF
s/[\\&%]/\\&/g
s%[\\$`]%\\&%g
s%#define \([A-Za-z_][A-Za-z0-9_]*\) *\(.*\)%${ac_dA}\1${ac_dB}\1${ac_dC}\2${ac_dD}%gp
s%ac_d%ac_u%gp
s%ac_u%ac_e%gp
EOF
sed -n -f conftest.hdr confdefs.h > conftest.vals
rm -f conftest.hdr

# This sed command replaces #undef with comments.  This is necessary, for
# example, in the case of _POSIX_SOURCE, which is predefined and required
# on some systems where configure will not decide to define it.
cat >> conftest.vals <<\EOF
s%^[ 	]*#[ 	]*undef[ 	][ 	]*[a-zA-Z_][a-zA-Z_0-9]*%/* & */%
EOF

# Break up conftest.vals because some shells have a limit on
# the size of here documents, and old seds have small limits too.

rm -f conftest.tail
while :
do
  ac_lines=`grep -c . conftest.vals`
  # grep -c gives empty output for an empty file on some AIX systems.
  if test -z "$ac_lines" || test "$ac_lines" -eq 0; then break; fi
  # Write a limited-size here document to conftest.frag.
  echo '  cat > conftest.frag <<CEOF' >> $CONFIG_STATUS
  sed ${ac_max_here_lines}q conftest.vals >> $CONFIG_STATUS
  echo 'CEOF
  sed -f conftest.frag conftest.in > conftest.out
  rm -f conftest.in
  mv conftest.out conftest.in
' >> $CONFIG_STATUS
  sed 1,${ac_max_here_lines}d conftest.vals > conftest.tail
  rm -f conftest.vals
  mv conftest.tail conftest.vals
done
rm -f conftest.vals

cat >> $CONFIG_STATUS <<\EOF
  rm -f conftest.frag conftest.h
  echo "/* $ac_file.  Generated automatically by configure.  */" > conftest.h
  cat conftest.in >> conftest.h
  rm -f conftest.in
  if cmp -s $ac_file conftest.h 2>/dev/null; then
    echo "$ac_file is unchanged"
    rm -f conftest.h
  else
    # Remove last slash and all that follows it.  Not all systems have dirname.
      ac_dir=`echo $ac_file|sed 's%/[^/][^/]*$%%'`
      if test "$ac_dir" != "$ac_file" && test "$ac_dir" != .; then
      # The file is in a subdirectory.
      test ! -d "$ac_dir" && mkdir "$ac_dir"
    fi
    rm -f $ac_file
    mv conftest.h $ac_file
  fi
fi; done

EOF
cat >> $CONFIG_STATUS <<EOF

EOF
cat >> $CONFIG_STATUS <<\EOF

exit 0
EOF
chmod +x $CONFIG_STATUS
rm -fr confdefs* $ac_clean_files
test "$no_create" = yes || ${CONFIG_SHELL-/bin/sh} $CONFIG_STATUS || exit 1

echo "Extensions to basic version: use configure --with-opt1 --with-opt2"
echo "  Option:"
echo "    --with-complex     incorporate complex functions"
echo "    --with-sparse      incorporate sparse matrix functions"
echo "    --with-all         both of the above"
echo "    --with-unroll      unroll low level loops on vectors"
echo "    --with-munroll     unroll low level loops on matrices"
echo "    --with-float       single precision"
echo "    --with-double      double precision (default)"
echo "Re-run configure with these options if you want them"
# configure.in copyright (C) Brook Milligan and David Stewart, 1993, 2002
