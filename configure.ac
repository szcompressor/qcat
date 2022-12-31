#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(QCAT, 1.5, sdi1@anl.gov)
AM_INIT_AUTOMAKE([foreign -Wall -Werror subdir-objects])
AC_CONFIG_SRCDIR([qcat/src/rw.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
: ${CFLAGS=-O3 -std=gnu99 -Wall}

AC_LANG(C++)
AC_PROG_CXX
AC_PROG_CC

# Checks for programs.
#: ${CFLAGS=-O3 -std=c99 -Wall}
#: ${CFLAGS=-O3 -std=gnu99 -Wall}
: ${CFLAGS=-O3 -std=gnu99 -Wall}

AC_PROG_RANLIB
AM_PROG_AR
AC_PROG_LIBTOOL
#m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AM_MAINTAINER_MODE([disable])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h limits.h malloc.h stddef.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memchr memset strchr strtol])

# Setup CUDA paths
#default: 64bit
AC_ARG_ENABLE([32bit],
	AC_HELP_STRING([--enable-32bit],[force compilation in 32 bit mode.]),
	[ARCH="386"; CUDA_LIB_DIR="lib"; PLATFORM="-m32";
		AC_MSG_NOTICE([Enabling 32bit compilation])
	],
	[ARCH="x86_64"; CUDA_LIB_DIR="lib64"; PLATFORM="-m64";
		AC_MSG_NOTICE([Enabling 64bit compilation])
	]
)

AC_ARG_WITH([cuda],
   [AC_HELP_STRING([--with-cuda=PATH],[path where cuda is installed (default not using cuda])],
   [CUDA_HOME=$withval],
   [CUDA_HOME=notusing])

#AC_ARG_WITH([dummy],
#	[AC_HELP_STRING([--without-cuda],[disable CUDA compilation])]
#	)
#AC_MSG_RESULT($CUDA_HOME)

AM_CONDITIONAL([COMPILE_CUDA], [test "x$CUDA_HOME" != "xnotusing"])

AS_IF([test "x$CUDA_HOME" != "xnotusing"], [

	CUDA_CFLAGS="-I$CUDA_HOME/include"
	CUDA_LIBS="-L$CUDA_HOME/$CUDA_LIB_DIR"
	CUDA_LDFLAGS="-L$CUDA_HOME/$CUDA_LIB_DIR"

	AC_CHECK_FILE($CUDA_HOME/$CUDA_LIB_DIR/libcudart.so, [],
		[AC_MSG_ERROR([libcudart.so library not found in path: $CUDA_HOME/$CUDA_LIB_DIR. Try running 
	                  configure again with the './configure --with-cuda=PATH' argument.])]
	)
	AC_ARG_WITH([nvcc],
	   [AC_HELP_STRING([--with-nvcc=PATH],[path where nvcc is installed (e.g. /usr/local/cuda/bin)])],
	   [NVCC_PATH=$withval],
	   [NVCC_PATH=$CUDA_HOME/bin]
	)
	AC_PATH_PROG([NVCC],[nvcc],[],[$NVCC_PATH])
	if test -z $NVCC; then
		AC_MSG_ERROR([CUDA nvcc compiler not found in path: $NVCC_PATH. Try running 
	                  configure again with the './configure --with-nvcc=PATH' argument.])
	fi
	
	#NVCCFLAGS="--compiler-bindir $NVCC_COMPILER_BINDIR $CUDA_CFLAGS"
	
	#
	# Check for CUDA libraries
	#
	save_LDFLAGS="$LDFLAGS"
	if [[ $ARCH == "386" ]]; then
		LDFLAGS="$LDFLAGS $CUDA_LDFLAGS -m32"
		echo $LDFLAGS
		
		AC_CHECK_LIB([z], [deflate],[],[AC_MSG_ERROR([
	                  libz.so:i386 error. Install the lib32z1-dev package or
	                  try running configure again with the './configure --enable-64bit' argument.
	                  You can also supply a ./configure LDFLAGS='-L<DIR>', where <DIR>
	                  is the directory where the 32bit version of libz.so resides.])]
		)
		
		AC_CHECK_LIB([stdc++], [memcpy],[],[AC_MSG_ERROR([
	                  libstdc++.so:i386 error. Install the g++-multilib package or
	                  try running configure again with the './configure --enable-64bit' argument.
	                  You can also supply a ./configure LDFLAGS='-L<DIR>', where <DIR>
	                  is the directory where the 32bit version of libstdc++.so resides.])]
		)
		
		AC_CHECK_LIB([cuda], [cuMemcpy],[],[AC_MSG_ERROR([
	                  libcuda.so:i386 error. Install the CUDA 32bit package or
	                  try running configure again with the './configure --enable-64bit' argument.
	                  You can also supply a ./configure LDFLAGS='-L<DIR>', where <DIR>
	                  is the directory where the 32bit version of libcuda.so resides.])]
		)
	else
		LDFLAGS="$LDFLAGS $CUDA_LDFLAGS"
		
		AC_CHECK_LIB([cuda], [cuMemcpy],[],[AC_MSG_ERROR([
	                  libcuda.so (cuMemcpy) was not found.])
		])
	fi
	LDFLAGS="$save_LDFLAGS"

	save_CPPFLAGS="$LDFLAGS"
	CPPFLAGS+=$CUDA_CFLAGS
	# Checks for header files.
	AC_CHECK_HEADERS([cuda_runtime.h],,
		AC_MSG_ERROR([Missing required header files.]))
	CPPFLAGS="$save_CPPFLAGS"


	AC_DEFINE([WITH_CUDA],[1],[Defined if CUDA should be used])
	NVCCFLAGS="$CUDA_CFLAGS $NVCCFLAGS $NVCC_COMPILER_BINDIR $NVCC_CUDA_ARCH $NVCC_CUDA_DEBUG $CUDA_REGISTERS"

	AC_SUBST(CUDA_CFLAGS)
	AC_SUBST(CUDA_LIBS)
	AC_SUBST(NVCC)
	AC_SUBST(NVCCFLAGS)
])


AC_CONFIG_FILES([Makefile
		 zstd/Makefile
                 qcat/Makefile
		 examples/Makefile])
AC_OUTPUT
