#AM_CFLAGS = -I./include -I../zlib
#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
include_HEADERS=include/rw.h include/ByteToolkit.h include/DynamicFloatArray.h include/DynamicDoubleArray.h include/qca_gnuplot.h

lib_LTLIBRARIES=libqca.la
libqca_la_CFLAGS=-I./include
libqca_la_LDFLAGS = -version-info  0:1:0
libqca_la_LIDADD=
libqca_la_SOURCES=src/rw.c src/ByteToolkit.c src/DynamicFloatArray.c src/DynamicDoubleArray.c src/qca_gnuplot.c

libqca_la_LINK= $(AM_V_CC)$(LIBTOOL) --tag=CC --mode=link $(CCLD) $(libqca_la_CFLAGS) -O3 $(libqca_la_LDFLAGS) -o $(lib_LTLIBRARIES)
