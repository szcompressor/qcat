#AM_CFLAGS = -I./include -I../zlib
#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
include_HEADERS=include/rw.h include/ByteToolkit.h include/DynamicFloatArray.h include/DynamicDoubleArray.h include/qcat_gnuplot.h include/qcat_dataAnalysis.h

lib_LTLIBRARIES=libqcat.la
libqcat_la_CFLAGS=-I./include
libqcat_la_LDFLAGS = -version-info  0:1:0
libqcat_la_LIDADD=
libqcat_la_SOURCES=src/rw.c src/ByteToolkit.c src/DynamicFloatArray.c src/DynamicDoubleArray.c src/qcat_gnuplot.c src/qcat_dataAnalysis.c

libqcat_la_LINK= $(AM_V_CC)$(LIBTOOL) --tag=CC --mode=link $(CCLD) $(libqcat_la_CFLAGS) -O3 $(libqcat_la_LDFLAGS) -o $(lib_LTLIBRARIES)
