#AM_CFLAGS = -I./include -I../zlib
#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
include_HEADERS=include/DynamicIntArray.h include/DynamicByteArray.h include/Huffman.h include/rw.h include/ByteToolkit.h include/DynamicFloatArray.h include/DynamicDoubleArray.h include/qcat_gnuplot.h include/qcat_dataAnalysis.h include/qcat_compressionAnalysis.h include/sz_dummy_compression.h include/sz_utility.h include/qcat.h

lib_LTLIBRARIES=libqcat.la
libqcat_la_CFLAGS=-I./include
libqcat_la_LDFLAGS = -version-info  0:1:0
libqcat_la_LIDADD=
libqcat_la_SOURCES=src/DynamicIntArray.c src/DynamicByteArray.c src/Huffman.c src/rw.c src/ByteToolkit.c src/DynamicFloatArray.c src/DynamicDoubleArray.c src/qcat_gnuplot.c src/qcat_dataAnalysis.c src/qcat_compressionAnalysis.c src/sz_dummy_compression.c src/sz_utility.c

libqcat_la_LINK= $(AM_V_CC)$(LIBTOOL) --tag=CC --mode=link $(CCLD) $(libqcat_la_CFLAGS) -O3 $(libqcat_la_LDFLAGS) -o $(lib_LTLIBRARIES)
