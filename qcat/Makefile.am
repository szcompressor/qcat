#AM_CFLAGS = -I./include -I../zlib
#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
include_HEADERS=include/qcat_hashtable.h include/DynamicIntArray.h include/DynamicByteArray.h include/Huffman.h include/rw.h include/ByteToolkit.h include/DynamicFloatArray.h include/DynamicDoubleArray.h include/qcat_gnuplot.h include/qcat_dataAnalysis.h include/qcat_compressionAnalysis.h include/sz_dummy_compression.h include/sz_utility.h include/qcat.h include/qcat_ssim.h include/der.hpp include/idx.hpp include/matrix.hpp  include/m_basic.hpp  include/m_der.hpp  include/m_der_metrics.hpp include/m_norm.hpp  include/norm.hpp  include/others.hpp  include/ssim.hpp

lib_LTLIBRARIES=libqcat.la
libqcat_la_CPPFLAGS=-I./include -I../zstd/
libqcat_la_LDFLAGS = -version-info  0:1:0
libqcat_la_LIDADD=../zstd/.libs/libzstd.a
libqcat_la_SOURCES=src/qcat_hashtable.c src/DynamicIntArray.c src/DynamicByteArray.c src/Huffman.c src/rw.c src/ByteToolkit.c src/DynamicFloatArray.c src/DynamicDoubleArray.c src/qcat_gnuplot.c src/qcat_dataAnalysis.c src/qcat_compressionAnalysis.c src/sz_dummy_compression.c src/sz_utility.c src/qcat_ssim.c src/qcat_derivativeMetrics.cpp

libqcat_la_LINK= $(AM_V_CC)$(LIBTOOL) --tag=CC --mode=link $(CCLD) $(libqcat_la_CPPFLAGS) -O3 $(libqcat_la_LDFLAGS) -o $(lib_LTLIBRARIES)
