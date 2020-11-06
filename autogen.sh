#!/bin/bash

aclocal
libtoolize -f -c
autoconf
autoheader
automake -a

#make dist

#cp -rf example/README example/Makefile.bk sz-1.0/example/
#cp -rf doc sz-1.1/

#After generating sz/Makefile, then modify sz/Makefile as follows, by adding "-rm -f *.mod"
#mostlyclean-compile:
#        -rm -f *.$(OBJEXT)
#        -rm -f *.mod
