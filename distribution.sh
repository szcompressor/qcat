#!/bin/bash

ver=1.6

rm -rf qcat-${ver}.tar.gz
make dist
tar -xzvf qcat-${ver}.tar.gz
cp README.md qcat-${ver}
cp CMakeLists.txt qcat-${ver}
cp qcat/CMakeLists.txt qcat-${ver}/qcat

tar -czvf qcat-${ver}.tar.gz qcat-${ver}
