cmake_minimum_required (VERSION 3.12)
project (qcat VERSION 1.5.0)
include(GNUInstallDirs)
find_package(PkgConfig)

option(QCAT_USE_BUNDLES "prefer bundled versions of sofware" ON)
option(BUILD_SHARED_LIBS "build shared libraies instead of static ones" ON)
if(NOT QCAT_USE_BUNDLES)
  pkg_search_module(ZSTD IMPORTED_TARGET GLOBAL libzstd)
endif()

if(NOT ZSTD_FOUND)
  add_subdirectory(zstd)
endif()

add_subdirectory (qcat)
add_subdirectory (examples)
