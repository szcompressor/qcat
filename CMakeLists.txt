cmake_minimum_required (VERSION 3.1.3)
project (qcat)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -lm")

set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

include_directories (
	${CMAKE_CURRENT_SOURCE_DIR}/qcat/include 
	)

add_subdirectory (qcat)
add_subdirectory (examples)
