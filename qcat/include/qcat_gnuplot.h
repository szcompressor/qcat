/**
 *  @file QCAT_GNUPLOT.h
 *  @author Sheng Di
 *  @date July, 2019
 *  @brief Header file for the SDA_GNUPLOT.c.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#ifndef _QCAT_GNUPLOT_H
#define _QCAT_GNUPLOT_H

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

char** genGnuplotScript_sliceImage(char* dataFileName, size_t r2, size_t r1, char* imageFileName);
float* computeSlice2DLog_float(size_t r2, size_t r1, float* data);
float* computeSlice3DLog_float(int sliceNumber, size_t r3, size_t r2, size_t r1, float* data);
double* computeSlice2DLog_double(size_t r2, size_t r1, double* data);
double* computeSlice3DLog_double(int sliceNumber, size_t r3, size_t r2, size_t r1, double* data);
float* transformData_float(int plotDim, size_t r3, size_t r2, size_t r1, float* oriData);
double* transformData_double(int plotDim, size_t r3, size_t r2, size_t r1, double* oriData);
float* generateSliceData_float(int plotDim, int sliceNumber, size_t r3, size_t r2, size_t r1, float* data, int domain);
double* generateSliceData_double(int plotDim, int sliceNumber, size_t r3, size_t r2, size_t r1, double* data, int domain);
float* generateSliceDiff_float(int plotDim, int sliceNumber, size_t r3, size_t r2, size_t r1, float* oriData, float* decData, int domain);
double* generateSliceDiff_double(int plotDim, int sliceNumber, size_t r3, size_t r2, size_t r1, double* oriData, double* decData, int domain);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _QCAT_GNUPLOT_H  ----- */

