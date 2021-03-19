/**
 *  @file qcat_ssim.h
 *  @author Ali Gok, Sheng Di
 *  @date April, 2021
 *  @brief Header file for the whole io interface.
 *  (C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifndef _QCAT_SSIM_H
#define _QCAT_SSIM_H

#define K1 0.01
#define K2 0.03

#ifdef __cplusplus
extern "C" {
#endif

double SSIM_1d_calcWindow_float(float* data, float* other, int offset0, int windowSize0);
double SSIM_1d_calcWindow_double(double* data, double* other, int offset0, int windowSize0);
double SSIM_1d_windowed_float(float* oriData, float* decData, size_t size0, int windowSize0, int windowShift0);
double SSIM_1d_windowed_double(double* oriData, double* decData, size_t size0, int windowSize0, int windowShift0);

double SSIM_2d_calcWindow_float(float* data, float *other, size_t size0, int offset0, int offset1, int windowSize0, int windowSize1);
double SSIM_2d_calcWindow_double(double* data, double *other, size_t size0, int offset0, int offset1, int windowSize0, int windowSize1);
double SSIM_2d_windowed_float(float* oriData, float* decData, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowShift0, int windowShift1);double SSIM_2d_windowed_double(double* oriData, double* decData, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowShift0, int windowShift1);

double SSIM_3d_calcWindow_float(float* data, float* other, size_t size1, size_t size0, int offset0, int offset1, int offset2, int windowSize0, int windowSize1, int windowSize2);
double SSIM_3d_calcWindow_double(double* data, double* other, size_t size1, size_t size0, int offset0, int offset1, int offset2, int windowSize0, int windowSize1, int windowSize2);
double SSIM_3d_windowed_float(float* oriData, float* decData, size_t size2, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowSize2, int windowShift0, int windowShift1, int windowShift2);
double SSIM_3d_windowed_double(double* oriData, double* decData, size_t size2, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowSize2, int windowShift0, int windowShift1, int windowShift2);

double SSIM_4d_calcWindow_float(float* data, float* other, size_t size2, size_t size1, size_t size0, int offset0, int offset1, int offset2, int offset3,int windowSize0, int windowSize1, int windowSize2, int windowSize3);
double SSIM_4d_calcWindow_double(double* data, double* other, size_t size2, size_t size1, size_t size0, int offset0, int offset1, int offset2, int offset3,int windowSize0, int windowSize1, int windowSize2, int windowSize3);
double SSIM_4d_windowed_float(float* oriData, float* decData, size_t size3, size_t size2, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowSize2, int windowSize3, int windowShift0, int windowShift1, int windowShift2, int windowShift3);
double SSIM_4d_windowed_double(double* oriData, double* decData, size_t size3, size_t size2, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowSize2, int windowSize3, int windowShift0, int windowShift1, int windowShift2, int windowShift3);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _QCAT_SSIM_H  ----- */
