
/**
 *  @file qcat_dataAnalysis.h
 *  @author Sheng Di
 *  @date Feb, 2021
 *  @brief Header file.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdint.h>
#include <qcat_dataAnalysis.h>
#include <Huffman.h>
#include <zstd.h>

#ifndef _QCAT_COMPRESSIONANALYSIS_H
#define _QCAT_COMPRESSIONANALYSIS_H

#define PDF_INTERVALS 10000

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QCAT_CompressionResult
{
	QCAT_DataProperty* property;
	float compressionSize;
	float compressionRatio;
	float mse;
	float nrmse;
	float psnr;
	float normErr;
	float normErr_norm;
	float maxABSError;
	float maxRELError;
	float maxPWRError;
	float pearsonEff;
	size_t unpredictableCount;
	float unpredictablePercent;
	double ssim;
} QCAT_CompressionResult;

double* computePDF(int dataType, void* oriData, void* decData, size_t numOfElem, double* min_diff, double* err_interval);
double calculateSSIM(void* oriData, void* decData, int dataType, size_t r4, size_t r3, size_t r2, size_t r1);
QCAT_CompressionResult* compareData(int dataType, size_t nbEle, void* data, void* dec);
QCAT_CompressionResult* getCompressionResult(int dataType, float errBound, int quantBinCapacity, void* origData, void* predData, QCAT_DataProperty* property);
void printCompressionResult(QCAT_CompressionResult* result);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _QCAT_COMPRESSIONANALYSIS_H  ----- */
