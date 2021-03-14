
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
} QCAT_CompressionResult;

QCAT_CompressionResult* compareData(int dataType, size_t nbEle, void* data, void* dec);
QCAT_CompressionResult* getCompressionResult(int dataType, float errBound, int quantBinCapacity, void* origData, void* predData, QCAT_DataProperty* property);
void printCompressionResult(QCAT_CompressionResult* result);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _QCAT_COMPRESSIONANALYSIS_H  ----- */
