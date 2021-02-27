
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

#ifndef _QCAT_DATAANALYSIS_H
#define _QCAT_DATAANALYSIS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QCAT_DataProperty
{
	int dataType; /*DA_DOUBLE or DA_FLOAT*/
	size_t r5;
	size_t r4;
	size_t r3;
	size_t r2;
	size_t r1;
	
	long numOfElem;
	double minValue;
	double maxValue;
	double valueRange;
	double avgValue;
	double entropy;
	double zeromean_variance;
	size_t totalByteSize;
} QCAT_DataProperty;


double computeEntropy(int dataType, void* data, size_t nbEle);
QCAT_DataProperty* computeProperty(int dataType, void* data, size_t nbEle);
void printProperty(QCAT_DataProperty* property);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _QCAT_DATAANALYSIS_H  ----- */
