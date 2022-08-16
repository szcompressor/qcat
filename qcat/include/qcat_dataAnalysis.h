
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

#define QCAT_FLOAT 0
#define QCAT_DOUBLE 1
#define QCAT_INT32 2
#define QCAT_INT16 3
#define QCAT_UINT32 4
#define QCAT_UINT16 5

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
	double entropy_8bits;
	double entropy_32bits;
	double entropy_64bits;
	double zeromean_variance;
	size_t totalByteSize;
} QCAT_DataProperty;


typedef struct QCAT__ELEMENT
{
	double value;
	unsigned int counter;
} QCAT_ELEMENT;

double computeLosslessEntropy(int dataType, void* data, size_t nbEle);
double computeLosslessEntropy_32bits_fast(void* data, size_t nbEle);
double computeLosslessEntropy_32bits(void* data, size_t nbEle);
double computeLosslessEntropy_64bits(void* data, size_t nbEle);
QCAT_DataProperty* computeProperty(int dataType, void* data, size_t nbEle, int entropyType);
void printProperty(QCAT_DataProperty* property);
double* computeDataPDF_int32(void* data, size_t numOfElem, int* min, int* intervals);
double* computeDataPDF_float(float* data, size_t numOfElem, int intervals, float* min, double* unit, float mint, float maxt);

void computeLaplacian_float(float *data, float *lap, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
void computeLaplacian_double(double *data, double *lap, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
int computeLaplacian(void *data, void *lap, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _QCAT_DATAANALYSIS_H  ----- */
