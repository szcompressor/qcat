/**
 *  @file sz_utility.h
 *  @author Sheng Di
 *  @date March, 2021
 *  @brief Header file for sz_utility.c.
 *  (C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <DynamicByteArray.h>
#include <DynamicIntArray.h>
#include <rw.h>

#ifndef _SZ_UTILITY_H
#define _SZ_UTILITY_H

#ifdef __cplusplus
extern "C" {
#endif

#define DynArrayInitLen 1024

#define LORENZO_1D_1LAYER 0
#define LORENZO_1D_2LAYER 1
#define LORENZO_1D_3LAYER 2
#define LORENZO_2D_1LAYER 3
#define LORENZO_3D_1LAYER 4

#define QUANT_CODE_ORIGINAL 0
#define QUANT_CODE_NORMALIZE 1
#define QUANT_CODE_SHIFT 2

typedef struct DoubleValueCompressElement
{
	double data;
	long curValue;
	unsigned char curBytes[8]; //big_endian
	int reqBytesLength;
	int resiBitsLength;
} DoubleValueCompressElement;

typedef struct FloatValueCompressElement
{
	float data;
	int curValue;
	unsigned char curBytes[4]; //big_endian
	int reqBytesLength;
	int resiBitsLength;
} FloatValueCompressElement;

typedef struct LossyCompressionElement
{
	int leadingZeroBytes; //0,1,2,or 3
	unsigned char integerMidBytes[8];
	int integerMidBytes_Length; //they are mid_bits actually
	//char curBytes[8];
	//int curBytes_Length; //4 for single_precision or 8 for double_precision       
	int resMidBitsLength;
	int residualMidBits;
} LossyCompressionElement;


int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
void computeReqLength_float(double realPrecision, short radExpo, int* reqLength, float* medianValue);
void compressSingleFloatValue(FloatValueCompressElement *vce, float tgtValue, float precision, float medianValue, 
		int reqLength, int reqBytesLength, int resiBitsLength);
void updateLossyCompElement_Float(unsigned char* curBytes, unsigned char* preBytes, 
		int reqBytesLength, int resiBitsLength,  LossyCompressionElement *lce);
void addExactData(DynamicByteArray *exactMidByteArray, DynamicIntArray *exactLeadNumArray, 
		DynamicIntArray *resiBitArray, LossyCompressionElement *lce);
int compIdenticalLeadingBytesCount_double(unsigned char* preBytes, unsigned char* curBytes);
int compIdenticalLeadingBytesCount_float(unsigned char* preBytes, unsigned char* curBytes);
size_t convertIntArray2ByteArray_fast_2b(unsigned char* timeStepType, size_t timeStepTypeLength, unsigned char **result);
void convertByteArray2IntArray_fast_2b(size_t stepLength, unsigned char* byteArray, size_t byteArrayLength, unsigned char **intArray);
size_t convertIntArray2ByteArray_fast_dynamic(unsigned char* timeStepType, unsigned char resiBitLength, size_t nbEle, unsigned char **bytes);
int getLeftMovingSteps(size_t k, unsigned char resiBitLength);
float computeRangeSize_float(float* oriData, size_t size, float* valueRangeSize, float* medianValue);
double computeRangeSize_double(double* oriData, size_t size, double* valueRangeSize, double* medianValue);
size_t bytesToSize(unsigned char* bytes);
void sizeToBytes(unsigned char* outBytes, size_t size);

int lorenzoPredictorQuant_Cmpr_NoOutlier_float(float* data, int lorenzoMode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, int* out);
int lorenzoPredictorQuant_Cmpr_NoOutlier_double(double* data, int lorenzoMode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, int* out);

int lorenzoPredictorQuant_Decmpr_NoOutlier_float(int* diffQuantData, int lorenzoMode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, float* result);
int lorenzoPredictorQuant_Decmpr_NoOutlier_double(int* diffQuantData, int lorenzoMode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, double* result);

int lorenzoPredictorQuant_Cmpr_NoOutlier(void* data, int dataType, int lorenzoMode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, int* out);
int lorenzoPredictorQuant_Decmpr_NoOutlier(int* diffQuantData, int dataType, int lorenzoMode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, void* result);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _SZ_UTILITY_H  ----- */
