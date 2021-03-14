/**
 *  @file sz_dummy_compression.h
 *  @author Sheng Di
 *  @date March, 2021
 *  @brief Header file for sz_dummy_compression.c.
 *  (C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _SZ_DUMMY_COMPRESSION_H
#define _SZ_DUMMY_COMPRESSION_H

#ifdef __cplusplus
extern "C" {
#endif

QCAT_CompressionResult* huffmanAndZstd(int dataType, int* type, int quantBinCapacity, size_t nbEle, void* origData, void* decData);
unsigned char* SZ_fast_compress_args_unpredictable_float(int dataType, float *data, size_t *outSize, float absErrBound, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
void SZ_fast_decompress_args_unpredictable_float(float** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, 
size_t cmpSize);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _SZ_DUMMY_COMPRESSION_H  ----- */
