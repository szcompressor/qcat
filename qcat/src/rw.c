/**
 *  @file rw.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief io interface for fortrance
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ByteToolkit.h>
#include <DynamicFloatArray.h>
#include <DynamicDoubleArray.h>
#include <libgen.h>
#include "rw.h"
#include "qcat.h"

int dataEndianType = 0;
int sysEndianType = 0; //0 means little endian, 1 means big endian

char *removeFileExtension(char* myStr) 
{
    char *retStr;
    char *lastExt;
    if (myStr == NULL) return NULL;
    if ((retStr = malloc (strlen (myStr) + 1)) == NULL) return NULL;
    strcpy (retStr, myStr);
    lastExt = strrchr (retStr, '.');
    if (lastExt != NULL)
        *lastExt = '\0';
    return retStr;
}

void symTransform_8bytes(unsigned char data[8])
{
        unsigned char tmp = data[0];
        data[0] = data[7];
        data[7] = tmp;

        tmp = data[1];
        data[1] = data[6];
        data[6] = tmp;

        tmp = data[2];
        data[2] = data[5];
        data[5] = tmp;

        tmp = data[3];
        data[3] = data[4];
        data[4] = tmp;
}

inline void symTransform_2bytes(unsigned char data[2])
{
        unsigned char tmp = data[0];
        data[0] = data[1];
        data[1] = tmp;
}

inline void symTransform_4bytes(unsigned char data[4])
{
        unsigned char tmp = data[0];
        data[0] = data[3];
        data[3] = tmp;

        tmp = data[1];
        data[1] = data[2];
        data[2] = tmp;
}


int checkFileExistance(char* filePath)
{
	if( access( filePath, F_OK ) != -1 ) {
		// file exists
		return 1;
	} else {
		// file doesn't exist
		return 0;
	}	
}

float** create2DArray_float(size_t m, size_t n)
{
	size_t i=0;
	float **data = (float**)malloc(sizeof(float*)*m);
	for(i=0;i<m;i++)
		data[i] = (float*)malloc(sizeof(float)*n);
	return data;
}

void free2DArray_float(float** data, size_t m)
{
	size_t i = 0;
	for(i=0;i<m;i++)
		free(data[i]);
	free(data);	
}

float*** create3DArray_float(size_t p, size_t m, size_t n)
{
	size_t i = 0, j = 0;
	float ***data = (float***)malloc(sizeof(float**)*m);
	for(i=0;i<p;i++)
	{
		data[i] = (float**)malloc(sizeof(float*)*n);
		for(j=0;j<m;j++)
			data[i][j] = (float*)malloc(sizeof(float)*n);
	}
	return data;
}

void free3DArray_float(float*** data, size_t p, size_t m)
{
	size_t i,j;
	for(i=0;i<p;i++)
	{
		for(j=0;j<m;j++)
			free(data[i][j]);
		free(data[i]);
	}
	free(data);	
}

double** create2DArray_double(size_t m, size_t n)
{
	size_t i=0;
	double **data = (double**)malloc(sizeof(double*)*m);
	for(i=0;i<m;i++)
			data[i] = (double*)malloc(sizeof(double)*n);
			
	return data;
}

void free2DArray_double(double** data, size_t m)
{
	size_t i;
	for(i=0;i<m;i++)
		free(data[i]);
	free(data);	
}

double*** create3DArray_double(size_t p, size_t m, size_t n)
{
	size_t i = 0, j = 0;
	double ***data = (double***)malloc(sizeof(double**)*m);
	for(i=0;i<p;i++)
	{
		data[i] = (double**)malloc(sizeof(double*)*n);
		for(j=0;j<m;j++)
			data[i][j] = (double*)malloc(sizeof(double)*n);
	}
	return data;
}

void free3DArray_double(double*** data, size_t p, size_t m)
{
	size_t i,j;
	for(i=0;i<p;i++)
	{
		for(j=0;j<m;j++)
			free(data[i][j]);
		free(data[i]);
	}
	free(data);	
}

size_t checkFileSize(char *srcFilePath, int *status)
{
	size_t filesize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = RW_FERR;
		return -1;
	}
	fseek(pFile, 0, SEEK_END);
    filesize = ftell(pFile);
    fclose(pFile);
    *status = RW_SCES;
    return filesize;
}

unsigned char *readByteData(char *srcFilePath, size_t *byteLength, int *status)
{
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        *status = RW_FERR;
        return 0;
    }
	fseek(pFile, 0, SEEK_END);
    *byteLength = ftell(pFile);
    fclose(pFile);
    
    unsigned char *byteBuf = ( unsigned char *)malloc((*byteLength)*sizeof(unsigned char)); //sizeof(char)==1
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        *status = RW_FERR;
        return 0;
    }
    fread(byteBuf, 1, *byteLength, pFile);
    fclose(pFile);
    *status = RW_SCES;
    return byteBuf;
}

double *readDoubleData(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = RW_SCES;
	if(dataEndianType==sysEndianType)
	{
		double *daBuf = readDoubleData_systemEndian(srcFilePath, nbEle,&state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;
		
		size_t byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state==RW_FERR)
		{
			*status = RW_FERR;
			return NULL;
		}
		double *daBuf = (double *)malloc(byteLength);
		*nbEle = byteLength/8;
		
		ldouble buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*8;
			memcpy(buf.byte, bytes+j, 8);
			symTransform_8bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		return daBuf;
	}
}


int8_t *readInt8Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = RW_SCES;
	int8_t *daBuf = readInt8Data_systemEndian(srcFilePath, nbEle, &state);
	*status = state;
	return daBuf;
}

int16_t *readInt16Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = RW_SCES;
	if(dataEndianType==sysEndianType)
	{
		int16_t *daBuf = readInt16Data_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state == RW_FERR)
		{
			*status = RW_FERR;
			return NULL;
		}
		int16_t *daBuf = (int16_t *)malloc(byteLength);
		*nbEle = byteLength/2;

		lint16 buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i << 1;//*2
			memcpy(buf.byte, bytes+j, 2);
			symTransform_2bytes(buf.byte);
			daBuf[i] = buf.svalue;
		}
		free(bytes);
		return daBuf;
	}
}

uint16_t *readUInt16Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = RW_SCES;
	if(dataEndianType==sysEndianType)
	{
		uint16_t *daBuf = readUInt16Data_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state == RW_FERR)
		{
			*status = RW_FERR;
			return NULL;
		}
		uint16_t *daBuf = (uint16_t *)malloc(byteLength);
		*nbEle = byteLength/2;

		lint16 buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i << 1;//*2
			memcpy(buf.byte, bytes+j, 2);
			symTransform_2bytes(buf.byte);
			daBuf[i] = buf.usvalue;
		}
		free(bytes);
		return daBuf;
	}
}

int32_t *readInt32Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = RW_SCES;
	if(dataEndianType==sysEndianType)
	{
		int32_t *daBuf = readInt32Data_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state == RW_FERR)
		{
			*status = RW_FERR;
			return NULL;
		}
		int32_t *daBuf = (int32_t *)malloc(byteLength);
		*nbEle = byteLength/4;

		lint32 buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*4;
			memcpy(buf.byte, bytes+j, 4);
			symTransform_4bytes(buf.byte);
			daBuf[i] = buf.ivalue;
		}
		free(bytes);
		return daBuf;
	}
}

uint32_t *readUInt32Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = RW_SCES;
	if(dataEndianType==sysEndianType)
	{
		uint32_t *daBuf = readUInt32Data_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state == RW_FERR)
		{
			*status = RW_FERR;
			return NULL;
		}
		uint32_t *daBuf = (uint32_t *)malloc(byteLength);
		*nbEle = byteLength/4;

		lint32 buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i << 2; //*4
			memcpy(buf.byte, bytes+j, 4);
			symTransform_4bytes(buf.byte);
			daBuf[i] = buf.uivalue;
		}
		free(bytes);
		return daBuf;
	}
}

int64_t *readInt64Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = RW_SCES;
	if(dataEndianType==sysEndianType)
	{
		int64_t *daBuf = readInt64Data_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state == RW_FERR)
		{
			*status = RW_FERR;
			return NULL;
		}
		int64_t *daBuf = (int64_t *)malloc(byteLength);
		*nbEle = byteLength/8;

		lint64 buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i << 3; //*8
			memcpy(buf.byte, bytes+j, 8);
			symTransform_8bytes(buf.byte);
			daBuf[i] = buf.lvalue;
		}
		free(bytes);
		return daBuf;
	}
}

uint64_t *readUInt64Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = RW_SCES;
	if(dataEndianType==sysEndianType)
	{
		uint64_t *daBuf = readUInt64Data_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state == RW_FERR)
		{
			*status = RW_FERR;
			return NULL;
		}
		uint64_t *daBuf = (uint64_t *)malloc(byteLength);
		*nbEle = byteLength/8;

		lint64 buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i << 3; //*8
			memcpy(buf.byte, bytes+j, 8);
			symTransform_8bytes(buf.byte);
			daBuf[i] = buf.ulvalue;
		}
		free(bytes);
		return daBuf;
	}
}


float *readFloatData(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = RW_SCES;
	if(dataEndianType==sysEndianType)
	{
		float *daBuf = readFloatData_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;
		
		size_t byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state == RW_FERR)
		{
			*status = RW_FERR;
			return NULL;
		}
		float *daBuf = (float *)malloc(byteLength);
		*nbEle = byteLength/4;
		
		lfloat buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*4;
			memcpy(buf.byte, bytes+j, 4);
			symTransform_4bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		return daBuf;
	}
}

double *readDoubleData_systemEndian_k(char *srcFilePath, size_t nbEle, int *status)
{
	size_t inSize = sizeof(double)*nbEle;
    
    double *daBuf = (double *)malloc(inSize);
    
	int state = RW_SCES;
	if(dataEndianType==sysEndianType)
	{
		FILE* pFile = fopen(srcFilePath, "rb");
		if (pFile == NULL)
		{
			printf("Failed to open input file. 2\n");
			*status = RW_FERR;
			return NULL;
		}
		fread(daBuf, sizeof(double), nbEle, pFile);
		fclose(pFile);
		*status = RW_SCES;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state == RW_FERR)
		{
			*status = RW_FERR;
			return NULL;
		}

		ldouble buf;
		for(i = 0;i<nbEle;i++)
		{
			j = i << 3; //*8
			memcpy(buf.byte, bytes+j, 8);
			symTransform_8bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		*status = RW_SCES;
		return daBuf;
	}    
}

float *readFloatData_systemEndian_k(char *srcFilePath, size_t nbEle, int *status)
{
	size_t inSize = sizeof(float)*nbEle;
    
    float *daBuf = (float *)malloc(inSize);


	int state = RW_SCES;
	if(dataEndianType==sysEndianType)
	{
		FILE* pFile = fopen(srcFilePath, "rb");
		if (pFile == NULL)
		{
			printf("Failed to open input file. 2\n");
			*status = RW_FERR;
			return NULL;
		}
		fread(daBuf, sizeof(float), nbEle, pFile);
		fclose(pFile);
		*status = RW_SCES;
		return daBuf;
	}
	else
	{
		size_t i,j;
		
		size_t byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state == RW_FERR)
		{
			*status = RW_FERR;
			return NULL;
		}
		
		lfloat buf;
		for(i = 0;i<nbEle;i++)
		{
			j = i*4;
			memcpy(buf.byte, bytes+j, 4);
			symTransform_4bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		*status = RW_SCES;
		return daBuf;
	}
}

float *readFloatData_systemEndian_sk(char *srcFilePath, size_t startLocation, size_t nbEle, int *status)
{
	size_t inSize = sizeof(float)*(startLocation+nbEle);
    
    float *daBuf = (float *)malloc(inSize);


	int state = RW_SCES;
	if(dataEndianType==sysEndianType)
	{
		FILE* pFile = fopen(srcFilePath, "rb");
		if (pFile == NULL)
		{
			printf("Failed to open input file. 2\n");
			*status = RW_FERR;
			return NULL;
		}
		fread(daBuf, sizeof(float), (startLocation+nbEle), pFile);
		fclose(pFile);
		*status = RW_SCES;
		
		size_t i;
		for(i=0;i<nbEle;i++)
		{
			daBuf[i] = daBuf[i+startLocation];
		}
	}
	else
	{
		size_t i,j;
		
		size_t byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state == RW_FERR)
		{
			*status = RW_FERR;
			return NULL;
		}
		
		lfloat buf;
		for(i = startLocation;i<nbEle;i++)
		{
			j = i*4;
			memcpy(buf.byte, bytes+j, 4);
			symTransform_4bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		*status = RW_SCES;
	}
	return daBuf;
}

double *readDoubleData_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        *status = RW_FERR;
        return NULL;
    }
	fseek(pFile, 0, SEEK_END);
    inSize = ftell(pFile);
    *nbEle = inSize/8; //only support double in this version
    fclose(pFile);
    
    double *daBuf = (double *)malloc(inSize);
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        *status = RW_FERR;
        return NULL;
    }
    fread(daBuf, 8, *nbEle, pFile);
    fclose(pFile);
    *status = RW_SCES;
    return daBuf;
}


int8_t *readInt8Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = RW_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize;
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = RW_FERR;
	}

	int8_t *daBuf = (int8_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = RW_FERR;
		return NULL;
	}
	fread(daBuf, 1, *nbEle, pFile);
	fclose(pFile);
	*status = RW_SCES;
	return daBuf;
}


int16_t *readInt16Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = RW_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize/2; 
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = RW_FERR;
	}

	int16_t *daBuf = (int16_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = RW_FERR;
		return NULL;
	}
	fread(daBuf, 2, *nbEle, pFile);
	fclose(pFile);
	*status = RW_SCES;
	return daBuf;	
}

uint16_t *readUInt16Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = RW_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize/2; 
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = RW_FERR;
	}

	uint16_t *daBuf = (uint16_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = RW_FERR;
		return NULL;
	}
	fread(daBuf, 2, *nbEle, pFile);
	fclose(pFile);
	*status = RW_SCES;
	return daBuf;	
}

int32_t *readInt32Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = RW_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize/4; 
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = RW_FERR;
	}

	int32_t *daBuf = (int32_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = RW_FERR;
		return NULL;
	}
	fread(daBuf, 4, *nbEle, pFile);
	fclose(pFile);
	*status = RW_SCES;
	return daBuf;	
}

uint32_t *readUInt32Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = RW_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize/4; 
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = RW_FERR;
	}

	uint32_t *daBuf = (uint32_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = RW_FERR;
		return NULL;
	}
	fread(daBuf, 4, *nbEle, pFile);
	fclose(pFile);
	*status = RW_SCES;
	return daBuf;	
}

int64_t *readInt64Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = RW_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize/8; 
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = RW_FERR;
	}

	int64_t *daBuf = (int64_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = RW_FERR;
		return NULL;
	}
	fread(daBuf, 8, *nbEle, pFile);
	fclose(pFile);
	*status = RW_SCES;
	return daBuf;
}

uint64_t *readUInt64Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = RW_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize/8; 
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = RW_FERR;
	}

	uint64_t *daBuf = (uint64_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = RW_FERR;
		return NULL;
	}
	fread(daBuf, 8, *nbEle, pFile);
	fclose(pFile);
	*status = RW_SCES;
	return daBuf;
}

float *readFloatData_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        *status = RW_FERR;
        return NULL;
    }
	fseek(pFile, 0, SEEK_END);
    inSize = ftell(pFile);
    *nbEle = inSize/4; 
    fclose(pFile);
    
    if(inSize<=0)
    {
		printf("Error: input file is wrong!\n");
		*status = RW_FERR;
	}
    
    float *daBuf = (float *)malloc(inSize);
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        *status = RW_FERR;
        return NULL;
    }
    fread(daBuf, 4, *nbEle, pFile);
    fclose(pFile);
    *status = RW_SCES;
    return daBuf;
}

void writeByteData(unsigned char *bytes, size_t byteLength, char *tgtFilePath, int *status)
{
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        *status = RW_FERR;
        return;
    }
    
    fwrite(bytes, 1, byteLength, pFile); //write outSize bytes
    fclose(pFile);
    *status = RW_SCES;
}

void writeDoubleData(double *data, size_t nbEle, char *tgtFilePath, int *status)
{
	size_t i = 0;
	char s[64];
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        *status = RW_FERR;
        return;
    }
    
    for(i = 0;i<nbEle;i++)
	{
		sprintf(s,"%.20G\n",data[i]);
		fputs(s, pFile);
	}
    
    fclose(pFile);
    *status = RW_SCES;
}

void writeFloatData(float *data, size_t nbEle, char *tgtFilePath, int *status)
{
	size_t i = 0;
	char s[64];
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        *status = RW_FERR;
        return;
    }
   
    for(i = 0;i<nbEle;i++)
	{
		//printf("i=%d\n",i);
		//printf("data[i]=%f\n",data[i]);
		sprintf(s,"%.30G\n",data[i]);
		fputs(s, pFile);
	}
    
    fclose(pFile);
    *status = RW_SCES;
}

void writeData(void *data, int dataType, size_t nbEle, char *tgtFilePath, int *status)
{
	int state = RW_SCES;
	if(dataType == QCAT_FLOAT)
	{
		float* dataArray = (float *)data;
		writeFloatData(dataArray, nbEle, tgtFilePath, &state);
	}
	else if(dataType == QCAT_DOUBLE)
	{
		double* dataArray = (double *)data;
		writeDoubleData(dataArray, nbEle, tgtFilePath, &state);	
	}
	else
	{
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		*status = RW_TERR; //wrong type
		return;
	}
	*status = state;
}

void writeFloatData_inBytes(float *data, size_t nbEle, char* tgtFilePath, int *status)
{
	size_t i = 0; 
	int state = RW_SCES;
	lfloat buf;
	unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(float));
	for(i=0;i<nbEle;i++)
	{
		buf.value = data[i];
		bytes[i*4+0] = buf.byte[0];
		bytes[i*4+1] = buf.byte[1];
		bytes[i*4+2] = buf.byte[2];
		bytes[i*4+3] = buf.byte[3];					
	}

	size_t byteLength = nbEle*sizeof(float);
	writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void writeDoubleData_inBytes(double *data, size_t nbEle, char* tgtFilePath, int *status)
{
	size_t i = 0, index = 0; 
	int state = RW_SCES;
	ldouble buf;
	unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(double));
	for(i=0;i<nbEle;i++)
	{
		index = i*8;
		buf.value = data[i];
		bytes[index+0] = buf.byte[0];
		bytes[index+1] = buf.byte[1];
		bytes[index+2] = buf.byte[2];
		bytes[index+3] = buf.byte[3];
		bytes[index+4] = buf.byte[4];
		bytes[index+5] = buf.byte[5];
		bytes[index+6] = buf.byte[6];
		bytes[index+7] = buf.byte[7];
	}

	size_t byteLength = nbEle*sizeof(double);
	writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void writeShortData_inBytes(short *states, size_t stateLength, char *tgtFilePath, int *status)
{
	int state = RW_SCES;
	size_t byteLength = stateLength*2;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertShortArrayToBytes(states, stateLength, bytes);
	writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void writeUShortData_inBytes(unsigned short *states, size_t stateLength, char *tgtFilePath, int *status)
{
	int state = RW_SCES;
	size_t byteLength = stateLength*2;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertUShortArrayToBytes(states, stateLength, bytes);
	writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void writeIntData_inBytes(int *states, size_t stateLength, char *tgtFilePath, int *status)
{
	int state = RW_SCES;
	size_t byteLength = stateLength*4;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertIntArrayToBytes(states, stateLength, bytes);
	writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void writeUIntData_inBytes(unsigned int *states, size_t stateLength, char *tgtFilePath, int *status)
{
	int state = RW_SCES;
	size_t byteLength = stateLength*4;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertUIntArrayToBytes(states, stateLength, bytes);
	writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void writeLongData_inBytes(int64_t *states, size_t stateLength, char *tgtFilePath, int *status)
{
	int state = RW_SCES;
	size_t byteLength = stateLength*8;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertLongArrayToBytes(states, stateLength, bytes);
	writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void writeULongData_inBytes(uint64_t *states, size_t stateLength, char *tgtFilePath, int *status)
{
	int state = RW_SCES;
	size_t byteLength = stateLength*8;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertULongArrayToBytes(states, stateLength, bytes);
	writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

unsigned short* readShortData(char *srcFilePath, size_t *dataLength, int *status)
{
	size_t byteLength = 0; 
	int state = RW_SCES;
	unsigned char * bytes = readByteData(srcFilePath, &byteLength, &state);
	*dataLength = byteLength/2;
	unsigned short* states = convertByteDataToUShortArray(bytes, byteLength);
	free(bytes);
	*status = state;
	return states;
}

void writeStrings(int nbStr, char *str[], char *tgtFilePath, int *status)
{
	size_t i = 0;
	char s[256];
	FILE *pFile = fopen(tgtFilePath, "wb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 3\n");
		*status = RW_FERR;
		return;
	}

	for(i = 0;i<nbStr;i++)
	{
		sprintf(s,"%s\n",str[i]);
		fputs(s, pFile);
	}

	fclose(pFile);
	*status = RW_SCES;
}

double *readDoubleData_inTxt(char* filePath, size_t* nbEle) 
{
	DynamicDoubleArray* dda = NULL;
	new_DDA(&dda, 1024);
	char buffer[255];
	FILE *fp = fopen(filePath, "r");
    if (fp == NULL)
	{
		printf("Error: please check if the file path: %s is correct.\n", filePath);
		exit(0);
	}

	while(fgets(buffer, 255, (FILE*) fp)) {
		double value = atof(buffer);
		addDDA_Data(dda, value);
	}

	*nbEle = dda->size;
	double* result = (double*)malloc(sizeof(double)*(*nbEle));
	memcpy(result, dda->array, sizeof(double)*(*nbEle));
	free_DDA(dda);
	
    fclose(fp);
	return result;
}

float *readFloatData_inTxt(char* filePath, size_t* nbEle) 
{
	DynamicFloatArray* dfa = NULL;
	new_DFA(&dfa, 1024);
	char buffer[255];
	FILE *fp = fopen(filePath, "r");
    if (fp == NULL)
	{
		printf("Error: please check if the file path: %s is correct.\n", filePath);
		exit(0);
	}

	while(fgets(buffer, 255, (FILE*) fp)) {
		float value = (float)atof(buffer);
		addDFA_Data(dfa, value);
	}

	*nbEle = dfa->size;
	float* result = (float*)malloc(sizeof(float)*(*nbEle));
	memcpy(result, dfa->array, sizeof(float)*(*nbEle));
	free_DFA(dfa);
	
    fclose(fp);
	return result;
}

void RW_writeFloatData_gnuplotImage(float *data, size_t r2, size_t r1, char *tgtFilePath)
{
	char s[64];
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        exit(1);
    }
   
    size_t i = 0, j = 0;
    for(i = 0;i<r2;i++)
	{
		for(j=0;j<r1;j++)
		{
			sprintf(s, "%.10G\n", data[i*r1+j]);
			fputs(s, pFile);			
		}
		if(i<r2-1)
			fputs("\n",pFile);
	}
    
    fclose(pFile);	
}

void RW_writeDoubleData_gnuplotImage(double *data, size_t r2, size_t r1, char *tgtFilePath)
{
	char s[64];
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        exit(1);
    }
   
    size_t i = 0, j = 0;
    for(i = 0;i<r2;i++)
	{
		for(j=0;j<r1;j++)
		{
			sprintf(s, "%.10G\n", data[i*r1+j]);
			fputs(s, pFile);			
		}
		if(i<r2-1)
			fputs("\n",pFile);
	}
    
    fclose(pFile);	
}

void RW_writeData_genuplotImage(void *data, int dataType, size_t r2, size_t r1, char *tgtFilePath)
{
	if(dataType == 0)
	{
		float* dataArray = (float *)data;	
		RW_writeFloatData_gnuplotImage(dataArray, r2, r1, tgtFilePath);
	}
	else if(dataType == 1)
	{
		double* dataArray = (double *)data;
		RW_writeDoubleData_gnuplotImage(dataArray, r2, r1, tgtFilePath);
	}
	else
	{
		printf("Error: data type cannot be the types other than ZC_FLOAT or ZC_DOUBLE\n");
		exit(0);
	}
}

int RW_writeStrings(int string_size, char **string, char *tgtFilePath)
{
	size_t i = 0;
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        exit(1);
    }
   
    for(i = 0;i<string_size;i++)
	{
		if(string[i]==0)
			break;
		fputs(string[i], pFile);
	}
    
    fclose(pFile);
    
    return i;	
   
}

char* extractDirFromPath(char* filePath)
{
	if(strstr(filePath, "/")==NULL) 
		return NULL;
	char* ts1 = strdup(filePath);
	char* dir = dirname(ts1);
	return dir;
}

char *extractFileNameFromPath(char *filePath)
{
    char ch = '/';
	if(strstr(filePath, "/")==NULL) 
		return filePath;    
    char *q = strrchr(filePath,ch) + 1;
    return q;
}

void writePDFData(char* tgtFilePath, double err_minValue, double err_interval, int pdf_intervals, double* pdfData)
{
	size_t i = 0;
	if(err_interval==0)
	{
		char *ss[2];
		ss[0] = (char*)malloc(sizeof(char)*QCAT_BUFS);
		sprintf(ss[0], "x errpdf\n");
		ss[1] = (char*)malloc(sizeof(char)*QCAT_BUFS);
		strcpy(ss[1],"0 1\n");
		RW_writeStrings(2, ss, tgtFilePath);
		free(ss[0]);
		free(ss[1]);
	}
	else
	{
		char *ss[pdf_intervals+1];
		ss[0] = (char*)malloc(sizeof(char)*QCAT_BUFS);
		sprintf(ss[0], "x errpdf\n");
		for(i=0;i<pdf_intervals;i++)
		{
			//printf("%d\n", i);
			ss[i+1] = (char*)malloc(sizeof(char)*QCAT_BUFS);
			double x = err_minValue+i*err_interval;
			sprintf(ss[i+1], "%.10G %.10G\n", x, pdfData[i]);
		}
		RW_writeStrings(pdf_intervals+1, ss, tgtFilePath);
		for(i=0;i<pdf_intervals+1;i++)
			free(ss[i]);
	}
	
}

void writePDFData_int32(char* tgtFilePath, double min, int intervals, double* pdfData)
{
	size_t i = 0;
	char** ss = (char**)malloc((intervals+1)*sizeof(char*));
	ss[0] = (char*)malloc(sizeof(char)*QCAT_BUFS);
	sprintf(ss[0], "x errpdf\n");
	for(i=0;i<intervals;i++)
	{
		//printf("%d\n", i);
		ss[i+1] = (char*)malloc(sizeof(char)*QCAT_BUFS);
		double x = min+i;
		if(pdfData[i]!=0)
			sprintf(ss[i+1], "%.10G %.10G\n", x, pdfData[i]);
	}
	RW_writeStrings(intervals+1, ss, tgtFilePath);
	for(i=0;i<intervals+1;i++)
		free(ss[i]);
		
	free(ss);
	
}
