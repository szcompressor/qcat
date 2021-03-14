/**
 *  @file predCR.c
 *  @author Sheng Di
 *  @date March, 2021
 *  @brief This is an example of using compression interface
 *  (C) 2020 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ByteToolkit.h"
#include "rw.h"
#include "qcat.h"

int main(int argc, char * argv[])
{
	int status = 0;
	char dataType[4];
	char oriFilePath[640], preFilePath[640];

	if(argc < 3)
	{
		printf("Test case: predCR [datatype (-f or -d)] [quantBinCapacity] [errorBound] [original data file] [predcted data file]\n");
		printf("			-f means single precision; -d means double precision\n");
		printf("Example: predCR -f 1024 1E-1 original.dat predicted.dat \n");
		exit(0);
	}

	sprintf(dataType, "%s", argv[1]);
	int quantBinCapacity = atoi(argv[2]);
	float errBound = atof(argv[3]);
	sprintf(oriFilePath, "%s", argv[4]);
	sprintf(preFilePath, "%s", argv[5]);

	size_t nbEle;
	size_t nbEle2;
	void* oriData;
	void* predData;

	int x = 1;
	char *y = (char*)&x;

	if(*y==1)
	{     
		sysEndianType = 0; //LITTLE_ENDIAN_SYSTEM;
		printf("This is little-endian system.\n");
	}
	else //=0
	{
		sysEndianType = 1; //BIG_ENDIAN_SYSTEM;
		printf("This is big-endian system.\n");
	}
   
	int dType = 0;
	if(strcmp(dataType, "-f")==0)
	{
		dType = QCAT_FLOAT;
		oriData = readFloatData(oriFilePath, &nbEle, &status);
		predData = readFloatData(preFilePath, &nbEle2, &status);		
	}
	else
	{
		dType = QCAT_DOUBLE;
		oriData = readDoubleData(oriFilePath, &nbEle, &status);
		predData = readDoubleData(preFilePath, &nbEle2, &status);		
	}
	
	QCAT_DataProperty* property = computeProperty(dType, oriData, nbEle);	
	printProperty(property);
	
	QCAT_CompressionResult* result = getCompressionResult(dType, errBound, quantBinCapacity, oriData, predData, property);
	printCompressionResult(result);
   
	free(oriData);
	free(predData);
	free(property);
	free(result);
    
	return 0;
}
