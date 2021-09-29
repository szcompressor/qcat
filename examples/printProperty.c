/**
 *  @file printProperty.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include "ByteToolkit.h"
#include "rw.h"
#include "qcat_dataAnalysis.h"
#include "math.h"

int main(int argc, char * argv[])
{
	int status = 0;
	char oriFilePath[640];

	if(argc < 3)
	{
		printf("Usage: printDataProperty [dataType (-f or -d)] [tgtFilePath] [entropyType]\n");
		printf("As for [entropyType]:\n");
		printf("\t\t0 means don't compute entropy at all.\n");	
		printf("\t\t1 means compute only 8bit byte entropy.\n");
		printf("\t\t2 means compute only floating-point entropy.\n");
		printf("\t\t3 means compute both 8bit and floating-point entropy.\n");
		printf("Example: printDataProperty -f testfloat_8_8_128.dat 1\n");
		exit(0);
	}

	int dataType = strcmp(argv[1],"-f") == 0 ? QCAT_FLOAT : QCAT_DOUBLE; //0: float , 1: double
	sprintf(oriFilePath, "%s", argv[2]);
	int entropyType = atoi(argv[3]);

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

	unsigned char* data  = NULL;
	size_t nbEle = 0;
	if(dataType == QCAT_FLOAT)
		data = (unsigned char*)readFloatData(oriFilePath, &nbEle, &status);
	else
		data = (unsigned char*)readDoubleData(oriFilePath, &nbEle, &status);

	int i = 0;
	printf("The first 10 values are: \n");
	if(dataType == QCAT_FLOAT)
		for(i=0;i<10;i++)
			printf("%f ", ((float*)data)[i]);
	else
		for(i=0;i<10;i++)
			printf("%f ", ((double*)data)[i]);		
	printf("....\n------------------------\n");
	if(status != RW_SCES)
	{
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
	}
	
	QCAT_DataProperty* property = computeProperty(dataType, data, nbEle, entropyType);
	
	printProperty(property);
	
	free(property);
	free(data);

	return 0;
}
