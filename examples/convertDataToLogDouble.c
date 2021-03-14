/**
 *  @file test_compress.c
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
#include <math.h>

int main(int argc, char * argv[])
{
    int status = 0;
    char oriFilePath[640], outputFilePath[640];
    
    if(argc < 3)
    {
		printf("Test case: convertDataToLogDouble [srcFilePath] newsize [tgtFilePath]\n");
		printf("Example: convertDataToLogDouble testfloat_8_8_128.dat -1 testfloat_8_8_128.xls\n");
		exit(0);
    }
   
    sprintf(oriFilePath, "%s", argv[1]);
    int newSize = atoi(argv[2]);
    sprintf(outputFilePath, "%s", argv[3]);
   
    double *data = NULL;
    size_t nbEle = 0;
    if(newSize==-1)
	data = readDoubleData(oriFilePath, &nbEle, &status);
    else
    {
	nbEle = newSize;
        data = readDoubleData_systemEndian_k(oriFilePath, nbEle, &status);
    }

    double* newData = (double*)malloc(sizeof(double)*nbEle);

    int i = 0;
    for(i=0;i<nbEle;i++)
    {
	if(data[i]<0) 
		newData[i] = log2(-data[i]);
	else if(data[i]==0) 
		newData[i] = -100;
	else
    		newData[i] = log2(data[i]);
    }
    for(i=0;i<10;i++)
	    printf("data[%d]=%.20G --> %.20G\n", i, data[i], newData[i]);
    if(status != RW_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
    writeDoubleData_inBytes(newData, nbEle, outputFilePath, &status);
    printf("output file: %s\n", outputFilePath);
    free(data);
    free(newData);
    
    return 0;
}
