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

int main(int argc, char * argv[])
{
    int status = 0;
    char oriFilePath[640], realDataPath[640], imagDataPath[640];
    
    if(argc < 2)
    {
		printf("Test case: splitComplexDataDouble [srcFilePath] newsize\n");
		printf("Example: splitComplexDataDouble testdouble_8_8_128.dat 200\n");
		exit(0);
    }
   
    sprintf(oriFilePath, "%s", argv[1]);
    sprintf(realDataPath, "%s.real", argv[1]);
    sprintf(imagDataPath, "%s.imag", argv[1]);
    int newSize = atoi(argv[2]);
   
    size_t nbEle;
    double *data = readDoubleData(oriFilePath, &nbEle, &status);
    int i = 0;
    for(i=0;i<10;i++)
	    printf("data[%d]=%f\n", i, data[i]);
    if(status != RW_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
    size_t j = 0;
    if(newSize>0)
	nbEle = newSize;
    else
	nbEle/=2;
    double *realData = (double*)malloc(sizeof(double)*nbEle);
    double *imagData = (double*)malloc(sizeof(double)*nbEle);
    for(j=0;j<nbEle;j++)
    {
	realData[j] = data[j*2];
	imagData[j] = data[j*2+1];
    }   

    writeDoubleData_inBytes(realData, newSize, realDataPath, &status);
    writeDoubleData_inBytes(imagData, newSize, imagDataPath, &status);
    
    free(data);
    free(realData);
    free(imagData);
    
    return 0;
}

