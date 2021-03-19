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
    char oriFilePath[640];
    
    if(argc < 2)
    {
		printf("Test case: changeDataSize [srcFilePath] newsize\n");
		printf("Example: changeDataSize testfloat_8_8_128.dat 200\n");
		exit(0);
    }
   
    sprintf(oriFilePath, "%s", argv[1]);
    int newSize = atoi(argv[2]);
   
    size_t nbEle;
    float *data = readFloatData(oriFilePath, &nbEle, &status);
    int i = 0;
    for(i=0;i<10;i++)
	    printf("data[%d]=%f\n", i, data[i]);
    if(status != RW_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
   
    char outputFile[644];
    sprintf(outputFile, "%s.%d", oriFilePath, newSize);
    writeFloatData_inBytes(data, newSize, outputFile, &status);
    printf("outputFile= %s\n", outputFile);
    free(data);
    
    return 0;
}
