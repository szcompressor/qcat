/**
 *  @file changeDataSize.c
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
    char type[3];
    if(argc < 2)
    {
		printf("Test case: changeDataSize [srcFilePath] [type(-f/-d)] newsize\n");
		printf("Example: changeDataSize testfloat_8_8_128.dat -f 200\n");
		exit(0);
    }
   
    int newSize = 0;
    if(argc==2)
    {
    	sprintf(oriFilePath, "%s", argv[1]);
    	newSize = atoi(argv[2]);
    	strcpy(type, "-f");
    }
    else
    {	    
    	sprintf(oriFilePath, "%s", argv[1]);
    	sprintf(type, "%s", argv[2]);
    	newSize = atoi(argv[3]);
    }
    size_t nbEle;

    if(strcmp(type, "-f")==0)
    {
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
    }
    else
    {
	double *data = readDoubleData(oriFilePath, &nbEle, &status);
	int i = 0;
	for(i=0;i<10;i++)
		printf("data[%d]=%.20G\n", i, data[i]);
	if(status != RW_SCES)
	{
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
	}
	char outputFile[644];
    	sprintf(outputFile, "%s.%d", oriFilePath, newSize);
    	writeDoubleData_inBytes(data, newSize, outputFile, &status);
    	printf("outputFile= %s\n", outputFile);
    	free(data);
    }
	
    return 0;
}
