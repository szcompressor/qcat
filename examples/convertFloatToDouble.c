/**
 *  @file Convert.c
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
    char oriFilePath[640], outputFilePath[640];
    
    if(argc < 3)
    {
		printf("Test case: convertFloatToDouble [srcFilePath] [tgtFilePath]\n");
		printf("Example: convertFloatToDouble testfloat_8_8_128.dat testdouble_8_8_128.f64 \n");
		exit(0);
    }
   
    sprintf(oriFilePath, "%s", argv[1]);
    sprintf(outputFilePath, "%s", argv[2]);
   
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
   
   
    float *floatData = NULL;
    size_t nbEle = 0;
    printf("reading data from %s \n", oriFilePath);
    floatData = readFloatData(oriFilePath, &nbEle, &status);

    printf("converting...\n");
    double* doubleData = (double*)malloc(nbEle*sizeof(double));

    size_t i = 0;
    float min = floatData[0];
    float max = floatData[0];
    for(i=0;i<nbEle;i++)
    {	
	float value = floatData[i];
	doubleData[i] = value;
	if(min>value)
		min = value;
	if(max<value)
		max = value;
    }

    printf("min=%f, max=%f\n", min, max);
    for(i=0;i<10;i++)
	    printf("doubleData[%zu]=%f\n", i, doubleData[i]);
    if(status != DA_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
    
    printf("writing data to %s\n", outputFilePath);
    writeDoubleData_inBytes(doubleData, nbEle, outputFilePath, &status);
    free(floatData);
    free(doubleData);
    
    return 0;
}
