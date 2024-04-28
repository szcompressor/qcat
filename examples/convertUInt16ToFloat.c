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
		printf("Test case: convertDoubleToFloat [srcFilePath] [tgtFilePath]\n");
		printf("Example: convertDoubleToFloat testdouble_8_8_128.dat testfloat_8_8_128.dat \n");
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
   
   
    uint16_t *uint16Data = NULL;
    size_t nbEle = 0;
    printf("reading data from %s \n", oriFilePath);

    uint16Data = readUInt16Data(oriFilePath, &nbEle, &status);

    printf("converting...\n");
    float* floatData = (float*)malloc(nbEle*sizeof(float));

    size_t i = 0;
    uint16_t min = uint16Data[0];
    uint16_t max = uint16Data[0];
    for(i=0;i<nbEle;i++)
    {	
	uint16_t value = uint16Data[i];
	floatData[i] = value;
	if(min>value)
		min = value;
	if(max<value)
		max = value;
    }

    printf("min=%hu, max=%hu\n", min, max);
    for(i=0;i<10;i++)
	    printf("floatData[%zu]=%f\n", i, floatData[i]);
    if(status != RW_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
    
    printf("writing data to %s\n", outputFilePath);
    writeFloatData_inBytes(floatData, nbEle, outputFilePath, &status);
    free(floatData);
    free(uint16Data);
    
    return 0;
}
