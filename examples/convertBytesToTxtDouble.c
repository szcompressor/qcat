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
    char oriFilePath[640], outputFilePath[640];
    
    if(argc < 3)
    {
		printf("Test case: convertBytesToTxtDouble [srcFilePath] newsize [tgtFilePath]\n");
		printf("Example: convertBytesToTxtDouble testfloat_8_8_128.dat 200 testfloat_8_8_128.xls\n");
		exit(0);
    }
   
    sprintf(oriFilePath, "%s", argv[1]);
    int newSize = atoi(argv[2]);
    sprintf(outputFilePath, "%s", argv[3]);
   
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
	 printf("This is little-endian system.\n");
    }
   
   
    double *data = NULL;
    size_t nbEle = 0;
    if(newSize==-1)
		data = readDoubleData(oriFilePath, &nbEle, &status);
    else
    {
	nbEle = newSize;
        data = readDoubleData_systemEndian_k(oriFilePath, nbEle, &status);
    }
    int i = 0;
    for(i=0;i<10;i++)
	    printf("data[%d]=%f\n", i, data[i]);
    if(status != RW_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
    writeDoubleData(data, nbEle, outputFilePath, &status);
    free(data);
    
    return 0;
}
