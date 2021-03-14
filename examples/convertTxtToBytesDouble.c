/**
 *  @file convertTxtToBytesDouble.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief 
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
		printf("Test case: convertTxtToBytesDouble [srcFilePath] [tgtFilePath]\n");
		printf("Example: convertTxtToBytesDouble testfloat_8_8_128.txt testfloat_8_8_128.dat\n");
		exit(0);
    }
   
    sprintf(oriFilePath, "%s", argv[1]);
    sprintf(outputFilePath, "%s", argv[2]);
   
    double *data = NULL;
    size_t nbEle = 0;

	data = readDoubleData_inTxt(oriFilePath, &nbEle);

    int i = 0;
    for(i=0;i<10;i++)
	    printf("data[%d]=%f\n", i, data[i]);
    if(status != RW_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
    writeDoubleData_inBytes(data, nbEle, outputFilePath, &status);
    free(data);
    
    return 0;
}
