/**
 *  @file convertEndianType.c
 *  @author Sheng Di
 *  @date April, 2021
 *  @brief This is an example of using compression interface
 *  (C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include "ByteToolkit.h"
#include "rw.h"

int main(int argc, char * argv[])
{
    int status = 0;
    char datatype[6]; 
    char oriFilePath[640], outputFilePath[640];
    
    if(argc < 4)
    {
		printf("Usage: convertEndianType -f/-d input output (-f means single-precision; -d means double-precision)\n");
		printf("Example: convertEndianType -f CLOUD_48_bigendian.dat CLOUD_48_littleendian.dat \n");
		exit(0);
    }
   
    sprintf(datatype, "%s", argv[1]);
    sprintf(oriFilePath, "%s", argv[2]);
    sprintf(outputFilePath, "%s", argv[3]);
   
   
    unsigned char *byteData = NULL;
    size_t nbBytes = 0;
    printf("reading data from %s \n", oriFilePath);
    byteData = readByteData(oriFilePath, &nbBytes, &status);

    printf("converting...\n");

    size_t i = 0;
    unsigned char* p = byteData;

    if(strcmp(datatype, "-f")==0)
    {
    	for(i=0;i<nbBytes;i+=4) //only for float
    	{
    		symTransform_4bytes(p);
		p += 4;
    	}
    }
    else //-d
    {
    	for(i=0;i<nbBytes;i+=8)
	{
		symTransform_8bytes(p);
		p += 8;
	}
    }

    printf("writing data to %s\n", outputFilePath);
    writeByteData(byteData, nbBytes, outputFilePath, &status);
    free(byteData);
    
    return 0;
}
