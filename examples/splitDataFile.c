/**
 *  @file splitDataFile.c
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
    char outFilePath[660];
    if(argc < 2)
    {
		printf("Test case: splitDataFile [srcFilePath] [# bytes per split file]\n");
		printf("Example: splitDataFile testfloat_8_8_128.dat 4096\n");
		exit(0);
    }
   
    sprintf(oriFilePath, "%s", argv[1]);
    size_t bytesPerFile = atoll(argv[2]);
    size_t nbBytes;

    unsigned char* data = readByteData(oriFilePath, &nbBytes, &status);
    
    size_t nbFiles = nbBytes/bytesPerFile;

    int i = 0;
    for(i = 0;i<nbFiles;i++)
    {
    	size_t offsetStart = i*bytesPerFile;
	sprintf(outFilePath, "%s-%03d.bin", oriFilePath, i);

       	printf("writing data from %zu to %zu into %s\n", offsetStart, offsetStart+bytesPerFile, outFilePath);
	writeByteData(&(data[offsetStart]), bytesPerFile, outFilePath, &status);
    }

    unsigned long remainder = nbBytes%bytesPerFile;
    if(remainder!=0)
    {
    	size_t offsetStart = i*bytesPerFile;
	sprintf(outFilePath, "%s-%03d.bin", oriFilePath, i);
       	printf("writing data from %zu to %zu into %s\n", offsetStart, offsetStart+ remainder, outFilePath);
	writeByteData(&(data[offsetStart]), remainder, outFilePath, &status);
    }

    free(data);
    return 0;
}
