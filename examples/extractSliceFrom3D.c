/**
 *  @file extractSliceFrom3D.c
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
    char oriFilePath[640];
    char outFilePath[644];
    
    int type = 0; //float:0, double:1
    
    if(argc != 7)
    {
		printf("Test case: extractSliceFrom3D [datatype (-f or -d)] [srcFilePath] [slice number] r1 r2 r3\n");
		printf("Example: extractSliceFrom3D -f CLOUDf48.bin.dat 51 500 500 100\n");
		exit(0);
    }
   
    if(strcmp(argv[1],"-f")==0)
		type = 0;
	else if(strcmp(argv[2],"-d")==0)
		type = 1;
	else
	{
		printf("Error: missing data type\n");
		exit(0);
	}	
	
    sprintf(oriFilePath, "%s", argv[2]);
	size_t r1 = 0, r2 = 0, r3 = 0;
    int slice = atoi(argv[3]);
	r1 = atoi(argv[4]); //8
	r2 = atoi(argv[5]); //8
	r3 = atoi(argv[6]); //128
	if(slice>=r3)
	{
		printf("Error: slice is greater than r3\n");
		exit(0);
	}   
	size_t startIndex = slice*r1*r2;
	size_t length = r1*r2;
	sprintf(outFilePath, "%s.ext", oriFilePath);
    size_t nbEle;
    if(type==0)
    {
		float *data = readFloatData(oriFilePath, &nbEle, &status);
		if(status != RW_SCES)
		{
			printf("Error: data file %s cannot be read!\n", oriFilePath);
			exit(0);
		}
		writeFloatData_inBytes(&(data[startIndex]), length, outFilePath, &status);		
		free(data);		
	}
	else if(type==1)
	{
		double *data = readDoubleData(oriFilePath, &nbEle, &status);
		if(status != RW_SCES)
		{
			printf("Error: data file %s cannot be read!\n", oriFilePath);
			exit(0);
		}

		writeDoubleData_inBytes(&(data[startIndex]), length, outFilePath, &status);		
		free(data);
	}

    printf("outputFile= %s\n", outFilePath);
    
    return 0;
}
