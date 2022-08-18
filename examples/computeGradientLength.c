/**
 *  @file computeGradientLength.c
 *  @author Sheng Di
 *  @date Aug, 2022
 *  @brief This is an example of using compression interface
 *  (C) 2020 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ByteToolkit.h"
#include "rw.h"
#include "qcat.h"

int main(int argc, char * argv[])
{
	size_t r4 = 0, r3 = 0, r2 = 0, r1 = 0;	
    int status = 0;
    char dataType[4];
    char oriFilePath[640], gradLenFilePath[640];
    
    if(argc < 4)
    {
		printf("Usage: computeGradientLength [datatype (-f or -d)] [original data file (input)] [laplacian matrix file (output)]\n");
		printf("			-f means single precision; -d means double precision\n");
		printf("Example: computeGradientLength -f CLOUD_100x500x500.dat CLOUD_100x500x500.dat.gradl 500 500 100\n");
		exit(0);
    }
   
    sprintf(dataType, "%s", argv[1]);
    sprintf(oriFilePath, "%s", argv[2]);
    sprintf(gradLenFilePath, "%s", argv[3]);
    
  
    if(argc>=5)
		r1 = atoi(argv[4]); //8
    if(argc>=6)
		r2 = atoi(argv[5]); //8
    if(argc>=7)
		r3 = atoi(argv[6]); //128
    if(argc>=8)
		r4 = atoi(argv[7]);    
  
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
   
   int dType = 0;
   size_t nbEle = 0;
    if(strcmp(dataType, "-f")==0) //single precision
	{
		dType = QCAT_FLOAT;
		float *ori_data = NULL;
		printf("reading data from %s \n", oriFilePath);
		ori_data = readFloatData(oriFilePath, &nbEle, &status);
		
		float* gradLen_data = (float *)malloc(nbEle*sizeof(float));
		memset(gradLen_data, 0, nbEle*sizeof(float));
		
		printf("calcaulting....\n");
		computeGradientLength(ori_data, gradLen_data, dType, 0, r4, r3, r2, r1);
		
		//writing results
		printf("Writing resulting data to %s\n", gradLenFilePath);		
		writeData_inBytes(gradLen_data, dType, nbEle, gradLenFilePath, &status);
		
		free(ori_data);
		free(gradLen_data);			
	} 
	else
	{
		dType = QCAT_DOUBLE;
		double *ori_data = NULL;
		printf("reading data from %s \n", oriFilePath);
		ori_data = readDoubleData(oriFilePath, &nbEle, &status);

		double* gradLen_data = (double *)malloc(nbEle*sizeof(double));

		printf("calcaulting....\n");
		computeGradientLength(ori_data, gradLen_data, dType, 0, r4, r3, r2, r1);
	
		//writing results
		printf("Writing resulting data to %s\n", gradLenFilePath);
		writeData_inBytes(gradLen_data, dType, nbEle, gradLenFilePath, &status);
		
		free(ori_data);
		free(gradLen_data);	
	}
	
    return 0;
}
