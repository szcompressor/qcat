/**
 *  @file calculateSSIM.c
 *  @author Sheng Di
 *  @date Sept, 2021
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
    char oriFilePath[640], decFilePath[640];
    
    if(argc < 3)
    {
		printf("Usage: calculateSSIM [datatype (-f or -d)] [original data file] [decompressed data file] [dimesions... (from fast to slow)]\n");
		printf("			-f means single precision; -d means double precision\n");
		printf("Example: calculateSSIM -f CLOUD_100x500x500.dat CLOUD_100x500x500.dat.sz.out 500 500 100\n");
		exit(0);
    }
   
    sprintf(dataType, "%s", argv[1]);
    sprintf(oriFilePath, "%s", argv[2]);
    sprintf(decFilePath, "%s", argv[3]);
  
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
   
   
    if(strcmp(dataType, "-f")==0) //single precision
	{
		float *ori_data = NULL, *dec_data = NULL;
		size_t nbEle = 0, nbEle2 = 0;
		printf("reading data from %s \n", oriFilePath);
		ori_data = readFloatData(oriFilePath, &nbEle, &status);
		
		dec_data = readFloatData(decFilePath, &nbEle2, &status);
		
		if(nbEle != nbEle2)
		{
			printf("Error: number of elements is not consistent\n");
			exit(0);
		}

		printf("calcaulting....\n");
		double ssim = calculateSSIM(ori_data, dec_data, QCAT_FLOAT, r4, r3, r2, r1);
		printf("ssim = %f\n", ssim);
		
		free(ori_data);
		free(dec_data);		
	} 
	else
	{
		double *ori_data = NULL, *dec_data = NULL;
		size_t nbEle = 0, nbEle2 = 0;
		printf("reading data from %s \n", oriFilePath);
		ori_data = readDoubleData(oriFilePath, &nbEle, &status);
		
		dec_data = readDoubleData(decFilePath, &nbEle2, &status);
				
		if(nbEle != nbEle2)
		{
			printf("Error: number of elements is not consistent\n");
			exit(0);
		}		

		printf("calcaulting....\n");
		double ssim = calculateSSIM(ori_data, dec_data, QCAT_DOUBLE, r4, r3, r2, r1);
		printf("ssim = %f\n", ssim);
		
		free(ori_data);
		free(dec_data);
	}

    return 0;
}
