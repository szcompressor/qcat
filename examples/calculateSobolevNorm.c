/**
 *  @file calculateSobolevNorm.c
 *  @author Sheng Di
 *  @date Sept, 2022
 *  @brief This is an example of using compression interface
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
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
    char dataFilePath[640];
    
    if(argc < 3)
    {
		printf("Usage: calculateSobolevNorm [datatype (-f or -d)] [data file] [order(0,1,2,-1)] [dimesions... (from fast to slow)]\n");
		printf("			-f means single precision; -d means double precision\n");
		printf("			order: -1 means doing the calculation for all orders (0, 1, and 2)\n");
		printf("			The type of norm is always 2-norm, i.e., p=2\n");
		printf("Example: calculateSobolevNorm -f CLOUD_100x500x500.dat -1 500 500 100\n");
		exit(0);
    }
   
    sprintf(dataType, "%s", argv[1]);
    sprintf(dataFilePath, "%s", argv[2]);
    int order = atoi(argv[3]);
  
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
   
    double result = 0;
    if(strcmp(dataType, "-f")==0) //single precision
	{
		int dType = QCAT_FLOAT;
		float *data = NULL;
		size_t nbEle = 0;
		printf("reading data from %s \n", dataFilePath);
		data = readFloatData(dataFilePath, &nbEle, &status);
		
		printf("calcaulting....\n");
		if(order==0)
		{
			result = calculateSobolevNorm_p2(data, dType, 0, 0, r4, r3, r2, r1);
			printf("sobolev norm (order=%d): %.10G\n", 0,result);
		}
		else if(order==1)
		{
			result = calculateSobolevNorm_p2(data, dType, 1, 0, r4, r3, r2, r1);
			printf("sobolev norm (order=%d): %.10G\n", 1,result);			
		}
		else if(order==2)
		{
			result = calculateSobolevNorm_p2(data, dType, 2, 0, r4, r3, r2, r1);
			printf("sobolev norm (order=%d): %.10G\n", 2,result);				
		}
		else if(order==-1)
		{
			result = calculateSobolevNorm_p2(data, dType, 0, 0, r4, r3, r2, r1);
			printf("sobolev norm (order=%d): %.10G\n", 0,result);	
			result = calculateSobolevNorm_p2(data, dType, 1, 0, r4, r3, r2, r1);
			printf("sobolev norm (order=%d): %.10G\n", 1,result);	
			result = calculateSobolevNorm_p2(data, dType, 2, 0, r4, r3, r2, r1);
			printf("sobolev norm (order=%d): %.10G\n", 2,result);											
		}
		
		free(data);
	} 
	else
	{
		int dType = QCAT_DOUBLE;
		double *data = NULL;
		size_t nbEle = 0;
		printf("reading data from %s \n", dataFilePath);
		data = readDoubleData(dataFilePath, &nbEle, &status);
		
		printf("calcaulting....\n");
		if(order==0)
		{
			result = calculateSobolevNorm_p2(data, dType, 0, 0, r4, r3, r2, r1);
			printf("sobolev norm (order=%d): %.10G\n", 0,result);
		}
		else if(order==1)
		{
			result = calculateSobolevNorm_p2(data, dType, 1, 0, r4, r3, r2, r1);
			printf("sobolev norm (order=%d): %.10G\n", 1,result);			
		}
		else if(order==2)
		{
			result = calculateSobolevNorm_p2(data, dType, 2, 0, r4, r3, r2, r1);
			printf("sobolev norm (order=%d): %.10G\n", 2,result);				
		}
		else if(order==-1)
		{
			result = calculateSobolevNorm_p2(data, dType, 0, 0, r4, r3, r2, r1);
			printf("sobolev norm (order=%d): %.10G\n", 0,result);	
			result = calculateSobolevNorm_p2(data, dType, 1, 0, r4, r3, r2, r1);
			printf("sobolev norm (order=%d): %.10G\n", 1,result);	
			result = calculateSobolevNorm_p2(data, dType, 2, 0, r4, r3, r2, r1);
			printf("sobolev norm (order=%d): %.10G\n", 2,result);											
		}
		
		free(data);
	}

    return 0;
}
