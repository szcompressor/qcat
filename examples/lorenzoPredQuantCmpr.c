/**
 *  @file lorenzoPredQuant.c
 *  @author Sheng Di
 *  @date April, 2022
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include "rw.h"
#include "sz_utility.h"

#include <sys/time.h>      /* For gettimeofday(), in microseconds */
#include <time.h>          /* For time(), in seconds */

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


void cost_start()
{
    totalCost = 0;
    gettimeofday(&costStart, NULL);
}

void cost_end()
{
    double elapsed;
    struct timeval costEnd;
    gettimeofday(&costEnd, NULL);
    elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
    totalCost += elapsed;
}


int main(int argc, char * argv[])
{
    size_t r1 = 0, r2 = 0, r3 = 0;
    int status = 0;
    char oriFilePath[640], outFilePath[645];
    char lmode[30];
    char type[3];
    if(argc < 3)
    {
	printf("Test case: lorenzoPredQuantCmpr [type(-f/-d)] [dataFilePath] [lorenzo_mode] [error_bound] [dims...]\n");
	printf("Example: lorenzoPredQuantCmpr -f Hurricane.dat LORENZO_1D_1LAYER 1E-2 500 500 100\n");
	exit(0);
    }
  
    sprintf(type, "%s", argv[1]);
    sprintf(oriFilePath, "%s", argv[2]);
    sprintf(outFilePath, "%s.i32", oriFilePath);

    sprintf(lmode, "%s", argv[3]);
    double errorBound = atof(argv[4]);

    if(argc>=6)
	r1 = atoi(argv[5]); //8
    if(argc>=7)
	r2 = atoi(argv[6]); //8
    if(argc>=8)
	r3 = atoi(argv[7]); //128
    
    int mode = 0;
    if(strcmp(lmode, "LORENZO_1D_1LAYER")==0)
	mode = LORENZO_1D_1LAYER;
    else if(strcmp(lmode, "LORENZO_1D_2LAYER")==0)
	mode = LORENZO_1D_2LAYER;
    else if(strcmp(lmode, "LORENZO_1D_3LAYER")==0)
	mode = LORENZO_1D_3LAYER;
    else if(strcmp(lmode, "LORENZO_2D_1LAYER")==0)
	mode = LORENZO_2D_1LAYER;
    else if(strcmp(lmode, "LORENZO_3D_1LAYER")==0)
	mode = LORENZO_3D_1LAYER;
    else
    {
	printf("Error: wrong lorenzo_mode\n");
	exit(0);
    }

    int dataType = 0; 
    size_t nbEle = computeDataLength(0, 0, r3, r2, r1);
    int* out = (int*)malloc(sizeof(int)*nbEle);
    if(strcmp(type, "-f")==0)
    {
	dataType = QCAT_FLOAT;
	float *data = readFloatData(oriFilePath, &nbEle, &status);
	cost_start();
	status = lorenzoPredictorQuant_Cmpr_NoOutlier(data, dataType, mode, errorBound, r3, r2, r1, out);
	cost_end();
	printf("total time cost = %f\n", totalCost);
    }
    else
    {
	dataType = QCAT_DOUBLE;
	double *data = readDoubleData(oriFilePath, &nbEle, &status);
	cost_start();
	status = lorenzoPredictorQuant_Cmpr_NoOutlier(data, dataType, mode, errorBound, r3, r2, r1, out);
	cost_end();
	printf("total time cost = %f\n", totalCost);
    }
   
    if(status!=0)
    {
    	printf("Error state returned by lorenzoPredictorQuant function.\n");
	exit(0);
    } 
    writeIntData_inBytes(out, nbEle, outFilePath, &status);
    printf("quantization codes are stored in %s\n", outFilePath);

    return status;
}
