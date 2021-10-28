/**
 *  @file compute auto correlation for compression error
 *  @author Sheng Di
 *  @date Nov., 2021
 *  @brief 
 *  (C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ByteToolkit.h"
#include "rw.h"
#include "qcat.h"
#include "qcat_gnuplot.h"

int main(int argc, char * argv[])
{
    int status = 0;
    int dataType_ = 0;
    char dataType[4];
    char oriFilePath[640];
    char decFilePath[640];
    char outputFilePath[640];

    if(argc < 3)
    {
	printf("Usage: computeErrAutoCorrelation [datatype (-f/-d)] (optional: -o [data file])\n");
	printf("Example 1: computeErrAutoCorrelation -f test.f32 dec.f32 \n");
	printf("Example 2: computeErrAutoCorrelation -f test.f32 dec.f32 -o test-err.f32.ac\n");	
	exit(0);
    }
   
    sprintf(dataType, "%s", argv[1]);
    sprintf(oriFilePath, "%s", argv[2]);
    sprintf(decFilePath, "%s", argv[3]);
    if(argc > 5)
    	sprintf(outputFilePath, "%s", argv[5]);
  
    int x = 1;
    char *y = (char*)&x;

    if(*y==1)
    {     
	sysEndianType = 0; //LITTLE_ENDIAN_SYSTEM;
	//printf("This is little-endian system.\n");
    }
    else //=0
    {
         sysEndianType = 1; //BIG_ENDIAN_SYSTEM;
         //printf("This is big-endian system.\n");
    }
    
    double* acEffs = NULL;

    size_t i = 0, nbEle = 0;

    if(strcmp(dataType, "-f")==0)
    {
	dataType_ = QCAT_FLOAT;
	float *data = readFloatData(oriFilePath, &nbEle, &status);
	float *dec = readFloatData(decFilePath, &nbEle, &status);
	float *diff = (float*)malloc(sizeof(float)*nbEle);
	for(i=0;i<nbEle;i++)
		diff[i] = data[i] - dec[i];
	acEffs = ZC_compute_autocorrelation1D(diff, dataType_, nbEle);
	free(data);
	free(dec);
	free(diff);
    }
    else if(strcmp(dataType, "-d")==0)
    {
	dataType_ = QCAT_DOUBLE;
	double *data = readDoubleData(oriFilePath, &nbEle, &status);
	double *dec = readDoubleData(decFilePath, &nbEle, &status);
	double *diff = (double*)malloc(sizeof(double)*nbEle);
	for(i=0;i<nbEle;i++)
		diff[i] = data[i] - dec[i];
	acEffs = ZC_compute_autocorrelation1D(diff, dataType_, nbEle);	
	free(data);
	free(dec);
	free(diff);
    }
    else
    {
	printf("Error: wrong data type!\n");
	exit(0);
    }
    printf("Lag-1 auto correlation: %.20G\n", acEffs[1]);

    if(argc > 4)
    {
    	writeDoubleData(acEffs, AUTOCORR_SIZE+1, outputFilePath, &status);

    	printf("The full result has been written to %s\n", outputFilePath);
    }
	
    return 0;
}
