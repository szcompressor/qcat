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
    char outputFilePath[640];

    if(argc < 3)
    {
	printf("Usage: computeErrAutoCorrelation [datatype (-f/-d)] (optional: -o [data file])\n");
	printf("Example 1: computeAutoCorrelation -f test.f32\n");
	printf("Example 2: computeAutoCorrelation -f test.f32 -o test.f32.ac\n");	
	exit(0);
    }
   
    sprintf(dataType, "%s", argv[1]);
    sprintf(oriFilePath, "%s", argv[2]);
    if(argc > 4)
    	sprintf(outputFilePath, "%s", argv[4]);
  
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
	acEffs = ZC_compute_autocorrelation1D(data, dataType_, nbEle);
	free(data);
    }
    else if(strcmp(dataType, "-d")==0)
    {
	dataType_ = QCAT_DOUBLE;
	double *data = readDoubleData(oriFilePath, &nbEle, &status);
	acEffs = ZC_compute_autocorrelation1D(data, dataType_, nbEle);	
	free(data);
    }
    else
    {
	printf("Error: wrong data type!\n");
	exit(0);
    }
    for(i = 0;i<512;i++)
    	printf("Lag-%zu auto correlation: %.20G\n", i, acEffs[i]);

    if(argc > 4)
    {
    	writeDoubleData(acEffs, AUTOCORR_SIZE+1, outputFilePath, &status);

    	printf("The full result has been written to %s\n", outputFilePath);
    }

    free(acEffs); 
    return 0;
}
