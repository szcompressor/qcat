/**
 *  @file generateRandomData.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ByteToolkit.h"
#include "rw.h"

int main(int argc, char * argv[])
{
    int status = 0;
    char outputFilePath[640];
    
    if(argc < 3)
    {
		printf("Test case: generateRandomData [tgtFilePath] outSize nonzeroSize\n");
		printf("Example: generateRandomData test.dat 200 4\n");
		exit(0);
    }
   
    sprintf(outputFilePath, "%s", argv[1]);
    long newSize = atol(argv[2]);
    long nonzeroSize = atol(argv[3]);
   
    size_t nbEle = newSize;
    double* data = (double*)malloc(sizeof(double)*nbEle);
    memset(data, 0, sizeof(double)*nbEle);
    int i = 0;
    for(i=0;i<nonzeroSize;i++)
	data[i] = (double)rand()/(double)(RAND_MAX);

    for(i=0;i<10;i++)
	    printf("data[%d]=%f\n", i, data[i]);
    if(status != DA_SCES)
    {
		printf("Error: data file %s cannot be read!\n", outputFilePath);
		exit(0);
    }
  
    printf("writing %ld data into %s\n", nbEle, outputFilePath);
    writeDoubleData_inBytes(data, nbEle, outputFilePath, &status);
    free(data);
    
    return 0;
}
