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
#include "qcat_gnuplot.h"

int main(int argc, char * argv[])
{
    int status = 0;
    char dataType[4];
    char oriFilePath[640];
    
    if(argc < 3)
    {
		printf("Usage: computeDataPDF [datatype (-i] [data file]\n");
		printf("Example 1: computeDataPDF -i FREQSH_55_1e-3-lorenzo.q\n");		
		exit(0);
    }
   
    sprintf(dataType, "%s", argv[1]);
    sprintf(oriFilePath, "%s", argv[2]);
  
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
    
	double* pdf = NULL;
	double min = 0;
	int intervals = 0;
	
    if(strcmp(dataType, "-i")==0)
    {
		int* data;
		size_t nbEle = 0;
		printf("reading data from %s \n", oriFilePath);
		data = readInt32Data(oriFilePath, &nbEle, &status);
		printf("calcaulting....\n");
		int intMin = 0;
		pdf = computeDataPDF_int32(QCAT_INT32, data, nbEle, &intMin, &intervals);
		min = intMin;	
		
		free(data);
	}
	
	//write pdf to a file
	char* tgtPDFDataPath = (char*)malloc(QCAT_BUFS_LONG);
	sprintf(tgtPDFDataPath, "%s.dis", oriFilePath);
	writePDFData_int32(tgtPDFDataPath, min, intervals, pdf);
	free(tgtPDFDataPath);

	//generate gnuplot script file
	char* tgtPDFScriptPath = (char*)malloc(QCAT_BUFS_LONG);
	sprintf(tgtPDFScriptPath, "%s.p", oriFilePath);
	char** gnuplotScriptLines = genGnuplotScript_fillsteps(oriFilePath, "dis", 2, "Data", "Probability Density Function (PDF)");
	RW_writeStrings(19, gnuplotScriptLines, tgtPDFScriptPath);	
	
	printf("Done.\n");
	printf("You need to run following lines to generate distrituion img files:\n");
	char* dir = extractDirFromPath(tgtPDFScriptPath);
	char* filename = extractFileNameFromPath(tgtPDFScriptPath);
	printf("cd %s; gnuplot %s\n", dir, filename);
	free(tgtPDFScriptPath);	
	free(pdf);	
	
    return 0;
}
