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
    char oriFilePath[640], decFilePath[640];
    
    if(argc < 4)
    {
		printf("Usage: computeErrPDF [datatype (-f or -d)] [fix_interval] [original data file] [decompressed data file]\n");
		printf("			-f means single precision; -d means double precision\n");
		printf("			fix_interval: a positive number means fix_interval, -1 means fix_number_of_intevals\n");
		printf("Example: computeErrPDF -f -1 CLOUD_100x500x500.dat CLOUD_100x500x500.dat.sz.out\n");
		printf("Example: computeErrPDF -f 0.001 CLOUD_100x500x500.dat CLOUD_100x500x500.dat.sz.out\n");
		exit(0);
    }
   
    sprintf(dataType, "%s", argv[1]);
    double fix_interval = atof(argv[2]);
    sprintf(oriFilePath, "%s", argv[3]);
    sprintf(decFilePath, "%s", argv[4]);
  
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
	double min_diff = 0;
	double err_interval = 0;
	int intervals = 0;
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
		pdf = computeErrPDF(QCAT_FLOAT, ori_data, dec_data, nbEle, fix_interval, &min_diff, &err_interval, &intervals);	
		
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
		pdf = computeErrPDF(QCAT_DOUBLE, ori_data, dec_data, nbEle, fix_interval, &min_diff, &err_interval, &intervals);
		
		free(ori_data);
		free(dec_data);
	}
	
	//write pdf to a file
	char* tgtPDFDataPath = (char*)malloc(QCAT_BUFS_LONG);
	sprintf(tgtPDFDataPath, "%s.dis", decFilePath);
	writePDFData_error(tgtPDFDataPath, min_diff, err_interval, intervals, pdf);
	free(tgtPDFDataPath);

	//generate gnuplot script file
	char* tgtPDFScriptPath = (char*)malloc(QCAT_BUFS_LONG);
	sprintf(tgtPDFScriptPath, "%s.p", decFilePath);
	char** gnuplotScriptLines = genGnuplotScript_lines(decFilePath, "dis", 2, "Compression Error", "Probability Density Function (PDF)");
	RW_writeStrings(24, gnuplotScriptLines, tgtPDFScriptPath);	
	
	printf("Done.\n");
	printf("You need to run following lines to generate distrituion img files:\n");
	char* dir = extractDirFromPath(tgtPDFScriptPath);
	char* filename = extractFileNameFromPath(tgtPDFScriptPath);
	printf("cd %s; gnuplot %s\n", dir, filename);
	free(tgtPDFScriptPath);	
	free(pdf);	
	
    return 0;
}
