/**
 *  @file convertDataType.c
 *  @author Sheng Di
 *  @date April, 2021
 *  @brief This is an executable to do data type convertion.
 *  (C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include "ByteToolkit.h"
#include "rw.h"

void usage()
{
	printf("Usage: convertDataType <options>\n");
	printf("Options:\n");
	printf("*DataTypeChanging Mode:\n");
	printf("	-I : input data type (FLOAT, DOUBLE, UINT32, UINT16, INT32, INT16, ....)\n");
	printf("	-O : output data type (FLOAT, DOUBLE, UINT32, UINT16, INT32, INT16, ....)\n");
	printf("* input data file:\n");
	printf("	-i <original data file> : specify original data file\n");
	printf("* output data file:\n");
	printf("	-o <output data file> : specify output data file\n");
	printf("Examples:\n");
	printf("	convertDataFile -I UINT32 -O FLOAT -i data.f32 -o date.u32\n");	
}

int main(int argc, char * argv[])
{	
    int status = 0;
    char inputFilePath[100], outputFilePath[100];
    char inputType[100], outputType[100];
    
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
    
	if(argc==1)
	{
		usage();
		return 0;
	}

	size_t i = 0;
	for(i=1;i<argc;i++)
	{
		if (argv[i][0] != '-' || argv[i][2])
			usage();
		switch (argv[i][1])
		{
		case 'i':
			if (++i == argc)
				usage();
			strcpy(inputFilePath, argv[i]);
			break;
		case 'o':
			if (++i == argc)
				usage();
			strcpy(outputFilePath, argv[i]);	
			break;
		case 'I':
			if (++i == argc)
				usage();
			strcpy(inputType, argv[i]);
			break;
		case 'O':
			if (++i == argc)
				usage();
			strcpy(outputType, argv[i]);
			break;			
		default: 
			usage();
			break;
		}
	}

	int inType = 0, outType = 0;
	if(strcmp(inputType, "FLOAT")==0)
		inType = QCAT_FLOAT;
	else if(strcmp(inputType, "DOUBLE")==0)
		inType = QCAT_DOUBLE;
	else if(strcmp(inputType, "UINT32")==0)
		inType = QCAT_UINT32;
	else if(strcmp(inputType, "UINT16")==0)
		inType = QCAT_UINT16;
	else if(strcmp(inputType, "INT32")==0)
		inType = QCAT_INT32;
	else if(strcmp(inputType, "INT16")==0)
		inType = QCAT_INT16;
	else
	{
		printf("Error: wrong type - %s for input\n", inputType);
		exit(0);
	}
	
	if(strcmp(outputType, "FLOAT")==0)
		outType = QCAT_FLOAT;
	else if(strcmp(outputType, "DOUBLE")==0)
		outType = QCAT_DOUBLE;
	else if(strcmp(outputType, "UINT32")==0)
		outType = QCAT_UINT32;
	else if(strcmp(outputType, "UINT16")==0)
		outType = QCAT_UINT16;
	else if(strcmp(outputType, "INT32")==0)
		outType = QCAT_INT32;
	else if(strcmp(outputType, "INT16")==0)
		outType = QCAT_INT16;	
	else
	{
		printf("Error: wrong type - %s for output\n", outputType);
		exit(0);
	}	

	if(checkFileExistance(inputFilePath)==0)
	{
		printf("Error: the input file %s doesn't exist.\n", inputFilePath);
		exit(0);
	}    
    
    size_t nbEle = 0;
	if(inType == QCAT_FLOAT && outType == QCAT_DOUBLE)
	{
		float* in_data = readFloatData_systemEndian(inputFilePath, &nbEle, &status);
		double* out_data = (double*)malloc(sizeof(double)*nbEle);
		for(i = 0;i < nbEle;i++)
			out_data[i] = in_data[i];
		writeDoubleData_inBytes(out_data, nbEle, outputFilePath, &status);
		free(in_data);
		free(out_data);
	}
	else if(inType == QCAT_FLOAT && outType == QCAT_UINT32)
    {
		float* in_data = readFloatData_systemEndian(inputFilePath, &nbEle, &status);
		unsigned int* out_data = (unsigned int*)malloc(sizeof(unsigned int)*nbEle);
		for(i = 0;i < nbEle;i++)
			out_data[i] = (unsigned int)in_data[i];
		writeUIntData_inBytes(out_data, nbEle, outputFilePath, &status);
		free(in_data);
		free(out_data);		
	}
	else if(inType == QCAT_FLOAT && outType == QCAT_UINT16)
    {
		float* in_data = readFloatData_systemEndian(inputFilePath, &nbEle, &status);
		unsigned short* out_data = (unsigned short*)malloc(sizeof(unsigned short)*nbEle);
		for(i = 0;i < nbEle;i++)
			out_data[i] = (unsigned short)in_data[i];
		writeUShortData_inBytes(out_data, nbEle, outputFilePath, &status);
		free(in_data);
		free(out_data);			
	}
	else if(inType == QCAT_FLOAT && outType == QCAT_INT32)
    {
		float* in_data = readFloatData_systemEndian(inputFilePath, &nbEle, &status);
		int* out_data = (int*)malloc(sizeof(int)*nbEle);
		for(i = 0;i < nbEle;i++)
			out_data[i] = (int)in_data[i];
		writeIntData_inBytes(out_data, nbEle, outputFilePath, &status);
		free(in_data);
		free(out_data);			
	}
	else if(inType == QCAT_FLOAT && outType == QCAT_INT16)
    {
		float* in_data = readFloatData_systemEndian(inputFilePath, &nbEle, &status);
		short* out_data = (short*)malloc(sizeof(short)*nbEle);
		for(i = 0;i < nbEle;i++)
			out_data[i] = (short)in_data[i];
		writeShortData_inBytes(out_data, nbEle, outputFilePath, &status);
		free(in_data);
		free(out_data);			
	}
	else if(inType == QCAT_DOUBLE && outType == QCAT_FLOAT)
	{
		double* in_data = readDoubleData_systemEndian(inputFilePath, &nbEle, &status);
		float* out_data = (float*)malloc(sizeof(float)*nbEle);
		for(i = 0;i < nbEle;i++)
			out_data[i] = (float)in_data[i];
		writeFloatData_inBytes(out_data, nbEle, outputFilePath, &status);
		free(in_data);
		free(out_data);			
	}
	else if(inType == QCAT_UINT32 && outType == QCAT_FLOAT)
	{
		unsigned int* in_data = readUInt32Data_systemEndian(inputFilePath, &nbEle, &status);
		float* out_data = (float*)malloc(sizeof(float)*nbEle);
		for(i = 0;i < nbEle;i++)
			out_data[i] = (float)in_data[i];
		writeFloatData_inBytes(out_data, nbEle, outputFilePath, &status);
		free(in_data);
		free(out_data);			
	}
	else if(inType == QCAT_UINT16 && outType == QCAT_FLOAT)
	{
		unsigned short* in_data = readUInt16Data_systemEndian(inputFilePath, &nbEle, &status);
		float* out_data = (float*)malloc(sizeof(float)*nbEle);
		for(i = 0;i < nbEle;i++)
			out_data[i] = (float)in_data[i];
		writeFloatData_inBytes(out_data, nbEle, outputFilePath, &status);
		free(in_data);
		free(out_data);			
	}
	else if(inType == QCAT_INT32 && outType == QCAT_FLOAT)
	{
		int* in_data = readInt32Data_systemEndian(inputFilePath, &nbEle, &status);
		float* out_data = (float*)malloc(sizeof(float)*nbEle);
		for(i = 0;i < nbEle;i++)
			out_data[i] = (float)in_data[i];
		writeFloatData_inBytes(out_data, nbEle, outputFilePath, &status);
		free(in_data);
		free(out_data);			
	}
	else if(inType == QCAT_INT16 && outType == QCAT_FLOAT)
	{
		short* in_data = readInt16Data_systemEndian(inputFilePath, &nbEle, &status);
		float* out_data = (float*)malloc(sizeof(float)*nbEle);
		for(i = 0;i < nbEle;i++)
			out_data[i] = (float)in_data[i];
		writeFloatData_inBytes(out_data, nbEle, outputFilePath, &status);
		free(in_data);
		free(out_data);			
	}
       
    printf("writing data to %s\n", outputFilePath);    
    return 0;
}
