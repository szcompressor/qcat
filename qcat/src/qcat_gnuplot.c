/**
 *  @file qcat_gnuplot.c
 *  @author Sheng Di
 *  @date July, 2019
 *  @brief gnuplot related operations
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "qcat_gnuplot.h"
#include <stdio.h> 
#include <math.h>
#include <string.h> 
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>

char** genGnuplotScript_sliceImage(char* dataFileName, size_t r2, size_t r1, char* imageFileName)
{
	char** lines = (char**)malloc(10*sizeof(char*));
	size_t s1, s2;
	s1 = r1<900?900:r1;
	s2 = r2<900?900:r2;
	
	int i = 0;
	for(i=0;i<10;i++)
	{
		lines[i] = (char*)malloc(250);
		memset(lines[i], 0, 250);
	}
	s1=900;
	s2=900;
	sprintf(lines[0], "#!/usr/bin/gnuplot\n");
	sprintf(lines[1], "set term png size %zu, %zu enhanced font \"Arial,24\"\n", s1, s2);
	sprintf(lines[2], "set pm3d map\n");
	sprintf(lines[3], "#set cbrange [-4:4]\n");
	sprintf(lines[4], "set output \"%s\"\n", imageFileName);
	sprintf(lines[5], "#set size square\n");
	sprintf(lines[6], "set xrange [0:%zu]\n", r1);
	sprintf(lines[7], "set yrange [%zu:0]\n", r2);
	sprintf(lines[8], "set palette rgbformulae 33,13,10\n");
	sprintf(lines[9], "splot \"%s\"", dataFileName);

	return lines;
}

float* computeSlice2DLog_float(size_t r2, size_t r1, float* data)
{
	float* logSliceData = (float*)malloc(sizeof(float)*r2*r1);
	size_t i = 0, j = 0;
	for(i=0;i<r2;i++)
		for(j=0;j<r1;j++)
		{
			size_t index = i*r1+j;
			logSliceData[index] = log10f(fabsf(data[index]));
		}
	return logSliceData;
}

float* computeSlice3DLog_float(int sliceNumber, size_t r3, size_t r2, size_t r1, float* data)
{
	float* oriSliceData = &data[sliceNumber*r2*r1];
	float* logSliceData = (float*)malloc(sizeof(float)*r2*r1);
	size_t i = 0, j = 0;
	for(i=0;i<r2;i++)
		for(j=0;j<r1;j++)
		{
			size_t index = i*r1+j;
			logSliceData[index] = log10f(fabsf(oriSliceData[index]));
		}
	return logSliceData;
}

double* computeSlice2DLog_double(size_t r2, size_t r1, double* data)
{
	double* logSliceData = (double*)malloc(sizeof(double)*r2*r1);
	size_t i = 0, j = 0;
	for(i=0;i<r2;i++)
		for(j=0;j<r1;j++)
		{
			size_t index = i*r1+j;
			logSliceData[index] = log10(fabs(data[index]));
		}
	return logSliceData;
}

double* computeSlice3DLog_double(int sliceNumber, size_t r3, size_t r2, size_t r1, double* data)
{
	double* oriSliceData = &data[sliceNumber*r2*r1];
	double* logSliceData = (double*)malloc(sizeof(double)*r2*r1);
	size_t i = 0, j = 0;
	for(i=0;i<r2;i++)
		for(j=0;j<r1;j++)
		{
			size_t index = i*r1+j;
			logSliceData[index] = log10(fabs(oriSliceData[index]));
		}
	return logSliceData;
}

float* transformData_float(int plotDim, size_t r3, size_t r2, size_t r1, float* oriData)
{
	float* result = (float*)malloc(r3*r2*r1*sizeof(float));
	
	size_t i, j, k;
	
	switch(plotDim)
	{
		case 3:
			memcpy(result, oriData, r3*r2*r1*sizeof(float));
			break;		
		case 2:
		for(k=0;k<r3;k++)
			for(j=0;j<r2;j++)
				for(i=0;i<r1;i++)
				{
					result[(r2-1-j)*r1*r3+k*r1+i] = oriData[k*r2*r1+j*r1+i];
				}
		break;
		case 1:
		for(k=0;k<r3;k++)
			for(j=0;j<r2;j++)
				for(i=0;i<r1;i++)
				{
					result[(r1-1-i)*r2*r3+k*r2+j] = oriData[k*r2*r1+j*r1+i];
				}
		break;			
	}
	
	return result;

}

double* transformData_double(int plotDim, size_t r3, size_t r2, size_t r1, double* oriData)
{
	double* result = (double*)malloc(r3*r2*r1*sizeof(double));
	
	size_t i, j, k;
	
	switch(plotDim)
	{
		case 3:
			memcpy(result, oriData, r3*r2*r1*sizeof(double));	
			break;	
		case 2:
		for(k=0;k<r3;k++)
			for(j=0;j<r2;j++)
				for(i=0;i<r1;i++)
				{
					result[(r2-1-j)*r1*r3+k*r1+i] = oriData[k*r2*r1+j*r1+i];
				}
		break;
		case 1:
		for(k=0;k<r3;k++)
			for(j=0;j<r2;j++)
				for(i=0;i<r1;i++)
				{
					result[(r1-1-i)*r2*r3+k*r2+j] = oriData[k*r2*r1+j*r1+i];
				}
		break;		
	}
	
	return result;

}

float* generateSliceData_float(int plotDim, int sliceNumber, size_t r3, size_t r2, size_t r1, float* data, int domain)
{	
	float* sliceData = NULL;
	if(r3 == 0)
	{
		if(domain==1)	
			sliceData = computeSlice2DLog_float(r2, r1, data);
		else
		{
			sliceData = (float*)malloc(sizeof(float)*r2*r1);
			memcpy(sliceData, data, sizeof(float)*r2*r1);			
		}
	}
	else //r3 >=1
	{
		//transform the data in place
		float* tdata = transformData_float(plotDim, r3, r2, r1, data);	
			
		if(domain == 1)
		{
			if(plotDim==3)
				sliceData = computeSlice3DLog_float(sliceNumber, r3, r2, r1, tdata);
			else if(plotDim==2)
				sliceData = computeSlice3DLog_float(sliceNumber, r2, r3, r1, tdata);
			else 
				sliceData = computeSlice3DLog_float(sliceNumber, r1, r3, r2, tdata);
		}
		else if(domain == 0)
		{
			size_t ss = 0;
			switch(plotDim)
			{
			case 3:
				ss = r1*r2;
				break;
			case 2:
				ss = r1*r3;
				break;
			case 1:
				ss = r2*r3;
				break;
			}
			sliceData = (float*)malloc(sizeof(float)*ss);
			memcpy(sliceData, &tdata[sliceNumber*ss], sizeof(float)*ss);		
		}
		free(tdata);		
	}
	return sliceData;
}

double* generateSliceData_double(int plotDim, int sliceNumber, size_t r3, size_t r2, size_t r1, double* data, int domain)
{
	double* sliceData = NULL;
	if(r3 == 0)
	{
		if(domain==1)	
			sliceData = computeSlice2DLog_double(r2, r1, data);
		else
		{
			sliceData = (double*)malloc(sizeof(double)*r2*r1);
			memcpy(sliceData, data, sizeof(double)*r2*r1);			
		}
	}
	else //r3 >=1
	{
		//transform the data in place
		double* tdata = transformData_double(plotDim, r3, r2, r1, data);	
			
		if(domain == 1)
		{
			if(plotDim==3)
				sliceData = computeSlice3DLog_double(sliceNumber, r3, r2, r1, tdata);
			else if(plotDim==2)
				sliceData = computeSlice3DLog_double(sliceNumber, r2, r3, r1, tdata);
			else 
				sliceData = computeSlice3DLog_double(sliceNumber, r1, r3, r2, tdata);
		}
		else if(domain == 0)
		{
			size_t ss = 0;
			switch(plotDim)
			{
			case 3:
				ss = r1*r2;
				break;
			case 2:
				ss = r1*r3;
				break;
			case 1:
				ss = r2*r3;
				break;
			}
			sliceData = (double*)malloc(sizeof(double)*ss);
			memcpy(sliceData, &tdata[sliceNumber*ss], sizeof(double)*ss);		
		}
		free(tdata);		
	}
	return sliceData;
}

float* generateSliceDiff_float(int plotDim, int sliceNumber, size_t r3, size_t r2, size_t r1, float* oriData, float* decData, int domain)
{
	size_t i = 0, j = 0, index = 0;
	
	float* sliceOriData = NULL;
	float* sliceDecData = NULL;
	float* t_oriData = NULL;
	float* t_decData = NULL;
	
	if(r3 == 0)
	{
		if(domain == 1)
		{
			sliceOriData = computeSlice2DLog_float(r2, r1, oriData);
			sliceDecData = computeSlice2DLog_float(r2, r1, decData);
		}
		else
		{
			sliceOriData = oriData;
			sliceDecData = decData;			
		}
	}
	else
	{
		t_oriData = transformData_float(plotDim, r3, r2, r1, oriData);	
		t_decData = transformData_float(plotDim, r3, r2, r1, decData);
		if(domain == 1)
		{
			sliceOriData = computeSlice3DLog_float(sliceNumber, r3, r2, r1, t_oriData);
			sliceDecData = computeSlice3DLog_float(sliceNumber, r3, r2, r1, t_decData);				
		}
		else
		{
			sliceOriData = &t_oriData[sliceNumber*r2*r1];
			sliceDecData = &t_decData[sliceNumber*r2*r1];
		}
	}

	float* sliceDiff = (float*)malloc(sizeof(float)*r2*r1);	
	for(i=0;i<r2;i++)
		for(j=0;j<r1;j++)
		{
			index = i*r1+j;
			sliceDiff[index] = sliceOriData[index] - sliceDecData[index];
		}
	
	if(domain==1)
	{
		free(sliceOriData);
		free(sliceDecData);
	}
	
	if(t_oriData!=NULL)
		free(t_oriData);
	if(t_decData!=NULL)
		free(t_decData);
	
	return sliceDiff;
}

double* generateSliceDiff_double(int plotDim, int sliceNumber, size_t r3, size_t r2, size_t r1, double* oriData, double* decData, int domain)
{
	size_t i = 0, j = 0, index = 0;
	
	double* sliceOriData = NULL;
	double* sliceDecData = NULL;
	double* t_oriData = NULL;
	double* t_decData = NULL;
	
	if(r3 == 0)
	{
		if(domain == 1)
		{
			sliceOriData = computeSlice2DLog_double(r2, r1, oriData);
			sliceDecData = computeSlice2DLog_double(r2, r1, decData);
		}
		else
		{
			sliceOriData = oriData;
			sliceDecData = decData;			
		}
	}
	else
	{
		t_oriData = transformData_double(plotDim, r3, r2, r1, oriData);	
		t_decData = transformData_double(plotDim, r3, r2, r1, decData);
		if(domain == 1)
		{
			sliceOriData = computeSlice3DLog_double(sliceNumber, r3, r2, r1, t_oriData);
			sliceDecData = computeSlice3DLog_double(sliceNumber, r3, r2, r1, t_decData);				
		}
		else
		{
			sliceOriData = &t_oriData[sliceNumber*r2*r1];
			sliceDecData = &t_decData[sliceNumber*r2*r1];
		}
	}

	double* sliceDiff = (double*)malloc(sizeof(double)*r2*r1);	
	for(i=0;i<r2;i++)
		for(j=0;j<r1;j++)
		{
			index = i*r1+j;
			sliceDiff[index] = sliceOriData[index] - sliceDecData[index];
		}
	
	if(domain==1)
	{
		free(sliceOriData);
		free(sliceDecData);
	}
	
	if(t_oriData!=NULL)
		free(t_oriData);
	if(t_decData!=NULL)
		free(t_decData);
	
	return sliceDiff;
}
