/**
 *  @file qcat_dataAnalysis.c
 *  @author Sheng Di
 *  @date Feb, 2021
 *  @brief data anlaysis
 *  (C) 2015-2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "qcat_dataAnalysis.h"
#include <stdio.h> 
#include <math.h>
#include <string.h> 
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <rw.h>
#include <qcat_hashtable.h>

double computeLosslessEntropy_8bits(int dataType, void* data, size_t nbEle)
{
	size_t i = 0;
	unsigned char* bytes = (unsigned char*)data;
	size_t totalLen = dataType == QCAT_FLOAT? nbEle*sizeof(float): nbEle*sizeof(double);

	double entVal = 0.0;
	unsigned char index = 0;
	size_t table_size = 256;
	long *table = (long*)malloc(table_size*sizeof(long));
	memset(table, 0, table_size*sizeof(long));
				
	for(i=0;i<totalLen;i++)
	{
		index = bytes[i];
		table[index]++;
	}
		
	size_t sum = dataType==QCAT_FLOAT?nbEle*sizeof(float):nbEle*sizeof(double);
	for (i = 0; i<table_size; i++)
		if (table[i] != 0)
		{
			double prob = (double)table[i]/sum;
			entVal -= prob*log(prob)/log(2);
		}

	free(table);
	return entVal;
}

double computeLosslessEntropy_32bits(void* data, size_t nbEle)
{
	size_t i = 0;
	char vs[256];	
	
	hashtable_t* entropy_table = ht_create(nbEle);
	float* value = (float*)data;
	for(i=0;i<nbEle;i++)
	{
        	lfloat buf;
		buf.value = value[i];
 		unsigned int v = buf.ivalue;
		
		sprintf(vs, "%d", v);	
		QCAT_ELEMENT* qe = ht_get(entropy_table, vs);
		if(qe==NULL)
		{
			qe = (QCAT_ELEMENT*) malloc(sizeof(QCAT_ELEMENT));
			memset(qe, 0, sizeof(QCAT_ELEMENT));
			qe->value = value[i];
			ht_set(entropy_table, vs, qe);
		}
		qe->counter ++;
	}		
	
	size_t sum = nbEle;
	size_t j = 0;
	double entVal = 0;
	for(i=0;i<entropy_table->capacity&&j<entropy_table->count;i++)
	{
		entry_t* t = entropy_table->table[i];
		while(t!=NULL)
		{
			QCAT_ELEMENT* qe = (QCAT_ELEMENT*)t->value;
			double prob = ((double)qe->counter)/sum;
			entVal -= prob*log(prob)/log(2);
			free(qe);
			t = t->next;
		}
	}	
	
	ht_freeTable(entropy_table);
	return entVal;
}

double computeLosslessEntropy_64bits(void* data, size_t nbEle)
{
	char vs[256];
	size_t i = 0;
	hashtable_t* entropy_table = ht_create(nbEle);
	double* value = (double*)data;	
	for(i=0;i<nbEle;i++)
	{
        ldouble buf;
        buf.value = value[i];
        unsigned long v = buf.lvalue;
		
		sprintf(vs, "%lu", v);	
		QCAT_ELEMENT* qe = ht_get(entropy_table, vs);
		if(qe==NULL)
		{
			qe = (QCAT_ELEMENT*) malloc(sizeof(QCAT_ELEMENT));
			memset(qe, 0, sizeof(QCAT_ELEMENT));
			qe->value = value[i];
			ht_set(entropy_table, vs, qe);
		}
		qe->counter ++;
	}
	
	size_t sum = nbEle;
	size_t j = 0;
	double entVal = 0;
	for(i=0;i<entropy_table->capacity&&j<entropy_table->count;i++)
	{
		entry_t* t = entropy_table->table[i];
		while(t!=NULL)
		{
			QCAT_ELEMENT* qe = (QCAT_ELEMENT*)t->value;
			double prob = ((double)qe->counter)/sum;
			entVal -= prob*log(prob)/log(2);
			free(qe);
			t = t->next;
		}
	}	
	
	ht_freeTable(entropy_table);
	return entVal;	
}

QCAT_DataProperty* computeProperty(int dataType, void* data, size_t nbEle, int entropyType)
{
	QCAT_DataProperty* property = (QCAT_DataProperty*)malloc(sizeof(QCAT_DataProperty));
	memset(property, 0, sizeof(QCAT_DataProperty));
	
	property->dataType = dataType;
	property->numOfElem = nbEle;
	property->entropy_8bits = -1;
	property->entropy_32bits = -1;
	property->entropy_64bits = -1;

	size_t i = 0;
	if(dataType == QCAT_FLOAT)
	{
		float* data_ = (float*)data;
		double min=data_[0],max=data_[0],sum=0;		
		for(i=0;i<nbEle;i++)
		{
			if(min>data_[i]) min = data_[i];
			if(max<data_[i]) max = data_[i];
			sum += data_[i];
		}

		double med = min+(max-min)/2;
		double sum_of_square = 0;
		for(i=0;i<nbEle;i++)
			sum_of_square += (data_[i] - med)*(data_[i] - med);
		property->zeromean_variance = sum_of_square/nbEle;

		property->minValue = min;
		property->maxValue = max;

		property->avgValue = sum/nbEle;
		property->valueRange = max - min;
		property->totalByteSize = nbEle*sizeof(float);
				
	}
	else //QCAT_DOUBLE
	{
		double* data_ = (double*)data;
		double min=data_[0],max=data_[0],sum=0;		
		for(i=0;i<nbEle;i++)
		{
			if(min>data_[i]) min = data_[i];
			if(max<data_[i]) max = data_[i];
			sum += data_[i];
		}

		double med = min+(max-min)/2;
		double sum_of_square = 0;
		for(i=0;i<nbEle;i++)
			sum_of_square += (data_[i] - med)*(data_[i] - med);
		property->zeromean_variance = sum_of_square/nbEle;

		property->minValue = min;
		property->maxValue = max;

		property->avgValue = sum/nbEle;
		property->valueRange = max - min;	
		property->totalByteSize = nbEle*sizeof(double);	
	}
	
	//entropyType == 0: do not compute any entropy
	//entropyType == 1: compute lossless entropy 8 bits
	//entropyType == 2: compute lossless entropy 32 bits , 64 bits
	//entropyType == 3: compute lossless entropy both 8 bits and 32 bits/64 bits
	if(entropyType == 1 || entropyType == 3)
		property->entropy_8bits = computeLosslessEntropy_8bits(dataType, data, nbEle);
	if(entropyType >= 2)
	{
		if(property->dataType == QCAT_FLOAT)
			property->entropy_32bits = computeLosslessEntropy_32bits(data, nbEle);
		else if(property->dataType == QCAT_DOUBLE)
			property->entropy_64bits = computeLosslessEntropy_64bits(data, nbEle);
	}
	return property;
}

void printProperty(QCAT_DataProperty* property)
{
	printf("numOfElem = %zu\n", property->numOfElem);
	printf("totalDataSize = %zu bytes (%f MB)\n", property->totalByteSize, property->totalByteSize/(1024.0f*1024.0f));
	printf("min = %f\n", property->minValue);
	printf("max = %f\n", property->maxValue);
	printf("valueRange = %f\n", property->valueRange);
	printf("avgValue = %f\n", property->avgValue);
	if(property->entropy_8bits>=0)
		printf("entropy_8bits = %f\n", property->entropy_8bits);
	if(property->entropy_32bits>=0)
		printf("entropy_32bits = %f\n", property->entropy_32bits);
	if(property->entropy_64bits>=0)
		printf("entropy_64bits = %f\n", property->entropy_64bits);		
	printf("zeromean_variance = %f\n", property->zeromean_variance);
	
}

double* computeDataPDF_int32(void* data, size_t numOfElem, int* min, int* intervals)
{
	size_t i = 0;
	int index = 0;
	
	double* dataPDF = NULL;
	int* intData = (int*)data;
	int minData = intData[0];
	int maxData = intData[0];
	for (i = 0; i < numOfElem; i++)
	{
		if(minData > intData[i])
			minData = intData[i];
		if(maxData < intData[i])
			maxData = intData[i];
	}
	int range = maxData - minData;
	
	*min = minData;
	int pdf_intervals = range + 1;
	
	dataPDF = (double*)malloc(sizeof(double)*1000000);
	memset(dataPDF, 0, 1000000*sizeof(double));			
	for (i = 0; i < numOfElem; i++)
	{
		index = intData[i] - minData;
		dataPDF[index] += 1;
	}

	for (i = 0; i < pdf_intervals; i++)
		dataPDF[i]/=numOfElem;	
	
	*intervals = pdf_intervals;
	return dataPDF;		

}

double* computeDataPDF_float(float* data, size_t numOfElem, int intervals, float* min, double* unit, float mint, float maxt)
{
	size_t i = 0;
	int index = 0;
	int threshold = 0;
	double* dataPDF = NULL;
	dataPDF = (double*)malloc(sizeof(double)*intervals);
	memset(dataPDF, 0, intervals*sizeof(double));	
	
	if(mint<maxt)
		threshold = 1;
		
	if(threshold==0) //use all data points
	{
		float minData = data[0];
		float maxData = data[0];
		for (i = 0; i < numOfElem; i++)
		{
			if(minData > data[i])
				minData = data[i];
			if(maxData < data[i])
				maxData = data[i];
		}
		float range = maxData - minData;
		
		*min = minData;
		*unit= range/intervals;
		
		for (i = 0; i < numOfElem; i++)
		{
			index = (int)((data[i] - minData)/(*unit));
			if(index==intervals)
				index--;
			dataPDF[index] += 1;
		}		
	}
	else //focus on only selected value range [min]~[max]
	{
		float minData = 3E38f;
		float maxData = -3E38f;
		for (i = 0; i < numOfElem; i++)
		{
			if(data[i]>=mint && data[i]<=maxt)
			{
				if(minData > data[i])
					minData = data[i];
				if(maxData < data[i])
					maxData = data[i];				
			}
		}
		
		if(maxData==-3E38f || minData==3E38f || maxData==minData)
		{
			printf("Error: no data points or only 1 data point in ([min],[max]).\n");
			exit(0);
		}
		
		float range = maxData - minData;
		
		*min = minData;
		*unit= range/intervals;
		
		for (i = 0; i < numOfElem; i++)
		{
			if(data[i]>=mint && data[i]<=maxt)
			{
				index = (int)((data[i] - minData)/(*unit));
				if(index==intervals)
					index--;
				dataPDF[index] += 1;				
			}
		}				
	}

	for (i = 0; i < intervals; i++)
		dataPDF[i]/=numOfElem;	
	
	return dataPDF;		

}
