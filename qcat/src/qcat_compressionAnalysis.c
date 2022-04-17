/**
 *  @file qcat_dataAnalysis.c
 *  @author Sheng Di
 *  @date Feb, 2021
 *  @brief data anlaysis
 *  (C) 2015-2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "qcat_dataAnalysis.h"
#include "qcat_compressionAnalysis.h"
#include <stdio.h> 
#include <math.h>
#include <string.h> 
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <Huffman.h>
#include <zstd.h>
#include <qcat_ssim.h>
#include <sz_utility.h>
#include <sz_dummy_compression.h>

double* ZC_compute_autocorrelation1D_double(double* data, size_t numOfElem, double avg)
{
	double *autocorr = (double*)malloc((AUTOCORR_SIZE+1)*sizeof(double));

	size_t i = 0;
	int delta = 0;

	if (numOfElem > 4096)
	{
		double cov = 0;
		for (i = 0; i < numOfElem; i++)
			cov += (data[i] - avg)*(data[i] - avg);

		cov = cov/numOfElem;

		if (cov == 0)
		{
			for (delta = 1; delta <= AUTOCORR_SIZE; delta++)
				autocorr[delta] = 0;
		}
		else
		{
			for(delta = 1; delta <= AUTOCORR_SIZE; delta++)
			{
				double sum = 0;

				for (i = 0; i < numOfElem-delta; i++)
					sum += (data[i]-avg)*(data[i+delta]-avg);

				autocorr[delta] = sum/(numOfElem-delta)/cov;
			}
		}
	}
	else
	{
		for (delta = 1; delta <= AUTOCORR_SIZE; delta++)
		{
			double avg_0 = 0;
			double avg_1 = 0;

			for (i = 0; i < numOfElem-delta; i++)
			{
				avg_0 += data[i];
				avg_1 += data[i+delta];
			}

			avg_0 = avg_0 / (numOfElem-delta);
			avg_1 = avg_1 / (numOfElem-delta);

			double cov_0 = 0;
			double cov_1 = 0;

			for (i = 0; i < numOfElem-delta; i++)
			{
				cov_0 += (data[i]-avg_0)*(data[i]-avg_0);
				cov_1 += (data[i+delta]-avg_1)*(data[i+delta]-avg_1);
			}

			cov_0 = cov_0/(numOfElem-delta);
			cov_1 = cov_1/(numOfElem-delta);

			cov_0 = sqrt(cov_0);
			cov_1 = sqrt(cov_1);

			if (cov_0*cov_1 == 0)
			{
				for (delta = 1; delta <= AUTOCORR_SIZE; delta++)
					autocorr[delta] = 1;
			}
			else
			{
				double sum = 0;

				for (i = 0; i < numOfElem-delta; i++)
					sum += (data[i]-avg_0)*(data[i+delta]-avg_1);

				autocorr[delta] = sum/(numOfElem-delta)/(cov_0*cov_1);
			}
		}
	}

	autocorr[0] = 1;	
	return autocorr;
}

double* ZC_compute_autocorrelation1D_float(float* data, size_t numOfElem, double avg)
{
	double *autocorr = (double*)malloc((AUTOCORR_SIZE+1)*sizeof(double));

	size_t i = 0;
	int delta = 0;

	if (numOfElem > 4096)
	{
		double cov = 0;
		for (i = 0; i < numOfElem; i++)
			cov += (data[i] - avg)*(data[i] - avg);

		cov = cov/numOfElem;

		if (cov == 0)
		{
			for (delta = 1; delta <= AUTOCORR_SIZE; delta++)
				autocorr[delta] = 0;
		}
		else
		{
			for(delta = 1; delta <= AUTOCORR_SIZE; delta++)
			{
				double sum = 0;

				for (i = 0; i < numOfElem-delta; i++)
					sum += (data[i]-avg)*(data[i+delta]-avg);

				autocorr[delta] = sum/(numOfElem-delta)/cov;
			}
		}
	}
	else
	{
		for (delta = 1; delta <= AUTOCORR_SIZE; delta++)
		{
			double avg_0 = 0;
			double avg_1 = 0;

			for (i = 0; i < numOfElem-delta; i++)
			{
				avg_0 += data[i];
				avg_1 += data[i+delta];
			}

			avg_0 = avg_0 / (numOfElem-delta);
			avg_1 = avg_1 / (numOfElem-delta);

			double cov_0 = 0;
			double cov_1 = 0;

			for (i = 0; i < numOfElem-delta; i++)
			{
				cov_0 += (data[i]-avg_0)*(data[i]-avg_0);
				cov_1 += (data[i+delta]-avg_1)*(data[i+delta]-avg_1);
			}

			cov_0 = cov_0/(numOfElem-delta);
			cov_1 = cov_1/(numOfElem-delta);

			cov_0 = sqrt(cov_0);
			cov_1 = sqrt(cov_1);

			if (cov_0*cov_1 == 0)
			{
				for (delta = 1; delta <= AUTOCORR_SIZE; delta++)
					autocorr[delta] = 1;
			}
			else
			{
				double sum = 0;

				for (i = 0; i < numOfElem-delta; i++)
					sum += (data[i]-avg_0)*(data[i+delta]-avg_1);

				autocorr[delta] = sum/(numOfElem-delta)/(cov_0*cov_1);
			}
		}
	}

	autocorr[0] = 1;	
	return autocorr;
}

double* ZC_compute_autocorrelation1D(void* data, int dataType, size_t numOfElem)
{
	double* result = NULL;
	size_t i = 0;
	double avg = 0;
	if(dataType == QCAT_FLOAT)
	{
		float* data_ = (float*)data;
		for(i=0;i<numOfElem;i++)
			avg += data_[i];
		avg /= numOfElem;
		result = ZC_compute_autocorrelation1D_float(data, numOfElem, avg);
	}
	else if(dataType == QCAT_DOUBLE)
	{
		double* data_ = (double*)data;
		for(i=0;i<numOfElem;i++)
			avg += data_[i];
		avg /= numOfElem;		
		result = ZC_compute_autocorrelation1D_double(data, numOfElem, avg);
	}
	else
		return NULL;
	
	return result;
}

double calculateSSIM(void* oriData, void* decData, int dataType, size_t r4, size_t r3, size_t r2, size_t r1)
{
	int dim = computeDimension(0, r4, r3, r2, r1);
	
	int windowSize0 = 7;
	int windowSize1 = 7;
	int windowSize2 = 7;
	int windowSize3 = 7;
	
	int windowShift0 = 2;
	int windowShift1 = 2;
	int windowShift2 = 2;
	int windowShift3 = 2;
	
	double result = -1;
	if(dataType==QCAT_FLOAT) //float type
	{
		switch(dim)
		{
		case 1:
			result = SSIM_1d_windowed_float(oriData, decData, r1, windowSize0, windowShift0);
			break;
		case 2:
			result = SSIM_2d_windowed_float(oriData, decData, r2, r1, windowSize0, windowSize1, windowShift0, windowShift1);
			break;
		case 3:
			result = SSIM_3d_windowed_float(oriData, decData, r3, r2, r1, windowSize0, windowSize1, windowSize2, windowShift0, windowShift1, windowShift2);
			break;
		case 4:
			result = SSIM_4d_windowed_float(oriData, decData, r4, r3, r2, r1, windowSize0, windowSize1, windowSize2, windowSize3, windowShift0, windowShift1, windowShift2, windowShift3);
			break;
		}
	}
	else //double type
	{
		switch(dim)
		{
		case 1:
			result = SSIM_1d_windowed_double(oriData, decData, r1, windowSize0, windowShift0);
			break;
		case 2:
			result = SSIM_2d_windowed_double(oriData, decData, r2, r1, windowSize0, windowSize1, windowShift0, windowShift1);
			break;
		case 3:
			result = SSIM_3d_windowed_double(oriData, decData, r3, r2, r1, windowSize0, windowSize1, windowSize2, windowShift0, windowShift1, windowShift2);
			break;
		case 4:
			result = SSIM_4d_windowed_double(oriData, decData, r4, r3, r2, r1, windowSize0, windowSize1, windowSize2, windowSize3, windowShift0, windowShift1, windowShift2, windowShift3);
			break;
		}		
	}
	return result;
}


double* computeErrPDF(int dataType, void* oriData, void* decData, size_t numOfElem, double fix_interval, double* min_diff, double* err_interval, int* intervals)
{
	size_t i = 0;
	int index = 0;
	double *absErrPDF = NULL;
	if(dataType==QCAT_FLOAT)
	{
		float* data = (float*)oriData;
		float* dec = (float*)decData;
		float* diff = (float*)malloc(sizeof(float)*numOfElem);
		float minDiff = data[0] - dec[0];
		float maxDiff = minDiff;
		for (i = 0; i < numOfElem; i++)
		{
			diff[i] = data[i] - dec[i];
			if(minDiff > diff[i])
				minDiff = diff[i];
			if(maxDiff < diff[i])
				maxDiff = diff[i];
		}
		double diffRange = maxDiff - minDiff;
		double interval = fix_interval;
		int pdf_intervals = PDF_INTERVALS;
		if(interval<=0)
		{
			interval = diffRange/PDF_INTERVALS;
			pdf_intervals = PDF_INTERVALS;
		}
		else
		{
			interval = fix_interval;
			pdf_intervals = (int)(diffRange/interval)+1;
		}
		*min_diff = minDiff;
		*err_interval = interval;
		*intervals = pdf_intervals;

		if(interval==0)
		{
			absErrPDF = (double*)malloc(sizeof(double));
			*absErrPDF = 0;
		}
		else
		{
			absErrPDF = (double*)malloc(sizeof(double)*pdf_intervals);
			memset(absErrPDF, 0, pdf_intervals*sizeof(double));			
			for (i = 0; i < numOfElem; i++)
			{
				index = (int)((diff[i]-minDiff)/interval);
				if(index==pdf_intervals)
					index = pdf_intervals-1;
				absErrPDF[index] += 1;
			}

			for (i = 0; i < pdf_intervals; i++)
				absErrPDF[i]/=numOfElem;	
				
			free(diff);			
		}
	}
	else
	{
		double* data = (double*)oriData;
		double* dec = (double*)decData;
		double* diff = (double*)malloc(sizeof(double)*numOfElem);
		double minDiff = data[0] - dec[0];
		double maxDiff = minDiff;
		for (i = 0; i < numOfElem; i++)
		{
			diff[i] = data[i] - dec[i];
			if(minDiff > diff[i])
				minDiff = diff[i];
			if(maxDiff < diff[i])
				maxDiff = diff[i];
		}
		double diffRange = maxDiff - minDiff;
		double interval = fix_interval;
		int pdf_intervals = PDF_INTERVALS;
		if(interval<=0)
		{
			interval = diffRange/PDF_INTERVALS;
			pdf_intervals = PDF_INTERVALS;
		}
		else
		{
			interval = fix_interval;
			pdf_intervals = (int)(diffRange/interval)+1;
		}
		*min_diff = minDiff;
		*err_interval = interval;
		*intervals = pdf_intervals;

		if(interval==0)
		{
			absErrPDF = (double*)malloc(sizeof(double));
			*absErrPDF = 0;
		}
		else
		{
			absErrPDF = (double*)malloc(sizeof(double)*pdf_intervals);
			memset(absErrPDF, 0, pdf_intervals*sizeof(double));			
			for (i = 0; i < numOfElem; i++)
			{
				index = (int)((diff[i]-minDiff)/interval);
				if(index==pdf_intervals)
					index = pdf_intervals-1;
				absErrPDF[index] += 1;
			}

			for (i = 0; i < pdf_intervals; i++)
				absErrPDF[i]/=numOfElem;	
				
			free(diff);			
		}

	}	
	return absErrPDF;			
}

QCAT_CompressionResult* compareData(int dataType, size_t nbEle, void* data, void* dec)
{
	QCAT_CompressionResult* compressionResult = NULL;
    if(dataType==QCAT_FLOAT) //single precision
	{
		float *ori_data = (float*)data; 
		float *dec_data = (float*)dec;
		
		size_t i = 0;
		double Max = 0, Min = 0, diffMax = 0;
		Max = ori_data[0];
		Min = ori_data[0];
		diffMax = dec_data[0]>ori_data[0]?dec_data[0]-ori_data[0]:ori_data[0]-dec_data[0];

		//diffMax = fabs(data[0] - ori_data[0]);
		double sum1 = 0, sum2 = 0, sum22 = 0;

		for (i = 0; i < nbEle; i++)
		{
			sum1 += ori_data[i];
			sum2 += dec_data[i];
			sum22 += dec_data[i]*dec_data[i];
		}
		double mean1 = sum1/nbEle;
		double mean2 = sum2/nbEle;

		double sum3 = 0, sum4 = 0;
		double sum = 0, prodSum = 0, relerr = 0;

		double maxpw_relerr = 0; 
		for (i = 0; i < nbEle; i++)
		{
			if (Max < ori_data[i]) Max = ori_data[i];
			if (Min > ori_data[i]) Min = ori_data[i];

			float err = fabs(dec_data[i] - ori_data[i]);
			if(ori_data[i]!=0)
			{
				relerr = err/fabs(ori_data[i]);
				if(maxpw_relerr<relerr)
				  maxpw_relerr = relerr;
			}

			if (diffMax < err)
			  diffMax = err;
			prodSum += (ori_data[i]-mean1)*(dec_data[i]-mean2);
			sum3 += (ori_data[i] - mean1)*(ori_data[i]-mean1);
			sum4 += (dec_data[i] - mean2)*(dec_data[i]-mean2);
			sum += err*err;	
		}
		double std1 = sqrt(sum3/nbEle);
		double std2 = sqrt(sum4/nbEle);
		double ee = prodSum/nbEle;
		double pearsonEff = ee/std1/std2;

		double mse = sum/nbEle;
		double range = Max - Min;
		double psnr = 20*log10(range)-10*log10(mse);
		double normErr = sqrt(sum);
		double normErr_norm = normErr/sqrt(sum22);
		double nrmse = sqrt(mse)/range;

		compressionResult = (QCAT_CompressionResult*)malloc(sizeof(QCAT_CompressionResult));
		
		compressionResult->mse = mse;
		compressionResult->psnr = psnr;
		compressionResult->nrmse = nrmse;
		compressionResult->maxABSError = diffMax;
		compressionResult->maxRELError = diffMax/(Max-Min);
		compressionResult->maxPWRError = maxpw_relerr;
		compressionResult->normErr = normErr;
		compressionResult->normErr_norm = normErr_norm;
		compressionResult->pearsonEff = pearsonEff;
		
	} 
	else
	{
		double *ori_data = (double*)data;
		double *dec_data = (double*)dec;
		
		size_t i = 0;
		double Max = 0, Min = 0, diffMax = 0;
		Max = ori_data[0];
		Min = ori_data[0];
		diffMax = dec_data[0]>ori_data[0]?dec_data[0]-ori_data[0]:ori_data[0]-dec_data[0];

		//diffMax = fabs(data[0] - ori_data[0]);
		double sum1 = 0, sum2 = 0, sum22 = 0;

		for (i = 0; i < nbEle; i++)
		{
			sum1 += ori_data[i];
			sum2 += dec_data[i];
			sum22 += dec_data[i]*dec_data[i];
		}
		double mean1 = sum1/nbEle;
		double mean2 = sum2/nbEle;

		double sum3 = 0, sum4 = 0;
		double sum = 0, prodSum = 0, relerr = 0;

		double maxpw_relerr = 0; 
		for (i = 0; i < nbEle; i++)
		{
			if (Max < ori_data[i]) Max = ori_data[i];
			if (Min > ori_data[i]) Min = ori_data[i];

			float err = fabs(dec_data[i] - ori_data[i]);
			if(ori_data[i]!=0)
			{
				relerr = err/fabs(ori_data[i]);
				if(maxpw_relerr<relerr)
				  maxpw_relerr = relerr;
			}

			if (diffMax < err)
			  diffMax = err;
			prodSum += (ori_data[i]-mean1)*(dec_data[i]-mean2);
			sum3 += (ori_data[i] - mean1)*(ori_data[i]-mean1);
			sum4 += (dec_data[i] - mean2)*(dec_data[i]-mean2);
			sum += err*err;	
		}
		double std1 = sqrt(sum3/nbEle);
		double std2 = sqrt(sum4/nbEle);
		double ee = prodSum/nbEle;
		double pearsonEff = ee/std1/std2;

		double mse = sum/nbEle;
		double range = Max - Min;
		double psnr = 20*log10(range)-10*log10(mse);
		double normErr = sqrt(sum);
		double normErr_norm = normErr/sqrt(sum22);
		double nrmse = sqrt(mse)/range;

		compressionResult = (QCAT_CompressionResult*)malloc(sizeof(QCAT_CompressionResult));
		
		compressionResult->mse = mse;
		compressionResult->psnr = psnr;
		compressionResult->nrmse = nrmse;
		compressionResult->maxABSError = diffMax;
		compressionResult->maxRELError = diffMax/(Max-Min);
		compressionResult->maxPWRError = maxpw_relerr;
		compressionResult->normErr = normErr;
		compressionResult->normErr_norm = normErr_norm;
		compressionResult->pearsonEff = pearsonEff;
	}
	
	return compressionResult;
}

QCAT_CompressionResult* getCompressionResult(int dataType, float errBound, int quantBinCapacity, void* origData, void* predData, QCAT_DataProperty* property)
{
	size_t i = 0;
	size_t n = property->numOfElem;

	QCAT_CompressionResult* result = NULL;
	int* type = (int*)malloc(sizeof(int)*n);
	memset(type, 0, sizeof(int)*n);
	
	float checkRadius = (quantBinCapacity-1)*errBound;
	int intvRadius = quantBinCapacity/2;
	float interval = 2*errBound;
	
	size_t unpredictableCount = 0;
	
	if(dataType == QCAT_FLOAT)
	{
		//float* delta = (float*)malloc(sizeof(float)*n); //created for possible post analysis of delta
		float* data = (float*)origData;
		float* pred = (float*)predData;
		float* dec = (float*)malloc(sizeof(float)*n);
		float* unpred = (float*)malloc(sizeof(float)*n);
		for(i = 0;i < n;i++)
		{
			float d = fabsf(data[i] - pred[i]);
			if(d <= checkRadius)
			{
				int state = ((int)(d/errBound+1))>>1;
				float decValue = 0;
				if(data[i] >= pred[i])
				{
					type[i] = intvRadius+state;
					decValue = pred[i] + state*interval;
				}
				else
				{
					type[i] = intvRadius-state;
					decValue = pred[i] - state*interval;
				}
				dec[i] = decValue;				
			}
			else
			{
				type[i] = 0;
				//Collect unpredictable data data[i]
				unpred[unpredictableCount++] = data[i];
			}
		}
		//compress unpredictable data
		size_t unpredOutSize = 0;
		unsigned char* unpredCompBytes = SZ_fast_compress_args_unpredictable_float(dataType, unpred, &unpredOutSize, errBound, 0, 0, 0, 0, unpredictableCount);
		
		float* dec_unpred = NULL;
		SZ_fast_decompress_args_unpredictable_float(&dec_unpred, 0, 0, 0, 0, unpredictableCount, unpredCompBytes, unpredOutSize);
		
		size_t j = 0;
		for(i = 0; i < n;i++)
		{
			if(type[i]==0)
				dec[i] = dec_unpred[j++];
		}
		
		result = huffmanAndZstd(dataType, type, quantBinCapacity, n, data, dec);
		size_t zstdSize = result->compressionSize;
		result->compressionSize = zstdSize + unpredOutSize;
		result->compressionRatio = 1.0f*n*sizeof(float)/result->compressionSize;
		
		//free(delta);
		free(unpred);
		free(dec);
	}
	else if(dataType == QCAT_DOUBLE)
	{
		//double* delta = (double*)malloc(sizeof(double)*n); //for possible post analysis of delta
		double* data = (double*)origData;
		double* pred = (double*)predData;		
		double* dec = (double*)malloc(sizeof(double)*n);
		for(i = 0;i < n;i++)
		{
			double d = fabs(data[i] - pred[i]);
			if(d <=  checkRadius)
			{
				int state = ((int)(d/errBound+1))>>1;
				double decValue = 0;
				if(data[i] >= pred[i])
				{
					type[i] = intvRadius+state;
					decValue = pred[i] + state*interval;
				}
				else
				{
					type[i] = intvRadius-state;
					decValue = pred[i] - state*interval;
				}
				dec[i] = decValue;	
			}
			else
			{
				type[i] = 0;
				//TODO: compress unpredictable data data[i]
				unpredictableCount ++;				
			}
		}
		
		result = huffmanAndZstd(dataType, type, quantBinCapacity, n, data, dec);
		//free(delta);
		free(dec);
	}
	
	
	result->unpredictableCount = unpredictableCount;
	result->unpredictablePercent = 1.0f*unpredictableCount/n;
	free(type);
	return result;
}

void printCompressionResult(QCAT_CompressionResult* result)
{	
	printf("compressionRatio = %f\n", result->compressionRatio);
	printf("max absolute error = %f\n", result->maxABSError);
	printf("max relative error = %f\n", result->maxRELError);
	printf("max point-wise relative error  = %f\n", result->maxPWRError);	
	printf("mean squared error (MSE) = %.20G\n", result->mse);	
	printf("normalized root mean squared error (NRMSE) = %.20G\n", result->nrmse);
	printf("peak signal-to-noise ratio (PSNR) = %f\n", result->psnr);
	printf("norm Error = %.20G\n", result->normErr);
	printf("norm Error (norm) = %.20G\n", result->normErr_norm);
	printf("pearson correlation coefficient = %f\n", result->pearsonEff);
	printf("unpredictableCount = %zu\n", result->unpredictableCount);
	printf("unpredictablePercent = %f\n", result->unpredictablePercent);
}
