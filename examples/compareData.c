/**
 *  @file compareData.c
 *  @author Sheng Di
 *  @date Sept, 2020
 *  @brief This is an example of using compression interface
 *  (C) 2020 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ByteToolkit.h"
#include "rw.h"

int main(int argc, char * argv[])
{
    int status = 0;
    char dataType[4];
    char oriFilePath[640], decFilePath[640];
    if(argc < 4)
    {
		printf("Test case: compareData [datatype (-f or -d)] [original data file] [decompressed data file] <-p>\n");
		printf("			-f means single precision; -d means double precision\n");
		printf("Example1: compareData -f testfloat_8_8_128.dat testfloat_8_8_128.dat.sz.out \n");
		printf("Example2: compareData -f testfloat_8_8_128.dat testfloat_8_8_128.dat.sz.out -p\n");
		printf("Note: -p is optional. It indicates printing the compression errors in a text file.\n");
		exit(0);
    }
    int printError = 0;
    sprintf(dataType, "%s", argv[1]);
    sprintf(oriFilePath, "%s", argv[2]);
    sprintf(decFilePath, "%s", argv[3]);
    if(argc == 5)
	    printError = 1;
    int data_type = strcmp(dataType, "-f")==0? sizeof(float): sizeof(double);
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
   
   
    if(strcmp(dataType, "-f")==0) //single precision
	{
		float* errors = NULL;
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
		
    		if(printError)
			errors = (float*)malloc(data_type*nbEle);
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
			float err_ = dec_data[i] - ori_data[i];
			if(printError)
				errors[i] = err_;
			float err = fabs(err_);
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
		double acEff = ee/std1/std2;

		double mse = sum/nbEle;
		double range = Max - Min;
		double psnr = 20*log10(range)-10*log10(mse);
		double normErr = sqrt(sum);
		double normErr_norm = normErr/sqrt(sum22);
		double nrmse = sqrt(mse)/range;


		printf ("Min = %.20G, Max = %.20G, range = %.20G\n", Min, Max, range);
		printf ("Max absolute error = %.10f\n", diffMax);
		printf ("Max relative error = %f\n", diffMax/(Max-Min));
		printf ("Max pw relative error = %f\n", maxpw_relerr);
		printf ("PSNR = %f, NRMSE = %.20G\n", psnr,nrmse);
		printf ("normErr = %f, normErr_norm = %f\n", normErr, normErr_norm);
		printf ("pearson coeff = %f\n", acEff);

		if(printError)
		{
			char errFile[644];
			sprintf(errFile, "%s.err", decFilePath);
			printf ("Printing the errors into a file: %s ....\n", errFile);
			writeFloatData(errors, nbEle, errFile, &status);
			printf ("Done\n");
		}
		
		free(ori_data);
		free(dec_data);	
		free(errors);
	} 
	else
	{
		double* errors = NULL;
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

    		if(printError)
			errors = (double*)malloc(data_type*nbEle);
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

			float err_ = dec_data[i] - ori_data[i];
			if(printError)
				errors[i] = err_;
			float err = fabs(err_);
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
		double acEff = ee/std1/std2;

		double mse = sum/nbEle;
		double range = Max - Min;
		double psnr = 20*log10(range)-10*log10(mse);
		double normErr = sqrt(sum);
		double normErr_norm = normErr/sqrt(sum22);
		double nrmse = sqrt(mse)/range;


		printf ("Min = %.20G, Max = %.20G, range = %.20G\n", Min, Max, range);
		printf ("Max absolute error = %.10f\n", diffMax);
		printf ("Max relative error = %f\n", diffMax/(Max-Min));
		printf ("Max pw relative error = %f\n", maxpw_relerr);
		printf ("PSNR = %f, NRMSE = %.20G\n", psnr,nrmse);
		printf ("normErr = %f, normErr_norm = %f\n", normErr, normErr_norm);
		printf ("pearson coeff = %f\n", acEff);

		if(printError)
		{
			char errFile[644];
			sprintf(errFile, "%s.err", decFilePath);
			printf ("Printing the errors into a file: %s ....\n", errFile);
			writeDoubleData(errors, nbEle, errFile, &status);
			printf ("Done\n");
		}

		free(ori_data);
		free(dec_data);
		free(errors);
	}

    
    return 0;
}
