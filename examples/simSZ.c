/**
 *  @file simSZ.c
 *  @author Sheng Di
 *  @date June, 2021
 *  @brief This is an example of using compression interface
 *  (C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ByteToolkit.h"
#include "rw.h"
#include "qcat.h"

int main(int argc, char * argv[])
{
	int status = 0;
	char dataType[4];
	char huffmanQuantBinFile[640], unpredictableDataFile[640];

	if(argc < 3)
	{
		printf("Test case: simSZ [datatype (-f or -d)] [quantBinCapacity] [quantization bin file] [error bound] [unpredictable data file]\n");
		printf("			-f means single precision; -d means double precision (support only -f now)\n");
		printf("Example: simSZ -f 65536 1E-3 huffmanQuantBin.int32 unpredictableData.f32\n");
		exit(0);
	}

	sprintf(dataType, "%s", argv[1]);
	int quantBinCapacity = atoi(argv[2]);
	float errBound = atof(argv[3]);
	sprintf(huffmanQuantBinFile, "%s", argv[4]);
	sprintf(unpredictableDataFile, "%s", argv[5]);

	size_t nbEle;
	int* type = NULL;
	float* unpredData;
	size_t unpredCount = 0;

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
   
	if(strcmp(dataType, "-f")==0)
	{
		//dType = QCAT_FLOAT;
		type = readInt32Data(huffmanQuantBinFile, &nbEle, &status);
		unpredData = readFloatData(unpredictableDataFile, &unpredCount, &status);

		printf("the number of elements is %zu.\n", nbEle);	
		printf("the number of unpred data is %zu (%f).\n", unpredCount, ((float)unpredCount)/((float)nbEle));
		
		printf("Performing compression ....\n");
		size_t outSize = 0;
		unsigned char* bytes = huffmanAndZstdCompress(type, nbEle, quantBinCapacity, errBound, unpredData, unpredCount, &outSize);
		
		printf("compressed size = %zu, CR = %f\n", outSize, 1.0f*(sizeof(float)*nbEle)/outSize);
		printf("Performing decompression ....\n");
		int *dec_type = NULL;
		float *dec_unpData = NULL;
		unsigned int dec_unpredictableCount = 0;
		zstdAndHuffmanDecompress(quantBinCapacity, bytes, outSize, nbEle, &dec_type, &dec_unpData, &dec_unpredictableCount);
		
		printf("the number of unpred data is %d.\n", dec_unpredictableCount);
		
		//verification
		size_t i = 0;
		for(i=0;i<nbEle;i++)
		{
			if(type[i] != dec_type[i])
			{
				printf("Error: type[%zu]==%d, dec_type[%zu]==%d.\n", i, type[i], i, dec_type[i]);
				exit(0);
			}
		}
		printf("Verification of quantization bin PASS: type is identical to dec_type\n");
				
		if(unpredCount!=dec_unpredictableCount)
		{
			printf("Error: unpredictableCount!=dec_unpredictableCount\n");
			exit(0);
		}
		float max = 0;
		for(i=0;i<dec_unpredictableCount;i++)
		{
			float error = fabsf(unpredData[i] - dec_unpData[i]);
			if(max < error)
				max = error;
		}
		if(max < errBound)
			printf("Unpredictable data: ERROR BOUND check PASS\n");
		else
			printf("Unpredictable data: ERROR BOUND check FAIL\n");
		printf("Max Error = %f\n", max);
	}
	else
	{
		printf("Error: doesn't support -d\n");
		exit(0);
	}
	
	free(unpredData);
    
	return 0;
}
