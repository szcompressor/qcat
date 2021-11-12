/**
 *  @file compressQuantBins.c
 *  @author Sheng Di
 *  @date Nov, 2021
 *  @brief This is an example of using compression interface
 *  (C) 2021 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include "ByteToolkit.h"
#include "rw.h"
#include "qcat_dataAnalysis.h"
#include "sz_dummy_compression.h"
#include "math.h"

int main(int argc, char * argv[])
{
	int status = 0;
	char oriFilePath[640];

	if(argc < 2)
	{
		printf("Usage: compressQuantBins [inputFilePath]\n");
		exit(0);
	}

	sprintf(oriFilePath, "%s", argv[1]);

	int x = 1;
	char *y = (char*)&x;

	if(*y==1)
	{     
		sysEndianType = 0; //LITTLE_ENDIAN_SYSTEM;
	}
	else //=0
	{
		sysEndianType = 1; //BIG_ENDIAN_SYSTEM;
	}
	
	size_t nbEle = 0;
	int * type = readInt32Data(oriFilePath, &nbEle, &status);

	int quantBinCapacity = 65536;

        int stateNum = 2*quantBinCapacity;
        HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
        unsigned char* huffmanOutput = NULL;
        size_t huffmanOutSize = 0;
        encode_withTree(huffmanTree, type, nbEle, &huffmanOutput, &huffmanOutSize);
        SZ_ReleaseHuffman(huffmanTree);

        //compress huffman output by zstd
        size_t estimatedCompressedSize = 0;
        if(nbEle < 100)
                        estimatedCompressedSize = 200;
        else
                        estimatedCompressedSize = nbEle*sizeof(int)*1.2;
        unsigned char* compressBytes = (unsigned char*)malloc(estimatedCompressedSize);
        size_t zstdOutSize = ZSTD_compress(compressBytes, estimatedCompressedSize, huffmanOutput, huffmanOutSize, 3);

        free(compressBytes);
        free(huffmanOutput);

	printf("huffmanOutSize = %zu (%f), zstdOutSize = %zu (%f)\n", huffmanOutSize, 1.0f*nbEle/huffmanOutSize, zstdOutSize, 1.0f*nbEle/zstdOutSize);

	free(type);

	return 0;
}
