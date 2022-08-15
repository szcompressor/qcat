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

#include <sys/time.h>      /* For gettimeofday(), in microseconds */
#include <time.h>          /* For time(), in seconds */

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


void cost_start()
{
    totalCost = 0;
    gettimeofday(&costStart, NULL);
}

void cost_end()
{
    double elapsed;
    struct timeval costEnd;
    gettimeofday(&costEnd, NULL);
    elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
    totalCost += elapsed;
}


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

    cost_start();
    int quantBinCapacity = 65536;

    int stateNum = 2*quantBinCapacity;
    HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
    unsigned char* huffmanOutput = NULL;
    size_t huffmanOutSize = 0;
    encode_withTree(huffmanTree, type, nbEle, &huffmanOutput, &huffmanOutSize);
    //SZ_ReleaseHuffman(huffmanTree);
	cost_end();
	double huffman_compress_time = totalCost;

    //compress huffman output by zstd
    cost_start();
    size_t estimatedCompressedSize = 0;
    if(nbEle < 100)
        estimatedCompressedSize = 200;
    else
        estimatedCompressedSize = nbEle*sizeof(int)*1.2;
    unsigned char* compressBytes = (unsigned char*)malloc(estimatedCompressedSize);
    size_t zstdOutSize = ZSTD_compress(compressBytes, estimatedCompressedSize, huffmanOutput, huffmanOutSize, 3);
	cost_end();
	double zstd_compress_time = totalCost;

    printf("huffmanOutSize = %zu (%f), zstdOutSize = %zu (%f), overallCompressionRatio = %f\n", huffmanOutSize, 1.0*nbEle*sizeof(int)/huffmanOutSize, zstdOutSize, 1.0*huffmanOutSize/zstdOutSize, 1.0f*nbEle*sizeof(int)/zstdOutSize);
    double huffman_compress_throughput = 1.0*nbEle*sizeof(int)/1024/1024/huffman_compress_time;
    double zstd_compress_throughput = huffmanOutSize/1024/1024/zstd_compress_time;
    double overall_compression_throughput = 1.0*nbEle*sizeof(int)/1024/1024/(huffman_compress_time+zstd_compress_time);
	printf("huffman compression time = %f (%f MB/s), zstd compression time = %f (%f MB/s), overall compression time = %f (%f MB/s)\n", huffman_compress_time,  huffman_compress_throughput, zstd_compress_time, zstd_compress_throughput, huffman_compress_time+zstd_compress_time, overall_compression_throughput);

	//check zstd's decompression time
	unsigned char* decData = (unsigned char*)malloc(nbEle*sizeof(int));
	cost_start();
	ZSTD_decompress(decData, nbEle*sizeof(int), compressBytes, zstdOutSize);
	cost_end();
	double zstd_decompression_time = totalCost;

	//check huffman decoding time
	cost_start();

	HuffmanTree* huffmanTree2 = createHuffmanTree(stateNum);
	decode_withTree(huffmanTree2, decData, nbEle, type);
	SZ_ReleaseHuffman(huffmanTree2);
	cost_end();
	double huffman_decompression_time = totalCost;
	
	double overall_decompression_throughput = 1.0*nbEle*sizeof(int)/1024/1024/(huffman_decompression_time + zstd_decompression_time);
	printf("huffman decompression time = %f (%f MB/s), zstd decompression time = %f (%f MB/s), overall decompression time = %f (%f MB/s)\n", huffman_decompression_time, 1.0*nbEle*sizeof(int)/1024/1024/huffman_decompression_time, zstd_decompression_time, huffmanOutSize/1024/1024/zstd_decompression_time, huffman_decompression_time + zstd_decompression_time, overall_decompression_throughput);
	
    free(compressBytes);
    free(huffmanOutput);
    free(decData);

    free(type);

    return 0;
}
