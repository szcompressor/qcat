#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <qcat.h>
#include <DynamicByteArray.h>
#include <DynamicIntArray.h>
#include <sz_utility.h>

unsigned char* SZ_fast_compress_args_unpredictable_float(int dataType, float *data, size_t *outSize, float absErrBound, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	float valueRangeSize = 0, medianValue = 0;
	
	size_t dataLength = computeDataLength(r5, r4, r3, r2, r1);
	computeRangeSize_float(data, dataLength, &valueRangeSize, &medianValue);	

	size_t i;
	int reqLength;
	short radExpo = getExponent_float(valueRangeSize/2);
	
	computeReqLength_float(absErrBound, radExpo, &reqLength, &medianValue);	

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);
	
	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);
	
	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);
	
	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
	
	//total: 1.6 s
	for(i=0;i<dataLength;i++)
	{
		compressSingleFloatValue(vce, data[i], absErrBound, medianValue, reqLength, reqBytesLength, resiBitsLength); //0.6 s
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce); //0.8 s
		memcpy(preDataBytes,vce->curBytes,sizeof(float));
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);		
	}
	
	//construct output bytes	
	size_t exactDataNum = exactLeadNumArray->size;
	unsigned char* leadNumIntArray = exactLeadNumArray->array;
	size_t residualMidBytes_size = exactMidByteArray->size;
	
	unsigned char* leadNumArray;
	unsigned char* residualMidBytes = exactMidByteArray->array;
	unsigned char* residualMidBits;
	
	size_t leadNumArray_size = convertIntArray2ByteArray_fast_2b(leadNumIntArray, exactDataNum, &leadNumArray);
	size_t residualMidBits_size = convertIntArray2ByteArray_fast_dynamic(resiBitArray->array, resiBitsLength, exactDataNum, &residualMidBits);

	size_t totalSize = 1+sizeof(float)+2*sizeof(size_t)+leadNumArray_size+residualMidBytes_size+residualMidBits_size;
	unsigned char* bytes = (unsigned char*)malloc(totalSize);
	memset(bytes, 0, totalSize);
	size_t k = 0;
	
	unsigned char reqLengthB = (unsigned char)reqLength;
	bytes[k] = reqLengthB;
	k++;
	floatToBytes(&(bytes[k]), medianValue);
	k+=sizeof(float);
	//floatToBytes(&(bytes[k]), absErrBound);
	//k+=sizeof(float);
	sizeToBytes(&(bytes[k]), leadNumArray_size);
	k+=sizeof(size_t);
	sizeToBytes(&(bytes[k]), residualMidBytes_size);
	k+=sizeof(size_t);
	//sizeToBytes(&(bytes[k]), residualMidBits_size);
	//k+=sizeof(size_t);
		
	memcpy(&(bytes[k]), leadNumArray, leadNumArray_size);
	k+=leadNumArray_size;
	memcpy(&(bytes[k]), residualMidBytes, residualMidBytes_size);
	k+=residualMidBytes_size;
	memcpy(&(bytes[k]), residualMidBits, residualMidBits_size);

	free_DIA(exactLeadNumArray);
	free_DBA(exactMidByteArray);
	free_DIA(resiBitArray);
	free(vce);
	free(lce);
	
	free(leadNumArray);
	
	*outSize = totalSize;
	return bytes;
}

void SZ_fast_decompress_args_unpredictable_float(float** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, 
size_t cmpSize)
{
	size_t nbEle = computeDataLength(r5, r4, r3, r2, r1);
	*newData = (float*)malloc(sizeof(float)*nbEle);	
	
	float medianValue, exactData;
	size_t leadNumArray_size, residualMidBytes_size;
	
	size_t k = 0;
	int reqLength = (int)cmpBytes[k];
	k++;
	medianValue = bytesToFloat(&(cmpBytes[k]));
	k+=sizeof(float);
	//absErrBound = bytesToFloat(&(cmpBytes[k]));
	//k+=sizeof(sizeof(float));
	leadNumArray_size = bytesToSize(&(cmpBytes[k]));
	k+=sizeof(size_t);
	residualMidBytes_size = bytesToSize(&(cmpBytes[k]));
	k+=sizeof(size_t);
	//residualMidBits_size = bytesToSize(&(cmpBytes[k]));
	//k+=sizeof(sizeof(size_t));	
	
	unsigned char* leadNumArray = &(cmpBytes[k]);
	k += leadNumArray_size;
	unsigned char* residualMidBytes = &(cmpBytes[k]);	
	k += residualMidBytes_size;
	unsigned char* residualMidBits = &(cmpBytes[k]);
	
	unsigned char* leadNum = NULL;
	convertByteArray2IntArray_fast_2b(nbEle, leadNumArray, leadNumArray_size, &leadNum);		
		
	size_t i = 0;
	size_t j, p = 0, l = 0; 
	k = 0; // k is to track the location of residual_bit
	
	unsigned char preBytes[4];
	unsigned char curBytes[4];

	memset(preBytes, 0, 4);

	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits; 
	unsigned char leadingNum;

	reqBytesLength = reqLength/8;
	resiBitsLength = reqLength%8;
	
	for(i=0;i < nbEle;i++)
	{
		resiBits = 0;
		if (resiBitsLength != 0) {
			int kMod8 = k % 8;
			int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
			if (rightMovSteps > 0) {
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (residualMidBits[p] & code) >> rightMovSteps;
			} else if (rightMovSteps < 0) {
				int code1 = getLeftMovingCode(kMod8);
				int code2 = getRightMovingCode(kMod8, resiBitsLength);
				int leftMovSteps = -rightMovSteps;
				rightMovSteps = 8 - leftMovSteps;
				resiBits = (residualMidBits[p] & code1) << leftMovSteps;
				p++;
				resiBits = resiBits
						| ((residualMidBits[p] & code2) >> rightMovSteps);
			} 
			else // rightMovSteps == 0
			{
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (residualMidBits[p] & code);
				p++;
			}
			k += resiBitsLength;
		}

		// recover the exact data	
		memset(curBytes, 0, 4);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = residualMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}

		exactData = bytesToFloat(curBytes);
		(*newData)[i] = exactData + medianValue;
		memcpy(preBytes,curBytes,4);		
	}
	
}


QCAT_CompressionResult* huffmanAndZstd(int dataType, int* type, int quantBinCapacity, size_t nbEle, void* origData, void* decData)
{
	QCAT_CompressionResult* result= NULL;
	//compress type[] by Huffman encoding
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
			estimatedCompressedSize = nbEle*1.2;
	unsigned char* compressBytes = (unsigned char*)malloc(estimatedCompressedSize);
	size_t zstdOutSize = ZSTD_compress(compressBytes, estimatedCompressedSize, huffmanOutput, huffmanOutSize, 3);	
	
	free(compressBytes);
	free(huffmanOutput);
	
	//analyze compression results
	result = compareData(dataType, nbEle, origData, decData);
	int dataSize = 4+dataType*4;
	result->compressionSize = zstdOutSize;
	//float compressionRatio = ((float)dataSize*nbEle)/zstdOutSize;
	//result->compressionRatio = compressionRatio;	
	
	return result;
}
