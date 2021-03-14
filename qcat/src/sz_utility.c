#include <stdio.h>
#include <sz_utility.h>

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t dataLength;
	if(r1==0) 
	{
		dataLength = 0;
	}
	else if(r2==0) 
	{
		dataLength = r1;
	}
	else if(r3==0) 
	{
		dataLength = r1*r2;
	}
	else if(r4==0) 
	{
		dataLength = r1*r2*r3;
	}
	else if(r5==0) 
	{
		dataLength = r1*r2*r3*r4;
	}
	else 
	{
		dataLength = r1*r2*r3*r4*r5;
	}
	return dataLength;
}

void computeReqLength_float(double realPrecision, short radExpo, int* reqLength, float* medianValue)

{

	short reqExpo = getPrecisionReqLength_double(realPrecision);

	*reqLength = 9+radExpo - reqExpo+1; //radExpo-reqExpo == reqMantiLength

	if(*reqLength<9)

		*reqLength = 9;

	if(*reqLength>32)

	{	

		*reqLength = 32;

		*medianValue = 0;

	}			

}

inline void compressSingleFloatValue(FloatValueCompressElement *vce, float tgtValue, float precision, float medianValue, 
		int reqLength, int reqBytesLength, int resiBitsLength)
{		
	float normValue = tgtValue - medianValue;

	lfloat lfBuf;
	lfBuf.value = normValue;
			
	int ignBytesLength = 32 - reqLength;
	if(ignBytesLength<0)
		ignBytesLength = 0;
	
	int tmp_int = lfBuf.ivalue;
	intToBytes_bigEndian(vce->curBytes, tmp_int);
		
	lfBuf.ivalue = (lfBuf.ivalue >> ignBytesLength) << ignBytesLength;
	
	//float tmpValue = lfBuf.value;
	
	vce->data = lfBuf.value+medianValue;
	vce->curValue = tmp_int;
	vce->reqBytesLength = reqBytesLength;
	vce->resiBitsLength = resiBitsLength;
}

inline void updateLossyCompElement_Float(unsigned char* curBytes, unsigned char* preBytes, 
		int reqBytesLength, int resiBitsLength,  LossyCompressionElement *lce)
{
	int resiIndex, intMidBytes_Length = 0;
	int leadingNum = compIdenticalLeadingBytesCount_float(preBytes, curBytes); //in fact, float is enough for both single-precision and double-precisiond ata.
	int fromByteIndex = leadingNum;
	int toByteIndex = reqBytesLength; //later on: should use "< toByteIndex" to tarverse....
	if(fromByteIndex < toByteIndex)
	{
		intMidBytes_Length = reqBytesLength - leadingNum;
		memcpy(lce->integerMidBytes, &(curBytes[fromByteIndex]), intMidBytes_Length);
	}
	int resiBits = 0;
	if(resiBitsLength!=0)
	{
		resiIndex = reqBytesLength;
		if(resiIndex < 8)
			resiBits = (curBytes[resiIndex] & 0xFF) >> (8-resiBitsLength);
	}
	lce->leadingZeroBytes = leadingNum;
	lce->integerMidBytes_Length = intMidBytes_Length;
	lce->resMidBitsLength = resiBitsLength;
	lce->residualMidBits = resiBits;
}

//TODO double-check the correctness...
inline void addExactData(DynamicByteArray *exactMidByteArray, DynamicIntArray *exactLeadNumArray, 
		DynamicIntArray *resiBitArray, LossyCompressionElement *lce)
{
	int i;
	int leadByteLength = lce->leadingZeroBytes;
	addDIA_Data(exactLeadNumArray, leadByteLength);
	unsigned char* intMidBytes = lce->integerMidBytes;
	int integerMidBytesLength = lce->integerMidBytes_Length;
	int resMidBitsLength = lce->resMidBitsLength;
	if(intMidBytes!=NULL||resMidBitsLength!=0)
	{
		if(intMidBytes!=NULL)
			for(i = 0;i<integerMidBytesLength;i++)
				addDBA_Data(exactMidByteArray, intMidBytes[i]);
		if(resMidBitsLength!=0)
			addDIA_Data(resiBitArray, lce->residualMidBits);
	}
}

int compIdenticalLeadingBytesCount_double(unsigned char* preBytes, unsigned char* curBytes)
{
	int i, n = 0;
	for(i=0;i<8;i++)
		if(preBytes[i]==curBytes[i])
			n++;
		else
			break;
	if(n>3) n = 3;
	return n;
}

inline int compIdenticalLeadingBytesCount_float(unsigned char* preBytes, unsigned char* curBytes)
{
	int i, n = 0;
	for(i=0;i<4;i++)
		if(preBytes[i]==curBytes[i])
			n++;
		else
			break;
	if(n>3) n = 3;
	return n;
}

/**
 * little endian
 * [01|10|11|00|....]-->[01|10|11|00][....]
 * @param timeStepType
 * @return
 */
size_t convertIntArray2ByteArray_fast_2b(unsigned char* timeStepType, size_t timeStepTypeLength, unsigned char **result)
{
	size_t i, j, byteLength = 0;
	if(timeStepTypeLength%4==0)
		byteLength = timeStepTypeLength*2/8;
	else
		byteLength = timeStepTypeLength*2/8+1;
	if(byteLength>0)
		*result = (unsigned char*)malloc(byteLength*sizeof(unsigned char));
	else
		*result = NULL;
	size_t n = 0;
	for(i = 0;i<byteLength;i++)
	{
		int tmp = 0;
		for(j = 0;j<4&&n<timeStepTypeLength;j++)
		{
			int type = timeStepType[n];
			switch(type)
			{
			case 0: 
				
				break;
			case 1:
				tmp = (tmp | (1 << (6-j*2)));
				break;
			case 2:
				tmp = (tmp | (2 << (6-j*2)));
				break;
			case 3:
				tmp = (tmp | (3 << (6-j*2)));
				break;
			default:
				printf("Error: wrong timestep type...: type[%zu]=%d\n", n, type);
				exit(0);
			}
			n++;
		}
		(*result)[i] = (unsigned char)tmp;
	}
	return byteLength;
}

void convertByteArray2IntArray_fast_2b(size_t stepLength, unsigned char* byteArray, size_t byteArrayLength, unsigned char **intArray)
{
	if(stepLength > byteArrayLength*4)
	{
		printf("Error: stepLength > byteArray.length*4\n");
		printf("stepLength=%zu, byteArray.length=%zu\n", stepLength, byteArrayLength);
		exit(0);
	}
	if(stepLength>0)
		*intArray = (unsigned char*)malloc(stepLength*sizeof(unsigned char));
	else
		*intArray = NULL;
	size_t i, n = 0;

	int mod4 = stepLength%4;
	if(mod4==0)
	{
		for (i = 0; i < byteArrayLength; i++) {
			unsigned char tmp = byteArray[i];
			(*intArray)[n++] = (tmp & 0xC0) >> 6;
			(*intArray)[n++] = (tmp & 0x30) >> 4;
			(*intArray)[n++] = (tmp & 0x0C) >> 2;
			(*intArray)[n++] = tmp & 0x03;
		}	
	}
	else
	{
		size_t t = byteArrayLength - mod4;
		for (i = 0; i < t; i++) {
			unsigned char tmp = byteArray[i];
			(*intArray)[n++] = (tmp & 0xC0) >> 6;
			(*intArray)[n++] = (tmp & 0x30) >> 4;
			(*intArray)[n++] = (tmp & 0x0C) >> 2;
			(*intArray)[n++] = tmp & 0x03;
		}
		unsigned char tmp = byteArray[i];				
		switch(mod4)
		{
		case 1:
			(*intArray)[n++] = (tmp & 0xC0) >> 6;
			break;
		case 2:
			(*intArray)[n++] = (tmp & 0xC0) >> 6;
			(*intArray)[n++] = (tmp & 0x30) >> 4;			
			break;
		case 3:	
			(*intArray)[n++] = (tmp & 0xC0) >> 6;
			(*intArray)[n++] = (tmp & 0x30) >> 4;
			(*intArray)[n++] = (tmp & 0x0C) >> 2;		
			break;
		}
	}


}

/**
 * 
 * @param timeStepType is the resiMidBits
 * @param resiBitLength is the length of resiMidBits for each element, (the number of resiBitLength == the # of unpredictable elements
 * @return
 */
size_t convertIntArray2ByteArray_fast_dynamic(unsigned char* timeStepType, unsigned char resiBitLength, size_t nbEle, unsigned char **bytes)
{
	size_t i = 0, j = 0, k = 0; 
	int value;
	DynamicByteArray* dba;
	new_DBA(&dba, 1024);
	int tmp = 0, leftMovSteps = 0;
	for(j = 0;j<nbEle;j++)
	{
		if(resiBitLength==0)
			continue;
		value = timeStepType[i];
		leftMovSteps = getLeftMovingSteps(k, resiBitLength);
		if(leftMovSteps < 0)
		{
			tmp = tmp | (value >> (-leftMovSteps));
			addDBA_Data(dba, (unsigned char)tmp);
			tmp = 0 | (value << (8+leftMovSteps));
		}
		else if(leftMovSteps > 0)
		{
			tmp = tmp | (value << leftMovSteps);
		}
		else //==0
		{
			tmp = tmp | value;
			addDBA_Data(dba, (unsigned char)tmp);
			tmp = 0;
		}
		i++;
		k += resiBitLength;
	}
	if(leftMovSteps != 0)
		addDBA_Data(dba, (unsigned char)tmp);
	convertDBAtoBytes(dba, bytes);
	size_t size = dba->size;
	free_DBA(dba);
	return size;
}

inline int getLeftMovingSteps(size_t k, unsigned char resiBitLength)
{
	return 8 - k%8 - resiBitLength;
}


float computeRangeSize_float(float* oriData, size_t size, float* valueRangeSize, float* medianValue)
{
	size_t i = 0;
	float min = oriData[0];
	float max = min;
	for(i=1;i<size;i++)
	{
		float data = oriData[i];
		if(min>data)
			min = data;
		else if(max<data)
			max = data;
	}

	*valueRangeSize = max - min;
	*medianValue = min + *valueRangeSize/2;
	return min;
}

double computeRangeSize_double(double* oriData, size_t size, double* valueRangeSize, double* medianValue)
{
	size_t i = 0;
	double min = oriData[0];
	double max = min;
	for(i=1;i<size;i++)
	{
		double data = oriData[i];
		if(min>data)
			min = data;
		else if(max<data)
			max = data;
	}
	
	*valueRangeSize = max - min;
	*medianValue = min + *valueRangeSize/2;
	return min;
}

inline size_t bytesToSize(unsigned char* bytes)
{
	size_t result = bytesToLong_bigEndian(bytes);//8	
	return result;
}

inline void sizeToBytes(unsigned char* outBytes, size_t size)
{
	longToBytes_bigEndian(outBytes, size);//8
}
