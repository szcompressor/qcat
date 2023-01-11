#include <stdio.h>
#include <sz_utility.h>

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1) {
    size_t dataLength;
    if (r1 == 0) {
        dataLength = 0;
    } else if (r2 == 0) {
        dataLength = r1;
    } else if (r3 == 0) {
        dataLength = r1 * r2;
    } else if (r4 == 0) {
        dataLength = r1 * r2 * r3;
    } else if (r5 == 0) {
        dataLength = r1 * r2 * r3 * r4;
    } else {
        dataLength = r1 * r2 * r3 * r4 * r5;
    }
    return dataLength;
}

int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1) {
    if (r1 == 0)
        return 0;
    else if (r2 == 0)
        return 1;
    else if (r3 == 0)
        return 2;
    else if (r4 == 0)
        return 3;
    else if (r5 == 0)
        return 4;
    return 5;
}

void computeReqLength_float(double realPrecision, short radExpo, int *reqLength, float *medianValue) {

    short reqExpo = getPrecisionReqLength_double(realPrecision);

    *reqLength = 9 + radExpo - reqExpo + 1; //radExpo-reqExpo == reqMantiLength

    if (*reqLength < 9)

        *reqLength = 9;

    if (*reqLength > 32) {

        *reqLength = 32;

        *medianValue = 0;

    }

}

inline void compressSingleFloatValue(FloatValueCompressElement *vce, float tgtValue, float precision, float medianValue,
                                     int reqLength, int reqBytesLength, int resiBitsLength) {
    float normValue = tgtValue - medianValue;

    lfloat lfBuf;
    lfBuf.value = normValue;

    int ignBytesLength = 32 - reqLength;
    if (ignBytesLength < 0)
        ignBytesLength = 0;

    int tmp_int = lfBuf.ivalue;
    intToBytes_bigEndian(vce->curBytes, tmp_int);

    lfBuf.ivalue = (lfBuf.ivalue >> ignBytesLength) << ignBytesLength;

    //float tmpValue = lfBuf.value;

    vce->data = lfBuf.value + medianValue;
    vce->curValue = tmp_int;
    vce->reqBytesLength = reqBytesLength;
    vce->resiBitsLength = resiBitsLength;
}

inline void updateLossyCompElement_Float(unsigned char *curBytes, unsigned char *preBytes,
                                         int reqBytesLength, int resiBitsLength, LossyCompressionElement *lce) {
    int resiIndex, intMidBytes_Length = 0;
    int leadingNum = compIdenticalLeadingBytesCount_float(preBytes,
                                                          curBytes); //in fact, float is enough for both single-precision and double-precisiond ata.
    int fromByteIndex = leadingNum;
    int toByteIndex = reqBytesLength; //later on: should use "< toByteIndex" to tarverse....
    if (fromByteIndex < toByteIndex) {
        intMidBytes_Length = reqBytesLength - leadingNum;
        memcpy(lce->integerMidBytes, &(curBytes[fromByteIndex]), intMidBytes_Length);
    }
    int resiBits = 0;
    if (resiBitsLength != 0) {
        resiIndex = reqBytesLength;
        if (resiIndex < 8)
            resiBits = (curBytes[resiIndex] & 0xFF) >> (8 - resiBitsLength);
    }
    lce->leadingZeroBytes = leadingNum;
    lce->integerMidBytes_Length = intMidBytes_Length;
    lce->resMidBitsLength = resiBitsLength;
    lce->residualMidBits = resiBits;
}

//TODO double-check the correctness...
inline void addExactData(DynamicByteArray *exactMidByteArray, DynamicIntArray *exactLeadNumArray,
                         DynamicIntArray *resiBitArray, LossyCompressionElement *lce) {
    int i;
    int leadByteLength = lce->leadingZeroBytes;
    addDIA_Data(exactLeadNumArray, leadByteLength);
    unsigned char *intMidBytes = lce->integerMidBytes;
    int integerMidBytesLength = lce->integerMidBytes_Length;
    int resMidBitsLength = lce->resMidBitsLength;
    if (intMidBytes != NULL || resMidBitsLength != 0) {
        if (intMidBytes != NULL)
            for (i = 0; i < integerMidBytesLength; i++)
                addDBA_Data(exactMidByteArray, intMidBytes[i]);
        if (resMidBitsLength != 0)
            addDIA_Data(resiBitArray, lce->residualMidBits);
    }
}

int compIdenticalLeadingBytesCount_double(unsigned char *preBytes, unsigned char *curBytes) {
    int i, n = 0;
    for (i = 0; i < 8; i++)
        if (preBytes[i] == curBytes[i])
            n++;
        else
            break;
    if (n > 3) n = 3;
    return n;
}

inline int compIdenticalLeadingBytesCount_float(unsigned char *preBytes, unsigned char *curBytes) {
    int i, n = 0;
    for (i = 0; i < 4; i++)
        if (preBytes[i] == curBytes[i])
            n++;
        else
            break;
    if (n > 3) n = 3;
    return n;
}

/**
 * little endian
 * [01|10|11|00|....]-->[01|10|11|00][....]
 * @param timeStepType
 * @return
 */
size_t convertIntArray2ByteArray_fast_2b(unsigned char *timeStepType, size_t timeStepTypeLength, unsigned char **result) {
    size_t i, j, byteLength = 0;
    if (timeStepTypeLength % 4 == 0)
        byteLength = timeStepTypeLength * 2 / 8;
    else
        byteLength = timeStepTypeLength * 2 / 8 + 1;
    if (byteLength > 0)
        *result = (unsigned char *) malloc(byteLength * sizeof(unsigned char));
    else
        *result = NULL;
    size_t n = 0;
    for (i = 0; i < byteLength; i++) {
        int tmp = 0;
        for (j = 0; j < 4 && n < timeStepTypeLength; j++) {
            int type = timeStepType[n];
            switch (type) {
                case 0:

                    break;
                case 1:
                    tmp = (tmp | (1 << (6 - j * 2)));
                    break;
                case 2:
                    tmp = (tmp | (2 << (6 - j * 2)));
                    break;
                case 3:
                    tmp = (tmp | (3 << (6 - j * 2)));
                    break;
                default:
                    printf("Error: wrong timestep type...: type[%zu]=%d\n", n, type);
                    exit(0);
            }
            n++;
        }
        (*result)[i] = (unsigned char) tmp;
    }
    return byteLength;
}

void convertByteArray2IntArray_fast_2b(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned char **intArray) {
    if (stepLength > byteArrayLength * 4) {
        printf("Error: stepLength > byteArray.length*4\n");
        printf("stepLength=%zu, byteArray.length=%zu\n", stepLength, byteArrayLength);
        exit(0);
    }
    if (stepLength > 0)
        *intArray = (unsigned char *) malloc(stepLength * sizeof(unsigned char));
    else
        *intArray = NULL;
    size_t i, n = 0;

    int mod4 = stepLength % 4;
    if (mod4 == 0) {
        for (i = 0; i < byteArrayLength; i++) {
            unsigned char tmp = byteArray[i];
            (*intArray)[n++] = (tmp & 0xC0) >> 6;
            (*intArray)[n++] = (tmp & 0x30) >> 4;
            (*intArray)[n++] = (tmp & 0x0C) >> 2;
            (*intArray)[n++] = tmp & 0x03;
        }
    } else {
        size_t t = byteArrayLength - mod4;
        for (i = 0; i < t; i++) {
            unsigned char tmp = byteArray[i];
            (*intArray)[n++] = (tmp & 0xC0) >> 6;
            (*intArray)[n++] = (tmp & 0x30) >> 4;
            (*intArray)[n++] = (tmp & 0x0C) >> 2;
            (*intArray)[n++] = tmp & 0x03;
        }
        unsigned char tmp = byteArray[i];
        switch (mod4) {
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
size_t convertIntArray2ByteArray_fast_dynamic(unsigned char *timeStepType, unsigned char resiBitLength, size_t nbEle, unsigned char **bytes) {
    size_t i = 0, j = 0, k = 0;
    int value;
    DynamicByteArray *dba;
    new_DBA(&dba, 1024);
    int tmp = 0, leftMovSteps = 0;
    for (j = 0; j < nbEle; j++) {
        if (resiBitLength == 0)
            continue;
        value = timeStepType[i];
        leftMovSteps = getLeftMovingSteps(k, resiBitLength);
        if (leftMovSteps < 0) {
            tmp = tmp | (value >> (-leftMovSteps));
            addDBA_Data(dba, (unsigned char) tmp);
            tmp = 0 | (value << (8 + leftMovSteps));
        } else if (leftMovSteps > 0) {
            tmp = tmp | (value << leftMovSteps);
        } else //==0
        {
            tmp = tmp | value;
            addDBA_Data(dba, (unsigned char) tmp);
            tmp = 0;
        }
        i++;
        k += resiBitLength;
    }
    if (leftMovSteps != 0)
        addDBA_Data(dba, (unsigned char) tmp);
    convertDBAtoBytes(dba, bytes);
    size_t size = dba->size;
    free_DBA(dba);
    return size;
}

inline int getLeftMovingSteps(size_t k, unsigned char resiBitLength) {
    return 8 - k % 8 - resiBitLength;
}


float computeRangeSize_float(float *oriData, size_t size, float *valueRangeSize, float *medianValue) {
    size_t i = 0;
    float min = oriData[0];
    float max = min;
    for (i = 1; i < size; i++) {
        float data = oriData[i];
        if (min > data)
            min = data;
        else if (max < data)
            max = data;
    }

    *valueRangeSize = max - min;
    *medianValue = min + *valueRangeSize / 2;
    return min;
}

double computeRangeSize_double(double *oriData, size_t size, double *valueRangeSize, double *medianValue) {
    size_t i = 0;
    double min = oriData[0];
    double max = min;
    for (i = 1; i < size; i++) {
        double data = oriData[i];
        if (min > data)
            min = data;
        else if (max < data)
            max = data;
    }

    *valueRangeSize = max - min;
    *medianValue = min + *valueRangeSize / 2;
    return min;
}

inline size_t bytesToSize(unsigned char *bytes) {
    size_t result = bytesToLong_bigEndian(bytes);//8
    return result;
}

inline void sizeToBytes(unsigned char *outBytes, size_t size) {
    longToBytes_bigEndian(outBytes, size);//8
}

int lorenzoPredictorQuant_Cmpr_NoOutlier_float(float *data, int mode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, int *out) {
    size_t i = 0, j = 0, k = 0;
    size_t nbEle = computeDataLength(0, 0, n3, n2, n1);
    int dim = computeDimension(0, 0, n3, n2, n1);

    register float data_recip = 0;
    register float recip_Precision = 0.5f / errorBound;

    if (mode == LORENZO_1D_1LAYER) {
        register int s = 0;
        register int prev_quant_value = 0;
        register int curr_quant_value = 0;
        register int x = 0;

        data_recip = data[0] * recip_Precision;
        s = data_recip >= -0.5f ? 0 : 1;
        out[0] = (int) (data_recip + 0.5f) - s;
        prev_quant_value = out[0];

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (i = 1; i < nbEle; i++) {
                data_recip = data[i] * recip_Precision;
                s = data_recip >= -0.5f ? 0 : 1;
                curr_quant_value = (int) (data_recip + 0.5f) - s;

                out[i] = curr_quant_value - prev_quant_value;
                prev_quant_value = curr_quant_value;
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (i = 1; i < nbEle; i++) {
                data_recip = data[i] * recip_Precision;
                s = data_recip >= -0.5f ? 0 : 1;
                curr_quant_value = (int) (data_recip + 0.5f) - s;

                x = curr_quant_value - prev_quant_value;
                out[i] = (x << 1) ^ (x >> 31);
                prev_quant_value = curr_quant_value;
            }
        } else //QUANT_CODE_SHIFT
        {

        }

    } else if (mode == LORENZO_1D_2LAYER) {
        register int s = 0;
        register int curr_quant_value = 0;
        register int x = 0;

        data_recip = data[0] * recip_Precision;
        s = data_recip >= -0.5f ? 0 : 1;
        out[0] = (int) (data_recip + 0.5f) - s;

        data_recip = data[1] * recip_Precision;
        s = data_recip >= -0.5f ? 0 : 1;
        out[1] = (int) (data_recip + 0.5f) - s;

        register int prev2_quant_value = out[0];
        register int prev1_quant_value = out[1];

        //quantization  & 2-layer 1D prediction
        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (i = 2; i < nbEle; i++) {
                data_recip = data[i] * recip_Precision;
                s = data_recip >= -0.5f ? 0 : 1;
                curr_quant_value = (int) (data_recip + 0.5f) - s;

                out[i] = curr_quant_value - (2 * prev1_quant_value - prev2_quant_value);
                prev2_quant_value = prev1_quant_value;
                prev1_quant_value = curr_quant_value;
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (i = 2; i < nbEle; i++) {
                data_recip = data[i] * recip_Precision;
                s = data_recip >= -0.5f ? 0 : 1;
                curr_quant_value = (int) (data_recip + 0.5f) - s;

                x = curr_quant_value - (2 * prev1_quant_value - prev2_quant_value);
                out[i] = (x << 1) ^ (x >> 31);
                prev2_quant_value = prev1_quant_value;
                prev1_quant_value = curr_quant_value;
            }
        } else //QUANT_CODE_SHIFT
        {

        }
    } else if (mode == LORENZO_1D_3LAYER) {
        register int s = 0;
        register int curr_quant_value = 0;
        register int x = 0;

        data_recip = data[0] * recip_Precision;
        s = data_recip >= -0.5f ? 0 : 1;
        out[0] = (int) (data_recip + 0.5f) - s;

        data_recip = data[1] * recip_Precision;
        s = data_recip >= -0.5f ? 0 : 1;
        out[1] = (int) (data_recip + 0.5f) - s;

        data_recip = data[2] * recip_Precision;
        s = data_recip >= -0.5f ? 0 : 1;
        out[2] = (int) (data_recip + 0.5f) - s;

        register int prev3_quant_value = out[0];
        register int prev2_quant_value = out[1];
        register int prev1_quant_value = out[2];

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (i = 3; i < nbEle; i++) {
                data_recip = data[i] * recip_Precision;
                s = data_recip >= -0.5f ? 0 : 1;
                curr_quant_value = (int) (data_recip + 0.5f) - s;

                out[i] = curr_quant_value - (3 * prev1_quant_value - 3 * prev2_quant_value + prev3_quant_value);
                prev3_quant_value = prev2_quant_value;
                prev2_quant_value = prev1_quant_value;
                prev1_quant_value = curr_quant_value;
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (i = 3; i < nbEle; i++) {
                data_recip = data[i] * recip_Precision;
                s = data_recip >= -0.5f ? 0 : 1;
                curr_quant_value = (int) (data_recip + 0.5f) - s;

                x = curr_quant_value - (3 * prev1_quant_value - 3 * prev2_quant_value + prev3_quant_value);
                out[i] = (x << 1) ^ (x >> 31);
                prev3_quant_value = prev2_quant_value;
                prev2_quant_value = prev1_quant_value;
                prev1_quant_value = curr_quant_value;
            }
        } else //QUANT_CODE_SHIFT
        {

        }
    } else if (mode == LORENZO_2D_1LAYER) {
        if (dim == 1) {
            return 1; //error! The dimension is 1 but the mode is required to be LORENZO_2D_1LAYER
        }
        register int x = 0, s = 0;
        size_t dim2_offset = n1 + 1;

        int *quant_with_buffer = calloc(2 * dim2_offset, sizeof(int));
        int *quant = quant_with_buffer + dim2_offset + 1;

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (j = 0; j < n3 * n2; j++) {
                for (i = 0; i < n1; i++) {
                    data_recip = data[j * n1 + i] * recip_Precision;
                    s = (data_recip >= -0.5f ? 0 : 1);
                    quant[i] = (int) (data_recip + 0.5f) - s;

                    out[j * n1 + i] = quant[i] - (quant[i - 1] + quant[i - dim2_offset] - quant[i - 1 - dim2_offset]);
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim2_offset, dim2_offset);
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (j = 0; j < n3 * n2; j++) {
                for (i = 0; i < n1; i++) {
                    data_recip = data[j * n1 + i] * recip_Precision;
                    s = (data_recip >= -0.5f ? 0 : 1);
                    quant[i] = (int) (data_recip + 0.5f) - s;

                    x = quant[i] - (quant[i - 1] + quant[i - dim2_offset] - quant[i - 1 - dim2_offset]);
                    out[j * n1 + i] = (x << 1) ^ (x >> 31);
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim2_offset, dim2_offset);
            }
        } else {//QUANT_CODE_SHIFT

        }

    } else if (mode == LORENZO_3D_1LAYER) {
        if (dim == 1 || dim == 2) {
            return 1;
        }
        register int x = 0, s = 0;
        size_t dim2_offset = n1 + 1;
        size_t dim3_offset = (n1 + 1) * (n2 + 1);

        int *quant_with_buffer = calloc(2 * dim3_offset, sizeof(int));
        int *quant = quant_with_buffer + dim3_offset + dim2_offset + 1;

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (k = 0; k < n3; k++) {
                int *quant_p = quant;
                for (j = 0; j < n2; j++) {
                    for (i = 0; i < n1; i++) {
                        data_recip = data[k * n2 * n1 + j * n1 + i] * recip_Precision;
                        s = (data_recip >= -0.5f ? 0 : 1);
                        *quant_p = (int) (data_recip + 0.5f) - s;

                        out[k * n1 * n2 + j * n1 + i] = *quant_p -
                                                        (quant_p[-1] + quant_p[-dim2_offset] + quant_p[-dim3_offset]
                                                         - quant_p[-dim2_offset - 1] - quant_p[-dim3_offset - 1] -
                                                         quant_p[-dim3_offset - dim2_offset] +
                                                         quant_p[-dim3_offset - dim2_offset - 1]);
                        quant_p++;
                    }
                    quant_p++;
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim3_offset, dim3_offset);
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (k = 0; k < n3; k++) {
                int *quant_p = quant;
                for (j = 0; j < n2; j++) {
                    for (i = 0; i < n1; i++) {
                        data_recip = data[k * n2 * n1 + j * n1 + i] * recip_Precision;
                        s = (data_recip >= -0.5f ? 0 : 1);
                        *quant_p = (int) (data_recip + 0.5f) - s;

                        x = *quant_p - (quant_p[-1] + quant_p[-dim2_offset] + quant_p[-dim3_offset]
                                        - quant_p[-dim2_offset - 1] - quant_p[-dim3_offset - 1] - quant_p[-dim3_offset - dim2_offset] +
                                        quant_p[-dim3_offset - dim2_offset - 1]);
                        out[k * n1 * n2 + j * n1 + i] = (x << 1) ^ (x >> 31);
                        quant_p++;
                    }
                    quant_p++;
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim3_offset, dim3_offset);
            }
        } else {//QUANT_CODE_SHIFT

        }
    }

    return 0;

}

int
lorenzoPredictorQuant_Cmpr_NoOutlier_double(double *data, int mode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, int *out) {
    size_t i = 0, j = 0, k = 0;
    size_t nbEle = computeDataLength(0, 0, n3, n2, n1);
    int dim = computeDimension(0, 0, n3, n2, n1);

    register double data_recip = 0;
    register double recip_Precision = 0.5 / errorBound;

    if (mode == LORENZO_1D_1LAYER) {
        register int s = 0;
        register int prev_quant_value = 0;
        register int curr_quant_value = 0;
        register int x = 0;

        data_recip = data[0] * recip_Precision;
        s = data_recip >= -0.5 ? 0 : 1;
        out[0] = (int) (data_recip + 0.5) - s;
        prev_quant_value = out[0];
        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (i = 1; i < nbEle; i++) {
                data_recip = data[i] * recip_Precision;
                s = data_recip >= -0.5 ? 0 : 1;
                curr_quant_value = (int) (data_recip + 0.5) - s;

                out[i] = curr_quant_value - prev_quant_value;
                prev_quant_value = curr_quant_value;
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (i = 1; i < nbEle; i++) {
                data_recip = data[i] * recip_Precision;
                s = data_recip >= -0.5 ? 0 : 1;
                curr_quant_value = (int) (data_recip + 0.5) - s;

                x = curr_quant_value - prev_quant_value;
                out[i] = (x << 1) ^ (x >> 31);
                prev_quant_value = curr_quant_value;
            }
        } else //codeFormat == QUANT_CODE_SHIFT
        {

        }
    } else if (mode == LORENZO_1D_2LAYER) {
        register int s = 0;
        register int curr_quant_value = 0;
        register int x = 0;

        data_recip = data[0] * recip_Precision;
        s = data_recip >= -0.5 ? 0 : 1;
        out[0] = (int) (data_recip + 0.5) - s;

        data_recip = data[1] * recip_Precision;
        s = data_recip >= -0.5 ? 0 : 1;
        out[1] = (int) (data_recip + 0.5) - s;

        register int prev2_quant_value = out[0];
        register int prev1_quant_value = out[1];

        //quantization  & 2-layer 1D prediction
        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (i = 2; i < nbEle; i++) {
                data_recip = data[i] * recip_Precision;
                s = data_recip >= -0.5 ? 0 : 1;
                curr_quant_value = (int) (data_recip + 0.5) - s;

                out[i] = curr_quant_value - (2 * prev1_quant_value - prev2_quant_value);
                prev2_quant_value = prev1_quant_value;
                prev1_quant_value = curr_quant_value;
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            data_recip = data[i] * recip_Precision;
            s = data_recip >= -0.5 ? 0 : 1;
            curr_quant_value = (int) (data_recip + 0.5) - s;

            x = curr_quant_value - (2 * prev1_quant_value - prev2_quant_value);
            out[i] = (x << 1) ^ (x >> 31);
            prev2_quant_value = prev1_quant_value;
            prev1_quant_value = curr_quant_value;
        } else //codeFormat == QUANT_CODE_SHIFT
        {

        }
    } else if (mode == LORENZO_1D_3LAYER) {
        register int s = 0;
        register int curr_quant_value = 0;
        register int x = 0;

        data_recip = data[0] * recip_Precision;
        s = data_recip >= -0.5 ? 0 : 1;
        out[0] = (int) (data_recip + 0.5) - s;

        data_recip = data[1] * recip_Precision;
        s = data_recip >= -0.5 ? 0 : 1;
        out[1] = (int) (data_recip + 0.5) - s;

        data_recip = data[2] * recip_Precision;
        s = data_recip >= -0.5 ? 0 : 1;
        out[2] = (int) (data_recip + 0.5) - s;

        register int prev3_quant_value = out[0];
        register int prev2_quant_value = out[1];
        register int prev1_quant_value = out[2];

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (i = 3; i < nbEle; i++) {
                data_recip = data[i] * recip_Precision;
                s = data_recip >= -0.5 ? 0 : 1;
                curr_quant_value = (int) (data_recip + 0.5) - s;

                out[i] = curr_quant_value - (3 * prev1_quant_value - 3 * prev2_quant_value + prev3_quant_value);
                prev3_quant_value = prev2_quant_value;
                prev2_quant_value = prev1_quant_value;
                prev1_quant_value = curr_quant_value;
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (i = 3; i < nbEle; i++) {
                data_recip = data[i] * recip_Precision;
                s = data_recip >= -0.5 ? 0 : 1;
                curr_quant_value = (int) (data_recip + 0.5) - s;

                x = curr_quant_value - (3 * prev1_quant_value - 3 * prev2_quant_value + prev3_quant_value);
                out[i] = (x << 1) ^ (x >> 31);
                prev3_quant_value = prev2_quant_value;
                prev2_quant_value = prev1_quant_value;
                prev1_quant_value = curr_quant_value;
            }
        } else //codeFormat == QUANT_CODE_SHIFT
        {

        }
    } else if (mode == LORENZO_2D_1LAYER) {
        if (dim == 1) {
            return 1; //error! The dimension is 1 but the mode is required to be LORENZO_2D_1LAYER
        }
        register int x = 0, s = 0;
        size_t dim2_offset = n1 + 1;

        int *quant_with_buffer = calloc(2 * dim2_offset, sizeof(int));
        int *quant = quant_with_buffer + dim2_offset + 1;

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (j = 0; j < n3 * n2; j++) {
                for (i = 0; i < n1; i++) {
                    data_recip = data[j * n1 + i] * recip_Precision;
                    s = (data_recip >= -0.5f ? 0 : 1);
                    quant[i] = (int) (data_recip + 0.5f) - s;

                    out[j * n1 + i] = quant[i] - (quant[i - 1] + quant[i - dim2_offset] - quant[i - 1 - dim2_offset]);
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim2_offset, dim2_offset);
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (j = 0; j < n3 * n2; j++) {
                for (i = 0; i < n1; i++) {
                    data_recip = data[j * n1 + i] * recip_Precision;
                    s = (data_recip >= -0.5f ? 0 : 1);
                    quant[i] = (int) (data_recip + 0.5f) - s;

                    x = quant[i] - (quant[i - 1] + quant[i - dim2_offset] - quant[i - 1 - dim2_offset]);
                    out[j * n1 + i] = (x << 1) ^ (x >> 31);
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim2_offset, dim2_offset);
            }
        } else {//QUANT_CODE_SHIFT

        }

    } else if (mode == LORENZO_3D_1LAYER) {
        if (dim == 1 || dim == 2) {
            return 1;
        }
        register int x = 0, s = 0;
        size_t dim2_offset = n1 + 1;
        size_t dim3_offset = (n1 + 1) * (n2 + 1);

        int *quant_with_buffer = calloc(2 * dim3_offset, sizeof(int));
        int *quant = quant_with_buffer + dim3_offset + dim2_offset + 1;

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (k = 0; k < n3; k++) {
                int *quant_p = quant;
                for (j = 0; j < n2; j++) {
                    for (i = 0; i < n1; i++) {
                        data_recip = data[k * n2 * n1 + j * n1 + i] * recip_Precision;
                        s = (data_recip >= -0.5f ? 0 : 1);
                        *quant_p = (int) (data_recip + 0.5f) - s;

                        out[k * n1 * n2 + j * n1 + i] = *quant_p -
                                                        (quant_p[-1] + quant_p[-dim2_offset] + quant_p[-dim3_offset]
                                                         - quant_p[-dim2_offset - 1] - quant_p[-dim3_offset - 1] -
                                                         quant_p[-dim3_offset - dim2_offset] +
                                                         quant_p[-dim3_offset - dim2_offset - 1]);
                        quant_p++;
                    }
                    quant_p++;
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim3_offset, dim3_offset);
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (k = 0; k < n3; k++) {
                int *quant_p = quant;
                for (j = 0; j < n2; j++) {
                    for (i = 0; i < n1; i++) {
                        data_recip = data[k * n2 * n1 + j * n1 + i] * recip_Precision;
                        s = (data_recip >= -0.5f ? 0 : 1);
                        *quant_p = (int) (data_recip + 0.5f) - s;

                        x = *quant_p - (quant_p[-1] + quant_p[-dim2_offset] + quant_p[-dim3_offset]
                                        - quant_p[-dim2_offset - 1] - quant_p[-dim3_offset - 1] - quant_p[-dim3_offset - dim2_offset] +
                                        quant_p[-dim3_offset - dim2_offset - 1]);
                        out[k * n1 * n2 + j * n1 + i] = (x << 1) ^ (x >> 31);
                        quant_p++;
                    }
                    quant_p++;
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim3_offset, dim3_offset);
            }
        } else {//QUANT_CODE_SHIFT

        }
    }

    return 0;
}

<<<<<<< HEAD
int lorenzoPredictorQuant_Decmpr_NoOutlier_float(int* diffQuantData, int mode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, float* result)
{
	size_t i = 0;
	
	size_t nbEle = computeDataLength(0, 0, n3, n2, n1);
	int dim = computeDimension(0, 0, n3, n2, n1);	
	
	register int curQuantValue = 0;
	register float e2 = errorBound*2;
	register int x = 0;
	
	if(mode == LORENZO_1D_1LAYER)
	{
		int preQuantValue = diffQuantData[0];
		result[0] = e2*preQuantValue;
		
		if(codeFormat == QUANT_CODE_ORIGINAL)
		{
			for(i=1;i<nbEle;i++)
			{
				curQuantValue = preQuantValue + diffQuantData[i];
				result[i] = e2*curQuantValue;			
				preQuantValue = curQuantValue;
			}			
		}
		else if(codeFormat == QUANT_CODE_NORMALIZE)
		{
			for(i=1;i<nbEle;i++)
			{
				x = diffQuantData[i];
				curQuantValue = preQuantValue + (((x >> 1) ^ ((x&1)*0xFFFFFFFF)));
				result[i] = e2*curQuantValue;			
				preQuantValue = curQuantValue;
			}				
		}
		else //QUANT_CODE_SHIFT
		{
			
		}
	}
	else if(mode == LORENZO_1D_2LAYER)
	{
		result[0] = e2*diffQuantData[0];
		result[1] = e2*diffQuantData[1];
		register int pred2_quant_value = diffQuantData[0];
		register int pred1_quant_value = diffQuantData[1];		
		if(codeFormat == QUANT_CODE_ORIGINAL)
		{						
			for(i=2;i<nbEle;i++)
			{
				curQuantValue = (2*pred1_quant_value - pred2_quant_value) + diffQuantData[i];
				result[i] = e2*curQuantValue;
				pred2_quant_value = pred1_quant_value;
				pred1_quant_value = curQuantValue;
			}
		}
		else if(codeFormat == QUANT_CODE_NORMALIZE)
		{
			for(i=2;i<nbEle;i++)
			{
				x = diffQuantData[i];
				curQuantValue = (2*pred1_quant_value - pred2_quant_value) + (((x >> 1) ^ ((x&1)*0xFFFFFFFF)));
				result[i] = e2*curQuantValue;
				pred2_quant_value = pred1_quant_value;
				pred1_quant_value = curQuantValue;
			}			
		}	
		else //QUANT_CODE_SHIFT
		{
			
		}
	}	
	else if(mode == LORENZO_1D_3LAYER)
	{
		result[0] = e2*diffQuantData[0];
		result[1] = e2*diffQuantData[1];
		result[2] = e2*diffQuantData[2];		
		register int pred3_quant_value = diffQuantData[0];
		register int pred2_quant_value = diffQuantData[1];		
		register int pred1_quant_value = diffQuantData[2];				

		if(codeFormat == QUANT_CODE_ORIGINAL)
		{							
			for(i=3;i<nbEle;i++)
			{
				curQuantValue = (3*pred1_quant_value - 3*pred2_quant_value + pred3_quant_value) + diffQuantData[i];
				result[i] = e2*curQuantValue;
				pred3_quant_value = pred2_quant_value;
				pred2_quant_value = pred1_quant_value;
				pred1_quant_value = curQuantValue;
			}
		}
		else if(codeFormat == QUANT_CODE_NORMALIZE)				
		{
			for(i=3;i<nbEle;i++)
			{
				x = diffQuantData[i];
				curQuantValue = (3*pred1_quant_value - 3*pred2_quant_value + pred3_quant_value) + (((x >> 1) ^ ((x&1)*0xFFFFFFFF)));
				result[i] = e2*curQuantValue;
				pred3_quant_value = pred2_quant_value;
				pred2_quant_value = pred1_quant_value;
				pred1_quant_value = curQuantValue;
			}			
		}
		else //QUANT_CODE_SHIFT
		{
			
		}
	}
	else if(mode == LORENZO_2D_1LAYER)
	{
		if(dim==1)
			return 1; //error! The dimension is 1 but the mode is required to be LORENZO_2D_1LAYER
		//TODO implement 2D lorenzo			
	}
	else if(mode == LORENZO_3D_1LAYER)
	{
		if(dim==1 || dim==2)
			return 1;
		//TODO implement 3D lorenzo		
	}
	
	return 0;	
}

int lorenzoPredictorQuant_Decmpr_NoOutlier_double(int* diffQuantData, int mode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, double* result)
{
	size_t i = 0;
	
	size_t nbEle = computeDataLength(0, 0, n3, n2, n1);
	int dim = computeDimension(0, 0, n3, n2, n1);	
	
	register int curQuantValue = 0;
	register double e2 = errorBound*2;
	register int x = 0;
		
	if(mode == LORENZO_1D_1LAYER)
	{
		int preQuantValue = diffQuantData[0];
		result[0] = e2*preQuantValue;
		
		if(codeFormat == QUANT_CODE_ORIGINAL)
		{		
			for(i=1;i<nbEle;i++)
			{
				curQuantValue = preQuantValue + diffQuantData[i];
				result[i] = e2*curQuantValue;			
				preQuantValue = curQuantValue;
			}
		}
		else if(codeFormat == QUANT_CODE_NORMALIZE)
		{
			for(i=1;i<nbEle;i++)
			{
				x = diffQuantData[i];
				curQuantValue = preQuantValue + (((x >> 1) ^ ((x&1)*0xFFFFFFFFFFFFFFFF)));
				result[i] = e2*curQuantValue;			
				preQuantValue = curQuantValue;
			}			
		}
		else
		{
			
		}
	}
	else if(mode == LORENZO_1D_2LAYER)
	{
		result[0] = e2*diffQuantData[0];
		result[1] = e2*diffQuantData[1];
		register int pred2_quant_value = diffQuantData[0];
		register int pred1_quant_value = diffQuantData[1];		

		if(codeFormat == QUANT_CODE_ORIGINAL)
		{							
			for(i=2;i<nbEle;i++)
			{	
				curQuantValue = (2*pred1_quant_value - pred2_quant_value) + diffQuantData[i];
				result[i] = e2*curQuantValue;
				pred2_quant_value = pred1_quant_value;
				pred1_quant_value = curQuantValue;
			}		
		}
		else if(codeFormat == QUANT_CODE_NORMALIZE)
		{
			for(i=2;i<nbEle;i++)
			{	
				x = diffQuantData[i];
				curQuantValue = (2*pred1_quant_value - pred2_quant_value) + (((x >> 1) ^ ((x&1)*0xFFFFFFFFFFFFFFFF)));
				result[i] = e2*curQuantValue;
				pred2_quant_value = pred1_quant_value;
				pred1_quant_value = curQuantValue;
			}					
		}
		else if(codeFormat == QUANT_CODE_SHIFT)
		{
			
		}
	}	
	else if(mode == LORENZO_1D_3LAYER)
	{
		result[0] = e2*diffQuantData[0];
		result[1] = e2*diffQuantData[1];
		result[2] = e2*diffQuantData[2];		
		register int pred3_quant_value = diffQuantData[0];
		register int pred2_quant_value = diffQuantData[1];		
		register int pred1_quant_value = diffQuantData[2];				
					
		if(codeFormat == QUANT_CODE_ORIGINAL)
		{								
			for(i=3;i<nbEle;i++)
			{
				curQuantValue = (3*pred1_quant_value - 3*pred2_quant_value + pred3_quant_value) + diffQuantData[i];
				result[i] = e2*curQuantValue;
				pred3_quant_value = pred2_quant_value;
				pred2_quant_value = pred1_quant_value;
				pred1_quant_value = curQuantValue;
			}				
		}
		else if(codeFormat == QUANT_CODE_NORMALIZE)
		{
			for(i=3;i<nbEle;i++)
			{
				x = diffQuantData[i];
				curQuantValue = (3*pred1_quant_value - 3*pred2_quant_value + pred3_quant_value) + (((x >> 1) ^ ((x&1)*0xFFFFFFFFFFFFFFFF)));
				result[i] = e2*curQuantValue;
				pred3_quant_value = pred2_quant_value;
				pred2_quant_value = pred1_quant_value;
				pred1_quant_value = curQuantValue;
			}						
		}
		else if(codeFormat == QUANT_CODE_SHIFT)
		{
			
		}
	}
	else if(mode == LORENZO_2D_1LAYER)
	{
		if(dim==1)
			return 1; //error! The dimension is 1 but the mode is required to be LORENZO_2D_1LAYER
		//TODO implement 2D lorenzo			
	}
	else if(mode == LORENZO_3D_1LAYER)
	{
		if(dim==1 || dim==2)
			return 1;
		//TODO implement 3D lorenzo		
	}
	return 0;	
=======
int lorenzoPredictorQuant_Decmpr_NoOutlier_float(int *diffQuantData, int mode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1,
                                                 float *result) {
    size_t i = 0, j = 0, k = 0;

    size_t nbEle = computeDataLength(0, 0, n3, n2, n1);
    int dim = computeDimension(0, 0, n3, n2, n1);

    register int curQuantValue = 0;
    register float e2 = errorBound * 2;
    register int x = 0;

    if (mode == LORENZO_1D_1LAYER) {
        int preQuantValue = diffQuantData[0];
        result[0] = e2 * preQuantValue;

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (i = 1; i < nbEle; i++) {
                curQuantValue = preQuantValue + diffQuantData[i];
                result[i] = e2 * curQuantValue;
                preQuantValue = curQuantValue;
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (i = 1; i < nbEle; i++) {
                x = diffQuantData[i];
                curQuantValue = preQuantValue + ((x >> 1) ^ (-(x & 1)));
                result[i] = e2 * curQuantValue;
                preQuantValue = curQuantValue;
            }
        } else //QUANT_CODE_SHIFT
        {

        }
    } else if (mode == LORENZO_1D_2LAYER) {
        result[0] = e2 * diffQuantData[0];
        result[1] = e2 * diffQuantData[1];
        register int pred2_quant_value = diffQuantData[0];
        register int pred1_quant_value = diffQuantData[1];
        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (i = 2; i < nbEle; i++) {
                curQuantValue = (2 * pred1_quant_value - pred2_quant_value) + diffQuantData[i];
                result[i] = e2 * curQuantValue;
                pred2_quant_value = pred1_quant_value;
                pred1_quant_value = curQuantValue;
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (i = 2; i < nbEle; i++) {
                x = diffQuantData[i];
                curQuantValue = (2 * pred1_quant_value - pred2_quant_value) + ((x >> 1) ^ (-(x & 1)));
                result[i] = e2 * curQuantValue;
                pred2_quant_value = pred1_quant_value;
                pred1_quant_value = curQuantValue;
            }
        } else //QUANT_CODE_SHIFT
        {

        }
    } else if (mode == LORENZO_1D_3LAYER) {
        result[0] = e2 * diffQuantData[0];
        result[1] = e2 * diffQuantData[1];
        result[2] = e2 * diffQuantData[2];
        register int pred3_quant_value = diffQuantData[0];
        register int pred2_quant_value = diffQuantData[1];
        register int pred1_quant_value = diffQuantData[2];

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (i = 3; i < nbEle; i++) {
                curQuantValue = (3 * pred1_quant_value - 3 * pred2_quant_value + pred3_quant_value) + diffQuantData[i];
                result[i] = e2 * curQuantValue;
                pred3_quant_value = pred2_quant_value;
                pred2_quant_value = pred1_quant_value;
                pred1_quant_value = curQuantValue;
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (i = 3; i < nbEle; i++) {
                x = diffQuantData[i];
                curQuantValue = (3 * pred1_quant_value - 3 * pred2_quant_value + pred3_quant_value) + ((x >> 1) ^ (-(x & 1)));
                result[i] = e2 * curQuantValue;
                pred3_quant_value = pred2_quant_value;
                pred2_quant_value = pred1_quant_value;
                pred1_quant_value = curQuantValue;
            }
        } else //QUANT_CODE_SHIFT
        {

        }
    } else if (mode == LORENZO_2D_1LAYER) {
        if (dim == 1)
            return 1; //error! The dimension is 1 but the mode is required to be LORENZO_2D_1LAYER
        size_t dim2_offset = n1 + 1;

        int *quant_with_buffer = calloc(2 * dim2_offset, sizeof(int));
        int *quant = quant_with_buffer + dim2_offset + 1;

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (j = 0; j < n3 * n2; j++) {
                for (i = 0; i < n1; i++) {
                    quant[i] = diffQuantData[j * n1 + i] + (quant[i - 1] + quant[i - dim2_offset] - quant[i - 1 - dim2_offset]);
                    result[j * n1 + i] = e2 * quant[i];
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim2_offset, dim2_offset);
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (j = 0; j < n3 * n2; j++) {
                for (i = 0; i < n1; i++) {
                    x = diffQuantData[j * n1 + i];
                    quant[i] = ((x >> 1) ^ (-(x & 1))) + (quant[i - 1] + quant[i - dim2_offset] - quant[i - 1 - dim2_offset]);
                    result[j * n1 + i] = e2 * quant[i];
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim2_offset, dim2_offset);
            }
        } else {//QUANT_CODE_SHIFT

        }
    } else if (mode == LORENZO_3D_1LAYER) {
        if (dim == 1 || dim == 2) {
            return 1;
        }
        size_t dim2_offset = n1 + 1;
        size_t dim3_offset = (n1 + 1) * (n2 + 1);

        int *quant_with_buffer = calloc(2 * dim3_offset, sizeof(int));
        int *quant = quant_with_buffer + dim3_offset + dim2_offset + 1;

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (k = 0; k < n3; k++) {
                int *quant_p = quant;
                for (j = 0; j < n2; j++) {
                    for (i = 0; i < n1; i++) {
                        *quant_p = diffQuantData[k * n2 * n1 + j * n1 + i] + (quant_p[-1] + quant_p[-dim2_offset] + quant_p[-dim3_offset]
                                                                              - quant_p[-dim2_offset - 1] - quant_p[-dim3_offset - 1] -
                                                                              quant_p[-dim3_offset - dim2_offset] +
                                                                              quant_p[-dim3_offset - dim2_offset - 1]);
                        result[k * n2 * n1 + j * n1 + i] = e2 * (*quant_p);

                        quant_p++;
                    }
                    quant_p++;
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim3_offset, dim3_offset);
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (k = 0; k < n3; k++) {
                int *quant_p = quant;
                for (j = 0; j < n2; j++) {
                    for (i = 0; i < n1; i++) {
                        x = diffQuantData[k * n2 * n1 + j * n1 + i];
                        *quant_p = ((x >> 1) ^ (-(x & 1))) + (quant_p[-1] + quant_p[-dim2_offset] + quant_p[-dim3_offset]
                                                              - quant_p[-dim2_offset - 1] - quant_p[-dim3_offset - 1] -
                                                              quant_p[-dim3_offset - dim2_offset] +
                                                              quant_p[-dim3_offset - dim2_offset - 1]);
                        result[k * n2 * n1 + j * n1 + i] = e2 * (*quant_p);
                        quant_p++;
                    }
                    quant_p++;
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim3_offset, dim3_offset);
            }
        } else {//QUANT_CODE_SHIFT

        }
    }

    return 0;
}

int lorenzoPredictorQuant_Decmpr_NoOutlier_double(int *diffQuantData, int mode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1,
                                                  double *result) {
    size_t i = 0, j = 0, k = 0;

    size_t nbEle = computeDataLength(0, 0, n3, n2, n1);
    int dim = computeDimension(0, 0, n3, n2, n1);

    register int curQuantValue = 0;
    register double e2 = errorBound * 2;
    register int x = 0;

    if (mode == LORENZO_1D_1LAYER) {
        int preQuantValue = diffQuantData[0];
        result[0] = e2 * preQuantValue;

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (i = 1; i < nbEle; i++) {
                curQuantValue = preQuantValue + diffQuantData[i];
                result[i] = e2 * curQuantValue;
                preQuantValue = curQuantValue;
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (i = 1; i < nbEle; i++) {
                x = diffQuantData[i];
                curQuantValue = preQuantValue + ((x >> 1) ^ (-(x & 1)));
                result[i] = e2 * curQuantValue;
                preQuantValue = curQuantValue;
            }
        } else {

        }
    } else if (mode == LORENZO_1D_2LAYER) {
        result[0] = e2 * diffQuantData[0];
        result[1] = e2 * diffQuantData[1];
        register int pred2_quant_value = diffQuantData[0];
        register int pred1_quant_value = diffQuantData[1];

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (i = 2; i < nbEle; i++) {
                curQuantValue = (2 * pred1_quant_value - pred2_quant_value) + diffQuantData[i];
                result[i] = e2 * curQuantValue;
                pred2_quant_value = pred1_quant_value;
                pred1_quant_value = curQuantValue;
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (i = 2; i < nbEle; i++) {
                x = diffQuantData[i];
                curQuantValue = (2 * pred1_quant_value - pred2_quant_value) + ((x >> 1) ^ (-(x & 1)));
                result[i] = e2 * curQuantValue;
                pred2_quant_value = pred1_quant_value;
                pred1_quant_value = curQuantValue;
            }
        } else if (codeFormat == QUANT_CODE_SHIFT) {

        }
    } else if (mode == LORENZO_1D_3LAYER) {
        result[0] = e2 * diffQuantData[0];
        result[1] = e2 * diffQuantData[1];
        result[2] = e2 * diffQuantData[2];
        register int pred3_quant_value = diffQuantData[0];
        register int pred2_quant_value = diffQuantData[1];
        register int pred1_quant_value = diffQuantData[2];

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (i = 3; i < nbEle; i++) {
                curQuantValue = (3 * pred1_quant_value - 3 * pred2_quant_value + pred3_quant_value) + diffQuantData[i];
                result[i] = e2 * curQuantValue;
                pred3_quant_value = pred2_quant_value;
                pred2_quant_value = pred1_quant_value;
                pred1_quant_value = curQuantValue;
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (i = 3; i < nbEle; i++) {
                x = diffQuantData[i];
                curQuantValue = (3 * pred1_quant_value - 3 * pred2_quant_value + pred3_quant_value) + ((x >> 1) ^ (-(x & 1)));
                result[i] = e2 * curQuantValue;
                pred3_quant_value = pred2_quant_value;
                pred2_quant_value = pred1_quant_value;
                pred1_quant_value = curQuantValue;
            }
        } else if (codeFormat == QUANT_CODE_SHIFT) {

        }
    } else if (mode == LORENZO_2D_1LAYER) {
        if (dim == 1)
            return 1; //error! The dimension is 1 but the mode is required to be LORENZO_2D_1LAYER
        size_t dim2_offset = n1 + 1;

        int *quant_with_buffer = calloc(2 * dim2_offset, sizeof(int));
        int *quant = quant_with_buffer + dim2_offset + 1;

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (j = 0; j < n3 * n2; j++) {
                for (i = 0; i < n1; i++) {
                    quant[i] = diffQuantData[j * n1 + i] + (quant[i - 1] + quant[i - dim2_offset] - quant[i - 1 - dim2_offset]);
                    result[j * n1 + i] = e2 * quant[i];
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim2_offset, dim2_offset);
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (j = 0; j < n3 * n2; j++) {
                for (i = 0; i < n1; i++) {
                    x = diffQuantData[j * n1 + i];
                    quant[i] = ((x >> 1) ^ (-(x & 1))) + (quant[i - 1] + quant[i - dim2_offset] - quant[i - 1 - dim2_offset]);
                    result[j * n1 + i] = e2 * quant[i];
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim2_offset, dim2_offset);
            }
        } else {//QUANT_CODE_SHIFT

        }
    } else if (mode == LORENZO_3D_1LAYER) {
        if (dim == 1 || dim == 2) {
            return 1;
        }
        size_t dim2_offset = n1 + 1;
        size_t dim3_offset = (n1 + 1) * (n2 + 1);

        int *quant_with_buffer = calloc(2 * dim3_offset, sizeof(int));
        int *quant = quant_with_buffer + dim3_offset + dim2_offset + 1;

        if (codeFormat == QUANT_CODE_ORIGINAL) {
            for (k = 0; k < n3; k++) {
                int *quant_p = quant;
                for (j = 0; j < n2; j++) {
                    for (i = 0; i < n1; i++) {
                        *quant_p = diffQuantData[k * n2 * n1 + j * n1 + i] + (quant_p[-1] + quant_p[-dim2_offset] + quant_p[-dim3_offset]
                                                                              - quant_p[-dim2_offset - 1] - quant_p[-dim3_offset - 1] -
                                                                              quant_p[-dim3_offset - dim2_offset] +
                                                                              quant_p[-dim3_offset - dim2_offset - 1]);
                        result[k * n2 * n1 + j * n1 + i] = e2 * (*quant_p);

                        quant_p++;
                    }
                    quant_p++;
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim3_offset, dim3_offset);
            }
        } else if (codeFormat == QUANT_CODE_NORMALIZE) {
            for (k = 0; k < n3; k++) {
                int *quant_p = quant;
                for (j = 0; j < n2; j++) {
                    for (i = 0; i < n1; i++) {
                        x = diffQuantData[k * n2 * n1 + j * n1 + i];
                        *quant_p = ((x >> 1) ^ (-(x & 1))) + (quant_p[-1] + quant_p[-dim2_offset] + quant_p[-dim3_offset]
                                                              - quant_p[-dim2_offset - 1] - quant_p[-dim3_offset - 1] -
                                                              quant_p[-dim3_offset - dim2_offset] +
                                                              quant_p[-dim3_offset - dim2_offset - 1]);
                        result[k * n2 * n1 + j * n1 + i] = e2 * (*quant_p);
                        quant_p++;
                    }
                    quant_p++;
                }
                memcpy(quant_with_buffer, quant_with_buffer + dim3_offset, dim3_offset);
            }
        } else {//QUANT_CODE_SHIFT

        }
    }

    return 0;
>>>>>>> 93817303910da322ecce3e7dea9ee32f4fa9713a
}


int
lorenzoPredictorQuant_Cmpr_NoOutlier(void *data, int dataType, int lorenzoMode, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1,
                                     int *out) {
    int status = 0;
    if (dataType == QCAT_FLOAT) {
        status = lorenzoPredictorQuant_Cmpr_NoOutlier_float(data, lorenzoMode, codeFormat, errorBound, n3, n2, n1, out);
    } else if (dataType == QCAT_DOUBLE) {
        status = lorenzoPredictorQuant_Cmpr_NoOutlier_double(data, lorenzoMode, codeFormat, errorBound, n3, n2, n1, out);
    }
    return status;
}

int lorenzoPredictorQuant_Decmpr_NoOutlier(int *diffQuantData, int dataType, int lorenzoMode, int codeFormat, double errorBound, size_t n3, size_t n2,
                                           size_t n1, void *result) {
    int status = 0;
    if (dataType == QCAT_FLOAT) {
        status = lorenzoPredictorQuant_Decmpr_NoOutlier_float(diffQuantData, lorenzoMode, codeFormat, errorBound, n3, n2, n1, result);
    } else if (dataType == QCAT_DOUBLE) {
        status = lorenzoPredictorQuant_Decmpr_NoOutlier_double(diffQuantData, lorenzoMode, codeFormat, errorBound, n3, n2, n1, result);
    }

    return status;
}
