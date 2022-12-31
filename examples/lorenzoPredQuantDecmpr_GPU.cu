#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "rw.h"
#include "sz_utility.h"

struct PrivateTimingGPU {
    cudaEvent_t start;
    cudaEvent_t stop;
};

class TimingGPU
{
    private:
        PrivateTimingGPU *privateTimingGPU;

    public:

        TimingGPU();

        ~TimingGPU();

        void StartCounter();
        void StartCounterFlags();

        float GetCounter();

}; // TimingGPU class

// default constructor
TimingGPU::TimingGPU() { privateTimingGPU = new PrivateTimingGPU;  }

// default destructor
TimingGPU::~TimingGPU() { }

void TimingGPU::StartCounter()
{
    cudaEventCreate(&((*privateTimingGPU).start));
    cudaEventCreate(&((*privateTimingGPU).stop));
    cudaEventRecord((*privateTimingGPU).start,0);
}

void TimingGPU::StartCounterFlags()
{
    int eventflags = cudaEventBlockingSync;

    cudaEventCreateWithFlags(&((*privateTimingGPU).start),eventflags);
    cudaEventCreateWithFlags(&((*privateTimingGPU).stop),eventflags);
    cudaEventRecord((*privateTimingGPU).start,0);
}

// Gets the counter in ms
float TimingGPU::GetCounter()
{
    float time;
    cudaEventRecord((*privateTimingGPU).stop, 0);
    cudaEventSynchronize((*privateTimingGPU).stop);
    cudaEventElapsedTime(&time,(*privateTimingGPU).start,(*privateTimingGPU).stop);
    return time;
}


__global__ void kernel_quant_Decmpr_float(float* decData, int* deLorenzoArray, float e2, int bunch)
{
    int index = (threadIdx.x + blockIdx.x * blockDim.x) * bunch;
    int tempIdx;

    for(int i=0; i<bunch; i++)
    {
        tempIdx = index + i;
        decData[tempIdx] = e2 * deLorenzoArray[tempIdx];
    }
}

__global__ void kernel_quant_Decmpr_double(double* decData, int* deLorenzoArray, double e2, int bunch) // yafan reaches here
{
    int index = (threadIdx.x + blockIdx.x * blockDim.x) * bunch;
    int tempIdx;

    for(int i=0; i<bunch; i++)
    {
        tempIdx = index + i;
        decData[tempIdx] = e2 * deLorenzoArray[tempIdx];
    }
}

int lorenzoPredictorQuant_Decmpr_NoOutlier_GPU_float(int* quantData, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, float* result)
{
    TimingGPU timer_GPU;
    size_t nbEle = computeDataLength(0, 0, n3, n2, n1);
    float e2 = errorBound*2;
    int bsize = 256, bunch = 2;
    int gsize = nbEle / (bsize * bunch) + (nbEle % (bsize * bunch) ==0 ? 0 : 1);
    int pad_nbEle = gsize * bsize * bunch;
    int* deLorenzoArray = (int*)malloc(sizeof(int)*nbEle);

    register int curQuantValue = 0;
    register int preQuantValue = 0;
    if(codeFormat == QUANT_CODE_ORIGINAL)
    {
        deLorenzoArray[0] = 0;
        for(size_t i=1; i<nbEle; i++)
        {
            curQuantValue = preQuantValue + quantData[i];
            deLorenzoArray[i] = curQuantValue;
            preQuantValue = curQuantValue;
        }
    }
    else if(codeFormat == QUANT_CODE_NORMALIZE)
    {
        int x = 0;
        deLorenzoArray[0] = 0;
        for(size_t i=1; i<nbEle; i++)
        {
            x = quantData[i];
            curQuantValue = preQuantValue + ((x>>1)^(-(x&1)));
            deLorenzoArray[i] = curQuantValue;
            preQuantValue = curQuantValue;
        }
    }
    else
    {
        printf("Error: wrong quantization_code_mode\n");
        exit(0);
    }

    int* d_deLorenzoArray;
    float* d_result;

    cudaMalloc((void**)&d_deLorenzoArray, sizeof(int)*pad_nbEle);
    cudaMalloc((void**)&d_result, sizeof(float)*pad_nbEle);
    cudaMemcpy(d_deLorenzoArray, deLorenzoArray, sizeof(int)*nbEle, cudaMemcpyHostToDevice);

    dim3 blockSize(bsize);
    dim3 gridSize(gsize);

    timer_GPU.StartCounter(); // set timer
    kernel_quant_Decmpr_float<<<gridSize, blockSize>>>(d_result, d_deLorenzoArray, e2, bunch);
    printf("lorenzoPredQuantDecmpr-float speed: %f GB/s\n", (nbEle*sizeof(float)/1024.0/1024.0)/timer_GPU.GetCounter()); // print speed

    cudaMemcpy(result, d_result, sizeof(float)*nbEle, cudaMemcpyDeviceToHost);

    free(deLorenzoArray);
    cudaFree(d_deLorenzoArray);
    cudaFree(d_result);

    return 0;
}

int lorenzoPredictorQuant_Decmpr_NoOutlier_GPU_double(int* quantData, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, double* result)
{
    TimingGPU timer_GPU;
    size_t nbEle = computeDataLength(0, 0, n3, n2, n1);
    double e2 = errorBound*2;
    int bsize = 256, bunch = 2;
    int gsize = nbEle / (bsize * bunch) + (nbEle % (bsize * bunch) ==0 ? 0 : 1);
    int pad_nbEle = gsize * bsize * bunch;
    int* deLorenzoArray = (int*)malloc(sizeof(int)*nbEle);

    register int curQuantValue = 0;
    register int preQuantValue = 0;
    if(codeFormat == QUANT_CODE_ORIGINAL)
    {
        deLorenzoArray[0] = 0;
        for(size_t i=1; i<nbEle; i++)
        {
            curQuantValue = preQuantValue + quantData[i];
            deLorenzoArray[i] = curQuantValue;
            preQuantValue = curQuantValue;
        }
    }
    else if(codeFormat == QUANT_CODE_NORMALIZE)
    {
        int x = 0;
        deLorenzoArray[0] = 0;
        for(size_t i=1; i<nbEle; i++)
        {
            x = quantData[i];
            curQuantValue = preQuantValue + ((x>>1)^(-(x&1)));
            deLorenzoArray[i] = curQuantValue;
            preQuantValue = curQuantValue;
        }
    }
    else
    {
        printf("Error: wrong quantization_code_mode\n");
        exit(0);
    }

    int* d_deLorenzoArray;
    double* d_result;

    cudaMalloc((void**)&d_deLorenzoArray, sizeof(int)*pad_nbEle);
    cudaMalloc((void**)&d_result, sizeof(double)*pad_nbEle);
    cudaMemcpy(d_deLorenzoArray, deLorenzoArray, sizeof(int)*nbEle, cudaMemcpyHostToDevice);

    dim3 blockSize(bsize);
    dim3 gridSize(gsize);

    timer_GPU.StartCounter(); // set timer
    kernel_quant_Decmpr_double<<<gridSize, blockSize>>>(d_result, d_deLorenzoArray, e2, bunch);
    printf("lorenzoPredQuantDecmpr-double speed: %f GB/s\n", (nbEle*sizeof(double)/1024.0/1024.0)/timer_GPU.GetCounter()); // print speed

    cudaMemcpy(result, d_result, sizeof(double)*nbEle, cudaMemcpyDeviceToHost);

    free(deLorenzoArray);
    cudaFree(d_deLorenzoArray);
    cudaFree(d_result);

    return 0;
}


int main(int argc, char*argv[])
{
    size_t r1 = 0, r2 = 0, r3 = 0;
    int status = 0;
    char oriFilePath[640], outFilePath[645];
    char qmode[30];
    char type[3];
    if(argc < 3)
    {
        printf("Test case: lorenzoPredQuantDecmpr_GPU [type(-f/-d)] [dataFilePath] [quantization_mode] [error_bound] [dims...]\n");
        printf("Example: lorenzoPredQuantDecmpr_GPU -f Hurricane.dat.bin QUANT_CODE_ORIGINAL 1E-2 500 500 100\n");
	    exit(0);
    }

    sprintf(type, "%s", argv[1]);
    sprintf(oriFilePath, "%s", argv[2]);
    sprintf(outFilePath, "%s.f32", oriFilePath);
    sprintf(qmode, "%s", argv[3]);
    double errorBound = atof(argv[4]);

    if(argc>=6)
        r1 = atoi(argv[5]);
    if(argc>=7)
        r2 = atoi(argv[6]);
    if(argc>=8)
        r3 = atoi(argv[7]);

    int mode = 0;
    if(strcmp(qmode, "QUANT_CODE_ORIGINAL")==0)
        mode = QUANT_CODE_ORIGINAL;
    else if(strcmp(qmode, "QUANT_CODE_NORMALIZE")==0)
        mode = QUANT_CODE_NORMALIZE;
    else
    {
        printf("Error: wrong quantization_code_mode\n");
        exit(0);
    }

    size_t nbEle = 0;
    int *quantData = readInt32Data(oriFilePath, &nbEle, &status);
    nbEle = computeDataLength(0, 0, r3, r2, r1); //update the number of elements
    if(strcmp(type, "-f")==0)
    {
        float* out = (float*)malloc(sizeof(float)*nbEle);
        status = lorenzoPredictorQuant_Decmpr_NoOutlier_GPU_float(quantData, mode, errorBound, r3, r2, r1, out);
    	writeFloatData_inBytes(out, nbEle, outFilePath, &status);
        free(out);
    }
    else
    {
        double* out = (double*)malloc(sizeof(double)*nbEle);
        status = lorenzoPredictorQuant_Decmpr_NoOutlier_GPU_double(quantData, mode, errorBound, r3, r2, r1, out);
    	writeDoubleData_inBytes(out, nbEle, outFilePath, &status);
        free(out);
    }
   
    printf("decompressed data are stored in %s\n", outFilePath);

    free(quantData);
    return status;
}