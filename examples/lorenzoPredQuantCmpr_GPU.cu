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


__device__ inline int kernel_quantization_float(float data, float recipPrecision)
{
    float dataRecip = data*recipPrecision;
    int s = dataRecip>=-0.5f?0:1;
    return (int)(dataRecip+0.5f) - s;
}

__device__ inline int kernel_quantization_double(double data, double recipPrecision)
{
    double dataRecip = data*recipPrecision;
    int s = dataRecip>=-0.5?0:1;
    return (int)(dataRecip+0.5) - s;
}


__global__ void kernel_lorenzoPredictorQuant_Cmpr_original_float(float* oriData, int* quantArray, float recipPrecision, int bunch)
{
    int index = (threadIdx.x + blockIdx.x * blockDim.x) * bunch;
    int currQuant, pre1Quant;
    int tempIdx;

    // Quantization and 1-layer Lorenzo, original mode.
    pre1Quant = index==0 ? 0 : kernel_quantization_float(oriData[index-1], recipPrecision);
    for(int i=0; i<bunch; i++)
    {
        tempIdx = index + i;
        currQuant = kernel_quantization_float(oriData[tempIdx], recipPrecision);
        quantArray[tempIdx] = currQuant - pre1Quant;
        pre1Quant = currQuant;
    }
}


__global__ void kernel_lorenzoPredictorQuant_Cmpr_original_double(double* oriData, int* quantArray, double recipPrecision, int bunch)
{
    int index = (threadIdx.x + blockIdx.x * blockDim.x) * bunch;
    int currQuant, pre1Quant;
    int tempIdx;

    // Quantization and 1-layer Lorenzo, original mode.
    pre1Quant = index==0 ? 0 : kernel_quantization_double(oriData[index-1], recipPrecision);
    for(int i=0; i<bunch; i++)
    {
        tempIdx = index + i;
        currQuant = kernel_quantization_double(oriData[tempIdx], recipPrecision);
        quantArray[tempIdx] = currQuant - pre1Quant;
        pre1Quant = currQuant;
    }
}


__global__ void kernel_lorenzoPredictorQuant_Cmpr_normalized_float(float* oriData, int* quantArray, float recipPrecision, int bunch)
{
    int index = (threadIdx.x + blockIdx.x * blockDim.x) * bunch;
    int currQuant, pre1Quant, x;
    int tempIdx;

    // Quantization and 1-layer Lorenzo, normalized mode              .
    pre1Quant = index==0 ? 0 : kernel_quantization_float(oriData[index-1], recipPrecision);
    for(int i=0; i<bunch; i++)
    {
        tempIdx = index + i;
        currQuant = kernel_quantization_float(oriData[tempIdx], recipPrecision);
        x = currQuant - pre1Quant;
        quantArray[tempIdx] = (x<<1)^(x>>31);
        pre1Quant = currQuant;
    }
}


__global__ void kernel_lorenzoPredictorQuant_Cmpr_normalized_double(double* oriData, int* quantArray, double recipPrecision, int bunch)
{
    int index = (threadIdx.x + blockIdx.x * blockDim.x) * bunch;
    int currQuant, pre1Quant, x;
    int tempIdx;

    // Quantization and 1-layer Lorenzo, normalized mode              .
    pre1Quant = index==0 ? 0 : kernel_quantization_double(oriData[index-1], recipPrecision);
    for(int i=0; i<bunch; i++)
    {
        tempIdx = index + i;
        currQuant = kernel_quantization_double(oriData[tempIdx], recipPrecision);
        x = currQuant - pre1Quant;
        quantArray[tempIdx] = (x<<1)^(x>>31);
        pre1Quant = currQuant;
    }
}

int lorenzoPredictorQuant_Cmpr_NoOutlier_GPU_float(float* data, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, int* out)
{
    TimingGPU timer_GPU;
    size_t nbEle = computeDataLength(0, 0, n3, n2, n1);
    float recipPrecision = 0.5f/errorBound;
    int bsize = 256, bunch = 2;
    int gsize = nbEle / (bsize * bunch) + (nbEle % (bsize * bunch) ==0 ? 0 : 1);
    int pad_nbEle = gsize * bsize * bunch;
    float* oriData = (float*)malloc(sizeof(float)*pad_nbEle);

    memcpy(oriData, data, sizeof(float)*pad_nbEle);

    float* d_oriData;
    int* d_quantArray;

    cudaMalloc((void**)&d_oriData, sizeof(float)*pad_nbEle);
    cudaMalloc((void**)&d_quantArray, sizeof(int)*pad_nbEle);
    cudaMemcpy(d_oriData, oriData, sizeof(float)*pad_nbEle, cudaMemcpyHostToDevice);

    dim3 blockSize(bsize);
    dim3 gridSize(gsize);

    timer_GPU.StartCounter(); // set timer
    if(codeFormat == QUANT_CODE_ORIGINAL)
        kernel_lorenzoPredictorQuant_Cmpr_original_float<<<gridSize, blockSize>>>(d_oriData, d_quantArray, recipPrecision, bunch);
    else if(codeFormat == QUANT_CODE_NORMALIZE)
        kernel_lorenzoPredictorQuant_Cmpr_normalized_float<<<gridSize, blockSize>>>(d_oriData, d_quantArray, recipPrecision, bunch);
    else
    {   
        // Temporarily leaving blank for other implementation.
    }
    printf("lorenzoPredQuantCmpr-float kernel speed: %f GB/s\n", (nbEle*sizeof(float)/1024.0/1024.0)/timer_GPU.GetCounter()); // print speed

    cudaMemcpy(out, d_quantArray, sizeof(int)*nbEle, cudaMemcpyDeviceToHost);

    free(oriData);
    cudaFree(d_oriData);
    cudaFree(d_quantArray);
    return 0;
}

int lorenzoPredictorQuant_Cmpr_NoOutlier_GPU_double(double* data, int codeFormat, double errorBound, size_t n3, size_t n2, size_t n1, int* out)
{
    TimingGPU timer_GPU;
    size_t nbEle = computeDataLength(0, 0, n3, n2, n1);
    double recipPrecision = 0.5/errorBound;
    int bsize = 256, bunch = 2;
    int gsize = nbEle / (bsize * bunch) + (nbEle % (bsize * bunch) ==0 ? 0 : 1);
    int pad_nbEle = gsize * bsize * bunch;
    double* oriData = (double*)malloc(sizeof(double)*pad_nbEle);

    memcpy(oriData, data, sizeof(double)*pad_nbEle);

    double* d_oriData;
    int* d_quantArray;

    cudaMalloc((void**)&d_oriData, sizeof(double)*pad_nbEle);
    cudaMalloc((void**)&d_quantArray, sizeof(int)*pad_nbEle);
    cudaMemcpy(d_oriData, oriData, sizeof(double)*pad_nbEle, cudaMemcpyHostToDevice);

    dim3 blockSize(bsize);
    dim3 gridSize(gsize);

    timer_GPU.StartCounter(); // set timer
    if(codeFormat == QUANT_CODE_ORIGINAL)
        kernel_lorenzoPredictorQuant_Cmpr_original_double<<<gridSize, blockSize>>>(d_oriData, d_quantArray, recipPrecision, bunch);
    else if(codeFormat == QUANT_CODE_NORMALIZE)
        kernel_lorenzoPredictorQuant_Cmpr_normalized_double<<<gridSize, blockSize>>>(d_oriData, d_quantArray, recipPrecision, bunch);
    else
    {   
        // Temporarily leaving blank for other implementation.
    }
    printf("lorenzoPredQuantCmpr-double kernel speed: %f GB/s\n", (nbEle*sizeof(double)/1024.0/1024.0)/timer_GPU.GetCounter()); // print speed

    cudaMemcpy(out, d_quantArray, sizeof(int)*nbEle, cudaMemcpyDeviceToHost);

    free(oriData);
    cudaFree(d_oriData);
    cudaFree(d_quantArray);
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
        printf("Test case: lorenzoPredQuantCmpr_GPU [type(-f/-d)] [dataFilePath] [quant_mode(QUANT_CODE_ORIGINAL/QUANT_CODE_NORMALIZE] [error_bound] [dims...]\n");
        printf("Example: lorenzoPredQuantCmpr_GPU -f Hurricane.dat QUANT_CODE_ORIGINAL 1E-2 500 500 100\n");
        exit(0);
    }

    sprintf(type, "%s", argv[1]);
    sprintf(oriFilePath, "%s", argv[2]);
    sprintf(outFilePath, "%s.i32", oriFilePath);
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

    size_t nbEle = computeDataLength(0, 0, r3, r2, r1);
    int* out = (int*)malloc(sizeof(int)*nbEle);
    if(strcmp(type, "-f")==0)
    {
        float* data = readFloatData(oriFilePath, &nbEle, &status);
        status = lorenzoPredictorQuant_Cmpr_NoOutlier_GPU_float(data, mode, errorBound, r3, r2, r1, out);
        free(data);
    }
    else if(strcmp(type, "-d")==0)
    {
        double* data = readDoubleData(oriFilePath, &nbEle, &status);
        status = lorenzoPredictorQuant_Cmpr_NoOutlier_GPU_double(data, mode, errorBound, r3, r2, r1, out);
        free(data);
    }
    else
    {
        printf("Error: wrong data type\n");
        exit(0);
    }

    if(status!=0)
    {
        printf("Error state returned by lorenzoPredictorQuant_GPU function.\n");
	    exit(0);
    }


    writeIntData_inBytes(out, nbEle, outFilePath, &status);
    printf("quantization codes are stored in %s\n", outFilePath);
    free(out);

    return status;
}