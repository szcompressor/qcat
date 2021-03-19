#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <qcat_ssim.h>

/////////////////// 1D

double SSIM_1d_calcWindow_float(float* data, float* other, int offset0, int windowSize0)
{
	int i0;
	int np=0; //Number of points
	
	float xMin=data[0];
	float xMax=data[0];
	float yMin=data[0];
	float yMax=data[0];
	double xSum=0;
	double x2Sum=0;
	double ySum=0;
	double y2Sum=0;
	double xySum=0;

	for(i0=offset0;i0<offset0+windowSize0;i0++){
		np++;
		if(xMin>data[i0])
			xMin=data[i0];
		if(xMax<data[i0])
			xMax=data[i0];
		if(yMin>other[i0])
			yMin=other[i0];
		if(yMax<other[i0])
			yMax=other[i0];
		xSum+=data[i0];
		x2Sum+=(data[i0]*data[i0]);
		ySum+=other[i0];
		y2Sum+=(other[i0]*other[i0]);
		xySum+=(data[i0]*other[i0]);
	}


	double xMean=xSum/np;
	double yMean=ySum/np;
	double xSigma=sqrt((x2Sum/np)-(xMean*xMean));
	double ySigma=sqrt((y2Sum/np)-(yMean*yMean));
	double xyCov=(xySum/np)-(xMean*yMean);

	double c1,c2;
	if(xMax-xMin==0){
		c1=K1*K1;
		c2=K2*K2;
	}else{
		c1=K1*K1*(xMax-xMin)*(xMax-xMin);
		c2=K2*K2*(xMax-xMin)*(xMax-xMin);
	}
	double c3=c2/2;

	double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
	double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
	double structure=(xyCov+c3)/(xSigma*ySigma+c3);
	double ssim=luminance*contrast*structure;
	return ssim;	
}

double SSIM_1d_calcWindow_double(double* oriData, double* decData, int offset0, int windowSize0)
{
	int i0;
	int np=0; //Number of points
	
	double* data = oriData;
	double* other = decData;
	
	double xMin=data[0];
	double xMax=data[0];
	double yMin=data[0];
	double yMax=data[0];
	double xSum=0;
	double x2Sum=0;
	double ySum=0;
	double y2Sum=0;
	double xySum=0;

	for(i0=offset0;i0<offset0+windowSize0;i0++){
		np++;
		if(xMin>data[i0])
			xMin=data[i0];
		if(xMax<data[i0])
			xMax=data[i0];
		if(yMin>other[i0])
			yMin=other[i0];
		if(yMax<other[i0])
			yMax=other[i0];
		xSum+=data[i0];
		x2Sum+=(data[i0]*data[i0]);
		ySum+=other[i0];
		y2Sum+=(other[i0]*other[i0]);
		xySum+=(data[i0]*other[i0]);
	}


	double xMean=xSum/np;
	double yMean=ySum/np;
	double xSigma=sqrt((x2Sum/np)-(xMean*xMean));
	double ySigma=sqrt((y2Sum/np)-(yMean*yMean));
	double xyCov=(xySum/np)-(xMean*yMean);

	double c1,c2;
	if(xMax-xMin==0){
		c1=K1*K1;
		c2=K2*K2;
	}else{
		c1=K1*K1*(xMax-xMin)*(xMax-xMin);
		c2=K2*K2*(xMax-xMin)*(xMax-xMin);
	}
	double c3=c2/2;

	double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
	double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
	double structure=(xyCov+c3)/(xSigma*ySigma+c3);
	double ssim=luminance*contrast*structure;
	return ssim;	
}

double SSIM_1d_windowed_float(float* oriData, float* decData, size_t size0, int windowSize0, int windowShift0){
	int offset0;
	int nw=0; //Number of windows
	double ssimSum=0;
	int offsetInc0;

	if(windowSize0>size0){printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);}

	//offsetInc0=windowSize0/2;
	offsetInc0=windowShift0;

  
	for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0){ //MOVING WINDOW
		nw++;
		ssimSum+=SSIM_1d_calcWindow_float(oriData, decData, offset0, windowSize0);
	}
  
	return ssimSum/nw;
}

double SSIM_1d_windowed_double(double* oriData, double* decData, size_t size0, int windowSize0, int windowShift0)
{
	int offset0;
	int nw=0; //Number of windows
	double ssimSum=0;
	int offsetInc0;

	if(windowSize0>size0){printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);}

	//offsetInc0=windowSize0/2;
	offsetInc0=windowShift0;

  
	for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0){ //MOVING WINDOW
		nw++;
		ssimSum+=SSIM_1d_calcWindow_double(oriData, decData, offset0, windowSize0);
	}
  
	return ssimSum/nw;
}

//////////////////// 2D

double SSIM_2d_windowed_float(float* oriData, float* decData, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowShift0, int windowShift1)
{
	int offset0,offset1;
	int nw=0; //Number of windows
	double ssimSum=0;
	int offsetInc0,offsetInc1;

	float* data = oriData;
	float* other = decData;

	if(windowSize0>size0){printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);}
	if(windowSize1>size1){printf("ERROR: windowSize1 = %d > %zu\n", windowSize1, size1);}

	//offsetInc0=windowSize0/2;
	//offsetInc1=windowSize1/2;
	offsetInc0=windowShift0;
	offsetInc1=windowShift1;
  
	for(offset1=0; offset1+windowSize1<=size1; offset1+=offsetInc1){ //MOVING WINDOW

		for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0){ //MOVING WINDOW
			nw++;
			ssimSum+=SSIM_2d_calcWindow_float(data, other, size0, offset0, offset1, windowSize0, windowSize1);
		}
	}

	return ssimSum/nw;
}

double SSIM_2d_calcWindow_float(float* data, float *other, size_t size0, int offset0, int offset1, int windowSize0, int windowSize1)
{
	int i0,i1,index;
	int np=0; //Number of points
	float xMin=data[0];
	float xMax=data[0];
	float yMin=data[0];
	float yMax=data[0];
	double xSum=0;
	double x2Sum=0;
	double ySum=0;
	double y2Sum=0;
	double xySum=0;

	for(i1=offset1;i1<offset1+windowSize1;i1++){
	for(i0=offset0;i0<offset0+windowSize0;i0++){
	  np++;
	  index=i0+size0*i1;
	  if(xMin>data[index])
		xMin=data[index];
	  if(xMax<data[index])
		xMax=data[index];
	  if(yMin>other[index])
		yMin=other[index];
	  if(yMax<other[index])
		yMax=other[index];
	  xSum+=data[index];
	  x2Sum+=(data[index]*data[index]);
	  ySum+=other[index];
	  y2Sum+=(other[index]*other[index]);
	  xySum+=(data[index]*other[index]);
	}
	}

	double xMean=xSum/np;
	double yMean=ySum/np;
	double xSigma=sqrt((x2Sum/np)-(xMean*xMean));
	double ySigma=sqrt((y2Sum/np)-(yMean*yMean));
	double xyCov=(xySum/np)-(xMean*yMean);

	double c1,c2;
	if(xMax-xMin==0){
	c1=K1*K1;
	c2=K2*K2;
	}else{
	c1=K1*K1*(xMax-xMin)*(xMax-xMin);
	c2=K2*K2*(xMax-xMin)*(xMax-xMin);
	}
	double c3=c2/2;

	double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
	double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
	double structure=(xyCov+c3)/(xSigma*ySigma+c3);
	double ssim=luminance*contrast*structure;
	return ssim;
}

double SSIM_2d_windowed_double(double* oriData, double* decData, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowShift0, int windowShift1)
{
	int offset0,offset1;
	int nw=0; //Number of windows
	double ssimSum=0;
	int offsetInc0,offsetInc1;

	double* data = oriData;
	double* other = decData;

	if(windowSize0>size0){printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);}
	if(windowSize1>size1){printf("ERROR: windowSize1 = %d > %zu\n", windowSize1, size1);}

	//offsetInc0=windowSize0/2;
	//offsetInc1=windowSize1/2;
	offsetInc0=windowShift0;
	offsetInc1=windowShift1;
  
	for(offset1=0; offset1+windowSize1<=size1; offset1+=offsetInc1){ //MOVING WINDOW

		for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0){ //MOVING WINDOW
			nw++;
			ssimSum+=SSIM_2d_calcWindow_double(data, other, size0, offset0, offset1, windowSize0, windowSize1);
		}
	}

	return ssimSum/nw;
}

double SSIM_2d_calcWindow_double(double* data, double *other, size_t size0, int offset0, int offset1, int windowSize0, int windowSize1){
	int i0,i1,index;
	int np=0; //Number of points
	double xMin=data[0];
	double xMax=data[0];
	double yMin=data[0];
	double yMax=data[0];
	double xSum=0;
	double x2Sum=0;
	double ySum=0;
	double y2Sum=0;
	double xySum=0;

	for(i1=offset1;i1<offset1+windowSize1;i1++){
	for(i0=offset0;i0<offset0+windowSize0;i0++){
	  np++;
	  index=i0+size0*i1;
	  if(xMin>data[index])
		xMin=data[index];
	  if(xMax<data[index])
		xMax=data[index];
	  if(yMin>other[index])
		yMin=other[index];
	  if(yMax<other[index])
		yMax=other[index];
	  xSum+=data[index];
	  x2Sum+=(data[index]*data[index]);
	  ySum+=other[index];
	  y2Sum+=(other[index]*other[index]);
	  xySum+=(data[index]*other[index]);
	}
	}

	double xMean=xSum/np;
	double yMean=ySum/np;
	double xSigma=sqrt((x2Sum/np)-(xMean*xMean));
	double ySigma=sqrt((y2Sum/np)-(yMean*yMean));
	double xyCov=(xySum/np)-(xMean*yMean);

	double c1,c2;
	if(xMax-xMin==0){
	c1=K1*K1;
	c2=K2*K2;
	}else{
	c1=K1*K1*(xMax-xMin)*(xMax-xMin);
	c2=K2*K2*(xMax-xMin)*(xMax-xMin);
	}
	double c3=c2/2;

	double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
	double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
	double structure=(xyCov+c3)/(xSigma*ySigma+c3);
	double ssim=luminance*contrast*structure;
	return ssim;
}


//////////////////////// 3D

double SSIM_3d_windowed_float(float* oriData, float* decData, size_t size2, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowSize2, int windowShift0, int windowShift1, int windowShift2)
{
	int offset0,offset1,offset2;
	int nw=0; //Number of windows
	double ssimSum=0;
	int offsetInc0,offsetInc1,offsetInc2;

	if(windowSize0>size0){printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);}
	if(windowSize1>size1){printf("ERROR: windowSize1 = %d > %zu\n", windowSize1, size1);}
	if(windowSize2>size2){printf("ERROR: windowSize2 = %d > %zu\n", windowSize2, size2);}

	//offsetInc0=windowSize0/2;
	//offsetInc1=windowSize1/2;
	//offsetInc2=windowSize2/2;
	offsetInc0=windowShift0;
	offsetInc1=windowShift1;
	offsetInc2=windowShift2;

	for(offset2=0; offset2+windowSize2<=size2; offset2+=offsetInc2){ //MOVING WINDOW
	  
	for(offset1=0; offset1+windowSize1<=size1; offset1+=offsetInc1){ //MOVING WINDOW
	  
	  for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0){ //MOVING WINDOW
		nw++;
		ssimSum+=SSIM_3d_calcWindow_float(oriData, decData, size1, size0, offset0, offset1, offset2, windowSize0, windowSize1, windowSize2);
		
	  }
	}
	}

	return ssimSum/nw;
}

double SSIM_3d_calcWindow_float(float* data, float* other, size_t size1, size_t size0, int offset0, int offset1, int offset2, int windowSize0, int windowSize1, int windowSize2){
	int i0,i1,i2,index;
	int np=0; //Number of points
	float xMin=data[0];
	float xMax=data[0];
	float yMin=data[0];
	float yMax=data[0];
	double xSum=0;
	double x2Sum=0;
	double ySum=0;
	double y2Sum=0;
	double xySum=0;

	for(i2=offset2;i2<offset2+windowSize2;i2++){
	for(i1=offset1;i1<offset1+windowSize1;i1++){
	  for(i0=offset0;i0<offset0+windowSize0;i0++){
		np++;
		index=i0+size0*(i1+size1*i2);
		if(xMin>data[index])
		  xMin=data[index];
		if(xMax<data[index])
		  xMax=data[index];
		if(yMin>other[index])
		  yMin=other[index];
		if(yMax<other[index])
		  yMax=other[index];
		xSum+=data[index];
		x2Sum+=(data[index]*data[index]);
		ySum+=other[index];
		y2Sum+=(other[index]*other[index]);
		xySum+=(data[index]*other[index]);
	  }
	}
	}


	double xMean=xSum/np;
	double yMean=ySum/np;
	double xSigma=sqrt(fabs((x2Sum/np)-(xMean*xMean)));
	double ySigma=sqrt(fabs((y2Sum/np)-(yMean*yMean)));
	double xyCov=(xySum/np)-(xMean*yMean);


	double c1,c2;
	if(xMax-xMin==0){
	c1=K1*K1;
	c2=K2*K2;
	}else{
	c1=K1*K1*(xMax-xMin)*(xMax-xMin);
	c2=K2*K2*(xMax-xMin)*(xMax-xMin);
	}
	double c3=c2/2;

	double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
	double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
	double structure=(xyCov+c3)/(xSigma*ySigma+c3);
	double ssim=luminance*contrast*structure;

	return ssim;
}


double SSIM_3d_windowed_double(double* oriData, double* decData, size_t size2, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowSize2, int windowShift0, int windowShift1, int windowShift2)
{
	int offset0,offset1,offset2;
	int nw=0; //Number of windows
	double ssimSum=0;
	int offsetInc0,offsetInc1,offsetInc2;

	if(windowSize0>size0){printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);}
	if(windowSize1>size1){printf("ERROR: windowSize1 = %d > %zu\n", windowSize1, size1);}
	if(windowSize2>size2){printf("ERROR: windowSize2 = %d > %zu\n", windowSize2, size2);}

	//offsetInc0=windowSize0/2;
	//offsetInc1=windowSize1/2;
	//offsetInc2=windowSize2/2;
	offsetInc0=windowShift0;
	offsetInc1=windowShift1;
	offsetInc2=windowShift2;

	for(offset2=0; offset2+windowSize2<=size2; offset2+=offsetInc2){ //MOVING WINDOW
	  
	for(offset1=0; offset1+windowSize1<=size1; offset1+=offsetInc1){ //MOVING WINDOW
	  
	  for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0){ //MOVING WINDOW
		nw++;
		ssimSum+=SSIM_3d_calcWindow_double(oriData, decData, size1, size0, offset0, offset1, offset2, windowSize0, windowSize1, windowSize2);
		
	  }
	}
	}

	return ssimSum/nw;
}

double SSIM_3d_calcWindow_double(double* data, double* other, size_t size1, size_t size0, int offset0, int offset1, int offset2, int windowSize0, int windowSize1, int windowSize2)
{
	int i0,i1,i2,index;
	int np=0; //Number of points
	double xMin=data[0];
	double xMax=data[0];
	double yMin=data[0];
	double yMax=data[0];
	double xSum=0;
	double x2Sum=0;
	double ySum=0;
	double y2Sum=0;
	double xySum=0;

	for(i2=offset2;i2<offset2+windowSize2;i2++){
	for(i1=offset1;i1<offset1+windowSize1;i1++){
	  for(i0=offset0;i0<offset0+windowSize0;i0++){
		np++;
		index=i0+size0*(i1+size1*i2);
		if(xMin>data[index])
		  xMin=data[index];
		if(xMax<data[index])
		  xMax=data[index];
		if(yMin>other[index])
		  yMin=other[index];
		if(yMax<other[index])
		  yMax=other[index];
		xSum+=data[index];
		x2Sum+=(data[index]*data[index]);
		ySum+=other[index];
		y2Sum+=(other[index]*other[index]);
		xySum+=(data[index]*other[index]);
	  }
	}
	}


	double xMean=xSum/np;
	double yMean=ySum/np;
	double xSigma=sqrt(fabs((x2Sum/np)-(xMean*xMean)));
	double ySigma=sqrt(fabs((y2Sum/np)-(yMean*yMean)));
	double xyCov=(xySum/np)-(xMean*yMean);


	double c1,c2;
	if(xMax-xMin==0){
	c1=K1*K1;
	c2=K2*K2;
	}else{
	c1=K1*K1*(xMax-xMin)*(xMax-xMin);
	c2=K2*K2*(xMax-xMin)*(xMax-xMin);
	}
	double c3=c2/2;

	double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
	double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
	double structure=(xyCov+c3)/(xSigma*ySigma+c3);
	double ssim=luminance*contrast*structure;

	return ssim;
}

double SSIM_4d_windowed_float(float* oriData, float* decData, size_t size3, size_t size2, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowSize2, int windowSize3, int windowShift0, int windowShift1, int windowShift2, int windowShift3)
{
	int offset0,offset1,offset2,offset3;
	int nw=0; //Number of windows
	double ssimSum=0;
	int offsetInc0,offsetInc1,offsetInc2,offsetInc3;

	if(windowSize0>size0){printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);}
	if(windowSize1>size1){printf("ERROR: windowSize1 = %d > %zu\n", windowSize1, size1);}
	if(windowSize2>size2){printf("ERROR: windowSize2 = %d > %zu\n", windowSize2, size2);}
	if(windowSize3>size3){printf("ERROR: windowSize3 = %d > %zu\n", windowSize3, size3);}

	//offsetInc0=windowSize0/2;
	//offsetInc1=windowSize1/2;
	//offsetInc2=windowSize2/2;
	//offsetInc3=windowSize3/2;
	offsetInc0=windowShift0;
	offsetInc1=windowShift1;
	offsetInc2=windowShift2;
	offsetInc3=windowShift3;

	for(offset3=0; offset3+windowSize3<=size3; offset3+=offsetInc3){ //MOVING WINDOW

		for(offset2=0; offset2+windowSize2<=size2; offset2+=offsetInc2){ //MOVING WINDOW
			
		  for(offset1=0; offset1+windowSize1<=size1; offset1+=offsetInc1){ //MOVING WINDOW
			
			for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0){ //MOVING WINDOW
			  nw++;
			  ssimSum+=SSIM_4d_calcWindow_float(oriData, decData, size2, size1, size0, offset0, offset1, offset2, offset3, windowSize0, windowSize1, windowSize2, windowSize3);
			}
		  }
		}
	}

	return ssimSum/nw;
	return 0;
}

double SSIM_4d_calcWindow_float(float* data, float* other, size_t size2, size_t size1, size_t size0, int offset0, int offset1, int offset2, int offset3,int windowSize0, int windowSize1, int windowSize2, int windowSize3)
{
	int i0,i1,i2,i3,index;
	int np=0; //Number of points
	float xMin=data[0];
	float xMax=data[0];
	float yMin=data[0];
	float yMax=data[0];
	double xSum=0;
	double x2Sum=0;
	double ySum=0;
	double y2Sum=0;
	double xySum=0;
	for(i3=offset3;i3<offset3+windowSize3;i3++){
	for(i2=offset2;i2<offset2+windowSize2;i2++){
	  for(i1=offset1;i1<offset1+windowSize1;i1++){
		for(i0=offset0;i0<offset0+windowSize0;i0++){
		  np++;
		  index=i0+size0*(i1+size1*(i2+size2*i3));
		  if(xMin>data[index])
			xMin=data[index];
		  if(xMax<data[index])
			xMax=data[index];
		  if(yMin>other[index])
			yMin=other[index];
		  if(yMax<other[index])
			yMax=other[index];
		  xSum+=data[index];
		  x2Sum+=(data[index]*data[index]);
		  ySum+=other[index];
		  y2Sum+=(other[index]*other[index]);
		  xySum+=(data[index]*other[index]);
		}
	  }
	}
	}

	double xMean=xSum/np;
	double yMean=ySum/np;
	double xSigma=sqrt((x2Sum/np)-(xMean*xMean));
	double ySigma=sqrt((y2Sum/np)-(yMean*yMean));
	double xyCov=(xySum/np)-(xMean*yMean);

	double c1,c2;
	if(xMax-xMin==0){
	c1=K1*K1;
	c2=K2*K2;
	}else{
	c1=K1*K1*(xMax-xMin)*(xMax-xMin);
	c2=K2*K2*(xMax-xMin)*(xMax-xMin);
	}
	double c3=c2/2;

	double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
	double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
	double structure=(xyCov+c3)/(xSigma*ySigma+c3);
	double ssim=luminance*contrast*structure;
	return ssim;
}

double SSIM_4d_windowed_double(double* oriData, double* decData, size_t size3, size_t size2, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowSize2, int windowSize3, int windowShift0, int windowShift1, int windowShift2, int windowShift3)
{
	int offset0,offset1,offset2,offset3;
	int nw=0; //Number of windows
	double ssimSum=0;
	int offsetInc0,offsetInc1,offsetInc2,offsetInc3;

	if(windowSize0>size0){printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);}
	if(windowSize1>size1){printf("ERROR: windowSize1 = %d > %zu\n", windowSize1, size1);}
	if(windowSize2>size2){printf("ERROR: windowSize2 = %d > %zu\n", windowSize2, size2);}
	if(windowSize3>size3){printf("ERROR: windowSize3 = %d > %zu\n", windowSize3, size3);}

	//offsetInc0=windowSize0/2;
	//offsetInc1=windowSize1/2;
	//offsetInc2=windowSize2/2;
	//offsetInc3=windowSize3/2;
	offsetInc0=windowShift0;
	offsetInc1=windowShift1;
	offsetInc2=windowShift2;
	offsetInc3=windowShift3;

	for(offset3=0; offset3+windowSize3<=size3; offset3+=offsetInc3){ //MOVING WINDOW

	for(offset2=0; offset2+windowSize2<=size2; offset2+=offsetInc2){ //MOVING WINDOW
		
	  for(offset1=0; offset1+windowSize1<=size1; offset1+=offsetInc1){ //MOVING WINDOW
		
		for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0){ //MOVING WINDOW
		  nw++;
		  ssimSum+=SSIM_4d_calcWindow_double(oriData, decData, size2, size1, size0, offset0, offset1, offset2, offset3, windowSize0, windowSize1, windowSize2, windowSize3);
		}
      }
    }
  }
  
  return ssimSum/nw;
  return 0;
}

double SSIM_4d_calcWindow_double(double* data, double* other, size_t size2, size_t size1, size_t size0, int offset0, int offset1, int offset2, int offset3,int windowSize0, int windowSize1, int windowSize2, int windowSize3)
{
  int i0,i1,i2,i3,index;
  int np=0; //Number of points
  double xMin=data[0];
  double xMax=data[0];
  double yMin=data[0];
  double yMax=data[0];
  double xSum=0;
  double x2Sum=0;
  double ySum=0;
  double y2Sum=0;
  double xySum=0;
  for(i3=offset3;i3<offset3+windowSize3;i3++){
    for(i2=offset2;i2<offset2+windowSize2;i2++){
      for(i1=offset1;i1<offset1+windowSize1;i1++){
        for(i0=offset0;i0<offset0+windowSize0;i0++){
          np++;
          index=i0+size0*(i1+size1*(i2+size2*i3));
          if(xMin>data[index])
            xMin=data[index];
          if(xMax<data[index])
            xMax=data[index];
          if(yMin>other[index])
            yMin=other[index];
          if(yMax<other[index])
            yMax=other[index];
          xSum+=data[index];
          x2Sum+=(data[index]*data[index]);
          ySum+=other[index];
          y2Sum+=(other[index]*other[index]);
          xySum+=(data[index]*other[index]);
        }
      }
    }
  }

  double xMean=xSum/np;
  double yMean=ySum/np;
  double xSigma=sqrt((x2Sum/np)-(xMean*xMean));
  double ySigma=sqrt((y2Sum/np)-(yMean*yMean));
  double xyCov=(xySum/np)-(xMean*yMean);
  
  double c1,c2;
  if(xMax-xMin==0){
    c1=K1*K1;
    c2=K2*K2;
  }else{
    c1=K1*K1*(xMax-xMin)*(xMax-xMin);
    c2=K2*K2*(xMax-xMin)*(xMax-xMin);
  }
  double c3=c2/2;
    
  double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
  double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
  double structure=(xyCov+c3)/(xSigma*ySigma+c3);
  double ssim=luminance*contrast*structure;
  return ssim;
}
