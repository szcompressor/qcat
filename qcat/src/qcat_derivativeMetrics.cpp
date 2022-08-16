
#include <cmath>
#include <matrix.hpp>
#include "qcat.h"

double qcat_calc_metric_der_order1_sobolev_p2_float(float *data1, float *data2,
  size_t dim, size_t r4, size_t r3, size_t r2, size_t r1, int *status){
  *status=-1;
  
  if(dim<1){
    cout<<"ERROR: Dimension less than 1!"<<endl;
    *status=1;
    return 0;
  }else if(dim==1 && r1<9){
    cout<<"ERROR: dim=1, but r1<9."<<endl;
    cout<<"REMINDER: Derivative order 1 SSIM metrics require every valid dimension to be at least 9."<<endl;
    *status=1;
    return 0;
  }else if(dim==2 && (r1<9 || r2<9)){
    cout<<"ERROR: dim=1, but { r1<9 or r2<9 }."<<endl;
    cout<<"REMINDER: Derivative order 1 SSIM metrics require every valid dimension to be at least 9."<<endl;
    *status=1;
    return 0;
  }else if(dim==3 && (r1<9 || r2<9 || r3<9)){
    cout<<"ERROR: dim=1, but { r1<9 or r2<9 or r3<9 }."<<endl;
    cout<<"REMINDER: Derivative order 1 SSIM metrics require every valid dimension to be at least 9."<<endl;
    *status=1;
    return 0;
  }else if(dim==4 && (r1<9 || r2<9 || r3<9 || r4<9)){
    cout<<"ERROR: dim=1, but { r1<9 or r2<9 or r3<9 or r4<9 }."<<endl;
    cout<<"REMINDER: Derivative order 1 SSIM metrics require every valid dimension to be at least 9."<<endl;
    *status=1;
    return 0;
  }else if(dim>4){
    cout<<"ERROR: dim > 4 is not supported."<<endl;
    *status=1;
    return 0;
  }
  
  if(dim==1){ r2=1;r3=1;r4=1; }
  if(dim==2){ r3=1;r4=1; }
  if(dim==3){ r4=1; }
  
  matrix<float> orig;  orig.nDim=dim; orig.size0=r1; orig.size1=r2; orig.size2=r3; orig.size3=r4;
  matrix<float> lossy;  lossy.nDim=dim; lossy.size0=r1; lossy.size1=r2; lossy.size2=r3; lossy.size3=r4;
  matrix<float> temp1;
  matrix<float> temp2;
  
  orig.data=data1;
  lossy.data=data2;
  
  double result=orig.sobolevNorm_s1_p2(lossy);
  
  orig.data=NULL; //To prevent deallocation from ~matrix()
  lossy.data=NULL; //To prevent deallocation from ~matrix()
  
  *status=0;
  return result;  

}

double qcat_calc_metric_der_order1_sobolev_p2_double(double *data1, double *data2,
  size_t dim, size_t r4, size_t r3, size_t r2, size_t r1, int *status){
  *status=-1;
  
  if(dim<1){
    cout<<"ERROR: Dimension less than 1!"<<endl;
    *status=1;
    return 0;
  }else if(dim==1 && r1<9){
    cout<<"ERROR: dim=1, but r1<9."<<endl;
    cout<<"REMINDER: Derivative order 1 SSIM metrics require every valid dimension to be at least 9."<<endl;
    *status=1;
    return 0;
  }else if(dim==2 && (r1<9 || r2<9)){
    cout<<"ERROR: dim=1, but { r1<9 or r2<9 }."<<endl;
    cout<<"REMINDER: Derivative order 1 SSIM metrics require every valid dimension to be at least 9."<<endl;
    *status=1;
    return 0;
  }else if(dim==3 && (r1<9 || r2<9 || r3<9)){
    cout<<"ERROR: dim=1, but { r1<9 or r2<9 or r3<9 }."<<endl;
    cout<<"REMINDER: Derivative order 1 SSIM metrics require every valid dimension to be at least 9."<<endl;
    *status=1;
    return 0;
  }else if(dim==4 && (r1<9 || r2<9 || r3<9 || r4<9)){
    cout<<"ERROR: dim=1, but { r1<9 or r2<9 or r3<9 or r4<9 }."<<endl;
    cout<<"REMINDER: Derivative order 1 SSIM metrics require every valid dimension to be at least 9."<<endl;
    *status=1;
    return 0;
  }else if(dim>4){
    cout<<"ERROR: dim > 4 is not supported."<<endl;
    *status=1;
    return 0;
  }
  
  if(dim==1){ r2=1;r3=1;r4=1; }
  if(dim==2){ r3=1;r4=1; }
  if(dim==3){ r4=1; }
  
  matrix<double> orig;  orig.nDim=dim; orig.size0=r1; orig.size1=r2; orig.size2=r3; orig.size3=r4;
  matrix<double> lossy;  lossy.nDim=dim; lossy.size0=r1; lossy.size1=r2; lossy.size2=r3; lossy.size3=r4;
  matrix<double> temp1;
  matrix<double> temp2;
  
  orig.data=data1;
  lossy.data=data2;
  
  double result=orig.sobolevNorm_s1_p2(lossy);
  
  orig.data=NULL; //To prevent deallocation from ~matrix()
  lossy.data=NULL; //To prevent deallocation from ~matrix()
  
  *status=0;
  return result;  

}

