#include <cmath>
using namespace std;

#include "debug.h"
#include "fcns.h"
#include <assert.h>

////////////////////////////////////////////////////////////////////////////////
// SPATIALLY DEPENDENT PDE FUNCTIONS

// Coefficients
static const double PI = 3.141592653589793238462643383279502884197169399375105820974944;
//int n=1, m=0, l=0;

double coefA(double x, double y, double z) {
  return 0;
}

void coefB(double x, double y, double z, Tensor1& b) {
  b.set(0,0,0);
}

void coefC(double x, double y, double z, Tensor1& c) {
  c.set(0,0,0);
}

void coefD(double x, double y, double z, Tensor2& d) {
  //d.set(0,0,0,0,0,0,0,0,0);
  d.set(1,0,0,0,1,0,0,0,1);
}

void coefD_inv(double x, double y, double z, Tensor2& d) {
  //d.set(0,0,0,0,0,0,0,0,0);
  d.set(1,0,0,0,1,0,0,0,1);
}

Tensor2 trueHessianSoln(double x, double y, double z);
// Source f
double sourceVal(double x, double y, double z) {
  //Id D=I, a=b=0
  Tensor2 h = trueHessianSoln(x,y,z);
  return - h(0,0) - h(1,1) - h(2,2);
  //return trueSoln(x,y,z);
}

// BC values g
double bcVal(double x, double y, double z) {
  return trueSoln(x,y,z);
}

////////////////////////////////////////////////////////////////////////////////
// TRUE SOLUTION (IF KNOWN)

bool trueSolnKnown() { return true; }

// Real solution
double trueSoln(double x, double y, double z) {
  //return pow(x,n)*pow(y,m)*pow(z,l);
  return sin(PI*x)*sin(PI*y)*sin(PI*z);
}

Tensor1 trueGradSoln(double x, double y, double z) {
  //double dx = (n==0)? 0 : n*pow(x,n-1)*pow(y,m)*pow(z,l);
  //double dy = (m==0)? 0 : m*pow(x,n)*pow(y,m-1)*pow(z,l);
  //double dz = (l==0)? 0 : l*pow(x,n)*pow(y,m)*pow(z,l-1);
  //return Tensor1(dx,dy,dz);
  return Tensor1(PI*cos(PI*x)*sin(PI*y)*sin(PI*z),
                 PI*sin(PI*x)*cos(PI*y)*sin(PI*z),
                 PI*sin(PI*x)*sin(PI*y)*cos(PI*z));
}

Tensor2 trueHessianSoln(double x, double y, double z) {
  /*
  double dx2 = (n<=1)? 0 : n*(n-1)*pow(x,n-2)*pow(y,m)*pow(z,l);
  double dy2 = (m<=1)? 0 : m*(m-1)*pow(x,n)*pow(y,m-2)*pow(z,l);
  double dz2 = (l<=1)? 0 : l*(l-1)*pow(x,n)*pow(y,m)*pow(z,l-2);
  double dxdy = (n==0 || m==0)? 0 : n*m*pow(x,n-1)*pow(y,m-1)*pow(z,l);
  double dxdz = (n==0 || l==0)? 0 : n*l*pow(x,n-1)*pow(y,m)*pow(z,l-1);
  double dydz = (m==0 || l==0)? 0 : m*l*pow(x,n)*pow(y,m-1)*pow(z,l-1);

  return Tensor2(dx2,dxdy,dxdz,
                 dxdy,dy2,dydz,
                 dxdz,dydz,dz2);  
  */


  return Tensor2(-PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z),  PI*PI*cos(PI*x)*cos(PI*y)*sin(PI*z),  PI*PI*cos(PI*x)*sin(PI*y)*cos(PI*z),
                  PI*PI*cos(PI*x)*cos(PI*y)*sin(PI*z), -PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z),  PI*PI*sin(PI*x)*cos(PI*y)*cos(PI*z),
                  PI*PI*cos(PI*x)*sin(PI*y)*cos(PI*z),  PI*PI*sin(PI*x)*cos(PI*y)*cos(PI*z), -PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z));
}
