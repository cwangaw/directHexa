#include <cmath>
using namespace std;

#include "debug.h"
#include "fcns.h"
#include <assert.h>

////////////////////////////////////////////////////////////////////////////////
// SPATIALLY DEPENDENT PDE FUNCTIONS

// Coefficients
static const double PI = 3.141592653589793238462643383279502884197169399375105820974944;

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
  d.set(1,0,0,0,1,0,0,0,1);
}

void coefD_inv(double x, double y, double z, Tensor2& d) {
  d.set(1,0,0,0,1,0,0,0,1);
}

Tensor2 trueHessianSoln(double x, double y, double z);
// Source f
double sourceVal(double x, double y, double z) {
  //Id D=I, a=b=0
  Tensor2 h = trueHessianSoln(x,y,z);
  return - h(0,0) - h(1,1) - h(2,2);
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
  return sin(PI*x)*sin(PI*y)*sin(PI*z);
}

Tensor1 trueGradSoln(double x, double y, double z) {
  return Tensor1(PI*cos(PI*x)*sin(PI*y)*sin(PI*z),
                 PI*sin(PI*x)*cos(PI*y)*sin(PI*z),
                 PI*sin(PI*x)*sin(PI*y)*cos(PI*z));
}

Tensor2 trueHessianSoln(double x, double y, double z) {
  return Tensor2(-PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z),  PI*PI*cos(PI*x)*cos(PI*y)*sin(PI*z),  PI*PI*cos(PI*x)*sin(PI*y)*cos(PI*z),
                  PI*PI*cos(PI*x)*cos(PI*y)*sin(PI*z), -PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z),  PI*PI*sin(PI*x)*cos(PI*y)*cos(PI*z),
                  PI*PI*cos(PI*x)*sin(PI*y)*cos(PI*z),  PI*PI*sin(PI*x)*cos(PI*y)*cos(PI*z), -PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z));
}
