#ifndef __fcns_h_included__
#define __fcns_h_included__

#include "Mesh/baseObjects.h"

////////////////////////////////////////////////////////////////////////////////
// SPATIALLY DEPENDENT PDE FUNCTIONS

double coefA(double x, double y, double z);
void coefB(double x, double y, double z, Tensor1& b);
void coefC(double x, double y, double z, Tensor1& c);
void coefD(double x, double y, double z, Tensor2& d);
void coefD_inv(double x, double y, double z, Tensor2& d);

double sourceVal(double x, double y, double z);
double bcVal(double x, double y, double z);

////////////////////////////////////////////////////////////////////////////////
// TRUE SOLUTION (IF KNOWN)

bool trueSolnKnown();
double trueSoln(double x, double y, double z);
Tensor1 trueGradSoln(double x, double y, double z);
Tensor2 trueHessianSoln(double x, double y, double z);

#endif
