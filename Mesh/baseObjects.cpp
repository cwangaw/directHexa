#include "baseObjects.h"
#include "debug.h"


////////////////////////////////////////////////////////////////////////////////
// Class Point

Point operator+(const Point& p0, const Point& p1) {
  Point p(p0); p += p1; return p;
};

Point operator-(const Point& p0, const Point& p1) {
  Point p(p0); p -= p1; return p;
};

Point operator*(double scalar, const Point& p0) {
  Point p(p0); p *= scalar; return p;
};

Point operator*(const Point& p0, double scalar) {
  Point p(p0); p *= scalar; return p;
};

Point operator/(const Point& p0, double scalar) {
  Point p(p0); p /= scalar; return p;
};

////////////////////////////////////////////////////////////////////////////////
// Class Tensor2

static const double pi = 3.141592653589793238463;

void Tensor2::mult(double scalar, Tensor2& result) const {
  result.the_tensor[0][0] = scalar*the_tensor[0][0]; result.the_tensor[0][1] = scalar*the_tensor[0][1]; result.the_tensor[0][2] = scalar*the_tensor[0][2]; 
  result.the_tensor[1][0] = scalar*the_tensor[1][0]; result.the_tensor[1][1] = scalar*the_tensor[1][1]; result.the_tensor[1][2] = scalar*the_tensor[1][2];
  result.the_tensor[2][0] = scalar*the_tensor[2][0]; result.the_tensor[2][1] = scalar*the_tensor[2][1]; result.the_tensor[2][2] = scalar*the_tensor[2][2];
}

void Tensor2::mult(const Tensor1& v, Tensor1& result) const {
  result[0] = the_tensor[0][0]*v[0] + the_tensor[0][1]*v[1] + the_tensor[0][2]*v[2];
  result[1] = the_tensor[1][0]*v[0] + the_tensor[1][1]*v[1] + the_tensor[1][2]*v[2];
  result[2] = the_tensor[2][0]*v[0] + the_tensor[2][1]*v[1] + the_tensor[2][2]*v[2];
}

void Tensor2::mult(const Point& v, Point& result) const {
  result[0] = the_tensor[0][0]*v[0] + the_tensor[0][1]*v[1] + the_tensor[0][2]*v[2];
  result[1] = the_tensor[1][0]*v[0] + the_tensor[1][1]*v[1] + the_tensor[1][2]*v[2];
  result[2] = the_tensor[2][0]*v[0] + the_tensor[2][1]*v[1] + the_tensor[2][2]*v[2];
}

void Tensor2::mult(const Tensor2& t, Tensor2& result) const {
  result.the_tensor[0][0] = the_tensor[0][0]*t.the_tensor[0][0] + the_tensor[0][1]*t.the_tensor[1][0] + the_tensor[0][2]*t.the_tensor[2][0];
  result.the_tensor[0][1] = the_tensor[0][0]*t.the_tensor[0][1] + the_tensor[0][1]*t.the_tensor[1][1] + the_tensor[0][2]*t.the_tensor[2][1];
  result.the_tensor[0][2] = the_tensor[0][0]*t.the_tensor[0][2] + the_tensor[0][1]*t.the_tensor[1][2] + the_tensor[0][2]*t.the_tensor[2][2];
  result.the_tensor[1][0] = the_tensor[1][0]*t.the_tensor[0][0] + the_tensor[1][1]*t.the_tensor[1][0] + the_tensor[1][2]*t.the_tensor[2][0];
  result.the_tensor[1][1] = the_tensor[1][0]*t.the_tensor[0][1] + the_tensor[1][1]*t.the_tensor[1][1] + the_tensor[1][2]*t.the_tensor[2][1];
  result.the_tensor[1][2] = the_tensor[1][0]*t.the_tensor[0][2] + the_tensor[1][1]*t.the_tensor[1][2] + the_tensor[1][2]*t.the_tensor[2][2];
  result.the_tensor[2][0] = the_tensor[2][0]*t.the_tensor[0][0] + the_tensor[2][1]*t.the_tensor[1][0] + the_tensor[2][2]*t.the_tensor[2][0];
  result.the_tensor[2][1] = the_tensor[2][0]*t.the_tensor[0][1] + the_tensor[2][1]*t.the_tensor[1][1] + the_tensor[2][2]*t.the_tensor[2][1];
  result.the_tensor[2][2] = the_tensor[2][0]*t.the_tensor[0][2] + the_tensor[2][1]*t.the_tensor[1][2] + the_tensor[2][2]*t.the_tensor[2][2];
}

Tensor2& Tensor2::operator*=(const Tensor2& t) {
  double t00 = the_tensor[0][0], t01 = the_tensor[0][1], t02 = the_tensor[0][2]; 
  the_tensor[0][0] = t00*t.the_tensor[0][0] + t01*t.the_tensor[1][0] + t02*t.the_tensor[2][0];
  the_tensor[0][1] = t00*t.the_tensor[0][1] + t01*t.the_tensor[1][1] + t02*t.the_tensor[2][1];
  the_tensor[0][2] = t00*t.the_tensor[0][2] + t01*t.the_tensor[1][2] + t02*t.the_tensor[2][2];
  double t10 = the_tensor[1][0], t11 = the_tensor[1][1], t12 = the_tensor[1][2];
  the_tensor[1][0] = t10*t.the_tensor[0][0] + t11*t.the_tensor[1][0] + t12*t.the_tensor[2][0];
  the_tensor[1][1] = t10*t.the_tensor[0][1] + t11*t.the_tensor[1][1] + t12*t.the_tensor[2][1];
  the_tensor[1][2] = t10*t.the_tensor[0][2] + t11*t.the_tensor[1][2] + t12*t.the_tensor[2][2];
  double t20 = the_tensor[2][0], t21 = the_tensor[2][1], t22 = the_tensor[2][2];
  the_tensor[2][0] = t20*t.the_tensor[0][0] + t21*t.the_tensor[1][0] + t22*t.the_tensor[2][0];
  the_tensor[2][1] = t20*t.the_tensor[0][1] + t21*t.the_tensor[1][1] + t22*t.the_tensor[2][1];
  the_tensor[2][2] = t20*t.the_tensor[0][2] + t21*t.the_tensor[1][2] + t22*t.the_tensor[2][2];
  return *this;
}

double Tensor2::determinant() const {
  return the_tensor[0][0] * ( the_tensor[1][1] * the_tensor[2][2] - the_tensor[2][1] * the_tensor[1][2] )
       - the_tensor[1][0] * ( the_tensor[0][1] * the_tensor[2][2] - the_tensor[2][1] * the_tensor[0][2] )
       + the_tensor[2][0] * ( the_tensor[0][1] * the_tensor[1][2] - the_tensor[1][1] * the_tensor[0][2] );
}

bool Tensor2::inverse(Tensor2& invT) const {
  double det = determinant();
  if(det == 0) return false;

  invT(0,0) =  (the_tensor[1][1] * the_tensor[2][2] - the_tensor[1][2] * the_tensor[2][1])/det;
  invT(0,1) = -(the_tensor[0][1] * the_tensor[2][2] - the_tensor[2][1] * the_tensor[0][2])/det;
  invT(0,2) =  (the_tensor[0][1] * the_tensor[1][2] - the_tensor[0][2] * the_tensor[1][1])/det;
  invT(1,0) = -(the_tensor[1][0] * the_tensor[2][2] - the_tensor[1][2] * the_tensor[2][0])/det;
  invT(1,1) =  (the_tensor[0][0] * the_tensor[2][2] - the_tensor[0][2] * the_tensor[2][0])/det;
  invT(1,2) = -(the_tensor[0][0] * the_tensor[1][2] - the_tensor[0][2] * the_tensor[1][0])/det;
  invT(2,0) =  (the_tensor[1][0] * the_tensor[2][1] - the_tensor[1][1] * the_tensor[2][0])/det;
  invT(2,1) = -(the_tensor[0][0] * the_tensor[2][1] - the_tensor[0][1] * the_tensor[2][0])/det;
  invT(2,2) =  (the_tensor[0][0] * the_tensor[1][1] - the_tensor[0][1] * the_tensor[1][0])/det;
  return true;
}
