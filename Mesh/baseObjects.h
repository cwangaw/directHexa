#ifndef __baseobjects_h_included__
#define __baseobjects_h_included__

#include <cmath>
#include <iostream>
#include "debug.h"

////////////////////////////////////////////////////////////////////////////////
// Base Objects
// Class Point and Tensor1
////////////////////////////////////////////////////////////////////////////////

static const double baseobjects_eps = 1e-8;

class Point {
 protected:
  double the_point[3];
  void set_point(double p0, double p1, double p2) { the_point[0] = p0; the_point[1] = p1; the_point[2] = p2; };

 public:
  Point() { set_point(0,0,0); };
  Point(double p0, double p1, double p2) { set_point(p0,p1,p2); };
  Point(const Point& p) { set_point(p.the_point[0],p.the_point[1],p.the_point[2]); };

  Point& operator=(const Point& p1) {
    the_point[0]=p1.the_point[0]; the_point[1]=p1.the_point[1]; the_point[2]=p1.the_point[2]; return *this; }

  void set() { set_point(0,0,0); };
  void set(double p0=0, double p1=0, double p2=0) { set_point(p0,p1,p2); };
  void set(const Point& p) { set_point(p.the_point[0],p.the_point[1],p.the_point[2]); };

  double& operator() (int i)       { return the_point[i % 3]; }
  double  operator() (int i) const { return the_point[i % 3]; }
  double& operator[] (int i)       { return the_point[i % 3]; }
  double  operator[] (int i) const { return the_point[i % 3]; }
  double& val(int i)               { return the_point[i % 3]; }
  double  val(int i)         const { return the_point[i % 3]; }

  bool operator==(const Point& p1) const {
  return the_point[0]==p1.the_point[0] && the_point[1]==p1.the_point[1] && the_point[2]==p1.the_point[2]; }
  bool operator!=(const Point& p1) const {
  return the_point[0]!=p1.the_point[0] || the_point[1]!=p1.the_point[1] || the_point[2]!=p1.the_point[2]; }

 
  Point& operator+=(const Point& p1) {
    the_point[0] += p1.the_point[0]; the_point[1] += p1.the_point[1]; the_point[2] += p1.the_point[2]; return *this; }
  Point& operator-=(const Point& p1) {
    the_point[0] -= p1.the_point[0]; the_point[1] -= p1.the_point[1]; the_point[2] -= p1.the_point[2]; return *this; }
  Point& operator*=(double scalar) {
    the_point[0] *= scalar; the_point[1] *= scalar; the_point[2] *= scalar; return *this; }
  Point& operator/=(double scalar) {
    the_point[0] /= scalar; the_point[1] /= scalar; the_point[2] /= scalar; return *this; }

  friend class Tensor1;
};

inline bool eq_eps(const Point& p0, const Point& p1) {
  return std::abs(p0[0]-p1[0]) < baseobjects_eps && std::abs(p0[1]-p1[1]) < baseobjects_eps && std::abs(p0[2]-p1[2]) < baseobjects_eps; }

inline bool ne_eps(const Point& p0, const Point& p1) {
  return std::abs(p0[0]-p1[0]) > baseobjects_eps || std::abs(p0[1]-p1[1]) > baseobjects_eps || std::abs(p0[2]-p1[2]) > baseobjects_eps; }

Point operator+(const Point& p0, const Point& p1);
Point operator-(const Point& p0, const Point& p1);
Point operator*(double scalar, const Point& p0);
Point operator*(const Point& p0, double scalar);
Point operator/(const Point& p0, double scalar);

inline std::istream& operator>> (std::istream& is, Point& p) { is >> p[0] >> p[1] >> p[2]; return is; };
inline std::ostream& operator<< (std::ostream& os, const Point& p) { 
  os << p[0] << "," << p[1] << "," << p[2]; return os; 
};

class Tensor1 : public Point {
 public:
  Tensor1() { set_point(0,0,0); };
  Tensor1(const Point& p) { set_point(p.the_point[0],p.the_point[1],p.the_point[2]); };
  Tensor1(const Tensor1& p) { set_point(p.the_point[0],p.the_point[1],p.the_point[2]); };
  Tensor1(double p0, double p1, double p2) { set_point(p0,p1,p2); };

  double norm() const { return std::sqrt(the_point[0]*the_point[0] + the_point[1]*the_point[1] + the_point[2]*the_point[2]); }

  double operator*(const Tensor1& t1) { // scalar dot product
    return the_point[0]*t1.the_point[0] + the_point[1]*t1.the_point[1] + the_point[2]*t1.the_point[2]; };
};

inline bool eq_eps(const Tensor1& t0, const Tensor1& t1) {
  return std::abs(t0[0]-t1[0]) < baseobjects_eps && std::abs(t0[1]-t1[1]) < baseobjects_eps && std::abs(t0[2]-t1[2]) < baseobjects_eps; }

inline bool ne_eps(const Tensor1& t0, const Tensor1& t1) {
  return std::abs(t0[0]-t1[0]) > baseobjects_eps || std::abs(t0[1]-t1[1]) > baseobjects_eps || std::abs(t0[2]-t1[2]) > baseobjects_eps; }

inline Tensor1 cross(Tensor1 t0, Tensor1 t1) {
  return Tensor1(t0[1]*t1[2]-t0[2]*t1[1], t0[2]*t1[0]-t0[0]*t1[2], t0[0]*t1[1]-t0[1]*t1[0]);
};

class Tensor2 {
 private:
  double the_tensor[3][3];
  void set_tensor2(double t00, double t01, double t02, double t10, double t11, double t12, double t20, double t21, double t22) { 
    // By row
    the_tensor[0][0]=t00; the_tensor[0][1]=t01; the_tensor[0][2]=t02; 
    the_tensor[1][0]=t10; the_tensor[1][1]=t11; the_tensor[1][2]=t12;
    the_tensor[2][0]=t20; the_tensor[2][1]=t21; the_tensor[2][2]=t22;}

 public:
  Tensor2() { set_tensor2(1,0,0,0,1,0,0,0,1); };
  Tensor2(const Tensor2& t) { set_tensor2(t.the_tensor[0][0],t.the_tensor[0][1],t.the_tensor[0][2],
					  t.the_tensor[1][0],t.the_tensor[1][1],t.the_tensor[1][2],
            t.the_tensor[2][0],t.the_tensor[2][1],t.the_tensor[2][2]); }
  Tensor2(double t00, double t01, double t02, double t10, double t11, double t12, double t20, double t21, double t22) {
    set_tensor2(t00, t01, t02, t10, t11, t12, t20, t21, t22); };

  Tensor2& operator=(const Tensor2& t) {
    the_tensor[0][0]=t.the_tensor[0][0]; the_tensor[0][1]=t.the_tensor[0][1]; the_tensor[0][2]=t.the_tensor[0][2];
    the_tensor[1][0]=t.the_tensor[1][0]; the_tensor[1][1]=t.the_tensor[1][1]; the_tensor[1][2]=t.the_tensor[1][2];
    the_tensor[2][0]=t.the_tensor[2][0]; the_tensor[2][1]=t.the_tensor[2][1]; the_tensor[2][2]=t.the_tensor[2][2];
    return *this; }

  void set() { set_tensor2(1,0,0,0,1,0,0,0,1); };
  void set(double t00, double t01, double t02, double t10, double t11, double t12, double t20, double t21, double t22) {
    set_tensor2(t00, t01, t02, t10, t11, t12, t20, t21, t22); };
  void set(const Tensor2& t) { set_tensor2(t.the_tensor[0][0],t.the_tensor[0][1],t.the_tensor[0][2],
					   t.the_tensor[1][0],t.the_tensor[1][1],t.the_tensor[1][2],
             t.the_tensor[2][0],t.the_tensor[2][1],t.the_tensor[2][2]); }

  double& operator() (int i, int j)       { return the_tensor[i % 3][j % 3]; }
  double  operator() (int i, int j) const { return the_tensor[i % 3][j % 3]; }
  double& val(int i, int j)               { return the_tensor[i % 3][j % 3]; }
  double  val(int i, int j)         const { return the_tensor[i % 3][j % 3]; }
  Tensor1 col(int j)                      { return Tensor1(the_tensor[0][j % 3], the_tensor[1][j % 3], the_tensor[2][j % 3]); };

  bool operator==(const Tensor2& t) const { return
      the_tensor[0][0]==t.the_tensor[0][0] && the_tensor[0][1]==t.the_tensor[0][1] && the_tensor[0][2]==t.the_tensor[0][2] &&
      the_tensor[1][0]==t.the_tensor[1][0] && the_tensor[1][1]==t.the_tensor[1][1] && the_tensor[1][2]==t.the_tensor[1][2] &&
      the_tensor[2][0]==t.the_tensor[2][0] && the_tensor[2][1]==t.the_tensor[2][1] && the_tensor[1][2]==t.the_tensor[2][2] ; };
  bool operator!=(const Tensor2& t) const { return
      the_tensor[0][0]!=t.the_tensor[0][0] || the_tensor[0][1]!=t.the_tensor[0][1] || the_tensor[0][2]!=t.the_tensor[0][2] ||
      the_tensor[1][0]!=t.the_tensor[1][0] || the_tensor[1][1]!=t.the_tensor[1][1] || the_tensor[1][2]!=t.the_tensor[1][2] ||
      the_tensor[2][0]!=t.the_tensor[2][0] || the_tensor[2][1]!=t.the_tensor[2][1] || the_tensor[2][2]!=t.the_tensor[2][2]; };

  Tensor2& operator+=(const Tensor2& t) {
    the_tensor[0][0] += t.the_tensor[0][0]; the_tensor[0][1] += t.the_tensor[0][1]; the_tensor[0][2] += t.the_tensor[0][2];
    the_tensor[1][0] += t.the_tensor[1][0]; the_tensor[1][1] += t.the_tensor[1][1]; the_tensor[1][2] += t.the_tensor[1][2];
    the_tensor[2][0] += t.the_tensor[2][0]; the_tensor[2][1] += t.the_tensor[2][1]; the_tensor[2][2] += t.the_tensor[2][2];
    return *this; }
  Tensor2& operator-=(const Tensor2& t) {
    the_tensor[0][0] -= t.the_tensor[0][0]; the_tensor[0][1] -= t.the_tensor[0][1];
    the_tensor[1][0] -= t.the_tensor[1][0]; the_tensor[1][1] -= t.the_tensor[1][1];
    return *this; }
  Tensor2& operator*=(double scalar) {
    the_tensor[0][0] *= scalar; the_tensor[0][1] *= scalar; the_tensor[0][2] *= scalar;
    the_tensor[1][0] *= scalar; the_tensor[1][1] *= scalar; the_tensor[1][2] *= scalar;
    the_tensor[2][0] *= scalar; the_tensor[2][1] *= scalar; the_tensor[2][2] *= scalar;
    return *this; }
  Tensor2& operator/=(double scalar) {
    the_tensor[0][0] /= scalar; the_tensor[0][1] /= scalar; the_tensor[0][2] /= scalar;
    the_tensor[1][0] /= scalar; the_tensor[1][1] /= scalar; the_tensor[1][2] /= scalar;
    the_tensor[2][0] /= scalar; the_tensor[2][1] /= scalar; the_tensor[2][2] /= scalar;
    return *this; }

  void mult(double scalar, Tensor2& result) const;
  void mult(const Tensor1& v, Tensor1& result) const;
  void mult(const Point& v, Point& result) const;
  void mult(const Tensor2& t, Tensor2& result) const;
  Tensor2& operator*=(const Tensor2& t);

  double determinant() const;
  bool inverse(Tensor2& invT) const;

  Tensor1 operator*(const Tensor1& t1) { // Ab
    return Tensor1(the_tensor[0][0]*t1.val(0) + the_tensor[0][1]*t1.val(1) + the_tensor[0][2]*t1.val(2),
                   the_tensor[1][0]*t1.val(0) + the_tensor[1][1]*t1.val(1) + the_tensor[1][2]*t1.val(2),
                   the_tensor[2][0]*t1.val(0) + the_tensor[2][1]*t1.val(1) + the_tensor[2][2]*t1.val(2)); };

  friend class Tensor1;
};

inline std::istream& operator>> (std::istream& is, Tensor2& t) {
  is >> t(0,0) >> t(0,1) >> t(0,2) >> t(1,0) >> t(1,1) >> t(1,2) >> t(2,0) >> t(2,1) >> t(2,2); return is; }; // enter by rows
inline std::ostream& operator<< (std::ostream& os, const Tensor2& t) {
  os << t(0,0) << "," << t(0,1) << "," << t(0,2) << ";" << t(1,0) << "," << t(1,1) << "," << t(1,2) << ";" << t(2,0) << "," << t(2,1) << "," << t(2,2); return os; };

inline bool eq_eps(const Tensor2& t0, const Tensor2& t1) {
  return std::abs(t0(0,0)-t1(0,0)) < baseobjects_eps && std::abs(t0(0,1)-t1(0,1)) < baseobjects_eps && std::abs(t0(0,2)-t1(0,2)) < baseobjects_eps
      && std::abs(t0(1,0)-t1(1,0)) < baseobjects_eps && std::abs(t0(1,1)-t1(1,1)) < baseobjects_eps && std::abs(t0(1,2)-t1(1,2)) < baseobjects_eps
      && std::abs(t0(2,0)-t1(2,0)) < baseobjects_eps && std::abs(t0(2,1)-t1(2,1)) < baseobjects_eps && std::abs(t0(2,2)-t1(2,2)) < baseobjects_eps; }

inline bool ne_eps(const Tensor2& t0, const Tensor2& t1) {
  return std::abs(t0(0,0)-t1(0,0)) > baseobjects_eps || std::abs(t0(0,1)-t1(0,1)) > baseobjects_eps || std::abs(t0(0,2)-t1(0,2)) > baseobjects_eps
      || std::abs(t0(1,0)-t1(1,0)) > baseobjects_eps || std::abs(t0(1,1)-t1(1,1)) > baseobjects_eps || std::abs(t0(1,2)-t1(1,2)) > baseobjects_eps
      || std::abs(t0(2,0)-t1(2,0)) > baseobjects_eps || std::abs(t0(2,1)-t1(2,1)) > baseobjects_eps || std::abs(t0(2,2)-t1(2,2)) > baseobjects_eps; }


#endif