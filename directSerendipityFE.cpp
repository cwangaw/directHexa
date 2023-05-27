#include <cmath>
#include <iostream>
#include <vector>

#include <complex.h>
#include "lapacke.h"
#include <stdio.h>
#include <assert.h>
#include <cblas.h>

using namespace std;

#include "Utilities/debug.h"
#include "Mesh/hexaMesh.h"
using namespace hexamesh;
#include "directSerendipity.h"
using namespace directserendipity;

////////////////////////////////////////////////////////////////////////////////
// Class DirectSerendipityFE Usage Functions

// Matrix inversion
lapack_int mat_inv(double *A, int n)
{
  // inplace inverse n x n matrix A.
  // matrix A is arranged as: first row, second row ... 
  // returns:
  //   ret = 0 on success
  //   ret < 0 illegal argument value
  //   ret > 0 singular matrix
  int ipiv[n+1];
  lapack_int ret;

  ret =  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n,n,A,n,ipiv);

  if (ret !=0) return ret;
  ret = LAPACKE_dgetri(LAPACK_ROW_MAJOR,n,A,n,ipiv);
  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// Class DirectSerendipityFE

void DirectSerendipityFE::set_directserendipityfe(DirectSerendipity* dsSpace, Element* element) {
  my_ds_space = dsSpace;
  my_element = element;

  polynomial_degree = my_ds_space->degPolyn();
  deg_face_poly = polynomial_degree-4;
  num_face_dofs = (deg_face_poly<0)? 0:3*(polynomial_degree-2)*(polynomial_degree-3);
  deg_cell_poly = polynomial_degree-6;
  num_cell_dofs = (deg_cell_poly<0)? 0:(polynomial_degree-3)*(polynomial_degree-4)*(polynomial_degree-5)/6;
  num_dofs = 8 + 12*(polynomial_degree-1) + num_face_dofs + num_cell_dofs; // Redundant with above

  // Set up higher order element (if nescessary)
  if(polynomial_degree<3) {
    if(one_element_mesh) delete one_element_mesh;
    one_element_mesh = new hexamesh::HexaMesh(my_element);
    if(high_order_ds_space) delete high_order_ds_space;
    high_order_ds_space = new DirectSerendipity(3,my_ds_space->supplementType(),one_element_mesh);
  }

  // Set up face_basis_coefficients if num_face_dofs>0
  if (num_face_dofs>0) {
    if(face_basis_coefficients) delete face_basis_coefficients;
    face_basis_coefficients = new double[num_face_dofs*num_face_dofs/6];

    for (int nFace=0; nFace<6; nFace++) {
      int face_global_index = my_element->faceGlobal(nFace);

      std::vector<double> A_vec(num_face_dofs*num_face_dofs/36,0);
      double* A = A_vec.data();
      for (int iRow=0; iRow<num_face_dofs/6; iRow++) {
        Point* pt = my_ds_space->faceDoFPtr(face_global_index,iRow);
        for (int jCol=0; jCol<num_face_dofs/6; jCol++) {
          double value; Tensor1 gradvalue;
          faceVarphi(nFace, jCol, *pt, value, gradvalue);
          A[iRow*num_face_dofs/6+jCol] = value;
        }
      }
      (void) mat_inv(A,num_face_dofs/6);
      for (int iFunc=0; iFunc<num_face_dofs/6; iFunc++) { 
        for (int jCoeff=0; jCoeff<num_face_dofs/6; jCoeff++) {
          face_basis_coefficients[nFace*num_face_dofs*num_face_dofs/36+iFunc*num_face_dofs/6+jCoeff] = A[jCoeff*num_face_dofs/6+iFunc];
        }
      }
    }
  }

  // Set up edge_basis_coefficients if polynomial_degree>=3
  if (polynomial_degree>=3) {
    if (edge_basis_coefficients) delete edge_basis_coefficients;
    edge_basis_coefficients = new double[12*(polynomial_degree-1)*(polynomial_degree-1)];
    if (edgecheby_eval_mat_inv) delete edgecheby_eval_mat_inv;
    edgecheby_eval_mat_inv = new std::vector<double>[12];    
    if (edge_nodes) delete[] edge_nodes;
    edge_nodes = new Point[12*(polynomial_degree-1)];

    for (int nEdge=0; nEdge<12; nEdge++) {
      std::vector<double> A_vec((polynomial_degree-1)*(polynomial_degree-1),0);
      double* A = A_vec.data();

      std::vector<double> B_vec((polynomial_degree-1)*(polynomial_degree-1),0);
      double* B = B_vec.data();

      std::vector<double> coeff_vec((polynomial_degree-1)*(polynomial_degree-1),0);
      double* coeff = coeff_vec.data();

      Vertex* v0 = my_element->edgeVertexPtr(0,nEdge);
      Vertex* v1 = my_element->edgeVertexPtr(1,nEdge);
      Tensor1 tangent(*v1-*v0);
      double seg = tangent.norm()/polynomial_degree;
      tangent /= tangent.norm();
      for (int iRow=0; iRow<polynomial_degree-1; iRow++) {
        Point pt = *v0 + (iRow+1)*seg*tangent;
        edge_nodes[nEdge*(polynomial_degree-1)+iRow] = pt;
        
        for (int jCol=0; jCol<polynomial_degree-1; jCol++) {
          double value; Tensor1 gradvalue;
          edgeVarphi(nEdge, jCol, pt, value, gradvalue);
          A[iRow*(polynomial_degree-1)+jCol] = value;
          B[iRow*(polynomial_degree-1)+jCol] = edgeCheby(nEdge, iRow, jCol);
        }        
      }

      // coeff = A^{-1}*B
      (void) mat_inv(A,polynomial_degree-1);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, polynomial_degree-1, polynomial_degree-1,
              polynomial_degree-1, 1, A, polynomial_degree-1, B, polynomial_degree-1,
              0, coeff, polynomial_degree-1);
      for (int iFunc=0; iFunc<polynomial_degree-1; iFunc++) { 
        for (int jCoeff=0; jCoeff<polynomial_degree-1; jCoeff++) {
          edge_basis_coefficients[nEdge*(polynomial_degree-1)*(polynomial_degree-1)+iFunc*(polynomial_degree-1)+jCoeff] = coeff[jCoeff*(polynomial_degree-1)+iFunc];
        }
      }

      (void) mat_inv(B,polynomial_degree-1);
      edgecheby_eval_mat_inv[nEdge]=B_vec;
    }
    
  }
}


DirectSerendipityFE::~DirectSerendipityFE() {
  if (high_order_ds_space) delete high_order_ds_space;
  if (one_element_mesh) delete one_element_mesh;
  if (face_basis_coefficients) delete[] face_basis_coefficients;
  if (edge_basis_coefficients) delete[] edge_basis_coefficients;
  if (edge_nodes) delete[] edge_nodes;
  if (edgecheby_eval_mat_inv) delete[] edgecheby_eval_mat_inv;
  if (value_n) delete[] value_n;
  if (gradvalue_n) delete[] gradvalue_n;
}

void DirectSerendipityFE::getAB(int n, int m, int iEdge, double& A, double& B) {
  Vertex* v0 = my_element->edgeVertexPtr(0,iEdge);
  Vertex* v1 = my_element->edgeVertexPtr(1,iEdge);
  B = (my_element->lambda(n,*v1)-my_element->lambda(n,*v0)) / (my_element->lambda(m,*v1)-my_element->lambda(m,*v0));
  A = my_element->lambda(n,*v0) - B*my_element->lambda(m,*v0);
  return;
}

void DirectSerendipityFE::faceVarphi(int n, int s, const Point& pt, double& value, Tensor1& gradvalue) {
  int s0 = ceil(0.5*(2*deg_face_poly+3-sqrt((2*deg_face_poly+3)*(2*deg_face_poly+3)-8*(s+1))))-1;
  int s1 = s - 0.5*(2*deg_face_poly+3-s0)*s0;

  n = n%6;
  int n_opposite=0, dir=0, n0=0, n1=0;
  if (n==0 || n==1) { dir=0; n_opposite=1-n; n0=2; n1=4; }
  if (n==2 || n==3) { dir=1; n_opposite=5-n; n0=0; n1=4; }
  if (n==4 || n==5) { dir=2; n_opposite=9-n; n0=0; n1=2; }

  
  double common = pow(my_element->lambda(n0, pt),s0) * pow(my_element->lambda(n1, pt),s1);
  gradvalue = (s0==0)? Tensor1(0,0,0) : s0 * pow(my_element->lambda(n0, pt),s0-1) * pow(my_element->lambda(n1, pt),s1) * my_element->dLambda(n0);
  gradvalue += (s1==0)? Tensor1(0,0,0) : s1 * pow(my_element->lambda(n0, pt),s0) * pow(my_element->lambda(n1, pt),s1-1) * my_element->dLambda(n1);
  for (int nFace=0; nFace<6; nFace++) {
    if (nFace==n || nFace==n_opposite) continue;
    common *= my_element->lambda(nFace, pt);
    gradvalue *= my_element->lambda(nFace, pt);
  }

  for (int nFace=0; nFace<6; nFace++) {
    if (nFace==n || nFace==n_opposite) continue;
    Tensor1 dL = pow(my_element->lambda(n0, pt),s0) * pow(my_element->lambda(n1, pt),s1) * my_element->dLambda(nFace);
    for (int mFace=0; mFace<6; mFace++) {
      if (mFace==n || mFace==n_opposite || mFace==nFace) continue;
      dL *= my_element->lambda(mFace, pt);
    }
    gradvalue += dL;
  }
  
  if (s0+s1<deg_face_poly) {
    value = common*my_element->lambda(n_opposite, pt);
    gradvalue *= my_element->lambda(n_opposite, pt);
    gradvalue += common*my_element->dLambda(n_opposite);
  } else {
    double specfunc_value;
    Tensor1 specfunc_gradvalue;
    specFunc(dir,pt,specfunc_value,specfunc_gradvalue);
    value = (n%2==0)? 0.5*common*(1-specfunc_value) : 0.5*common*(1+specfunc_value);
    gradvalue *= (n%2==0)? 0.5*(1-specfunc_value) : 0.5*(1+specfunc_value);
    gradvalue += (n%2==0)? -0.5*common*specfunc_gradvalue : 0.5*common*specfunc_gradvalue;
  }
  return;
}

void DirectSerendipityFE::edgeVarphi(int n, int s, const Point& pt, double& value, Tensor1& gradvalue) {
  assert(s<=polynomial_degree-2);
  int f0, f1;
  int f_neg = 4-(n/4)*2;
  my_element -> edgeFace(n,f0,f1);

  if (s<=polynomial_degree-4) {
    value = pow(my_element->lambda(f_neg,pt),s);
    gradvalue = (s==0)? Tensor1(0,0,0) : s*pow(my_element->lambda(f_neg,pt),s-1)*my_element->dLambda(f_neg);
    for (int nFace=0; nFace<6; nFace++) {
      if (nFace==f0 || nFace==f1) continue;
      value *= my_element->lambda(nFace,pt);
      gradvalue *= my_element->lambda(nFace,pt);
    }

    Tensor1 dL;
    for (int nFace=0; nFace<6; nFace++) {
      if (nFace==f0 || nFace==f1) continue;
      dL = pow(my_element->lambda(f_neg,pt),s)*my_element->dLambda(nFace);
      for (int mFace=0; mFace<6; mFace++) {
        if (mFace==f0 || mFace==f1 || mFace==nFace) continue;
        dL *= my_element->lambda(mFace,pt);
      }
      gradvalue += dL;
    }
  } else if (s==polynomial_degree-3) {
    double value_psi;
    Tensor1 gradvalue_psi;
    edgePsi(f_neg/2,pt,value_psi,gradvalue_psi);
    if (n%4==3) {
      // +,+
      value = pow(my_element->lambda(f_neg,pt),s+1)*my_element->lambda(f_neg+1,pt)*value_psi;
      gradvalue =((s+1)*pow(my_element->lambda(f_neg,pt),s)*my_element->lambda(f_neg+1,pt)*my_element->dLambda(f_neg)
                + pow(my_element->lambda(f_neg,pt),s+1)*my_element->dLambda(f_neg+1))*value_psi
                + pow(my_element->lambda(f_neg,pt),s+1)*my_element->lambda(f_neg+1,pt)*gradvalue_psi;
    } else if (n%4>0) {
      // +,- & -,+
      double A_this, A_pp, B_this, B_pp;
      if (f0%2==1 && f1%2==0) {
        getAB(f0-1,f_neg,n,A_this,B_this);
        getAB(f0-1,f_neg,4*(n/4)+3,A_pp, B_pp);
      } else {
        getAB(f1-1,f_neg,n,A_this,B_this);
        getAB(f1-1,f_neg,4*(n/4)+3,A_pp, B_pp);
      }

      value = pow(my_element->lambda(f_neg,pt),s+1)*my_element->lambda(f_neg+1,pt);
      gradvalue = (s+1)*pow(my_element->lambda(f_neg,pt),s)*my_element->lambda(f_neg+1,pt)*my_element->dLambda(f_neg)
                + pow(my_element->lambda(f_neg,pt),s+1)*my_element->dLambda(f_neg+1);
      value *= (f0%2==1 && f1%2==0)? my_element->lambda(f0-1,pt) : my_element->lambda(f1-1,pt);
      gradvalue *= (f0%2==1 && f1%2==0)? my_element->lambda(f0-1,pt) : my_element->lambda(f1-1,pt);
      gradvalue +=  pow(my_element->lambda(f_neg,pt),s+1)*my_element->lambda(f_neg+1,pt)*((f0%2==1 && f1%2==0)? my_element->dLambda(f0-1) : my_element->dLambda(f1-1));
      
      double value_varphi;
      Tensor1 gradvalue_varphi;
      edgeVarphi(4*(n/4)+3,s,pt,value_varphi,gradvalue_varphi);
      value -= A_pp*value_varphi;
      gradvalue -= A_pp*gradvalue_varphi;
      edgeVarphi(4*(n/4)+3,s+1,pt,value_varphi,gradvalue_varphi);
      value -= B_pp*value_varphi;
      gradvalue -= B_pp*gradvalue_varphi;
      edgeVarphi(n,s+1,pt,value_varphi,gradvalue_varphi);
      value -= B_this*value_varphi;
      gradvalue -= B_this*gradvalue_varphi;

      value /= A_this;
      gradvalue /= A_this;
    } else {
      // -,-
      value = pow(my_element->lambda(f_neg,pt),s+1)*my_element->lambda(f_neg+1,pt);
      gradvalue = (s+1)*pow(my_element->lambda(f_neg,pt),s)*my_element->lambda(f_neg+1,pt)*my_element->dLambda(f_neg)
                + pow(my_element->lambda(f_neg,pt),s+1)*my_element->dLambda(f_neg+1);
      double value_varphi;
      Tensor1 gradvalue_varphi;
      for (int i=1; i<=3; i++) {
        edgeVarphi(n+i,s,pt,value_varphi,gradvalue_varphi);
        value -= value_varphi;
        gradvalue -= gradvalue_varphi;
      }
    }
  } else {
    // s=polynomial_degree-2
    double value_spec0, value_spec1;
    Tensor1 gradvalue_spec0, gradvalue_spec1;
    specFunc(f0/2,pt,value_spec0, gradvalue_spec0);
    specFunc(f1/2,pt,value_spec1, gradvalue_spec1);
    value =  0.25*pow(my_element->lambda(f_neg,pt),s+1)*my_element->lambda(f_neg+1,pt) 
                * (1+(2*(f0%2)-1)*value_spec0)*(1+(2*(f1%2)-1)*value_spec1);
    gradvalue = ((s+1)*pow(my_element->lambda(f_neg,pt),s)*my_element->lambda(f_neg+1,pt)*my_element->dLambda(f_neg)
                + pow(my_element->lambda(f_neg,pt),s+1)*my_element->dLambda(f_neg+1))
                * (1+(2*(f0%2)-1)*value_spec0)*(1+(2*(f1%2)-1)*value_spec1);
    gradvalue += pow(my_element->lambda(f_neg,pt),s+1)*my_element->lambda(f_neg+1,pt)
                * (((2*(f0%2)-1)*gradvalue_spec0)*(1+(2*(f1%2)-1)*value_spec1) + (1+(2*(f0%2)-1)*value_spec0)*((2*(f1%2)-1)*gradvalue_spec1));
    gradvalue *= 0.25;
  }
}

double DirectSerendipityFE::edgeCheby(int iEdge, int nPt, int s) {
  return my_ds_space->edgeCheby(my_element->edgeGlobal(iEdge), nPt, s);
}

void DirectSerendipityFE::specFunc(int dir, const Point& pt, double& value, Tensor1& gradvalue) {
  assert(dir>=0 && dir <=2);
  switch(my_ds_space->supplementType()) {
  case 0: case 1: {
    // piecewise polynomial with T_D or T_M
    std::vector<std::vector<int>>* subtetra = my_element->subtetrahedra(my_ds_space->supplementType());
    int tetra = my_element->inSubtetrahedron(*subtetra, pt);
    // pos_vert store the vertex indices with value 1
    std::vector<int> pos_vert; pos_vert.clear();
    // neg_vert store the vertex indices with value -1
    std::vector<int> neg_vert; neg_vert.clear();
    for (int i=0; i<4; i++) {
      int sgn[3];
      my_element->vertexPos((*subtetra)[tetra][i], sgn[0], sgn[1], sgn[2]);
      if (sgn[dir] > 0) { 
        pos_vert.push_back((*subtetra)[tetra][i]); 
      } else {
        neg_vert.push_back((*subtetra)[tetra][i]);
      }
    }
    switch (pos_vert.size()) {
      case 1: {
        Tensor1 normal = cross(Tensor1(*my_element->vertexPtr(neg_vert[1]) - *my_element->vertexPtr(neg_vert[0])),
                               Tensor1(*my_element->vertexPtr(neg_vert[2]) - *my_element->vertexPtr(neg_vert[0])));
        normal /= normal.norm();
        double scaling = Tensor1(*my_element->vertexPtr(pos_vert[0])-*my_element->vertexPtr(neg_vert[0]))*normal;
        if ( scaling<0 ) { normal *= -1; scaling *= -1; }
        value = 2 * (Tensor1(pt-*my_element->vertexPtr(neg_vert[0]))*normal) / scaling - 1;
        gradvalue = 2 * normal / scaling;
        break;
      }
      case 2: {
        Tensor1 normal = cross(Tensor1(*my_element->vertexPtr(pos_vert[1]) - *my_element->vertexPtr(pos_vert[0])),
                               Tensor1(*my_element->vertexPtr(neg_vert[1]) - *my_element->vertexPtr(neg_vert[0])));
        normal /= normal.norm();
        double scaling = Tensor1(*my_element->vertexPtr(pos_vert[0])-*my_element->vertexPtr(neg_vert[0]))*normal;
        if ( scaling<0 ) { normal *= -1; scaling *= -1; }
        value = 2 * (Tensor1(pt-*my_element->vertexPtr(neg_vert[0]))*normal) / scaling - 1;
        gradvalue = 2 * normal / scaling;
        break;
      }
      case 3: {
        Tensor1 normal = cross(Tensor1(*my_element->vertexPtr(pos_vert[1]) - *my_element->vertexPtr(pos_vert[0])),
                               Tensor1(*my_element->vertexPtr(pos_vert[2]) - *my_element->vertexPtr(pos_vert[0])));
        normal /= normal.norm();
        double scaling = Tensor1(*my_element->vertexPtr(neg_vert[0])-*my_element->vertexPtr(pos_vert[0]))*normal;
        if ( scaling<0 ) { normal *= -1; scaling *= -1; }
        value = -2 * (Tensor1(pt-*my_element->vertexPtr(pos_vert[0]))*normal) / scaling + 1;
        gradvalue = -2 * normal / scaling;
        break;
      }
      default: std::cerr << "Error: Problematic subtetrahedral division" << std::endl;
    }
    break;
  }
  case 2: {
    // pullback map
    Point p_orig = my_element->backwardMap(pt);
    Tensor2 df; Tensor2 df_inv; double jac;
    my_element->dForwardMap(p_orig, df, jac);
    df.inverse(df_inv);
    value = p_orig[dir];
    gradvalue.set(df_inv.val(dir,0), df_inv.val(dir,1), df_inv.val(dir,2));
    break;
  } 
  default: return;
  }	   
  return;
}

void DirectSerendipityFE::edgePsi(int dir, const Point& pt, double& value, Tensor1& gradvalue) {
  assert(dir>=0 && dir <=2);
  switch (my_ds_space->supplementType()) {
  case 0: case 1: {
    // piecewise polynomial with T_D and T_M
    std::vector<std::vector<int>>* subtetra = my_element->subtetrahedra(my_ds_space->supplementType());
    int tetra = my_element->inSubtetrahedron(*subtetra, pt);
    std::vector<double> eval; eval.clear();

    int iEdge = 4*(2-dir)+3;
    double R0, R1; Tensor1 gradR0, gradR1;
    double A0, B0, A1, B1;
    getAB((dir+2)%3*2,dir*2,iEdge,A0,B0);
    getAB((dir+1)%3*2,dir*2,iEdge,A1,B1);

    for (int iV=0; iV<4; iV++) {
      int sgn[3];
      my_element->vertexPos((*subtetra)[tetra][iV],sgn[0],sgn[1],sgn[2]);
      if (sgn[(dir+1)%3]==1 && sgn[(dir+2)%3]==1) {
        eval.push_back(1);
      } else {
        eval.push_back(0);
      }
    }

    int iSgn[3], jSgn[3];
    for (int iV=0; iV<4; iV++) {
      my_element->vertexPos((*subtetra)[tetra][iV],iSgn[0],iSgn[1],iSgn[2]);
      for(int jV=iV+1; jV<4; jV++) {
        my_element->vertexPos((*subtetra)[tetra][jV],jSgn[0],jSgn[1],jSgn[2]);
        if (iSgn[(dir+1)%3]==1 && jSgn[(dir+1)%3]==1) {
          Point pt=(*my_element->vertexPtr((*subtetra)[tetra][iV])+*my_element->vertexPtr((*subtetra)[tetra][jV]))/2;
          specFunc((dir+2)%3,pt,R0,gradR0);
          eval.push_back((my_element->lambda((dir+2)%3*2,pt)-0.5*B0*my_element->lambda(dir*2,pt)*(1+R0))/A0);
        } else if (iSgn[(dir+2)%3]==1 && jSgn[(dir+2)%3]==1) {
          Point pt=(*my_element->vertexPtr((*subtetra)[tetra][iV])+*my_element->vertexPtr((*subtetra)[tetra][jV]))/2;
          specFunc((dir+1)%3,pt,R1,gradR1);
          eval.push_back((my_element->lambda((dir+1)%3*2,pt)-0.5*B1*my_element->lambda(dir*2,pt)*(1+R1))/A1);
        } else {
          eval.push_back((eval[iV]+eval[jV])/2);
        }
      }
    }

    value=0; gradvalue.set(0,0);
    for (int iV=0; iV<4; iV++) {
      Tensor1 normal = cross(Tensor1(*my_element->vertexPtr((*subtetra)[tetra][(iV+2)%4]) - *my_element->vertexPtr((*subtetra)[tetra][(iV+1)%4])),
                              Tensor1(*my_element->vertexPtr((*subtetra)[tetra][(iV+3)%4]) - *my_element->vertexPtr((*subtetra)[tetra][(iV+1)%4])));
      normal /= normal.norm();
      double scaling = Tensor1(*my_element->vertexPtr((*subtetra)[tetra][iV])-*my_element->vertexPtr((*subtetra)[tetra][(iV+1)%4]))*normal;
      if ( scaling<0 ) { normal *= -1; scaling *= -1; }
      scaling *= 0.5*scaling;
      Point midpt = (*my_element->vertexPtr((*subtetra)[tetra][iV])+*my_element->vertexPtr((*subtetra)[tetra][(iV+1)%4]))/2;
      value += eval[iV] * (Tensor1(pt-*my_element->vertexPtr((*subtetra)[tetra][(iV+1)%4]))*normal) 
                        * (Tensor1(pt-midpt)*normal) / scaling;
      gradvalue += eval[iV] * ((Tensor1(pt-*my_element->vertexPtr((*subtetra)[tetra][(iV+1)%4]))*normal) * normal
                              + normal*(Tensor1(pt-midpt)*normal)) / scaling;
    }

    int ind=3;
    for (int iV=0; iV<4; iV++) {
      for (int jV=iV+1; jV<4; jV++) {
        ind++;
        int common_vertices[2]; int vind=0;
        for (int kV=0; kV<4; kV++) {
          if (kV!=iV && kV!=jV) {
            common_vertices[vind]=(*subtetra)[tetra][kV];
            vind++;
          }
        }
        Point midpt = (*my_element->vertexPtr((*subtetra)[tetra][iV])+*my_element->vertexPtr((*subtetra)[tetra][jV]))/2;
        Tensor1 normal0 = cross(Tensor1(*my_element->vertexPtr((*subtetra)[tetra][iV]) - *my_element->vertexPtr(common_vertices[0])),
                                Tensor1(*my_element->vertexPtr(common_vertices[1]) - *my_element->vertexPtr(common_vertices[0])));
        normal0 /= normal0.norm();
        double scaling0 = Tensor1(midpt-*my_element->vertexPtr((*subtetra)[tetra][iV]))*normal0;
        if ( scaling0<0 ) { normal0 *= -1; scaling0 *= -1; }
        Tensor1 normal1 = cross(Tensor1(*my_element->vertexPtr((*subtetra)[tetra][jV]) - *my_element->vertexPtr(common_vertices[0])),
                                Tensor1(*my_element->vertexPtr(common_vertices[1]) - *my_element->vertexPtr(common_vertices[0])));
        normal1 /= normal1.norm();
        double scaling1 = Tensor1(midpt-*my_element->vertexPtr((*subtetra)[tetra][jV]))*normal1;
        if ( scaling1<0 ) { normal1 *= -1; scaling1 *= -1; }

        value += eval[ind] * (Tensor1(pt-*my_element->vertexPtr((*subtetra)[tetra][iV]))*normal0)
                           * (Tensor1(pt-*my_element->vertexPtr((*subtetra)[tetra][jV]))*normal1) / (scaling0*scaling1);
        gradvalue += eval[ind] * ((Tensor1(pt-*my_element->vertexPtr((*subtetra)[tetra][iV]))*normal0)*normal1 
                                + normal0*(Tensor1(pt-*my_element->vertexPtr((*subtetra)[tetra][jV]))*normal1)) / (scaling0*scaling1);
      }
    }
    break;
  } 
  case 2: {
    // pullback map
    Point p_orig = my_element->backwardMap(pt);
    Tensor2 df; Tensor2 df_inv; double jac;
    my_element->dForwardMap(p_orig, df, jac);
    df.inverse(df_inv);

    int iEdge = 4*(2-dir)+3;

    Point eval_p0, eval_p1;
    double R0, R1; Tensor1 gradR0, gradR1;
    double A0, B0, A1, B1;
    Tensor2 chain0=df_inv, chain1=df_inv, dftemp; double jactemp;

    switch (dir) {
      case 0: {
        eval_p0 = my_element->forwardMap(Point(p_orig[0],1,p_orig[2])); //y
        my_element->dForwardMap(Point(p_orig[0],1,p_orig[2]),dftemp,jactemp);
        chain0 *= Tensor2(1,0,0,0,0,0,0,0,1); chain0 *= dftemp;
        eval_p1 = my_element->forwardMap(Point(p_orig[0],p_orig[1],1)); //z
        my_element->dForwardMap(Point(p_orig[0],p_orig[1],1),dftemp,jactemp);
        chain1 *= Tensor2(1,0,0,0,1,0,0,0,0); chain1 *= dftemp;
        break;
      }
      case 1: {
        eval_p0 = my_element->forwardMap(Point(p_orig[0],p_orig[1],1)); //z
        my_element->dForwardMap(Point(p_orig[0],p_orig[1],1),dftemp,jactemp);
        chain0 *= Tensor2(1,0,0,0,1,0,0,0,0); chain0 *= dftemp;
        eval_p1 = my_element->forwardMap(Point(1,p_orig[1],p_orig[2])); //x
        my_element->dForwardMap(Point(1,p_orig[1],p_orig[2]),dftemp,jactemp);
        chain1 *= Tensor2(0,0,0,0,1,0,0,0,1); chain1 *= dftemp;
        break;
      }
      case 2: {
        eval_p0 = my_element->forwardMap(Point(1,p_orig[1],p_orig[2])); //x
        my_element->dForwardMap(Point(1,p_orig[1],p_orig[2]),dftemp,jactemp);
        chain0 *= Tensor2(0,0,0,0,1,0,0,0,1); chain0 *= dftemp;
        eval_p1 = my_element->forwardMap(Point(p_orig[0],1,p_orig[2])); //y
        my_element->dForwardMap(Point(p_orig[0],1,p_orig[2]),dftemp,jactemp);
        chain1 *= Tensor2(1,0,0,0,0,0,0,0,1); chain1 *= dftemp;
        break;
      }
    }

    specFunc((dir+2)%3,pt,R0,gradR0);
    specFunc((dir+1)%3,pt,R1,gradR1);

    getAB((dir+2)%3*2,dir*2,iEdge,A0,B0);
    getAB((dir+1)%3*2,dir*2,iEdge,A1,B1);

    double psi0 = (my_element->lambda((dir+2)%3*2,eval_p0)-0.5*B0*my_element->lambda(dir*2,eval_p0)*(1+R0))/A0;
    double psi1 = (my_element->lambda((dir+1)%3*2,eval_p1)-0.5*B1*my_element->lambda(dir*2,eval_p1)*(1+R1))/A1;

    value = psi0 * psi1;
    gradvalue = (chain0 * (my_element->dLambda((dir+2)%3*2)-0.5*B0*(my_element->dLambda(dir*2)*(1+R0)+my_element->lambda(dir*2,eval_p0)*gradR0))/A0) * psi1
              + (chain1 * (my_element->dLambda((dir+1)%3*2)-0.5*B1*(my_element->dLambda(dir*2)*(1+R1)+my_element->lambda(dir*2,eval_p1)*gradR1))/A1) * psi0;

    break;
  }
  default: return;
  }
  return;
}

void DirectSerendipityFE::initBasis(const Point* pt, int num_pts) {
  if(num_pts <= 0) return;
  num_eval_pts = num_pts;


  // Allocate space for the resulting values
  int sizeOfArray = num_pts * num_dofs;

  if(value_n) delete[] value_n;
  value_n = new double[sizeOfArray];
  if (gradvalue_n) delete[] gradvalue_n;
  gradvalue_n = new Tensor1[sizeOfArray];  

  if (polynomial_degree < 3) {
    // Define DS_1 and DS_2 as subspaces of DS_3 
    high_order_ds_space->finiteElementPtr(0)->initBasis(pt, num_pts);

    // Update phi_{v,i}
    for (int i = 0; i < 8; i++) {
      // We directly take higher order vertex basis functions that are linear on each edge
      for (int pt_index = 0; pt_index < num_pts; pt_index++) {
        value_n[pt_index*nDoFs()+i] = high_order_ds_space->finiteElementPtr(0)->vertexBasis(i,pt_index);
        gradvalue_n[pt_index*nDoFs()+i] = high_order_ds_space->finiteElementPtr(0)->gradVertexBasis(i,pt_index);
      }
    }

    // Update phi_{e,nEdge,jNode}
    if (polynomial_degree > 1) {
      for (int nEdge = 0; nEdge < 12; nEdge++) {
          for (int pt_index = 0; pt_index < num_pts; pt_index++) {
            // We directly take higher order edge basis functions that are quadratic on each edge
            value_n[pt_index*nDoFs()+nVertexDoFs()+nEdge] = high_order_ds_space->finiteElementPtr(0)->edgeBasis(2*nEdge,pt_index);
            gradvalue_n[pt_index*nDoFs()+nVertexDoFs()+nEdge] = high_order_ds_space->finiteElementPtr(0)->gradEdgeBasis(2*nEdge,pt_index);
          }
      }  
    }
  } else {
    // Cell Basis Functions
    if (num_cell_dofs>0) {
      int iCellDoF = 0;
      for (int n=0; n<=deg_cell_poly; n++) {
        for (int m=0; m<=deg_cell_poly-n; m++) {
          for (int l=0; l<=deg_cell_poly-n-m; l++) {
            for (int pt_index = 0; pt_index < num_pts; pt_index++) {
              double value = pow(my_element->lambda(0,pt[pt_index]),n) * pow(my_element->lambda(2,pt[pt_index]),m) * pow(my_element->lambda(4,pt[pt_index]),l);
              for (int iFace=0; iFace<6; iFace++) {
                value *= my_element->lambda(iFace,pt[pt_index]);
              }
              value_n[pt_index*nDoFs()+nVertexDoFs()+nEdgeDoFs()+nFaceDoFs()+iCellDoF] = value;

              Tensor1 gradvalue = my_element->lambda(3,pt[pt_index])*my_element->lambda(5,pt[pt_index]) * my_element->dLambda(1)
                                + my_element->lambda(1,pt[pt_index])*my_element->lambda(5,pt[pt_index]) * my_element->dLambda(3)
                                + my_element->lambda(1,pt[pt_index])*my_element->lambda(3,pt[pt_index]) * my_element->dLambda(5);
              gradvalue *= pow(my_element->lambda(0,pt[pt_index]),n+1) * pow(my_element->lambda(2,pt[pt_index]),m+1) * pow(my_element->lambda(4,pt[pt_index]),l+1);
              gradvalue += pow(my_element->lambda(0,pt[pt_index]),n) * my_element->lambda(1,pt[pt_index])
                         * pow(my_element->lambda(2,pt[pt_index]),m) * my_element->lambda(3,pt[pt_index])
                         * pow(my_element->lambda(4,pt[pt_index]),l) * my_element->lambda(5,pt[pt_index])
                         * ((n+1)*my_element->lambda(2,pt[pt_index])*my_element->lambda(4,pt[pt_index])*my_element->dLambda(0)
                           +(m+1)*my_element->lambda(0,pt[pt_index])*my_element->lambda(4,pt[pt_index])*my_element->dLambda(2)
                           +(l+1)*my_element->lambda(0,pt[pt_index])*my_element->lambda(2,pt[pt_index])*my_element->dLambda(4));
              gradvalue_n[pt_index*nDoFs()+nVertexDoFs()+nEdgeDoFs()+nFaceDoFs()+iCellDoF] = gradvalue;
              iCellDoF++;
            }
          }
        }
      }
    }

    // Face Basis Functions
    if (num_face_dofs>0) {
      for (int nFace=0; nFace<6; nFace++) {
        for (int iFunc=0; iFunc<num_face_dofs/6; iFunc++) {
          for (int pt_index = 0; pt_index < num_pts; pt_index++) {
            double value=0; 
            Tensor1 gradvalue(0,0,0);
            
            double value_varphi; Tensor1 gradvalue_varphi;
            for (int jCoeff=0; jCoeff<num_face_dofs/6; jCoeff++) {
              faceVarphi(nFace, jCoeff, pt[pt_index], value_varphi, gradvalue_varphi);
              value += face_basis_coefficients[nFace*num_face_dofs*num_face_dofs/36+iFunc*num_face_dofs/6+jCoeff] * value_varphi;
              gradvalue += face_basis_coefficients[nFace*num_face_dofs*num_face_dofs/36+iFunc*num_face_dofs/6+jCoeff] * gradvalue_varphi;
            }
            value_n[pt_index*nDoFs()+nVertexDoFs()+nEdgeDoFs()+ nFace*num_face_dofs/6+iFunc] = value;
            gradvalue_n[pt_index*nDoFs()+nVertexDoFs()+nEdgeDoFs()+ nFace*num_face_dofs/6+iFunc] = gradvalue;
          }
        }
      }
    }
    // Edge Basis Functions
    for (int nEdge=0; nEdge<12; nEdge++) {
      // Find out the two faces containing this edge
      int f_index[2];
      my_element->edgeFace(nEdge, f_index[0], f_index[1]);

      for (int iFunc=0; iFunc<polynomial_degree-1; iFunc++) {
        // If num_face_dofs>0, first evaluate the edge basis functions at 
        // all the face nodes on the two faces containing this edge
        std::vector<double> eval_face_nodes;
        eval_face_nodes.clear();

        if (num_face_dofs>0) {
          for (int n=0; n<2; n++) {
            for (int iDoF=0; iDoF<num_face_dofs/6; iDoF++) {
              Point* node=my_ds_space->faceDoFPtr(my_element->faceGlobal(f_index[n]),iDoF);
              double value=0; 
              Tensor1 gradvalue(0,0,0);
            
              double value_varphi; Tensor1 gradvalue_varphi;
              for (int jCoeff=0; jCoeff<polynomial_degree-1; jCoeff++) {
                edgeVarphi(nEdge, jCoeff, *node, value_varphi, gradvalue_varphi);
                value += edge_basis_coefficients[nEdge*(polynomial_degree-1)*(polynomial_degree-1)+iFunc*(polynomial_degree-1)+jCoeff] * value_varphi;
              }
              eval_face_nodes.push_back(value);
            }          
          }
        }
        for (int pt_index = 0; pt_index < num_pts; pt_index++) {
          double value=0; 
          Tensor1 gradvalue(0,0,0);
          
          double value_varphi; Tensor1 gradvalue_varphi;
          for (int jCoeff=0; jCoeff<polynomial_degree-1; jCoeff++) {
            edgeVarphi(nEdge, jCoeff, pt[pt_index], value_varphi, gradvalue_varphi);
            value += edge_basis_coefficients[nEdge*(polynomial_degree-1)*(polynomial_degree-1)+iFunc*(polynomial_degree-1)+jCoeff] * value_varphi;
            gradvalue += edge_basis_coefficients[nEdge*(polynomial_degree-1)*(polynomial_degree-1)+iFunc*(polynomial_degree-1)+jCoeff] * gradvalue_varphi;
          }

          if (num_face_dofs>0) {
            for (int n=0; n<2; n++) {
              for (int iDoF=0; iDoF<num_face_dofs/6; iDoF++) {
                value -= eval_face_nodes[n*num_face_dofs/6+iDoF] * value_n[pt_index*nDoFs()+nVertexDoFs()+nEdgeDoFs()+f_index[n]*num_face_dofs/6+iDoF];
                gradvalue -= eval_face_nodes[n*num_face_dofs/6+iDoF] * gradvalue_n[pt_index*nDoFs()+nVertexDoFs()+nEdgeDoFs()+f_index[n]*num_face_dofs/6+iDoF];
              }
            }   
          }

          value_n[pt_index*nDoFs()+nVertexDoFs()+nEdge*(polynomial_degree-1)+iFunc] = value;
          gradvalue_n[pt_index*nDoFs()+nVertexDoFs()+nEdge*(polynomial_degree-1)+iFunc] = gradvalue;   
        }
      }
    }
    
    // Vertex Basis Functions
    for (int nVertex=0; nVertex<8; nVertex++) {
      Vertex* v = my_element->vertexPtr(nVertex);
      int sgnx, sgny, sgnz;
      my_element->vertexPos(nVertex,sgnx,sgny,sgnz);
      int v_index[3] = { sgnx, sgny, sgnz };
      // Find faces and edges that connect this vertex
      int f_connect[3] = { my_element->faceIndex(sgnx,0,0),my_element->faceIndex(0,sgny,0),my_element->faceIndex(0,0,sgnz) };
      int e_connect[3] = { my_element->edgeIndex(sgnx,sgny,0),my_element->edgeIndex(sgnx,0,sgnz),my_element->edgeIndex(0,sgny,sgnz) };
      // Find faces that do not connect this vertex
      int fx_neg = my_element->faceIndex(-1*sgnx,0,0);
      int fy_neg = my_element->faceIndex(0,-1*sgny,0);
      int fz_neg = my_element->faceIndex(0,0,-1*sgnz);

      double scaling_factor = my_element->lambda(fx_neg,*v) * my_element->lambda(fy_neg,*v) * my_element->lambda(fz_neg,*v);

      // If num_face_dofs>0, first evaluate the face basis function at 
      // all the face nodes on the three faces containing this vertex
      std::vector<double> eval_face_nodes;
      eval_face_nodes.clear();

      if (num_face_dofs>0) {
        for (int n=0; n<3; n++) {
          for (int iDoF=0; iDoF<num_face_dofs/6; iDoF++) {
            Point* node=my_ds_space->faceDoFPtr(my_element->faceGlobal(f_connect[n]),iDoF);
            eval_face_nodes.push_back( my_element->lambda(fx_neg,*node)*my_element->lambda(fy_neg,*node)*my_element->lambda(fz_neg,*node) / scaling_factor );
          }          
        }
      }

      // Calculate and store the coefficients corresponding to edge basis functions
      // such that the sum of vertex varphi and the 'weighted' edge basis functions
      // are linear on each edge
      std::vector<double> normalization_coeffs[3];  
      for (int n=0; n<3; n++) {
        int nEdge = e_connect[n];

        std::vector<double> expected_value_vec(polynomial_degree-1,0);
        double* expected_value = expected_value_vec.data();
        std::vector<double> coeff_vec(polynomial_degree-1,0);
        double* coeff =coeff_vec.data();

        for (int iRow=0; iRow<polynomial_degree-1; iRow++) {
          expected_value[iRow] = (v_index[2-n]<0)? 1-(double)(iRow+1)/polynomial_degree : (double)(iRow+1)/polynomial_degree;        
          expected_value[iRow] -= my_element->lambda(fx_neg,edge_nodes[nEdge*(polynomial_degree-1)+iRow])*my_element->lambda(fy_neg,edge_nodes[nEdge*(polynomial_degree-1)+iRow])*my_element->lambda(fz_neg,edge_nodes[nEdge*(polynomial_degree-1)+iRow]) / scaling_factor;
        }

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, polynomial_degree-1, 1,
              polynomial_degree-1, 1, edgecheby_eval_mat_inv[nEdge].data(), polynomial_degree-1, expected_value, 1,
              0, coeff, 1);
        normalization_coeffs[n] = coeff_vec;
      }

      for (int pt_index = 0; pt_index < num_pts; pt_index++) { 
        double value = my_element->lambda(fx_neg,pt[pt_index])*my_element->lambda(fy_neg,pt[pt_index])*my_element->lambda(fz_neg,pt[pt_index]) / scaling_factor;
        Tensor1 gradvalue = my_element->dLambda(fx_neg) * my_element->lambda(fy_neg,pt[pt_index]) * my_element->lambda(fz_neg,pt[pt_index])
                          + my_element->lambda(fx_neg,pt[pt_index]) * my_element->dLambda(fy_neg) * my_element->lambda(fz_neg,pt[pt_index])
                          + my_element->lambda(fx_neg,pt[pt_index]) * my_element->lambda(fy_neg,pt[pt_index]) * my_element->dLambda(fz_neg);
        gradvalue /= scaling_factor;

        // eliminate face dofs
        if (num_face_dofs>0) {
          for (int n=0; n<3; n++) {
            for (int iDoF=0; iDoF<num_face_dofs/6; iDoF++) {
              value -= eval_face_nodes[n*num_face_dofs/6+iDoF] * value_n[pt_index*nDoFs()+nVertexDoFs()+nEdgeDoFs()+f_connect[n]*num_face_dofs/6+iDoF];
              gradvalue -= eval_face_nodes[n*num_face_dofs/6+iDoF] * gradvalue_n[pt_index*nDoFs()+nVertexDoFs()+nEdgeDoFs()+f_connect[n]*num_face_dofs/6+iDoF];
            }
          }   
        }

        // linearize on the edges
        for (int n=0; n<3; n++) { 
          for (int iDoF=0; iDoF<polynomial_degree-1; iDoF++) {
            value += normalization_coeffs[n][iDoF] * value_n[pt_index*nDoFs()+nVertexDoFs()+e_connect[n]*(polynomial_degree-1)+iDoF];
            gradvalue += normalization_coeffs[n][iDoF] * gradvalue_n[pt_index*nDoFs()+nVertexDoFs()+e_connect[n]*(polynomial_degree-1)+iDoF];
          }
        }
        value_n[pt_index*nDoFs()+nVertex] = value;
        gradvalue_n[pt_index*nDoFs()+nVertex] = gradvalue;
      }   
    }
  }
}


// Eval for DirectSerendipityFE

void DirectSerendipityFE::eval(const Point* pt, double* result, Tensor1* gradResult, int num_pts,
			       double* vertex_dofs, double* edge_dofs, double* face_dofs, double* cell_dofs) {
  initBasis(pt,num_pts);

  for(int n=0; n<num_pts; n++) {
    result[n] = 0;
    gradResult[n].set(0,0);

    for(int i=0; i<8; i++) {
      result[n] += vertex_dofs[i]*basis(i,n);
      gradResult[n] += vertex_dofs[i]*gradVertexBasis(i,n);
    } 

    for(int k=0; k<12*(polynomial_degree-1); k++) {
      result[n] += edge_dofs[k]*edgeBasis(k,n);
      gradResult[n] += edge_dofs[k]*gradEdgeBasis(k,n);
    }

    for(int k=0; k<nFaceDoFs(); k++) {
      result[n] += face_dofs[k]*faceBasis(k,n);
      gradResult[n] += face_dofs[k]*gradFaceBasis(k,n);
    }

    for(int k=0; k<nCellDoFs(); k++) {
      result[n] += cell_dofs[k]*cellBasis(k,n);
      gradResult[n] += cell_dofs[k]*gradCellBasis(k,n);
    }
  }
}

void DirectSerendipityFE::eval(const Point* pt, double* result, int num_pts, 
				 double* vertex_dofs, double* edge_dofs, double* face_dofs, double* cell_dofs) {
  initBasis(pt,num_pts);

  for(int n=0; n<num_pts; n++) {
    result[n] = 0;

    for(int i=0; i<8; i++) {
      result[n] += vertex_dofs[i]*vertexBasis(i,n);
    } 

    for(int k=0; k<12*(polynomial_degree-1); k++) {
      result[n] += edge_dofs[k]*edgeBasis(k,n);
    }

    for(int k=0; k<nFaceDoFs(); k++) {
      result[n] += face_dofs[k]*faceBasis(k,n);
    }

    for(int k=0; k<nCellDoFs(); k++) {
      result[n] += cell_dofs[k]*cellBasis(k,n);
    }
  }
}

// Output for DirectSerendipityFE

void DirectSerendipityFE::write_raw(std::ofstream& fout) const {
  fout << "    DIRECT SERENDIPITY FE\n";
  fout << "    my_ds_space       = " << my_ds_space << "\n";
  fout << "    polynomial_degree = " <<  polynomial_degree << "\n";
  fout << "    my_element: " << my_element << "\n";
  my_element -> write_raw(fout);
  
}

int DirectSerendipityFE::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

void DirectSerendipityFE::write_tecplot(std::ofstream* fout, std::ofstream* fout_grad,
        int num_pts_x, int num_pts_y, int num_pts_z,
        double* vertex_dofs, double* edge_dofs, double* face_dofs, double* cell_dofs) {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;
  if(num_pts_z <= 1) num_pts_z = 2;

  // Determine mesh of points
  double x_min=my_element->vertexPtr(0)->val(0);
  double x_max=x_min;
  double y_min=my_element->vertexPtr(0)->val(1);
  double y_max=y_min;
  double z_min=my_element->vertexPtr(0)->val(2); 
  double z_max=z_min;

  for (int iVertex=1; iVertex<8; iVertex++) {
    if (my_element->vertexPtr(iVertex)->val(0) < x_min) x_min = my_element->vertexPtr(iVertex)->val(0);
    if (my_element->vertexPtr(iVertex)->val(1) < y_min) y_min = my_element->vertexPtr(iVertex)->val(1);
    if (my_element->vertexPtr(iVertex)->val(2) < z_min) z_min = my_element->vertexPtr(iVertex)->val(2);
    if (my_element->vertexPtr(iVertex)->val(0) > x_max) x_max = my_element->vertexPtr(iVertex)->val(0);
    if (my_element->vertexPtr(iVertex)->val(1) > y_max) y_max = my_element->vertexPtr(iVertex)->val(1);
    if (my_element->vertexPtr(iVertex)->val(2) > z_max) z_max = my_element->vertexPtr(iVertex)->val(2);
  }

  double dx = (x_max - x_min)/(num_pts_x-1);
  double dy = (y_max - y_min)/(num_pts_y-1);
  double dz = (z_max - z_min)/(num_pts_z-1);

  Point* pts = new Point[num_pts_x*num_pts_y*num_pts_z];
  Point* vertices = new Point[8];

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      for (int k=0; k<num_pts_z; k++) {
        pts[k + num_pts_z*j + num_pts_z*num_pts_y*i].set(x_min+i*dx, y_min+j*dy, z_min+k*dz);
      }
    }
  }

  for(int i=0; i<8; i++) {
    vertices[i] = *my_element->vertexPtr(i);
  }

  // Evaluate
  double result[num_pts_x*num_pts_y*num_pts_z];
  double vertexResult[8];
  Tensor1* gradResult = new Tensor1[num_pts_x*num_pts_y*num_pts_z];
  Tensor1* vertexGradResult = new Tensor1[8];
  eval(pts, result, gradResult, num_pts_x*num_pts_y*num_pts_z, vertex_dofs, edge_dofs, face_dofs, cell_dofs);
  eval(vertices, vertexResult, vertexGradResult, 8, vertex_dofs, edge_dofs, face_dofs, cell_dofs);

  // Modify the evaluation of the points outside the element to be zero
  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      for (int k=0; k<num_pts_z; k++) {
        if (!my_element->isInElement(pts[k + num_pts_z*j + num_pts_z*num_pts_y*i])) {
          result[k + num_pts_z*j + num_pts_z*num_pts_y*i] = 0;
          gradResult[k + num_pts_z*j + num_pts_z*num_pts_y*i].set(0,0,0);
        }
      }
    }
  }  

  // Write file  
  *fout << "TITLE = \"Direct serendipity array output\"\n";
  *fout << "VARIABLES = \"X\", \"Y\", \"Z\", \"value\"\n";
  *fout << "ZONE I=" << num_pts_x << ", J=" << num_pts_y << ", K=" << num_pts_z << ", DATAPACKING=POINT\n";

  if(fout_grad) {
    *fout_grad << "TITLE = \"Direct serendipity array output\"\n";
    *fout_grad << "VARIABLES = \"X\", \"Y\", \"Z\", \"gradvalue[1]\", \"gradvalue[2]\"\n";
    *fout_grad << "ZONE I=" << num_pts_x << ", J=" << num_pts_y << ", K=" << num_pts_z << ", DATAPACKING=POINT\n";
  }

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      for (int k=0; k<num_pts_z; k++) {
        int ind = k + num_pts_z*j + num_pts_z*num_pts_y*i;
        *fout << pts[ind][0] << " " << pts[ind][1] << " " << pts[ind][2] << " " << result[ind] << "\n";
        if(fout_grad) *fout_grad << pts[ind][0] << " " << pts[ind][1] << " " << pts[ind][2] << " " << gradResult[ind][0] << " " << gradResult[ind][1] << "\n";
      }
    }
  }


  *fout << "ZONE NODES=" << 8 << ", ELEMENTS=" << 1 << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n";
  for(int i=0; i<8; i++) {
    *fout <<my_element->vertexPtr(i)->val(0) << " ";
    if(fout_grad) *fout_grad << my_element->vertexPtr(i)->val(0) << " ";
  }
  *fout << "\n";
  if(fout_grad) *fout_grad << "\n";

  for(int i=0; i<8; i++) {
    *fout << my_element->vertexPtr(i)->val(1) << " ";
    if(fout_grad) *fout_grad << my_element->vertexPtr(i)->val(1) << " ";
  }
  *fout << "\n";
  if(fout_grad) *fout_grad << "\n";

  for(int i=0; i<8; i++) {
    *fout << my_element->vertexPtr(i)->val(2) << " ";
    if(fout_grad) *fout_grad << my_element->vertexPtr(i)->val(2) << " ";
  }
  *fout <<  "\n";
  if(fout_grad) *fout_grad << "\n";

  for(int i=0; i<8; i++) {
    *fout << vertexResult[i] << " ";
  }
  *fout <<  "\n";
  if(fout_grad) *fout_grad << "\n";


  *fout << my_element->vertexGlobal(4) << " " << my_element->vertexGlobal(6) << " " << my_element->vertexGlobal(2) << " " << my_element->vertexGlobal(0) << " "
        << my_element->vertexGlobal(5) << " " << my_element->vertexGlobal(7) << " " << my_element->vertexGlobal(3) << " " << my_element->vertexGlobal(1) << " ";
  *fout << "\n";
  if(fout_grad) {
  *fout_grad << my_element->vertexGlobal(4) << " " << my_element->vertexGlobal(6) << " " << my_element->vertexGlobal(2) << " " << my_element->vertexGlobal(0) << " "
              << my_element->vertexGlobal(5) << " " << my_element->vertexGlobal(7) << " " << my_element->vertexGlobal(3) << " " << my_element->vertexGlobal(1) << " ";
  *fout_grad << "\n";
  }


  delete[] gradResult;
  delete[] vertexGradResult;
  delete[] pts;
  delete[] vertices;
}

int DirectSerendipityFE::write_tecplot(std::string& filename, std::string& filename_grad,
        int num_pts_x, int num_pts_y, int num_pts_z,
        double* vertex_dofs, double* edge_dofs, double* face_dofs, double* cell_dofs) {
  std::ofstream fout(filename+".dat");
  if( !fout ) return 1;
  std::ofstream fout_grad(filename_grad+".dat");
  if( !fout_grad ) return 1;
  write_tecplot(&fout, &fout_grad, num_pts_x, num_pts_y, num_pts_z, vertex_dofs, edge_dofs, face_dofs, cell_dofs);
  return 0;
}