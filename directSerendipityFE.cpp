#include <cmath>
#include <iostream>
#include <vector>

#include <complex.h>
#include "lapacke.h"
#include <stdio.h>
#include <assert.h>

using namespace std;

#include "Utilities/debug.h"
#include "Mesh/hexaMesh.h"
using namespace hexamesh;
#include "directSerendipity.h"
using namespace directserendipity;

////////////////////////////////////////////////////////////////////////////////
// Class DirectSerendipityFE

void DirectSerendipityFE::set_directserendipityfe(DirectSerendipity* dsSpace, Element* element) {
  my_ds_space = dsSpace;
  my_element = element;

  polynomial_degree = my_ds_space->degPolyn();
  deg_face_poly = std::max(0,polynomial_degree-4);
  num_face_dofs = (deg_face_poly=0)? 0:3*(polynomial_degree-2)*(polynomial_degree-3);
  deg_cell_poly = std::max(0,polynomial_degree-6);
  num_cell_dofs = (deg_cell_poly=0)? 0:(polynomial_degree-3)*(polynomial_degree-4)*(polynomial_degree-5)/6;
  num_dofs = 8 + 12*(polynomial_degree-1) + num_face_dofs + num_cell_dofs; // Redundant with above

  // Set up higher order element (if nescessary)
  if(polynomial_degree<3) {
    if(one_element_mesh) delete one_element_mesh;
    one_element_mesh = new hexamesh::HexaMesh(my_element);
    if(high_order_ds_space) delete high_order_ds_space;
    high_order_ds_space = new DirectSerendipity(3,my_ds_space->supplementType(),one_element_mesh);
  }
}

DirectSerendipityFE::~DirectSerendipityFE() {
  if(high_order_ds_space) delete high_order_ds_space;
  if(one_element_mesh) delete one_element_mesh;
  if(value_n) delete[] value_n;
  if (gradvalue_n) delete[] gradvalue_n;
}

void DirectSerendipityFE::initBasis(const Point* pt, int num_pts) {
  return;
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