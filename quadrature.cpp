#include <cmath>

#include "quadrature.h"
#include "debug.h"
#include <assert.h>

using namespace quadrature;

////////////////////////////////////////////////////////////////////////////////
// class Quadrature

void Quadrature::set_rule(int desired_dop) {
  my_desired_dop = desired_dop;

  my_rule = ruleForTetra.size()-1;
  for(unsigned long int i=0; i<ruleForTetra.size(); i++) {
    if(ruleForTetra[i].dop >= desired_dop) { my_rule = i; break; }
  }
  my_dop  = ruleForTetra[my_rule].dop;
  num_pts_ref = ruleForTetra[my_rule].num;

  if(my_pts_ref) delete[] my_pts_ref;
  if(my_wts_ref) delete[] my_wts_ref;
  my_pts_ref = new Point[num_pts_ref];
  my_wts_ref = new double[num_pts_ref];

  // set them
  for (int iPt = 0; iPt < num_pts_ref; iPt++) {
    my_pts_ref[iPt] = ruleForTetra[my_rule].pts[iPt];
    my_wts_ref[iPt] = ruleForTetra[my_rule].wts[iPt];
  }
}

void Quadrature::set_element(int supplement_type, int refinement_level, hexamesh::Element* element) {
  my_element = element;
  if(!element) return;

  std::vector<std::vector<Point>> tetrahedra = refine(refinement_level, *my_element->subtetrahedra(std::min(1,supplement_type)));

  int num_tetrahedra = tetrahedra.size();
  num_pts = num_tetrahedra*num_pts_ref;

  // Quadrature points and weights

  if(my_pts) delete[] my_pts;
  my_pts = new Point[num_pts];

  if(my_wts) delete[] my_wts;
  my_wts = new double[num_pts];

  for(int i=0; i<num_tetrahedra; i++) {
    // v0 = (0,0,0)
    Point v1(tetrahedra[i][1]); v1 -= tetrahedra[i][0];
    Point v2(tetrahedra[i][2]); v2 -= tetrahedra[i][0];
    Point v3(tetrahedra[i][3]); v3 -= tetrahedra[i][0];


    Tensor2 mappingMatrix(v1[0], v2[0], v3[0], v1[1], v2[1], v3[1], v1[2], v2[2], v3[2]);
    double jacobian = std::abs(mappingMatrix.determinant())/6;

    int kk = i*num_pts_ref;
    for(int j=0; j<num_pts_ref; j++) {
      mappingMatrix.mult(my_pts_ref[j], my_pts[kk+j]);
      my_pts[kk+j] += tetrahedra[i][0];
      my_wts[kk+j] = my_wts_ref[j]*jacobian;
    }
  }
} 


std::vector<std::vector<Point>> Quadrature::refine(int level, const std::vector<std::vector<int>>& vertices) {
  std::vector<std::vector<Point>> my_vertices; my_vertices.clear();

  for (unsigned int iTetra = 0; iTetra < vertices.size(); iTetra++) {
    std::vector<Point> this_tetra(4);
    for (int i = 0; i < 4; i++) {
      this_tetra[i].set(*my_element->vertexPtr(vertices[iTetra][i]));
    }
    my_vertices.push_back(this_tetra);
  }

  std::vector<std::vector<Point>> my_old_vertices; my_old_vertices.clear();

  for (int iLevel=0; iLevel<level; iLevel++) {
    // Store all the data of my_vertices in my_old_vertices
    my_old_vertices = my_vertices;

    // Set up my_vertices with vertices of next level
    my_vertices.clear();
    for (unsigned int iTetra = 0; iTetra < my_old_vertices.size(); iTetra++) {
      // Calculate all the midpoints at once
      Point v01 = (my_old_vertices[iTetra][0]+my_old_vertices[iTetra][1])/2;
      Point v02 = (my_old_vertices[iTetra][0]+my_old_vertices[iTetra][2])/2;
      Point v03 = (my_old_vertices[iTetra][0]+my_old_vertices[iTetra][3])/2;
      Point v12 = (my_old_vertices[iTetra][1]+my_old_vertices[iTetra][2])/2;
      Point v13 = (my_old_vertices[iTetra][1]+my_old_vertices[iTetra][3])/2;
      Point v23 = (my_old_vertices[iTetra][2]+my_old_vertices[iTetra][3])/2;

      // t0
      my_vertices.push_back(std::vector<Point>{my_old_vertices[iTetra][0], v01, v02, v03});
      // t1
      my_vertices.push_back(std::vector<Point>{v01, my_old_vertices[iTetra][1], v12, v13});
      // t2
      my_vertices.push_back(std::vector<Point>{v02, v12, my_old_vertices[iTetra][2], v23});
      // t3
      my_vertices.push_back(std::vector<Point>{v03, v13, v23, my_old_vertices[iTetra][3]});
      // t4 (shares a common face with t0)
      my_vertices.push_back(std::vector<Point>{v12, v01, v02, v03});
      // t5 (shares a common face with t1)
      my_vertices.push_back(std::vector<Point>{v01, v03, v12, v13});
      // t6 (shares a common face with t2)
      my_vertices.push_back(std::vector<Point>{v02, v12, v03, v23});
      // t7 (shares a common face with t3)
      my_vertices.push_back(std::vector<Point>{v03, v13, v23, v12});
    }
  }
  return my_vertices;
}

Quadrature::~Quadrature() {
  if(my_pts) delete[] my_pts;
  if(my_wts) delete[] my_wts;
  if(my_pts_ref) delete[] my_pts_ref;
  if(my_wts_ref) delete[] my_wts_ref;
}

// Test degree of precision on a mesh over the domain [0,1]^2
void quadrature::testQuadrature(hexamesh::HexaMesh* mesh, int supplement_type, int refinement_level, double eps, int toDOP) {
  auto f = [](Point& x, int i, int j, int k) { return std::pow(x[0],i)*std::pow(x[1],j)*std::pow(x[2],k); };
  auto trueIntegF = [](int i, int j, int k) { return 1/double(i+1)*1/double(j+1)*1/double(k+1); };

  Quadrature quadrature;

  std::cout << std::endl;
  for(int testDOP=2; testDOP<=toDOP; testDOP++) {
    quadrature.setRule(testDOP);
    if(testDOP != quadrature.dop() ) continue;
    std::cout << "DOP = " << testDOP << "\n";

    for (int i=0; i<=testDOP; i++) {
      for(int j=0; j<=testDOP-i; j++) {
        for(int k=0; k<=testDOP-i-j; k++) {
          double full_integ = 0;
          for(int iElem=0; iElem<mesh->nElements(); iElem++) {
            quadrature.setElement(supplement_type, refinement_level, mesh->elementPtr(iElem));

            double integ = 0;
            for(int l=0; l<quadrature.num(); l++) {
            integ += f(quadrature.pts()[l],i,j,k)*quadrature.wt(l);
            }
          full_integ += integ;
          }
          double true_integ = trueIntegF(i,j,k);
          double err = std::abs(true_integ - full_integ) / true_integ;
          if(err > eps) std::cout << "  i = " << i << ", j = " << j << ", k = " << k << ", err = " << err << "\n";
        }
      }
    }
  }
  std::cout << std::endl;
};
