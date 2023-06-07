#include <cmath>
#include <iostream>
#include <vector>

#include <stdio.h>
#include <assert.h>
#include "lapacke.h"

using namespace std;

#include "Utilities/debug.h"
#include "Mesh/hexaMesh.h"
using namespace hexamesh;
#include "directSerendipity.h"
using namespace directserendipity;
#include "quadrature.h"
using namespace quadrature;

// fact
int fact(int n) {
   if (n == 0 || n == 1)
   return 1;
   else
   return n * fact(n - 1);
}

////////////////////////////////////////////////////////////////////////////////
// class DirectSerendipityArray

void DirectSerendipityArray::set_directserendipityarray(DirectSerendipity* dsSpace) {
  my_ds_space = dsSpace;
  num_dofs = my_ds_space->nDoFs();
  
  if(the_array) delete[] the_array;
  the_array = new double[num_dofs];
}

DirectSerendipityArray::~DirectSerendipityArray() {
  if(the_array) delete[] the_array;
}

void DirectSerendipityArray::eval(const Point* pts, double* result,
				  Tensor1* gradResult, int num_pts) const {
  bool ptEvaluated[num_pts];
  for(int i=0; i<num_pts; i++) ptEvaluated[i] = false;


  // Loop through the elements
  std::vector<Point> elementPts;
  std::vector<int> elementPtsIndex;
  for(int iElement=0; iElement < my_ds_space->my_mesh->nElements(); iElement++) {
    DirectSerendipityFE* finiteElement = &(my_ds_space->the_ds_elements[iElement]);
    Element* element = finiteElement->elementPtr();

    // Set list of points in element
    elementPts.clear();
    elementPtsIndex.clear();
    for(int i=0; i<num_pts; i++) {
      if(ptEvaluated[i]) continue;
      
      if(element->isInElement(pts[i])) {
        elementPts.push_back(pts[i]);
        elementPtsIndex.push_back(i);
        ptEvaluated[i] = true;
      }
    }
    if(elementPts.size() == 0) continue;

    // Set DoFs for element
    double vertex_dofs[8];
    double edge_dofs[12*(finiteElement->degPolyn() - 1)];
    double face_dofs[finiteElement->nFaceDoFs()];
    double cell_dofs[finiteElement->nCellDoFs()];

    for (int i=0; i<8; i++) {
      int iVertexDoF = element->vertexGlobal(i);
      vertex_dofs[i] = the_array[iVertexDoF];
    }

    for (int i=0; i<12; i++) {
      int iEdge = element->edgeGlobal(i);
      int iEdgeDoF = my_ds_space->nVertexDoFs() + iEdge*(finiteElement->degPolyn() - 1);
      for(int j=0; j<(finiteElement->degPolyn() - 1); j++) {
	      edge_dofs[i*(finiteElement->degPolyn() - 1) + j] = the_array[iEdgeDoF+j];
      }
    }

    for (int i=0; i<6; i++) {
      int iFace = element->faceGlobal(i);
      int iFaceDoF = my_ds_space->nVertexDoFs() + my_ds_space->nEdgeDoFs() + iFace*my_ds_space->nSingleFaceDoFs();
      for (int j=0; j<my_ds_space->nSingleFaceDoFs(); j++) {
        face_dofs[i*my_ds_space->nSingleFaceDoFs() + j] = the_array[iFaceDoF+j];
      }
    }

    {
      int iCell = element->meshIndex();
      int iCellDoF = my_ds_space->nVertexDoFs() + my_ds_space->nEdgeDoFs() + my_ds_space->nFaceDoFs() + iCell*finiteElement->nCellDoFs();
      for (int j=0; j<finiteElement->nCellDoFs(); j++) {
        cell_dofs[j] = the_array[iCellDoF+j];
      }
    }

    // Evaluate array at points on element
    double elementResult[elementPts.size()];
    Tensor1* elementGradResult = new Tensor1[elementPts.size()];
    finiteElement->eval(elementPts.data(), elementResult, elementGradResult, elementPts.size(),
			vertex_dofs, edge_dofs, face_dofs, cell_dofs);
    // Place results in global array
    for(unsigned long int i=0; i<elementPts.size(); i++) {
      result[elementPtsIndex[i]] = elementResult[i];
      gradResult[elementPtsIndex[i]] = elementGradResult[i];
    }

    delete[] elementGradResult;
  }

  // Attend to unset points (outside the mesh)
  for(int i=0; i<num_pts; i++) {
    if(!ptEvaluated[i]) {
      result[i] = 0;
      gradResult[i].set(0,0);
    }
  }
  return;
}

void DirectSerendipityArray::eval(const Point& pt, double& result, Tensor1& gradResult) const {
  int iElement = my_ds_space->my_mesh->inElement(pt);
  if(iElement < 0) {
    result = 0;
    gradResult.set(0,0);
    return;
  }
  DirectSerendipityFE* elem = &(my_ds_space->the_ds_elements[iElement]);

  // Set DoFs for element
  double vertex_dofs[8];
  double edge_dofs[12*(elem->degPolyn() - 1)];
  double face_dofs[elem->nFaceDoFs()];
  double cell_dofs[elem->nCellDoFs()];

  for (int i=0; i<8; i++) {
    int iVertexDoF = elem->elementPtr()->vertexGlobal(i);
    vertex_dofs[i] = the_array[iVertexDoF];
  }

  for (int i=0; i<12; i++) {
    int iEdge = elem->elementPtr()->edgeGlobal(i);
    int iEdgeDoF = my_ds_space->nVertexDoFs() + iEdge*(elem->degPolyn() - 1);
    for(int j=0; j<(elem->degPolyn() - 1); j++) {
      edge_dofs[i*(elem->degPolyn() - 1) + j] = the_array[iEdgeDoF+j];
    }
  }

  for (int i=0; i<6; i++) {
    int iFace = elem->elementPtr()->faceGlobal(i);
    int iFaceDoF = my_ds_space->nVertexDoFs() + my_ds_space->nEdgeDoFs() + iFace*my_ds_space->nSingleFaceDoFs();
    for (int j=0; j<my_ds_space->nSingleFaceDoFs(); j++) {
      face_dofs[i*my_ds_space->nSingleFaceDoFs() + j] = the_array[iFaceDoF+j];
    }
  }

  {
    int iCell = elem->elementPtr()->meshIndex();
    int iCellDoF = my_ds_space->nVertexDoFs() + my_ds_space->nEdgeDoFs() + my_ds_space->nFaceDoFs() + iCell*elem->nCellDoFs();
    for (int j=0; j<elem->nCellDoFs(); j++) {
      cell_dofs[j] = the_array[iCellDoF+j];
    }
  }

  // Evaluate
  elem->eval(pt, result, gradResult, vertex_dofs, edge_dofs, face_dofs, cell_dofs);
  return;
 }

double DirectSerendipityArray::eval(const Point& pt) const {
  int iElement = my_ds_space->my_mesh->inElement(pt);
  if(iElement < 0) return 0;

  DirectSerendipityFE* elem = &(my_ds_space->the_ds_elements[iElement]);

  // Set DoFs for element
  double vertex_dofs[8];
  double edge_dofs[12*(elem->degPolyn() - 1)];
  double face_dofs[elem->nFaceDoFs()];
  double cell_dofs[elem->nCellDoFs()];

  for (int i=0; i<8; i++) {
    int iVertexDoF = elem->elementPtr()->vertexGlobal(i);
    vertex_dofs[i] = the_array[iVertexDoF];
  }

  for (int i=0; i<12; i++) {
    int iEdge = elem->elementPtr()->edgeGlobal(i);
    int iEdgeDoF = my_ds_space->nVertexDoFs() + iEdge*(elem->degPolyn() - 1);
    for(int j=0; j<(elem->degPolyn() - 1); j++) {
      edge_dofs[i*(elem->degPolyn() - 1) + j] = the_array[iEdgeDoF+j];
    }
  }

  for (int i=0; i<6; i++) {
    int iFace = elem->elementPtr()->faceGlobal(i);
    int iFaceDoF = my_ds_space->nVertexDoFs() + my_ds_space->nEdgeDoFs() + iFace*my_ds_space->nSingleFaceDoFs();
    for (int j=0; j<my_ds_space->nSingleFaceDoFs(); j++) {
      face_dofs[i*my_ds_space->nSingleFaceDoFs() + j] = the_array[iFaceDoF+j];
    }
  }

  {
    int iCell = elem->elementPtr()->meshIndex();
    int iCellDoF = my_ds_space->nVertexDoFs() + my_ds_space->nEdgeDoFs() + my_ds_space->nFaceDoFs() + iCell*elem->nCellDoFs();
    for (int j=0; j<elem->nCellDoFs(); j++) {
      cell_dofs[j] = the_array[iCellDoF+j];
    }
  }

  // Evaluate
  return elem->eval(pt, vertex_dofs, edge_dofs, face_dofs, cell_dofs);
}

void DirectSerendipityArray::l2normError(double& l2Error, double& l2GradError, double& l2Norm, double& l2GradNorm, int refinement_level,
					 double (*referenceFcn)(double,double,double),
					 Tensor1 (*referenceGradFcn)(double,double,double)) {
  l2Error = 0, l2GradError = 0, l2Norm = 0, l2GradNorm = 0;
  Quadrature quadRule(8);

  for(int iElement=0; iElement < my_ds_space->mesh()->nElements(); iElement++) {
    DirectSerendipityFE* fePtr = my_ds_space->finiteElementPtr(iElement);
    quadRule.setElement(my_ds_space->supplementType(), refinement_level, fePtr->elementPtr());

    double result[quadRule.num()];
    Tensor1* gradResult = new Tensor1[quadRule.num()];

    eval(quadRule.pts(), result, gradResult,quadRule.num());

    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt).val(0);
      double y = quadRule.pt(iPt).val(1);
      double z = quadRule.pt(iPt).val(2);
      
      double diff = (referenceFcn == nullptr) ? result[iPt] : (result[iPt] - referenceFcn(x,y,z));
      Tensor1 diffGrad = (referenceGradFcn == nullptr) ? gradResult[iPt] : (gradResult[iPt] - referenceGradFcn(x,y,z));
      
      l2Error += pow(diff,2) * quadRule.wt(iPt);
      l2GradError += diffGrad * diffGrad * quadRule.wt(iPt);

      l2Norm += pow(result[iPt],2) * quadRule.wt(iPt);
      l2GradNorm += gradResult[iPt] * gradResult[iPt] * quadRule.wt(iPt);
    }
    delete[] gradResult;
  }
  
  l2Error = sqrt(l2Error);
  l2GradError = sqrt(l2GradError);
  l2Norm = sqrt(l2Norm);
  l2GradNorm = sqrt(l2GradNorm);

  return;
}

void DirectSerendipityArray::write_raw(std::ofstream& fout) const {
  for(int i=0; i<num_dofs; i++) {
    fout << the_array[i] << "  ";
    if( !(i%10) ) fout << "\n";
  }
};

int DirectSerendipityArray::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

void DirectSerendipityArray::write_tecplot_mesh(std::ofstream* fout, std::ofstream* fout_grad,
					       int num_pts_x, int num_pts_y, int num_pts_z) const {
  if(num_pts_x <= 1) num_pts_x = 2;
  if(num_pts_y <= 1) num_pts_y = 2;
  if(num_pts_z <= 1) num_pts_z = 2;

  // Determine mesh of points
  double xMin = my_ds_space->my_mesh->xMin();
  double xMax = my_ds_space->my_mesh->xMax();
  double yMin = my_ds_space->my_mesh->yMin();
  double yMax = my_ds_space->my_mesh->yMax();
  double zMin = my_ds_space->my_mesh->zMin();
  double zMax = my_ds_space->my_mesh->zMax();

  double dx = (xMax - xMin)/(num_pts_x-1);
  double dy = (yMax - yMin)/(num_pts_y-1);
  double dz = (zMax - zMin)/(num_pts_z-1);

  Point* pts = new Point[num_pts_x*num_pts_y*num_pts_z];
  
  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      for (int k=0; k<num_pts_z; k++) {
        pts[k + num_pts_z*j + num_pts_z*num_pts_y*i].set(xMin+i*dx, yMin+j*dy, zMin+k*dz);
      }
    }
  }


  // Evaluate
  double result[num_pts_x*num_pts_y*num_pts_z];
  Tensor1* gradResult = new Tensor1[num_pts_x*num_pts_y*num_pts_z];

  eval(pts, result, gradResult, num_pts_x*num_pts_y*num_pts_z);
  // Write file  
  *fout << "TITLE = \"Direct serendipity array output\"\n";
  *fout << "VARIABLES = \"X\", \"Y\", \"Z\", \"value\"\n";
  *fout << "ZONE I=" << num_pts_x << ", J=" << num_pts_y << ", K=" << num_pts_z << ", DATAPACKING=POINT\n";

  if(fout_grad) {
    *fout_grad << "TITLE = \"Direct serendipity array output\"\n";
    *fout_grad << "VARIABLES = \"X\", \"Y\", \"Z\", \"gradvalue[x]\", \"gradvalue[y]\", \"gradvalue[z]\"\n";
    *fout_grad << "ZONE I=" << num_pts_x << ", J=" << num_pts_y << ", K=" << num_pts_z << ", DATAPACKING=POINT\n";
  }

  for(int i=0; i<num_pts_x; i++) {
    for(int j=0; j<num_pts_y; j++) {
      for (int k=0; k<num_pts_z; k++) {
        int ind = k + num_pts_z*j + num_pts_z*num_pts_y*i;
        *fout << pts[ind][0] << " " << pts[ind][1] << " " << pts[ind][2] << " " << result[ind] << "\n";
        if(fout_grad) *fout_grad << pts[ind][0] << " " << pts[ind][1] << " " << pts[ind][2] << " " << gradResult[ind][0] << " " << gradResult[ind][1] << " " << gradResult[ind][2] << "\n";
      }
    }
  }

 
  delete[] pts;
  delete[] gradResult;

  return;
};

int DirectSerendipityArray::write_tecplot_mesh(std::string& filename, std::string& filename_grad,
					      int num_pts_x, int num_pts_y, int num_pts_z) const {
  std::ofstream fout(filename+".dat");
  if( !fout ) return 1;
  std::ofstream fout_grad(filename_grad+".dat");
  if( !fout_grad ) return 1;
  write_tecplot_mesh(&fout, &fout_grad, num_pts_x, num_pts_y, num_pts_z);
  return 0;
}

int DirectSerendipityArray::write_tecplot_mesh(std::string& filename,
					      int num_pts_x, int num_pts_y, int num_pts_z) const {
  std::ofstream fout(filename+".dat");
  if( !fout ) return 1;
  write_tecplot_mesh(&fout, nullptr, num_pts_x, num_pts_y, num_pts_z);
  return 0;
}


void DirectSerendipity::set_directserendipity(int polyDeg, int suppType, HexaMesh* mesh) {
  polynomial_degree = polyDeg;
  supplement_type = suppType;
  my_mesh = mesh;

  num_dofs_per_face = (polynomial_degree < 4)? 0 : (polynomial_degree-2)*(polynomial_degree-3)/2; 
  num_dofs_per_cell = (polynomial_degree < 6)? 0 : (polynomial_degree-3)*(polynomial_degree-4)*(polynomial_degree-5)/6;
  
  num_dofs = my_mesh->nVertices() + my_mesh->nEdges() * (polynomial_degree-1) + my_mesh->nFaces() * num_dofs_per_face + my_mesh->nElements() * num_dofs_per_cell;

  // set edge nodes
  if (polynomial_degree > 1) {
    if (edge_nodes) delete[] edge_nodes; 
    edge_nodes = new Point[my_mesh->nEdges()*(polynomial_degree-1)];

    for (int iEdge=0; iEdge<my_mesh->nEdges(); iEdge++) {
      Vertex* v0 = my_mesh->edgeVertexPtr(0,iEdge);
      Vertex* v1 = my_mesh->edgeVertexPtr(1,iEdge);
      Tensor1 tangent(*v1-*v0);
      double length = tangent.norm();
      tangent /= length;
      for (int nPt=0; nPt<polynomial_degree-1; nPt++) {
        edge_nodes[iEdge*(polynomial_degree-1)+nPt].set(*v0 + (nPt+1)*length*tangent/polynomial_degree);    
      }
    } 
  }

  // set face dofs
  if (num_dofs_per_face>0) {
    if (face_dofs) delete[] face_dofs;
    face_dofs = new Point[num_dofs_per_face*my_mesh->nFaces()];
    int node_num = 0;
    int degFacePoly = polynomial_degree-4;

    for (int iFace=0; iFace<my_mesh->nFaces(); iFace++) {
      Point center(*my_mesh->faceCenterPtr(iFace));
      if (degFacePoly > 0) {
        double dist = 0.5*my_mesh->faceRadius(iFace);
        Tensor1 x_vec = *my_mesh->faceVertexPtr(1,0,iFace)-*my_mesh->faceVertexPtr(0,0,iFace);
        x_vec /= x_vec.norm();
        Tensor1 y_vec = cross(cross(x_vec,Tensor1(*my_mesh->faceVertexPtr(0,1,iFace)-*my_mesh->faceVertexPtr(0,0,iFace))), x_vec);
        y_vec /= y_vec.norm();

        Point dx(std::sqrt(3)*dist/degFacePoly*x_vec);
        Point dy(0.5*std::sqrt(3)*dist/degFacePoly*x_vec+1.5*dist/degFacePoly*y_vec);

        Point basePoint(center-0.5*std::sqrt(3)*dist*x_vec-0.5*dist*y_vec);
        for(int j=0; j<=degFacePoly; j++) {
          Point nodePoint(basePoint);
          for(int i=0; i<=degFacePoly-j; i++) {
            face_dofs[node_num] = nodePoint;
            nodePoint += dx;
            node_num++;
          }    
          basePoint += dy;
        }
      } else {
        // Just take the center to be the only face dof
        face_dofs[iFace] = center;
      }
    }
  }

  // Define and assemble edge_cheby
  if(edge_cheby) delete[] edge_cheby;
  edge_cheby = new double[my_mesh->nEdges() * (polynomial_degree-1) * (polynomial_degree-1)];

  for (int iEdge=0; iEdge<my_mesh->nEdges(); iEdge++) {
      Vertex* v0 = my_mesh->edgeVertexPtr(0,iEdge);
      Vertex* v1 = my_mesh->edgeVertexPtr(1,iEdge);
      Tensor1 tangent(*v1-*v0);
      double length = tangent.norm();
      tangent /= length;
    for (int nPt=0; nPt<polynomial_degree-1; nPt++) {  
      for (int s=0; s<polynomial_degree-1; s++) {
        double x = Tensor1(edge_nodes[iEdge*(polynomial_degree-1)+nPt]-(*v0+*v1)/2)*tangent;

        double value=0;
        for (int t=0; t<=floor(s/2); t++) {
          value += fact(s)/fact(2*t)/fact(s-2*t) * pow(4*x*x/length/length-1,t) * pow(2*x/length,s-2*t);
        }
        value *= Tensor1(edge_nodes[iEdge*(polynomial_degree-1)+nPt]-*v0)*Tensor1(*v1-edge_nodes[iEdge*(polynomial_degree-1)+nPt]);
        edge_cheby[iEdge*(polynomial_degree-1)*(polynomial_degree-1) + nPt*(polynomial_degree-1) + s] = value / pow(length/2,2);
      }
    }
  }

  // ALLOCATE ELEMENTS

  if(the_ds_elements) delete[] the_ds_elements;
  the_ds_elements = new DirectSerendipityFE[my_mesh->nElements()];
  for(int iElement=0; iElement<my_mesh->nElements(); iElement++) { 
    the_ds_elements[iElement].set(this, my_mesh->elementPtr(iElement));
  }

  // store if the dof is on the boundary
  if (is_interior) delete[] is_interior;
  is_interior = new bool[num_dofs];
  int index = 0;

  // store the edge indices on the boundary
  if (bc_edge_index) delete[] bc_edge_index;
  bc_edge_index = new int[2*(mesh->xElements()+mesh->yElements())*mesh->zElements() + 2*(mesh->xElements()+mesh->zElements())*mesh->yElements() + 2*(mesh->yElements()+mesh->zElements())*mesh->xElements()];
  int bc_index = 0;

  // vertex DoFs
  for (int i=0; i<mesh->nVertices(); i++) {
    int vertex_pos[3];
    mesh->vertexPos(i,vertex_pos[0],vertex_pos[1],vertex_pos[2]);
    is_interior[index] = true;
    if (vertex_pos[0] == 0 || vertex_pos[0] == mesh->xElements()) is_interior[index] = false;
    if (vertex_pos[1] == 0 || vertex_pos[1] == mesh->yElements()) is_interior[index] = false;
    if (vertex_pos[2] == 0 || vertex_pos[2] == mesh->zElements()) is_interior[index] = false;
    index++;
  }

  // edge DoFs
  if (polynomial_degree>1) {
    for (int i=0; i<mesh->nEdges(); i++) {
      bool edge_is_interior = true;
      if (i<(mesh->xElements()+1)*(mesh->yElements()+1)*mesh->zElements()) {
        // z-dir edges
        int ix = mesh->edgeVertexPtr(0,i)->meshPos(0);
        int iy = mesh->edgeVertexPtr(0,i)->meshPos(1);
        if (ix == 0 || ix == mesh->xElements()) edge_is_interior = false;
        if (iy == 0 || iy == mesh->yElements()) edge_is_interior = false;
      } else if (i<(mesh->xElements()+1)*(mesh->yElements()+1)*mesh->zElements()+(mesh->xElements()+1)*mesh->yElements()*(mesh->zElements()+1)) {
        // y-dir edges
        int ix = mesh->edgeVertexPtr(0,i)->meshPos(0);
        int iz = mesh->edgeVertexPtr(0,i)->meshPos(2);
        if (ix == 0 || ix == mesh->xElements()) edge_is_interior = false;
        if (iz == 0 || iz == mesh->zElements()) edge_is_interior = false;
      } else {
        // x-dir edges
        int iy = mesh->edgeVertexPtr(0,i)->meshPos(1);
        int iz = mesh->edgeVertexPtr(0,i)->meshPos(2);
        if (iy == 0 || iy == mesh->yElements()) edge_is_interior = false;
        if (iz == 0 || iz == mesh->zElements()) edge_is_interior = false;
      }

      if (!edge_is_interior) {
        bc_edge_index[bc_index] = i;
        bc_index++;
      }

      for (int jDoF=0; jDoF<polynomial_degree-1; jDoF++) {
        is_interior[index] = edge_is_interior;
        index++;
      }
    }
  }


  // face DoFs
  if (num_dofs_per_face>0) {
    for (int i=0; i<mesh->nFaces(); i++) {
      bool face_is_interior = true;
      if (i < (mesh->xElements()+1)*mesh->yElements()*mesh->zElements()) {
        // yz-dir faces
        int ix = mesh->faceVertexPtr(0,0,i)->meshPos(0);
        if (ix == 0 || ix == mesh->xElements()) face_is_interior = false;
      } else if (i < (mesh->xElements()+1)*mesh->yElements()*mesh->zElements()+mesh->xElements()*(mesh->yElements()+1)*mesh->zElements()) {
        // xz-dir faces
        int iy = mesh->faceVertexPtr(0,0,i)->meshPos(1);
        if (iy == 0 || iy == mesh->yElements()) face_is_interior = false;
      } else {
        // xy - dir faces
        int iz = mesh->faceVertexPtr(0,0,i)->meshPos(2);
        if (iz == 0 || iz == mesh->zElements()) face_is_interior = false;
      }
      for (int jDoF=0; jDoF<num_dofs_per_face; jDoF++) {
        is_interior[index] = face_is_interior;
        index++;
      }
    }
  }

  // cell DoFs
  if (num_dofs_per_cell>0) {
    for (int i=0; i<mesh->nElements(); i++) {
      is_interior[index] = true;
      index++;      
    }
  }
}

DirectSerendipity::~DirectSerendipity() {
  if (edge_nodes) delete[] edge_nodes; 
  if (face_dofs) delete[] face_dofs;
  if (the_ds_elements) delete[] the_ds_elements;
  if (is_interior) delete[] is_interior;
  if (edge_cheby) delete[] edge_cheby;
  if (bc_edge_index) delete[] bc_edge_index;
};

void DirectSerendipity::bcModification(double* bc_vals) {
  int num_bc_vertices = 2 + 2*my_mesh->xElements()*my_mesh->yElements() + 2*my_mesh->xElements()*my_mesh->zElements() + 2*my_mesh->yElements()*my_mesh->zElements();
  int num_bc_edges = 2*(my_mesh->xElements()+my_mesh->yElements())*my_mesh->zElements() + 2*(my_mesh->xElements()+my_mesh->zElements())*my_mesh->yElements() + 2*(my_mesh->yElements()+my_mesh->zElements())*my_mesh->xElements();

  // vertex index mapping
  int ind = 0;
  int index_mapping[my_mesh->nVertices()];

  for (int i=0; i<my_mesh->nVertices(); i++) {
    if (!isInterior(i)) {
      // on boundary
      index_mapping[i] = ind;
      ind++;
    } else {
      // interior
      index_mapping[i] = -1;
    }
  }

  for (int iEdge=0; iEdge<num_bc_edges; iEdge++) {
    int real_index = bc_edge_index[iEdge];
    int starting_dof = num_bc_vertices + iEdge * (degPolyn()-1);

    double eval_v0 = bc_vals[index_mapping[my_mesh->edgeVertexPtr(0,real_index)->meshIndex()]];
    double eval_v1 = bc_vals[index_mapping[my_mesh->edgeVertexPtr(1,real_index)->meshIndex()]];

    std::vector<double> A_vec((degPolyn()-1)*(degPolyn()-1),0);
    double* A = A_vec.data();

    std::vector<double> rhs_vec((degPolyn()-1)*(degPolyn()-1),0);
    double* rhs = rhs_vec.data();

    
    for (int iRow=0; iRow<degPolyn()-1; iRow++) {
      rhs[iRow] = bc_vals[starting_dof+iRow] - (eval_v0 * (1-(double)(iRow+1)/degPolyn()) + eval_v1 * (double)(iRow+1)/degPolyn());
      for (int jCol=0; jCol<degPolyn()-1; jCol++) {
        A[iRow*(degPolyn()-1)+jCol] = edgeCheby(real_index, iRow, jCol);
      }        
    }

    lapack_int* ipiv; char norm = 'I';
    int size = degPolyn()-1;
    ipiv = (lapack_int*)malloc(size * sizeof(lapack_int));
    double anorm = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, degPolyn()-1, degPolyn()-1, A, degPolyn()-1);
    int ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, degPolyn()-1, 1, A, degPolyn()-1, ipiv, rhs, 1); //mat updated to be LU
    if(ierr) { // ?? what should we do ???
      std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
    }
    double rcond = 0;
    ierr = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, degPolyn()-1, A, degPolyn()-1, anorm, &rcond);

    for (int i=0; i<degPolyn()-1; i++) {
      bc_vals[starting_dof+i] = rhs[i];
    }
    free(ipiv);
  }
  return;
}

void DirectSerendipity::nodeModification(double* bc_vals) {
  for (int iEdge=0; iEdge<my_mesh->nEdges(); iEdge++) {
    int starting_dof = my_mesh->nVertices() + iEdge * (degPolyn()-1);

    double eval_v0 = bc_vals[my_mesh->edgeVertexPtr(0,iEdge)->meshIndex()];
    double eval_v1 = bc_vals[my_mesh->edgeVertexPtr(1,iEdge)->meshIndex()];

    std::vector<double> A_vec((degPolyn()-1)*(degPolyn()-1),0);
    double* A = A_vec.data();

    std::vector<double> rhs_vec((degPolyn()-1)*(degPolyn()-1),0);
    double* rhs = rhs_vec.data();

    
    for (int iRow=0; iRow<degPolyn()-1; iRow++) {
      rhs[iRow] = bc_vals[starting_dof+iRow] - (eval_v0 * (1-(double)(iRow+1)/degPolyn()) + eval_v1 * (double)(iRow+1)/degPolyn());
      for (int jCol=0; jCol<degPolyn()-1; jCol++) {
        A[iRow*(degPolyn()-1)+jCol] = edgeCheby(iEdge, iRow, jCol);
      }        
    }

    lapack_int* ipiv; char norm = 'I';
    int size = degPolyn()-1;
    ipiv = (lapack_int*)malloc(size * sizeof(lapack_int));
    double anorm = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, degPolyn()-1, degPolyn()-1, A, degPolyn()-1);
    int ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, degPolyn()-1, 1, A, degPolyn()-1, ipiv, rhs, 1); //mat updated to be LU
    if(ierr) { // ?? what should we do ???
      std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
    }
    double rcond = 0;
    ierr = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, degPolyn()-1, A, degPolyn()-1, anorm, &rcond);

    for (int i=0; i<degPolyn()-1; i++) {
      bc_vals[starting_dof+i] = rhs[i];
    }
    free(ipiv);
  }
  return;
}

void DirectSerendipity::write_raw(std::ofstream& fout) const {
  fout << "DIRECT SERENDIPITY SPACE\n";
  fout << "polynomial_degree = " <<  polynomial_degree << "\n";
  fout << "supplement_type = " <<  supplement_type << "\n";
  fout << "my_mesh           = " << my_mesh << "\n";
  fout << "num_dofs         = " << num_dofs << "\n";

  fout << "num_vertex_dofs         = " << nVertexDoFs() << "\n";
  fout << "num_edge_dofs         = " << nEdgeDoFs() << "\n";
  fout << "num_face_dofs         = " << nFaceDoFs() << "\n";
  fout << "num_cell_dofs         = " << nCellDoFs() << "\n";

  fout << "\nface dofs:\n";
  for (int iFace=0; iFace<my_mesh->nFaces(); iFace++) {
    fout << "\tface[" << iFace << "]\n";
    for (int iDoF=0; iDoF<nSingleFaceDoFs(); iDoF++) {
      fout << "\t\tDoF[" << iDoF << "]: (" << face_dofs[iFace*nSingleFaceDoFs()+iDoF] << ")\n";
    }
  }
}

int DirectSerendipity::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

int DirectSerendipity::write_matlab(std::string& filename) const {
  std::ofstream fout(filename + ".m");
  if( !fout ) return 1;

  // MESH
  // Define the color of faces
  std::string color = "[0.6471 0.8706 0.8941]";

  fout << "clf;\n";

  int i = 1; // face index

  for (int ix=0; ix<my_mesh->xElements()+1; ix++) {
    for (int iy=0; iy<my_mesh->yElements()+1; iy++) {
      for (int iz=0; iz<my_mesh->zElements()+1; iz++) {
        if (iy < my_mesh->yElements() && iz < my_mesh->zElements()) {
          fout << "x" << i << " = [" << my_mesh->vertexPtr(ix,iy,iz)->val(0) << " "
                                     << my_mesh->vertexPtr(ix,iy+1,iz)->val(0) << " "
                                     << my_mesh->vertexPtr(ix,iy+1,iz+1)->val(0) << " "
                                     << my_mesh->vertexPtr(ix,iy,iz+1)->val(0) << "];\n";
          fout << "y" << i << " = [" << my_mesh->vertexPtr(ix,iy,iz)->val(1) << " "
                                     << my_mesh->vertexPtr(ix,iy+1,iz)->val(1) << " "
                                     << my_mesh->vertexPtr(ix,iy+1,iz+1)->val(1) << " "
                                     << my_mesh->vertexPtr(ix,iy,iz+1)->val(1) << "];\n";
          fout << "z" << i << " = [" << my_mesh->vertexPtr(ix,iy,iz)->val(2) << " "
                                     << my_mesh->vertexPtr(ix,iy+1,iz)->val(2) << " "
                                     << my_mesh->vertexPtr(ix,iy+1,iz+1)->val(2) << " "
                                     << my_mesh->vertexPtr(ix,iy,iz+1)->val(2) << "];\n";
          i += 1;
        }

        if (ix < my_mesh->xElements() && iz < my_mesh->zElements()) {
          fout << "x" << i << " = [" << my_mesh->vertexPtr(ix,iy,iz)->val(0) << " "
                                     << my_mesh->vertexPtr(ix+1,iy,iz)->val(0) << " "
                                     << my_mesh->vertexPtr(ix+1,iy,iz+1)->val(0) << " "
                                     << my_mesh->vertexPtr(ix,iy,iz+1)->val(0) << "];\n";
          fout << "y" << i << " = [" << my_mesh->vertexPtr(ix,iy,iz)->val(1) << " "
                                     << my_mesh->vertexPtr(ix+1,iy,iz)->val(1) << " "
                                     << my_mesh->vertexPtr(ix+1,iy,iz+1)->val(1) << " "
                                     << my_mesh->vertexPtr(ix,iy,iz+1)->val(1) << "];\n";
          fout << "z" << i << " = [" << my_mesh->vertexPtr(ix,iy,iz)->val(2) << " "
                                     << my_mesh->vertexPtr(ix+1,iy,iz)->val(2) << " "
                                     << my_mesh->vertexPtr(ix+1,iy,iz+1)->val(2) << " "
                                     << my_mesh->vertexPtr(ix,iy,iz+1)->val(2) << "];\n";
          i += 1;
        }

        if (ix < my_mesh->xElements() && iy < my_mesh->yElements()) {
          fout << "x" << i << " = [" << my_mesh->vertexPtr(ix,iy,iz)->val(0) << " "
                                     << my_mesh->vertexPtr(ix+1,iy,iz)->val(0) << " "
                                     << my_mesh->vertexPtr(ix+1,iy+1,iz)->val(0) << " "
                                     << my_mesh->vertexPtr(ix,iy+1,iz)->val(0) << "];\n";
          fout << "y" << i << " = [" << my_mesh->vertexPtr(ix,iy,iz)->val(1) << " "
                                     << my_mesh->vertexPtr(ix+1,iy,iz)->val(1) << " "
                                     << my_mesh->vertexPtr(ix+1,iy+1,iz)->val(1) << " "
                                     << my_mesh->vertexPtr(ix,iy+1,iz)->val(1) << "];\n";
          fout << "z" << i << " = [" << my_mesh->vertexPtr(ix,iy,iz)->val(2) << " "
                                     << my_mesh->vertexPtr(ix+1,iy,iz)->val(2) << " "
                                     << my_mesh->vertexPtr(ix+1,iy+1,iz)->val(2) << " "
                                     << my_mesh->vertexPtr(ix,iy+1,iz)->val(2) << "];\n";
          i += 1;
        }
      }
    }
  }

  fout << "p = fill3(";

  for (int k=1; k<i; k++) {
    fout << "x" << k << ",y" << k << ",z" << k << "," << color;
    if (k<i-1) fout << ","; 
  }
  fout << ");\n";
  fout << "hold on;\n";
  fout << "set(p,'facealpha',.2);\n";
  
  // NODAL DOFs

  // Face nodes
  if (num_dofs_per_face>0) {
    fout << "scatter3([ ";
    for(int i=0; i<nFaceDoFs(); i++) {
      fout << face_dofs[i][0] << " ";
    }
    fout << "],[ ";
    for(int i=0; i<nFaceDoFs(); i++) {
      fout << face_dofs[i][1] << " ";
    }
    fout << "],[ ";
    for(int i=0; i<nFaceDoFs(); i++) {
      fout << face_dofs[i][2] << " ";
    }
    fout << "],'sb','filled')\n";

  }

  // CENTERS

  // Face centers
  fout << "scatter3([ ";
  for(int i=0; i<my_mesh->nFaces(); i++) {
    fout << my_mesh->faceCenterPtr(i)->val(0) << " ";
  }
  fout << "],[ ";
  for(int i=0; i<my_mesh->nFaces(); i++) {
    fout << my_mesh->faceCenterPtr(i)->val(1) << " ";
  }
  fout << "],[ ";
  for(int i=0; i<my_mesh->nFaces(); i++) {
    fout << my_mesh->faceCenterPtr(i)->val(2) << " ";
  }
  fout << "],'ok')\n";

  fout << "hold off;\n";
  return 0;
}