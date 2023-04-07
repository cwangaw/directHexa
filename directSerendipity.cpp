#include <cmath>
#include <iostream>
#include <vector>

#include <stdio.h>
#include <assert.h>


using namespace std;

#include "Utilities/debug.h"
#include "Mesh/hexaMesh.h"
using namespace hexamesh;
#include "directSerendipity.h"
using namespace directserendipity;
#include "quadrature.h"
using namespace quadrature;

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
  Quadrature quadRule(8,refinement_level);

  for(int iElement=0; iElement < my_ds_space->mesh()->nElements(); iElement++) {
    DirectSerendipityFE* fePtr = my_ds_space->finiteElementPtr(iElement);
    quadRule.setElement(my_ds_space->supplementType(), refinement_level, fePtr->elementPtr());
    
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt).val(0);
      double y = quadRule.pt(iPt).val(1);
      double z = quadRule.pt(iPt).val(2);
      
      double result; Tensor1 gradResult;
      eval(quadRule.pt(iPt), result, gradResult);
      
      double diff = (referenceFcn == nullptr) ? result : (result - referenceFcn(x,y,z));
      Tensor1 diffGrad = (referenceGradFcn == nullptr) ? gradResult : (gradResult - referenceGradFcn(x,y,z));
      
      l2Error += pow(diff,2) * quadRule.wt(iPt);
      l2GradError += diffGrad * diffGrad * quadRule.wt(iPt);

      l2Norm += pow(result,2) * quadRule.wt(iPt);
      l2GradNorm += gradResult * gradResult * quadRule.wt(iPt);
    }
  }
  
  l2Error = sqrt(l2Error);
  l2GradError = sqrt(l2GradError);
  l2Norm = sqrt(l2Norm);
  l2GradNorm = sqrt(l2GradNorm);
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

void DirectSerendipity::set_directserendipity(int polyDeg, int suppType, HexaMesh* mesh) {
  polynomial_degree = polyDeg;
  supplement_type = suppType;
  my_mesh = mesh;

  num_dofs_per_face = (polynomial_degree < 4)? 0 : (polynomial_degree-2)*(polynomial_degree-3)/2; 
  num_dofs_per_cell = (polynomial_degree < 6)? 0 : (polynomial_degree-3)*(polynomial_degree-4)*(polynomial_degree-5)/6;
  
  num_dofs = my_mesh->nVertices() + my_mesh->nEdges() * (polynomial_degree-1) + my_mesh->nFaces() * num_dofs_per_face + my_mesh->nElements() * num_dofs_per_cell;

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
/*
  // set cell dofs
  if (num_dofs_per_cell>0) {
    if (cell_dofs) delete[] cell_dofs;
    cell_dofs = new Point[num_dofs_per_cell*my_mesh->nElements()];
    int node_num = 0;
    int degCellPoly = polynomial_degree-6;

    for (int iElement=0; iElement<my_mesh->nElements(); iElement++) {
      Point center(*my_mesh->elementCenterPtr(iElement));
      if (degCellPoly > 0) {
        double dist = 0.5*my_mesh->elementRadius(iElement);

        Point dx(std::sqrt(double(8)/3)*dist/degCellPoly,0,0);
        Point dy(std::sqrt(double(2)/3)*dist/degCellPoly,std::sqrt(2)*dist/degCellPoly,0);
        Point dz(std::sqrt(double(2)/3)*dist/degCellPoly,std::sqrt(2)/3*dist/degCellPoly,(double(4)/3)*dist/degCellPoly);

        Point zBasePoint(-std::sqrt(double(2)/3)*dist,-std::sqrt(2)/3*dist,-(double(1)/3)*dist);
        zBasePoint += center;

        for (int k=0; k<=degCellPoly; k++) {
          Point yBasePoint(zBasePoint);
          for (int j=0; j<=degCellPoly-k; j++) {
            Point nodePoint(yBasePoint);
            for (int i=0; i<=degCellPoly-k-j; i++) {
              cell_dofs[node_num] = nodePoint;
              nodePoint += dx;
              node_num++;              
            }
            yBasePoint += dy;
          }
          zBasePoint += dz;
        }
      } else {
        // Just take the center to be the only face dof
        cell_dofs[iElement] = center;
      }
    }
  }
*/  
  // ALLOCATE ELEMENTS
  if(the_ds_elements) delete[] the_ds_elements;
  the_ds_elements = new DirectSerendipityFE[my_mesh->nElements()];
  for(int iElement=0; iElement<my_mesh->nElements(); iElement++) { 
    the_ds_elements[iElement].set(this, my_mesh->elementPtr(iElement));
  }
}

DirectSerendipity::~DirectSerendipity() {
  if (face_dofs) delete[] face_dofs;
//  if (cell_dofs) delete[] cell_dofs;
  if (the_ds_elements) delete[] the_ds_elements;
};

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
/*
  fout << "\ncell dofs:\n";
  for (int iElement=0; iElement<my_mesh->nElements(); iElement++) {
    fout << "\telement[" << iElement << "]\n";
    for (int iDoF=0; iDoF<nSingleCellDoFs(); iDoF++) {
      fout << "\t\tDoF[" << iDoF << "]: (" << cell_dofs[iElement*nSingleCellDoFs()+iDoF] << ")\n";
    }
  }
*/
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
/*
  // Cell nodes
  if (num_dofs_per_cell>0) {
    fout << "scatter3([ ";
    for(int i=0; i<num_dofs_per_cell*my_mesh->nElements(); i++) {
      fout << cell_dofs[i][0] << " ";
    }
    fout << "],[ ";
    for(int i=0; i<num_dofs_per_cell*my_mesh->nElements(); i++) {
      fout << cell_dofs[i][1] << " ";
    }
    fout << "],[ ";
    for(int i=0; i<num_dofs_per_cell*my_mesh->nElements(); i++) {
      fout << cell_dofs[i][2] << " ";
    }
    fout << "],'^r','filled')\n";
  }
*/

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

  // Element centers
  fout << "scatter3([ ";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << my_mesh->elementPtr(i)->center()[0] << " ";
  }
  fout << "],[ ";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << my_mesh->elementPtr(i)->center()[1] << " ";
  }
  fout << "],[ ";
  for(int i=0; i<my_mesh->nElements(); i++) {
    fout << my_mesh->elementPtr(i)->center()[2] << " ";
  }
  fout << "],'*g')\n";

  fout << "hold off;\n";
  return 0;
}