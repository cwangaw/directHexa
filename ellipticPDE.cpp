#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "ellipticPDE.h"
#include "fcns.h"
#include "Mesh/baseObjects.h"
#include "Mesh/hexaMesh.h"
#include "directSerendipity.h"
#include "parameterData.h"
#include "quadrature.h"
#include "Utilities/monitor.h"
#include "Utilities/debug.h"
#include <complex.h>
#include "lapacke.h"
#include <cblas.h>

using namespace directserendipity;
using namespace hexamesh;
using namespace quadrature;

int EllipticPDE::solve(Monitor& monitor) {
  monitor(0,"Polynomial Degree = ", parameterDataPtr()->dsSpace.degPolyn());

  ParameterData& param = *parameterDataPtr();

  // TEST QUADRATURE ///////////////////////////////////////////////////////

  if(false) {
    monitor(1,"Test Quadrature Rules");
    testQuadrature(&(param.mesh), param.supplement_type , param.refinement_level,1e-14);
  }

  // TEST BASIS FUNCTION ///////////////////////////////////////////////////////

  if (false) {
    DirectSerendipityArray testingarray(&(param.dsSpace));
    for(int i=0; i<param.mesh.nVertices(); i++) {
      testingarray[i] = 0;  
    }

    testingarray[13] = 1;

    std::string fileName(param.directory_name);
    fileName += "test_mesh";
    std::string fileNameGrad(param.directory_name);
    fileNameGrad += "test_grad_mesh";
    testingarray.write_tecplot_mesh(fileName,fileNameGrad,
       param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y,param.output_mesh_numPts_DS_z);
  }

  // DEBUG ///////////////////////////////////////////////////////////

  if (false) {
    DirectSerendipityArray testingarray(&(param.dsSpace));
    for(int i=0; i<param.mesh.nVertices(); i++) {
      testingarray[i] = 0;  
    }
  }
    

  // SOLVE THE PDE ///////////////////////////////////////////////////////////
  
  monitor(0,"\nSolve the PDE\n");
  
  DirectSerendipityArray solution(&(param.dsSpace));

  // Correct for BCSs: If node i is an interior node, i-th element of index_correction 
  // would give its index in our vector without BC nodes. If node i is a BC node, the  
  // i-th element would give -1.
  std::vector<int> index_correction(param.dsSpace.nDoFs());
  int nn = 0;

  // Correct for BCSs: If node i is a BC node, i-th element of index_correction 
  // would give its index in our vector without interior nodes. If node i is an interior node, the  
  // i-th element would give -1.
  std::vector<int> bc_correction(param.dsSpace.nDoFs());
  int mm = 0;

  std::vector<double> bc_vals; bc_vals.clear();

  for (int i = 0; i < param.dsSpace.nDoFs(); i++) {
    if (param.dsSpace.isInterior(i)) {
      bc_correction[i] = -1;
      index_correction[i] = nn;
      nn++;
    } else {
      index_correction[i] = -1;
      bc_correction[i] = mm;
      mm++;
      // evaluate the bc point
      if (i<param.dsSpace.nVertexDoFs()) {
        bc_vals.push_back(bcVal(param.mesh.vertexPtr(i)->val(0), param.mesh.vertexPtr(i)->val(1), param.mesh.vertexPtr(i)->val(2)));
      } else if (i<param.dsSpace.nVertexDoFs()+param.dsSpace.nEdgeDoFs()) {
        int iEdge = (i-param.dsSpace.nVertexDoFs()) / (param.dsSpace.degPolyn()-1);
        int nPt = (i-param.dsSpace.nVertexDoFs()) % (param.dsSpace.degPolyn()-1);
        Point pt(*param.mesh.edgeVertexPtr(0,iEdge) + (double)(nPt+1)/param.dsSpace.degPolyn() * (*param.mesh.edgeVertexPtr(1,iEdge)-*param.mesh.edgeVertexPtr(0,iEdge)));
        bc_vals.push_back(bcVal(pt[0],pt[1],pt[2]));
      } else {
        int iFace = (i-param.dsSpace.nVertexDoFs()-param.dsSpace.nEdgeDoFs()) / param.dsSpace.nSingleFaceDoFs();
        int nPt = (i-param.dsSpace.nVertexDoFs()-param.dsSpace.nEdgeDoFs()) % param.dsSpace.nSingleFaceDoFs();
        bc_vals.push_back(bcVal(param.dsSpace.faceDoFPtr(iFace,nPt)->val(0), param.dsSpace.faceDoFPtr(iFace,nPt)->val(1), param.dsSpace.faceDoFPtr(iFace,nPt)->val(2)));
      }
    }
  }

  if(param.dsSpace.degPolyn() >= 2) param.dsSpace.bcModification(bc_vals.data());

  std::vector<double> mat_vector(nn*nn,0);
  double* mat = mat_vector.data();
  std::vector<double> rhs_vector(nn,0);
  double* rhs = rhs_vector.data();

  // quadrature points
  quadrature::Quadrature quadRule(10);

  monitor(1,"Matrix and RHS Assembly"); ////////////////////////////////////////

  for(int iElement=0; iElement<param.mesh.nElements(); iElement++) {
    DirectSerendipityFE* fePtr = param.dsSpace.finiteElementPtr(iElement);
    quadRule.setElement(param.supplement_type,param.refinement_level,fePtr->elementPtr());
    fePtr->initBasis(quadRule.pts(), quadRule.num());
    // Local matrix and rhs
    int nn_loc = fePtr->nDoFs();
    std::vector<double> mat_loc(nn_loc*nn_loc,0), rhs_loc(nn_loc,0);

    // Determine local to global map
    int node_loc_to_gbl[nn_loc];

    // Determine local to global bc map
    int node_loc_to_bc[nn_loc];

    for(int i=0; i<8; i++) {
      node_loc_to_gbl[i] = index_correction[fePtr->elementPtr()->vertexGlobal(i)];
      node_loc_to_bc[i] = bc_correction[fePtr->elementPtr()->vertexGlobal(i)];
    }

    for (int iEdge=0; iEdge<12; iEdge++) {
      int iEdge_global = fePtr->elementPtr()->edgeGlobal(iEdge);
      for (int jDoF=0; jDoF<param.dsSpace.degPolyn()-1; jDoF++) {
        node_loc_to_gbl[8 + iEdge*(param.dsSpace.degPolyn()-1) + jDoF] = index_correction[param.dsSpace.nVertexDoFs() + iEdge_global*(param.dsSpace.degPolyn()-1) + jDoF];
        node_loc_to_bc[8 + iEdge*(param.dsSpace.degPolyn()-1) + jDoF] = bc_correction[param.dsSpace.nVertexDoFs() + iEdge_global*(param.dsSpace.degPolyn()-1) + jDoF];
      }
    }

    if (fePtr->nFaceDoFs()>0) {
      for (int iFace=0; iFace<6; iFace++) {
        int iFace_global = fePtr->elementPtr()->faceGlobal(iFace);
        for (int jDoF=0; jDoF<fePtr->nFaceDoFs()/6; jDoF++) {
          node_loc_to_gbl[8 + 12*(param.dsSpace.degPolyn()-1) + iFace*fePtr->nFaceDoFs()/6 + jDoF] = index_correction[param.dsSpace.nVertexDoFs() + param.dsSpace.nEdgeDoFs() + iFace_global*fePtr->nFaceDoFs()/6 + jDoF];
          node_loc_to_bc[8 + 12*(param.dsSpace.degPolyn()-1) + iFace*fePtr->nFaceDoFs()/6 + jDoF] = bc_correction[param.dsSpace.nVertexDoFs() + param.dsSpace.nEdgeDoFs() + iFace_global*fePtr->nFaceDoFs()/6 + jDoF];
        }
      }
    }

    if (fePtr->nCellDoFs()>0) {
      int iCell_global = fePtr->elementPtr()->meshIndex();
      for (int jDoF=0; jDoF<fePtr->nCellDoFs(); jDoF++) {
        node_loc_to_gbl[8 + 12*(param.dsSpace.degPolyn()-1) + fePtr->nFaceDoFs() + jDoF] = index_correction[param.dsSpace.nVertexDoFs() + param.dsSpace.nEdgeDoFs() + param.dsSpace.nFaceDoFs() + iCell_global*fePtr->nCellDoFs() + jDoF];
        node_loc_to_gbl[8 + 12*(param.dsSpace.degPolyn()-1) + fePtr->nFaceDoFs() + jDoF] = -1;
      }
    }

    // Matrix and rhs assembly over elements
    for(int iPt=0; iPt<quadRule.num(); iPt++) {
      double x = quadRule.pt(iPt)[0];
      double y = quadRule.pt(iPt)[1];
      double z = quadRule.pt(iPt)[2];

      double valA = coefA(x,y,z);
      Tensor1 valB, valC;
      Tensor2 valD;
      coefB(x,y,z,valB); coefC(x,y,z,valC); coefD(x,y,z,valD);

      // Local interactions	      
      for(int jNode=0; jNode<nn_loc; jNode++) {
        if (node_loc_to_gbl[jNode] == -1) continue;
	      double valj = fePtr->basis(jNode, iPt); 
        Tensor1 gradValj = fePtr->basisGrad(jNode, iPt); 

        for(int iNode=0; iNode<nn_loc; iNode++) {
          //In case "shape functions" are not delta_{i,j} for BC nodes one day
          //if on BC, we use nodal basis function, otherwise we use shape functions

          double vali = fePtr->basis(iNode, iPt);
          Tensor1 gradVali = fePtr->basisGrad(iNode, iPt);

          if (node_loc_to_gbl[iNode] != -1) {
            // +(a N_i,N_j);
	          mat_loc[iNode + nn_loc*jNode] += valA*vali*valj*quadRule.wt(iPt);
            // +(c dot grad(N_i), N_j)
            mat_loc[iNode + nn_loc*jNode] += (valC*gradVali)*valj*quadRule.wt(iPt);
            // +(b N_i, grad(N_j))
            mat_loc[iNode + nn_loc*jNode] += vali*(valB*gradValj)*quadRule.wt(iPt);
            // +(D grad(N_i), grad(N_j))
            mat_loc[iNode + nn_loc*jNode] += ((valD*gradVali)*gradValj)*quadRule.wt(iPt);
          }
          
          //rhs

          if (node_loc_to_gbl[iNode] == -1) {
	          double bcVali = bc_vals[node_loc_to_bc[iNode]];
            // -(a g(x_i)*N^{BC}_i, N_j)
            rhs_loc[jNode] -= valA*vali*valj*bcVali*quadRule.wt(iPt); 
            // -(c dot g(x_i)*grad(N^{BC}_i), N_j)
            rhs_loc[jNode] -= (valC*gradVali)*valj*bcVali*quadRule.wt(iPt);
            // -(b g(x_i)*N^{BC}_i, grad(N_j))
            rhs_loc[jNode] -= vali*(valB*gradValj)*bcVali*quadRule.wt(iPt);
            // -(D g(x_i)*grad(N^{BC}_i), grad(N_j))
            rhs_loc[jNode] -= ((valD*gradVali)*gradValj)*bcVali*quadRule.wt(iPt);
          }
	      }
        // +(f, N_j)
	      rhs_loc[jNode] += sourceVal(quadRule.pt(iPt)[0],quadRule.pt(iPt)[1],quadRule.pt(iPt)[2])*valj*quadRule.wt(iPt);
      }
    }

    // Map local matrix and rhs to global
    
    for(int jNode=0; jNode<nn_loc; jNode++) {
      if (node_loc_to_gbl[jNode] != -1) {
        for(int iNode=0; iNode<nn_loc; iNode++) {
          if (node_loc_to_gbl[iNode] != -1) {
            mat[node_loc_to_gbl[iNode] + nn*node_loc_to_gbl[jNode]]
              += mat_loc[iNode + nn_loc*jNode];
          }
        }
      }
    }

    for(int iNode=0; iNode<nn_loc; iNode++) {
      if (node_loc_to_gbl[iNode] != -1) {
        rhs[node_loc_to_gbl[iNode]] += rhs_loc[iNode];
      }
    }
  }

  std::ofstream fout("test/matrix.txt");
  fout.precision(24);
  for(int j=0; j<nn; j++) {
    for(int i=0; i<nn; i++) {
      fout << mat[i + nn*j] << "\t";
    }
    if (j < nn - 1) fout << "\n";
  }


  std::ofstream rout("test/rhs.txt");
  rout.precision(24);
  for(int i=0; i<nn; i++) {
    rout << rhs[i];
    if (i < nn - 1) rout << "\n";
  }



  monitor(1,"===============================================");
  monitor(1,"===Solution of linear system (DIRECT SOLVER)==="); ////////////////////////////////////////
  monitor(1,"===============================================");
  
  //Solve the matrix, result would be stored in rhs
  lapack_int* ipiv; char norm = 'I'; 
  ipiv = (lapack_int*)malloc(nn * sizeof(lapack_int));
  double anorm = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, nn, nn, mat, nn);
  int ierr = LAPACKE_dgesv(LAPACK_ROW_MAJOR, nn, 1, mat, nn, ipiv, rhs, 1); //mat updated to be LU
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  double rcond = 0;
  ierr = LAPACKE_dgecon(LAPACK_ROW_MAJOR, norm, nn, mat, nn, anorm, &rcond);
  if(ierr) { // ?? what should we do ???
    std::cerr << "ERROR: Lapack failed with code " << ierr << std::endl; 
  }
  rcond = 1/rcond;
  free(ipiv);

  std::ofstream sout("test/solution.txt");
  sout.precision(24);
  for(int i=0; i<nn; i++) {
    sout << rhs[i];
    if (i < nn - 1) sout << "\n";
  }

  //Calculate inf condition number
 
  std::cout << "\tNorm Format:\t" << norm << std::endl;
  std::cout << "\tNorm of mat:\t" << anorm << std::endl;
  std::cout << "\tCond number:\t" << rcond << std::endl;

  for(int i=0; i<solution.size(); i++) {
    if(index_correction[i] == -1) {
      solution[i] = bc_vals[bc_correction[i]];
    } else {
      solution[i] = rhs[index_correction[i]];
    }
  }

  if(param.output_soln_DS_format > 0) {
    monitor(1,"Write Solution"); //////////////////////////////////////////////////

    switch(param.output_soln_DS_format) {
    case 1: {
      std::string fileName(param.directory_name);
      fileName += "solution_raw";
      solution.write_raw(fileName);
      break;
    }
    case 2: {
      std::string fileName(param.directory_name);
      fileName += "solution_mesh";
      std::string fileNameGrad(param.directory_name);
      fileNameGrad += "solution_grad_mesh";
      solution.write_tecplot_mesh(fileName,fileNameGrad,
				 param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y,param.output_mesh_numPts_DS_z);
      break;
    }
    }
  }

  if(trueSolnKnown()) {
    monitor(0,"\nError estimate\n"); ///////////////////////////////////////////////

    double l2Error = 0, l2GradError = 0, l2Norm = 0, l2GradNorm = 0;
    solution.l2normError(l2Error, l2GradError, l2Norm, l2GradNorm, param.refinement_level, trueSoln, trueGradSoln);

    std::cout << "  L_2 Error:      " << l2Error << std::endl;
    std::cout << "  L_2 Grad Error: " << l2GradError << std::endl;
    std::cout << std::endl;
    std::cout << "  Relative L_2 Error:      " << l2Error/l2Norm << std::endl;
    std::cout << "  Relative L_2 Grad Error: " << l2GradError/l2GradNorm << std::endl;
    std::cout << std::endl;

    if(param.output_soln_DS_format > 0) {

      monitor(1,"Write Interpolation Solution, only work for polynomial degree r < 6\n"); ////////////////////////////////////////////

      DirectSerendipityArray u(&(param.dsSpace));

      for(int i=0; i<u.size(); i++) {
        if (i<param.dsSpace.nVertexDoFs()) {
          Vertex* v = param.mesh.vertexPtr(i);
          u[i] = trueSoln(v->val(0),v->val(1),v->val(2));
        } else if (i<param.dsSpace.nVertexDoFs()+param.dsSpace.nEdgeDoFs()) {
          Point* pt = param.dsSpace.edgeNodePtr(i-param.dsSpace.nVertexDoFs());
          u[i] = trueSoln(pt->val(0),pt->val(1),pt->val(2));
        } else {
          Point* pt = param.dsSpace.faceDoFPtr(i-param.dsSpace.nVertexDoFs()-param.dsSpace.nEdgeDoFs());
          u[i] = trueSoln(pt->val(0),pt->val(1),pt->val(2));
        }
      }

      param.dsSpace.nodeModification(u.theArray());

      double l2InterpError = 0, l2InterpGradError = 0, l2InterpNorm = 0, l2InterpGradNorm = 0;
      u.l2normError(l2InterpError, l2InterpGradError, l2InterpNorm, l2InterpGradNorm, param.refinement_level, trueSoln, trueGradSoln);

      std::cout << "  L_2 Interpolation Error:      " << l2InterpError << std::endl;
      std::cout << "  L_2 Interpolation Grad Error: " << l2InterpGradError << std::endl;
      std::cout << std::endl;
      std::cout << "  Relative L_2 Interpolation Error:      " << l2InterpError/l2InterpNorm << std::endl;
      std::cout << "  Relative L_2 Interpolation Grad Error: " << l2InterpGradError/l2InterpGradNorm << std::endl;
      std::cout << std::endl;

      switch(param.output_soln_DS_format) {
      case 1: {
        std::string fileName(param.directory_name);
        fileName += "interpolation_raw";
        u.write_raw(fileName);
        break;
      }
      case 2: {
        std::string fileName(param.directory_name);
        fileName += "interpolation_mesh";
        std::string fileNameGrad(param.directory_name);
        fileNameGrad += "interpolation_grad_mesh";
        u.write_tecplot_mesh(fileName,fileNameGrad,
        param.output_mesh_numPts_DS_x,param.output_mesh_numPts_DS_y,param.output_mesh_numPts_DS_z);
        break;
      }
      }
    }
  }
  return 0;
}