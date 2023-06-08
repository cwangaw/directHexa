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
#include "umfpack.h"

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
    std::vector<std::array<int,2>>* vec_ptr = param.dsSpace.indicesPtr();
    for (unsigned int iArray=0; iArray<(*vec_ptr).size(); iArray++) {
      std::cout << iArray << ":  (" <<  (*vec_ptr)[iArray][0] << "," << (*vec_ptr)[iArray][1] << ")" << std::endl;
    }
  }
    

  // SOLVE THE PDE ///////////////////////////////////////////////////////////
  
  monitor(0,"\nSolve the PDE\n");
  
  DirectSerendipityArray solution(&(param.dsSpace));

  std::vector<double> bc_vals; bc_vals.clear();

  for (int i = 0; i < param.dsSpace.nDoFs(); i++) {
    if (!param.dsSpace.isInterior(i)) {
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

  int nn = param.dsSpace.nInteriorDoFs();
  std::vector<std::array<int,2>>* vec_ptr = param.dsSpace.indicesPtr();
  std::vector<int32_t> Ap(nn+1); 
  std::vector<int32_t> Ai((*vec_ptr).size());
  std::vector<double> Ax((*vec_ptr).size(),0);
  std::vector<double> b((*vec_ptr).size(),0);
  std::vector<double> x(nn);

  Ap[0] = 0;
  for (unsigned int iArray=0; iArray<(*vec_ptr).size(); iArray++) {
    Ai[iArray] = (*vec_ptr)[iArray][0];
    if (iArray>0 && (*vec_ptr)[iArray][1]!=(*vec_ptr)[iArray-1][1]) {
      Ap[(*vec_ptr)[iArray][1]] = iArray;
    }
  }
  Ap[nn] = (*vec_ptr).size();

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

    for (int i=0; i<nn_loc; i++) {
      node_loc_to_gbl[i] = param.dsSpace.indexCorrection(fePtr->globalDoF(i));
      node_loc_to_bc[i] = param.dsSpace.bcCorrection(fePtr->globalDoF(i));
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
            // get the global index in Ax
            int global_index = 0;
            for (int possible_index=Ap[node_loc_to_gbl[iNode]]; possible_index<Ap[node_loc_to_gbl[iNode]+1]; possible_index++) {
              if (node_loc_to_gbl[jNode] == Ai[possible_index]) {
                global_index = possible_index;
                break;
              }
            }
            Ax[global_index] += mat_loc[iNode + nn_loc*jNode];
          }
        }
      }
    }

    for(int iNode=0; iNode<nn_loc; iNode++) {
      if (node_loc_to_gbl[iNode] != -1) {
        b[node_loc_to_gbl[iNode]] += rhs_loc[iNode];
      }
    }
  }

  monitor(1,"===============================================");
  monitor(1,"===Solution of linear system (SPARSE SOLVER)==="); ////////////////////////////////////////
  monitor(1,"===============================================");
  
  double *null = (double *) NULL ;
  void *Symbolic, *Numeric ;
  
  int error_sym = umfpack_di_symbolic(nn, nn, Ap.data(), Ai.data(), Ax.data(), &Symbolic, null, null);
  if (error_sym != 0) {
    std::cout << "ERROR: umfpack_di_symbolic with error code " << error_sym << std::endl;
  }
  int error_num =  umfpack_di_numeric(Ap.data(), Ai.data(), Ax.data(), Symbolic, &Numeric, null, null);
  if (error_num != 0) {
    std::cout << "ERROR: umfpack_di_numeric with error code " << error_num << std::endl;
  }
  umfpack_di_free_symbolic(&Symbolic);
  
  int error_sol = umfpack_di_solve(UMFPACK_A, Ap.data(), Ai.data(), Ax.data(), x.data(), b.data(), Numeric, null, null) ;
  if (error_sol != 0) {
    std::cout << "ERROR: umfpack_di_solve with error code " << error_sol << std::endl;
  }
  umfpack_di_free_numeric(&Numeric);

  std::ofstream sout("test/solution.txt");
  sout.precision(24);
  for(int i=0; i<nn; i++) {
    sout << x[i];
    if (i < nn - 1) sout << "\n";
  }

  for(int i=0; i<solution.size(); i++) {
    if(param.dsSpace.indexCorrection(i) == -1) {
      solution[i] = bc_vals[param.dsSpace.bcCorrection(i)];
    } else {
      solution[i] = x[param.dsSpace.indexCorrection(i)];
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
  }
  return 0;
}