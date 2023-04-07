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
  //monitor(0,"Polynomial Degree = ", parameterDataPtr()->dsSpace.degPolyn());

  ParameterData& param = *parameterDataPtr();

  // TEST QUADRATURE ///////////////////////////////////////////////////////

  if(false) {
    monitor(1,"Test Quadrature Rules");
    testQuadrature(&(param.mesh), param.supplement_type , param.refinement_level,1e-15);
  }

  return 0;
}