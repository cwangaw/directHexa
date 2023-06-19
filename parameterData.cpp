#include <fstream>
#include <algorithm>

#include "debug.h"
#include "parameterData.h"
#include "Reader/reader.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// PARAMETERDATA

ParameterData::ParameterData()
  : infile_name("infile"),
    echofile_name("echo"),
    echo(0) {}

ParameterData::ParameterData(const std::string& infileName)
  : infile_name(infileName),
    echofile_name("echo"),
    echo(0) {}

ParameterData::ParameterData(const std::string& infileName,
				  const std::string& echofileName, bool echo0)
  : infile_name(infileName),
    echofile_name(echofileName),
    echo(echo0) {}

// CASES =========================================================================

void ParameterData::print_cases(ostream& fout) {
  int n;

  // Soln DS output
  fout<< "Soln output formats (" << n_case_soln_DS_output << " total):\n";
  n=-1;
  fout<< "  " << case_soln_DS_output_omit
      << " = do not compute DS soln\n";
  n++;
  n=0;
  fout<< "  " << case_soln_DS_output_none
      << " = no soln output\n";
  n++;
  fout<< "  " << case_soln_DS_output_raw
      << " = raw text file\n";
  n++;
  fout<< "  " << case_soln_DS_output_matlab
      << " = matlab file\n";
  n++;
  if(n != n_case_soln_DS_output) {
    fout<< "  " << "other cases exist" << endl;
  }

  // Mesh output
  fout<< "Mesh output formats (" << n_case_mesh_output << " total):\n";
  n=0;
  fout<< "  " << case_mesh_output_none
      << " = no mesh output\n";
  n++;
  fout<< "  " << case_mesh_output_raw
      << " = raw text file\n";
  n++;
  fout<< "  " << case_mesh_output_matlab
      << " = matlab file\n";
  n++;
  if(n != n_case_mesh_output) {
    fout<< "  " << "other cases exist" << endl;
  }

  // DSSpace output
  fout<< "DS space output formats (" << n_case_dsSpace_output << " total):\n";
  n=0;
  fout<< "  " << case_dsSpace_output_none
      << " = no DS space output\n";
  n++;
  fout<< "  " << case_dsSpace_output_raw
      << " = raw text file\n";
  n++;
  fout<< "  " << case_dsSpace_output_matlab
      << " = matlab file\n";
  n++;
  if(n != n_case_dsSpace_output) {
    fout<< "  " << "other cases exist" << endl;
  }
}

// ===================================================================================

#define ERRCHK(s) {error=s; if(error){processReaderError(error);return error;}}
#define ERRRET(s) {error=s; if(error) return error;}

// MESH AND FINITE ELEMENTS ==========================================================

// Output mesh
int ParameterData::writeMesh() {
  if(output_mesh_format == case_mesh_output_none) return 0;

  string fileName = directory_name;
  fileName += "meshOut";

  switch(output_mesh_format) {
  case case_mesh_output_raw: {
    fileName += ".txt";
    if(mesh.write_raw(fileName)) return ERR_FILESYSTEM;
    break;
  }
  case case_mesh_output_matlab: {
    if(mesh.write_matlab(fileName)) return ERR_FILESYSTEM;
    break;
  }
  default:
    return ERR_UNSUPPORTED_OPTION;
  }
  return 0;
}

// Output DS space
int ParameterData::writeDSSpace() {
  if(output_dsSpace_format == case_dsSpace_output_none) return 0;

  string fileName = directory_name;
  fileName += "dsSpaceOut";

  switch(output_dsSpace_format) {
  case case_dsSpace_output_raw: {
    fileName += ".txt";
    if(dsSpace.write_raw(fileName)) return ERR_FILESYSTEM;
    break;
  }
  case case_dsSpace_output_matlab: {
    if(dsSpace.write_matlab(fileName)) return ERR_FILESYSTEM;
    break;
  }  
  default:
    return ERR_UNSUPPORTED_OPTION;
  }

  return 0;
}

// READ ==========================================================================

int ParameterData::read() {
  int error = 0;

  ERRCHK(openReader(infile_name, echofile_name, 121, 1, echo));

  // OUTPUT DIRECTORY

  ERRCHK(readString(directory_name));

  // MESH PARAMETERS

  ERRCHK(readString(mesh_type));
  char meshTypeC;
  switch(mesh_type[0]) {
  case 'r': case 'R': meshTypeC = 'r'; break;
  case 'o': case 'O': meshTypeC = 'o'; break;
  case 't': case 'T': meshTypeC = 't'; break;
  default: return processReaderError(ERR_UNSUPPORTED_OPTION);
  }	       

  switch(meshTypeC) {
  case 'r':
  case 'o':
  case 't': {
    double xMin, xMax;
    ERRCHK(readScalar(xMin));
    ERRCHK(readScalar(xMax));
    if(xMin >= xMax) return processReaderError(ERR_BAD_DATA);
    
    double yMin, yMax;
    ERRCHK(readScalar(yMin));
    ERRCHK(readScalar(yMax));
    if(yMin >= yMax) return processReaderError(ERR_BAD_DATA);

    double zMin, zMax;
    ERRCHK(readScalar(zMin));
    ERRCHK(readScalar(zMax));
    if(zMin >= zMax) return processReaderError(ERR_BAD_DATA);

    int nx, ny, nz;
    ERRCHK(readScalar(nx));
    ERRCHK(readScalar(ny));
    ERRCHK(readScalar(nz));
    if(nx < 1) nx = 1;
    if(ny < 1) ny = 1;
    if(nz < 1) nz = 1;

    ERRCHK(readScalar(distortion_factor));
    if(distortion_factor >= 0.5) return processReaderError(ERR_BAD_DATA);

    mesh.createMesh(meshTypeC, nx,ny,nz, xMin,xMax, yMin,yMax, zMin,zMax, distortion_factor);
    break;
  }
  default:
    return processReaderError(ERR_UNSUPPORTED_OPTION);
  }

  // FINITE ELEMENTS

  ERRCHK(readScalar(polynomial_degree));
  if (polynomial_degree < 0) processReaderError(ERR_BAD_DATA);

  ERRCHK(readScalar(supplement_type));
  if (supplement_type < 0  || supplement_type > 2) processReaderError(ERR_BAD_DATA);

  ERRCHK(readScalar(refinement_level));
  if (refinement_level < 0) processReaderError(ERR_BAD_DATA);
  dsSpace.set(polynomial_degree, supplement_type, &mesh);

  // OUTPUT PARAMETERS

  // DS
  ERRCHK(readScalar(output_soln_DS_format));
  if(output_soln_DS_format < -1 || output_soln_DS_format >= n_case_soln_DS_output)
    return processReaderError(ERR_OUT_OF_RANGE);

  ERRCHK(readScalar(output_mesh_numPts_DS_x));
  if(output_mesh_numPts_DS_x < 2) output_mesh_numPts_DS_x = 2;
  ERRCHK(readScalar(output_mesh_numPts_DS_y));
  if(output_mesh_numPts_DS_y < 2) output_mesh_numPts_DS_y = 2;
  ERRCHK(readScalar(output_mesh_numPts_DS_z));
  if(output_mesh_numPts_DS_z < 2) output_mesh_numPts_DS_z = 2;


  // Mesh and Spaces
  ERRCHK(readScalar(output_mesh_format));
  if(output_mesh_format < 0 || output_mesh_format >= n_case_mesh_output)
    return processReaderError(ERR_OUT_OF_RANGE);
  ERRCHK(writeMesh());

  ERRCHK(readScalar(output_dsSpace_format));
  if(output_dsSpace_format < 0 || output_dsSpace_format >= n_case_dsSpace_output)
    return processReaderError(ERR_OUT_OF_RANGE);
  ERRCHK(writeDSSpace());

  ERRCHK(readScalar(monitor_to_level));

  return error;
}
