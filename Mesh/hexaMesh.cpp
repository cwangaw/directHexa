#include <cmath>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <assert.h>
#include <random>
#include <string>
#include "debug.h"
#include "hexaMesh.h"
using namespace hexamesh;

////////////////////////////////////////////////////////////////////////////////
// Class Vertex

void Vertex::set_vertex(double px, double py, double pz, int ix, int iy, int iz,  int i, HexaMesh* myMesh) {
  the_point[0] = px; 
  the_point[1] = py; 
  the_point[2] = pz; 
  my_mesh_pos[0] = ix; my_mesh_pos[1] = iy; my_mesh_pos[2] = iz; 
  my_mesh_index = i;
  my_mesh = myMesh; 
}

void Vertex::write_raw(std::ofstream& fout) const {
  fout << "      VERTEX (" << val(0) << "," <<val(1) << "," << val(2) <<")\n";
  fout << "      my_mesh       = " << my_mesh << "\n";
  fout << "      my_mesh_pos = (" << my_mesh_pos[0] << "," << my_mesh_pos[1] << "," << my_mesh_pos[2] <<")\n";
  fout << "      my_mesh_index = " << my_mesh_index <<"\n";
}

int Vertex::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// Class Element
static const double element_eps = 1e-8;

void Element::set_element(int ix, int iy, int iz, int i, HexaMesh* myMesh) {
  my_mesh_pos[0] = ix;
  my_mesh_pos[1] = iy;
  my_mesh_pos[2] = iz;
  my_mesh_index = i;
  my_mesh = myMesh;
  good_element = true;

  if(vertex_global_index) delete[] vertex_global_index;
  vertex_global_index = new int[8];
  if(!vertex_global_index) return;

  for (int ix = 0; ix < 2; ix++) {
    for (int iy = 0; iy < 2; iy++) {
      for (int iz = 0; iz < 2; iz++) {
        vertex_global_index[ix*4+iy*2+iz] = my_mesh -> vertexIndex(meshPos(0)+ix, meshPos(1)+iy, meshPos(2)+iz);
      }
    }
  }

  for (int iFace = 0; iFace < 6; iFace++) {
    // If v11 is not in the space spanned by {v01-v00, v10-v00}, we set good_element to false
    if ( Tensor1(*faceVertexPtr(1,1,iFace)-*faceVertexPtr(0,0,iFace))*normal(iFace)  > element_eps ) good_element = false;
  }

  subtetra_d.clear();
  if ((my_mesh_pos[0]+my_mesh_pos[1]+my_mesh_pos[2]) % 2 == 0) {
    // boundary tetra with a vertex at v_{-1,-2,-3}
    subtetra_d.push_back({vertexIndex(-1,-1,-1), vertexIndex(1,-1,-1), vertexIndex(-1,1,-1), vertexIndex(-1,-1,1)});
    // boundary tetra with a vertex at v_{-1,2,3}
    subtetra_d.push_back({vertexIndex(-1,1,1), vertexIndex(1,1,1), vertexIndex(-1,-1,1), vertexIndex(-1,1,-1)});
    // boundary tetra with a vertex at v_{1,-2,3}
    subtetra_d.push_back({vertexIndex(1,-1,1), vertexIndex(-1,-1,1), vertexIndex(1,1,1), vertexIndex(1,-1,-1)});
    // boundary tetra with a vertex at v_{1,2,-3}
    subtetra_d.push_back({vertexIndex(1,1,-1), vertexIndex(-1,1,-1), vertexIndex(1,-1,-1), vertexIndex(1,1,1)});
    // center tetra
    subtetra_d.push_back({vertexIndex(1,1,1), vertexIndex(1,-1,-1), vertexIndex(-1,1,-1),  vertexIndex(-1,-1,1)});
  } else {
    // boundary tetra with a vertex at v_{1,2,3}
    subtetra_d.push_back({vertexIndex(1,1,1), vertexIndex(-1,1,1), vertexIndex(1,-1,1), vertexIndex(1,1,-1)});
    // boundary tetra with a vertex at v_{1,-2,-3}
    subtetra_d.push_back({vertexIndex(1,-1,-1), vertexIndex(-1,-1,-1), vertexIndex(1,1,-1), vertexIndex(1,-1,1)});
    // boundary tetra with a vertex at v_{-1,2,-3}
    subtetra_d.push_back({vertexIndex(-1,1,-1), vertexIndex(1,1,-1), vertexIndex(-1,-1,-1), vertexIndex(-1,1,1)});
    // boundary tetra with a vertex at v_{-1,-2,3}
    subtetra_d.push_back({vertexIndex(-1,-1,1), vertexIndex(1,-1,1), vertexIndex(-1,1,1), vertexIndex(-1,-1,-1)});
    // center tetra
    subtetra_d.push_back({vertexIndex(-1,-1,-1), vertexIndex(-1,1,1), vertexIndex(1,-1,1), vertexIndex(1,1,-1)});
  }

  subtetra_m.clear();
  // T_{1,-2}
  subtetra_m.push_back({vertexIndex(1,1,1), vertexIndex(-1,-1,-1), vertexIndex(1,-1,-1), vertexIndex(1,-1,1)});
  // T_{1,-3}
  subtetra_m.push_back({vertexIndex(1,1,1), vertexIndex(-1,-1,-1), vertexIndex(1,-1,-1), vertexIndex(1,1,-1)});
  // T_{2,-1}
  subtetra_m.push_back({vertexIndex(1,1,1), vertexIndex(-1,-1,-1), vertexIndex(-1,1,-1), vertexIndex(-1,1,1)});
  // T_{2,-3}
  subtetra_m.push_back({vertexIndex(1,1,1), vertexIndex(-1,-1,-1), vertexIndex(-1,1,-1), vertexIndex(1,1,-1)});        
  // T_{3,-1}
  subtetra_m.push_back({vertexIndex(1,1,1), vertexIndex(-1,-1,-1), vertexIndex(-1,-1,1), vertexIndex(-1,1,1)});
  // T_{3,-2}
  subtetra_m.push_back({vertexIndex(1,1,1), vertexIndex(-1,-1,-1), vertexIndex(-1,-1,1), vertexIndex(1,-1,1)});     
}

Element::~Element() {
  if(vertex_global_index) delete[] vertex_global_index;
}

void Element::vertexPos(int i, int& sgnx, int& sgny, int& sgnz) {
  i = i % 8;
  sgnx = (i/4 < 1)? -1:1;
  sgny = (i%4 < 2)? -1:1;
  sgnz = (i%2 < 1)? -1:1;
  return;
}

int Element::vertexIndex(int sgnx, int sgny, int sgnz) {
  return std::max(sgnx,0) * 4 + std::max(sgny,0) * 2 + std::max(sgnz,0);
}

int Element::vertexGlobal(int sgnx, int sgny, int sgnz) {
  return vertex_global_index[ vertexIndex(sgnx, sgny, sgnz) ];
}

int Element::vertexGlobal(int i) {
  int sgnx, sgny, sgnz;
  vertexPos(i, sgnx, sgny, sgnz);
  return vertexGlobal(sgnx, sgny, sgnz);
}

Vertex* Element::vertexPtr(int sgnx, int sgny, int sgnz) {
  return my_mesh->vertexPtr(vertexGlobal(sgnx, sgny, sgnz));
}

Vertex* Element::vertexPtr(int i) {
  return my_mesh->vertexPtr(vertexGlobal(i));
}

void Element::edgePos(int i, int& sgnx, int& sgny, int& sgnz) {
  i = i % 12;
  sgnx = (i%4 < 2)? -1:1;
  sgnz = (i%2 < 1)? -1:1;

  if (i < 4) {
    sgnz = 0;
    sgny = (i%2 < 1)? -1:1;
  } else if (i > 7) {
    sgnx = 0;
    sgny = (i%4 < 2)? -1:1;
  } else {
    sgny = 0;
  }

  return;
}

int Element::edgeIndex(int sgnx, int sgny, int sgnz) {
  if (sgnx == 0) {
    return 8 + std::max(sgny,0)*2 + std::max(sgnz,0);
  } else if ( sgny == 0 ) {
    return 4 + std::max(sgnx,0)*2 + std::max(sgnz,0);
  } else if ( sgnz == 0 ) {
    return std::max(sgnx,0)*2 + std::max(sgny,0);
  } else {
    std::cout << "ERROR: Problematic element edge position index";
    return -1;
  }
}

Vertex* Element::edgeVertexPtr(bool pm, int sgnx, int sgny, int sgnz) {
  if (sgnx == 0) {
    return vertexPtr(2*int(pm)-1, sgny, sgnz);
  } else if (sgny == 0) {
    return vertexPtr(sgnx, 2*int(pm)-1, sgnz);
  } else if (sgnz == 0) {
    return vertexPtr(sgnx, sgny, 2*int(pm)-1);
  } else {
    std::cout << "ERROR: Problematic element edge position index";
    return nullptr;
  }
}

Vertex* Element::edgeVertexPtr(bool pm, int i) {
  int sgnx, sgny, sgnz;
  edgePos(i, sgnx, sgny, sgnz);
  return edgeVertexPtr(pm, sgnx, sgny, sgnz);
}

int Element::edgeGlobal(int sgnx, int sgny, int sgnz) {
  Vertex* v0 = edgeVertexPtr(0, sgnx, sgny, sgnz);
  if (sgnx == 0) {
    int base = (mesh()->xElements()+1) * (mesh()->yElements()+1) * mesh()->zElements() + (mesh()->xElements()+1) * mesh()->yElements() * (mesh()->zElements()+1);
    return base + v0->meshPos(1) * mesh()->xElements() * (mesh()->zElements()+1) +  v0->meshPos(2) * mesh()->xElements() + v0->meshPos(0);
  } else if (sgny == 0) {
    int base = (mesh()->xElements()+1) * (mesh()->yElements()+1) * mesh()->zElements();
    return base + v0->meshPos(0) * mesh()->yElements() * (mesh()->zElements()+1) + v0->meshPos(2) * mesh()->yElements() + v0->meshPos(1);
  } else if (sgnz == 0) {
    return v0->meshPos(0) * mesh()->zElements() * (mesh()->yElements()+1) + v0->meshPos(1) * mesh()->zElements() +  v0->meshPos(2);
  } else {
    std::cout << "ERROR: Problematic element edge position index";
    return -1;
  }
}

void Element::edgeFace(int nEdge, int& f0, int& f1) {
  int sgnx, sgny, sgnz;
  edgePos(nEdge, sgnx, sgny, sgnz);
  if (sgnx==0) {
    f0 = faceIndex(0,sgny,0);
    f1 = faceIndex(0,0,sgnz);
  } else if (sgny==0) {
    f0 = faceIndex(sgnx,0,0);
    f1 = faceIndex(0,0,sgnz);
  } else {
    f0 = faceIndex(sgnx,0,0);
    f1 = faceIndex(0,sgny,0);        
  }
  return;
}

int Element::edgeGlobal(int i) {
  int sgnx, sgny, sgnz;
  edgePos(i, sgnx, sgny, sgnz);
  return edgeGlobal(sgnx, sgny, sgnz);
}

void Element::facePos(int i, int& sgnx, int& sgny, int& sgnz) {
  sgnx = 0; sgny = 0; sgnz = 0;
  i = i % 6;
  int nz_entry = 2*(i%2)-1; //nonzero entry

  if (i < 2) {
    sgnx = nz_entry;
  } else if (i < 4) {
    sgny = nz_entry;
  } else {
    sgnz = nz_entry;
  }

  return;
}

int Element::faceIndex(int sgnx, int sgny, int sgnz) {
  if (sgny == 0 && sgnz == 0) {
    return (sgnx+1)/2;
  } else if (sgnx == 0 && sgnz == 0) {
    return 2 + (sgny+1)/2;
  } else if (sgnx == 0 && sgny == 0) {
    return 4 + (sgnz+1)/2;
  } else {
    std::cout << "ERROR: Problematic element face position index";
    return -1;
  }
}

Vertex* Element::faceVertexPtr(bool pm0, bool pm1, int sgnx, int sgny, int sgnz) {
  if (sgny == 0 && sgnz == 0) {
    return vertexPtr(sgnx, 2*int(pm0)-1, 2*int(pm1)-1);
  } else if (sgnx == 0 && sgnz == 0) {
    return vertexPtr(2*int(pm0)-1, sgny, 2*int(pm1)-1);
  } else if (sgnx == 0 && sgny == 0) {
    return vertexPtr(2*int(pm0)-1, 2*int(pm1)-1, sgnz);
  } else {
    std::cout << "ERROR: Problematic element face position index";
    return nullptr;
  }  
}

Vertex* Element::faceVertexPtr(bool pm0, bool pm1, int i) {
  int sgnx, sgny, sgnz;
  facePos(i, sgnx, sgny, sgnz);
  return faceVertexPtr(pm0, pm1, sgnx, sgny, sgnz);
}

int Element::faceGlobal(int sgnx, int sgny, int sgnz) {
  Vertex* v0 = faceVertexPtr(0, 0, sgnx, sgny, sgnz);
  if (sgny == 0 && sgnz == 0) {
    return v0->meshPos(0) * mesh()->yElements() * mesh()->zElements() + v0->meshPos(1) * mesh()->zElements() + v0->meshPos(2);
  } else if (sgnx == 0 && sgnz == 0) {
    int base = (mesh()->xElements()+1) * mesh()->yElements() * mesh()->zElements();
    return base + v0->meshPos(1) * mesh()->xElements() * mesh()->zElements() + v0->meshPos(0) * mesh()->zElements() + v0->meshPos(2);
  } else if (sgnx == 0 && sgny == 0) {
    int base = (mesh()->xElements()+1) * mesh()->yElements() * mesh()->zElements() + mesh()->xElements() * (mesh()->yElements()+1) * mesh()->zElements();
    return base + v0->meshPos(2) * mesh()->xElements() * mesh()->yElements() + v0->meshPos(0) * mesh()->yElements() + v0->meshPos(1);
  } else {
    std::cout << "ERROR: Problematic element face position index";
    return -1;
  }
}

int Element::faceGlobal(int i) {
  int sgnx, sgny, sgnz;
  facePos(i, sgnx, sgny, sgnz);
  return faceGlobal(sgnx, sgny, sgnz);
}

std::vector<std::vector<int>>* Element::subtetrahedra(int type) {
  // type=0 : T_D, type=1 : T_M
  switch(type) {
    case 0: {
      return &subtetra_d;
    }
    case 1: {
      return &subtetra_m;
    }
    default: return &subtetra_m;
  }
}

bool Element::isInSubtetrahedron(const std::vector<int>& subtetrahedron, const Point& pt) {
  for (int iVertex=0; iVertex<4; iVertex++) {
    Tensor1 normal = cross(Tensor1(*vertexPtr(subtetrahedron[(iVertex+2)%4])-*vertexPtr(subtetrahedron[(iVertex+1)%4])), 
                           Tensor1(*vertexPtr(subtetrahedron[(iVertex+3)%4])-*vertexPtr(subtetrahedron[(iVertex+1)%4])));
      if((Tensor1(*vertexPtr(subtetrahedron[iVertex])-*vertexPtr(subtetrahedron[(iVertex+1)%4]))*normal) 
       * (Tensor1(pt-*vertexPtr(subtetrahedron[(iVertex+1)%4]))*normal) < -element_eps) return false;
  }
  return true;
}

int Element::inSubtetrahedron(const std::vector<std::vector<int>>& subtetrahedra, const Point& pt) {
  static int previousSubTetraIndex = 0;

  for(unsigned int i=previousSubTetraIndex; i<subtetrahedra.size(); i++) {
    if(isInSubtetrahedron(subtetrahedra[i],pt)) {
      previousSubTetraIndex = i;
      return i;
    }
  }
  for(int i=0; i<previousSubTetraIndex; i++) {
    if(isInSubtetrahedron(subtetrahedra[i],pt)) {
      previousSubTetraIndex = i;
      return i;
    }
  }
  return -1;
}

Tensor1 Element::normal(int i) {
  i = i % 6;
  Tensor1 e01(*faceVertexPtr(0,1,i)-*faceVertexPtr(0,0,i));
  Tensor1 e10(*faceVertexPtr(1,0,i)-*faceVertexPtr(0,0,i));
  if (i == 0 || i == 3 || i == 4) {
    Tensor1 normal_vector = cross(e01,e10);
    return normal_vector/normal_vector.norm();
  } else {
    Tensor1 normal_vector = cross(e10,e01);
    return normal_vector/normal_vector.norm();
  }
}

double Element::lambda(int i, double x, double y, double z) {
  Tensor1 normal_vector = normal(i);
  return (faceVertexPtr(0,0,i)->val(0) - x) * normal_vector[0]
        +(faceVertexPtr(0,0,i)->val(1) - y) * normal_vector[1]
        +(faceVertexPtr(0,0,i)->val(2) - z) * normal_vector[2];
}

Point Element::forwardMap(const Point& p) {
  return ((double)1/8)*( 
      *vertexPtr(0)*(1-p(0))*(1-p(1))*(1-p(2)) + *vertexPtr(4)*(1+p(0))*(1-p(1))*(1-p(2)) 
    + *vertexPtr(2)*(1-p(0))*(1+p(1))*(1-p(2)) + *vertexPtr(6)*(1+p(0))*(1+p(1))*(1-p(2)) 
    + *vertexPtr(1)*(1-p(0))*(1-p(1))*(1+p(2)) + *vertexPtr(5)*(1+p(0))*(1-p(1))*(1+p(2)) 
    + *vertexPtr(3)*(1-p(0))*(1+p(1))*(1+p(2)) + *vertexPtr(7)*(1+p(0))*(1+p(1))*(1+p(2))
    );
}

Point Element::backwardMap(const Point& p) {
  double tol = 1e-12;
  int max_iter = 20;

  Point p_orig((*vertexPtr(0)+*vertexPtr(7))/2); // initial guess
  Tensor1 resid;
  for (int k=0; k<max_iter; k++) {
    resid = forwardMap(p_orig) - p;
    if (resid.norm() < tol) return p_orig;
    Tensor2 df; double det;
    dForwardMap(p_orig,df,det);
    Tensor2 df_inv; df.inverse(df_inv);
    p_orig = p_orig - df_inv*resid;
  }
  return p_orig;
}

void Element::dForwardMap(const Point& p, Tensor2& df, double& jac) {
  Tensor1 col[3];
  col[0] = ((double)1/8)*((-1)*
    (*vertexPtr(0))*(1-p(1))*(1-p(2)) + *vertexPtr(4)*(1-p(1))*(1-p(2)) 
    - *vertexPtr(2)*(1+p(1))*(1-p(2)) + *vertexPtr(6)*(1+p(1))*(1-p(2)) 
    - *vertexPtr(1)*(1-p(1))*(1+p(2)) + *vertexPtr(5)*(1-p(1))*(1+p(2)) 
    - *vertexPtr(3)*(1+p(1))*(1+p(2)) + *vertexPtr(7)*(1+p(1))*(1+p(2))
    );
  col[1] = ((double)1/8)*( (-1)*
    (*vertexPtr(0))*(1-p(0))*(1-p(2)) - *vertexPtr(4)*(1+p(0))*(1-p(2)) 
    + *vertexPtr(2)*(1-p(0))*(1-p(2)) + *vertexPtr(6)*(1+p(0))*(1-p(2)) 
    - *vertexPtr(1)*(1-p(0))*(1+p(2)) - *vertexPtr(5)*(1+p(0))*(1+p(2)) 
    + *vertexPtr(3)*(1-p(0))*(1+p(2)) + *vertexPtr(7)*(1+p(0))*(1+p(2))
    );
  col[2] = ((double)1/8)*( (-1)*
    (*vertexPtr(0))*(1-p(0))*(1-p(1)) - *vertexPtr(4)*(1+p(0))*(1-p(1))
    - *vertexPtr(2)*(1-p(0))*(1+p(1)) - *vertexPtr(6)*(1+p(0))*(1+p(1))
    + *vertexPtr(1)*(1-p(0))*(1-p(1)) + *vertexPtr(5)*(1+p(0))*(1-p(1)) 
    + *vertexPtr(3)*(1-p(0))*(1+p(1)) + *vertexPtr(7)*(1+p(0))*(1+p(1))
    );

  df.set(col[0][0],col[1][0],col[2][0],col[0][1],col[1][1],col[2][1],col[0][2],col[1][2],col[2][2]);
  jac = abs(df.determinant());
}

bool Element::isInElement(const Point& pt) {
  for(int i=0; i<6; i++) {
    if(lambda(i, pt) < -element_eps) return false;
  }
  return true;
}

void Element::write_raw(std::ofstream& fout) {
  fout << "  ELEMENT\n";
  fout << "  good_element  = " << good_element << "\n";
  fout << "  my_mesh       = " << my_mesh << "\n";
  fout << "  my_mesh_index       = " << my_mesh_index << "\n";
  fout << "  my_mesh_pos = (" << my_mesh_pos[0] << "," << my_mesh_pos[1] << "," << my_mesh_pos[2] <<")\n";
  fout << "\nvertices:\n";
  for (int i=0; i<8; i++) {
    int sgnx, sgny, sgnz;
    vertexPos(i, sgnx, sgny, sgnz);
    int i_check=vertexIndex(sgnx, sgny, sgnz);
    fout << "  my_v[" << i_check << "]_pos = (" << sgnx << "," << sgny << "," << sgnz << ")\n";
    vertexPtr(i_check)->write_raw(fout);
  }
  fout << "\nedges:\n";
  for (int i=0; i<12; i++) {
    int sgnx, sgny, sgnz;
    edgePos(i, sgnx, sgny, sgnz);
    int i_check=edgeIndex(sgnx, sgny, sgnz);
    fout << "  edge local index: " << i_check << "\n";
    fout << "    position = (" << sgnx << "," << sgny << "," << sgnz << ")\n";
    fout << "    global index = " << edgeGlobal(i_check) << "\n";
    fout << "    v0:\n\t\t";
    edgeVertexPtr(0, i_check)->write_raw(fout);
    fout << "    v1:\n\t\t";
    edgeVertexPtr(1, i_check)->write_raw(fout);
  }
  fout << "\nfaces:\n";
  for (int i=0; i<6; i++) {
    int sgnx, sgny, sgnz;
    facePos(i, sgnx, sgny, sgnz);
    int i_check=faceIndex(sgnx, sgny, sgnz);    
    fout << "  face local index: " << i_check << "\n";
    fout << "    position = (" << sgnx << "," << sgny << "," << sgnz << ")\n";
    fout << "    global index = " << faceGlobal(i_check) << "\n";   
    fout << "    v00:\n\t\t";
    faceVertexPtr(0, 0, i_check)->write_raw(fout);
    fout << "    v01:\n\t\t";
    faceVertexPtr(0, 1, i_check)->write_raw(fout);
    fout << "    v10:\n\t\t";
    faceVertexPtr(1, 0, i_check)->write_raw(fout);
    fout << "    v11:\n\t\t";
    faceVertexPtr(1, 1, i_check)->write_raw(fout);
  }
}

int Element::write_raw(std::string& filename) {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Class HexaMesh

void HexaMesh::set_hexamesh(int nx, int ny, int nz, Point* vertices) {
  num_x = nx; num_y= ny; num_z = nz;
  num_elements = nx*ny*nz;
  num_vertices = (nx+1)*(ny+1)*(nz+1);
  min_x = 0, min_y = 0, min_z = 0;
  max_x = 1, max_y = 1, max_z = 1;
  good_mesh = false;

  if(the_vertices) delete[] the_vertices;
  the_vertices = new Vertex[num_vertices];
  if(!the_vertices) return;

  if(the_elements) delete[] the_elements;
  the_elements = new Element[num_elements];
  if(!the_elements) return;

  if(face_center) delete[] face_center;
  face_center = new Point[nFaces()];

  if(face_radius) delete[] face_radius;
  face_radius = new double[nFaces()];

  for (int ix=0; ix<nx+1; ix++) {
    for (int iy=0; iy<ny+1; iy++) {
      for (int iz=0; iz<nz+1; iz++) {
        int vertex_index = vertexIndex(ix,iy,iz);
        the_vertices[vertex_index].set(vertices[vertex_index],ix,iy,iz,vertex_index,this);
        if (vertices[vertex_index].val(0)<min_x) min_x=vertices[vertex_index].val(0);
        if (vertices[vertex_index].val(0)>max_x) max_x=vertices[vertex_index].val(0);
        if (vertices[vertex_index].val(1)<min_y) min_y=vertices[vertex_index].val(1);
        if (vertices[vertex_index].val(1)>max_y) max_y=vertices[vertex_index].val(1);
        if (vertices[vertex_index].val(2)<min_z) min_z=vertices[vertex_index].val(2);
        if (vertices[vertex_index].val(2)>max_z) max_z=vertices[vertex_index].val(2);
      }
    }
  }
  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
        int element_index = elementIndex(ix,iy,iz);
        the_elements[element_index].set(ix,iy,iz,element_index,this);
        if (!the_elements[element_index].isGood()) return; 
      }
    }
  }

  Point center; double radius;
  for (int i=0; i<nFaces(); i++) {
    faceCenter(i, center, radius);
    face_center[i] = center;
    face_radius[i] = radius;
  }
  good_mesh = true;
}

HexaMesh::HexaMesh(Element* single_element){
  Point vertices[8];
  for (int i=0; i<8; i++) {
    vertices[i] = Point(single_element->vertexPtr(i)->val(0), single_element->vertexPtr(i)->val(1), single_element->vertexPtr(i)->val(2));
  }
  set(1,1,1,vertices);
}

int HexaMesh::createMesh(char meshTypeC, int nx, int ny, int nz, 
                      double xMin, double xMax, double yMin, double yMax, double zMin, double zMax, double distortionFactor) {
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-1,1); //doubles from -1 to 1

  Point vertices[(nx+1)*(ny+1)*(nz+1)];
  double h_x = (xMax-xMin)/nx, h_y = (yMax-yMin)/ny, h_z = (zMax-zMin)/nz;
    
  for (int ix=0; ix<nx+1; ix++) {
    for (int iy=0; iy<ny+1; iy++) {
      for (int iz=0; iz<nz+1; iz++) {

        double px = 0, py = 0, pz = 0;
        // set vertices
        if (meshTypeC != 'r') {
          // deviated mesh
          px = xMin+ix*h_x, py = yMin+iy*h_y, pz = zMin+iz*h_z;

          if(0 < ix && ix < nx) {
            int parity_xy = ( 2*(iy % 2) - 1 ) * ( 2*(ix % 2) - 1 );
            px += parity_xy*distortionFactor*h_x;
          }

          if (0 < iz && iz < nz && meshTypeC == 't') {
            int parity_yz = ( 2*(iy % 2) - 1 ) * ( 2*(iz % 2) - 1 );
            pz += parity_yz*distortionFactor*h_z;
          }

        } else {
          // randomly deviated

          // Firstly, we move all the points on the 'lower' boundary
          if (0 < ix && ix < nx && ( iy == 0 || iz == 0 )) {
            double rand_x = distribution(generator);
            px = xMin+ix*h_x + distortionFactor*rand_x*h_x;
          }
          if (0 < iy && iy < ny && ( ix == 0 || iz == 0 )) {
            double rand_y = distribution(generator);
            py = yMin+iy*h_y + distortionFactor*rand_y*h_y;
          }
          if (0 < iz && iz < nz && ( ix == 0 || iy == 0)) {
            double rand_z = distribution(generator);
            pz = zMin+iz*h_z + distortionFactor*rand_z*h_z;              
          }
          
          // All the other points are uniquely fixed by previous points
          if (0 < ix && 0 < iy && 0 < iz) {
            Point p = flatHex(vertices[ (ix-1)*(ny+1)*(nz+1)+(iy-1)*(nz+1)+iz-1 ], vertices[ (ix-1)*(ny+1)*(nz+1)+(iy-1)*(nz+1)+iz ],
                              vertices[ (ix-1)*(ny+1)*(nz+1)+(iy)*(nz+1)+iz-1 ], vertices[ (ix-1)*(ny+1)*(nz+1)+(iy)*(nz+1)+iz ],
                              vertices[ (ix)*(ny+1)*(nz+1)+(iy-1)*(nz+1)+iz-1 ], vertices[ (ix)*(ny+1)*(nz+1)+(iy-1)*(nz+1)+iz ], vertices[ (ix)*(ny+1)*(nz+1)+(iy)*(nz+1)+iz-1 ]);
            px = p.val(0); py = p.val(1); pz = p.val(2);
          }

          // To eliminate the machine error on the boundary, 
          // we force all the points on the boundary to keep in the boundary          
          if (ix == 0) { px = xMin; }
          if (iy == 0) { py = yMin; }
          if (iz == 0) { pz = zMin; }
          if (ix == nx) { px = xMax; }
          if (iy == ny) { py = yMax; }
          if (iz == nz) { pz = zMax; }
        }
        vertices[ ix*(ny+1)*(nz+1)+iy*(nz+1)+iz ].set(px,py,pz);
      }
    }
  }
  set(nx,ny,nz,vertices);

  if (!isGood()) return -1;

  return 0;
}

HexaMesh::~HexaMesh() {
  if(the_vertices) delete[] the_vertices;
  if(the_elements) delete[] the_elements;
  if(face_center) delete[] face_center;
  if(face_radius) delete[] face_radius;

}

Point HexaMesh::flatHex(const Point& v000, const Point& v001, 
                        const Point& v010, const Point& v011, 
                        const Point& v100, const Point& v101, const Point& v110) {
  Tensor1 nu_f100 = cross(Tensor1(v110-v100), Tensor1(v101-v100)); //nu2
  Tensor1 nu_f010 = cross(Tensor1(v011-v010), Tensor1(v110-v010)); //nu4
  Tensor1 nu_f001 = cross(Tensor1(v101-v001), Tensor1(v011-v001)); //nu6

  Tensor2 A(nu_f001.val(0),nu_f001.val(1),nu_f001.val(2), 
            nu_f010.val(0),nu_f010.val(1),nu_f010.val(2),
            nu_f100.val(0),nu_f100.val(1),nu_f100.val(2));

  Tensor2 invA;

  (void) A.inverse(invA);
  Tensor1 v111 = invA * Tensor1(nu_f001*v001, nu_f010*v010, nu_f100*v100);
  return v111;
}

void HexaMesh::potentialFaceCenter(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3, 
                        Point& center_candidate, double& radius, bool& inside) {
  Tensor1 t10(*v0-*v1); t10 /= t10.norm();
  Tensor1 t12(*v2-*v1); t12 /= t12.norm();
  Tensor1 t1 = t10+t12; t1 /= t1.norm();

  Tensor1 t21(*v1-*v2); t21 /= t21.norm();
  Tensor1 t23(*v3-*v2); t23 /= t23.norm();
  Tensor1 t2 = t21+t23; t2 /= t2.norm();

  // Pick the two coordinates we are going to use to solve for the center
  int i0,i1;
  if (v0->meshPos(0) == v1-> meshPos(0) && v1->meshPos(0) == v2-> meshPos(0)) { i0 = 1; i1 = 2; }
  if (v0->meshPos(1) == v1-> meshPos(1) && v1->meshPos(1) == v2-> meshPos(1)) { i0 = 0; i1 = 2; }
  if (v0->meshPos(2) == v1-> meshPos(2) && v1->meshPos(2) == v2-> meshPos(2)) { i0 = 0; i1 = 1; }
  
  double c = (t2[i1]*(v2->val(i0)) + t2[i0]*(v1->val(i1)) - t2[i1]*(v1->val(i0)) - t2[i0]*(v2->val(i1)))
           / (t1[i0]*t2[i1]-t1[i1]*t2[i0]);

  center_candidate = *v1+c*t1;
  radius = c * cross(t1,t10).norm();

  Tensor1 inner_normal = cross(cross(-1*t21,t23),*v0-*v3);  inner_normal /= inner_normal.norm();
  double dist = Tensor1(center_candidate-*v0) * inner_normal;
  if (dist + element_eps < radius) { 
    inside = false; 
  } else {
    inside = true;
  } 

  return;
}

void HexaMesh::vertexPos(int i, int& ix, int& iy, int& iz) {
  i = i % num_vertices;
  ix = i / ( (num_y+1)*(num_z+1) );
  iy = (i - ix*(num_y+1)*(num_z+1)) / (num_z+1);
  iz = i - ix*(num_y+1)*(num_z+1) - iy*(num_z+1);
  return;
}

void HexaMesh::elementPos(int i, int& ix, int& iy, int& iz) {
  i = i % num_elements;
  ix = i / num_y*num_z;
  iy = (i - ix*num_y*num_z) / num_z;
  iz = i - ix*num_y*num_z - iy*num_z;
  return;
}

// Map the position index (ix, iy, iz) to the array index i
int HexaMesh::vertexIndex(int ix, int iy, int iz) const {
  ix = ix % (num_x+1);
  iy = iy % (num_y+1);
  iz = iz % (num_z+1);
  return ix*(num_y+1)*(num_z+1) + iy*(num_z+1) + iz;
}

int HexaMesh::elementIndex(int ix, int iy, int iz) const {
  ix = ix % num_x;
  iy = iy % num_y;
  iz = iz % num_z;
  return ix*num_y*num_z + iy*num_z + iz;
}


// Edges: The indexing of the edges are started with z-dir edges (with index priority z -> y -> x), 
//                                       followed by y-dir edges (with index priority y -> z -> x),
//                                       followed by x-dir edges (with index priority x -> z -> y)
//
//        z-dir edges (the edges intersecting both f_{\pm 3}) with v_0 (the intersection with f_{-3}):
//           (0,0,0) -> (0,0,1) -> ... -> (0,0,nz-1) 
//        -> (0,1,0) -> (0,1,1) -> ... -> (0,1,nz-1) -> ... -> (0,ny,nz-1)
//        -> (1,0,0) -> (1,0,1) -> ... -> (1,0,nz-1)
//        -> (1,1,0) -> (1,1,1) -> ... -> (1,1,nz-1) -> ... -> (1,ny,nz-1) -> ... -> (nx,ny,nz-1)
//
//        y-dir edges (the edges intersecting both f_{\pm 3}) with v_0 (the intersection with f_{-2}):
//           (0,0,0) -> (0,1,0) -> ... -> (0,ny-1,0)
//        -> (0,0,1) -> (0,1,1) -> ... -> (0,ny-1,1) -> ... -> (0,ny-1,nz)
//        -> (1,0,0) -> (1,1,0) -> ... -> (1,ny-1,0)
//        -> (1,0,1) -> (1,1,1) -> ... -> (1,ny-1,1) -> ... -> (1,ny-1,nz) -> ... ->  (nx,ny-1,nz)
//
//        x-dir edges (the edges intersecting both f_{\pm 1}) with v_0 (the intersection with f_{-1}):
//           (0,0,0) -> (1,0,0) -> ... -> (nx-1,0,0)
//        -> (0,0,1) -> (1,0,1) -> ... -> (nx-1,0,1) -> ... -> (nx-1,0,nz)
//        -> (0,1,0) -> (1,1,0) -> ... -> (nx-1,1,0)
//        -> (0,1,1) -> (1,1,1) -> ... -> (nx-1,1,1) -> ... -> (nx-1,1,nz) -> ... ->  (nx-1,ny,nz)

Vertex* HexaMesh::edgeVertexPtr(bool pm, int i) {
  i = i % nEdges();
  int base_x, base_y, base_z;
  if (i<(num_x+1)*(num_y+1)*num_z) {
    // z-dir edges
    base_x = i / ((num_y+1)*num_z);
    base_y = (i - base_x*(num_y+1)*num_z) / num_z;
    base_z = i % num_z;
    return vertexPtr(base_x, base_y, base_z+pm);
  } else if (i<(num_x+1)*(num_y+1)*num_z+(num_x+1)*num_y*(num_z+1)) {
    // y-dir edges
    i -= (num_x+1)*(num_y+1)*num_z;
    base_x = i / (num_y*(num_z+1));
    base_y = i % num_y;
    base_z = (i - base_x*num_y*(num_z+1)) / num_y;
    return vertexPtr(base_x, base_y+pm, base_z);
  } else {
    // x-dir edges
    i -= (num_x+1)*(num_y+1)*num_z+(num_x+1)*num_y*(num_z+1);
    base_x = i % num_x;
    base_y = i / (num_x*(num_z+1));
    base_z = (i - base_y*num_x*(num_z+1)) / num_x;
    return vertexPtr(base_x+pm, base_y, base_z);
  }
}

// Faces: The indexing of the edges are started with yz-dir faces (with index priority z -> y -> x),
//                                       followed by xz-dir faces (with index priority z -> x -> y),
//                                       followed by xy-dir faces (with index priority y -> x -> z)
//
//        yz-dir faces (locally f_{\pm 1}) with v_00 (the left-bottom-back corner vertex):
//           (0,0,0) -> (0,0,1) -> ... -> (0,0,nz-1) 
//        -> (0,1,0) -> (0,1,1) -> ... -> (0,1,nz-1) -> ... -> (0,ny-1,nz-1)
//        -> (1,0,0) -> (1,0,1) -> ... -> (1,0,nz-1)
//        -> (1,1,0) -> (1,1,1) -> ... -> (1,1,nz-1) -> ... -> (1,ny-1,nz-1) -> ... -> (nx,ny-1,nz-1)
//
//        xz-dir faces (locally f_{\pm 2}) with v_00 (the left-bottom-back corner vertex):
//           (0,0,0) -> (0,0,1) -> ... -> (0,0,nz-1) 
//        -> (1,0,0) -> (1,0,1) -> ... -> (1,0,nz-1) -> ... -> (nx-1,0,nz-1)
//        -> (0,1,0) -> (0,1,1) -> ... -> (0,1,nz-1)
//        -> (1,1,0) -> (1,1,1) -> ... -> (1,1,nz-1) -> ... -> (nx-1,1,nz-1) -> ... ->  (nx-1,ny,nz-1)
//
//        xy-dir faces (locally f_{\pm 3}) with v_00 (the left-bottom-back corner vertex):
//           (0,0,0) -> (0,1,0) -> ... -> (0,ny-1,0)
//        -> (1,0,0) -> (1,1,0) -> ... -> (1,ny-1,0) -> ... -> (nx-1,ny-1,0)
//           (0,0,1) -> (0,1,1) -> ... -> (0,ny-1,1)
//        -> (1,0,1) -> (1,1,1) -> ... -> (1,ny-1,1) -> ... -> (nx-1,ny-1,1) -> ... ->  (nx-1,ny-1,nz)

Vertex* HexaMesh::faceVertexPtr(bool pm0, bool pm1, int i) {
  i = i % nFaces();
  int base_x, base_y, base_z;
  if (i < (num_x+1)*num_y*num_z) {
    // yz-dir faces
    base_x = i / (num_y*num_z);
    base_y = (i - base_x*num_y*num_z) / num_z;
    base_z = i % num_z;
    return vertexPtr(base_x, base_y+pm0, base_z+pm1);
  } else if (i < (num_x+1)*num_y*num_z+num_x*(num_y+1)*num_z) {
    i -= (num_x+1)*num_y*num_z;
    base_y = i / (num_x*num_z);
    base_x = (i - base_y*num_x*num_z) / num_z;
    base_z = i % num_z;
    return vertexPtr(base_x+pm0, base_y, base_z+pm1);
  } else {
    i -= (num_x+1)*num_y*num_z+num_x*(num_y+1)*num_z;
    base_z = i / (num_x*num_y);
    base_x = (i - base_z*num_x*num_y) / num_y;
    base_y = i % num_y;
    return vertexPtr(base_x+pm0, base_y+pm1, base_z);
  }
}

void HexaMesh::faceCenter(int i, Point& center, double& radius) {
  // Compute center of largest inscribed circle
  double max_radius = 0;
  std::vector<Point> center_candidates;
  center_candidates.clear();

  Point candidate; double new_radius; bool inside;
  Vertex* vertex_list[4] = { faceVertexPtr(0,0,i), faceVertexPtr(1,0,i), faceVertexPtr(1,1,i), faceVertexPtr(0,1,i) };

  for (int iVertex=0; iVertex<4; iVertex++) {
    potentialFaceCenter(vertex_list[iVertex], vertex_list[(iVertex+1)%4], vertex_list[(iVertex+2)%4], vertex_list[(iVertex+3)%4],
                    candidate, new_radius, inside);
  	if(inside) {
	    if(max_radius < new_radius-element_eps) {
	      max_radius = new_radius;
	      center_candidates.clear();
	      center_candidates.push_back(candidate);
	    } else if (std::abs(max_radius-new_radius) < element_eps) {
	      center_candidates.push_back(candidate);
	    }
	  }  
  }

  center.set(0,0,0);
  if(center_candidates.size() == 0) {
    std::cerr << "ERROR: the face f[" << i << "] has no center of the largest inscribed circle" << std::endl;
  } else {
    for(unsigned long int i=0; i<center_candidates.size(); i++) {
      center += center_candidates[i];
    }
    center /= center_candidates.size();	  
  }
  radius = max_radius;
  return;
}

int HexaMesh::inElement(const Point& pt) {
  static int previousElementIndex = 0;

  for(int i=previousElementIndex; i<num_elements; i++) {
    if(the_elements[i].isInElement(pt)) {
      previousElementIndex = i;
      return i;
    }
  }
  for(int i=0; i<previousElementIndex; i++) {
    if(the_elements[i].isInElement(pt)) {
      previousElementIndex = i;
      return i;
    }
  }
  return -1;
}

void HexaMesh::write_raw(std::ofstream& fout) {
  fout << "HEXAMESH\n";
  fout << "good_mesh  = " << good_mesh << "\n";
  fout << "(num_x, num_y, num_z) = (" << num_x << "," << num_y << "," << num_z << ")\n";
  fout << "[min_x, max_x] x [min_y, max_y] x [min_z, max_z] = [" << min_x << "," << max_x
       << "]x[" << min_y <<"," << max_y << "]x[" << min_z <<"," << max_z << "]\n";
  
  fout << "\nmesh vertices:\n";
  for(int i=0; i<num_vertices; i++) {
    the_vertices[i].write_raw(fout);
  }
  fout << "\nmesh elements:\n";
  for(int i=0; i<num_elements; i++) {
    the_elements[i].write_raw(fout);
    fout << "\n";
  }
  fout << "\nmesh faces:\n";
  for (int i=0; i<nFaces(); i++) {
    fout << " face_index: " << i << "\n";
    fout << " face_vertices: " <<"\n";
    faceVertexPtr(0,0,i)->write_raw(fout);
    faceVertexPtr(0,1,i)->write_raw(fout);
    faceVertexPtr(1,0,i)->write_raw(fout);
    faceVertexPtr(1,1,i)->write_raw(fout);
    fout << " face_center: " << face_center[i] << "\n";
    fout << " face_radius: " << face_radius[i] << "\n";
  }
}

int HexaMesh::write_raw(std::string& filename) {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

int HexaMesh::write_matlab(std::string& filename) const {
  std::ofstream fout(filename + ".m");
  if( !fout ) return 1;

  // Define the color of faces
  std::string color = "[0.6471 0.8706 0.8941]";

  fout << "clf;\n";

  int i = 1; // face index

  for (int ix=0; ix<num_x+1; ix++) {
    for (int iy=0; iy<num_y+1; iy++) {
      for (int iz=0; iz<num_z+1; iz++) {
        if (iy < num_y && iz < num_z) {
          fout << "x" << i << " = [" << the_vertices[ vertexIndex(ix,iy,iz) ].val(0) << " "
                                     << the_vertices[ vertexIndex(ix,iy+1,iz) ].val(0) << " "
                                     << the_vertices[ vertexIndex(ix,iy+1,iz+1) ].val(0) << " "
                                     << the_vertices[ vertexIndex(ix,iy,iz+1) ].val(0) << "];\n";
          fout << "y" << i << " = [" << the_vertices[ vertexIndex(ix,iy,iz) ].val(1) << " "
                                     << the_vertices[ vertexIndex(ix,iy+1,iz) ].val(1) << " "
                                     << the_vertices[ vertexIndex(ix,iy+1,iz+1) ].val(1) << " "
                                     << the_vertices[ vertexIndex(ix,iy,iz+1) ].val(1) << "];\n";
          fout << "z" << i << " = [" << the_vertices[ vertexIndex(ix,iy,iz) ].val(2) << " "
                                     << the_vertices[ vertexIndex(ix,iy+1,iz) ].val(2) << " "
                                     << the_vertices[ vertexIndex(ix,iy+1,iz+1) ].val(2) << " "
                                     << the_vertices[ vertexIndex(ix,iy,iz+1) ].val(2) << "];\n";
          i += 1;
        }

        if (ix < num_x && iz < num_z) {
          fout << "x" << i << " = [" << the_vertices[ vertexIndex(ix,iy,iz) ].val(0) << " "
                                     << the_vertices[ vertexIndex(ix+1,iy,iz) ].val(0) << " "
                                     << the_vertices[ vertexIndex(ix+1,iy,iz+1) ].val(0) << " "
                                     << the_vertices[ vertexIndex(ix,iy,iz+1) ].val(0) << "];\n";
          fout << "y" << i << " = [" << the_vertices[ vertexIndex(ix,iy,iz) ].val(1) << " "
                                     << the_vertices[ vertexIndex(ix+1,iy,iz) ].val(1) << " "
                                     << the_vertices[ vertexIndex(ix+1,iy,iz+1) ].val(1) << " "
                                     << the_vertices[ vertexIndex(ix,iy,iz+1) ].val(1) << "];\n";
          fout << "z" << i << " = [" << the_vertices[ vertexIndex(ix,iy,iz) ].val(2) << " "
                                     << the_vertices[ vertexIndex(ix+1,iy,iz) ].val(2) << " "
                                     << the_vertices[ vertexIndex(ix+1,iy,iz+1) ].val(2) << " "
                                     << the_vertices[ vertexIndex(ix,iy,iz+1) ].val(2) << "];\n";
          i += 1;
        }

        if (ix < num_x && iy < num_y) {
          fout << "x" << i << " = [" << the_vertices[ vertexIndex(ix,iy,iz) ].val(0) << " "
                                     << the_vertices[ vertexIndex(ix+1,iy,iz) ].val(0) << " "
                                     << the_vertices[ vertexIndex(ix+1,iy+1,iz) ].val(0) << " "
                                     << the_vertices[ vertexIndex(ix,iy+1,iz) ].val(0) << "];\n";
          fout << "y" << i << " = [" << the_vertices[ vertexIndex(ix,iy,iz) ].val(1) << " "
                                     << the_vertices[ vertexIndex(ix+1,iy,iz) ].val(1) << " "
                                     << the_vertices[ vertexIndex(ix+1,iy+1,iz) ].val(1) << " "
                                     << the_vertices[ vertexIndex(ix,iy+1,iz) ].val(1) << "];\n";
          fout << "z" << i << " = [" << the_vertices[ vertexIndex(ix,iy,iz) ].val(2) << " "
                                     << the_vertices[ vertexIndex(ix+1,iy,iz) ].val(2) << " "
                                     << the_vertices[ vertexIndex(ix+1,iy+1,iz) ].val(2) << " "
                                     << the_vertices[ vertexIndex(ix,iy+1,iz) ].val(2) << "];\n";
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
  fout << "hold off;\n";
  return 0;
}
