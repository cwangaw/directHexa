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

void Vertex::set_vertex(double px, double py, double pz, int ix, int iy, int iz,  HexaMesh* myMesh) {
  the_point[0] = px; 
  the_point[1] = py; 
  the_point[2] = pz; 
  my_mesh_index[0] = ix; my_mesh_index[1] = iy; my_mesh_index[2] = iz; 
  my_mesh = myMesh; 
};

void Vertex::write_raw(std::ofstream& fout) const {
  fout << "      VERTEX (" << val(0) << "," <<val(1) << "," << val(2) <<")\n";
  fout << "      my_mesh       = " << my_mesh << "\n";
  fout << "      my_mesh_index = (" << my_mesh_index[0] << "," << my_mesh_index[1] << "," << my_mesh_index[2] <<")\n";
}

int Vertex::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// Class Element

void Element::set_element(int ix, int iy, int iz, HexaMesh* myMesh) {
  my_mesh_index[0] = ix;
  my_mesh_index[1] = iy;
  my_mesh_index[2] = iz;
  my_mesh = myMesh;

  // the_volume;
  // my_center;  // of largest inscribed ball
  // my_centroid; // center of mass
  // my_diameter;
  // max_radius;
}

void Element::write_raw(std::ofstream& fout) const {
  fout << "  ELEMENT\n";
  // fout << "  the_volume    = " << the_volume << "\n";
  // fout << "  my_center     = " << my_center << "\n";
  // fout << "  my_centroid   = " << my_centroid << "\n";
  fout << "  my_mesh       = " << my_mesh << "\n";
  fout << "  my_mesh_index = (" << my_mesh_index[0] << "," << my_mesh_index[1] << "," << my_mesh_index[2] <<")\n";

  /*
  for(int i=0; i<num_vertices; i++) {
    fout << "  the_oriented_edge "<< i << ":\n";
    the_oriented_edge[i].write_raw(fout);
    fout << "\n";
  }
  */
  
}

int Element::write_raw(std::string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Class HexaMesh

void HexaMesh::set_hexamesh(char meshTypeC, int nx, int ny, int nz, 
                      double xMin, double xMax, double yMin, double yMax, double zMin, double zMax, double distortionFactor) {
  min_x = xMin; max_x = xMax; min_y = yMin; max_y = yMax; min_z = zMin; max_z = zMax;
  num_x = nx; num_y= ny; num_z = nz;
  num_elements = nx * ny * nz;
  num_vertices = (nx+1) * (ny+1) * (nz+1);

  if(the_vertices) delete[] the_vertices;
  the_vertices = new Vertex[num_vertices];
  if(!the_vertices) return;

  if(the_elements) delete[] the_elements;
  the_elements = new Element[num_elements];
  if(!the_elements) return;

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-1,1); //doubles from -1 to 1

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
            Point p = flatHex(the_vertices[ vertexIndex(ix-1,iy-1,iz-1) ], the_vertices[ vertexIndex(ix-1,iy-1,iz) ],
                              the_vertices[ vertexIndex(ix-1,iy,iz-1) ], the_vertices[ vertexIndex(ix-1,iy,iz) ],
                              the_vertices[ vertexIndex(ix,iy-1,iz-1) ], the_vertices[ vertexIndex(ix,iy-1,iz) ], the_vertices[ vertexIndex(ix,iy,iz-1) ]);
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
    
        the_vertices[ vertexIndex(ix,iy,iz) ].set(px,py,pz,ix,iy,iz,this);
        // set elements
        the_elements[ elementIndex(ix,iy,iz) ].set(ix,iy,iz,this);
      }
    }
  }
}

HexaMesh::~HexaMesh() {
  if(the_vertices) delete[] the_vertices;
  if(the_elements) delete[] the_elements;
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
};

void HexaMesh::write_raw(std::ofstream& fout) const {
  fout << "HEXAMESH\n";
  fout << "(num_x, num_y, num_z) = (" << num_x << "," << num_y << "," << num_z << ")\n";
  fout << "[min_x, max_x] x [min_y, max_y] x [min_z, max_z] = [" << min_x << "," << max_x
       << "]x[" << min_y <<"," << max_y << "]x[" << min_z <<"," << max_z << "]\n";
  
  fout << "\nmesh vertices:\n";
  for(int i=0; i<num_vertices; i++) {
    the_vertices[i].write_raw(fout);
  }

  for(int i=0; i<num_elements; i++) {
    fout << "\nmesh elements "<<i<<":\n";
    the_elements[i].write_raw(fout);
    fout << "\n";
  }
}

int HexaMesh::write_raw(std::string& filename) const {
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
