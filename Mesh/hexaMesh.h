#ifndef __hexamesh_h_included__
#define __hexamesh_h_included__

#include <string>
#include <fstream>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
// HexaMesh class
//   Includes classes for Edge, Face and HexaElement
//   Vertices are base objects called Point
// Mesh of convex hexagons in 3D
//
// Assume base objects: Point and Tensor1 (both three numbers)
////////////////////////////////////////////////////////////////////////////////

#include "baseObjects.h"

namespace hexamesh {

  class HexaMesh;

  ////////////////////////////////////////////////////////////////////////////////
  // class Vertex
  ////////////////////////////////////////////////////////////////////////////////

  class Vertex : public Point
  {
  private:
    int my_mesh_index[3];
    HexaMesh* my_mesh;
    
    void set_vertex(double px, double py, double pz, int ix=-1, int iy=-1, int iz=-1,  HexaMesh* myMesh=nullptr);

  public:
    Vertex() { set_vertex(0,0,0); };
    Vertex(double px, double py, double pz, int ix=-1, int iy=-1, int iz=-1, HexaMesh* myMesh=nullptr) {
      set_vertex(px, py, pz, ix, iy, iz, myMesh); };
    Vertex(const Point& p, int ix=-1, int iy=-1, int iz=-1, HexaMesh* myMesh=nullptr) {
      set_vertex(p[0], p[1], p[2], ix, iy, iz, myMesh); };
    Vertex(const Point* p, int ix=-1, int iy=-1, int iz=-1, HexaMesh* myMesh=nullptr) {
      set_vertex((*p)[0], (*p)[1], (*p)[2], ix, iy, iz, myMesh); };

    void set() { set_vertex(0,0,0); };
    void set(double px=0, double py=0, double pz=0, int ix=-1, int iy=-1, int iz=-1, HexaMesh* myMesh=nullptr) {
      set_vertex(px, py, pz, ix, iy, iz, myMesh); };
    void set(const Point& p, int ix=-1, int iy=-1, int iz=-1, HexaMesh* myMesh=nullptr) {
      set_vertex(p[0], p[1], p[2], ix, iy, iz, myMesh); };
    void set(const Point* p, int ix=-1, int iy=-1, int iz=-1, HexaMesh* myMesh=nullptr) {
      set_vertex((*p)[0], (*p)[1], (*p)[2], ix, iy, iz, myMesh); };

    HexaMesh* mesh() const { return my_mesh; };
    int meshIndex(int i) const { return my_mesh_index[i % 3]; };

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };


  class Element
  {
  private:
    int my_mesh_index[3];

    // double the_volume;
    // Point my_center;  // of largest inscribed ball
    // Point my_centroid; // center of mass
    // double my_diameter;
    // double max_radius;

    HexaMesh* my_mesh = nullptr;

    //Vertex* my_vertices;
    // double computeAreaTetra(const Point* px, const Point* py, const Point* pz, const Point* p3);
    // double computeArea();
    // void computeCenterOfTetra(const Point* px, const Point* py, const Point* pz, const Point* p3, Point& center, bool& unique);

    void set_element(int ix, int iy, int iz, HexaMesh* myMesh=nullptr);

  public:
    Element(int ix, int iy, int iz, HexaMesh* myMesh=nullptr) {
      set_element(ix, iy, iz, myMesh); };
    Element() :
		    // the_volume(0), my_center(), my_centroid(), my_diameter(0), max_radius(0),
		    my_mesh_index(), my_mesh(nullptr) {};

    void set(int ix, int iy, int iz, HexaMesh* myMesh=nullptr) {
      set_element(ix, iy, iz, myMesh); };
    
    //Vertex vertexPtr(int ix, int iy, int iz) { return the_vertices[i % 8]; };

    HexaMesh* mesh() const { return my_mesh; };
    int meshIndex(int i) const { return my_mesh_index[i % 3]; };
    
    // double volume() const { return the_volume; };
    // Point center()   const { return my_center; };
    // Point centroid() const { return my_centroid; };
    // double maxRadius() const { return max_radius; }
    // double diameter() const { return my_diameter; };

    // bool isInElement(const Point& pt) const;
    // bool isOnElementBoundary(const Point& pt) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    friend class HexaMesh;
  };

  class HexaMesh
  {
  private:
    char mesh_type;

    int num_x, num_y, num_z;
    double min_x, max_x, min_y, max_y, min_z, max_z;
    
    int num_elements;
    int num_vertices;
    
    double distortion_factor;
    
    
    Vertex* the_vertices = nullptr;
    Element* the_elements = nullptr;

    // It return the last point of a hexahedron, given all the other 7 points
    Point flatHex(const Point& v000, const Point& v001, 
                  const Point& v010, const Point& v011, 
                  const Point& v100, const Point& v101, const Point& v110);

    void set_hexamesh(char meshType, int nx, int ny, int nz, 
                      double xMin, double xMax, double yMin, double yMax, double zMin, double zMax, double distortionFactor=0);

  public:
    HexaMesh(char meshType, int nx, int ny, int nz, 
             double xMin, double xMax, double yMin, double yMax, double zMin, double zMax, double distortionFactor=0) {
      set_hexamesh(meshType, nx,ny,nz, xMin, xMax, yMin, yMax, zMin, zMax, distortionFactor); };
    HexaMesh() : mesh_type('d'), num_x(1), num_y(1), num_z(1), 
                 min_x(0), max_x(1), min_y(0), max_y(1), min_z(0), max_z(1), 
                 num_elements(0), num_vertices(0), distortion_factor(0) {};
    //HexaMesh(Element* single_element); // Create a one element mesh from an element
    ~HexaMesh();

    void set(char mesh_type, int nx=1, int ny=1, int nz=1, 
             double xMin=0, double xMax=1, double yMin=0, double yMax=1, double zMin=0, double zMax=1, double distortionFactor=0) {
      set_hexamesh(mesh_type, nx,ny,nz, xMin, xMax, yMin, yMax, zMin, zMax, distortionFactor); };


    // Map the array index i to the position index (ix, iy, iz)
    void vertexPos(int i, int& ix, int& iy, int& iz);
    void elementPos(int i, int& ix, int& iy, int& iz);

    // Map the position index (ix, iy, iz) to the array index i
    int vertexIndex(int ix, int iy, int iz) const;
    int elementIndex(int ix, int iy, int iz) const;

    // Access functions

    // Access the element by its position in the array
    Element* elementPtr(int i) { return &the_elements[i % num_elements]; };
    // Access the element by its geometric position
    Element* elementPtr(int ix, int iy, int iz) { return elementPtr( elementIndex(ix,iy,iz) ); };

    // Access the vertex by its position in the array
    Vertex* vertexPtr(int i) { return &the_vertices[i % num_elements]; };
    // Access the vertex by its geometric position
    Vertex* vertexPtr(int ix, int iy, int iz) { return vertexPtr( vertexIndex(ix,iy,iz) ); };

    // Output
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    int write_matlab(std::string& filename) const;

    friend class Vertex;
    friend class Element;
  };
}
#endif