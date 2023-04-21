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
    int my_mesh_pos[3];
    int my_mesh_index;
    HexaMesh* my_mesh;
    
    void set_vertex(double px, double py, double pz, int ix=-1, int iy=-1, int iz=-1, int i=-1, HexaMesh* myMesh=nullptr);

  public:
    Vertex() { set_vertex(0,0,0); };
    Vertex(double px, double py, double pz, int ix=-1, int iy=-1, int iz=-1, int i=-1, HexaMesh* myMesh=nullptr) {
      set_vertex(px, py, pz, ix, iy, iz, i, myMesh); };
    Vertex(const Point& p, int ix=-1, int iy=-1, int iz=-1, int i=-1, HexaMesh* myMesh=nullptr) {
      set_vertex(p[0], p[1], p[2], ix, iy, iz, i, myMesh); };
    Vertex(const Point* p, int ix=-1, int iy=-1, int iz=-1, int i=-1, HexaMesh* myMesh=nullptr) {
      set_vertex((*p)[0], (*p)[1], (*p)[2], ix, iy, iz, i, myMesh); };

    void set() { set_vertex(0,0,0); };
    void set(double px=0, double py=0, double pz=0, int ix=-1, int iy=-1, int iz=-1, int i=-1, HexaMesh* myMesh=nullptr) {
      set_vertex(px, py, pz, ix, iy, iz, i, myMesh); };
    void set(const Point& p, int ix=-1, int iy=-1, int iz=-1, int i=-1, HexaMesh* myMesh=nullptr) {
      set_vertex(p[0], p[1], p[2], ix, iy, iz, i, myMesh); };
    void set(const Point* p, int ix=-1, int iy=-1, int iz=-1, int i=-1, HexaMesh* myMesh=nullptr) {
      set_vertex((*p)[0], (*p)[1], (*p)[2], ix, iy, iz, i, myMesh); };

    HexaMesh* mesh() const { return my_mesh; };
    int meshPos(int i) const { return my_mesh_pos[i % 3]; };
    int meshIndex() const { return my_mesh_index; };

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };

 ////////////////////////////////////////////////////////////////////////////////
  // class Element
  //
  // Local indexing
  //
  //      1------3           o------o
  //     /|     /|          /|     /|      5  0
  //    / |    / |         / |    / |      | /
  //   5------7  |        o------o  |      |/
  //   |  |   |  |        |  |   |  |  2---*---3
  //   |  0---|--2        |  o---|--o     /|
  //   | /    | /         | /    | /     / |
  //   |/     |/          |/     |/     1  4
  //   4------6           o------o
  //    Vertices                  Faces
  //
  //           o--5----o
  //          9|     11|      
  //         / |     / |      
  //        o----7--o  1     
  //        |  0    |  |  
  //        2  |    3  |
  //        |  o--4-|--o    
  //        | 8     | 10     
  //        |/      |/   
  //        o----6--o
  //          Edges
  //
  // The reference cube is [-1,1]^3
  ////////////////////////////////////////////////////////////////////////////////


  class Element
  {
  private:
    int my_mesh_pos[3];
    int my_mesh_index;
    bool good_element;
    HexaMesh* my_mesh = nullptr;
    int* vertex_global_index = nullptr;
    std::vector<std::vector<int>> subtetra_d;
    std::vector<std::vector<int>> subtetra_m;

    void set_element(int ix, int iy, int iz, int i, HexaMesh* myMesh=nullptr);

  public:
    Element(int ix, int iy, int iz, int i, HexaMesh* myMesh=nullptr) {
      set_element(ix, iy, iz, i, myMesh); };
    Element() :
		    my_mesh_pos(), my_mesh_index(), good_element(0), subtetra_d(), subtetra_m() {};
    ~Element();

    void set(int ix, int iy, int iz, int i, HexaMesh* myMesh=nullptr) {
      set_element(ix, iy, iz, i, myMesh); };
    
    HexaMesh* mesh() const { return my_mesh; };
    int meshPos(int i) const { return my_mesh_pos[i % 3]; };
    int meshIndex() const { return my_mesh_index; };
    bool isGood() const { return good_element; };

    // Access functions and index mapping functions

    // Vertices
    // (sgnx, sgny, sgnz) = (-1,-1,-1) -> f_{\hat{x}=-1} intersect f_{\hat{y}=-1} intersect f_{\hat{z}=-1}
    void vertexPos(int i, int& sgnx, int& sgny, int& sgnz); // Fetch position information by local index in 1d array
    int vertexIndex(int sgnx, int sgny, int sgnz); // Fetch the local index in 1d array by position information
    int vertexGlobal(int sgnx, int sgny, int sgnz);
    int vertexGlobal(int i);
    Vertex* vertexPtr(int sgnx, int sgny, int sgnz);
    Vertex* vertexPtr(int i);

    // Edges
    // (sgnx, sgny, sgnz) = (0,-1,-1) -> f_{\hat{y}=-1} intersect f_{\hat{z}=-1}
    void edgePos(int i, int& sgnx, int& sgny, int& sgnz); // Fetch position information by local index in 1d array
    int edgeIndex(int sgnx, int sgny, int sgnz); // Fetch the local index in 1d array by position information
    Vertex* edgeVertexPtr(bool pm, int sgnx, int sgny, int sgnz); // pm = 0 or 1, means v0 and v1
    Vertex* edgeVertexPtr(bool pm, int i);
    int edgeGlobal(int sgnx, int sgny, int sgnz);
    int edgeGlobal(int i);
    void edgeFace(int nEdge, int& f0, int& f1);

    // Faces
    // (sgnx, sgny, sgnz) = (0,0,-1) -> f_{\hat{z}=-1}
    void facePos(int i, int& sgnx, int& sgny, int& sgnz); // Fetch position information by local index in 1d array
    int faceIndex(int sgnx, int sgny, int sgnz); // Fetch the local index in 1d array by position information
    Vertex* faceVertexPtr(bool pm0, bool pm1, int sgnx, int sgny, int sgnz); // For pm0 -> pm1, the priority follows xdim -> ydim -> zdim
    Vertex* faceVertexPtr(bool pm0, bool pm1, int i);
    int faceGlobal(int sgnx, int sgny, int sgnz);
    int faceGlobal(int i);

    // Partition the element into subtetrahedra
    std::vector<std::vector<int>>* subtetrahedra(int type);
    bool isInSubtetrahedron(const std::vector<int>& subtetrahedron, const Point& pt);
    int inSubtetrahedron(const std::vector<std::vector<int>>& subtetrahedra, const Point& pt);

    Tensor1 normal(int i); // normal vector of face i
    Tensor1 normal(int sgnx, int sgny, int sgnz) { return normal( faceIndex(sgnx, sgny, sgnz) ); };

    double lambda(int i, double x, double y, double z);
    double lambda(int sgnx, int sgny, int sgnz, double x, double y, double z) {
      return lambda(faceIndex(sgnx, sgny, sgnz), x,y,z); };
    double lambda(int i, const Point& p) { return lambda(i, p[0], p[1], p[2]); };
    double lambda(int sgnx, int sgny, int sgnz, const Point& p) { return lambda(sgnx, sgny, sgnz, p[0], p[1], p[2]); };

    Tensor1 dLambda(int i) { return -1*normal(i); };
    Tensor1 dLambda(int sgnx, int sgny, int sgnz) { return dLambda(faceIndex(sgnx, sgny, sgnz)); };
    
    Point forwardMap(const Point& p);
    Point backwardMap(const Point& p);
    void dForwardMap(const Point& p, Tensor2& df, double& jac);

    bool isInElement(const Point& pt);
    

    void write_raw(std::ofstream& fout);
    int write_raw(std::string& filename);

    friend class HexaMesh;
  };

 ////////////////////////////////////////////////////////////////////////////////
  // class HexaMesh
  //
  // Vertices: The position indices of the (nx+1) x (ny+1) x (nz+1) vertices
  //           are arranged in the order with priority z -> y -> x
  //
  //           (0,0,0) -> (0,0,1) -> (0,0,2) -> ... -> (0,0,nz)
  //        -> (0,1,0) -> (0,1,1) -> (0,1,2) -> ... -> (0,1,nz) -> ... -> (0,ny,nz)
  //        -> (1,0,0) -> (1,0,1) -> (1,0,2) -> ... -> (1,0,nz) 
  //        -> (1,1,0) -> (1,1,1) -> (1,1,2) -> ... -> (1,1,nz) -> ... -> (1,ny,nz) -> ... -> (nx,ny,nz)
  // 
  // Elements: The position indices of the nx x ny x nz vertices
  //           are arranged in the order with priority z -> y -> x
  //
  //           (0,0,0) -> (0,0,1) -> (0,0,2) -> ... -> (0,0,nz-1)
  //        -> (0,1,0) -> (0,1,1) -> (0,1,2) -> ... -> (0,1,nz-1) -> ... -> (0,ny-1,nz-1)
  //        -> (1,0,0) -> (1,0,1) -> (1,0,2) -> ... -> (1,0,nz-1) 
  //        -> (1,1,0) -> (1,1,1) -> (1,1,2) -> ... -> (1,1,nz-1) -> ... -> (1,ny-1,nz-1) -> ... -> (nx-1,ny-1,nz-1)
  //
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
  //
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
  //
  ////////////////////////////////////////////////////////////////////////////////

  class HexaMesh
  {
  private:
    int num_x, num_y, num_z;
    double min_x, max_x, min_y, max_y, min_z, max_z;
    
    int num_elements;
    int num_vertices;
    
    bool good_mesh;

    Vertex* the_vertices = nullptr;
    Element* the_elements = nullptr;

    Point* face_center = nullptr;
    double* face_radius = nullptr;


    // It returns the last point of a hexahedron, given all the other 7 points
    Point flatHex(const Point& v000, const Point& v001, 
                  const Point& v010, const Point& v011, 
                  const Point& v100, const Point& v101, const Point& v110);

    // It returns the center of the circle touching three edges v0-v1, v1-v2, v2-v3 in a face
    void potentialFaceCenter(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3, 
                        Point& center_candidate, double& radius, bool& inside);

    void faceCenter(int i, Point& center, double& radius);

    void set_hexamesh(int nx, int ny, int nz, Point* vertices);

  public:
    HexaMesh(int nx, int ny, int nz, Point* vertices) {
      set_hexamesh(nx,ny,nz,vertices); };
    HexaMesh() : num_x(1), num_y(1), num_z(1), 
                 min_x(0), max_x(1), min_y(0), max_y(1), min_z(0), max_z(1), 
                 num_elements(0), num_vertices(0), good_mesh(0) {};
    HexaMesh(Element* single_element); // Create a one element mesh from an element
    ~HexaMesh();

    void set(int nx, int ny, int nz, Point* vertices) {
      set_hexamesh(nx,ny,nz,vertices); };

    // Create a mesh with rectangular cuboidal domain
    int createMesh(char meshType, int nx, int ny, int nz, 
                      double xMin, double xMax, double yMin, double yMax, double zMin, double zMax, double distortionFactor=0);

    // Map the array index i to the position index (ix, iy, iz)
    void vertexPos(int i, int& ix, int& iy, int& iz);
    void elementPos(int i, int& ix, int& iy, int& iz);

    // Map the position index (ix, iy, iz) to the array index i
    int vertexIndex(int ix, int iy, int iz) const;
    int elementIndex(int ix, int iy, int iz) const;

    // Access functions
    int nElements() const { return num_elements; };
    int nVertices() const { return num_vertices; };
    int nEdges() const { return num_x*(num_y+1)*(num_z+1) + num_y*(num_x+1)*(num_z+1) + num_z*(num_x+1)*(num_y+1); };
    int nFaces() const { return (num_x+1)*num_y*num_z + (num_y+1)*num_x*num_z + (num_z+1)*num_x*num_y; };
    // Number of elements in each direction
    int xElements() const { return num_x; };
    int yElements() const { return num_y; };
    int zElements() const { return num_z; };
    double xMin() const { return min_x; };
    double yMin() const { return min_y; };
    double zMin() const { return min_z; };
    double xMax() const { return max_x; };
    double yMax() const { return max_y; };
    double zMax() const { return max_z; };
    bool isGood() const { return good_mesh; };
    
    // Access the element by its position in the array
    Element* elementPtr(int i) { return &the_elements[i % num_elements]; };
    // Access the element by its geometric position
    Element* elementPtr(int ix, int iy, int iz) { return elementPtr( elementIndex(ix,iy,iz) ); };

    // Access the vertex by its position in the array
    Vertex* vertexPtr(int i) { return &the_vertices[i % num_vertices]; };
    // Access the vertex by its geometric position
    Vertex* vertexPtr(int ix, int iy, int iz) { return vertexPtr( vertexIndex(ix,iy,iz) ); };

    // Access the vertex of face by its position in the array
    Vertex* faceVertexPtr(bool pm0, bool pm1, int i);
    // Access the center of face by its position in the array
    Point* faceCenterPtr(int i) { return &face_center[i]; };
    double faceRadius(int i) { return face_radius[i]; };

    // Return the element containing the point pt
    int inElement(const Point& pt);

    // Output
    void write_raw(std::ofstream& fout);
    int write_raw(std::string& filename);

    int write_matlab(std::string& filename) const;

    friend class Vertex;
    friend class Element;
  };
}
#endif