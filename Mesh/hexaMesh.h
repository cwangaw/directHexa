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
    HexaMesh* my_mesh;
    int my_mesh_index;

    void set_vertex(double p0, double p1, double p2, int myIndex=-1, HexaMesh* myMesh=nullptr) {
      the_point[0] = p0; the_point[1] = p1; the_point[2] = p2; my_mesh_index = myIndex; my_mesh = myMesh; };

  public:
    Vertex() { set_vertex(0,0,0); };
    Vertex(double p0, double p1, double p2, int myIndex=-1, HexaMesh* myMesh=nullptr) {
      set_vertex(p0, p1, p2, myIndex, myMesh); };
    Vertex(const Point& p, int myIndex=-1, HexaMesh* myMesh=nullptr) {
      set_vertex(p[0], p[1], p[2], myIndex, myMesh); };
    Vertex(const Point* p, int myIndex=-1, HexaMesh* myMesh=nullptr) {
      set_vertex((*p)[0], (*p)[1], (*p)[2], myIndex, myMesh); };
    
    void set() { set_vertex(0,0,0); };
    void set(double p0=0, double p1=0, double p2=0, int myIndex=-1, HexaMesh* myMesh=nullptr) {
      set_vertex(p0, p1, p2, myIndex, myMesh); };
    void set(const Point& p, int myIndex=-1, HexaMesh* myMesh=nullptr) {
      set_vertex(p[0], p[1], p[2], myIndex, myMesh); };
    void set(const Point* p, int myIndex=-1, HexaMesh* myMesh=nullptr) {
      set_vertex((*p)[0], (*p)[1], (*p)[2], myIndex, myMesh); };

    HexaMesh* mesh() const { return my_mesh; };
    int meshIndex() const { return my_mesh_index; };
     
    void nbrElements(std::vector<int>& theNbrIndices) const;
    void nbrEdges(std::vector<int>& theNbrIndices) const;
    bool isOnBoundary() const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };


  class Element
  {
  private:
    Vertex* the_vertices;
    bool good_element;
    double the_volume;
    Point my_center;  // of largest inscribed ball
    Point my_centroid; // center of mass
    double my_diameter;
    double max_radius;

    HexaMesh* my_mesh;
    int my_mesh_index;
    
    bool isConnected();
    bool isConvexCounterclockwise(); //?
    double computeAreaTetra(const Point* p0, const Point* p1, const Point* p2, const Point* p3);
    double computeArea();
    void computeCenterOfTetra(const Point* p0, const Point* p1, const Point* p2, const Point* p3, Point& center, bool& unique);

    void set_element(Vertex** theVertices, int myIndex=-1, HexaMesh* myMesh=nullptr);

  public:
    Element(Vertex** theVertices, int myIndex=-1, HexaMesh* myMesh=nullptr) {
      set_element(theVertices, myIndex, myMesh); };
    Element() : good_element(false),
		    the_volume(0), my_center(), my_centroid(), my_diameter(0),
		    my_mesh(nullptr), my_mesh_index(-1) {};
    ~Element();

    void set(Vertex** theVertices, int myIndex=-1, HexaMesh* myMesh=nullptr) {
      set_element(theVertices, myIndex, myMesh); };
    
    bool isGood() const { return good_element; };

    Vertex vertexPtr(int i) { return the_vertices[i % 8]; };

    HexaMesh* mesh() const { return my_mesh; };
    int meshIndex() const { return my_mesh_index; };
    
    double volume() const { return the_volume; };
    Point center()   const { return my_center; };
    Point centroid() const { return my_centroid; };
    double maxRadius() const { return max_radius; }
    double diameter() const { return my_diameter; };
    double chunkParam();

    bool isInElement(const Point& pt) const;
    bool isOnElementBoundary(const Point& pt) const;

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    friend class HexaMesh;
  };
}
#endif