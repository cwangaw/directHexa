#ifndef __directserendipity_h_included__
#define __directserendipity_h_included__

////////////////////////////////////////////////////////////////////////////////
// DirectSerendipity class
//   Includes classes for DirectSerendipityFE
// Direct Serendipity Elements on convex polygons
//
// Uses the PolyMesh classes for the mesh and elements
// Assume base objects: Point and Tensor1 (both two numbers)
////////////////////////////////////////////////////////////////////////////////

#include "Mesh/baseObjects.h"
#include "Mesh/hexaMesh.h"
#include "fcns.h"

namespace directserendipity {

  class DirectSerendipity;

  ////////////////////////////////////////////////////////////////////////////////
  // class DirectSerendipityArray
  //    Gives values for each nodal point
  //    Also gives the gradient
  ////////////////////////////////////////////////////////////////////////////////
  
  class DirectSerendipityArray
  {
  private:
    int num_dofs;
    double* the_array = nullptr;

    DirectSerendipity* my_ds_space;

    void set_directserendipityarray(DirectSerendipity* dsSpace);

  public:
    DirectSerendipityArray() : num_dofs(0), the_array(nullptr), my_ds_space(nullptr) {};
    DirectSerendipityArray(DirectSerendipity* dsSpace) {
      set_directserendipityarray(dsSpace); };
    DirectSerendipityArray(const DirectSerendipityArray& a) : the_array(nullptr) {
      set_directserendipityarray(a.dsSpace()); };
    ~DirectSerendipityArray();
    
    void set(DirectSerendipity* dsSpace) { set_directserendipityarray(dsSpace); };

    DirectSerendipity* dsSpace() const { return my_ds_space; };
    int size() const { return num_dofs; };

    double& operator() (int i)       { return the_array[i]; }
    double  operator() (int i) const { return the_array[i]; }
    double& operator[] (int i)       { return the_array[i]; }
    double  operator[] (int i) const { return the_array[i]; }

    void eval(const Point* pts, double* result, Tensor1* gradResult, int num_pts) const;
    void eval(const Point& pt, double& result, Tensor1& gradResult) const;
    double eval(const Point& pt) const;

    void l2normError(double& l2Error, double& l2GradError, double& l2Norm, double& l2GradNorm, int refinement_level = 0,
		     double (*referenceFcn)(double,double,double) = nullptr, 
		     Tensor1 (*referenceGradFcn)(double,double,double) = nullptr);
    void l2norm(double& l2Norm, double& l2GradNorm, int refinement_level = 0) {
      double null1,null2; l2normError(l2Norm,l2GradNorm,null1,null2,refinement_level); };

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  // class DirectSerendipityFE
  //    Defined on an element (class Element)
  //    
  //    Gives basis functions
  //      First call initBasis to evaluate all basis functions corresponding to DoFs
  //      Then access the basis functions by node number, or type and number, and pt number
  //
  // Store the result in such a way that the evaluations of all basis functions at
  // one point are put together.
  ////////////////////////////////////////////////////////////////////////////////

 class DirectSerendipityFE 
  {
  private:
    DirectSerendipity* my_ds_space;
    hexamesh::Element* my_element;

    int polynomial_degree; // Redundant with my_ds_space
    int deg_face_poly; // Redundant with above
    int num_face_dofs; // Redundant with above
    int deg_cell_poly; // Redundant with above
    int num_cell_dofs; // Redundant with above
    int num_dofs; // Redundant with above

    // Evaluation storage
    int num_eval_pts;
    double* value_n = nullptr;
    Tensor1* gradvalue_n = nullptr;

    // If necessary (degPolyn small), the bigger space within which we construct the basis
    hexamesh::HexaMesh* one_element_mesh = nullptr;
    DirectSerendipity* high_order_ds_space = nullptr;

    void set_directserendipityfe(DirectSerendipity* dsSpace, hexamesh::Element* element);

  public:
    DirectSerendipityFE() : my_ds_space(nullptr), my_element(nullptr), 
			    polynomial_degree(0), deg_cell_poly(-1), num_cell_dofs(0), num_dofs(0) {};
    DirectSerendipityFE(DirectSerendipity* dsSpace, hexamesh::Element* element) {
      set_directserendipityfe(dsSpace, element); };
    ~DirectSerendipityFE();
      
    void set(DirectSerendipity* dsSpace, hexamesh::Element* element) {
      set_directserendipityfe(dsSpace, element); };

    // Access functions
    DirectSerendipity* dsSpace() const { return my_ds_space; };
    hexamesh::Element* elementPtr() const { return my_element; };

    int degPolyn() const { return polynomial_degree; };
    int degFacePolyn() const { return deg_face_poly; };
    int degCellPolyn() const { return deg_cell_poly; };
    
    int nVertexDoFs() const { return 8; };
    int nEdgeDoFs() const { return 12*(polynomial_degree-1); };
    int nFaceDoFs() const { return num_face_dofs; };
    int nCellDoFs() const { return num_cell_dofs; };
    int nDoFs() const { return num_dofs; };

    // FE basis functions

    void initBasis(const Point* pt, int num_pts); // (re)evaluate all basis fcns at the points

    double basis(int iNode, int iPt) const { return value_n[iNode + iPt*num_dofs]; };
    Tensor1 basisGrad(int iNode, int iPt) const { return gradvalue_n[iNode + iPt*num_dofs]; };
    
    double vertexBasis(int iVNode, int iPt) const { return value_n[iVNode + iPt*num_dofs]; };
    Tensor1 gradVertexBasis(int iVNode, int iPt) const { return gradvalue_n[iVNode + iPt*num_dofs]; };

    double edgeBasis(int iENode, int iPt) const {
      return value_n[iENode + nVertexDoFs() + iPt*num_dofs]; };
    Tensor1 gradEdgeBasis(int iENode, int iPt) const { return gradvalue_n[iENode + nVertexDoFs() + iPt*num_dofs]; }
    
    double faceBasis(int iFNode, int iPt) const {
      return value_n[iFNode + nVertexDoFs() + nEdgeDoFs() + iPt*num_dofs]; };
    Tensor1 gradFaceBasis(int iFNode, int iPt) const {
      return gradvalue_n[iFNode + nVertexDoFs() + nEdgeDoFs() + iPt*num_dofs]; };

    double cellBasis(int iCNode, int iPt) const {
      return value_n[iCNode + nVertexDoFs() + nEdgeDoFs() + nFaceDoFs() + iPt*num_dofs]; };
    Tensor1 gradCellBasis(int iCNode, int iPt) const {
      return gradvalue_n[iCNode + nVertexDoFs() + nEdgeDoFs() + nFaceDoFs() + iPt*num_dofs]; };  

    void eval(const Point* pt, double* result, Tensor1* gradResult, int num_pts,
	      double* vertex_dofs, double* edge_dofs=nullptr, double* face_dofs=nullptr, double* cell_dofs=nullptr);
    void eval(const Point& pt, double& result, Tensor1& gradResult, double* vertex_dofs,
              double* edge_dofs=nullptr, double* face_dofs=nullptr, double* cell_dofs=nullptr) {
      eval(&pt, &result, &gradResult, 1, vertex_dofs, edge_dofs, face_dofs, cell_dofs); };

    void eval(const Point* pt, double* result, int num_pts, double* vertex_dofs,
              double* edge_dofs=nullptr, double* face_dofs=nullptr, double* cell_dofs=nullptr);
    double eval(const Point& pt, double* vertex_dofs,
		double* edge_dofs=nullptr, double* face_dofs=nullptr, double* cell_dofs=nullptr) { double result;
      eval(&pt, &result, 1, vertex_dofs, edge_dofs, face_dofs, cell_dofs); return result; }

    // Output functions
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };
  ////////////////////////////////////////////////////////////////////////////////
  // class DirectSerendipity
  //   Defined on a logically rectangular mesh (class HexaMesh), consisting of direct serendipity
  //     finite elements (class DirectSerendipityFE)
  //    
  //   Gives function evaluation from global DoFs for the basis functions 
  //
  //   Vertex DoFs: evaluation
  //
  ////////////////////////////////////////////////////////////////////////////////

  
  class DirectSerendipity
  {
  private:
    int polynomial_degree;
    int supplement_type;
    hexamesh::HexaMesh* my_mesh;
    //DirectSerendipity* my_high_order_ds_space;

    int num_dofs;
    int num_dofs_per_face;
    int num_dofs_per_cell;

    Point* face_dofs = nullptr;
  //  Point* cell_dofs = nullptr;
    DirectSerendipityFE* the_ds_elements = nullptr;

    void set_directserendipity(int polyDeg, int suppType, hexamesh::HexaMesh* mesh);
    
  public:
    DirectSerendipity() : polynomial_degree(0), supplement_type(0), my_mesh(nullptr), num_dofs(0), num_dofs_per_face(0), num_dofs_per_cell(0), 
			  face_dofs(nullptr), the_ds_elements(nullptr) {};
    DirectSerendipity(int polyDeg, int suppType, hexamesh::HexaMesh* mesh) {
      set_directserendipity(polyDeg, suppType, mesh); };
    ~DirectSerendipity();

    void set(int polyDeg, int suppType, hexamesh::HexaMesh* mesh) { set_directserendipity(polyDeg, suppType, mesh); };
    
    int degPolyn() const { return polynomial_degree; };
    int supplementType() const { return supplement_type; }
    hexamesh::HexaMesh* mesh() const { return my_mesh; };
    
    int nDoFs() const { return num_dofs; };
    int nVertexDoFs() const { return my_mesh->nVertices(); };
    int nEdgeDoFs() const { return my_mesh->nEdges() * (polynomial_degree-1); };
    int nSingleFaceDoFs() const { return num_dofs_per_face; };
    int nFaceDoFs() const { return my_mesh->nFaces() * nSingleFaceDoFs(); };
    int nSingleCellDoFs() const { return num_dofs_per_cell; };
    int nCellDoFs() const { return my_mesh->nElements() * nSingleCellDoFs(); };

    DirectSerendipityFE* finiteElementPtr(int i) const { return &the_ds_elements[i]; }

    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;

    int write_matlab(std::string& filename) const;

    friend class DirectSerendipityFE;
    friend class DirectSerendipityArray;
  };
 };

#endif