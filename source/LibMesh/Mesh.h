#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <queue>

#include "Sorting.h"
#include "Edge.h"
#include "Triangle.h"
#include "map"
#include <set>
#include <list>
#include <assert.h>

/*
True if v is an endpoint of v1-v2
*/
#define is_endpoint(v,v1,v2) ( ( (v)==(v1) ) || ( (v)==(v2) ) )

using namespace std;

///A class representing a generic mesh parametrized by the type of top simplexes
class Mesh
{
protected:
    ///A private variable representing the vertices list of the mesh
    vector<Vertex3D> vertices;

    ///A private variable representing the top simplexes list of the mesh
    vector<Triangle> topSimplexes;

public:
    ///Contstructors
    Mesh();
    Mesh(const Mesh& orig);

    ///Descructors
    virtual ~Mesh();



    ///returns the number of vertices
    int getNumVertex();

    ///returns the number of top simplexes
    int getTopSimplexesNum();

    ///adds a vertex to the vertices list
    void addVertex(Vertex3D& v);

    ///returns the vertex at the i-th position in the mesh list
    Vertex3D& getVertex(int id);

    ///Adds a top simplex to the top simplexes list
    void addTopSimplex(Triangle& t);

    ///returns the triangle at the i-th position in the mesh list
    Triangle& getTopSimplex(int id);



    ///Initializes the space needed by the vertices and top simplexes lists
    void reserveVectorSpace(int numV, int numT);
    void reserveVectorSpace_TopSimplexes(int numT);
    void reserveVectorSpace_Vertices(int numV);

    ///builds the data structure
    void build();


    ///TOPOLOGICAL RELATIONS
    //retrieves simplexes incident into a vertex
    vector<Edge*> VE(int center);
    vector<int> VT(int center);
    //retrieves vertices adjacent to the given vertex
    vector<int> VV(int center);
    //true if the vertex is on the boundary of the mesh
    bool isBoundary(int center);

    //retrieves triangles incident into an edge
    vector<int> ET(Edge e);
    //retrieves edges adjacent to a given edge
    vector<Edge*> EE(Edge e);


    ///OUTPUT FUNCTIONS
    //output the mesh into VTK format with a given segmentation defined on its vertices or triangles (optional)
    void writeVTK(const char* filename, vector<int>* segm=NULL);

    //output the mesh into TRI format
    void writeTRI(const char* filename);

    //output the mesh into OFF format
    void writeOFF(const char* filename);


    ///FUNCTIONS FOR MODIFYING THE SCALAR FIELD
    void laplacianSmoothing(int iter);

protected:



    void find_adjacencies();
    void find_incidencies();
    void link_adj (int t1, int t2);
    int valid_vertex(int v);
    int valid_triangle(int t);

    template<class C> bool isIntoVector(C c, vector<C> &c_vect);
};


#endif // MESH_H
