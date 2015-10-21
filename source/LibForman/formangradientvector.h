#ifndef FORMANGRADIENTVECTOR_H
#define FORMANGRADIENTVECTOR_H

#include "../LibMesh/Reader.h"
#include "forman_arrow.h"
#include "aux_forman_gradient_vector.h"
#include <map>
#include <set>
#include <list>
#include "assert.h"

#include "../LibMesh/Mesh.h"
#include "../LibMesh/Vertex3D.h"

#include "../Performances/Timer.h"
#include "../Performances/Usage.h"

#define NON_VALID_CONF 1000


typedef std::vector<Arrows>				CompressedToExpandedLUT;
typedef std::map<Arrows, unsigned int>	ExpandedToCompressedLUT;

//incidence graph used for representing the lower star of a vertex
struct Simplex_Graph{
    vector<int> simplex;
    vector<Simplex_Graph*> coboundary;
};



class QElem{
private:
    float edgHeight;
    int tIndex;
    short int edgeIndexOnT;

public:
    QElem(){edgHeight=0; tIndex=-1; edgeIndexOnT=-1;}
    QElem(float eh, int ti, short int ei) : edgHeight(eh), tIndex(ti), edgeIndexOnT(ei){ }
    QElem(QElem const& el) : edgHeight(el.getEdgeHeight()), tIndex(el.getTriangle()), edgeIndexOnT(el.getEdgeOnT()){ }

    inline float getEdgeHeight() const { return edgHeight;}
    inline int getTriangle() const{ return tIndex;}
    inline short int getEdgeOnT() const {return edgeIndexOnT;}

    inline bool operator<(QElem const& edge) const{
        return edgHeight < edge.getEdgeHeight();
    }
};

struct CompareQElem
{
    bool operator()(QElem const& lhs, QElem const& rhs) const
    {
        return lhs.getEdgeHeight() < rhs.getEdgeHeight();
    }
};

struct CompareQElemOtherWay
{
    bool operator()(QElem const& lhs, QElem const& rhs) const
    {
        return lhs.getEdgeHeight() > rhs.getEdgeHeight();
    }
};



//Class used for representing a Forman gradient on a triangle mesh
class FormanGradientVector
{

    //the triangle mesh
    Mesh* mesh;

    //compact representation for the Forman gradient
    vector<unsigned int> forman_gradient;

    //lookup tables
    CompressedToExpandedLUT compressed_to_expanded;
    ExpandedToCompressedLUT expanded_to_compressed;

    map<int,int> maxToTri;

public:

    ///constructor
    FormanGradientVector(Mesh* mesh);

    ///Functions for computing the gradient vector field
        //initialize the forman gradient by homotopy expansion
    void homotopy_expansion();
        //initialize the forman gradient by using a watershed segmentation
    void watershedToForman(vector<int>&);



    ///functions for critical simplexes
        //true if the vertex is critical (i.e. not paired with any incident edge)
    bool is_vertex_critical(int v);

        //true if the edge e=(v,v1) is critical (i.e. not paired with an incident vertex or triangle)
    bool is_edge_critical(int v, int v1);

        //true if the triangle f is critical (i.e. not paired with any incident edge)
    bool is_face_critical(int f);


    ///functions for initializing a gradient pair
        //define a new pair between an edge (v,v1) and a vertex v
    void setVE(int v, int v1);

        //define a new pair between a face (v,v1,v2)=f and an edge (v1,v2)
    void setEF(int v, int f);


    ///functions for retrieveing gradient pairs
        //return the edge paired with v (NULL if v is critical)
    Edge* getVE(int v);

        //return the edge paired with v (NULL if v is critical)
    int getEF(Edge*);

        //return position (inside triangle tri) of the edge paired with tri (-1 if tri is critical)
    int getFE(int tri);


    ///functions for deleting a pair from the gradient
        //unpair an edge (v1,v2) and its paired vertex v1
    void freeVE(int v1, int v2);

        //unpair a triangle (v,v1,v2)=f and its paired edge (v1,v2)
    void freeEF(int v, int f);


    //true if the gradient local frame represents a valid configuration (i.e. each simplex is paired with at most one other simplex)
    bool isValidTetCase(Arrows arrow);

    //return a vector containing the number of critical simplexes found in the gradient (from minima to maxima).
    vector<int> count_critical_simplexes();

    ///functions for computing the cells of the asceding and descending Morse complexes
    ///print in a vtk file the Morse cells to be visualized in Paraview
        //descending Morse complex (cells associated with maxima and saddles)
    vector<int> descending_2cells_extraction(bool withGeometry);
    void descending_1cells_extraction(bool withGeometry);
        //ascending Morse complex (cells associated with minima and saddles)
    vector<int> ascending_2cells_extraction(bool withGeometry);
    void ascending_1cells_extraction(bool withGeometry);

    ///output functions for printing the Morse cells
        //write the descending 2-cells
    void writeVTK_2cells(char* nome_file_output, vector<int> triangles);

        //write the descending 1-cells
    void writeVTK_1cells(char*, set<int> const&, set<pair<int,int> > const&);

        //write the ascending 2-cells
    void writeVTK_2cells_on_vert(char* nome_file_output, vector<int> const& triangles);

        //write the ascending 1-cells
    void writeVTK_1cells(char* file_name, set<pair<Vertex3D,Vertex3D> > const& edges);

        //write the Forman gradient
    void writeVTK_gradient(char* nomeFile);

        //write a point cloud indicating the critical simplexes
    void writeVTK_criticalPoints(char* nomeFile);

private:

    //converter from the compact representation of a set of gradient arrows to the expanded representation
    TriGradient convert_compressed_to_expand(unsigned int ga);

    //converter from the expanded representaiton to its compact one
    unsigned int convert_expand_to_compressed(Arrows ga);

    //extract the lower star of a vertex
    map<vector<int>, Simplex_Graph*> compute_lower_star(int v, vector<int>* vt);

    //compute the gradient via homotopic expansion
    ExplicitGradient* compute_gradient(vector<int>& v, map<vector<int>, Simplex_Graph*>* lower_link_structure, map<vector<int>, int>* crtical_points);
    void push_in_coboundary(Simplex_Graph* sface, Simplex_Graph* co_sface);
    int num_unpared_faces(vector<int>& co_face, ExplicitGradient* gradient_v, map<vector<int>, int>* critical_points);
    vector<int> unique_pairable_face(vector<int>& s_face, ExplicitGradient* gradient_v, map<vector<int>, int>* critical_points);

    //translate a triangle (represented as the triple of vertices) into its index
    inline int vector_to_index(vector<int> simplex){

        vector<int> vt = mesh->VT(simplex.front());
        for(unsigned int i=0; i<vt.size(); i++){

            int t = vt.at(i);
            int founded=0;

            for(int j=0; j<3; j++){

                int find = simplex.at(j);
                for(int k=0; k<3; k++){
                    if(find == mesh->getTopSimplex(t).TV(k)){
                        founded++;
                        break;
                    }
                }

                if(founded != j+1) break;
            }

            if(founded == 3) return t;
        }

        return -1;
    }

    //sort the simplex based on the function values of its vertices
    inline vector<int> sort_simplex(vector<int>* simplex){
        list<int> ordered;
        for(int i=0; i<simplex->size(); i++){
            if(i==0) ordered.push_back(simplex->at(0));
            else{
                list<int>::iterator it = ordered.begin();
                for(; it!=ordered.end(); it++){
                    if(mesh->getVertex(simplex->at(i)).getF() > mesh->getVertex(*it).getF())
                        break;
                }

                if(it == ordered.begin()) ordered.push_front(simplex->at(i));
                else if(it == ordered.end()) ordered.push_back(simplex->at(i));
                else ordered.insert(it, simplex->at(i));
            }
        }

        return vector<int>(ordered.begin(), ordered.end());
    }

    //compare two simplexes based on the function values of their vertices
    inline bool is_lower(const vector<int>& simplex1,const vector<int>& simplex2, bool* same){

        if(simplex1.size() < simplex2.size()) return true;
        if(simplex1.size() > simplex2.size()) return false;

        int similar_pair=-1;
        for(int i=1; i<simplex1.size(); i++){
            if(mesh->getVertex(simplex1.at(i)).getF() == mesh->getVertex(simplex2.at(i)).getF()){
                    if(simplex1[i]!=simplex2[i])similar_pair=i;
                    continue;
            }
            else if(mesh->getVertex(simplex1.at(i)).getF() > mesh->getVertex(simplex2.at(i)).getF()) return false;
            else return true;
        }

        if(simplex1.back() == simplex2.back() && similar_pair == -1)
            *same=true;
        else
            *same=false;

        return (simplex1[similar_pair]<simplex2[similar_pair]);
    }

    inline void push_ordered(list<vector<int> >* coda, vector<int>& simplesso){

        bool same=false;

        list<vector<int> >::iterator it = coda->begin();
        for(; it != coda->end(); it++){
            if(!is_lower(*it, simplesso, &same)) break;
        }

        if(!same){
            if(it == coda->begin()) coda->push_front(simplesso);
            else if(it == coda->end()) coda->push_back(simplesso);
            else coda->insert(it, simplesso);
        }
    }

    //extract the missing vertex comparing two incident simplexes
    inline int extract_missing_id(vector<int> &smaller, vector<int> &bigger)
    {
        int ret=-1;
        for(vector<int>::iterator it=bigger.begin(); it!=bigger.end(); ++it)
        {
            bool found=false;
            for(vector<int>::iterator it2=smaller.begin(); it2!=smaller.end(); ++it2)
            {
                if(*it==*it2)
                {
                    found = true;
                    break;
                }
            }
            if(!found)
            {
                ret = *it;
                break;
            }
        }
        return ret;
    }

    ///Functions for computing a Forman gradient starting from a watershed segmentation
    void vertexEdgePairings(vector<int> const& segmentation);
    void vertexEdgePairingsPlateau(vector<int> const& vertexLabels);
    void edgeTrianglePairings(vector<int> const& vertexLabels, int);
    void boundaryPairings(vector<int> const& vertexLabels);
    void removeTerrainBoundaryPathsFast(vector<bool> const& isMaximum);

    void edgeTriangleFast(vector<int> const& vertexLabels);

    int pick_a_minimum(set<int> const&);
    int pick_best_maximum(int i);
    int pick_best_maximum(set<int> const&);

    void reversePath(int criticalT, vector<int> const&);
    int computePairingPerRegion(list<int> trianglesList, vector<int> const&);
    int computePairingPerRegionAlt(list<int> trianglesList, vector<int> const&);
    int readSegmentation(vector<int>&);
    void follow_descendingPath(int tri);

};

#endif // FORMANGRADIENTVECTOR_H
