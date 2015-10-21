#ifndef EDGE_H
#define EDGE_H
#include "Vertex3D.h"
#include <algorithm>

using namespace std;

    // Class for represeting edges explicitly
    class Edge{
    private:
         int vertexes[2]; // Indices of the endpoints
    public:

        // Constructor
        inline Edge(){}
        inline Edge(int vi1, int vi2){
            vertexes[0]=min(vi1,vi2);
            vertexes[1]=max(vi1,vi2);
        }

        // Smaller-index endpoint getter
        inline int minindex() const{return vertexes[0];}

        // Bigger-index endpoint getter
        inline int maxindex() const{return vertexes[1];}

        inline int EV(int vPos) const
        {
            return vertexes[vPos];
        }

        //index of vertex v inside the edge
        inline short vertexindex(int v){
            if(vertexes[0]==v) return 0;
            if(vertexes[1]==v) return 1;
            return -1;
        }

        // Operators --------
        inline bool operator==(Edge edge){
            if(minindex()!=edge.minindex()) return false;
            if(maxindex()!=edge.maxindex()) return false;
            return true;
        }

        inline bool operator!=(Edge edge){
            return !((*this)==(edge));
        }

        inline bool operator<(Edge edge) const{
            if(maxindex() < edge.maxindex()) return true;
            else if(maxindex() == edge.maxindex() && minindex() < edge.minindex()) return true;
            else return false;
        }

        inline bool operator>(Edge edge) const{
            if(maxindex() > edge.maxindex()) return true;
            else if(maxindex() == edge.maxindex() && minindex() > edge.minindex()) return true;
            else return false;
        }


    };

#endif // EDGE_H
