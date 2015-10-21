#ifndef Triangle_H
#define Triangle_H

#include "Edge.h"

class Triangle
{
public:

    ///Constructors
    Triangle();
    Triangle(int v1, int v2, int v3);

    ///return vertex index on position pos
    int TV(int pos);

    ///return Edge in position pos
    Edge* TE(int pos);

    ///return adjacent triangle
    int TT(int pos);

    ///setter for adjacencies
    void setTT(int pos, int adjId);

    ///retrieve the position of an edge inside the triangle
    short face_index(Edge* e);

    ///retrieve the position of a vertex inside the triangle
    int vertex_index(int v);

    ///operators -----
    inline bool operator== (const Triangle &p)
    {
        bool b[3];
        b[0] = false; b[1] = false; b[2] = false;
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(!b[j] && p.vertices[i]==this->vertices[j])
                {
                    b[j] = true;
                    break;
                }
            }
        }

        return b[0] && b[1] && b[2];
    }

    inline bool operator!= (const Triangle p)
    {
       return !((*this)==(p));
    }

    ///true if it contains the given vertex
    inline bool contains(int v){
        for(int i=0;i<3;i++)
        {
            if(vertices[i] == v) return true;
        }
        return false;
    }

    inline int getVerticesNum(){ return 3;}

private:
    int vertices[3]; /// indexes of its vertices
    int adj[3]; ///indexes of its adjacent triangles
};



#endif // Triangle_H
