#include "Triangle.h"

Triangle::Triangle()
{

    this->vertices[0] = -1;
    this->vertices[1] = -1;
    this->vertices[2] = -1;
}

Triangle::Triangle(int v1, int v2, int v3)
{
    this->vertices[0] = v1;
    this->vertices[1] = v2;
    this->vertices[2] = v3;
}

int Triangle::TV(int pos)
{
    return this->vertices[pos];
}

Edge* Triangle::TE(int pos)
{
    return new Edge(vertices[(pos+1)%3],vertices[(pos+2)%3]);
}

int Triangle::TT(int pos)
{
    return this->adj[pos];
}

void Triangle::setTT(int pos, int adjId)
{
    this->adj[pos]=adjId;
}

short Triangle::face_index(Edge* e){
    for(int i=0; i<3; i++){
        if(vertices[i] != e->EV(0) && vertices[i] != e->EV(1))
            return i;
    }

    return -1;
}


int Triangle::vertex_index(int v){
    for(int i=0; i<3; i++){
        if(vertices[i] == v)
            return i;
    }

    return -1;
}
