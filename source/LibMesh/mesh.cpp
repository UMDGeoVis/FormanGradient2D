#include "Mesh.h"


Mesh::Mesh() {
    vertices = vector<Vertex3D>();
    topSimplexes = vector<Triangle>();
}


Mesh::Mesh(const Mesh& orig) {
    this->vertices = orig.vertices;
    this->topSimplexes = orig.topSimplexes;
}


Mesh::~Mesh() {
    vertices.clear();
    topSimplexes.clear();
}


Vertex3D& Mesh::getVertex(int id){
    return this->vertices[id];
}


int Mesh::getNumVertex(){
    return this->vertices.size();
}


void Mesh::addVertex(Vertex3D& v){
    this->vertices.push_back(v);
}


Triangle& Mesh::getTopSimplex(int id){
    return this->topSimplexes.at(id);
}


int Mesh::getTopSimplexesNum(){
    return this->topSimplexes.size();
}


void Mesh::reserveVectorSpace(int numV, int numT){
    this->vertices.reserve(numV);
    this->topSimplexes.reserve(numT);
}


void Mesh::reserveVectorSpace_Vertices(int numV)
{
    this->vertices.reserve(numV);
}


void Mesh::reserveVectorSpace_TopSimplexes(int numT)
{
    this->topSimplexes.reserve(numT);
}


void Mesh::addTopSimplex(Triangle& t){
    this->topSimplexes.push_back(t);
}


void Mesh::build()
{
    this->find_adjacencies();
    this->find_incidencies();
}


void Mesh::find_adjacencies()
{
    aux *tr_vec ;
    int i, j, k;
    int t1, t2;
    int v1, v2;

    tr_vec = (aux*) calloc( this->getTopSimplexesNum()*3 , sizeof(aux) ) ;
    if (!tr_vec)
    {
        cerr << "malloc fallita in find_adj" <<endl;
        return;
    }

    k = 0;

    for (j=0; j<this->getTopSimplexesNum(); j++)
    {
        for (i=0;i<3;i++)
        {
            tr_vec[k].t = j;
            this->getTopSimplex(j).setTT(i,-1);
            v1 = this->getTopSimplex(j).TV(i);
            v2 = this->getTopSimplex(j).TV((i+1)%3);
            if (v1<v2) {  tr_vec[k].v1 = v1; tr_vec[k].v2 = v2;  }
            else {  tr_vec[k].v1 = v2; tr_vec[k].v2 = v1;  }
            k++;
        }
    }

    qsort(tr_vec,3*this->getTopSimplexesNum(),sizeof(aux),cmp_aux) ;

    for(k=0;k<3*this->getTopSimplexesNum()-1;k++)
    {
        if ( is_endpoint(tr_vec[k].v1,tr_vec[k+1].v1,tr_vec[k+1].v2) &&
             is_endpoint(tr_vec[k].v2,tr_vec[k+1].v1,tr_vec[k+1].v2) )
        {
            t1 = tr_vec[k].t;
            t2 = tr_vec[k+1].t;
            link_adj(t1,t2);
        }
    }
    free(tr_vec) ;
    return;
}


void Mesh::find_incidencies()
{
    int i, t;

    for (t=0;t<this->getTopSimplexesNum();t++)
    {
        for (i=0;i<3;i++)
        {
            if(this->getVertex(this->getTopSimplex(t).TV(i)).VTstar()==-1)
                this->getVertex(this->getTopSimplex(t).TV(i)).VTstar(t);
            else if(this->getTopSimplex(t).TT((i+2)%3)==-1)
                this->getVertex(this->getTopSimplex(t).TV(i)).VTstar(t);
        }
    }
    return;
}


void Mesh::link_adj(int t1, int t2)
{
    int i, j, k, pos1[2], pos2[2];
    if (valid_triangle(t1) && valid_triangle(t2))
    {
        k = 0;
        for (i=0; ((i<3)&&(k<2)); i++)
        {
            for (j=0; ((j<3)&&(k<2)); j++)
            {
                if(this->getTopSimplex(t1).TV(i) == this->getTopSimplex(t2).TV(j))
                {
                    pos1[k] = i;
                    pos2[k] = j;
                    k++;
                }
            }
        }
        if (k<2)
        {
            cerr << "error in link_adj" <<endl;
        }
        else
        {
            this->getTopSimplex(t1).setTT(3-pos1[0]-pos1[1],t2);
            this->getTopSimplex(t2).setTT(3-pos2[0]-pos2[1],t1);
        }
    }
    else
        cerr << "Error in link_adj: almeno uno dei triangoli e' nullo" << endl;
}


int Mesh::valid_vertex(int v)
{
    return ( (v>=0) && (v<this->getNumVertex()) );
}


int Mesh::valid_triangle(int t)
{
    return ( (t>=0) && (t<this->getTopSimplexesNum()) );
}


vector<int> Mesh::VT(int center)
{
    vector<int> triangles;
    int pred = -1;
    int current = this->getVertex(center).VTstar();
    triangles.push_back(this->getVertex(center).VTstar());

    int k=-1;

    for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
    {
        if(this->getTopSimplex(current).TV(i) == center)
        {
            k = i;
            break;
        }
    }

    pred = current;
    current = this->getTopSimplex(current).TT((k+1)%3);

    bool isBorder = false;
    while(1)
    {
        if(current == this->getVertex(center).VTstar())
            break;
        else if(current == -1)
        {
            isBorder = true;
            break;
        }
        else
            triangles.push_back(current);

        k=-1;
        for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
        {
            if(this->getTopSimplex(current).TV(i) == center)
            {
                k = i;
                break;
            }
        }

        if(this->getTopSimplex(current).TT((k+1)%3) == pred)
        {
            pred = current;
            current = this->getTopSimplex(current).TT((k+2)%3);
        }
        else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
        {
            pred = current;
            current = this->getTopSimplex(current).TT((k+1)%3);
        }
    }

    if(isBorder)
    {
        pred = this->getVertex(center).VTstar();
        for(int i=0;i<this->getTopSimplex(pred).getVerticesNum();i++)
        {
            if(this->getTopSimplex(pred).TV(i) == center)
            {
                k = i;
                break;
            }
        }
        current = this->getTopSimplex(pred).TT((k+2)%3);

        while(1)
        {
            if(current == -1)
                break;
            else
                triangles.push_back(current);

            k=-1;
            for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
            {
                if(this->getTopSimplex(current).TV(i) == center)
                {
                    k = i;
                    break;
                }
            }


            if(this->getTopSimplex(current).TT((k+1)%3) == pred)
            {
                pred = current;
                current = this->getTopSimplex(current).TT((k+2)%3);
            }
            else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
            {
                pred = current;
                current = this->getTopSimplex(current).TT((k+1)%3);
            }
        }
    }

    return triangles;
}



bool Mesh::isBoundary(int center)
{
    vector<int> triangles;
    int pred = -1;
    int current = this->getVertex(center).VTstar();
    triangles.push_back(this->getVertex(center).VTstar());

    int k=-1;

    for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
    {
        if(this->getTopSimplex(current).TV(i) == center)
        {
            k = i;
            break;
        }
    }

    pred = current;
    current = this->getTopSimplex(current).TT((k+1)%3);

    bool isBorder = false;
    while(1)
    {
        if(current == this->getVertex(center).VTstar())
            break;
        else if(current == -1)
        {
            return true;
        }
        else
            triangles.push_back(current);

        k=-1;
        for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
        {
            if(this->getTopSimplex(current).TV(i) == center)
            {
                k = i;
                break;
            }
        }

        if(this->getTopSimplex(current).TT((k+1)%3) == pred)
        {
            pred = current;
            current = this->getTopSimplex(current).TT((k+2)%3);
        }
        else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
        {
            pred = current;
            current = this->getTopSimplex(current).TT((k+1)%3);
        }
    }

    return false;
}


vector<Edge*> Mesh::VE(int center)
{
    vector<Edge*> edges;
    int pred = -1;
    int current = this->getVertex(center).VTstar();

    int k=-1;
    for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
    {
        if(this->getTopSimplex(current).TV(i) == center)
        {
            k = i;
            break;
        }
    }
    edges.push_back(this->getTopSimplex(current).TE((k+1)%3));

    pred = current;
    current = this->getTopSimplex(current).TT((k+1)%3);

    bool isBorder = false;

    while(1)
    {
        if(current == this->getVertex(center).VTstar())
            break;
        else if(current == -1)
        {
            isBorder = true;
            break;
        }

        k=-1;
        for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
        {
            if(this->getTopSimplex(current).TV(i) == center)
            {
                k = i;
                break;
            }
        }

        if(this->getTopSimplex(current).TT((k+1)%3) == pred)
        {
            edges.push_back(this->getTopSimplex(current).TE((k+2)%3));
            pred = current;
            current = this->getTopSimplex(current).TT((k+2)%3);
        }
        else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
        {
            edges.push_back(this->getTopSimplex(current).TE((k+1)%3));
            pred = current;
            current = this->getTopSimplex(current).TT((k+1)%3);
        }
    }

    if(isBorder)
    {
        pred = this->getVertex(center).VTstar();
        for(int i=0;i<this->getTopSimplex(pred).getVerticesNum();i++)
        {
            if(this->getTopSimplex(pred).TV(i) == center)
            {
                k = i;
                break;
            }
        }
        edges.push_back(this->getTopSimplex(pred).TE((k+2)%3));
        current = this->getTopSimplex(pred).TT((k+2)%3);

        while(1)
        {
            if(current == -1)
                break;

            k=-1;
            for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
            {
                if(this->getTopSimplex(current).TV(i) == center)
                {
                    k = i;
                    break;
                }
            }

            if(this->getTopSimplex(current).TT((k+1)%3) == pred)
            {
                edges.push_back(this->getTopSimplex(current).TE((k+2)%3));
                pred = current;
                current = this->getTopSimplex(current).TT((k+2)%3);
            }
            else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
            {
                edges.push_back(this->getTopSimplex(current).TE((k+1)%3));
                pred = current;
                current = this->getTopSimplex(current).TT((k+1)%3);
            }
        }
    }

    return edges;
}


vector<int> Mesh::VV(int center)
{
    vector<int> vertices;
    int pred = -1;
    int current = this->getVertex(center).VTstar();

    int k=-1;
    for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
    {
        if(this->getTopSimplex(current).TV(i) == center)
        {
            k = i;
            break;
        }
    }
    vertices.push_back(this->getTopSimplex(current).TV((k+1)%3));

    pred = current;
    current = this->getTopSimplex(current).TT((k+2)%3);

    bool isBorder = false;

    while(1)
    {
        if(current == this->getVertex(center).VTstar())
            break;
        else if(current == -1)
        {
            isBorder = true;
            break;
        }

        k=-1;
        for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
        {
            if(this->getTopSimplex(current).TV(i) == center)
            {
                k = i;
                break;
            }
        }

        if(this->getTopSimplex(current).TT((k+1)%3) == pred)
        {
            vertices.push_back(this->getTopSimplex(current).TV((k+1)%3));
            pred = current;
            current = this->getTopSimplex(current).TT((k+2)%3);
        }
        else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
        {
            vertices.push_back(this->getTopSimplex(current).TV((k+2)%3));
            pred = current;
            current = this->getTopSimplex(current).TT((k+1)%3);
        }
    }

    if(isBorder)
    {
        pred = this->getVertex(center).VTstar();
        for(int i=0;i<this->getTopSimplex(pred).getVerticesNum();i++)
        {
            if(this->getTopSimplex(pred).TV(i) == center)
            {
                k = i;
                break;
            }
        }
        vertices.push_back(this->getTopSimplex(pred).TV((k+2)%3));
        current = this->getTopSimplex(pred).TT((k+1)%3);

        while(1)
        {
            if(current == -1)
                break;

            k=-1;
            for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
            {
                if(this->getTopSimplex(current).TV(i) == center)
                {
                    k = i;
                    break;
                }
            }

            if(this->getTopSimplex(current).TT((k+1)%3) == pred)
            {
                vertices.push_back(this->getTopSimplex(current).TV((k+1)%3));
                pred = current;
                current = this->getTopSimplex(current).TT((k+2)%3);
            }
            else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
            {
                vertices.push_back(this->getTopSimplex(current).TV((k+2)%3));
                pred = current;
                current = this->getTopSimplex(current).TT((k+1)%3);
            }
        }
    }

    return vertices;
}


vector<int> Mesh::ET(Edge e)
{
    vector<int> triangles;
    vector<int> vtcomplete_triangles = this->VT(e.EV(0));
    for(unsigned int i=0;i<vtcomplete_triangles.size();i++)
    {
        for(int j=0;j<this->getTopSimplex(vtcomplete_triangles.at(i)).getVerticesNum();j++)
        {
            if(e.EV(1) == this->getTopSimplex(vtcomplete_triangles.at(i)).TV(j))
            {
                triangles.push_back(vtcomplete_triangles.at(i));
                break;
            }
        }
    }

    return triangles;
}


vector<Edge*> Mesh::EE(Edge e)
{
    vector<Edge*> edges;
    vector<Edge*> edge0 = this->VE(e.EV(0));
    vector<Edge*> edge1 = this->VE(e.EV(1));
    for(unsigned int i=0;i<edge0.size();i++)
    {
        if(e != *edge0.at(i))
            edges.push_back(edge0.at(i));
    }
    for(unsigned int i=0;i<edge1.size();i++)
    {
        if(e != *edge1.at(i))
            edges.push_back(edge1.at(i));
    }

    return edges;
}


template<class C> bool Mesh :: isIntoVector(C c, vector<C> &c_vect)
{
    for(unsigned int i=0; i<c_vect.size();i++)
        if(c == c_vect.at(i))
            return true;
    return false;
}


struct Comparer{
    inline bool operator()(pair<int,float> v1, pair<int, float> v2){
        return v1.second > v2.second;
    }
};



void Mesh::writeVTK(const char* filename,vector<int>* segm){

    FILE* file;
    file = fopen(filename, "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", getNumVertex());

    for(int i=0; i<getNumVertex(); i++)
        fprintf(file, "%f %f %f\n", getVertex(i).getX(), getVertex(i).getY(), getVertex(i).getZ());
    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", getTopSimplexesNum(), getTopSimplexesNum()*4);

    for(int i=0; i<getTopSimplexesNum(); i++)
        fprintf(file, "3 %d %d %d\n", getTopSimplex(i).TV(0), getTopSimplex(i).TV(1), getTopSimplex(i).TV(2));
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", getTopSimplexesNum());

    for(int i=0; i<getTopSimplexesNum(); i++)
        fprintf(file, "%d ", 5);
    fprintf(file, "\n\n");


    if(segm!=NULL)
    {
        if(segm->size() == getNumVertex()){
            fprintf(file, "POINT_DATA %d \n", getNumVertex());
        }
        else{
            assert(segm->size() == getTopSimplexesNum());
            fprintf(file, "CELL_DATA %d \n", getTopSimplexesNum());
        }

        fprintf(file, "FIELD FieldData 1\n");
        fprintf(file, "segmentation 1 %d int\n", segm->size());

        for(int i=0; i<segm->size(); i++)
            fprintf(file, "%d ", (*segm)[i]);

        fprintf(file, "\n\n");

    }

    fclose(file);

}




void Mesh::writeTRI(const char* filename){

    FILE* file;
    file = fopen(filename, "w");

    fprintf(file, "%d\n", getNumVertex());

    for(int i=0; i<getNumVertex(); i++)
        fprintf(file, "%f %f %f\n", getVertex(i).getX(), getVertex(i).getY(), getVertex(i).getZ());
    fprintf(file, "\n\n");

    fprintf(file, "%d\n", getTopSimplexesNum());

    for(int i=0; i<getTopSimplexesNum(); i++)
        fprintf(file, "%d %d %d\n", getTopSimplex(i).TV(0), getTopSimplex(i).TV(1), getTopSimplex(i).TV(2));
    fprintf(file, "\n");

    fclose(file);

}


void Mesh::writeOFF(const char* filename){

    FILE* file;
    file = fopen(filename, "w");

    fprintf(file, "OFF\n%d %d %D\n", getNumVertex(), getTopSimplexesNum(), 0);

    for(int i=0; i<getNumVertex(); i++)
        fprintf(file, "%f %f %f\n", getVertex(i).getX(), getVertex(i).getY(), getVertex(i).getZ());
    fprintf(file, "\n\n");

    for(int i=0; i<getTopSimplexesNum(); i++)
        fprintf(file, "3 %d %d %d\n", getTopSimplex(i).TV(0), getTopSimplex(i).TV(1), getTopSimplex(i).TV(2));
    fprintf(file, "\n");

    fclose(file);

}
