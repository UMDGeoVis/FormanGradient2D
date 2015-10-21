#include "Reader.h"
#include <sstream>
#include <algorithm>
#include <iostream>
#include <fstream>

Reader::Reader() {}

Reader::Reader(const Reader& ) {}

Reader::~Reader() {}

bool Reader::readMeshFile(Mesh &mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    int num_vertices;
    input >> num_vertices;

    if (num_vertices == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserveVectorSpace_Vertices(num_vertices);

    //legge i vertici aggiustando il dominio..
    for (int i = 0; i < num_vertices; i++) {
        double x, y, z;

        input >> x;
        input >> y;
        input >> z;
        if (input.eof())
            break;

        Vertex3D v = Vertex3D(x, y, z);
        v.setF(z);
        mesh.addVertex(v);
    }

    int num_topSimplexes;
    input >> num_topSimplexes;

    if(num_topSimplexes == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserveVectorSpace_TopSimplexes(num_topSimplexes);

    //legge i top simplessi
    for (int i = 0; i < num_topSimplexes; i++) {
        int v[3];
        for (int j = 0; j < 3; j++)
            input >> v[j];
        Triangle t = Triangle(v[0], v[1], v[2]);
        mesh.addTopSimplex(t);
    }

    return true;
}

bool Reader::readOFFMesh(Mesh &mesh, string path){

    FILE* file = fopen(path.c_str(), "r");

    int num_vertices, num_topsimplexes, num_edges; //edges is useless

    size_t len=0;
    char* offname=NULL;
    getline(&offname,&len,file);

    fscanf(file, "%d %d %d",&num_vertices,&num_topsimplexes,&num_edges);

    if(num_vertices == 0 || num_topsimplexes == 0){
        cerr<< "Number of simplexes is 0 Not a valid .off file "<<path<<endl;
        return false;
    }

    mesh.reserveVectorSpace_Vertices(num_vertices);

    //insert vertices
    for(unsigned long int i=0;i<num_vertices;i++){
        float x,y,z;
        fscanf(file,"%f %f %f", &x, &y, &z);


        Vertex3D v = Vertex3D(x,y,z);
        v.setF(z);
        mesh.addVertex(v);
    }

    mesh.reserveVectorSpace_TopSimplexes(num_topsimplexes);

    //insert top simplexes
    for(unsigned long int i=0;i<num_topsimplexes;i++){
        int v[4];
        fscanf(file,"%d %d %d %d", &v[0], &v[1], &v[2], &v[3]);

        Triangle t = Triangle(v[1],v[2],v[3]);
        mesh.addTopSimplex(t);
    }
    return true;
}

bool Reader::readScalarField(Mesh &mesh, string path){

    FILE* input = fopen(path.c_str(), "r");

    int vTot=0;
    fscanf(input, "%d", &vTot);

    if(vTot != mesh.getNumVertex()){
        cout << "[ERROR] - Vertices number in the input file mismatch with the number of vertices in the mesh ";
        cout << vTot << " vs " << mesh.getNumVertex() << endl;
        exit(0);
    }

    float val;
    for(int i=0; i<vTot; i++){
        fscanf(input, "%f ", &val);
        mesh.getVertex(i).setF(val);
    }

    fclose(input);

    return true;
}



bool Reader::readVTKFile(Mesh &mesh, string path){

    ifstream input(path);

    string a;
    for(int i=0; i<8; i++)
        input >> a;

    string name;
    int vertices;

    input >> name;
    input >> vertices;
    input >> name;

    mesh.reserveVectorSpace_Vertices(vertices);

    for(int i=0; i<vertices; i++){
        float x,y,z;
        input >> x;
        input >> y;
        input >> z;

        Vertex3D v = Vertex3D(x,y,z);
        v.setF(z);
        mesh.addVertex(v);
    }

    int triangles;

    input >> name;
    input >> triangles;
    input >> name;

    mesh.reserveVectorSpace_TopSimplexes(triangles);

    //insert top simplexes
    for(unsigned long int i=0;i<triangles;i++){
        int v[4];
        input >> v[0];
        input >> v[1];
        input >> v[2];
        input >> v[3];

        Triangle t = Triangle(v[1],v[2],v[3]);
        mesh.addTopSimplex(t);
    }


    input >> name;
    input >> triangles; //CELL_TYPES n_triangles
    int type=0;
    for(unsigned long int i=0;i<triangles;i++){
        input >> type;
    }

    for(int i=0; i<2; i++) //POINT_DATA vertices
        input >> a;

    int nFields;
    input >> name; //FIELD
    input >> name; //FieldData
    input >> nFields;

    input >> name; //field
    input >> name; //fieldValues

    int check_v;
    input >> check_v;
    input >> name;
    assert(vertices == check_v);

    float fValue;
    for(unsigned long int i=0;i<vertices;i++){
        input >> fValue;
        mesh.getVertex(i).setF(fValue);
    }

    input.close();

    return true;
}
