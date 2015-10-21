#ifndef _READER_H
#define	_READER_H

#include <vector>
#include <string>

#include "Vertex3D.h"
#include "Mesh.h"
#include "Triangle.h"

#include "assert.h"

using namespace std;
///A class that read file and initializite some relevant library structures
class Reader {
public:
    ///reader for .tri .off and .vtk files
    static bool readMeshFile(Mesh& mesh, string path);
    static bool readOFFMesh(Mesh& mesh, string path);
    static bool readVTKFile(Mesh& mesh, string path);
    static bool readScalarField(Mesh& mesh, string path);
private:
    ///A constructor method
    Reader();
    ///A constructor method
    Reader(const Reader& orig);
    ///A destructor method
    virtual ~Reader();
};

#endif	/* _READER_H */

