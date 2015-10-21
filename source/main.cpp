#include "./WatershedAlgs/simulatedimmersion.h"

#include "./LibForman/formangradientvector.h"

using namespace std;

void readInput(int argc, char** argv);
void changeMeshBorder(Mesh& mesh);

string meshFile("none");
string fieldFile("none");
string segmFile("none");


int main(int argc, char* argv[])
{

    MemoryUsage usage;
    Timer time;

    Mesh mesh = Mesh();

    readInput(argc, argv);


    if(meshFile.find(".tri") != string::npos)
        Reader::readMeshFile(mesh,meshFile);
    else if(meshFile.find(".off") != string::npos)
        Reader::readOFFMesh(mesh,meshFile);
    else if(meshFile.find(".vtk") != string::npos){
        Reader::readVTKFile(mesh, meshFile);
    }
    else{
        cout << "[Invalid input file] - supported formats: .off .tri .vtk" << endl;
        return 0;
    }
    mesh.build();
    cout << "Mesh with " << mesh.getNumVertex() << " vertices and " << mesh.getTopSimplexesNum() << "triangles [LOADED]" << endl;

    if(fieldFile != "none"){
        Reader::readScalarField(mesh,fieldFile);
    }


    FormanGradientVector gradient = FormanGradientVector(&mesh);

    vector<int> segm;
    if(segmFile == "immersion"){
        time.start();
        int nRegions = simulatedImmersionSegmentation(mesh,true,segm);
        cout << "Segmentation by Simulated immersion computed. " << endl;
        gradient.watershedToForman(segm);
    }
    else if(segmFile == "exp"){
        //Greedy discrimination of flat regions
        map<float, vector<int> > vert;
        for(int i=0; i<mesh.getNumVertex(); i++){
            vert[mesh.getVertex(i).getF()].push_back(i);
        }

        int count=0;
        for(auto m : vert){
            for(auto v : m.second){
                mesh.getVertex(v).setF(count++);
            }
        }
        time.start();
        gradient.homotopy_expansion();

    }
    else{

        assert(fieldFile == "none");
        cout << "Select a method for computing the Forman gradient. " << endl;
        cout << "   [immersion] Forman gradient combined with the watershed segmentaiton " << endl;
        cout << "   [exp] Forman gradient computed by homotopy expansion " << endl;
        return 1;
    }

    time.stop();
    cout << endl;
    cout << "Gradient computed in " << time.getElapsedTime() << " using " << usage.getValue_in_MB(false) << " MB" << endl;

    vector<int> cp = gradient.count_critical_simplexes();
    int numC = cp[0] + cp[1] + cp[2];
    cout << "Critical Points found (min sad max)" << cp[0] << " " << cp[1] << " " << cp[2] << " SUM: "<< numC << endl << endl;

///Uncomment for printing vtk files for the Forman gradient and the critical points
//    cout << "Writing Forman gradient on file: " << endl;
//    gradient.writeVTK_gradient("gradient.vtk");
//    cout << "- gradient.vtk" << endl;
//    gradient.writeVTK_criticalPoints("criticalPoints.vtk");
//    cout << "- criticalPoints.vtk" << endl << endl;

///Uncomment for printing vtk files with the Descending Morse complex
//    cout << "Writing Morse features on file: " << endl;
//    gradient.descending_2cells_extraction(true);
//    cout << "- descending2cells.vtk" << endl;
//    gradient.descending_1cells_extraction(true);
//    cout << "- descending1cells.vtk" << endl;

///Uncomment for printing vtk files with the Ascending Morse complex
//    gradient.ascending_2cells_extraction(true);
//    cout << "- ascending2cells.vtk" << endl;
//    gradient.ascending_1cells_extraction(true);
//    cout << "- ascending1cells.vtk" << endl;

    return 0;
}

void readInput(int argc, char** argv){

    //assert(strcmp("-i", argv[1])==0);
    int par=3;

    if(strcmp("-f", argv[1])==0){
        fieldFile = string(argv[3]);
        par=4;
    }

    meshFile = string(argv[2]);

    cout << endl;
    cout << "Mesh file        : " << meshFile << endl;
    cout << "Scalar field file: " << fieldFile << endl;

    if(argc > 3){


        segmFile = string(argv[par]);
        cout << "Segmentation method: " << segmFile << endl;
        assert(segmFile == "rain" || segmFile == "immersion" || segmFile == "exp" || segmFile == "filter");
    }
    cout << endl;
}

void changeMeshBorder(Mesh& mesh){

    vector<float> nValues(mesh.getNumVertex(),0);
    for(int i=0; i<mesh.getNumVertex(); i++){
        nValues[i]=mesh.getVertex(i).getF();

        if(mesh.isBoundary(i)){
            vector<int> vv = mesh.VV(i);
            int highest = -1;
            for(auto v : vv){
                if((highest == -1 || mesh.getVertex(highest).getF() <= mesh.getVertex(v).getF()) && !mesh.isBoundary(v)){
                    highest=v;
                }
            }

            if(highest != -1 && mesh.getVertex(highest).getF() <= mesh.getVertex(i).getF()){
                nValues[i]=(2*mesh.getVertex(highest).getF())-mesh.getVertex(i).getF();
            }

        }
    }

    for(int i=0; i<mesh.getNumVertex(); i++){
        mesh.getVertex(i).setF(nValues[i]);
    }

    for(int i=0; i<mesh.getNumVertex(); i++){
        nValues[i]=mesh.getVertex(i).getF();

        if(mesh.isBoundary(i)){
            vector<int> vv = mesh.VV(i);
            int highest = -1;
            for(auto v : vv){
                if((highest == -1 || mesh.getVertex(highest).getF() <= mesh.getVertex(v).getF()) && !mesh.isBoundary(v)){
                    highest=v;
                }
            }

            if(highest == -1){
                for(auto v : vv){
                    if((highest == -1 || mesh.getVertex(highest).getF() <= mesh.getVertex(v).getF()) ){
                        highest=v;
                    }
                }
                nValues[i]=(2*mesh.getVertex(highest).getF())-mesh.getVertex(i).getF();
            }

        }
    }

    for(int i=0; i<mesh.getNumVertex(); i++){
        mesh.getVertex(i).setF(nValues[i]);
    }

}
