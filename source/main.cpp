#include "./WatershedAlgs/simulatedimmersion.h"
#include "./LibForman/formangradientvector.h"

using namespace std;

void readInput(int argc, char** argv);

string meshFile("none");
string fieldFile("none");
string segmFile("none");


int main(int argc, char* argv[])
{

    MemoryUsage usage;
    Timer time;

    Mesh mesh = Mesh();

    //How to read an input mesh file
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

    mesh.build(); //build the mesh structure
    cout << "Mesh with " << mesh.getNumVertex() << " vertices and " << mesh.getTopSimplexesNum() << "triangles [LOADED]" << endl;


    //how to read a new scalar field
    if(fieldFile != "none"){
        Reader::readScalarField(mesh,fieldFile);
    }

    FormanGradientVector gradient = FormanGradientVector(&mesh);


    vector<int> segm;
    if(segmFile == "immersion"){
        time.start();
        //compute the watershed segmentation with a simulated immersion approach
        int nRegions = simulatedImmersionSegmentation(mesh,true,segm);
        cout << "Segmentation by Simulated immersion computed. " << nRegions << endl;
        //compute the Forman gradient using the watershe segmentation from [Computing a Forman Gradient From a Watershed Decomposition, Comic et al, 2016]
        gradient.watershedToForman(segm);
    }
    else if(segmFile == "exp"){
        //Greedy discrimination of flat regions with Simulation of Simplicity from [Simulation of Simplicity, Edelsbrunner et al., 1990]
        map<double, vector<int> > vert;
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
        //compute the Forman gradient via homotopty expansion
        gradient.homotopy_expansion();

    }
    else if(segmFile == "exp-improved"){
        //Improved simulation of simplicity from [Practical Considerations in Morse-Smale Complex Computation, Gyulassy et al. 2010]

        priority_queue<pair<int,float>, vector<pair<int,float> >, Comparer> simpleOrder;
        priority_queue<pair<int,float>, vector<pair<int,float> >, Comparer> bfsOrder;

        int count=0;
        vector<int> newF(mesh.getNumVertex(),-1);

        for(int i=0; i<mesh.getNumVertex(); i++)
            simpleOrder.push(pair<int,double>(i,mesh.getVertex(i).getF()));

        while(!simpleOrder.empty()){

            if(bfsOrder.empty() || bfsOrder.top().second > simpleOrder.top().second){

                int vertex = simpleOrder.top().first;
                simpleOrder.pop();

                if(newF[vertex]!=-1)
                    continue;

                newF[vertex]=count++;

                vector<int> vv = mesh.VV(vertex);
                for(auto adj : vv){
                    if(newF[adj]==-1){
                        bfsOrder.push(pair<int,double>(adj,mesh.getVertex(adj).getF()));
                    }
                }
            }

            else {

                int vertex = bfsOrder.top().first;
                bfsOrder.pop();

                if(newF[vertex]!=-1)
                    continue;

                newF[vertex]=count++;

                vector<int> vv = mesh.VV(vertex);
                for(auto adj : vv){
                    if(newF[adj]==-1){
                        bfsOrder.push(pair<int,double>(adj,mesh.getVertex(adj).getF()));
                    }
                }
            }
        }


        for(int i=0; i<newF.size(); i++){
            assert(newF[i] != -1);
            mesh.getVertex(i).setF(newF[i]);
        }

        time.start();
        //compute the Forman gradient via homotopty expansion
        gradient.homotopy_expansion();

    }
    else{
        assert(fieldFile == "none");
        cout << "Select a method for computing the Forman gradient. " << endl;
        cout << "   [immersion] Forman gradient combined with the watershed segmentaiton " << endl;
        cout << "   [exp] simulation of simplicity + homotopy expansion " << endl;
        cout << "   [exp-improved] improved simulation of simplicity + homotopy expansion " << endl;
        return 1;
    }


    time.stop();
    cout << endl;
    cout << "Gradient computed in " << time.getElapsedTime() << " using " << usage.getValue_in_MB(false) << " MB" << endl;

    vector<int> cp = gradient.count_critical_simplexes();
    int numC = cp[0] + cp[1] + cp[2];
    cout << "Critical Points found (min sad max)" << cp[0] << " " << cp[1] << " " << cp[2] << " SUM: "<< numC << endl << endl;

/////Uncomment for printing vtk files for the Forman gradient and the critical points
//    cout << "Writing Forman gradient on file: " << endl;
//    gradient.writeVTK_gradient("gradient.vtk");
//    cout << "- gradient.vtk" << endl;
//    gradient.writeVTK_criticalPoints("criticalPoints.vtk");
//    cout << "- criticalPoints.vtk" << endl << endl;

///////Uncomment for printing vtk files with the Descending Morse complex
//    cout << "Writing Morse features on file: " << endl;
//    gradient.descending_2cells_extraction(true);
//    cout << "- descending2cells.vtk" << endl;
//    gradient.descending_1cells_extraction(true);
//    cout << "- descending1cells.vtk" << endl;

///////Uncomment for printing vtk files with the Ascending Morse complex
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
        assert(segmFile == "immersion" || segmFile == "exp" || segmFile == "exp-improved" || segmFile == "filter");
    }
    cout << endl;
}


