#include "Mesh.h"


void Mesh::laplacianSmoothing(int iter){

    for(int k=0; k<iter; k++){
        vector<float> nV(getNumVertex());
        for(int i=0;i<getNumVertex(); i++){
            vector<int> vv = VV(i);
            float value=0;
            for(auto v : vv){
                value+=getVertex(v).getF();
            }
            nV[i]=value/(float)vv.size();
        }

        for(int i=0; i<getNumVertex(); i++)
            getVertex(i).setF(nV[i]);
    }

}


