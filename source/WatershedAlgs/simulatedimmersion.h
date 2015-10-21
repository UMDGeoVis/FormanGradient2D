#ifndef SIMULATEDIMMERSION
#define SIMULATEDIMMERSION

#include "../LibMesh/Mesh.h"


///Watershed by simulated immersion algorithm originally described in
///Vincent and Soille, Watersheds in Digital Spaces: An Efficient Algorithm Based on Immersion Simulations

/// Main function: simulatedImmersionSegmentation
/// INPUT
///     -triangle mesh
///     -boolean value:
///         -true for computing the ascending Morse cells (regions associated with minima)
///         -false for computing the descending Morse cells (regions associated with maxima)
///     -vector of integers where the segmentation is stored
/// RETURN:
///     -number of regions found




struct Comparer{
    inline bool operator()(pair<int,float> v1, pair<int, float> v2){
        return v1.second > v2.second;
    }
};

int simulatedImmersionSegmentation(Mesh& mesh, bool ascending, vector<int>& v_labels){

    int nRegions=0;

    v_labels = vector<int>(mesh.getNumVertex(), -1);
    vector<bool> visited = vector<bool>(mesh.getNumVertex(), false);

    priority_queue<pair<int,float>, vector<pair<int,float> >, Comparer> sorted_vertices;

    if(ascending)
        for(int i=0; i<mesh.getNumVertex(); i++)
            sorted_vertices.push(pair<int,double>(i,mesh.getVertex(i).getF()));
    else
        for(int i=0; i<mesh.getNumVertex(); i++)
            sorted_vertices.push(pair<int,double>(i,-mesh.getVertex(i).getF()));

    vector<int> vv;
    set<int> local_labels;
    queue<int> flat_region;
    list<int> touched_vertexes;
    bool watershed_in_area;

    while(!sorted_vertices.empty()){
        watershed_in_area = false;

        vv.clear();
        local_labels.clear();
        touched_vertexes.clear();

        int vertice = sorted_vertices.top().first;
        sorted_vertices.pop();


        if(visited[vertice]) continue;

        flat_region.push(vertice);
        visited[vertice] = true;
        touched_vertexes.push_back(vertice);

        while(!flat_region.empty()){
            int v = flat_region.front();
            flat_region.pop();
            vv = mesh.VV(v);

            for(int j=0; j<vv.size(); j++){
                if(mesh.getVertex(v).getF() == mesh.getVertex(vv[j]).getF() ){
                    if(!visited[vv[j]]){
                        visited[vv[j]] = true;
                        flat_region.push(vv[j]);
                        touched_vertexes.push_back(vv[j]);
                    }
                }
                else{
                    if(v_labels[vv[j]] > -1) local_labels.insert(v_labels[vv[j]]);
                    if(v_labels[vv[j]] == -2) watershed_in_area = true;
                }
            }

        }

        int label=0;
        if(local_labels.size() == 0 && !watershed_in_area){
            label = nRegions++;
        }
        else if(local_labels.size() == 1) label = *(local_labels.begin());
        else label = -2;

        for(list<int>::iterator it = touched_vertexes.begin(); it != touched_vertexes.end(); it++)
            v_labels[*it] = label;
    }

    return nRegions;
}



#endif // SIMULATEDIMMERSION

