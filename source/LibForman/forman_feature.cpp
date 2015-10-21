#include "formangradientvector.h"

vector<int> FormanGradientVector::descending_2cells_extraction(bool with_geometry){

    vector<int> labeling = vector<int>(mesh->getTopSimplexesNum(),-1);
    for(int i=0; i<mesh->getTopSimplexesNum(); i++){

        if(is_face_critical(i)){

            queue<Edge* > coda;
            labeling[i]=i;

            vector<int> simplexi;
            int v=-1;
            for(int j=0; j<3; j++){
                    simplexi.push_back(mesh->getTopSimplex(i).TV(j));

            }

            for(int j=0; j<3; j++){
                vector<int> simplex = simplexi;
                simplex.erase(simplex.begin()+j);
                coda.push(new Edge(simplex[0], simplex[1]));
            }
            simplexi.clear();

            while(!coda.empty()){
                Edge* edge = coda.front();
                coda.pop();

                int tri = getEF(edge);
                if(tri != -1){

                    labeling[tri]=i;

                    for(int j=0; j<3; j++){
                        Edge* new_e = mesh->getTopSimplex(tri).TE(j);
                        if(*new_e != *edge){
                            coda.push(new_e);
                        }
                        else delete new_e;
                    }

                    delete edge;
                }
            }
        }
    }


    if(with_geometry)
        writeVTK_2cells("descending2cells.vtk", labeling);

    return labeling;
}


void FormanGradientVector::descending_1cells_extraction(bool with_geometry){

    queue<int> coda;

    set<int> vertici;
    set<pair<int, bool> > critici;
    set<pair<int,int> > edges;

    vector<bool> visited = vector<bool>(mesh->getTopSimplexesNum(),false);
    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        visited[i]=true;
        for(int j=0; j<3; j++){
            if(mesh->getTopSimplex(i).TT(j) == -1 || !visited[mesh->getTopSimplex(i).TT(j)]){

                Edge* edge = mesh->getTopSimplex(i).TE(j);
                if(is_edge_critical(edge->EV(0), edge->EV(1))){


                    if(with_geometry){
                        mesh->getVertex(edge->EV(0)).getZ() > mesh->getVertex(edge->EV(1)).getZ() ? edges.insert(pair<int,int>(edge->EV(0),edge->EV(1))) : edges.insert(pair<int,int>(edge->EV(1),edge->EV(0)));
                        vertici.insert(edge->EV(0));
                        vertici.insert(edge->EV(1));
                    }

                    coda.push(edge->EV(0));
                    coda.push(edge->EV(1));
                    delete edge;

                    while(!coda.empty()){

                        int vert = coda.front();
                        coda.pop();

                        edge = getVE(vert);
                        if(edge != NULL){
                            if(with_geometry){
                                mesh->getVertex(edge->EV(0)).getZ() > mesh->getVertex(edge->EV(1)).getZ() ? edges.insert(pair<int,int>(edge->EV(0),edge->EV(1))) : edges.insert(pair<int,int>(edge->EV(1),edge->EV(0)));
                                vertici.insert(edge->EV(0));
                                vertici.insert(edge->EV(1));
                            }
                            int v2 = edge->EV(0) == vert ? edge->EV(1) : edge->EV(0);
                            delete edge;
                            coda.push(v2);
                        }
                    }
                }
            }
        }
    }


    if(with_geometry)
        writeVTK_1cells("descending1cells.vtk", vertici, edges);
}

vector<int> FormanGradientVector::ascending_2cells_extraction(bool with_geometry){

    queue<int> coda;
    vector<int> label = vector<int>(mesh->getNumVertex(), -1);

    int critical_edge =0;
    int minima=0;

    for(int i=0; i<mesh->getNumVertex(); i++){
        if(is_vertex_critical(i)){
            minima++;
            label[i]=i;

            vector<Edge*> ve = mesh->VE(i);
            for(unsigned int j=0; j<ve.size(); j++){
                int v2 = ve[j]->EV(0) == i ? ve[j]->EV(1) : ve[j]->EV(0);
                Edge* edge = getVE(v2);
                if(edge != NULL && *edge == *(ve[j])){
                    coda.push(v2);
                }

                if(edge != NULL) delete edge;
                delete ve[j];
            }


            while(!coda.empty()){

                int v = coda.front();
                coda.pop();
                label[v]=i;

                ve = mesh->VE(v);
                for(unsigned int j=0; j<ve.size(); j++){
                    int v2 = ve[j]->EV(0) == v ? ve[j]->EV(1) : ve[j]->EV(0);
                    Edge* edge = getVE(v2);
                    if(edge != NULL && *edge == *(ve[j])){
                        coda.push(v2);
                    }
                    if(edge != NULL) delete edge;
                    delete ve[j];
                }
            }

        }
    }


    if(with_geometry)
        writeVTK_2cells_on_vert("ascending2cells.vtk", label);

    return label;
}

void FormanGradientVector::ascending_1cells_extraction(bool with_geometry){


    vector<bool> visited = vector<bool>(mesh->getTopSimplexesNum(),false);
    map<int,int> visited_triangle;
    set<pair<Vertex3D,Vertex3D> > edges;

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        visited[i]=true;
        for(int j=0; j<3; j++){
            if(mesh->getTopSimplex(i).TT(j) == -1 || !visited[mesh->getTopSimplex(i).TT(j)]){

                queue<int> coda;
                Edge* edge = mesh->getTopSimplex(i).TE(j);
                if(is_edge_critical(edge->EV(0), edge->EV(1))){

                    vector<int> et = mesh->ET(*edge);
                    for(int j=0; j<et.size(); j++){
                        coda.push(et[j]);
                    }

                    while(!coda.empty()){

                        int t = coda.front();
                        coda.pop();

                        visited_triangle[t]=i;

                        if(is_face_critical(t)){
                            //here critical face
                        }
                        else{
                            int e = getFE(t);

                            int t_adj = mesh->getTopSimplex(t).TT(e);
                            if(t_adj != -1){
                                coda.push(t_adj);



                                Vertex3D v1 = mesh->getVertex(mesh->getTopSimplex(t_adj).TV(0));
                                v1 += mesh->getVertex(mesh->getTopSimplex(t_adj).TV(1));
                                v1 += mesh->getVertex(mesh->getTopSimplex(t_adj).TV(2));
                                v1 /= 3.0;

                                Edge* edg = mesh->getTopSimplex(t).TE(e);
                                Vertex3D ve = mesh->getVertex(edg->EV(0));
                                ve += mesh->getVertex(edg->EV(1));
                                ve /= 2.0;


                                Vertex3D v2 = mesh->getVertex(mesh->getTopSimplex(t).TV(0));
                                v2 += mesh->getVertex(mesh->getTopSimplex(t).TV(1));
                                v2 += mesh->getVertex(mesh->getTopSimplex(t).TV(2));
                                v2 /= 3.0;

                                edges.insert(pair<Vertex3D,Vertex3D>(v1,ve));
                                edges.insert(pair<Vertex3D,Vertex3D>(v2,ve));
                            }
                        }
                    }

                }
            }
        }
    }

    if(with_geometry)
        writeVTK_1cells("ascending1cells.vtk",edges);
}



vector<int> FormanGradientVector::count_critical_simplexes(){
    vector<int> cp(3,0);

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        if(is_face_critical(i))
            cp[2]++;
    }

    for(int i=0; i<mesh->getNumVertex(); i++){
        if(is_vertex_critical(i)){
            cp[0]++;
        }

        vector<int> vv = mesh->VV(i);
        for(auto v : vv){
            if(v < i){
                if(is_edge_critical(v,i))
                    cp[1]++;
            }
        }
    }

    return cp;
}
