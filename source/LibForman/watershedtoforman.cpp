#include "formangradientvector.h"

void FormanGradientVector::watershedToForman(vector<int>& segmentation){

    int nRegions = readSegmentation(segmentation);
    assert(segmentation.size() == mesh->getNumVertex());

    Timer time;
    MemoryUsage usage;
    time.start();
    cout << "Starting the gradient construction " << endl;
    vertexEdgePairings(segmentation);
    vertexEdgePairingsPlateau(segmentation);
    time.stop();
    cout << "Vertex-Edge [PAIRED] in " << time.getElapsedTime() << " occupying " << usage.getValue_in_MB(false) << " MB" << endl;

    time.start();
    edgeTrianglePairings(segmentation, nRegions);
    time.stop();
    cout << "Edge-Triangle inside regions [PAIRED] in " << time.getElapsedTime() << " occupying " << usage.getValue_in_MB(false) << " MB" << endl;

    time.start();
    boundaryPairings(segmentation);
    time.stop();
    cout << "Edge-Triangle on the boundaries [PAIRED] in " << time.getElapsedTime() << " occupying " << usage.getValue_in_MB(false) << " MB" << endl;

    time.start();

    cout << "Boundary paths [PAIRED] in " << time.getElapsedTime() << " occupying " << usage.getValue_in_MB(false) << " MB" << endl;
}

int FormanGradientVector::readSegmentation(vector<int>& segmentation){

    map<int,int> nRegions;
    int r=0;

    set<int> watershed = set<int>();

    for(int i=0; i<segmentation.size(); i++){
        assert(segmentation[i] != -1);
        if(segmentation[i] == -2)
            watershed.insert(i);
        else if(nRegions.find(segmentation[i]) == nRegions.end())
            nRegions[segmentation[i]]=r++;
    }

    if(watershed.size() > 0){
        queue<int> expandLabeling;
        for(auto v : watershed){
            vector<int> adjs = mesh->VV(v);
            for(auto v1 : adjs){
                if(mesh->getVertex(v).getF() >= mesh->getVertex(v1).getF() && segmentation[v1] != -2){
                    expandLabeling.push(v);
                }
            }
        }

        while(!expandLabeling.empty()){

            int v = expandLabeling.front();
            expandLabeling.pop();

            if(segmentation[v] >= 0)
                continue;

            vector<int> adjs = mesh->VV(v);
            int lowestLabeled=-1;
            for(auto v1 : adjs){
                if(mesh->getVertex(v).getF() >= mesh->getVertex(v1).getF() && segmentation[v1] != -2){
                    if(lowestLabeled==-1 || mesh->getVertex(lowestLabeled).getF() > mesh->getVertex(v1).getF())
                        lowestLabeled=v1;
                }

            }

            if(lowestLabeled == -1)
                expandLabeling.push(v);
            else{
                segmentation[v]=segmentation[lowestLabeled];
                for(auto v1 : adjs){
                    if(segmentation[v1] == -2)
                        expandLabeling.push(v1);
                }
            }
        }
    }


    return nRegions.size();
}


void FormanGradientVector::vertexEdgePairings(vector<int> const& vertexLabels){

    for(int i=0; i<vertexLabels.size(); i++){
        assert(vertexLabels[i]>=0);

        int other=i;
        vector<int> lowerVert = mesh->VV(i);
        for(int j=0; j<lowerVert.size(); j++){
            if(mesh->getVertex(lowerVert[j]).getF() < mesh->getVertex(other).getF() && vertexLabels[other] == vertexLabels[lowerVert[j]]){
                other = lowerVert[j];
            }
        }

        if(other != i){
            setVE(i,other);

        }
        else{
            //add vertex as critical
        }
    }

}


void FormanGradientVector::vertexEdgePairingsPlateau(vector<int> const& vertexLabels){

    vector<set<int> > plateau(mesh->getNumVertex(),set<int>());
    vector<bool> isMinimum(mesh->getNumVertex(),true);

    vector<bool> visitati = vector<bool>(mesh->getNumVertex(),false);

    int label=0;
    for(int i=0; i<mesh->getNumVertex(); i++){

        if(visitati[i]) continue;

        queue<int> cluster; //coda
        cluster.push(i);

        plateau[i].insert(i);
        visitati[i]=true;

        while(!cluster.empty()){

            int vertex = cluster.front();
            cluster.pop();

            vector<int> link = mesh->VV(vertex); //estraggo i vertici nell'intorno
            for(int j=0; j<link.size(); j++){

                if(mesh->getVertex(link[j]).getF() == mesh->getVertex(i).getF() && !visitati[link[j]]){
                    plateau[i].insert(link[j]);
                    visitati[link[j]] = true;
                    cluster.push(link[j]);

                }

                if(mesh->getVertex(link[j]).getF() < mesh->getVertex(i).getF() )
                    isMinimum[i]=false;

            }
        }

        label++;
    }

    for(int i=0; i<plateau.size(); i++){

        queue<int> next;
        if(plateau[i].size() > 1 && isMinimum[i]){
            set<int> plat_size = plateau[i];

            int minimum = pick_a_minimum(plat_size);

            vector<int> link = mesh->VV(minimum); //estraggo i vertici nell'intorno
            for(auto v : link){
                if(plat_size.find(v) != plat_size.end() && is_edge_critical(minimum,v)){
                    assert(vertexLabels[minimum] == vertexLabels[v]);
                    assert(getVE(v) == NULL);
                    setVE(v,minimum);
                    next.push(v);
                }
            }
        }
        else if(plateau[i].size() > 1){
            set<int> plat_size = plateau[i];
            for(auto v : plat_size){
                if(getVE(v) != NULL){
                    next.push(v);
                }
            }
        }


        while(!next.empty()){

            int v = next.front();
            next.pop();

            vector<int> link = mesh->VV(v);



            for(auto v1 : link){
                if(is_edge_critical(v,v1) && getVE(v1) == NULL && mesh->getVertex(v1).getF() >= mesh->getVertex(v).getF()
                        && vertexLabels[v1] == vertexLabels[v]){
                    setVE(v1,v);
                    next.push(v1);
                }
            }
        }
    }
}



void FormanGradientVector::edgeTrianglePairings(vector<int> const& vertexLabels, int nRegions){


    int pairings=0;
    vector<list<int> > trianglesPerRegion = vector<list<int> >(nRegions, list<int>());

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        if(vertexLabels[mesh->getTopSimplex(i).TV(0)] == vertexLabels[mesh->getTopSimplex(i).TV(1)] &&
           vertexLabels[mesh->getTopSimplex(i).TV(1)] == vertexLabels[mesh->getTopSimplex(i).TV(2)])
            trianglesPerRegion[vertexLabels[mesh->getTopSimplex(i).TV(0)]].push_back(i);
    }

    //here we can do a parallel implementation
    int tot_triangles=0;
    for(int i=0; i<nRegions; i++){
        list<int> triangles = trianglesPerRegion[i];
        tot_triangles += triangles.size();
        pairings += computePairingPerRegion(triangles, vertexLabels);
    }

}


void FormanGradientVector::boundaryPairings(vector<int> const& vertexLabels){

    vector<set<int> > plateau(mesh->getNumVertex(),set<int>());
    vector<bool> isMaximum(mesh->getNumVertex(),true);

    vector<bool> visitati = vector<bool>(mesh->getNumVertex(),false);
    vector<bool> maxTriangles = vector<bool>(mesh->getTopSimplexesNum(),false);

    for(int i=0; i<mesh->getNumVertex(); i++){
        //searching for maxima

        if(visitati[i]) continue;

        queue<int> cluster; //coda
        cluster.push(i);

        plateau[i].insert(i);
        visitati[i]=true;

        while(!cluster.empty()){

            int vertex = cluster.front();
            cluster.pop();

            vector<int> link = mesh->VV(vertex); //estraggo i vertici nell'intorno
            for(int j=0; j<link.size(); j++){

                if(mesh->getVertex(link[j]).getF() == mesh->getVertex(i).getF() && !visitati[link[j]]){
                    plateau[i].insert(link[j]);
                    visitati[link[j]] = true;
                    isMaximum[link[j]] = false;
                    cluster.push(link[j]);

                }

                if(mesh->getVertex(link[j]).getF() > mesh->getVertex(i).getF() )
                    isMaximum[i]=false;

            }
        }
    }

    for(int i=0; i<mesh->getNumVertex(); i++){
        if(isMaximum[i]){
            int max;
            if(plateau[i].size() > 1){
                max = pick_best_maximum(plateau[i]);
            }
            else{
                max = pick_best_maximum(i);
            }

            maxTriangles[max]=true;
        }
    }

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){

        if(maxTriangles[i]){

            //if the chosen triangle belongs to a region (of minima) lets create a path to the boundary
            if(vertexLabels[mesh->getTopSimplex(i).TV(0)] == vertexLabels[mesh->getTopSimplex(i).TV(1)] &&
               vertexLabels[mesh->getTopSimplex(i).TV(0)] == vertexLabels[mesh->getTopSimplex(i).TV(2)]){

                reversePath(i, vertexLabels);
            }
        }
    }


    for(int i=0; i<mesh->getTopSimplexesNum(); i++){

        if(maxTriangles[i]){

            if(vertexLabels[mesh->getTopSimplex(i).TV(0)] != vertexLabels[mesh->getTopSimplex(i).TV(1)] ||
               vertexLabels[mesh->getTopSimplex(i).TV(0)] != vertexLabels[mesh->getTopSimplex(i).TV(2)]){
                follow_descendingPath(i);
            }
        }
    }

    removeTerrainBoundaryPathsFast(maxTriangles);

}


void FormanGradientVector::removeTerrainBoundaryPathsFast(vector<bool> const& isMaximum){


    list<int> triangles;
    for(int i=0; i<mesh->getTopSimplexesNum(); i++){

        if(is_face_critical(i) && !isMaximum[i]){

            int count=0;
            Edge edge;
            for(int j=0; j<3; j++){
                Edge* e = mesh->getTopSimplex(i).TE(j);
                if(is_edge_critical(e->EV(0), e->EV(1)))
                {
                    count++;
                    edge = *e;
                }
                delete e;
            }

            assert(count != 0);

            if(count==1){
                setEF(mesh->getTopSimplex(i).TV(mesh->getTopSimplex(i).face_index(&edge)),i);
                follow_descendingPath(i);
            }
            else if(count == 3){
                int j=0;
                for(; j<3; j++){
                    int adj = mesh->getTopSimplex(i).TT(j);
                    if(adj == -1 || !is_face_critical(adj) || isMaximum[adj]){
                        setEF(mesh->getTopSimplex(i).TV(j),i);
                        follow_descendingPath(i);
                        break;
                    }
                }

                if(j==3)
                    triangles.push_back(i);
            }
            else
                triangles.push_back(i);

        }
    }


    //list<int> after;
    while(!triangles.empty()){

        list<int> pairing = triangles;
        triangles.clear();

        int minCount=4;
        int maxCount=0;
        int tMax=0;

        for( auto t : pairing){

            if(!is_face_critical(t))
                continue;

            int count=0;
            Edge edge;
            for(int j=0; j<3; j++){
                Edge* e = mesh->getTopSimplex(t).TE(j);
                if( is_edge_critical(e->EV(0), e->EV(1)) ){
                    count++;
                    edge = *e;
                }
                delete e;
            }

            if(minCount > count && count > 0)
                minCount = count;

            if(maxCount < count){
                maxCount = count;
                tMax=t;
            }

            if(count == 1){
                setEF(mesh->getTopSimplex(t).TV(mesh->getTopSimplex(t).face_index(&edge)), t);
                follow_descendingPath(t);
            }
            else if(count > 1){
                triangles.push_back(t);
            }

        }

        if(triangles.size() == 0)
            break;

        if(minCount >= 2){

            if(maxCount == 3){
                int j=0;
                for(; j<3; j++){
                    int adj = mesh->getTopSimplex(tMax).TT(j);
                    if(adj == -1 || !is_face_critical(adj) || isMaximum[adj]){
                        setEF(mesh->getTopSimplex(tMax).TV(j), tMax);
                        follow_descendingPath(tMax);
                        triangles.remove(tMax);
                        break;
                    }
                }

                if(j==3){
                    int t = triangles.front();
                    triangles.pop_front();
                    Edge* edge=NULL;
                    for(int i=0; i<3; i++){
                        Edge* e= mesh->getTopSimplex(t).TE(i);
                        if(is_edge_critical(e->EV(0), e->EV(1))){
                            if(!(mesh->getVertex(e->EV(0)).getF() <
                                 mesh->getVertex(mesh->getTopSimplex(t).TV(i)).getF() &&
                                 mesh->getVertex(e->EV(1)).getF() <
                                  mesh->getVertex(mesh->getTopSimplex(t).TV(i)).getF())

                                   )
                                edge = e;
                            else if(edge != NULL)
                                delete e;

                        }
                    }

                    assert(edge != NULL);
                    setEF(mesh->getTopSimplex(t).TV(mesh->getTopSimplex(t).face_index(edge)),t);
                    delete edge;
                    follow_descendingPath(t);
                }

            }
            else{
                int t = triangles.front();
                triangles.pop_front();
                Edge* edge=NULL;
                for(int i=0; i<3; i++){
                    Edge* e= mesh->getTopSimplex(t).TE(i);
                    if(is_edge_critical(e->EV(0), e->EV(1))){
                        if(!(mesh->getVertex(e->EV(0)).getF() <
                             mesh->getVertex(mesh->getTopSimplex(t).TV(i)).getF() &&
                             mesh->getVertex(e->EV(1)).getF() <
                              mesh->getVertex(mesh->getTopSimplex(t).TV(i)).getF())

                               )
                            edge = e;
                        else if(edge != NULL)
                            delete e;
                    }
                }

                assert(edge != NULL);
                setEF(mesh->getTopSimplex(t).TV(mesh->getTopSimplex(t).face_index(edge)),t);
                delete edge;
                follow_descendingPath(t);
            }
        }
    }

}

int FormanGradientVector::computePairingPerRegion(list<int> trianglesList, vector<int> const& vertexLabels){

    priority_queue<QElem, vector<QElem>, CompareQElem> edgeQueue= priority_queue<QElem, vector<QElem>, CompareQElem>();
    int pairings=0;

    for(auto t : trianglesList){
        for(int i=0; i<3; i++){
            int adjT = mesh->getTopSimplex(t).TT(i);
            if( adjT == -1 ||
               (vertexLabels[mesh->getTopSimplex(adjT).TV(0)] != vertexLabels[mesh->getTopSimplex(adjT).TV(1)] ||
                vertexLabels[mesh->getTopSimplex(adjT).TV(0)] != vertexLabels[mesh->getTopSimplex(adjT).TV(2)])){

                Edge* e = mesh->getTopSimplex(t).TE(i);
                assert(is_face_critical(t));

                if(is_face_critical(t) && is_edge_critical(e->EV(0), e->EV(1)))
                {
                    float height = (mesh->getVertex(e->EV(0)).getF() + mesh->getVertex(e->EV(1)).getF())/2.0;
                    edgeQueue.push(QElem(height,t,i));

                }
                delete e;
            }
        }
    }

    while(!edgeQueue.empty()){

        QElem elem = edgeQueue.top();
        edgeQueue.pop();

        int tIndex = elem.getTriangle();

        Edge* e = mesh->getTopSimplex(tIndex).TE(elem.getEdgeOnT());

        if(is_face_critical(tIndex) && is_edge_critical(e->EV(0),e->EV(1))){

            setEF(mesh->getTopSimplex(tIndex).TV(elem.getEdgeOnT()),tIndex);
            pairings++;

            for(int i=0; i<3; i++){

                if(i == elem.getEdgeOnT()) continue;

                int adjT = mesh->getTopSimplex(elem.getTriangle()).TT(i);
                if( adjT != -1 &&
                    vertexLabels[mesh->getTopSimplex(adjT).TV(0)] == vertexLabels[mesh->getTopSimplex(adjT).TV(1)] &&
                    vertexLabels[mesh->getTopSimplex(adjT).TV(0)] == vertexLabels[mesh->getTopSimplex(adjT).TV(2)]){

                    Edge* einside = mesh->getTopSimplex(elem.getTriangle()).TE(i);
                    if(is_face_critical(adjT) && is_edge_critical(einside->EV(0), einside->EV(1))){

                        float height = (mesh->getVertex(einside->EV(0)).getF() + mesh->getVertex(einside->EV(1)).getF()) /2.0;

                        int indexE = mesh->getTopSimplex(adjT).face_index(einside);
                        edgeQueue.push(QElem(height,adjT,indexE));
                    }
                    delete einside;
                }
            }
        }
        else if(is_edge_critical(e->EV(0),e->EV(1))){
            int eInd = getFE(tIndex);
            int adjT = mesh->getTopSimplex(tIndex).TT(eInd);
            if( adjT != -1 &&
                (vertexLabels[mesh->getTopSimplex(adjT).TV(0)] != vertexLabels[mesh->getTopSimplex(adjT).TV(1)] ||
                vertexLabels[mesh->getTopSimplex(adjT).TV(0)] != vertexLabels[mesh->getTopSimplex(adjT).TV(2)]) ){

                if(mesh->getVertex(mesh->getTopSimplex(tIndex).TV(eInd)).getF() >=
                   mesh->getVertex(mesh->getTopSimplex(tIndex).TV(elem.getEdgeOnT())).getF())
                {

                    freeEF(mesh->getTopSimplex(tIndex).TV(eInd),tIndex);
                    setEF(mesh->getTopSimplex(tIndex).TV(elem.getEdgeOnT()),tIndex);
                }
            }
        }

        delete e;
    }

    return pairings;
}




void FormanGradientVector::reversePath(int criticalT, vector<int> const& vertexLabels){

    //assert(!is_face_critical(criticalT));

    if(is_face_critical(criticalT))
        return;

    int tri = criticalT;
    int faceInT = getFE(tri);
    assert(faceInT != -1);
    freeEF(mesh->getTopSimplex(tri).TV(faceInT),tri);
    int adj = mesh->getTopSimplex(tri).TT(faceInT);

    while(adj != -1 && vertexLabels[mesh->getTopSimplex(adj).TV(0)] == vertexLabels[mesh->getTopSimplex(adj).TV(1)] &&
                       vertexLabels[mesh->getTopSimplex(adj).TV(0)] == vertexLabels[mesh->getTopSimplex(adj).TV(2)] ){

        Edge* e = mesh->getTopSimplex(tri).TE(faceInT);
        int faceInAdj = mesh->getTopSimplex(adj).face_index(e);

        faceInT = getFE(adj);

        if(mesh->getVertex(e->EV(0)).getF() < mesh->getVertex(mesh->getTopSimplex(adj).TV(faceInAdj)).getF() &&
           mesh->getVertex(e->EV(1)).getF() < mesh->getVertex(mesh->getTopSimplex(adj).TV(faceInAdj)).getF()){
            delete e;
            break;
        }
        delete e;


        freeEF(mesh->getTopSimplex(adj).TV(faceInT),adj);
        setEF(mesh->getTopSimplex(adj).TV(faceInAdj),adj);

        tri = adj;
        adj = mesh->getTopSimplex(tri).TT(faceInT);

    }

    assert(is_face_critical(criticalT));
}


int FormanGradientVector::pick_a_minimum(set<int> const& plateau){

    vector<float> barycenter(3,0);
    for(auto v : plateau){
        barycenter[0]+=mesh->getVertex(v).getX();
        barycenter[1]+=mesh->getVertex(v).getY();
        barycenter[2]+=mesh->getVertex(v).getZ();
    }

    barycenter[0]/=(float)plateau.size();
    barycenter[1]/=(float)plateau.size();
    barycenter[2]/=(float)plateau.size();

    int nearest=-1;
    float avg=0;
    for(auto v : plateau){

        float sum=sqrt(pow(mesh->getVertex(v).getX()-barycenter[0],2) +
                      pow(mesh->getVertex(v).getY()-barycenter[1],2) +
                      pow(mesh->getVertex(v).getZ()-barycenter[2],2) );


        if(nearest==-1 || avg > sum){
            avg = sum;
            nearest=v;
        }
    }

    assert(nearest != -1);

    return nearest;
}



int FormanGradientVector::pick_best_maximum(set<int> const& plateau){


    int best = -1;
    int count_best=0;

    for(auto i : plateau){
        vector<int> vt = mesh->VT(i);

        for(auto t : vt){

            if(is_face_critical(t)){
                int count=0;
                for(int e=0; e<3; e++){
                    Edge* edg = mesh->getTopSimplex(t).TE(e);
                    if(is_edge_critical(edg->EV(0), edg->EV(1)))
                        count++;
                }

                if(count > count_best || best == -1){
                    best = t;
                    count_best=count;
                }
            }
        }
    }

    if(best == -1){
        float max=0;
        for(auto i : plateau){
            vector<int> vt = mesh->VT(i);

            for(auto t : vt){
                float barycenter=0;
                int v1 = mesh->getTopSimplex(t).TV(0);
                int v2 = mesh->getTopSimplex(t).TV(1);
                int v3 = mesh->getTopSimplex(t).TV(2);

                barycenter = mesh->getVertex(v1).getF() + mesh->getVertex(v2).getF() + mesh->getVertex(v3).getF();
                barycenter /= 3.0;
                if(best == -1 || barycenter > max){
                    max = barycenter;
                    best = t;
                }
            }
        }
    }

    return best;
}



int FormanGradientVector::pick_best_maximum(int i){

    vector<int> vt = mesh->VT(i);
    int best = -1;
    int count_best=0;
    float max_bar;
    for(auto t : vt){

        if(is_face_critical(t)){
            int count=0;
            for(int e=0; e<3; e++){
                Edge* edg = mesh->getTopSimplex(t).TE(e);
                if(is_edge_critical(edg->EV(0), edg->EV(1)))
                    count++;
            }

            float barycenter=0;
            int v1 = mesh->getTopSimplex(t).TV(0);
            int v2 = mesh->getTopSimplex(t).TV(1);
            int v3 = mesh->getTopSimplex(t).TV(2);

            barycenter = mesh->getVertex(v1).getF() + mesh->getVertex(v2).getF() + mesh->getVertex(v3).getF();
            barycenter /= 3.0;


            if(count > count_best || best == -1){
                best = t;
                count_best=count;
                max_bar = barycenter;
            }
            else if(count == count_best && max_bar < barycenter){
                max_bar = barycenter;
                best = t;
            }
        }
    }

    if(best == -1){
        float max=0;
        for(auto t : vt){
            float barycenter=0;
            int v1 = mesh->getTopSimplex(t).TV(0);
            int v2 = mesh->getTopSimplex(t).TV(1);
            int v3 = mesh->getTopSimplex(t).TV(2);

            barycenter = mesh->getVertex(v1).getF() + mesh->getVertex(v2).getF() + mesh->getVertex(v3).getF();
            barycenter /= 3.0;
            if(best == -1 || barycenter > max){
                max = barycenter;
                best = t;
            }
        }
    }

    return best;
}

void FormanGradientVector::follow_descendingPath(int tri){

    queue<int> triangles;
    triangles.push(tri);

    while(!triangles.empty()){
        int tri = triangles.front();
        triangles.pop();


        for(int i=0; i<3; i++){
            Edge* edg = mesh->getTopSimplex(tri).TE(i);
            int adj=-1;
            if(is_edge_critical(edg->EV(0), edg->EV(1)) ){
                adj = mesh->getTopSimplex(tri).TT(i);
                if(adj != -1 && getFE(adj) == -1){

                    int indexE = mesh->getTopSimplex(adj).face_index(edg);

                    ///FOLLOWING CONDITION
                    if(!(mesh->getVertex(edg->EV(0)).getF() <
                        mesh->getVertex(mesh->getTopSimplex(adj).TV(indexE)).getF() &&
                        mesh->getVertex(edg->EV(1)).getF() <
                        mesh->getVertex(mesh->getTopSimplex(adj).TV(indexE)).getF()))

                    {

                        setEF(mesh->getTopSimplex(adj).TV(indexE),adj);
                        triangles.push(adj);
                    }
                }
                else{
                    //set edge critical
                }
            }
            delete edg;
        }
    }

}

