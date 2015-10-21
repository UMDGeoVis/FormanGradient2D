#include "formangradientvector.h"
#include <stdio.h>
#include "assert.h"

#include <climits>

void FormanGradientVector::homotopy_expansion(){

    for(int i=0; i< mesh->getNumVertex(); i++){

        vector<int> vt = mesh->VT(i);
        map<vector<int>, Simplex_Graph* > lower_star = compute_lower_star(i, &vt);

        vector<int> v;
        v.push_back(i);

        if(lower_star.size() == 0){
            //minimum found
        }
        else{

            map<vector<int>, int> critical_simplexes_lower = map<vector<int>, int>();
            ExplicitGradient* gradient_vector_in_lower = compute_gradient(v, &lower_star, &critical_simplexes_lower);

            for(simplices_map::iterator it = gradient_vector_in_lower->begin(vector<int>(1)); it != gradient_vector_in_lower->end(vector<int>(1)); it++){
                setVE(it->first[0], it->second);
            }

            for(simplices_map::iterator it = gradient_vector_in_lower->begin(vector<int>(2)); it != gradient_vector_in_lower->end(vector<int>(2)); it++){
                int v = 3;
                v -= mesh->getTopSimplex(it->second).vertex_index(it->first[0]);
                v -= mesh->getTopSimplex(it->second).vertex_index(it->first[1]);

                setEF(mesh->getTopSimplex(it->second).TV(v), it->second);
            }

            delete gradient_vector_in_lower;
        }


        for(map<vector<int>, Simplex_Graph* >::iterator it = lower_star.begin(); it != lower_star.end(); it++)
            delete it->second;
    }

}


map<vector<int>, Simplex_Graph*> FormanGradientVector::compute_lower_star(int v, vector<int>* vt){

    map<vector<int>, Simplex_Graph*> mappa_costruzione;

    for(int i=0; i<vt->size(); i++){
        int t = vt->at(i);
        vector<int> simplex;
        simplex.push_back(v);

        Triangle t1 = mesh->getTopSimplex(t);


        for(int j=0; j<3; j++){
            int v1 = t1.TV(j);
            if(mesh->getVertex(v).getF() > mesh->getVertex(v1).getF()){
                simplex.push_back(v1);
            }

            assert(mesh->getVertex(v).getF() != mesh->getVertex(v1).getF() || v == v1);
        }
        if(simplex.size() == 1) continue;

        simplex= sort_simplex(&simplex);

        Simplex_Graph* simplex_grafo = new Simplex_Graph();
        simplex_grafo->simplex = simplex;
        simplex_grafo->coboundary = vector<Simplex_Graph*>();
        if(mappa_costruzione.find(simplex) == mappa_costruzione.end())
            mappa_costruzione[simplex] = simplex_grafo;


        if(simplex.size()>2){
            for(int j=1; j<simplex.size(); j++){

                vector<int> sub_simplex = simplex;
                sub_simplex.erase(sub_simplex.begin()+j);


                    map<vector<int>, Simplex_Graph*>::iterator it = mappa_costruzione.find(sub_simplex);
                    Simplex_Graph* sopra;
                    if(it == mappa_costruzione.end()){
                        sopra = new Simplex_Graph();
                        sopra->simplex = sub_simplex;
                        sopra->coboundary = vector<Simplex_Graph*>();
                        mappa_costruzione[sub_simplex] = sopra;
                    }
                    else sopra = it->second;

                push_in_coboundary(sopra, mappa_costruzione[simplex]);
            }
        }
    }

    return mappa_costruzione;
}

ExplicitGradient* FormanGradientVector::compute_gradient(vector<int>& v, map<vector<int>, Simplex_Graph*>* lower_star_structure, map<vector<int>, int>* critical_points){

    ExplicitGradient* gradient_v = new ExplicitGradient();

    vector<int> minimal_edge = vector<int>(0);
    vector<vector<int> > other_edges;

    for(map< vector<int>, Simplex_Graph* >::iterator it = lower_star_structure->begin(); it != lower_star_structure->end(); it++){
        if(it->first.size() == 2){
            bool same=false;

            if(minimal_edge.size()==0){
                minimal_edge = it->first;
            }
            else if(is_lower(it->first, minimal_edge, &same)){
                    other_edges.push_back(minimal_edge);
                    minimal_edge = it->first;
            }
            else other_edges.push_back(it->first);
        }
    }

    list<vector<int> > pq_one = list<vector<int> >();
    list<vector<int> > pq_zero = list<vector<int> >();

    gradient_v->insert(v, minimal_edge[1]);
    for(int i=0; i<other_edges.size(); i++){
        push_ordered(&pq_zero, other_edges.at(i));
    }


    vector<Simplex_Graph*> coboundary_faces = (*lower_star_structure)[minimal_edge]->coboundary;
    for(int i=0; i<coboundary_faces.size(); i++){
        coboundary_faces.at(i)->simplex;
        if(num_unpared_faces(coboundary_faces.at(i)->simplex, gradient_v, critical_points) == 1){
            push_ordered(&pq_one, coboundary_faces.at(i)->simplex);

        }
    }

    while(pq_one.size() != 0 || pq_zero.size() != 0){
        while(pq_one.size() != 0){
            vector<int> popped = pq_one.front();
            pq_one.pop_front();
            if(num_unpared_faces(popped, gradient_v, critical_points) == 0)
                push_ordered(&pq_zero, popped);
            else{

                vector<int> pair = unique_pairable_face(popped, gradient_v, critical_points);

                for(int i=0; i<popped.size();i++){
                    int j=0;
                    for(j=0; j<pair.size();j++){
                        if(pair.at(j) == popped.at(i)) break;
                    }
                    if(j == pair.size()) break;
                }

                if(pair.size() == 2){
                    gradient_v->insert(pair, vector_to_index(popped));
                }
                else{
                    gradient_v->insert(pair, extract_missing_id(pair,popped));
                }

                pq_zero.remove(pair);

                vector<Simplex_Graph*> coboundary_faces_popped = (*lower_star_structure)[popped]->coboundary;
                for(int i=0; i<coboundary_faces_popped.size(); i++){

                    if(num_unpared_faces(coboundary_faces_popped.at(i)->simplex, gradient_v, critical_points) == 1){
                        push_ordered(&pq_one, coboundary_faces_popped.at(i)->simplex);
                     }
                }

                vector<Simplex_Graph*> coboundary_faces_pair= (*lower_star_structure)[pair]->coboundary;
                for(int i=0; i<coboundary_faces_pair.size(); i++){

                    if(num_unpared_faces(coboundary_faces_pair.at(i)->simplex, gradient_v, critical_points) == 1)
                        push_ordered(&pq_one, coboundary_faces_pair.at(i)->simplex);
                }
            }

        }

        if(pq_zero.size() != 0){
            vector<int> critico = pq_zero.front();
            pq_zero.pop_front();

            if(critico.size() == 2){
                (*critical_points)[critico] = critico.front();

            }
            if(critico.size() == 3){
                (*critical_points)[critico] = critico.front();
            }

            vector<Simplex_Graph*> coboundary_faces_critico= (*lower_star_structure)[critico]->coboundary;
            for(int i=0; i<coboundary_faces_critico.size(); i++){
                if(num_unpared_faces(coboundary_faces_critico.at(i)->simplex, gradient_v, critical_points) == 1)
                    push_ordered(&pq_one, coboundary_faces_critico.at(i)->simplex);
            }

        }
    }


    return gradient_v;
}

int FormanGradientVector::num_unpared_faces(vector<int>& co_face, ExplicitGradient* gradient_v, map<vector<int>, int>* critical_points){


    int num_unpaired=0;
    vector<int> simplex = co_face;

    for(int i=1; i<simplex.size(); i++){
        vector<int> simpl_face = simplex;
        simpl_face.erase(simpl_face.begin()+i);


        if(gradient_v->find(simpl_face) != gradient_v->end(simpl_face)){
            continue;
        }

        if(critical_points->find(simpl_face) != critical_points->end()){
            continue;
        }


        int j=-1;
        if(simpl_face.size() > 1){
            for(j=1; j<simpl_face.size();j++){
                vector<int> simpl_face_second = simpl_face;
                simpl_face_second.erase(simpl_face_second.begin()+j);

                simplices_map::iterator it = gradient_v->find(simpl_face_second);
                if(it != gradient_v->end(simpl_face_second)){
                    vector<int> simpl;

                    if(simpl_face_second.size() == 1){
                        simpl = simpl_face_second;
                        simpl.push_back(it->second);
                    }
                    else{
                        for(int pos=0; pos<3; pos++){
                            simpl.push_back(mesh->getTopSimplex(it->second).TV(pos));
                        }
                    }
                    simpl = sort_simplex(&simpl);

                    if(simpl == simpl_face)
                        break;
                }

            }
        }
        if(j==simpl_face.size()){num_unpaired++;}
    }

    return num_unpaired;
}

vector<int> FormanGradientVector::unique_pairable_face(vector<int>& s_face, ExplicitGradient* gradient_v, map<vector<int>, int>* critical_points){

    vector<int> simplex = s_face;
    for(int i=1; i<simplex.size(); i++){
        vector<int> simpl_face = simplex;
        simpl_face.erase(simpl_face.begin()+i);

        if(gradient_v->find(simpl_face) != gradient_v->end(simpl_face)) continue;
        if(critical_points->find(simpl_face) != critical_points->end()) continue;

        int j=-1;
        if(simpl_face.size() > 1){
            for(j=1; j<simpl_face.size();j++){
                vector<int> simpl_face_second = simpl_face;
                simpl_face_second.erase(simpl_face_second.begin()+j);

                simplices_map::iterator it = gradient_v->find(simpl_face_second);
                if(it != gradient_v->end(simpl_face_second) ){

                    vector<int> simpl;

                    if(simpl_face_second.size() == 1){
                        simpl = simpl_face_second;
                        simpl.push_back(it->second);
                    }
                    else{
                        for(int pos=0; pos<3; pos++){
                            simpl.push_back(mesh->getTopSimplex(it->second).TV(pos));
                        }
                    }
                    simpl = sort_simplex(&simpl);

                    if(simpl == simpl_face)
                        break;
                }
            }
        }
        if(j==simpl_face.size()) return simpl_face;
    }


}


void FormanGradientVector::push_in_coboundary(Simplex_Graph* sface, Simplex_Graph* co_sface){

    int i=0;
    for(; i<sface->coboundary.size(); i++){
        if(sface->coboundary.at(i)->simplex == co_sface->simplex){
            break;
        }
    }

    if(i==sface->coboundary.size()) sface->coboundary.push_back(co_sface);
}
