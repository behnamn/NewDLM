//
//  Loop.cpp
//  DLM_XC
//
//  Created by Behnam Najafi on 03/05/2017.
//  Copyright © 2017 Behnam Najafi. All rights reserved.
//

#include "Loop.h"

FaceManager::FaceManager(){}
FaceManager::FaceManager (Graph* g_) : g(g_) {
}

void Loop::calculate_weight(const Graph& g){
    weight = 0.;
    for (auto e = edges.begin(); e!= edges.end(); ++e){
        weight += g[*e].weight;
    }
}

void Loop::print(){
    cout << "Loop: ";
    for (auto e = edges.begin(); e!= edges.end(); ++e){
        cout << "(" << e->m_source << "->" << e->m_target << ")" << "\t";
    }
    cout << weight << endl;
}

void FaceManager::print(){
    cout << "---- Printing all loops ---\n";
    int i = 1;
    for (auto loop = faces.begin(); loop!= faces.end(); ++loop){
        cout << "Loop" << i << ": ";
        for (auto e = loop->edges.begin(); e!= loop->edges.end(); ++e){
            cout << "(" << e->m_source << "->" << e->m_target << ")" << "\t";
        }
        loop->calculate_weight(*g);
        cout << loop->weight << endl;
        i++;
    }
}







void print_embedding_storage(Graph& g, planar_embedding_storage_t& emb){
    vertex_t v;
    vertex_iterator_t vi, vi_end;
    vertex_map_t v_index = get(&VertexProperty::id, g);
    edge_map_t e_index = get(&EdgeProperty::id, g);
    for (boost::tie(vi,vi_end) = vertices(g); vi!=vi_end; ++vi){
        v = *vi;
        std::cout << v_index(v) << ":";
        for (auto e : emb[v_index(v)]){
            std::cout << e_index(e) << ",";
        }
        std::cout << "\t";
    }
    std::cout << std::endl;
}



traversal_visitor::traversal_visitor(const Graph &graph) :
        graph(graph), numFaces(0), total_weight(0), logsum(0){}

void traversal_visitor::begin_traversal() {
    numFaces = 0;
    logsum = 0;
    total_weight = 0;
    //std::cout<<"Beginning Travelsal"<<std::endl;
}

void traversal_visitor::begin_face() {
    face_weight = 0;
    //std::cout << "{ " << numFaces << " [";
}

void traversal_visitor::next_vertex(vertex_t v) {
    //std::cout << v << ",";
}

void traversal_visitor::next_edge(edge_t e) {
    //std::cout << e << " ";
    face_weight += graph[e].weight;
};

void traversal_visitor::end_face() {
    //std::cout << "] " << face_weight << " }";
    numFaces++;
    total_weight += face_weight;
    logsum += log(C_parameter/face_weight);
    //std::cout << "\t";
}
void traversal_visitor::end_traversal() {
    //std::cout << "Total Weight: " << total_weight;
    //std::cout << std::endl;
}
/*
template < typename Vertex >
void traversal_visitor::next_vertex(Vertex v) {}


template < typename Edge >
void traversal_visitor::next_edge(Edge e) {
    std::cout << e << " ";
    std::cout << graph[e].etype << " ";
}
 */

void set_face_data(Graph& g, planar_embedding_storage_t& emb, traversal_visitor& outvis){
    boost::planar_face_traversal(g, &emb[0], outvis, get(&EdgeProperty::id, g));
}
