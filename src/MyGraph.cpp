/*
 *
 *	MyGraph.cpp
 * 	Author: Behnam
 *
 */

#include "MyGraph.h"

void reset_edge_index_single(Graph& g) {
    edge_iterator_t ei,ef;
    tie(ei,ef) = edges(g);
    int i = 0;
    for ( ; ei != ef ; ++ei){
        g[*ei].id = i;
        i++;
    }
}


void remove_long_crossovers(Graph& graph, const Design* design){
    for (auto pool = design->staple_pools.begin(); pool!=design->staple_pools.end(); ++pool) {
        for (auto cross = pool->crossovers.begin(); cross!=pool->crossovers.end(); ++cross){
            if (cross->type == 'l') {
                if (cross->edge.second)
                    boost::remove_edge(cross->edge.first, graph);
            }
        }
    }
    reset_edge_index_single(graph);
}
void readd_long_crossovers(Graph& graph, Design* design){
    EdgeProperty EP;
    for (auto pool = design->staple_pools.begin(); pool!=design->staple_pools.end(); ++pool) {
        for (auto cross = pool->crossovers.begin(); cross!=pool->crossovers.end(); ++cross){
            if (cross->type == 'l') {
                if (cross->edge.second) {
                    EP.domain.second = false;
                    EP.crossover = make_pair(cross,true);
                    EP.weight = cs_hack;
                    EP.id = boost::num_edges(graph);
                    cross->edge.first = add_edge((cross->nodes.first)->id, (cross->nodes.second)->id, EP, graph).first;
                }
            }
        }
    }
}


bool sort_decending (int i,int j) { return (i>j); }


Comp_Graph::Comp_Graph(){}
Comp_Graph::Comp_Graph (int pool_id_, string pool_name_, int num_vertices_):
pool_id(pool_id_), pool_name(pool_name_), num_vertices(num_vertices_){
    this->component.resize(num_vertices);
    this->comp_sizes.resize(num_vertices);
    this->num_comps.resize(num_vertices);
    this->add_OPs();
}
void Comp_Graph::add_OPs(){
    OrderParameter OP;
    string name;
    name = "CLs"; if (pool_id > -1) {name+= "_"; name+= pool_name;}
    OP = OrderParameter(pool_id,NUM_CLUSTERS,name);
    this->OPs.push_back(OP);
    name = "LC1"; if (pool_id > -1) {name+= "_"; name+= pool_name;}
    OP = OrderParameter(pool_id,LARGE_CLUSTER,name);
    this->OPs.push_back(OP);
    name = "LC2"; if (pool_id > -1) {name+= "_"; name+= pool_name;}
    OP = OrderParameter(pool_id,LARGE_CLUSTER,name);
    this->OPs.push_back(OP);
    name = "LC3"; if (pool_id > -1) {name+= "_"; name+= pool_name;}
    OP = OrderParameter(pool_id,LARGE_CLUSTER,name);
    this->OPs.push_back(OP);
    name = "C0S"; if (pool_id > -1) {name+= "_"; name+= pool_name;}
    OP = OrderParameter(pool_id,LARGE_CLUSTER,name);
    this->OPs.push_back(OP);
}
void Comp_Graph::fill_components(){
    for (auto s=comp_sizes.begin(); s!=comp_sizes.end();++s){
        *s = 0;
    }
    num_components = connected_components(g, &component[0]);
    for (auto c = component.begin(); c!= component.end(); ++c){
        comp_sizes[*c]++;
    }
    this->OPs[4].state = this->comp_sizes[component[0]];
    std::sort (comp_sizes.begin(),comp_sizes.end()-1,sort_decending);
    for (auto n=num_comps.begin(); n!=num_comps.end();++n){
        *n = 0;
    }
    for (auto c = comp_sizes.begin(); c!=comp_sizes.end(); ++c){
        for (int i=0; i<num_comps.size(); i++){
            if (*c==i){
                num_comps[i]++;
            }
        }
    }
    this->OPs[0].state = this->num_components;
    this->OPs[1].state = this->comp_sizes[0];
    this->OPs[2].state = this->comp_sizes[1];
    this->OPs[3].state = this->comp_sizes[2];
}


//Constructor
MyGraph::MyGraph (Design* design_) : design(design_) {
    inputs = design->inputs;
    if (inputs->track_clusters){
        cg = Comp_Graph(-1,"",design->num_domains);
        if (design->num_staple_pools > 1){
            for (auto pool = design->staple_pools.begin(); pool!= design->staple_pools.end(); ++pool){
                p_cg.push_back(Comp_Graph(pool->id,pool->name,pool->num_domains));
            }
        }
    }

    add_vertices();
    if (inputs->make_movie) add_layouts();
    add_domains();
    initialise_state();
    
    distances.resize(design->num_domains);

    fill_components();
    //update_faces();
}
//Common Methods
void MyGraph::add_vertices(){
    // global 2
    embedding_storage.resize(design->num_domains);

    for(int i=0; i<design->num_domains; i++){
        add_vertex(g);
        if (inputs->track_clusters) add_vertex(cg.g);
    }
    if (inputs->track_clusters){
        for (auto cg = p_cg.begin(); cg!= p_cg.end(); ++cg){
            auto pool = design->staple_pools.begin()+cg->pool_id;
            for (auto dom = pool->domains.begin(); dom!= pool->domains.end(); ++dom){
                add_vertex(cg->g);
            }
        }
    }
    reset_vertex_index();
}
void MyGraph::add_layouts(){
    vertex_iterator_t vi,vf;
    for(tie(vi,vf) = vertices(g); vi != vf; ++vi) {
        g[*vi].pos = ",";
        g[*vi].pos += "!";
        g[*vi].height = 0.3;
        g[*vi].width = 0.3;
        g[*vi].other = "true";
    }
    for (int p = 0; p< p_cg.size(); p++) {
        for(tie(vi,vf) = vertices(p_cg[p].g); vi != vf; ++vi) {
            p_cg[p].g[*vi].pos = ",";
            p_cg[p].g[*vi].pos += "!";
            p_cg[p].g[*vi].height = 0.3;
            p_cg[p].g[*vi].width = 0.3;
            p_cg[p].g[*vi].other = "true";
        }
    }
}
void MyGraph::add_domains(){
    EdgeProperty EP;
    //EP.type = 's';
    EP.crossover.second = false;
    for (auto d = design->domains.begin(); d != design->domains.end(); ++d) {
        EP.id = std::distance(design->domains.begin(), d);
        //EP.domain = make_pair(d,true);
        EP.domain_s = make_pair(d,true);
        EP.weight = d->length * l_ss * lambda_ss;
        d->edge = add_edge(d->node_ids.first, d->node_ids.second, EP, g).first;

        embedding_storage[d->node_ids.first].push_back(d->edge);
    }

    vertex_t v1, v2;
    auto vs =  boost::make_iterator_range(vertices(g));
    for (auto vi = vs.begin(); vi!= vs.end(); ++vi){
        v1 = *vi;
        if (vi == std::prev(vs.end())) {v2 = *vs.begin();}
        else{v2 = *std::next(vi);}
        embedding_storage[g[v2].id].push_back(embedding_storage[g[v1].id][0]);
    }
}


//State update Methods
void MyGraph::initialise_state(){
    for (auto pool = design->staple_pools.begin(); pool!= design->staple_pools.end(); ++pool){
        for (auto cr = pool->crossovers.begin(); cr!= pool->crossovers.end(); ++cr){
            if (cr->state){ add_crossover(cr); }
        }
        for (auto dom = pool->domains.begin(); dom!= pool->domains.end(); ++dom){
            if (dom->state) { bind_domain(dom); }
        }
    }
}
void MyGraph::bind_domain(const PDOM domain){
    if (inputs->track_clusters && design->num_staple_pools > 1){
        EdgeProperty EP;
        EP.crossover.second = false;
        EP.id = domain->id;
        EP.domain_p = make_pair(domain,true);
        domain->edge_p = add_edge(domain->node_idps.first, domain->node_idps.second, EP, p_cg[domain->pool_id].g);
        reset_edge_index();
    }
    bind_domains(domain->children);
}
void MyGraph::unbind_domain(const PDOM domain){
    if (inputs->track_clusters && design->num_staple_pools > 1){
        remove_edge(domain->edge_p.first, p_cg[domain->pool_id].g);
        domain->edge_p.second = false;
        reset_edge_index();
    }
    unbind_domains(domain->children);
}
void MyGraph::bind_domains(vector<SDOM>& domains){
    for(auto domain = domains.begin(); domain!= domains.end(); ++domain){
        g[(*domain)->edge].weight = (*domain)->length * (*domain)->length * ds_hack;
        if (inputs->track_clusters){
            EdgeProperty EP;
            EP.crossover.second = false;
            EP.id = (*domain)->id;
            EP.domain_s = make_pair((*domain),true);
            (*domain)->edge2 = add_edge((*domain)->node_ids.first, (*domain)->node_ids.second, EP, cg.g);
            reset_edge_index();
        }
    }
}
void MyGraph::unbind_domains(vector<SDOM>& domains){
    for(auto domain = domains.begin(); domain!= domains.end(); ++domain){
        g[(*domain)->edge].weight = (*domain)->length * ss_hack;
        if (inputs->track_clusters){
            remove_edge((*domain)->edge2.first, cg.g);
            (*domain)->edge2.second = false;
            reset_edge_index();
        }
    }
}
void MyGraph::add_crossover(CR cross){
    EdgeProperty EP;
    EP.domain.second = false;
    EP.crossover = make_pair(cross,true);
    EP.weight = cs_hack;
    EP.id = boost::num_edges(g);
    cross->edge = add_edge((cross->nodes.first)->id, (cross->nodes.second)->id, EP, g);
    this->add2embedding(cross);

    if (inputs->track_clusters){
        EP.id = boost::num_edges(cg.g);
        cross->edge2 = add_edge((cross->nodes.first)->id, (cross->nodes.second)->id, EP, cg.g);
        if (design->num_staple_pools > 1)
            EP.id = boost::num_edges(p_cg[cross->pool_id].g);
        cross->edge_p = add_edge((cross->nodes.first)->idp, (cross->nodes.second)->idp, EP, p_cg[cross->pool_id].g);
    }
}
void MyGraph::remove_crossover(CR cross){
    this->removeFromEmbedding(cross);
    remove_edge(cross->edge.first, g);
    cross->edge.second = false;
    reset_edge_index_single(g);
    if (inputs->track_clusters){
        remove_edge(cross->edge2.first, cg.g);
        cross->edge2.second = false;
        if (design->num_staple_pools > 1){
            remove_edge(cross->edge_p.first, p_cg[cross->pool_id].g);
            cross->edge_p.second = false;
        }
        reset_edge_index();
    }

}

void MyGraph::reset_vertex_index() {
    vertex_iterator_t vi,vf;
    int i = 0;
    for(tie(vi,vf) = vertices(g); vi != vf; ++vi) {
        g[*vi].id = i;
        i++;
    }
    i = 0;
    for(tie(vi,vf) = vertices(cg.g); vi != vf; ++vi) {
        cg.g[*vi].id = i;
        i++;
    }
    for (int p = 0; p< p_cg.size(); p++){
        i = 0;
        for(tie(vi,vf) = vertices(p_cg[p].g); vi != vf; ++vi) {
            p_cg[p].g[*vi].id = i;
            i++;
        }
    }
}
void MyGraph::reset_edge_index() {
    edge_iterator_t ei,ef;
    int i = 0;
    tie(ei,ef) = edges(g);
    for ( ; ei != ef ; ++ei){
        g[*ei].id = i;
        i++;
    }
    i = 0;
    tie(ei,ef) = edges(cg.g);
    for ( ; ei != ef ; ++ei){
        cg.g[*ei].id = i;
        i++;
    }
    for (int p = 0; p< p_cg.size(); p++){
        i = 0;
        tie(ei,ef) = edges(p_cg[p].g);
        for ( ; ei != ef ; ++ei){
            (p_cg[p].g)[*ei].id = i;
            i++;
        }
    }
}


void MyGraph::add2embedding(const CR& cross){
    if (cross->type == 'l') return;
    int v1 = (cross->nodes.first)->id;
    int v2 = (cross->nodes.second)->id;
    std::size_t num_existing_edges = embedding_storage[v1].size();
    if (num_existing_edges == 2) {
        if (cross->type == 'i' || cross->type == 's') {
            embedding_storage[v1].insert(embedding_storage[v1].begin() + 1, cross->edge.first);
            embedding_storage[v2].insert(embedding_storage[v2].begin() + 1, cross->edge.first);
        }
        else if (cross->type =='o'){
            embedding_storage[v1].insert(embedding_storage[v1].begin() + 2, cross->edge.first);
            embedding_storage[v2].insert(embedding_storage[v2].begin() + 2, cross->edge.first);
        }
    } else if (num_existing_edges == 3) {
        if (cross->type == 'i' || cross->type == 's') {
            //int pos = adjacent_edge_position(v1,edge,emb_storage[u],g);
            embedding_storage[v1].insert(embedding_storage[v1].begin() + 2, cross->edge.first);
            embedding_storage[v2].insert(embedding_storage[v2].begin() + 1, cross->edge.first);
        }
        /*
        else if (cross->type =='o'){
            embedding_storage[v1].insert(embedding_storage[v1].begin() + 3, cross->edge.first);
            embedding_storage[v2].insert(embedding_storage[v2].begin() + 2, cross->edge.first);
        }
         */
        else if (cross->type =='o') {
            int d0_otherV = circ_add(v1, 1, design->num_domains);
            int d1_otherV = circ_substract(v1, 1, design->num_domains);
            int c0_otherV = boost::source(embedding_storage[v1][2], g);
            if (c0_otherV == v1) c0_otherV = boost::target(embedding_storage[v1][2], g);
            int c1_otherV = boost::source(cross->edge.first, g);
            if (c1_otherV == v1) c1_otherV = boost::target(cross->edge.first, g);

            if (c1_otherV == c0_otherV){
                embedding_storage[v1].insert(embedding_storage[v1].begin() + 3, cross->edge.first);
            }
            else{
                cout << "MUHAHAHAHAHAHAHAHAHHAHAHAHAHAHHA-----------------------------------------------------------" << endl;
                int a = abs(circ_substract(d0_otherV, c0_otherV, design->num_domains));
                int b = abs(circ_substract(d1_otherV, c0_otherV, design->num_domains));
                if (a > b)
                    embedding_storage[v1].insert(embedding_storage[v1].begin() + 3, cross->edge.first);
                else
                    embedding_storage[v1].insert(embedding_storage[v1].begin() + 2, cross->edge.first);
            }

            d0_otherV = circ_add(v2, 1, design->num_domains);
            d1_otherV = circ_substract(v2, 1, design->num_domains);
            c0_otherV = boost::source(embedding_storage[v2][2], g);
            if (c0_otherV == v2) c0_otherV = boost::target(embedding_storage[v2][2], g);
            c1_otherV = boost::source(cross->edge.first, g);
            if (c1_otherV == v2) c1_otherV = boost::target(cross->edge.first, g);
            if (c1_otherV == c0_otherV){
                embedding_storage[v2].insert(embedding_storage[v2].begin() + 2, cross->edge.first);
            }
            else{
                cout << "MUHAHAHAHAHAHAHAHAHHAHAHAHAHAHHA-----------------------------------------------------------" << endl;
                int a = abs(circ_substract(d0_otherV, c0_otherV, design->num_domains));
                int b = abs(circ_substract(d1_otherV, c0_otherV, design->num_domains));
                if (a > b)
                    embedding_storage[v2].insert(embedding_storage[v2].begin() + 3, cross->edge.first);
                else
                    embedding_storage[v2].insert(embedding_storage[v2].begin() + 2, cross->edge.first);
            }
        }
    } else {
        std::cout << "Add Crossover:: vertex can only have 2 or 3 edges." << std::endl;
    }
}
void MyGraph::removeFromEmbedding(const CR& cross){
    int v1 = (cross->nodes.first)->id;
    int v2 = (cross->nodes.second)->id;
    edge_t edge = cross->edge.first;
    embedding_storage[v1].erase(std::remove(embedding_storage[v1].begin(),embedding_storage[v1].end(),edge),embedding_storage[v1].end());
    embedding_storage[v2].erase(std::remove(embedding_storage[v2].begin(),embedding_storage[v2].end(),edge),embedding_storage[v2].end());

}
double MyGraph::faces_weight() {
    double result = 0;
    double C_parameter = 2.8 * pow(10,-18);
    remove_long_crossovers(g, design);

    traversal_visitor outvis(g);
    set_face_data(g, embedding_storage, outvis);
    result += outvis.logsum;
    for (auto pool = design->staple_pools.begin(); pool!=design->staple_pools.end(); ++pool) {
        for (auto cross = pool->crossovers.begin(); cross!=pool->crossovers.end(); ++cross){
            if (cross->type == 'l' && cross->edge.second) {
                result += log(C_parameter / (cs_hack + long_pathweight(cross)));
            }
        }
    }
    readd_long_crossovers(g,design);

    return result;
}
double MyGraph::long_pathweight(const CR crossover) {
    DistanceMap distanceMap(&distances[0], get(vertex_index, g));
    dijkstra_shortest_paths_no_color_map(g, crossover->node_ids.first,
                                         weight_map(get(&EdgeProperty::weight, g)).
                                                 distance_map(distanceMap));
    return distanceMap[crossover->node_ids.second];
}


//Local Method
/*
double MyGraph::total_weight(const CR crossover) {
    DistanceMap distanceMap(&distances[0], get(vertex_index, g));
    dijkstra_shortest_paths_no_color_map(g, crossover->node_ids.first,
                                         weight_map(get(&EdgeProperty::weight, g)).
                                         distance_map(distanceMap));
    return distanceMap[crossover->node_ids.second];
}
 */
double MyGraph::shortest_path(const CR crossover) {
    DistanceMap distanceMap(&distances[0], get(vertex_index, g));
    edge_t edge;
    if (crossover->edge.second){
        edge = crossover->edge.first;
        boost::remove_edge(edge, g);
        removeFromEmbedding(crossover);
        reset_edge_index();
        dijkstra_shortest_paths_no_color_map(g, crossover->node_ids.first,
                                         weight_map(get(&EdgeProperty::weight, g)).
                                         distance_map(distanceMap));
        EdgeProperty EP;
        EP.id = boost::num_edges(g);
        EP.domain.second = false;
        EP.crossover = make_pair(crossover,true);
        EP.weight = cs_hack;
        crossover->edge = add_edge((crossover->nodes.first)->id, (crossover->nodes.second)->id, EP, g);
        add2embedding(crossover);
    }
    else{
        dijkstra_shortest_paths_no_color_map(g, crossover->node_ids.first,
                                             weight_map(get(&EdgeProperty::weight, g)).
                                                     distance_map(distanceMap));
    }
    //if (crossover->is_seam == true) cout << "SEAM ";
    //cout << distanceMap[crossover->node_ids.second] << endl;
    return distanceMap[crossover->node_ids.second];
}

//Component methods
void MyGraph::fill_components(){
    if (inputs->track_clusters){
        cg.fill_components();
        for (auto it = p_cg.begin(); it!= p_cg.end(); ++it){
            it->fill_components();
        }
    }
}
//Printers
void MyGraph::print_edges(){
    cout << "\n\n---------------------- Printing Edges -------------------------\n";
    edge_iterator_t ei,ef;
    tie( ei, ef) = edges(g);
    for ( ; ei != ef ; ++ei ){
        cout << g[*ei].id << "\t";
        cout << *ei << "\t";
        //cout << "Type: " << g[*ei].type << "\t";
        //cout << "Edge Length: " << g[*ei].length << "\t";
        if(g[*ei].domain_s.second == true)
            cout << "Domain Length: " << g[*ei].domain_s.first->length << "\t";
        if(g[*ei].crossover.second == true)
            cout << "Crossover Length: " << g[*ei].crossover.first->length << "\t";
        cout << "\n";
    }
}
void MyGraph::print_embedding() {
    face_counter count_visitor;
    planar_face_traversal(g, &embedding_storage[0], count_visitor, get(&EdgeProperty::id, g));
    cout << "Traversing faces of planar embedding... the planar embedding of the graph has " << count_visitor.count << " faces." << endl;
    cout << endl << "Vertices on the faces: " << endl;
    vertex_output_visitor v_vis;
    planar_face_traversal(g, &embedding_storage[0], v_vis, get(&EdgeProperty::id, g));
    cout << endl << "Edges on the faces: " << endl;
    edge_output_visitor e_vis;
    planar_face_traversal(g, &embedding_storage[0], e_vis, get(&EdgeProperty::id, g));
}
//Writers
void MyGraph::write_single(const string& filename, Graph& g_, const string& mychar) {
    reset_vertex_index();
    reset_edge_index();
    vertex_iterator_t vi,vf;
    edge_iterator_t ei,ef;
    tie( ei, ef) = edges(g_);
    for ( ; ei != ef ; ++ei ){
        if (g_[*ei].domain_s.second){
            if      ((g_[*ei].domain_s.first)->state == -1){g_[*ei].colour = "black"; g_[*ei].thickness = "2.0";}
            else if ((g_[*ei].domain_s.first)->state == 1){g_[*ei].colour = "red"; g_[*ei].thickness = "6.0";}
            else if ((g_[*ei].domain_s.first)->state == 0){g_[*ei].colour = "blue"; g_[*ei].thickness = "6.0";}
            else                                         {cout << "domain type not recognised!\n";}
        }
        else if (g_[*ei].domain_p.second){
            if      ((g_[*ei].domain_p.first)->state == false){g_[*ei].colour = "black"; g_[*ei].thickness = "2.0";}
            else if ((g_[*ei].domain_p.first)->pool_id == 1){g_[*ei].colour = "red"; g_[*ei].thickness = "6.0";}
            else if ((g_[*ei].domain_p.first)->pool_id == 0){g_[*ei].colour = "blue"; g_[*ei].thickness = "6.0";}
            else                                         {cout << "domain type not recognised!\n";}
        }
        else if (g_[*ei].crossover.second){
            if      ((g_[*ei].crossover.first)->is_seam == false)    {g_[*ei].colour = "green"; g_[*ei].thickness = "3.0";}
            else                                                    {g_[*ei].colour = "purple"; g_[*ei].thickness = "3.0";}
        }
        else {cout << "Error while writing graph: Not domain nor crossover!\n";}
    }
    boost::dynamic_properties dp;
    //dp.property("len", get(&EdgeProperty::length, g_));
    //dp.property("type", get(&EdgeProperty::type, g_));
    dp.property("color", get(&EdgeProperty::colour, g_));
    dp.property("penwidth", get(&EdgeProperty::thickness, g_));
    dp.property("node_id", get(&VertexProperty::id, g_)); //must have!
    dp.property("pos", get(&VertexProperty::pos, g_));
    dp.property("height", get(&VertexProperty::height, g_));
    dp.property("width", get(&VertexProperty::width, g_));
    dp.property("fixedsize", get(&VertexProperty::other, g_));
    ofstream outfile;
    string str= filename + "_";   str+= mychar; str+= ".gv";
    outfile.open(str);
    write_graphviz_dp(outfile, g_, dp);
    outfile.close();
    //str= filename + "_";   str+= mychar; str+= ".gml";
    //outfile.open(str);
    //write_graphml(outfile, g_, dp, true);
    //outfile.close();
}
void MyGraph::write(const string& filename) {
    //write_single(filename,g,"M");
    //write_single(filename,cg.g,"C");
    //if (this->design->num_staple_pools == 1) write_single(filename,this->g,design->staple_pools[0].name);
    string str = "";
    for (const auto& pool : this->design->staple_pools){
        str+= pool.name;
    }
    write_single(filename,this->g,str);
    for (auto it = p_cg.begin(); it!= p_cg.end(); ++it){
        write_single(filename,it->g,design->staple_pools[it-p_cg.begin()].name);
    }
}










//Redundant:
void MyGraph::complete(){//Make full graph
    for (auto dom = design->domains.begin(); dom!=design->domains.end(); ++dom){
        bind_domain(dom);
    }
    for (auto pool = design->staple_pools.begin(); pool!= design->staple_pools.end(); ++pool){
        for (auto cr = pool->crossovers.begin(); cr!=pool->crossovers.end(); ++cr){
            if (cr->type!='l'){
                add_crossover(cr);
            }
        }
    }
}
edge_t MyGraph::id_to_edge(int i){
    edge_t result;
    edge_iterator_t ei,ef;
    tie( ei, ef) = edges(g);
    for ( ; ei != ef ; ++ei ){
        if (g[*ei].id == i){
            result = *ei;
            break;
        }
    }
    return result;
}

//Redundant: Only used in complete()
void MyGraph::bind_domain(const SDOM domain, bool pseudo){
    g[domain->edge].weight = domain->length * domain->length * ds_hack;
    //g[domain->edge].type = 'd';
    //New
    EdgeProperty EP;
    EP.crossover.second = false;
    EP.id = domain->id;
    //EP.domain = make_pair(domain,true);
    EP.domain_s = make_pair(domain,true);
    domain->edge2 = add_edge(domain->node_ids.first, domain->node_ids.second, EP, cg.g);
    reset_vertex_index();
    reset_edge_index();
}
void MyGraph::unbind_domain(const SDOM domain, bool pseudo){
    g[domain->edge].weight = domain->length * ss_hack;
    //g[domain->edge].type = 's';
    //New
    remove_edge(domain->edge2.first, cg.g);
    domain->edge2.second = false;
    reset_vertex_index();
    reset_edge_index();
}

