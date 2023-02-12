#include "../src/Constants.h"
#include "../src/Design.h"
#include "../src/MyGraph.h"
#include "../src/TempRamp.h"
#include "../src/Simulation.h"

void manual_apply_normal(TransitionManager *localManager, Uni& uni, MyGraph *G, TempRamp *ramp,
                         int pool_id, int domain_id, State_t initial_state, State_t final_state){
    int staple_id = localManager->design->staple_pools[pool_id].domains[domain_id].staple_id;
    string movie_file_name;
    localManager->fill_rates();
    localManager->select_transition(uni);
    cout << localManager->total_rate << "\n";
    localManager->manual_select_normal(pool_id, staple_id, domain_id, initial_state, final_state)->print();
    localManager->manual_select_normal(pool_id, staple_id, domain_id, initial_state, final_state)->apply(G);
    //if(!localManager->local) G->update_faces();
    localManager->step++;
    ramp->move_time(localManager->tau);
    localManager->reset_possibles();
    localManager->reset_recalculates();
    movie_file_name = "Output/Movie/";
    movie_file_name += to_string(localManager->step);
    movie_file_name += "_";
    movie_file_name += to_string_with_precision(centigrade(ramp->get_T()));
    G->write(movie_file_name);
}

void manual_apply_invasion(TransitionManager *localManager, Uni& uni, MyGraph *G, TempRamp *ramp,
                           int pool_id1, int domain_id1, State_t initial_state1, State_t final_state1,
                           int pool_id2, int domain_id2, State_t initial_state2, State_t final_state2){
    int staple_id1 = localManager->design->staple_pools[pool_id1].domains[domain_id1].staple_id;
    int staple_id2 = localManager->design->staple_pools[pool_id2].domains[domain_id2].staple_id;
    string movie_file_name;
    localManager->fill_rates();
    localManager->select_transition(uni);
    cout << localManager->total_rate << "\n";
    localManager->manual_select_invasion(pool_id1, staple_id1, domain_id1, initial_state1, final_state1,
                                         pool_id2, staple_id2, domain_id2, initial_state2, final_state2)->print();
    localManager->manual_select_invasion(pool_id1, staple_id1, domain_id1, initial_state1, final_state1,
                                         pool_id2, staple_id2, domain_id2, initial_state2, final_state2)->apply(G);
    //if(!localManager->local) G->update_faces();
    localManager->step++;
    ramp->move_time(localManager->tau);
    localManager->reset_possibles();
    localManager->reset_recalculates();
    movie_file_name = "Output/Movie/";
    movie_file_name += to_string(localManager->step);
    movie_file_name += "_";
    movie_file_name += to_string_with_precision(centigrade(ramp->get_T()));
    G->write(movie_file_name);
}

bool pairCompare(const std::pair<int, int*>& firstElem, const std::pair<int, int*>& secondElem) {
    return firstElem.first < secondElem.first;
}

void test_loops(Constants* constants, MyGraph *G,  Design *design){
    //cout << constants->C_parameter << "\t" << log(constants->C_parameter) << endl;
    design->staple_pools[0].print_crossovers();

    double W = 0;
    double Gshape = G->faces_weight();
    double Gshape2 = G->faces_weight();
    /*
    for (auto loop = G->Faces.faces.begin(); loop != G->Faces.faces.end(); ++loop){
        loop->calculate_weight(G->g);
        W+= loop->weight;
        Gshape+= log(constants->C_parameter / loop->weight);
    }*/
    //G->Faces.print();
    cout << "Total W = " << W << endl;
    cout << "Gshape = " << Gshape << endl;

    double Wlast, Glast, Wlocal;
    for (auto cross = design->staple_pools[0].crossovers.begin();
         cross != design->staple_pools[0].crossovers.end(); ++cross){
        std::cout << "--------------------------------------------------------------- adding cr ";
        std::cout << cross->id << std::endl;
        Glast = Gshape;
        G->bind_domain(cross->domains.first);
        G->bind_domain(cross->domains.second);
        G->add_crossover(cross);
        Gshape = G->faces_weight();
        Wlocal = G->shortest_path(cross);
        cout << "Gshape = " << Gshape << endl;
        cout << "Delta Gshape = " << Gshape - Glast << endl;
        cout << "Delta Glocal = " << log(constants->C_parameter / Wlocal) << endl;
    }

    for (auto cross = design->staple_pools[0].crossovers.begin();
         cross != design->staple_pools[0].crossovers.end(); ++cross){
        std::cout << "--------------------------------------------------------------- removing cr ";
        std::cout << cross->id << std::endl;
        Wlast = W;
        Glast = Gshape;
        G->unbind_domain(cross->domains.first);
        G->unbind_domain(cross->domains.second);
        G->remove_crossover(cross);
        Gshape = G->faces_weight();
        Wlocal = G->shortest_path(cross);
        cout << "Gshape = " << Gshape << endl;
        cout << "Delta Gshape = " << Gshape - Glast << endl;
        cout << "Delta Glocal = " << log(constants->C_parameter / Wlocal) << endl;
    }
}

void test_full(Simulation* sim){
    auto trManager = sim->trManager;

    /*
    std::cout << trManager->transitions.size() << endl;

    TransitionMap M;
    State_names state_names;
    for (auto &tp : M.FullMap){
        if (tp.second.allowed){
            std::cout << state_names.array(tp.first.first) << "\t";
            std::cout << state_names.array(tp.first.second) << std::endl;
        }
    }
    trManager->reset_possibles();
    trManager->fill_rates();


    for (auto &tpair : trManager->transition_pairs){
        std::cout << *(tpair.first) << "\t" << *(tpair.second) << std::endl;
    }

    for (auto &tr : trManager->transitions){
        tr.decide_possible();
        if (tr.possible){
            tr.print();
        }
    }
     */



    Design* design = trManager->design;
    StaplePool* pool = &design->staple_pools[0];
    vector<Transition*> path, reverse_path;
    for (auto &st: pool->staples){
        for (auto &trans : st.transitions){
            if (st.num_domains == 1 && trans->initial_state == State_t(s0) && trans->final_state == State_t(s1))
                path.push_back(trans);
            if (st.num_domains == 2 && trans->initial_state == State_t(s00) && trans->final_state == State_t(s01))
                path.push_back(trans);
            if (st.num_domains == 3 && trans->initial_state == State_t(s000) && trans->final_state == State_t(s001))
                path.push_back(trans);
        }
    }

    for (auto &st: pool->staples){
        for (auto &trans : st.transitions){
            if (st.num_domains == 2 && trans->initial_state == State_t(s01) && trans->final_state == State_t(s11))
                path.push_back(trans);
            if (st.num_domains == 3 && st.id <= 40 && trans->initial_state == State_t(s001) && trans->final_state == State_t(s101))
                path.push_back(trans);
            if (st.num_domains == 3 && st.id > 40 && trans->initial_state == State_t(s001) && trans->final_state == State_t(s011))
                path.push_back(trans);
        }
    }
    for (auto &st: pool->staples){
        for (auto &trans : st.transitions){
            if (st.num_domains == 3 && st.id > 40 && trans->initial_state == State_t(s011) && trans->final_state == State_t(s111))
                path.push_back(trans);
        }
    }
    for (auto &trans : path){
        reverse_path.push_back(&(*(trans->reverse)));
    }

    std::cout << "--------------------------------------------------------------- ";
    std::cout << trManager->G->faces_weight();
    std::cout << std::endl;

    double k_bi = trManager->inputs->k_bi;
    double k_uni = trManager->inputs->k_uni;
    double RT = gas_constant * 273.15+50;

    for (auto &tr : path){
        trManager->reset_possibles();
        trManager->next = trManager->transitions.begin() + tr->id;

        trManager->set_dG_duplex(trManager->next);
        trManager->set_dG_stack(trManager->next);
        //trManager->set_dG_shape(trManager->next);
        trManager->next->dG_shape = 0;

        if (tr->uni){
            if (tr->forward){
                trManager->set_dG_shape(trManager->next);
                tr->dG = tr->dG_shape;
                tr->rate = k_uni * exp( -tr->dG / RT );
            }
            else{
                tr->dG = tr->dG_duplex+tr->dG_stack;
                tr->rate = k_uni * exp( tr->dG / RT );
            }
        }
        else{
            if (tr->forward){
                tr->rate = k_bi * tr->staple->concentration;
            }
            else{
                if (!trManager->local){trManager->set_dG_shape(trManager->next);}
                tr->dG = tr->dG_duplex+tr->dG_stack-tr->dG_shape;
                tr->rate = k_uni * exp( tr->dG / RT );
            }
        }

        trManager->total_rate = 100;
        trManager->tau = 1/trManager->total_rate;
        trManager->next->tau = trManager->tau;

        std::cout << *trManager->next << "\t";
        std::cout << trManager->next->dG_duplex << "\t";
        std::cout << trManager->next->dG_stack << "\t";
        std::cout << trManager->next->dG_shape << "\t";
        if (trManager->next->properties.cross)
            std::cout << trManager->next->crossover.first->type;
        std::cout << std::endl;
        trManager->apply_next();
    }

    std::cout << "--------------------------------------------------------------- ";
    std::cout << trManager->G->faces_weight();
    std::cout << std::endl;

    for (auto &tr : reverse_path){
        trManager->reset_possibles();
        trManager->next = trManager->transitions.begin() + tr->id;

        trManager->set_dG_duplex(trManager->next);
        trManager->set_dG_stack(trManager->next);
        trManager->set_dG_shape(trManager->next);
        //trManager->next->dG_shape = 0;

        if (tr->uni){
            if (tr->forward){
                trManager->set_dG_shape(trManager->next);
                tr->dG = tr->dG_shape;
                tr->rate = k_uni * exp( -tr->dG / RT );
            }
            else{
                tr->dG = tr->dG_duplex+tr->dG_stack;
                tr->rate = k_uni * exp( tr->dG / RT );
            }
        }
        else{
            if (tr->forward){
                tr->rate = k_bi * tr->staple->concentration;
            }
            else{
                if (!trManager->local){trManager->set_dG_shape(trManager->next);}
                tr->dG = tr->dG_duplex+tr->dG_stack-tr->dG_shape;
                tr->rate = k_uni * exp( tr->dG / RT );
            }
        }

        trManager->total_rate = 100;
        trManager->tau = 1/trManager->total_rate;
        trManager->next->tau = trManager->tau;

        std::cout << *trManager->next << "\t";
        std::cout << trManager->next->dG_duplex << "\t";
        std::cout << trManager->next->dG_stack << "\t";
        std::cout << trManager->next->dG_shape << "\t";
        if (trManager->next->properties.cross)
            std::cout << trManager->next->crossover.first->type;
        std::cout << std::endl;
        trManager->apply_next();
    }

    std::cout << "--------------------------------------------------------------- ";
    std::cout << trManager->G->faces_weight();
    std::cout << std::endl;
}

int main(int argc, char * argv[]) {
    Inputs *inputs = new Inputs(argc,argv);
    Constants *constants = new Constants(inputs);
    Design *design = new Design(inputs);
    MyGraph *G = new MyGraph(design);
    TempRamp *ramp = new TempRamp(inputs);
    FileIO *ofiles = new FileIO(design);
    TransitionManager *trManager = new TransitionManager(constants, G, ramp, ofiles);
    StatManager *statManager = new StatManager(trManager);
    OPManager *opManager = new OPManager(statManager);
    Simulation *sim = new Simulation(opManager);
    if (inputs->test){
        //do_test();
        //design->print_domains();
        //design->print_OPs();
        //for (auto pool = design->staple_pools.begin(); pool!= design->staple_pools.end(); ++pool){
        //cout << pool->name << endl;
        //pool->print_domains();
        //for (auto st = pool->staples.begin(); st!=pool->staples.end(); ++st){
        //st->print();
        //}
        //pool->print_staples();
        //pool->print_helices();
        //pool->print_crossovers();
        //}

        //G->print_embedding();

        //test_loops(constants,G,design);
        test_full(sim);
        sim->ofiles->close_files();
    }
    else if (inputs->config_generator){
        sim->run();
    }
    else if (inputs->anneal || inputs->melt || inputs->isothermal){
        sim->run();
    }
    else{
        printf ("Please select sim type: test, anneal, melt, isothermal, config_generator.\n");
        exit (EXIT_FAILURE);
    }
    delete inputs;
    delete ofiles;
    delete constants;
    delete design;
    delete G;
    delete ramp;
    delete trManager;
    delete sim;
    return 0;
}
