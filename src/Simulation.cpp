/*
 *
 *	Simulation.cpp
 * 	Author: Behnam
 *
 */

#include "Simulation.h"

Simulation::Simulation(Inputs* inputs, Constants* constants) :
inputs(inputs), constants(constants){
    design = new Design(inputs);
    G = new MyGraph(design);
    ramp = new TempRamp(inputs);
    ofiles = new FileIO(design);
    trManager = new TransitionManager(constants, G, ramp, ofiles);
    statManager = new StatManager(trManager);
    opManager = new OPManager(statManager);
}

Simulation::~Simulation() {
    delete ofiles;
    delete design;
    delete G;
    delete ramp;
    delete trManager;
    delete statManager;
    delete opManager;
}

void Simulation::out_AM(){
    if (ramp->T_was_changed){
        cout << ramp->Temps[ramp->prev_idx].str_T << "\t";
        cout << ramp->Temps[ramp->prev_idx].num_tr[0] << "\t";
        cout << ramp->Temps[ramp->prev_idx].num_rejected << "\t";
        for (const auto& pool : design->staple_pools){
            cout << pool.OPs[0].stats[ramp->prev_idx].rolling_average << "\t";
        }
        cout << "\n";
    }
}
void Simulation::out_Iso(){
    if (trManager->step % 10000 == 0){
        cout << trManager->step << "\t";
        for (const auto& pool : design->staple_pools){
            cout << pool.OPs[1].state << "\t";
        }
        cout << "\n";
    }
}
void Simulation::write_movie(){
    str = "Output/Movie/";
    str += to_string(trManager->step);
    str += "_";
    str += to_string_with_precision(centigrade(ramp->get_T()));
    G->write(str);
}

void Simulation::prepare_config(char type){
    /* options:
     * 'r': randomly select initial config
     * 'l': select lowest possible opVal
     * 'h': select highest possible opVal
     * 'm': select middle opVal
    */

    int current_OPval, target_OPval;
    current_OPval = this->opManager->biased.second->state;
    vector<int> possible_OPvals;
    for (const auto &it: this->opManager->biased.second->weight) {
        possible_OPvals.push_back(it.first);
    }
    typedef boost::variate_generator<base_generator_type &, uniform_int<> > UniformInt_gen;
    base_generator_type generator(inputs->seed);
    int randIdx;

    if (type == 'r') {
        uniform_int<> uniformInt(0, possible_OPvals.size() - 1);
        UniformInt_gen uni(generator, uniformInt);
        randIdx = uni();
        target_OPval = possible_OPvals[randIdx];
    }
    else if (type == 'l') {
        target_OPval = *possible_OPvals.begin();
    }
    else if (type == 'h') {
        target_OPval = *possible_OPvals.rbegin();
    }
    else if (type == 'm') {
        target_OPval = (*possible_OPvals.begin() + *possible_OPvals.rbegin()) / 2;
    }
    else {
        cout << "Error: Select a type of initial configuration for Simulation::generate_config().\n";
        exit (EXIT_FAILURE);
    }

    vector<Transition*> possible_trs;
    while (target_OPval != current_OPval){
        //cout << current_OPval << "\t" << target_OPval << endl;
        trManager->reset_possibles();
        possible_trs.clear();
        for (auto& tr : this->trManager->transitions){
            if (tr.possible){
                if (target_OPval > current_OPval){ //binds (normal) and forwards (invasion)
                    if (tr.type == normal){
                        if (tr.pool == opManager->biased.second->pool_id){
                            if (tr.properties.bind) possible_trs.push_back(&tr);
                        }
                        else{
                            possible_trs.push_back(&tr);
                        }
                    }
                    else{
                        if (tr.forward) possible_trs.push_back(&tr);
                    }
                }
                else{ //unbinds (normal) and backwards (invasion)
                    if (tr.type == normal){
                        if (tr.pool == opManager->biased.second->pool_id){
                            if (!tr.properties.bind) possible_trs.push_back(&tr);
                        }
                        else{
                            possible_trs.push_back(&tr);
                        }
                    }
                    else{
                        if (!tr.forward) possible_trs.push_back(&tr);
                    }
                }
            }
        }
        uniform_int<> uniformInt2(0, possible_trs.size()-1);
        UniformInt_gen uni2(generator,uniformInt2);
        randIdx = uni2();
        possible_trs[randIdx]->apply(G);
        G->fill_components();
        current_OPval = this->opManager->biased.second->state;
    }
    trManager->reset_possibles();
}


void Simulation::print_debug() {
    std::cout << trManager->step << ":\n";
    for (auto st : design->staple_pools[0].staples){
        cout << st.state << "\t";
    }   cout << "\n";
    for (auto dom : design->staple_pools[0].domains){
        cout << dom.state << "\t";
    }   cout << "\n";
    for (auto cr : design->staple_pools[0].crossovers){
        cout << cr.state << "\t";
    }   cout << "\n";
    for (auto op : design->staple_pools[0].OPs){
        cout << op << "\n";
    }
    std::cout << G->embedding_storage << "\n";
}


void Simulation::run(){
    base_generator_type generator(inputs->seed);
    uniform_real<> uni_dist(0,1);
    Uni uni(generator, uni_dist);
	if (inputs->anneal || inputs->melt) {
        cout << "------ Starting Simulation ------\n";
		while (ramp->current_t < ramp->t_max){
            trManager->append_trajectory();
			trManager->fill_rates();
			trManager->select_transition(uni);
            trManager->write_transition();
            ramp->move_time(trManager->tau);
            statManager->update_times();
            statManager->update_counts();
            trManager->apply_next();
            statManager->write_all();

            this->out_AM();
            if (inputs->make_movie){
                this->write_movie();
            }
			trManager->reset_possibles();
			trManager->reset_recalculates();
		}
        cout << "------ Ending Simulation ------\n";
        ofiles->info_file << "N_transitions: " << trManager->step << "\n";
        opManager->write_last();
        statManager->write_in_times();
	}
	else if (inputs->isothermal) {
        cout << "------ Starting Iso Simulation ------\n";
        if (inputs->umbrella_sampling && inputs->initialise_as == "random") this->prepare_config();
        while (trManager->step < inputs->max_steps && ramp->current_t < ramp->t_max) {
            if (!design->target_reached){
                if (design->staple_pools[inputs->target_pool].OPs[2].state == design->staple_pools[inputs->target_pool].num_crossovers){
                    design->target_reached = true;
                    design->target_reached_at.first = trManager->step;
                    design->target_reached_at.second = ramp->current_t;
                    ofiles->info_file << "Target reached at: " << design->target_reached_at.first << "\t" << design->target_reached_at.second << "\n";
                }
            }
			trManager->fill_rates();
            if (inputs->umbrella_sampling)  opManager->fill_rates_w();
			trManager->select_transition(uni);
            trManager->write_transition();
            trManager->append_trajectory();
            statManager->write_all();
            ramp->move_time(trManager->tau);
            statManager->update_times();
            statManager->update_counts();
            opManager->update_times();
            if (trManager->step == inputs->burn_steps)
                opManager->write_burn();
            if ((trManager->step % inputs->write_hist_every == 0 && trManager->step > 0)
                || trManager->step == design->target_reached_at.first){
                ofiles->retrunc_hist_files();
                opManager->write();
                opManager->write_last();
                opManager->write_object_hist();
            }
            trManager->apply_next();
            opManager->set_values();

            this->out_Iso();
            if (inputs->make_movie) this->write_movie();
            trManager->reset_possibles();
            trManager->reset_recalculates();
            if(!inputs->umbrella_sampling && design->target_reached && statManager->all_in_times_done) {
                cout << "Target Reached!" << endl;
                break;
            }
		}
        statManager->write_all();
        trManager->append_trajectory();
        ofiles->info_file << "N_transitions: " << trManager->step << "\n";
        ofiles->retrunc_hist_files();
        opManager->write_last();
        opManager->write_object_hist();
        statManager->write_in_times();
        cout << "------ Ending Iso Simulation ------\n";
	}
    else if (inputs->weight_generator) {
        while (trManager->step < inputs->max_steps) {
            this->out_Iso();
            //this->print_debug();
            trManager->fill_rates(); 
            opManager->fill_rates_w(); 
            trManager->select_transition(uni); 
            ramp->move_time(trManager->tau); 
            statManager->update_times(); 
            statManager->update_counts(); 
            opManager->update_times(); 
            trManager->apply_next(); 
            opManager->set_values(); 
            trManager->reset_possibles(); 
            trManager->reset_recalculates();
        }
    }

    else if (inputs->config_generator) {
        cout << "------ Starting Configuration Generator ------\n";
        while (trManager->step < inputs->max_steps && !design->target_reached) {
            trManager->fill_rates();
            if (inputs->umbrella_sampling)  opManager->fill_rates_w();
            trManager->select_transition(uni);
            trManager->apply_next();
            opManager->set_values();
            ramp->move_time(trManager->tau);
            trManager->reset_possibles();
            trManager->reset_recalculates();
            opManager->generate_config();
        }
        cout << "------ Ending Config Generator ------\n";
    }
    else{
        printf ("Please select anneal, melt or isothermal. \n");
        exit (EXIT_FAILURE);
    }
    ofiles->close_files();
}

void Simulation::test(){
    for (auto pool=design->staple_pools.begin(); pool!= design->staple_pools.end(); ++pool){
        for (auto cr1 = pool->crossovers.begin(); cr1!= pool->crossovers.end(); ++cr1){
            for (auto cr2 = pool->crossovers.begin(); cr2!= pool->crossovers.end(); ++cr2){
                if (cr1->node_ids.first == cr2->node_ids.first && cr1->node_ids.second == cr2->node_ids.second && cr1->id != cr2->id)
                    cout << "A\n";
                if (cr1->node_ids.first == cr2->node_ids.second && cr1->node_ids.second == cr2->node_ids.first && cr1->id != cr2->id)
                    cout << "B\n";
                if (cr1->node_ids.second == cr2->node_ids.first && cr1->node_ids.first == cr2->node_ids.second && cr1->id != cr2->id)
                    cout << "C\n";
            }
        }
    }
}
