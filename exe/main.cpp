#include "../src/Constants.h"
#include "../src/Design.h"
#include "../src/MyGraph.h"
#include "../src/TempRamp.h"
#include "../src/Simulation.h"
#include "../src/Exact.h"


void test_loops(Constants* constants, MyGraph *G,  Design *design){
    //cout << constants->C_parameter << "\t" << log(constants->C_parameter) << endl;
    design->staple_pools[0].print_crossovers();

    double Wshape = G->faces_weight();
    double Gshape = ( (G->numFaces + G->numLong) * log(constants->C_parameter) ) + Wshape;
    /*
    for (auto loop = G->Faces.faces.begin(); loop != G->Faces.faces.end(); ++loop){
        loop->calculate_weight(G->g);
        W+= loop->weight;
        Gshape+= log(constants->C_parameter / loop->weight);
    }*/
    //G->Faces.print();
    cout << "Total W = " << Wshape << endl;
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
        Wshape = G->faces_weight();
        Gshape = ( (G->numFaces + G->numLong) * log(constants->C_parameter) ) + Wshape;
        Wlocal = G->shortest_path(cross);
        cout << "Gshape = " << Gshape << endl;
        cout << "Delta Gshape = " << Gshape - Glast << endl;
        cout << "Delta Glocal = " << log(constants->C_parameter / Wlocal) << endl;
    }

    for (auto cross = design->staple_pools[0].crossovers.begin();
         cross != design->staple_pools[0].crossovers.end(); ++cross){
        std::cout << "--------------------------------------------------------------- removing cr ";
        std::cout << cross->id << std::endl;
        Wlast = Wshape;
        Glast = Gshape;
        G->unbind_domain(cross->domains.first);
        G->unbind_domain(cross->domains.second);
        G->remove_crossover(cross);
        Wshape = G->faces_weight();
        Gshape = ( (G->numFaces + G->numLong) * log(constants->C_parameter) ) + Wshape;
        Wlocal = G->shortest_path(cross);
        cout << "Gshape = " << Gshape << endl;
        cout << "Delta Gshape = " << Gshape - Glast << endl;
        cout << "Delta Glocal = " << log(constants->C_parameter / Wlocal) << endl;
    }
}

void test_full(Simulation* sim){
    auto trManager = sim->trManager;
    MyGraph *G = trManager->G;

    print_embedding_storage(sim->G->g, sim->G->embedding_storage);

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
            if (st.num_domains == 3 && st.id <= 0 && trans->initial_state == State_t(s001) && trans->final_state == State_t(s101))
                path.push_back(trans);
            if (st.num_domains == 3 && st.id > 0 && trans->initial_state == State_t(s001) && trans->final_state == State_t(s011))
                path.push_back(trans);
        }
    }
    for (auto &st: pool->staples){
        for (auto &trans : st.transitions){
            if (st.num_domains == 3 && st.id > 0 && trans->initial_state == State_t(s011) && trans->final_state == State_t(s111))
                path.push_back(trans);
        }
    }
    for (auto tr = path.rbegin(); tr!=path.rend(); ++tr){
        reverse_path.push_back(&(*((*tr)->reverse)));
    }

    std::cout << "--------------------------------------------------------------- Binding";
    if (!trManager->local){
        double Wshape, Gshape;
        Wshape = G->faces_weight();
        Gshape = ( (G->numFaces + G->numLong) * log(trManager->constants->C_parameter) ) + Wshape;
        std::cout << Wshape << "\t" << Gshape;
    }
    std::cout << std::endl;

    double k_bi = trManager->inputs->k_bi;
    double k_uni = trManager->inputs->k_uni;
    double RT = gas_constant * trManager->ramp->get_T();

    for (auto &tr : path){
        trManager->reset_possibles();
        trManager->next = trManager->transitions.begin() + tr->id;

        trManager->set_dG_duplex(trManager->next);
        trManager->set_dG_stack(trManager->next);
        trManager->set_dG_shape(trManager->next);
        //trManager->next->dG_shape = 0;

        if (tr->uni){
            if (tr->forward){
                //trManager->set_dG_shape(trManager->next);
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
                //if (!trManager->local){trManager->set_dG_shape(trManager->next);}
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

    std::cout << "--------------------------------------------------------------- Unbinding";
    if (!trManager->local) {
        double Wshape, Gshape;
        Wshape = G->faces_weight();
        Gshape = ( (G->numFaces + G->numLong) * log(trManager->constants->C_parameter) ) + Wshape;
        std::cout << Wshape << "\t" << Gshape;
    }
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
                //trManager->set_dG_shape(trManager->next);
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
                //if (!trManager->local){trManager->set_dG_shape(trManager->next);}
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
    if (!trManager->local) {
        double Wshape, Gshape;
        Wshape = G->faces_weight();
        Gshape = ( (G->numFaces + G->numLong) * log(trManager->constants->C_parameter) ) + Wshape;
        std::cout << Wshape << "\t" << Gshape;
    }
    std::cout << std::endl;
}

void test_random(){
    int myidx;

    int inputrandom = 12341;
    base_generator_type generator(12341);
    typedef boost::variate_generator<base_generator_type&, uniform_int<> > UniformInt_gen;

    std::vector<int> myvec = {3,5,647,4352,7,4,2};
    uniform_int<> uniformInt(0, myvec.size()-1);
    UniformInt_gen uni(generator,uniformInt);

    for (int i = 0; i < 100; i++){
        myidx = uni();
        std::cout << i << "\t" << myidx << "\t" << myvec[myidx] << std::endl;
    }

    std::vector<int> myvec2 = {3,5,647,4352,7,4,2,6,12,65,2,1};
    uniform_int<> uniformInt2(0, myvec2.size()-1);
    UniformInt_gen uni2(generator,uniformInt2);

    for (int i = 0; i < 100; i++){
        myidx = uni2();
        std::cout << i << "\t" << myidx << "\t" << myvec2[myidx] << std::endl;
    }
}

void test_ramp(){
    Inputs *inputs = new Inputs();
    inputs->sim_type = "anneal";
    inputs->anneal = true;
    inputs->isothermal = inputs->exact = false;
    TempRamp *ramp = new TempRamp(inputs);


    std::cout << *ramp << std::endl;
    std::cout << ramp->T_was_changed << std::endl;
    ramp->print();
}

bool has_zero_value(const std::map<int, int>& my_map) {
    for (const auto& entry : my_map) {
        if (entry.second == 0) {
            return true;
        }
    }
    return false;
}


int main(int argc, char * argv[]) {
    Inputs inputs = Inputs(argc,argv);
    Constants constants = Constants(&inputs);
    if (inputs.test){
        Design design = Design(&inputs);
        //do_test();
        design.print_OPs();
        for (auto &pool : design.staple_pools){
            pool.print_OPs();
        }

        //G->print_embedding();

        //test_loops(constants,G,design);
        //test_full(sim);
        //double Tm = calculate_Tm(sim);
        //test_random();
        //test_ramp();
    }
    else if (inputs.exact){
        Simulation sim = Simulation(&inputs, &constants);
        double Tm = calculate_Tm(&sim);
        //std::cout << Tm << std::endl;
        //exact(&sim);
        //Er2_legacy(sim);
        sim.ofiles->close_files();
    }
    else if (inputs.anneal || inputs.melt || inputs.isothermal || inputs.config_generator){
        Simulation sim = Simulation(&inputs, &constants);
        sim.run();
    }
    else if (inputs.weight_generator){
        std::map<int,long double> times, weights;
        std::map<int,int> counts;
        std::string og_w_file_name = inputs.w_file_name;
        std::size_t dot_position = og_w_file_name.find('.');
        std::string name = og_w_file_name.substr(0, dot_position);
        std::string extension = og_w_file_name.substr(dot_position + 1);
        struct stat st;
        string idummy, dummy;
        string command = "mv Weights ";
        if (stat("Weights", &st) == 0) {
            for (int i = 1; i < 100; i++) {
                idummy = to_string(i);
                dummy = "Weights" + idummy;
                if (!(stat(dummy.c_str(), &st) == 0)) {
                    break;
                }
            }
            command += dummy;
            std::system(command.c_str());
            cout << "Renamed old Weights to " << dummy << "\n";
        }
        std::system("mkdir Weights");
        bool finished = false;
        int repeat = 0;
        long double StdOverMean = 9999999;
        ofstream info_file;
        std::string info_file_name = "Weights/Info.txt";
        info_file.open(info_file_name,std::ofstream::out | std::ofstream::trunc);
        info_file << "Repeat" << "\t" << "Seed" << "\t" << "wfile" << "\t";
        info_file << "Mean" << "\t" << "Variance" << "\t" << "Std" << "\t" << "Std/Mean";
        info_file << std::endl;
        std::cout << "------------------------------------------------ Beginning Generation Steps " << std::endl;
        while (!finished){
            std::cout << "------------------------ Weight Generation Step " << repeat << " of " << inputs.num_repeats << std::endl;
            inputs.seed = std::chrono::system_clock::now().time_since_epoch().count();
            info_file << repeat << "\t" << inputs.seed << "\t" << inputs.w_file_name << "\t";
            Simulation sim = Simulation(&inputs, &constants);
            sim.run();
            std::string newname = std::to_string(repeat) + "_" + name;
            std::string new_w_file_name = "Weights/"+ newname + "." + extension;
            std::string new_hist_file_name = "Weights/"+std::to_string(repeat)+"_Hist.txt";
            sim.opManager->write_new_weights(new_w_file_name);
            sim.opManager->write_weight_gen_hist(new_hist_file_name);
            // Calculate the mean, variance, and standard deviation of the w_new values
            long double sum = 0;
            long double sum_of_squares = 0;
            int count = 0;
            for (const auto &entry : sim.opManager->biased.second->stats[0].time) {
                sum += entry.second;
                sum_of_squares += entry.second * entry.second;
                count++;
            }
            long double mean = sum / count;
            long double variance = (sum_of_squares / count) - (mean * mean);
            long double std_dev = std::sqrt(variance);
            StdOverMean = std_dev / mean;
            info_file << mean << "\t" << variance << "\t" << std_dev << "\t" << StdOverMean << std::endl;

            inputs.w_file_name = new_w_file_name;
            repeat++;
            if (repeat == inputs.num_repeats || ( StdOverMean < 1 && !has_zero_value(sim.opManager->biased.second->stats[0].count) ) ){
                finished = true;
                times = sim.opManager->biased.second->stats[0].time;
                counts = sim.opManager->biased.second->stats[0].count;
                weights = sim.opManager->biased.second->weight;
            }
        }
        std::cout << "------------------------------------------------ Beginning Final Steps " << std::endl;
        info_file << "x" << inputs.num_repeats << "\t" << "NA" << "\t" << inputs.w_file_name << "\t";
        repeat = 0;
        while (repeat < inputs.num_repeats){
            std::cout << "------------------------ Final Step " << repeat << " of " << inputs.num_repeats << std::endl;
            inputs.seed = std::chrono::system_clock::now().time_since_epoch().count();
            Simulation sim = Simulation(&inputs, &constants);
            sim.run();
            for (auto& entry : times){
                entry.second += sim.opManager->biased.second->stats[0].time[entry.first];
                counts[entry.first] += sim.opManager->biased.second->stats[0].count[entry.first];
            }
            repeat++;
        }
        for (auto& entry : times){
            std::cout << entry.first << "\t" << entry.second << std::endl;
        }
        // Calculate the mean, variance, and standard deviation of the w_new values
        long double sum = 0;
        long double sum_of_squares = 0;
        int count = 0;
        for (const auto &entry : times) {
            sum += entry.second;
            sum_of_squares += entry.second * entry.second;
            count++;
        }
        long double mean = sum / count;
        long double variance = (sum_of_squares / count) - (mean * mean);
        long double std_dev = std::sqrt(variance);
        StdOverMean = std_dev / mean;
        info_file << mean << "\t" << variance << "\t" << std_dev << "\t" << StdOverMean << std::endl;
        ofstream hist_file;
        hist_file.open("Weights/Hist.txt", std::ofstream::out | std::ofstream::trunc);
        hist_file << "Val\t";
        hist_file << "Count\t";
        hist_file << "Time\t";
        hist_file << "Weight\t";
        hist_file << "\n";
        const int Ti = 0;
        for (const auto &entry: weights) {
            int key = entry.first;
            hist_file << key << "\t";
            hist_file << counts[key] << "\t";
            hist_file << times[key] << "\t";
            hist_file << weights[key] << "\t";
            hist_file << "\n";
        }
        std::map<int,long double> w_new = get_new_weights(weights,times);
        ofstream w_new_file;
        w_new_file.open("Weights/"+og_w_file_name,std::ofstream::out | std::ofstream::trunc);
        Simulation sim = Simulation(&inputs, &constants);
        w_new_file << sim.opManager->biased.second->name << std::endl;
        for (const auto &entry : w_new) {
            w_new_file << entry.first << "\t" << entry.second << std::endl;
        }
        w_new_file.close();
        hist_file.close();
        info_file.close();
    }
    else{
        printf ("Please select sim type: test, anneal, melt, isothermal, config_generator.\n");
        exit (EXIT_FAILURE);
    }

    return 0;
}
