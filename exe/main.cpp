#include "../src/Constants.h"
#include "../src/Design.h"
#include "../src/MyGraph.h"
#include "../src/TempRamp.h"
#include "../src/Simulation.h"
#include "../src/Exact.h"
#include "../src/WeighGen.h"


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

int main(int argc, char * argv[]) {
    Inputs inputs = Inputs(argc,argv);
    Constants constants = Constants(&inputs);

    if (inputs.weight_generator || inputs.isothermal) {
        std::cout << "Input Temperature: " << inputs.temp << endl;
        if (inputs.temp < 0) {
            string ogRateModel = inputs.rate_model;
            string ogSimType = inputs.sim_type;
            if (ogRateModel == "local") inputs.rate_model = "global";
            inputs.sim_type = "exact";
            double deltaT = -1. * (inputs.temp + 100.);
            Simulation sim = Simulation(&inputs, &constants);
            double Tm = calculate_Tm(&sim);
            inputs.temp = centigrade(Tm) + deltaT;
            inputs.rate_model = ogRateModel;
            inputs.sim_type = ogSimType;
        }
        std::cout << "New Temperature: " << inputs.temp << endl;
    }

    if (inputs.test){
        Design design = Design(&inputs);
        for (const auto& dom : design.domains){
            std::cout << dom.id << "\t" << dom.dH_seq << "\t" << dom.dH_seqave << std::endl;
        }
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
        make_weight_directory();
        ofstream infoFile;
        open_info_file(infoFile);
        string opName;
        std::map<int,long double> times, weights, Times, Weights;
        std::map<int,int> counts, Counts;
        vector<long double> mean_var_std;
        long double mean, variance, std_dev, StdOverMean;
        int repeat;
        bool finished;
        std::string wFileNameOG = inputs.w_file_name;
        std::vector<char> limits = {'l', 'm', 'h'};

        std::cout << "------------------------------------------------ Beginning LMH Steps " << std::endl;
        repeat = 1;
        finished = false;
        while (!finished){
            std::cout << "------------------------ LMH Step " << repeat << " of " << inputs.num_repeats << std::endl;
            for (char limit : limits) {
                std::cout << "--------- Step " << limit << std::endl;
                inputs.seed = std::chrono::system_clock::now().time_since_epoch().count();
                infoFile << limit << "\t" << inputs.seed << "\t" << inputs.w_file_name << "\t"; infoFile << std::flush;
                Simulation sim = Simulation(&inputs, &constants);
                sim.prepare_config(limit);
                sim.run();
                opName = sim.opManager->biased.second->name;
                weights = sim.opManager->biased.second->weight;
                times = sim.opManager->biased.second->stats[0].time;
                counts = sim.opManager->biased.second->stats[0].count;
                mean_var_std = get_mean_var_std(times);
                mean = mean_var_std[0]; variance = mean_var_std[1]; std_dev = mean_var_std[2];
                StdOverMean = std_dev / mean;
                infoFile << mean << "\t" << variance << "\t" << std_dev << "\t" << StdOverMean << std::endl; infoFile << std::flush;
                if (limit == 'l'){
                    for (auto& entry : weights) {
                        Weights[entry.first] = entry.second;
                        Times[entry.first] += 0;
                        Counts[entry.first] += 0;
                    }
                }
                for (auto& entry : Weights){
                    if (!(times.find(entry.first) == times.end()))
                        Times[entry.first] += times[entry.first];
                    if (!(counts.find(entry.first) == counts.end()))
                        Counts[entry.first] += counts[entry.first];
                }
            }
            std::string histFileNameNew = "Weights/LMH" + std::to_string(repeat) + "_Hist.txt";
            std::string wFileNameNew = "Weights/LMH" + std::to_string(repeat) + "_" + wFileNameOG;
            write_hist(histFileNameNew, Weights, Counts, Times);
            write_new_weights(wFileNameNew, opName, Weights, Times);
            mean_var_std = get_mean_var_std(Times);
            mean = mean_var_std[0]; variance = mean_var_std[1]; std_dev = mean_var_std[2];
            StdOverMean = std_dev / mean;
            infoFile << "LMH" << repeat << "\t" << "x3" << "\t" << inputs.w_file_name << "\t";
            infoFile << mean << "\t" << variance << "\t" << std_dev << "\t" << StdOverMean << std::endl;
            infoFile << std::flush;
            inputs.w_file_name = wFileNameNew;
            repeat++;
            if (!has_zero_value(Counts) || repeat == inputs.num_repeats){
                finished = true;
            }
        }

        std::cout << "------------------------------------------------ Beginning Random Steps " << std::endl;
        repeat = 1;
        StdOverMean = 9999999;
        finished = false;
        while (!finished){
            std::cout << "------------------------ Random Step " << repeat << " of " << inputs.num_repeats << std::endl;
            for (int limit : {1,2,3}) {
                std::cout << "--------- Step " << limit << std::endl;
                inputs.seed = std::chrono::system_clock::now().time_since_epoch().count();
                infoFile << limit << "\t" << inputs.seed << "\t" << inputs.w_file_name << "\t"; infoFile << std::flush;
                Simulation sim = Simulation(&inputs, &constants);
                sim.prepare_config('r');
                sim.run();
                opName = sim.opManager->biased.second->name;
                weights = sim.opManager->biased.second->weight;
                times = sim.opManager->biased.second->stats[0].time;
                counts = sim.opManager->biased.second->stats[0].count;
                mean_var_std = get_mean_var_std(times);
                mean = mean_var_std[0]; variance = mean_var_std[1]; std_dev = mean_var_std[2];
                StdOverMean = std_dev / mean;
                infoFile << mean << "\t" << variance << "\t" << std_dev << "\t" << StdOverMean << std::endl; infoFile << std::flush;
                if (limit == 1){
                    for (auto& entry : weights) {
                        Weights[entry.first] = entry.second;
                        Times[entry.first] += 0;
                        Counts[entry.first] += 0;
                    }
                }
                for (auto& entry : Weights){
                    if (!(times.find(entry.first) == times.end()))
                        Times[entry.first] += times[entry.first];
                    if (!(counts.find(entry.first) == counts.end()))
                        Counts[entry.first] += counts[entry.first];
                }
            }
            std::string histFileNameNew = "Weights/r" + std::to_string(repeat) + "_Hist.txt";
            std::string wFileNameNew = "Weights/r" + std::to_string(repeat) + "_" + wFileNameOG;
            write_hist(histFileNameNew, Weights, Counts, Times);
            write_new_weights(wFileNameNew, opName, Weights, Times);
            mean_var_std = get_mean_var_std(Times);
            mean = mean_var_std[0]; variance = mean_var_std[1]; std_dev = mean_var_std[2];
            StdOverMean = std_dev / mean;
            infoFile << "r" << repeat << "\t" << "x3" << "\t" << inputs.w_file_name << "\t";
            infoFile << mean << "\t" << variance << "\t" << std_dev << "\t" << StdOverMean << std::endl;
            infoFile << std::flush;
            inputs.w_file_name = wFileNameNew;
            repeat++;
            if ( (StdOverMean < 1 && !has_zero_value(Counts)) || repeat == inputs.num_repeats ){
                finished = true;
            }
        }

        std::cout << "------------------------------------------------ Beginning Final Steps a" << std::endl;
        repeat = 1;
        while (repeat < inputs.num_repeats){
            std::cout << "------------------------ Final Step " << repeat << " of " << inputs.num_repeats << std::endl;
            inputs.seed = std::chrono::system_clock::now().time_since_epoch().count();
            infoFile << repeat << "\t" << inputs.seed << "\t" << inputs.w_file_name << "\t"; infoFile << std::flush;
            Simulation sim = Simulation(&inputs, &constants);
            sim.prepare_config();
            sim.run();
            opName = sim.opManager->biased.second->name;
            Weights = weights = sim.opManager->biased.second->weight;
            times = sim.opManager->biased.second->stats[0].time;
            counts = sim.opManager->biased.second->stats[0].count;
            mean_var_std = get_mean_var_std(times);
            mean = mean_var_std[0]; variance = mean_var_std[1]; std_dev = mean_var_std[2];
            StdOverMean = std_dev / mean;
            infoFile << mean << "\t" << variance << "\t" << std_dev << "\t" << StdOverMean << std::endl; infoFile << std::flush;
            for (auto& entry : Weights){
                if (!(times.find(entry.first) == times.end()))
                    Times[entry.first] += times[entry.first];
                if (!(counts.find(entry.first) == counts.end()))
                    Counts[entry.first] += counts[entry.first];
            }
            repeat++;
        }
        for (auto& entry : Times){
            std::cout << entry.first << "\t" << entry.second << std::endl;
        }
        mean_var_std = get_mean_var_std(Times);
        mean = mean_var_std[0];
        variance = mean_var_std[1];
        std_dev = mean_var_std[2];
        StdOverMean = std_dev / mean;
        infoFile << "a_x" << inputs.num_repeats << "\t" << "NA" << "\t" << inputs.w_file_name << "\t";
        infoFile << mean << "\t" << variance << "\t" << std_dev << "\t" << StdOverMean << std::endl;
        infoFile << std::flush;
        write_hist("Weights/a_Hist.txt", Weights, Counts, Times);
        write_new_weights("Weights/a_"+wFileNameOG, opName, Weights, Times);
        inputs.w_file_name = "Weights/a_"+wFileNameOG;

        std::cout << "------------------------------------------------ Beginning Final Steps b" << std::endl;
        repeat = 1;
        while (repeat < inputs.num_repeats){
            std::cout << "------------------------ Final Step " << repeat << " of " << inputs.num_repeats << std::endl;
            inputs.seed = std::chrono::system_clock::now().time_since_epoch().count();
            infoFile << repeat << "\t" << inputs.seed << "\t" << inputs.w_file_name << "\t"; infoFile << std::flush;
            Simulation sim = Simulation(&inputs, &constants);
            sim.prepare_config();
            sim.run();
            opName = sim.opManager->biased.second->name;
            Weights = weights = sim.opManager->biased.second->weight;
            times = sim.opManager->biased.second->stats[0].time;
            counts = sim.opManager->biased.second->stats[0].count;
            mean_var_std = get_mean_var_std(times);
            mean = mean_var_std[0]; variance = mean_var_std[1]; std_dev = mean_var_std[2];
            StdOverMean = std_dev / mean;
            infoFile << mean << "\t" << variance << "\t" << std_dev << "\t" << StdOverMean << std::endl; infoFile << std::flush;
            for (auto& entry : Weights){
                if (!(times.find(entry.first) == times.end()))
                    Times[entry.first] += times[entry.first];
                if (!(counts.find(entry.first) == counts.end()))
                    Counts[entry.first] += counts[entry.first];
            }
            repeat++;
        }
        for (auto& entry : Times){
            std::cout << entry.first << "\t" << entry.second << std::endl;
        }
        mean_var_std = get_mean_var_std(Times);
        mean = mean_var_std[0];
        variance = mean_var_std[1];
        std_dev = mean_var_std[2];
        StdOverMean = std_dev / mean;
        infoFile << "b_x" << inputs.num_repeats << "\t" << "NA" << "\t" << inputs.w_file_name << "\t";
        infoFile << mean << "\t" << variance << "\t" << std_dev << "\t" << StdOverMean << std::endl;
        infoFile << std::flush;
        write_hist("Weights/b_Hist.txt", Weights, Counts, Times);
        write_new_weights("Weights/b_"+wFileNameOG, opName, Weights, Times);

        infoFile.close();
    }
    else{
        printf ("Please select sim type: test, anneal, melt, isothermal, config_generator.\n");
        exit (EXIT_FAILURE);
    }

    return 0;
}

/*
         std::cout << "------------------------------------------------ Beginning Generation Steps " << std::endl;
        finished = false;
        StdOverMean = 9999999;
        repeat = 0;
        while (!finished){
            std::cout << "------------------------ Weight Generation Step " << repeat << " of " << inputs.num_repeats << std::endl;
            inputs.seed = std::chrono::system_clock::now().time_since_epoch().count();
            infoFile << repeat << "\t" << inputs.seed << "\t" << inputs.w_file_name << "\t";
            Simulation sim = Simulation(&inputs, &constants);
            sim.prepare_config();
            sim.run();
            opName = sim.opManager->biased.second->name;
            weights = sim.opManager->biased.second->weight;
            times = sim.opManager->biased.second->stats[0].time;
            counts = sim.opManager->biased.second->stats[0].count;
            std::string wFileNameNew = "Weights/" + std::to_string(repeat) + "_" + wFileNameOG;
            std::string histFileNameNew = "Weights/" + std::to_string(repeat) + "_Hist.txt";
            write_new_weights(wFileNameNew, opName, weights, times);
            write_hist(histFileNameNew, weights, counts, times);
            mean_var_std = get_mean_var_std(times);
            mean = mean_var_std[0];
            variance = mean_var_std[1];
            std_dev = mean_var_std[2];
            StdOverMean = std_dev / mean;
            infoFile << mean << "\t" << variance << "\t" << std_dev << "\t" << StdOverMean << std::endl;
            inputs.w_file_name = wFileNameNew;
            repeat++;
            if (repeat == inputs.num_repeats || ( StdOverMean < 1 && !has_zero_value(counts) ) ){
                finished = true;
            }
        }
 */