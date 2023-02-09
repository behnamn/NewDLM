//
//  OPManager.h
//  DLM
//
//  Created by Behnam Najafi on 03/06/2017.
//  Copyright © 2017 Behnam Najafi. All rights reserved.
//

#ifndef OPManager_h
#define OPManager_h


#include "StatManager.h"
//#include "OrderParameter.h"

class OrderParameter2D{
public:
    OrderParameter2D(){}
    OrderParameter2D(OrderParameter*, OrderParameter*);
    virtual ~OrderParameter2D(){}
    
    OrderParameter* OP1;
    OrderParameter* OP2;
    
    string name;
    
    pair<int,int> prev_state;
    pair<int,int> state;
    pair<int,int> fut_state;
    
    std::set<pair<int,int> > explored_vals;
    std::map<pair<int,int>, int> count;
    std::map<pair<int,int>, double> time;
    std::map<pair<int,int>, double> weight;
    
    void set_new_value();
    void update(const double);
    
    void print();
    
    void write();
};


class OPManager{
public:
    OPManager(){}
    OPManager(StatManager*);
    virtual ~OPManager(){}
    
    Inputs* inputs;
    Design* design;
    MyGraph* G;
    TempRamp* ramp;
    TransitionManager* trManager;
    StatManager* statManager;
    FileIO* ofiles;

    vector<OrderParameter*> OPs;
    
    void initialise();
    void read_weights();
    void initialise_2D();
    
    void set_values();
    void update_times();
    
    std::set<int> all_vals;
    
    pair<bool,OrderParameter*> biased;
    
    vector<OrderParameter2D> OPs_2D;
    vector<ofstream> files_2D;
    vector<ofstream> burn_2D;
    
    void write();
    void write_last();
    void write_burn();
    
    void set_future_value(const TR& tr);
    void fill_rates_w();
    
    void generate_config();
    vector<int> done_configs;
    
    void print();
};


#endif
