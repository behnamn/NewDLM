//
//  OrderParameter.hpp
//  DLM
//
//  Created by Behnam Najafi on 03/06/2017.
//  Copyright © 2017 Behnam Najafi. All rights reserved.
//

#ifndef OrderParameter_h
#define OrderParameter_h

#include "Common.h"

class OrderParameter{
//private:
public:
    OrderParameter(){id = -1;}
    OrderParameter(int, OP_t, string);
    static int counter;
    int id;
    
    int pool_id;
    OP_t type;
    
    string name;
    
    int prev_prev_state;
    int prev_state;
    int state;
    
    bool biased;
    int fut_state;
    std::map<int,double> weight;
    
    
    std::set<int> group_ids;
    
    std::set<int> explored_vals;
    
    /*
     Temperature Dependent members:
     index of vector is current_index
     in TempRamp
     */
    vector<Stat<int>> stats;
};



 
#endif /* OrderParameter_h */
