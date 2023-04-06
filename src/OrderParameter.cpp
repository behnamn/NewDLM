//
//  OrderParameter.cpp
//  DLM
//
//  Created by Behnam Najafi on 03/06/2017.
//  Copyright © 2017 Behnam Najafi. All rights reserved.
//

#include "OrderParameter.h"

OrderParameter::OrderParameter(int pool_id_, OP_t type_, string name_):
pool_id(pool_id_), type(type_), name(name_) {
    this->biased = false;
    id = counter++;
}
int OrderParameter::counter = 0;

ostream &operator<<(ostream &os, const OrderParameter &orderParameter) {
    os << "name: " << orderParameter.name << " prev_prev_state: " << orderParameter.prev_prev_state << " prev_state: "
       << orderParameter.prev_state << " state: " << orderParameter.state << " biased: " << orderParameter.biased
       << " fut_state: " << orderParameter.fut_state;//  << " stats: " << orderParameter.stats;
    return os;
}


OrderParameter2D::OrderParameter2D(OrderParameter* OP1_, OrderParameter* OP2_) : OP1(OP1_) , OP2(OP2_){
    this->name = OP1->name;
    this->name+= "-";
    this->name+= OP2->name;
    state.first = OP1->state;
    state.second = OP2->state;
}
void OrderParameter2D::set_new_value(){
    prev_state = state;
    state.first = OP1->state;
    state.second = OP2->state;
}
void OrderParameter2D::update(const double tau){
    count[state]++;
    time[state]+= tau;
}