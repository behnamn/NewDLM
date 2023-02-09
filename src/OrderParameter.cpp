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
