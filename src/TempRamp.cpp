/*
 *
 *	TempRamp.cpp
 * 	Author: Behnam
 *
 */

#include "TempRamp.h"


T_Window::T_Window(){
    idx = -1;
}
T_Window::T_Window(double T_) : T(T_){
    idx = -1;
}
T_Window::T_Window(double T_, pair<double,double> t_limits_) :
T(T_), t_limits(t_limits_){
    idx = count++;
    str_T = to_string_with_dec(centigrade(T_));
    for(int i=0; i<Tr_MAX; i++){
        num_tr.push_back(0);
    }
}
int T_Window::count = 0;

TempJump::TempJump(double dt_, int max_idx_) : dt(dt_), max_idx(max_idx_){
    was_changed = false;
    backward_remainder = 0;
    forward_remainder = 0;
}
void TempJump::update(const int& idx0,
                      const int& idxf,
                      const double& t0,
                      const double& tf){
    this->idx_0 = idx0;
    this->idx_f = idxf;
    if (idx0 == idxf){
        was_changed = false;
    }
    else{
        was_changed = true;
        backward_remainder = (dt * (idx0+1) ) - t0;
        forward_remainder = tf - (idxf * dt);
    }
}



TempRamp::TempRamp(){}
TempRamp::TempRamp(Inputs * inputs_):
inputs(inputs_){
    this->anneal = inputs->anneal;
    this->melt = inputs->melt;
    this->isothermal = inputs->isothermal;
    
    if (inputs->test || inputs->exact || inputs->config_generator){
        this->isothermal = true;
    }
    
    if(anneal || melt) {
        T_high = kelvin(inputs->max_temp);
        T_low = kelvin(inputs->min_temp);
        cool_rate = Cpm2spC(inputs->cool_rate);
        dt = inputs->const_T_interval;
        t_max = (abs(T_high-T_low) * cool_rate) + dt;
        dT = dt / cool_rate;
    }
    else if(isothermal) {
        iso_T = kelvin(inputs->temp);
        t_max = inputs->max_time;
    }
    else{
        printf ("Please select anneal, melt or isothermal to fill ramp. \n");
        exit (EXIT_FAILURE);
    }
    this->initialise();
}
void TempRamp::initialise(){
    double t_low, t_high;
    int dummy_idx = 0;
    if (this->anneal){
        for (double T = T_high; T >= T_low-(dT/2); T-= dT){
            t_low = dummy_idx * dt;
            t_high = t_low + dt;
            Temps.push_back(T_Window(T,make_pair(t_low, t_high)));
            dummy_idx++;
        }
    }
    else if (this->melt){
        for (double T = T_low; T <= T_high+(dT/2); T+= dT){
            t_low = dummy_idx * dt;
            t_high = t_low + dt;
            Temps.push_back(T_Window(T,make_pair(t_low, t_high)));
            dummy_idx++;
        }
    }
    else if (this->isothermal){
        t_low = 0;
        t_high = t_max;
        Temps.push_back(T_Window(iso_T,make_pair(t_low, t_high)));
    }
    else{
        printf ("Please select anneal, melt or isothermal to fill ramp. \n");
        exit (EXIT_FAILURE);
    }
    previous_t = -1;
    current_t = 0.;
    prev_idx = -1;
    idx = 0;
    T_jump = TempJump(dt,Temps.size()-1);
    if (this->isothermal) {
        T_jump.idx_0 = idx;
        T_jump.idx_f = idx;
    }
}
/*
void TempRamp::move_time(const double& tau){
    prev_idx = idx;
    previous_t = current_t;
    if(anneal || melt) {
        if (tau>dt){
            current_t = current_t + dt;
            idx++;
            Temps[prev_idx].num_rejected++;
        }
        else{
         
            current_t = current_t + tau;
            idx = current_t / dt;
        }
        T_jump.update(prev_idx, idx, previous_t, current_t);
    }
    else if(isothermal) {
        current_t = current_t + tau;
    }
    else{
        printf ("Please select anneal, melt or isothermal to fill ramp. \n");
        exit (EXIT_FAILURE);
    }
}
 */
void TempRamp::move_time(const double& tau){
    prev_idx = idx;
    previous_t = current_t;
    if(anneal || melt) {
        //current_t = current_t + tau;
        //idx = current_t / dt;
        if (tau>dt){
            current_t = current_t + dt;
            idx++;
            Temps[prev_idx].num_rejected++;
        }
        else{
            current_t = current_t + tau;
            idx = current_t / dt;
        }
        T_jump.update(prev_idx, idx, previous_t, current_t);
    }
    else if(isothermal) {
        current_t = current_t + tau;
    }
    else{
        printf ("Please select anneal, melt or isothermal to fill ramp. \n");
        exit (EXIT_FAILURE);
    }
}


void TempRamp::print(){
    cout << "--------------\tTempRamp::print\t--------------\n";
    for (const auto& Temp : Temps){
        cout << Temp.idx << "\t";
        cout << Temp.str_T << "C\t";
        cout << Temp.t_limits.first << "s->";
        cout << Temp.t_limits.second << "s\n";
    }
}
double TempRamp::get_T(){
    return Temps[idx].T;
}





