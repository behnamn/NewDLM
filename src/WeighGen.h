//
// Created by Behnam Najafi on 31/03/2023.
//

#ifndef NEWDLM_WEIGHGEN_H
#define NEWDLM_WEIGHGEN_H

#include "../src/Simulation.h"


std::map<int,long double> get_new_weights(std::map<int,long double>& w_old, std::map<int,long double>& times){
    int key;
    const int Ti = 0;
    std::map<int,long double> probs;
    std::map<int,long double> w_new;
    for (const auto& entry : w_old){
        key = entry.first;
        probs[key] = times[key] / (w_old[key] * w_old[key]);
    }
    for (const auto &entry : probs) {
        w_new[entry.first] = std::sqrt(1 / entry.second);
    }
    // Find the minimum value of w_new
    long double min_w_new = std::numeric_limits<long double>::max();
    for (const auto &entry : w_new) {
        min_w_new = std::min(min_w_new, entry.second);
    }

    // Find the maximum value of w_new_dict excluding inf values
    long double max_w_new = std::numeric_limits<long double>::lowest();
    for (const auto &entry : w_new) {
        if (!std::isinf(entry.second)) {
            max_w_new = std::max(max_w_new, entry.second);
        }
    }

    // Find the first and last non-infinite values
    int first_non_inf_key = -1;
    int last_non_inf_key = -1;
    for (const auto& entry : w_new) {
        if (!std::isinf(entry.second)) {
            if (first_non_inf_key == -1) {
                first_non_inf_key = entry.first;
            }
            last_non_inf_key = entry.first;
        }
    }

    // Iteratively increase the weight based on the distance from the first and last non-infinite values
    long double increment_first, increment_last;
    //int min_key = w_new.begin()->first;
    //int max_key = w_new.rbegin()->first;;
    //if (first_non_inf_key != max_key)
    increment_first = w_new[first_non_inf_key] / w_new[first_non_inf_key+1];
    //if (last_non_inf_key != min_key)
    increment_last = w_new[last_non_inf_key] / w_new[last_non_inf_key-1];


    for (auto& entry : w_new) {
        if (std::isinf(entry.second)) {
            int distance_to_first = std::abs(entry.first - first_non_inf_key);
            int distance_to_last = std::abs(entry.first - last_non_inf_key);
            int min_distance = std::min(distance_to_first, distance_to_last);


            if (entry.first < first_non_inf_key){
                entry.second = pow(increment_first, distance_to_first) * w_new[first_non_inf_key];
            }
            else if (entry.first > last_non_inf_key){
                entry.second = pow(increment_last, distance_to_last) * w_new[last_non_inf_key];
            }
            else{
                std::cout << "Error: key has to be less than first non-inf key or greater than last non-inf key." << std::endl;
            }
        }
    }

    /*
    // Change inf values to double the max
    for (auto &entry : w_new) {
        if (std::isinf(entry.second)) {
            entry.second = max_w_new*2.;
        }
    }
     */
    // Normalize w_new_dict by dividing each value by the minimum value
    for (auto &entry : w_new) {
        entry.second /= min_w_new;
    }
    return w_new;
}

void make_weight_directory(){
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
}

void open_info_file(ofstream& info_file){
    std::string info_file_name = "Weights/Info.txt";
    info_file.open(info_file_name,std::ofstream::out | std::ofstream::trunc);
    info_file << "Repeat" << "\t" << "Seed" << "\t" << "wfile" << "\t";
    info_file << "Mean" << "\t" << "Variance" << "\t" << "Std" << "\t" << "Std/Mean";
    info_file << std::endl;
}

vector<long double> get_mean_var_std(const map<int,long double>& mymap){
    // Calculate the mean, variance, and standard deviation of the w_new values

    long double sum = 0;
    long double sum_of_squares = 0;
    int count = 0;
    for (const auto &entry : mymap) {
        sum += entry.second;
        sum_of_squares += entry.second * entry.second;
        count++;
    }
    long double mean = sum / count;
    long double variance = (sum_of_squares / count) - (mean * mean);
    long double std_dev = std::sqrt(variance);
    vector<long double> result = {mean, variance, std_dev};
    return result;
}

bool has_zero_value(const std::map<int, int>& my_map) {
    for (const auto& entry : my_map) {
        if (entry.second == 0) {
            return true;
        }
    }
    return false;
}

void write_hist(const string& path, map<int,long double>& weight, map<int,int>& count, map<int,long double>& time){
    int key;
    ofstream hist_file;
    hist_file.open(path, std::ofstream::out | std::ofstream::trunc);
    hist_file << "Val\t";
    hist_file << "Count\t";
    hist_file << "Time\t";
    hist_file << "Weight\t";
    hist_file << "\n";
    const int Ti = 0;
    for (const auto &entry: weight) {
        key = entry.first;
        hist_file << key << "\t";
        hist_file << count[key] << "\t";
        hist_file << time[key] << "\t";
        hist_file << weight[key] << "\t";
        hist_file << "\n";
    }
    hist_file.close();
}

void write_new_weights(const string& path, const string& opName, map<int,long double>& weight, map<int,long double>& time){
    std::map<int,long double> w_new = get_new_weights(weight,time);
    ofstream w_new_file;
    w_new_file.open(path,std::ofstream::out | std::ofstream::trunc);
    w_new_file << opName << std::endl;
    for (const auto &entry : w_new) {
        w_new_file << entry.first << "\t" << entry.second << std::endl;
    }
    w_new_file.close();
}


#endif //NEWDLM_WEIGHGEN_H
