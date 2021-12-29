#include <algorithm>
#include <chrono>
#include <functional>
#include<fstream>
#include<iostream>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stack>
#include <thread>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include "compute_graph_information.cpp"
#include "rootCounter-v2.cpp"


// Optimizations for speedup
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")


// read out fluxes
std::vector<std::vector<int>> read_fluxes(const int & file_number, const int & start, const int & end)
{
    
    // create file_name
    std::string file_name = "data_H1/fluxes_H1_" + std::to_string(file_number);
    
    // can be open the file?
    std::ifstream in(file_name.c_str());
    if(in.fail()){
        std::cout << "File " << file_name.c_str() << " not found \n";
    }
    
    // reserve a string
    std::string s;
    s.reserve(15);
        
    // skip as many lines as specified by variable start
    for(int j = 0; j < start; j++){
        std::getline(in, s);
    }
    
    // now read
    std::vector<std::vector<int>> fluxes;
    for(int i = start; i <= end; i++){
        
        // the get line
        std::getline(in,s);
        
        // cast the string s into a vector
        std::vector<int> flux;
        int pos0 = 0;
        int pos1 = s.find(",");
        int value;
        while (pos1 != -1) {
            std::stringstream ss;
            ss << s.substr(pos0, pos1 - pos0);
            ss >> value;
            flux.push_back(value);
            pos0 = pos1 + 1;
            pos1  = s.find(",", pos0);
        }
        std::stringstream ss;
        ss << s.substr(pos0, pos1 - pos0);
        ss >> value;
        flux.push_back(value);
        
        // save flux
        fluxes.push_back(flux);
        
    }
    
    // return the result
    return fluxes;
    
}



// determine root distribution for given outflux
void count_roots(const int & file_number, const int & start, const int & end)
{

    // (0) hard coded information for diagram 88
    int h0Max = 4;
    int root = 20;
    int genus = 6;
    std::vector<int> degrees = {42, 210, 84, 42, 42};
    std::vector<int> genera = {0,1,0,0,0};
    std::vector<std::vector<int>> edges = {{4,0},{0,3},{2,3},{2,4},{0,1},{1,4},{1,3},{1,2},{1,2}};
    
    // (1) compute additional information about this diagram
    std::vector<int> edge_numbers(degrees.size(),0);
    std::vector<std::vector<std::vector<int>>> graph_stratification;
    additional_graph_information(edges, edge_numbers, graph_stratification);
    
    // (2) read fluxes
    std::vector<std::vector<int>> fluxes = read_fluxes(file_number, start, end);
    
    // (3) for each flux, compute the distribution
    std::vector<std::vector<int>> non_trivial_fluxes;
    std::vector<std::vector<boost::multiprecision::int128_t>> non_trivial_distributions;
    for (int i = 0; i < fluxes.size(); i++){
        
        // (3.0) print status
        std::cout << "Status: " << i << '\r';
        
        // (3.1) compute the "reduced" degrees
        for (int j = 0; j < degrees.size(); j++){
            degrees[j] -= fluxes[i][j];
        }
        
        // (3.2) compute distribution on H1
        std::vector<boost::multiprecision::int128_t> dist(h0Max + 1, 0);
        int h0_min = (int)(std::accumulate(degrees.begin(),degrees.end(),0)/root) - genus + 1;
        if (h0_min < 0){
            h0_min = 0;
        }
        for (int j = h0_min; j <= h0Max; j++){
            dist[j] = NewRootDistributionCounter(degrees, genera, edges, root, graph_stratification, edge_numbers, j);
        }
        
        // (3.3) remember non-trivial results
        bool zeros = std::all_of(dist.begin(), dist.end(), [](boost::multiprecision::int128_t j) { return j==0; });
        if (!zeros){
            non_trivial_fluxes.push_back(fluxes[i]);
            non_trivial_distributions.push_back(dist);
        }
        
        // 3.5 flush line
        std::cout.flush();
        
    }
    
    // (4) print non-trivial fluxes
    std::ofstream ofile;
    ofile.open("results_H1/good_fluxes_H1_" + std::to_string(file_number), std::ios_base::app);
    for (int i = 0; i < non_trivial_fluxes.size(); i++){
        for (int j = 0; j < non_trivial_fluxes[i].size() -1; j ++){
            ofile << non_trivial_fluxes[i][j] << ",";
        }
        ofile << non_trivial_fluxes[i][non_trivial_fluxes[i].size()-1] << "\n";
    }
    ofile.close();

    // (5) print non-trivial distributions
    ofile.open("results_H1/distribution_H1_" + std::to_string(file_number), std::ios_base::app);
    for (int i = 0; i < non_trivial_distributions.size(); i++){
        for (int j = 0; j < non_trivial_distributions[i].size() -1; j ++){
            ofile << non_trivial_distributions[i][j] << ",";
        }
        ofile << non_trivial_distributions[i][non_trivial_distributions[i].size()-1] << "\n";
    }
    ofile.close();
    
}


// #################
// The main routine
// The main routine
// #################

int main(int argc, char* argv[]) {
    
    // check if we have the correct number of arguments
    if (argc != 2) {
        std::cout << "Error - number of arguments must be exactly 1 and not " << argc << "\n";
        std::cout << argv[ 0 ] << "\n";
        return 0;
    }
    
    // parse input
    std::string myString = argv[1];
    std::stringstream iss( myString );
    std::vector<int> input;
    int number;
    while ( iss >> number ){
        input.push_back( number );
    }
    
    // check input
    int file_number = input[0];
    int start = input[1];
    int end = input[2];
    if (start < 0 || end >= 728998 || file_number != 0){
        std::cout << "Invalid input.\n";
        return -1;
    }
    
    // compute distribution for given flux index
    std::cout << "Start: " << start << "\n";
    std::cout << "End: " << end << "\n\n";
    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    count_roots(file_number, start, end);
    std::chrono::steady_clock::time_point later = std::chrono::steady_clock::now();
    std::cout << "\nTime for run: " << std::chrono::duration_cast<std::chrono::seconds>(later - now).count() << "[s]\n";
    
    // return success
    return 0;
    
}
