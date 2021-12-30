// A program to compute the number of minimal limit roots on full blowups of nodal curves

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

// guard for thread-safe operations
boost::mutex myGuard;
#include "rootCounter-v2.cpp"

// Optimizations for speedup
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")

// Global variables
int thread_number = 8;
int display_details = false;

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
    
    // Original diagram 88
    int root = 20;
    int genus = 6;
    std::vector<int> degrees = {16, 80, 32, 16, 16};
    std::vector<int> flux_vector = {0,0,0,0,0};
    std::vector<int> genera = {0,1,0,0,0};
    std::vector<std::vector<int>> edges = {{4,0},{0,3},{2,3},{2,4},{0,1},{1,4},{1,3},{1,2},{1,2}};
    
    // Diagram 88
    /*int root = 20;
    int genus = 6;
    std::vector<int> degrees = {42, 210, 84, 42, 42};
    std::vector<int> flux_vector = {2,86,76,38,38};
    std::vector<int> genera = {0,1,0,0,0};
    std::vector<std::vector<int>> edges = {{4,0},{0,3},{2,3},{2,4},{0,1},{1,4},{1,3},{1,2},{1,2}};*/
    
    // Diagram 8
    /*int root = 12;
    int genus = 4;
    std::vector<int> degrees = {12,36,12,12};
    std::vector<int> flux_vector = {0,0,0,0};
    std::vector<int> genera = {0,1,0,0};
    std::vector<std::vector<int>> edges = {{3,0},{2,0},{2,3},{0,1},{1,3},{1,2}};*/
        
    // Hard coded information about example
    /*int root = 2;
    int genus = 1;
    std::vector<int> degrees = {4,4};
    std::vector<int> flux_vector = {0,0};
    std::vector<int> genera = {0,0};
    std::vector<std::vector<int>> edges = {{0,1},{0,1}};*/
    
    // Hard coded example
    /*int root = 8;
    int genus = 4;
    std::vector<int> degrees = {16,16,16};
    std::vector<int> flux_vector = {0,0,0};
    std::vector<int> genera = {0,0,0};
    std::vector<std::vector<int>> edges = {{0,1},{0,1},{0,1},{0,1},{0,2},{1,2}};*/
    
    // compute the "reduced" degrees
    for (int i = 0; i < degrees.size(); i++){
        degrees[i] -= flux_vector[i];
    }
    
    // compute additional information    
    std::vector<int> edge_numbers(degrees.size(),0);
    std::vector<std::vector<std::vector<int>>> graph_stratification;
    additional_graph_information(edges, edge_numbers, graph_stratification);
    
    // count roots
    boost::multiprecision::int128_t total = parallel_root_counter(degrees, genera, edges, root, graph_stratification, edge_numbers, input[0], thread_number, display_details);
    
    // return success
    return 0;
    
}
