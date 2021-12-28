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

// guard for thread-safe operations
boost::mutex myGuard;

//^#include "RootDistributionCounter.cpp"
//#include "ImprovedRootDistributionCounter.cpp"
#include "RootCounter-v2.cpp"


// Optimizations for speedup
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")


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
    
    // Hard coded information about diagram 88
    /*std::vector<int> degrees_H1 = {42, 210, 84, 42, 42};
    std::vector<int> flux_vector = {2,86,76,38,38};
    for (int i = 0; i < degrees_H1.size(); i++){
        degrees_H1[i] -= flux_vector[i];
    }
    std::vector<int> genera = {0,1,0,0,0};
    std::vector<std::vector<int>> edges = {{4,0},{0,3},{2,3},{2,4},{0,1},{1,4},{1,3},{1,2},{1,2}};
    int genus = 6;
    int root = 20;
    std::vector<std::vector<std::vector<int>>> graph_stratification = {{{1,3,4},{1,1,1},{4,2,2}},{{2,3,4},{2,1,1},{2,1,1}},{{3,4},{1,1},{0,0}}};
    std::vector<int> edge_numbers = {3,4,4,3,3};
    int h0_value = input[0];*/


    // Hard coded information about diagram 8
    /*std::vector<int> degrees_H1 = {12, 36, 12, 12};
    std::vector<int> genera = {0,1,0,0};
    std::vector<std::vector<int>> edges = {{3, 0}, {2, 0}, {2, 3}, {0, 1}, {1, 3}, {1, 2}};
    int genus = 4;
    int root = 12;
    std::vector<std::vector<std::vector<int>>> graph_stratification = {{{1,2,3},{1,1,1},{2,3,2}},{{2,3},{1,1},{1,1}},{{3},{1},{0}}};
    std::vector<int> edge_numbers = {3,3,3,3};
    int h0_value = input[0];*/

    // Hard coded information about diagram 8
    /*std::vector<int> degrees_H1 = {4, 4};
    std::vector<int> genera = {0,0};
    std::vector<std::vector<int>> edges = {{0, 1}, {0, 1}};
    int genus = 1;
    int root = 2;
    std::vector<std::vector<std::vector<int>>> graph_stratification = {{{1},{2},{0}}};
    std::vector<int> edge_numbers = {2,2};
    int h0_value = input[0];*/
    
    // Hard coded example
    std::vector<int> degrees_H1 = {16, 16, 16};
    std::vector<int> genera = {0, 0, 0};
    std::vector<std::vector<int>> edges = {{0,1}, {0,1}, {0,1}, {0,1}, {0,2}, {1,2}};
    int genus = 4;
    int root = 8;
    std::vector<std::vector<std::vector<int>>> graph_stratification = {{{1,2},{4,1},{1,1}},{{2},{1},{0}}};
    std::vector<int> edge_numbers = {5,5,2};
    int h0_value = input[0];
    
    // Count roots with new algorithm
    boost::multiprecision::int128_t total = NewRootDistributionCounter(degrees_H1, genera, edges, genus, root, graph_stratification, edge_numbers, h0_value );
    std::cout << "Total: " << total << "\n";
    
    //141264
    //142560
    // hi_old: [ 0, 0, 0, 2030, 70, 0, 0, 0, 0 ]
    // hi_new: [ 0, 0, 0, 14, 8, 0, 0, 0, 0 ]
    
    /*int N = 9;
    int n = 4;
    std::vector<int> minima = {1,1,1,1};
    std::vector<int> maxima = {7,7,7,7};
    std::vector<std::vector<int>> partitions;
    my_partitions(N, n, minima, maxima, partitions);
    std::cout << "\n";
    std::cout << partitions.size() << "\n\n";*/
    
    // return success
    return 0;
    
}
