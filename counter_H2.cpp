// A program to compute the number of minimal limit roots on full blowups of nodal curves

#include <algorithm>
#include <chrono>
#include <functional>
#include<fstream>
#include<iostream>
#include <mutex>
#include <sstream>
#include <stack>
#include <thread>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

// guard for thread-safe operations
boost::mutex myGuard;

#include "RootDistributionCounter.cpp"

// Optimizations for speedup
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")

// Number of threads
const int number_threads = 8;


// read the n-th line
std::string ReadNthLine(const std::string& filename, int N)
{
    std::ifstream in(filename.c_str());
    if(in.fail()){
        std::cout << "File " << filename.c_str() << " not found \n";
    }
    std::string s;
    s.reserve(15);
    
    // skip N lines
    for(int i = 0; i < N; ++i){
        std::getline(in, s);
    }
        
    // the get line
    std::getline(in,s);
    return s; 
}


// Determine root distribution for given outflux
void distributionOfFlux(int index)
{

    // (0) Hard coded information about diagram 88
    int numberVertices = 5;
    std::vector<int> vertices = {0,1,2,3,4};
    std::vector<int> degrees_H2 = {12, 60, 24, 12, 12};    
    std::vector<int> genera = {0,1,0,0,0};
    int numberEdges = 9;
    std::vector<std::vector<int>> edges = {{4,0},{0,3},{2,3},{2,4},{0,1},{1,4},{1,3},{1,2},{1,2}};
    std::vector<int> external_legs = {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4};
    std::vector<int> legs_per_component = {2, 10, 4, 2, 2};
    int genus = 6;
    int root = 20;
    int number_threads = 8;
    int h0Min = 0;
    int h0Max = 10;
    bool display_details = false;
    
    // (1) Read flux
    int file_number = (int) index / 1000000;
    int line_number = index % 1000000;
    std::string file_name = "data_H2/fluxes_H2_" + std::to_string(file_number);
    std::string flux = ReadNthLine( file_name, line_number );
    
    // (2) Cast the flux string into a vector
    std::vector<int> flux_vector;
    int start = 0;
    int end = flux.find(",");
    int value;
    while (end != -1) {
        std::stringstream ss;
        ss << flux.substr(start, end - start);
        ss >> value;
        flux_vector.push_back(value);
        start = end + 1;
        end = flux.find(",", start);
    }
    std::stringstream ss;
    ss << flux.substr(start, end - start);
    ss >> value;
    flux_vector.push_back(value);

    // (3) Compute distribution on H2
    std::vector<boost::multiprecision::int128_t> dist_H2(h0Max + 1, 0);
    int total_flux = std::accumulate( flux_vector.begin(), flux_vector.end(), 0 );
    if ( total_flux % root == 0 ){
        
        // (3.1) Generate weights for the outflux
        std::vector<int> external_weights(external_legs.size());
        int iterator = 0;
        for ( int j = 0; j < legs_per_component.size(); j++ ){
            int flux_for_component = flux_vector[j];
            int leg_number = legs_per_component[j];
            int remaining_legs = leg_number;
            for (int k = 0; k < leg_number; k++){
                external_weights[iterator] = flux_for_component / remaining_legs;
                flux_for_component = flux_for_component - (flux_for_component / remaining_legs);
                remaining_legs = remaining_legs - 1;
                iterator = iterator + 1;
            }
        }
        
        // (3.2) Compute distribution on H1
        WeightedDiagramWithExternalLegs dia_H2 = WeightedDiagramWithExternalLegs(vertices, degrees_H2, genera, edges, external_legs, external_weights, genus, root);
        int h0MinUsed_H2 = dia_H2.get_h0_min();
        std::vector<boost::multiprecision::int128_t> n_H2(h0Max - h0MinUsed_H2 + 1, 0);
        if (h0Max >= h0MinUsed_H2){
            countRootDistribution(dia_H2, number_threads, h0MinUsed_H2, h0Max, n_H2, display_details);
            for (int j = 0; j < n_H2.size(); j++){
                dist_H2[h0MinUsed_H2+j] = n_H2[j];
            }
        }
    
    }

    // (4) Save non-trivial results on H1
    bool zeros = std::all_of(dist_H2.begin(), dist_H2.end(), [](boost::multiprecision::int128_t j) { return j==0; });
    if (!zeros){
        
        // save this non-trivial flux
        std::ofstream ofile;
        ofile.open("results_H2/good_fluxes_H2_" + std::to_string(file_number), std::ios_base::app);
        for ( int i = 0; i < flux_vector.size() -1; i ++ ){
            ofile << flux_vector[i] << ",";
        }
        ofile << flux_vector[flux_vector.size()-1] << "\n";
        ofile.close();
        
        // save this non-trivial distribution
        ofile.open("results_H2/distribution_H2_" + std::to_string(file_number), std::ios_base::app);
        for ( int i = 0; i < dist_H2.size() -1; i ++ ){
            ofile << dist_H2[i] << ",";
        }
        ofile << dist_H2[dist_H2.size()-1] << "\n";
        ofile.close();
        
    }

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
    int start = input[0];
    int end = input[1];    
    if ( start < 0 || end >= 33463895 ){
        std::cout << "Flux index too small or too large.\n";
        return -1;
    }
    
    // compute distribution for given flux index
    std::cout << "Start: " << start << "\n";
    std::cout << "End: " << end << "\n\n";
    for (int i = start; i <= end; i++){
        distributionOfFlux(i);
        std::cout << i << "\n";
    }
    
    // return success
    return 0;
    
}
