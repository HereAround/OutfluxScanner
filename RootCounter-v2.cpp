#include "combinatorics.cpp"

// Count root bundle distribution on nodal curve
boost::multiprecision::int128_t NewRootDistributionCounter(
                                const std::vector<int> degrees,
                                const std::vector<int> genera,
                                const std::vector<std::vector<int>> edges,
                                const int genus,
                                const int root,
                                const std::vector<std::vector<std::vector<int>>> graph_stratification,
                                const std::vector<int> edge_numbers,
                                const int h0_value )
{
    
    // save total number of roots found
    boost::multiprecision::int128_t total = 0;
    
    // (1) Partition h0
    // (1) Partition h0
    std::vector<std::vector<int>> partitions;
    my_partitions(h0_value, degrees.size(), std::vector<int>(degrees.size(),0), std::vector<int>(degrees.size(),h0_value), partitions);
    
    
    // (2) Find corresponding fluxes
    // (2) Find corresponding fluxes
    struct flux_data{
        std::vector<int> flux;
        std::vector<int> partition;
    };
    std::vector<flux_data> outfluxes;
    for (int i = 0; i < partitions.size(); i++){
        
        // create stack
        std::stack<flux_data> snapshotStack;
        
        // add first snapshot
        flux_data currentSnapshot;
        currentSnapshot.flux = {};
        currentSnapshot.partition = partitions[i];
        snapshotStack.push(currentSnapshot);
        
        // Run...
        while(!snapshotStack.empty())
        {
        
            // pick the top snapshot and delete it from the stack
            currentSnapshot= snapshotStack.top();
            snapshotStack.pop();
            
            // any fluxes to be set?
            if (currentSnapshot.flux.size() < degrees.size()){
                
                // determine vertex for which we determine the outflux
                int j = currentSnapshot.flux.size();
                
                // non-trivial h0:
                if (currentSnapshot.partition[j] > 0){
                    int f = degrees[j] - root * currentSnapshot.partition[j];
                    if (genera[j] == 0){
                        f += root;
                    }
                    if ((edge_numbers[j] <= f) && (f <= edge_numbers[j] * (root-1)) && ((degrees[j] - f) % root == 0)){
                        std::vector<int> new_flux = currentSnapshot.flux;
                        new_flux.push_back(f);
                        flux_data newSnapshot;
                        newSnapshot.flux = new_flux;
                        newSnapshot.partition = currentSnapshot.partition;
                        snapshotStack.push(newSnapshot);
                    }
                }
                
                // trivial h0:
                if (currentSnapshot.partition[j] == 0){
                    int min_flux = degrees[j];
                    if (genera[j] == 0){
                        min_flux++;
                    }
                    if (min_flux < edge_numbers[j]){
                        min_flux = edge_numbers[j];
                    }
                    for (int k = min_flux; k <= edge_numbers[j]*(root-1); k++){
                        if ((degrees[j] - k) % root == 0){
                            std::vector<int> new_flux = currentSnapshot.flux;
                            new_flux.push_back(k);
                            flux_data newSnapshot;
                            newSnapshot.flux = new_flux;
                            newSnapshot.partition = currentSnapshot.partition;
                            snapshotStack.push(newSnapshot);
                        }
                    }
                }
            
            }
            // no more fluxes to be set --> add to list of fluxes if the sum of fluxes equals the number of edges * root (necessary and sufficient for non-zero number of weight assignments)
            else if (std::accumulate(currentSnapshot.flux.begin(),currentSnapshot.flux.end(),0) == root * edges.size()){
                outfluxes.push_back(currentSnapshot);
            }
            
        }
    
    }
    
    
    // (3) Count weight assignments
    // (3) Count weight assignments
    struct comb_data{
        std::vector<int> flux;
        int k;
        boost::multiprecision::int128_t mult;
    };
    for (int i = 0; i < outfluxes.size(); i++){
        
        // create stack
        std::stack<comb_data> snapshotStack;
        
        // add first snapshot
        comb_data currentSnapshot;
        currentSnapshot.flux = outfluxes[i].flux;
        currentSnapshot.k = 0;
        currentSnapshot.mult = (boost::multiprecision::int128_t) 1;
        snapshotStack.push(currentSnapshot);
        
        // Run...
        while(!snapshotStack.empty())
        {
        
            // pick the top snapshot and delete it from the stack
            currentSnapshot= snapshotStack.top();
            snapshotStack.pop();
            
            // action required...
            if (currentSnapshot.k < graph_stratification.size()){
                
                // gather data
                int N = currentSnapshot.flux[currentSnapshot.k];
                int n = graph_stratification[currentSnapshot.k][0].size();
                std::vector<int> minima, maxima;
                for (int j = 0; j < n; j++){
                    int vertex_number = graph_stratification[currentSnapshot.k][0][j];
                    int number_of_attached_edges = graph_stratification[currentSnapshot.k][1][j];
                    int remaining_edges = graph_stratification[currentSnapshot.k][2][j];
                    int min = number_of_attached_edges;
                    int f_other = currentSnapshot.flux[vertex_number];
                    if (min < number_of_attached_edges * root - (f_other - remaining_edges)){
                        min = number_of_attached_edges * root - (f_other - remaining_edges);
                    }
                    int max = number_of_attached_edges * (root-1);
                    minima.push_back(min);
                    maxima.push_back(max);
                }
                
                // compute partitions
                std::vector<std::vector<int>> partitions;
                my_partitions(N, n, minima, maxima, partitions);
                
                // create new snapshots
                std::vector<int> number_of_edges = graph_stratification[currentSnapshot.k][1];
                for(int j = 0; j < partitions.size(); j++){
                    
                    // create data of new snapshot (in particular the number of subpartitions)
                    boost::multiprecision::int128_t mult = currentSnapshot.mult;
                    std::vector<int> new_flux(currentSnapshot.flux.begin(), currentSnapshot.flux.end());
                    new_flux[currentSnapshot.k] = 0;
                    for (int a = 0; a < n; a++){
                        int index = graph_stratification[currentSnapshot.k][0][a];
                        new_flux[index] = new_flux[index] - (root * graph_stratification[currentSnapshot.k][1][a] - partitions[j][a]);
                        mult = mult * number_partitions(partitions[j][a], number_of_edges[a], root);
                    }
                    
                    // add snapshot
                    comb_data newSnapshot;
                    newSnapshot.flux = new_flux;
                    newSnapshot.mult = mult;
                    newSnapshot.k = currentSnapshot.k + 1;
                    snapshotStack.push(newSnapshot);
                    
                }
                
            }
            // no action required -> increase total
            else{
                boost::multiprecision::int128_t mult = currentSnapshot.mult;
                for (int j = 0; j < genera.size(); j++){
                    if ((genera[j] == 1) and (outfluxes[i].partition[j] == 0)){
                        mult = mult * (boost::multiprecision::int128_t) (root * root - 1);
                    }
                    if ((genera[j] == 1) and (outfluxes[i].partition[j] > 0)){
                        mult = mult * (boost::multiprecision::int128_t) (root * root);
                    }
                }
                total += mult;
            }
            
        }
        
    }
    
    // return the result
    return total;
    
}
