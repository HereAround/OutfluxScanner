#include "combinatorics.cpp"


// Thread-safe addition to the result
void UpdateCountThreadSafe(boost::multiprecision::int128_t & central, boost::multiprecision::int128_t & change, const bool & display_details)
{    
    boost::mutex::scoped_lock lock(myGuard);
    central = central + change;
    if (display_details){
        std::cout << "Worker complete: " << change << "\n";
    }
}


// Worker thread for parallel run
void worker(
                                const std::vector<int> degrees,
                                const std::vector<int> genera,
                                const std::vector<std::vector<int>> edges,
                                const int root,
                                const std::vector<std::vector<std::vector<int>>> graph_stratification,
                                const std::vector<int> edge_numbers,
                                const std::vector<std::vector<int>> partitions,
                                boost::multiprecision::int128_t & sum,
                                const bool & display_details )
{
    
    // save total number of roots found
    boost::multiprecision::int128_t total = 0;
    
    // (1) Find fluxes corresponding to partitions
    // (1) Find fluxes corresponding to partitions
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
    
    
    // (2) Count weight assignments
    // (2) Count weight assignments
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
                
                // compute flux_partitions
                std::vector<std::vector<int>> flux_partitions;
                comp_partitions(N, n, minima, maxima, flux_partitions);
                
                // create new snapshots
                std::vector<int> number_of_edges = graph_stratification[currentSnapshot.k][1];
                for(int j = 0; j < flux_partitions.size(); j++){
                    
                    // create data of new snapshot (in particular the number of subpartitions)
                    boost::multiprecision::int128_t mult = currentSnapshot.mult;
                    std::vector<int> new_flux(currentSnapshot.flux.begin(), currentSnapshot.flux.end());
                    new_flux[currentSnapshot.k] = 0;
                    for (int a = 0; a < n; a++){
                        int index = graph_stratification[currentSnapshot.k][0][a];
                        new_flux[index] = new_flux[index] - (root * graph_stratification[currentSnapshot.k][1][a] - flux_partitions[j][a]);
                        mult = mult * number_partitions(flux_partitions[j][a], number_of_edges[a], root);
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
    
    // Update result
    UpdateCountThreadSafe(sum, total, display_details);
    
}



// Count number of root bundles with prescribed number of sections
boost::multiprecision::int128_t parallel_root_counter(
                                const std::vector<int> degrees,
                                const std::vector<int> genera,
                                const std::vector<std::vector<int>> edges,
                                const int root,
                                const std::vector<std::vector<std::vector<int>>> graph_stratification,
                                const std::vector<int> edge_numbers,
                                const int & h0_value,
                                const int & thread_number,
                                const bool & display_details )
{
    
    // check input
    if (thread_number <= 0 or thread_number > 100){
        std::cout << "Corrupted input\n";
        return -1;
    }
    
    // (1) Partition h0
    // (1) Partition h0
    std::vector<std::vector<int>> partitions;
    comp_partitions(h0_value, degrees.size(), std::vector<int>(degrees.size(),0), std::vector<int>(degrees.size(),h0_value), partitions);
    
    // (2) Split the partitions into as many packages as determined by thread_number and start the threads
    // (2) Split the partitions into as many packages as determined by thread_number and start the threads
    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    boost::multiprecision::int128_t sum = (boost::multiprecision::int128_t) 0;
    if (thread_number > 1){
        boost::thread_group threadList;
        int package_size = (int) partitions.size()/thread_number;
        if (display_details){
            std::cout << "Computing in " << thread_number << " parallel threads (average load: " << package_size << ")...\n";
        }
        for (int i = 0; i < thread_number; i++)
        {
            if (i < thread_number - 1){
                std::vector<std::vector<int>> partial_partitions(partitions.begin() + i * package_size, partitions.begin() + (i+1) * package_size);
                boost::thread *t = new boost::thread(worker, degrees, genera, edges, root, graph_stratification, edge_numbers, partial_partitions, boost::ref(sum), boost::ref(display_details));
                threadList.add_thread(t);
            }
            else{
                std::vector<std::vector<int>> partial_partitions(partitions.begin() + i * package_size, partitions.end());
                boost::thread *t = new boost::thread(worker, degrees, genera, edges, root, graph_stratification, edge_numbers, partial_partitions, boost::ref(sum), boost::ref(display_details));
                threadList.add_thread(t);
            }
        }
        threadList.join_all();
    }
    else if (thread_number == 1){
        if (display_details){
            std::cout << "Computing in one thread...\n";
        }
        worker(degrees, genera, edges, root, graph_stratification, edge_numbers, partitions, boost::ref(sum), boost::ref(display_details));
    }
    std::chrono::steady_clock::time_point later = std::chrono::steady_clock::now();
    
    // inform about the result
    if (display_details){
        std::cout << "\nTime for run: " << std::chrono::duration_cast<std::chrono::seconds>(later - now).count() << "[s]\n";
    }
    std::cout << "Total: " << sum << "\n\n";
    return sum;
    
}
