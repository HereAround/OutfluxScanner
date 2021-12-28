// Task: Compute partitions of an integer N into a sum of n integers with specified minima and maxima.
void my_partitions(const int & N, const int & n, const std::vector<int> & minima, const std::vector<int> & maxima, std::vector<std::vector<int>> & partitions){
    
    // Only one value to set?
    if(n == 1) {
        if (minima[0] <= N and N <= maxima[0]){
            partitions.push_back({N});
        }
        return;
    }
    
    // Create stack
    struct SnapShotStruct {
        std::vector<int> p;
    };
    std::stack<SnapShotStruct> snapshotStack;
    
    // Add first snapshot
    SnapShotStruct currentSnapshot;
    currentSnapshot.p = {};
    snapshotStack.push(currentSnapshot);
    
    // Run...
    while(!snapshotStack.empty())
    {
        
        // pick the top snapshot and delete it from the stack
        currentSnapshot= snapshotStack.top();
        snapshotStack.pop();
        
        // current sum
        int current_sum = std::accumulate(currentSnapshot.p.begin(), currentSnapshot.p.end(), 0);

        // position at which we add
        int pos = currentSnapshot.p.size();            
        
        // there are at least two more values to be set
        if (pos < n-1){
        
            // Compute min and max for this placement
            int max = maxima[pos];
            if (N - current_sum < max ){
                max = N - current_sum;
            }
            
            // Iterate over these values
            for (int i = minima[pos]; i <= max; i++){
                std::vector<int> new_partition = currentSnapshot.p;
                new_partition.push_back(i);
                SnapShotStruct newSnapshot;
                newSnapshot.p = new_partition;
                snapshotStack.push(newSnapshot);
            }
            
        }
        
        // only one more value to be set
        else if ((pos == n-1) && (minima[pos] <= N - current_sum) && (N - current_sum <= maxima[pos])){
            std::vector<int> new_partition = currentSnapshot.p;
            new_partition.push_back(N - current_sum);
            partitions.push_back(new_partition);
        }
        
    }
    
}



// Task: Compute number of partitions of an integer f.
// Input: Integer f to be partitioned.
//           Integers r, n.
// Output: The number of partitions of f into a sum of exactly n integers w1, ... wn with 1 <= w1, ..., wn < r.
boost::multiprecision::int128_t number_partitions( const int & f, const int & n, const int & r ){
    boost::multiprecision::int128_t count = (boost::multiprecision::int128_t) 0;
    
    // Only one value to set?
    if( n == 1 ) {
        
        // Check if we have a partition
        if ( ( 1 <= f ) && ( f < r ) ){
            return (boost::multiprecision::int128_t) 1;
        }
        else{
            return (boost::multiprecision::int128_t) 0;
        }
    
    }
    
    // Pick values and make recursive call
    for ( int i = 1; i < r; i++ ){
        count = count + number_partitions( f - i, n-1, r );
    }
    
    // return final result
    return count;    
}



// Task: Compute partitions of an integer N into a sum of n integers with specified minima and maxima.
//           Each number in this partition represents a flux among two components, connected by n_i edges. There, we also have freedom to distribute.
//           This leads to further combinatorics prefactors.
void my_partitions_v2(
                      const int & N, 
                      const int & n, 
                      const std::vector<int> & edge_numbers, 
                      const std::vector<int> & minima, 
                      const std::vector<int> & maxima, 
                      const int & root, 
                      std::vector<std::vector<int>> & partitions,
                      std::vector<boost::multiprecision::int128_t> & subpartitions
                     )
{
    //std::cout << "Test1\n";
    
    // Only one value to set?
    if(n == 1) {
        if (minima[0] <= N and N <= maxima[0]){
            partitions.push_back({N});
        }
        return;
    }
    
    // Create stack
    struct SnapShotStruct {
        std::vector<int> p;
    };
    std::stack<SnapShotStruct> snapshotStack;
    
    // Add first snapshot
    SnapShotStruct currentSnapshot;
    currentSnapshot.p = {};
    snapshotStack.push(currentSnapshot);
    
    // Run...
    while(!snapshotStack.empty())
    {
        
        // pick the top snapshot and delete it from the stack
        currentSnapshot= snapshotStack.top();
        snapshotStack.pop();
        
        // current sum
        int current_sum = std::accumulate(currentSnapshot.p.begin(), currentSnapshot.p.end(), 0);

        // position at which we add
        int pos = currentSnapshot.p.size();            
        
        // there are at least two more values to be set
        if (pos < n-1){
        
            // Compute min and max for this placement
            int max = maxima[pos];
            if (N - current_sum < max ){
                max = N - current_sum;
            }
            
            // Iterate over these values
            for (int i = minima[pos]; i <= max; i++){
                std::vector<int> new_partition = currentSnapshot.p;
                new_partition.push_back(i);
                SnapShotStruct newSnapshot;
                newSnapshot.p = new_partition;
                snapshotStack.push(newSnapshot);
            }
            
        }
        
        // only one more value to be set
        else if ((pos == n-1) && (minima[pos] <= N - current_sum) && (N - current_sum <= maxima[pos])){
            // set final value for partition
            std::vector<int> new_partition = currentSnapshot.p;
            new_partition.push_back(N - current_sum);
            partitions.push_back(new_partition);
            
            // compute number of subpartitions
            /*boost::multiprecision::int128_t mult = 1;
            if (new_partition.size() != edge_numbers.size()){
                std::cout << "\n\nERROR:";
                std::cout << "Edge number: ";
                for (int k = 1; k < edge_numbers.size(); k++){
                    std::cout << edge_numbers[k] << ",";
                }
                std::cout << "\n";
                for (int k = 1; k < new_partition.size(); k++){
                    std::cout << new_partition[k] << ",";
                }
                std::cout << "\n\n";
            }
            for (int k = 1; k < edge_numbers.size(); k++){
                mult = mult * number_partitions(new_partition[k], edge_numbers[k], root);
            }
            subpartitions.push_back(mult);*/
        }
        
    }
    
    //std::cout << "Test2\n";
}
