#include <execution>
#include <vector>
#include <iostream>

// A simple function passed into std:execution::par
int ComputeVectorMultiple(int x) {
    return 2*x;
}  


int main() {
    std::vector<int> data = {1, 2, 3, 4, 5};
    std::vector<int> result(data.size());

    // Parallel execution using std::execution::par
    std::transform(std::execution::par, data.begin(), data.end(), result.begin(), ComputeVectorMultiple);

    // Print the result
    for (int i = 0; i < 5; i++) {
        std::cout << "Value of data[" << i << "] is " << result[i] << "\n"; 
    }
}