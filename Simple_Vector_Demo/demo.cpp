#include <execution>
#include <vector>
#include <iostream>
#include <thread>

// A simple function passed into std::execution::par
int ComputeVectorMultiple(int x) {
    return 2*x;
}  

// A slightly more complicated function passed into std::execution::par with two input vectors
int ComputeVectorDot(int x, int y) {
    return x*y;
}  

int main() {
    std::vector<int> data = {1, 2, 3, 4, 5};
    std::vector<int> more_data = {2, 4, 6, 8, 10};
    std::vector<int> result(data.size());

    unsigned int nThreads = std::thread::hardware_concurrency();
    std::cout << "Number of available threads: " << nThreads << "\n";

    // Parallel execution using std::execution::par with one input
    std::transform(std::execution::par, data.begin(), data.end(), result.begin(), ComputeVectorMultiple);

    // Print the result
    for (int i = 0; i < 5; i++) {
        std::cout << "Value of result[" << i << "] is " << result[i] << "\n"; 
    }

    // Parallel execution using std::execution::par with two inputs
    std::transform(std::execution::par, data.begin(), data.end(), more_data.begin(), result.begin(), ComputeVectorDot);

    // Print the result
    for (int i = 0; i < 5; i++) {
        std::cout << "Value of result[" << i << "] is " << result[i] << "\n"; 
    }

    return 0;
}
