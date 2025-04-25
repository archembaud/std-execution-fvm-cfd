#include <execution>
#include <vector>
#include <iostream>
#include "FVM.h"

#define N 100
#define CV (1/0.4)

void Init(std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2) {
    for (int i = 0; i < N; i++) {\
        p0[i] = 1.0;
        p1[i] = 0.5;
        p2[i] = 0.25;
    }
    // Debugging
    printf("Value for E = %g\n", p0[0]*(CV*p2[0] + 0.5*p1[0]*p1[0]));
}

int main() {
    
    // Set up large vectors for 1D Euler Computation
    size_t problemSize = N;
    // Primitives
    std::vector<float> p0(problemSize);
    std::vector<float> p1(problemSize);
    std::vector<float> p2(problemSize);
    // Conserved Quantities
    std::vector<float> u0(problemSize);
    std::vector<float> u1(problemSize);
    std::vector<float> u2(problemSize);

    // Call Init
    Init(p0, p1, p2);

    // Compute Conserved Quantities from Primitives
    ComputeConservedFromPrimitives(p0, p1, p2, u0, u1, u2);

    // Print a few of the first elements
    for (int i = 0; i < 5; i++) {
        std::cout << "---------------------------------------";
        std::cout << "u0[" << i << "]=" << u0[i] << ", p0 =" << p0[i] << "\n"; 
        std::cout << "u1, u2  [" << i << "] = " << u1[i] << "," << u2[i] << "\n"; 
    }

    return 0;
}