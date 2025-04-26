#include <execution>
#include <vector>
#include <iostream>
#include <fstream>
#include "FVM.h"

void Init(std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2) {
    for (int i = 0; i < N; i++) {
        if (i < 0.5*N) {
            p0[i] = 10.0;
            p1[i] = 0.0;
            p2[i] = 1.0;
        } else {
            p0[i] = 1.0;
            p1[i] = 0.0;
            p2[i] = 1.0;    
        }
    }
}

void Save_Results(const std::vector<float>&p0, const std::vector<float>&p1, const std::vector<float>&p2)   {
    std::ofstream ResultFile("results.txt");
    for (int i = 0; i < N; i++) {
        float cx = (i+0.5)*DX;
        ResultFile << cx << "\t" << p0[i] << "\t" << p1[i] << "\t" << p2[i] << "\n";
    }
    ResultFile.close();
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
    std::vector<float> du0(problemSize);
    std::vector<float> du1(problemSize);
    std::vector<float> du2(problemSize);
    // Fluxes
    std::vector<float> Fp0(problemSize);
    std::vector<float> Fp1(problemSize);
    std::vector<float> Fp2(problemSize);
    std::vector<float> Fm0(problemSize);
    std::vector<float> Fm1(problemSize);
    std::vector<float> Fm2(problemSize);
    float time = 0.0;
    int step = 0;
    // Call Init
    Init(p0, p1, p2);

    // Compute Conserved Quantities from Primitives
    ComputeConservedFromPrimitives(p0, p1, p2, u0, u1, u2);

    while (time < TOTAL_TIME) {
        ComputeFluxesFromPrimitives(p0, p1, p2, u0, u1, u2,
                                    Fp0, Fp1,Fp2, Fm0, Fm1, Fm2);

        ComputeConservedChangeFromFluxes(Fp0, Fp1, Fp2, Fm0, Fm1, Fm2,
                                    du0, du1, du2);
        
        UpdateConservedQuantitiesFromdU(du0, du1, du2, u0, u1, u2, p0, p1, p2);

        // Update time
        time+=DT;
        step++;
    }

    // Save the results
    Save_Results(p0, p1, p2);

    std::cout << "Simulation completed after " << step << " steps with " << N << " cells\n";
    return 0;
}