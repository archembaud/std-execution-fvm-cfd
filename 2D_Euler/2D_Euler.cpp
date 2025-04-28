#include <execution>
#include <vector>
#include <iostream>
#include <fstream>
#include "FVM.h"

void Init(std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2, std::vector<float>&p3,
          std::vector<float>&du0, std::vector<float>&du1, std::vector<float>&du2, std::vector<float>&du3) {
    for (int i = 0; i < N; i++) {
        int xcell = (int)i/NY;
        int ycell = i - xcell*NY;
        if (ycell < 0.5*NY) {
            p0[i] = 10.0;
            p1[i] = 0.0;
            p2[i] = 0.0;
            p3[i] = 1.0;
        } else {
            p0[i] = 1.0;
            p1[i] = 0.0;
            p2[i] = 0.0;
            p3[i] = 1.0;
        }
        du0[i] = 0.0;
        du1[i] = 0.0;
        du2[i] = 0.0;
        du3[i] = 0.0;
    }
}

void Save_Results(const std::vector<float>&p0, const std::vector<float>&p1, const std::vector<float>&p2, const std::vector<float>&p3)   {
    std::ofstream ResultFile("results.txt");
    for (int i = 0; i < N; i++) {
        int xcell = (int)i/NY;
        int ycell = i - xcell*NY;
        float cx = (xcell+0.5)*DX;
        float cy = (ycell+0.5)*DY;
        ResultFile << cx << "\t" << cy << "\t" << p0[i] << "\t" << p1[i] << "\t" << p2[i] << "\t" << p3[i] << "\n";
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
    std::vector<float> p3(problemSize);
    // Conserved Quantities
    std::vector<float> u0(problemSize);
    std::vector<float> u1(problemSize);
    std::vector<float> u2(problemSize);
    std::vector<float> u3(problemSize);
    std::vector<float> du0(problemSize);
    std::vector<float> du1(problemSize);
    std::vector<float> du2(problemSize);
    std::vector<float> du3(problemSize);
    // X Fluxes
    std::vector<float> Fp0(problemSize);
    std::vector<float> Fp1(problemSize);
    std::vector<float> Fp2(problemSize);
    std::vector<float> Fp3(problemSize);
    std::vector<float> Fm0(problemSize);
    std::vector<float> Fm1(problemSize);
    std::vector<float> Fm2(problemSize);
    std::vector<float> Fm3(problemSize);
    float time = 0.0;
    int step = 0;
    // Call Init
    Init(p0, p1, p2, p3, du0, du1, du2, du3);

    // Compute Conserved Quantities from Primitives
    ComputeConservedFromPrimitives(p0, p1, p2, p3, u0, u1, u2, u3);

    for (int step = 0; step < 200; step++) {
    //while (time < TOTAL_TIME) {
        printf("Taking time step %d\n", step);
        // X direction
        ComputeFluxesFromPrimitives(0, p0, p1, p2, p3, u0, u1, u2, u3,
                                    Fp0, Fp1, Fp2, Fp3, Fm0, Fm1, Fm2, Fm3);
        ComputeConservedChangeFromFluxes(0, Fp0, Fp1, Fp2, Fp3, Fm0, Fm1, Fm2, Fm3,
                                    du0, du1, du2, du3);
        // Y direction
        ComputeFluxesFromPrimitives(1, p0, p1, p2, p3, u0, u1, u2, u3,
             Fp0, Fp1, Fp2, Fp3, Fm0, Fm1, Fm2, Fm3);
        ComputeConservedChangeFromFluxes(1, Fp0, Fp1, Fp2, Fp3, Fm0, Fm1, Fm2, Fm3,
            du0, du1, du2, du3);

        // Update U, recompute all P, reset dU
        UpdateConservedQuantitiesFromdU(1, du0, du1, du2, du3, u0, u1, u2, u3, p0, p1, p2, p3);

        // Update time
        time+=DT;
        //step++;
    }

    // Save the results
    Save_Results(p0, p1, p2, p3);

    std::cout << "Simulation completed after " << step << " steps with " << N << " cells\n";
    return 0;
}