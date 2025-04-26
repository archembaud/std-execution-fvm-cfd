#include <execution>
#include <vector>
#include <iostream>
#include "FVM.h"

#define N 200
#define CV (1/0.4)
#define L 1.0
#define DX (L/N)
#define CFL 0.25
// CFL = 2.0*DT/DX
#define DT (CFL*DX/2.0)
#define NO_STEPS 100

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

    // Call Init
    Init(p0, p1, p2);

    // Compute Conserved Quantities from Primitives
    ComputeConservedFromPrimitives(p0, p1, p2, u0, u1, u2);

    for (int step = 0; step < NO_STEPS; step++) {
        std::cout << "Timestep " << step << " of " << NO_STEPS << "\n";
        ComputeFluxesFromPrimitives(p0, p1, p2, u0, u1, u2,
                                    Fp0, Fp1,Fp2, Fm0, Fm1, Fm2);

        ComputeConservedChangeFromFluxes(Fp0, Fp1, Fp2, Fm0, Fm1, Fm2,
                                    du0, du1, du2);
    }


    // Print a few of the first elements
    for (int i = 98; i < 104; i++) {
        std::cout << "---------------------------------------\n";
        std::cout << "du[" << i << "]=" << du0[i] << ", " << du1[i] << ", " << du2[i] << "\n"; 
        std::cout << "FP[" << i << "]=" << Fp0[i] << ", " << Fp1[i] << ", " << Fp2[i] << "\n"; 
        std::cout << "FM[" << i << "]=" << Fm0[i] << ", " << Fm1[i] << ", " << Fm2[i] << "\n"; 
    }

    return 0;
}