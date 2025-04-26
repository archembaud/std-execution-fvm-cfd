#include <execution>
#include <vector>
#include <cmath>
#include "FVM.h"

#define R 1.0
#define GAMMA 1.4
#define CV (R/(GAMMA-1.0))

void ComputeFluxesFromPrimitives(std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2,
                                 std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2,
                                 std::vector<float>&Fp0, std::vector<float>&Fp1, std::vector<float>&Fp2,
                                 std::vector<float>&Fm0, std::vector<float>&Fm1, std::vector<float>&Fm2) {

    /*
    Compute the split SHLL fluxes - forward (p) and backward (m) for each cell in parallel using std::execution
    instead of the more conventional (in HPC circles) OpenMP.
    */
    std::for_each(
        std::execution::par_unseq,
        p0.begin(),
        p0.end(),
        [&p0, &p1, &p2, &u0, &u1, &u2, &Fp0, &Fp1, &Fp2, &Fm0, &Fm1, &Fm2](float& elem) {
            ComputeFluxesFromP(elem, p0, p1, p2,  u0, u1, u2, Fp0, Fp1, Fp2, Fm0, Fm1, Fm2);
        }
    );
}

void ComputeConservedFromPrimitives(std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2, 
                                    std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2) { 

    /*
    // Update Conserved Mass in parallel
    std::transform(std::execution::par_unseq, p0.begin(), p0.end(), u0.begin(), ComputeMassFromP);
    // Update Conserved Momentum in parallel
    std::transform(std::execution::par_unseq, p0.begin(), p0.end(), p1.begin(), u1.begin(), ComputeMomFromP);
    // Update Conserved Energy in parallel
    // The std::transform is designed to support a max of two inputs.
    // For 3 inputs (such as energy) we have to use for_each.
    std::for_each(
        std::execution::par_unseq,
        p0.begin(),
        p0.end(),
        [&p0, &p1, &p2, &u2](float& elem) {
            ComputeEngFromP(elem, p0, p1, p2, u2);
        }
    );
    */
    // Instead of 2x std::transform and 1x std::for_each, let's do all the work in a single for_each call.
    // This should minimize the overhead and maximise the parallel work
    std::for_each(
        std::execution::par_unseq,
        p0.begin(),
        p0.end(),
        [&p0, &p1, &p2, &u0, &u1, &u2](float& elem) {
            ComputeAllUFromP(elem, p0, p1, p2, u0, u1, u2);
        }
    );
}


float ComputeMassFromP(float density) {
    // For density, this is trivial => p = u
    return density;
}

float ComputeMomFromP(float density, float xvel) {
    // Momentum is density x xvel
    return density*xvel;
}

void ComputeEngFromP(float& elem, const std::vector<float>& density, const std::vector<float>& xvel, const std::vector<float>& temp, std::vector<float>& eng) {
    size_t index = &elem - &density[0];    
    // Compute the energy based on density, velocity and temperature
    eng[index] = density[index] * (CV*temp[index] + 0.5f * xvel[index] * xvel[index]);
}

void ComputeAllUFromP(float& elem, const std::vector<float>& density, const std::vector<float>& xvel, const std::vector<float>& temp, std::vector<float>& mass, std::vector<float>& mom, std::vector<float>& eng) {
    size_t index = &elem - &density[0];    
    // Compute the energy based on density, velocity and temperature
    mass[index] = density[index];
    mom[index] = density[index]*xvel[index];
    eng[index] = density[index] * (CV*temp[index] + 0.5f * xvel[index] * xvel[index]);
}

void ComputeFluxesFromP(float& elem, const std::vector<float>& density, const std::vector<float>& xvel, const std::vector<float>& temp,
                        std::vector<float>& mass, std::vector<float>& mom, std::vector<float>& eng,
                        std::vector<float>& PmassFlux, std::vector<float>& PmomFlux, std::vector<float>& PengFlux,
                        std::vector<float>& MmassFlux, std::vector<float>& MmomFlux, std::vector<float>& MengFlux) {
    size_t index = &elem - &density[0];
    // Compute Mach number
    float a = sqrt(GAMMA*R*temp[index]);  // Speed of sound
    float Mach = xvel[index]/a;
    // Z invariants
    float Z1 = 0.5*(Mach + 1.0);
    float Z2 = 0.5*a*(1.0-Mach*Mach);
    float Z3 = 0.5*(Mach - 1.0);
    // Pressure using the ideal gas law
    float Pressure = density[index]*R*temp[index];
    // Fluxes of conserved quantities
    float massFlux = mom[index];
    float momFlux = massFlux*xvel[index] + Pressure;
    float engFlux = xvel[index]*(eng[index] + Pressure);

    // Split fluxes - positive (P)
    // FP[:,:,0] = F[:,:,0]*Z1 + U[:,:,0]*Z2
    PmassFlux[index] = massFlux*Z1 + mass[index]*Z2;
    PmomFlux[index] = momFlux*Z1 + mom[index]*Z2;
    PengFlux[index] = engFlux*Z1 + eng[index]*Z2;

    // Split fluxes - minus (M)
    // FM[:,:,0] = -F[:,:,0]*Z3 - U[:,:,0]*Z2
    MmassFlux[index] = -massFlux*Z3 - mass[index]*Z2;
    MmomFlux[index] = -momFlux*Z3 - mom[index]*Z2;
    MengFlux[index] = -engFlux*Z3 - eng[index]*Z2;
}
