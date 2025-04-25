#include <execution>
#include <vector>
#include "FVM.h"

#define R 1.0
#define GAMMA 1.4
#define CV (R/(GAMMA-1.0))

// A slightly more complicated function passed into std::execution::par with two input vectors
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
