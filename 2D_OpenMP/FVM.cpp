#include <execution>
#include <vector>
#include <cmath>
#include "FVM.h"

// Wrapping Functions called from the main loop

void ComputeFluxesFromPrimitives(const int direction, std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2, std::vector<float>&p3,
                                 std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2,  std::vector<float>&u3,
                                 std::vector<float>&Fp0, std::vector<float>&Fp1, std::vector<float>&Fp2, std::vector<float>&Fp3,
                                 std::vector<float>&Fm0, std::vector<float>&Fm1, std::vector<float>&Fm2, std::vector<float>&Fm3) {

    /*
    Compute the split SHLL fluxes - forward (p) and backward (m) for each cell in parallel using std::execution
    instead of the more conventional (in HPC circles) OpenMP.
    */
    #pragma omp parallel for
    for (int index = 0; index < N; index++) {
        ComputeFluxesFromP(index, direction, p0, p1, p2, p3, u0, u1, u2, u3, Fp0, Fp1, Fp2, Fp3, Fm0, Fm1, Fm2, Fm3);
    }
}



void ComputeConservedChangeFromFluxes(const int direction, std::vector<float>&Fp0, std::vector<float>&Fp1, std::vector<float>&Fp2, std::vector<float>&Fp3,
                                      std::vector<float>&Fm0, std::vector<float>&Fm1, std::vector<float>&Fm2, std::vector<float>&Fm3, 
                                      std::vector<float>&du0, std::vector<float>&du1, std::vector<float>&du2, std::vector<float>&du3) {
    /*
    Compute the change in conserved quantites based on the fluxes.
    */
    #pragma omp parallel for
    for (int index = 0; index < N; index++) {
        ComputeDeltaUFromP(index, direction, Fp0, Fp1, Fp2, Fp3, Fm0, Fm1, Fm2, Fm3, du0, du1, du2, du3);
    }

}


void UpdateConservedQuantitiesFromdU(const int update_primitives,
                                std::vector<float>&du0, std::vector<float>&du1, std::vector<float>&du2, std::vector<float>&du3,
                                std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2, std::vector<float>&u3,
                                std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2, std::vector<float>&p3) {
    /*
    Update the value of u based on dU, then update P values as well while we are there.
    */
    #pragma omp parallel for
    for (int index = 0; index < N; index++) {
        UpdateConservedQuantities(index, update_primitives, du0, du1, du2, du3, u0, u1, u2, u3, p0, p1, p2, p3);
    }
}



void ComputeConservedFromPrimitives(std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2, std::vector<float>&p3, 
                                    std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2, std::vector<float>&u3) { 

    #pragma omp parallel for
    for (int index = 0; index < N; index++) {
        ComputeAllUFromP(index, p0, p1, p2, p3, u0, u1, u2, u3);
    }
}

// Kernel functions called from within wrapping functions

void ComputeAllUFromP(int index, 
    const std::vector<float>& density, const std::vector<float>& xvel, const std::vector<float>& yvel, const std::vector<float>& temp, 
    std::vector<float>& mass, std::vector<float>& xmom, std::vector<float>& ymom, std::vector<float>& eng) {   
    // Compute the energy based on density, velocity and temperature
    mass[index] = density[index];
    xmom[index] = density[index]*xvel[index];
    ymom[index] = density[index]*yvel[index];
    eng[index] = density[index] * (CV*temp[index] + 0.5f * (xvel[index] * xvel[index] + yvel[index] * yvel[index]));
}


void ComputeDeltaUFromP(int index, const int direction,
                        const std::vector<float>& PmassFlux, const std::vector<float>& PxmomFlux, const std::vector<float>& PymomFlux, const std::vector<float>& PengFlux,
                        const std::vector<float>& MmassFlux, const std::vector<float>& MxmomFlux, const std::vector<float>& MymomFlux, const std::vector<float>& MengFlux,
                        std::vector<float>& dmass, std::vector<float>& dxmom, std::vector<float>& dymom, std::vector<float>& deng) {

    // Compute the change in mass by computing net fluxes
    // Compute the left net fluxes first
    float NetFluxLeft[4];
    float NetFluxRight[4];

    int xcell = (int)index/NY;
    int ycell = index - xcell*NY;

    int N_edge;
    float PHI;
    int direction_index;
    int neighbour_increment;
    if (direction == 0) {
        N_edge = NX;
        PHI = DT_ON_DX;
        direction_index = xcell;
        neighbour_increment = NY;
    } else {
        N_edge = NY;
        PHI = DT_ON_DY;
        direction_index = ycell;
        neighbour_increment = 1;
    }

    if (direction_index == 0) {
        // First cell, treat left hand boundary as outflow
        NetFluxLeft[0] = PmassFlux[index] + MmassFlux[index];
        NetFluxLeft[1] = PxmomFlux[index] + MxmomFlux[index];
        NetFluxLeft[2] = PymomFlux[index] + MymomFlux[index];
        NetFluxLeft[3] = PengFlux[index] + MengFlux[index];
        // Right Flux
        NetFluxRight[0] = MmassFlux[index+neighbour_increment] + PmassFlux[index];
        NetFluxRight[1] = MxmomFlux[index+neighbour_increment] + PxmomFlux[index];
        NetFluxRight[2] = MymomFlux[index+neighbour_increment] + PymomFlux[index];
        NetFluxRight[3] = MengFlux[index+neighbour_increment] + PengFlux[index];
    } else if (direction_index == (N_edge-1)) {
        // Last cell
        NetFluxLeft[0] = PmassFlux[index-neighbour_increment] + MmassFlux[index];
        NetFluxLeft[1] = PxmomFlux[index-neighbour_increment] + MxmomFlux[index];
        NetFluxLeft[2] = PymomFlux[index-neighbour_increment] + MymomFlux[index];
        NetFluxLeft[3] = PengFlux[index-neighbour_increment] + MengFlux[index];
        // Right Flux
        NetFluxRight[0] = MmassFlux[index] + PmassFlux[index];
        NetFluxRight[1] = MxmomFlux[index] + PxmomFlux[index];
        NetFluxRight[2] = MymomFlux[index] + PymomFlux[index];
        NetFluxRight[3] = MengFlux[index] + PengFlux[index];
    } else {
        // Interior Cell
        // Left
        NetFluxLeft[0] = PmassFlux[index-neighbour_increment] + MmassFlux[index];
        NetFluxLeft[1] = PxmomFlux[index-neighbour_increment] + MxmomFlux[index];
        NetFluxLeft[2] = PymomFlux[index-neighbour_increment] + MymomFlux[index];
        NetFluxLeft[3] = PengFlux[index-neighbour_increment] + MengFlux[index];
        // Right Flux
        NetFluxRight[0] = MmassFlux[index+neighbour_increment] + PmassFlux[index];
        NetFluxRight[1] = MxmomFlux[index+neighbour_increment] + PxmomFlux[index];
        NetFluxRight[2] = MymomFlux[index+neighbour_increment] + PymomFlux[index];
        NetFluxRight[3] = MengFlux[index+neighbour_increment] + PengFlux[index];
    }

    // Compute dU
    dmass[index] = dmass[index] -PHI*(NetFluxRight[0] - NetFluxLeft[0]);
    dxmom[index] = dxmom[index] -PHI*(NetFluxRight[1] - NetFluxLeft[1]);
    dymom[index] = dymom[index] -PHI*(NetFluxRight[2] - NetFluxLeft[2]);
    deng[index] = deng[index] -PHI*(NetFluxRight[3] - NetFluxLeft[3]);
}


void UpdateConservedQuantities(int index, const int update_primitives, std::vector<float>& dmass, std::vector<float>& dxmom, std::vector<float>& dymom, std::vector<float>& deng,
                               std::vector<float>& mass, std::vector<float>& xmom, std::vector<float>& ymom, std::vector<float>& eng, 
                               std::vector<float>& density, std::vector<float>& xvel, std::vector<float>& yvel, std::vector<float>& temp) {

    // Update conserved quantities
    mass[index] = mass[index] + dmass[index];
    xmom[index] = xmom[index] + dxmom[index];
    ymom[index] = ymom[index] + dymom[index];
    eng[index] = eng[index] + deng[index];
    if (update_primitives == 1) {
        density[index] = mass[index];
        xvel[index] = xmom[index]/density[index];
        yvel[index] = ymom[index]/density[index];
        temp[index] = ((eng[index]/density[index]) - 0.5*(xvel[index]*xvel[index]+yvel[index]*yvel[index]))/CV;
        // Reset dU
        dmass[index] = 0.0;
        dxmom[index] = 0.0;
        dymom[index] = 0.0;
        deng[index] = 0.0;
    }
}


void ComputeFluxesFromP(int index, const int direction, const std::vector<float>& density, const std::vector<float>& xvel, const std::vector<float>& yvel, const std::vector<float>& temp,
                        std::vector<float>& mass, std::vector<float>& xmom, std::vector<float>& ymom, std::vector<float>& eng,
                        std::vector<float>& PmassFlux, std::vector<float>& PxmomFlux, std::vector<float>& PymomFlux, std::vector<float>& PengFlux,
                        std::vector<float>& MmassFlux, std::vector<float>& MxmomFlux, std::vector<float>& MymomFlux, std::vector<float>& MengFlux) {

    // Compute Mach number
    float a = sqrt(GAMMA*R*temp[index]);  // Speed of sound
    float Pressure = density[index]*R*temp[index];
    float Mach, massFlux, normalMomFlux, perpMomFlux, engFlux; 
    if (direction == 0) {
        // X direction
        Mach = xvel[index]/a;
        // Z invariants
        float Z1 = 0.5*(Mach + 1.0);
        float Z2 = 0.5*a*(1.0-Mach*Mach);
        float Z3 = 0.5*(Mach - 1.0);

        // Fluxes of conserved quantities
        massFlux = xmom[index];
        normalMomFlux = massFlux*xvel[index] + Pressure;
        perpMomFlux = massFlux*yvel[index];
        engFlux = xvel[index]*(eng[index] + Pressure);

        // Split fluxes - positive (P)
        // FP[:,:,0] = F[:,:,0]*Z1 + U[:,:,0]*Z2
        PmassFlux[index] = massFlux*Z1 + mass[index]*Z2;
        PxmomFlux[index] = normalMomFlux*Z1 + xmom[index]*Z2;
        PymomFlux[index] = perpMomFlux*Z1 + ymom[index]*Z2;
        PengFlux[index] = engFlux*Z1 + eng[index]*Z2;
    
        // Split fluxes - minus (M)
        // FM[:,:,0] = -F[:,:,0]*Z3 - U[:,:,0]*Z2
        MmassFlux[index] = -massFlux*Z3 - mass[index]*Z2;
        MxmomFlux[index] = -normalMomFlux*Z3 - xmom[index]*Z2;
        MymomFlux[index] = -perpMomFlux*Z3 - ymom[index]*Z2;
        MengFlux[index] = -engFlux*Z3 - eng[index]*Z2;

    } else {
        // Y direction
        Mach = yvel[index]/a;

        // Z invariants
        float Z1 = 0.5*(Mach + 1.0);
        float Z2 = 0.5*a*(1.0-Mach*Mach);
        float Z3 = 0.5*(Mach - 1.0);

        // Fluxes of conserved quantities
        massFlux = ymom[index];
        normalMomFlux = massFlux*yvel[index] + Pressure;
        perpMomFlux = massFlux*xvel[index];
        engFlux = yvel[index]*(eng[index] + Pressure);

        // Split fluxes - positive (P)
        // FP[:,:,0] = F[:,:,0]*Z1 + U[:,:,0]*Z2
        PmassFlux[index] = massFlux*Z1 + mass[index]*Z2;
        PxmomFlux[index] = perpMomFlux*Z1 + xmom[index]*Z2;
        PymomFlux[index] = normalMomFlux*Z1 + ymom[index]*Z2;
        PengFlux[index] = engFlux*Z1 + eng[index]*Z2;
    
        // Split fluxes - minus (M)
        // FM[:,:,0] = -F[:,:,0]*Z3 - U[:,:,0]*Z2
        MmassFlux[index] = -massFlux*Z3 - mass[index]*Z2;
        MxmomFlux[index] = -perpMomFlux*Z3 - xmom[index]*Z2;
        MymomFlux[index] = -normalMomFlux*Z3 - ymom[index]*Z2;
        MengFlux[index] = -engFlux*Z3 - eng[index]*Z2;
    }
}
