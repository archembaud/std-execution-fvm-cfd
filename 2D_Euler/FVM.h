#define NX 256
#define NY 256
#define N (NX*NY)
#define R 1.0
#define GAMMA 1.4
#define CV (R/(GAMMA-1.0))
#define L 1.0
#define H 1.0
#define DX (L/NX)
#define DY (H/NY)
#define CFL 0.25
// CFL = 2.0*DT/DX
#define DT (CFL*DX/2.0)
#define DT_ON_DX (0.5*CFL)
#define DT_ON_DY (0.5*CFL)   // Assume equal; this is cheating
#define TOTAL_TIME 0.2

// Wrapping Functions called from the main loop

void ComputeConservedFromPrimitives(std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2, std::vector<float>&p3, 
    std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2, std::vector<float>&u3);

void ComputeFluxesFromPrimitives(const int direction, 
    std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2, std::vector<float>&p3,
    std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2,  std::vector<float>&u3,
    std::vector<float>&Fp0, std::vector<float>&Fp1, std::vector<float>&Fp2, std::vector<float>&Fp3,
    std::vector<float>&Fm0, std::vector<float>&Fm1, std::vector<float>&Fm2, std::vector<float>&Fm3);

void ComputeConservedChangeFromFluxes(const int direction, std::vector<float>&Fp0, std::vector<float>&Fp1, std::vector<float>&Fp2, std::vector<float>&Fp3,
    std::vector<float>&Fm0, std::vector<float>&Fm1, std::vector<float>&Fm2, std::vector<float>&Fm3, 
    std::vector<float>&du0, std::vector<float>&du1, std::vector<float>&du2, std::vector<float>&du3);

void UpdateConservedQuantitiesFromdU(const int update_primitives,
    std::vector<float>&du0, std::vector<float>&du1, std::vector<float>&du2, std::vector<float>&du3,
    std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2, std::vector<float>&u3,
    std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2, std::vector<float>&p3);


// Kernel functoins called from within wrapping functions

void ComputeAllUFromP(float& elem,
    const std::vector<float>& density, const std::vector<float>& xvel, const std::vector<float>& yvel,const std::vector<float>& temp,
    std::vector<float>& mass, std::vector<float>& xmom, std::vector<float>& ymom, std::vector<float>& eng);

void ComputeFluxesFromP(float& elem, int direction, const std::vector<float>& density, const std::vector<float>& xvel, const std::vector<float>& yvel, const std::vector<float>& temp,
    std::vector<float>& mass, std::vector<float>& xmom, std::vector<float>& ymom, std::vector<float>& eng,
    std::vector<float>& PmassFlux, std::vector<float>& PxmomFlux, std::vector<float>& PymomFlux, std::vector<float>& PengFlux,
    std::vector<float>& MmassFlux, std::vector<float>& MxmomFlux, std::vector<float>& MymomFlux, std::vector<float>& MengFlux);

void ComputeDeltaUFromP(float& elem, const int direction,
    const std::vector<float>& PmassFlux, const std::vector<float>& PxmomFlux, const std::vector<float>& PymomFlux, const std::vector<float>& PengFlux,
    const std::vector<float>& MmassFlux, const std::vector<float>& MxmomFlux, const std::vector<float>& MymomFlux, const std::vector<float>& MengFlux,
    std::vector<float>& dmass, std::vector<float>& dxmom, std::vector<float>& dymom, std::vector<float>& deng);
        
void UpdateConservedQuantities(float& elem, const int update_primitives,
    std::vector<float>& dmass, std::vector<float>& dxmom, std::vector<float>& dymom, std::vector<float>& deng,
    std::vector<float>& mass, std::vector<float>& xmom, std::vector<float>& ymom, std::vector<float>& eng, 
    std::vector<float>& density, std::vector<float>& xvel, std::vector<float>& yvel, std::vector<float>& temp);