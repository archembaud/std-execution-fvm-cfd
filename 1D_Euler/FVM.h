#define N 2048
#define R 1.0
#define GAMMA 1.4
#define CV (R/(GAMMA-1.0))
#define L 1.0
#define DX (L/N)
#define CFL 0.25
// CFL = 2.0*DT/DX
#define DT (CFL*DX/2.0)
#define DT_ON_DX (0.5*CFL)
#define TOTAL_TIME 0.2

float ComputeMassFromP(float density);

float ComputeMomFromP(float density, float xvel);

void ComputeEngFromP(float& elem, const std::vector<float>& density, 
                const std::vector<float>& xvel, 
                const std::vector<float>& temp,
                std::vector<float>& eng);

// Wrapping Functions called from the main loop

void ComputeConservedFromPrimitives(std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2, 
                                    std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2);

void ComputeFluxesFromPrimitives(std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2,
                                std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2,
                                std::vector<float>&Fp0, std::vector<float>&Fp1, std::vector<float>&Fp2,
                                std::vector<float>&Fm0, std::vector<float>&Fm1, std::vector<float>&Fm2);

void ComputeConservedChangeFromFluxes(std::vector<float>&Fp0, std::vector<float>&Fp1, std::vector<float>&Fp2,
                                std::vector<float>&Fm0, std::vector<float>&Fm1, std::vector<float>&Fm2,
                                std::vector<float>&du0, std::vector<float>&du1, std::vector<float>&du2);

void UpdateConservedQuantitiesFromdU(std::vector<float>&du0, std::vector<float>&du1, std::vector<float>&du2,
                                std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2,
                                std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2);


// Kernel functoins called from within wrapping functions

void ComputeAllUFromP(float& elem, const std::vector<float>& density, const std::vector<float>& xvel, 
                      const std::vector<float>& temp, std::vector<float>& mass, std::vector<float>& mom,
                      std::vector<float>& eng);

void ComputeFluxesFromP(float& elem, const std::vector<float>& density, const std::vector<float>& xvel, const std::vector<float>& temp,
                        std::vector<float>& mass, std::vector<float>& mom, std::vector<float>& eng,
                        std::vector<float>& PmassFlux, std::vector<float>& PmomFlux, std::vector<float>& PengFlux,
                        std::vector<float>& MmassFlux, std::vector<float>& MmomFlux, std::vector<float>& MengFlux);

void ComputeDeltaUFromP(float& elem,
                        const std::vector<float>& PmassFlux, const std::vector<float>& PmomFlux, const std::vector<float>& PengFlux,
                        const std::vector<float>& MmassFlux, const std::vector<float>& MmomFlux, const std::vector<float>& MengFlux,
                        std::vector<float>& dmass, std::vector<float>& dmom, std::vector<float>& deng);

void UpdateConservedQuantities(float& elem, const std::vector<float>& dmass, const std::vector<float>& dmom, const std::vector<float>& deng,
                                std::vector<float>& mass, std::vector<float>& mom, std::vector<float>& eng, 
                                std::vector<float>& density, std::vector<float>& xvel, std::vector<float>& temp);