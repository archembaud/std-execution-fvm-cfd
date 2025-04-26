
float ComputeMassFromP(float density);

float ComputeMomFromP(float density, float xvel);

void ComputeEngFromP(float& elem, const std::vector<float>& density, 
                const std::vector<float>& xvel, 
                const std::vector<float>& temp,
                std::vector<float>& eng);

void ComputeConservedFromPrimitives(std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2, 
                                    std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2);

void ComputeFluxesFromPrimitives(std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2,
                                std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2,
                                std::vector<float>&Fp0, std::vector<float>&Fp1, std::vector<float>&Fp2,
                                std::vector<float>&Fm0, std::vector<float>&Fm1, std::vector<float>&Fm2);


void ComputeAllUFromP(float& elem, const std::vector<float>& density, const std::vector<float>& xvel, 
                      const std::vector<float>& temp, std::vector<float>& mass, std::vector<float>& mom,
                      std::vector<float>& eng);

void ComputeFluxesFromP(float& elem, const std::vector<float>& density, const std::vector<float>& xvel, const std::vector<float>& temp,
                        std::vector<float>& mass, std::vector<float>& mom, std::vector<float>& eng,
                        std::vector<float>& PmassFlux, std::vector<float>& PmomFlux, std::vector<float>& PengFlux,
                        std::vector<float>& MmassFlux, std::vector<float>& MmomFlux, std::vector<float>& MengFlux);
