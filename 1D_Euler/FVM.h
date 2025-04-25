
float ComputeMassFromP(float density);

float ComputeMomFromP(float density, float xvel);

void ComputeEngFromP(float& elem, const std::vector<float>& density, 
                const std::vector<float>& xvel, 
                const std::vector<float>& temp,
                std::vector<float>& eng);

void ComputeConservedFromPrimitives(std::vector<float>&p0, std::vector<float>&p1, std::vector<float>&p2, 
                                    std::vector<float>&u0, std::vector<float>&u1, std::vector<float>&u2);

void ComputeAllUFromP(float& elem, const std::vector<float>& density, const std::vector<float>& xvel, 
                      const std::vector<float>& temp, std::vector<float>& mass, std::vector<float>& mom,
                      std::vector<float>& eng);
