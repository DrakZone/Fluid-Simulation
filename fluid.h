
class Fluid {
public:
    Fluid(float density, int numX, int numY, float h);

    void integrate(float dt, float gravity);
    void solveIncompressibility(int numIters, float dt);
    void extrapolate();
    float sampleField(float x, float y, int field);
    float avgU(int i, int j);
    float avgV(int i, int j);
    void advectVel(float dt);
    void advectSmoke(float dt);
    void simulate(float dt, float gravity, int numIters);

private:
    int numX, numY, numCells;
    float h;
    float density;
    std::vector<float> u, v, newU, newV, p, s, m, newM;

    static constexpr int U_FIELD = 0;
    static constexpr int V_FIELD = 1;
    static constexpr int S_FIELD = 2;
};