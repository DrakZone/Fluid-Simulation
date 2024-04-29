#include <cmath>
#include <vector>
#include <thread>
#include "physics.h"


// The IX function is a utility function to convert 2D grid coordinates to a linear array index, ensuring that the coordinates are within the valid range for the grid.
// This is a helper function - does not contribute to the physics 
int IX(int x, int y, int N) {
    if (x < 0) { x = 0; }
    if (x > N - 1) { x = N - 1; }

    if (y < 0) { y = 0; }
    if (y > N - 1) { y = N - 1; }

    return (y * N) + x;
}
// Constructor & Destructor
Physics::Physics() {}
Physics::~Physics() {}

// Set Boundaries
// NOTES: 
// b : Boundary condition parameter.
// x: The array to be updated/solved.
// N: Size of the simulation grid.
void Physics::SetBnd(int b, float x[], int N) {
    // Loops over Horizontal Boundaries ( excluding the first and last colums) 
    for (int i = 1; i < N - 1; i++) {
        x[IX(i, 0, N)] = b == 2 ? -x[IX(i, 1, N)] : x[IX(i, 1, N)];
        x[IX(i, N - 1, N)] = b == 2 ? -x[IX(i, N - 2, N)] : x[IX(i, N - 2, N)];
    }

    // Loops over Vertical Boundaries ( excluding the first and last rows)
    for (int j = 1; j < N - 1; j++) {
        x[IX(0, j, N)] = b == 1 ? -x[IX(1, j, N)] : x[IX(1, j, N)];
        x[IX(N - 1, j, N)] = b == 1 ? -x[IX(N - 2, j, N)] : x[IX(N - 2, j, N)];
    }

    // Set Corner Values
    // Top left corner
    x[IX(0, 0, N)] = 0.5f * (x[IX(1, 0, N)] + x[IX(0, 1, N)] + x[IX(0, 0, N)]);

    // Top right corner
    x[IX(0, N - 1, N)] = 0.5f * (x[IX(1, N - 1, N)] + x[IX(0, N - 2, N)] + x[IX(0, N - 1, N)]);

    // Bottom left corner
    x[IX(N - 1, 0, N)] = 0.5f * (x[IX(N - 2, 0, N)] + x[IX(N - 1, 1, N)] + x[IX(N - 1, 0, N)]);

    // Bottom right corner
    x[IX(N - 1, N - 1, N)] = 0.5f * (x[IX(N - 2, N - 1, N)] + x[IX(N - 1, N - 2, N)] + x[IX(N - 1, N - 1, N)]);
}

// To solve Linear equation
// NOTES: 
// b : Boundary condition parameter.
// x: The array to be updated/solved.
// x0: The input array.
// a: Weighting factor for the neighboring points.
// c: Scaling factor.
// iter: Number of iterations to perform.
// N: Size of the simulation grid.
void Physics::LinSolve(int b, float x[], float x0[], float a, float c, int iter, int N) {
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                x[IX(i, j, N)] = (x0[IX(i, j, N)] + a * (x[IX(i + 1, j, N)] + x[IX(i - 1, j, N)] + x[IX(i, j + 1, N)] + x[IX(i, j - 1, N)] + x[IX(i, j, N)] + x[IX(i, j, N)])) * cRecip;
            }
        }
        this->SetBnd(b, x, N);
    }
}

void Physics::LinSolveThreads(int b, float x[], float x0[], float a, float c, int iter, int N) {
    float cRecip = 1.0 / c;
    std::vector<std::thread> threads;
    for (int k = 0; k < iter; k++) {
        for (int j = 1; j < N - 1; j++) {
            threads.emplace_back([&, j]() {
                for (int i = 1; i < N - 1; i++) {
                    int index = IX(i, j, N);
                    float currentValue = x[index];
                    x[index] = (x0[index] + a * (x[IX(i + 1, j, N)] + x[IX(i - 1, j, N)] + x[IX(i, j + 1, N)] + x[IX(i, j - 1, N)] + 2 * currentValue)) * cRecip;
                }
                });
        }
        for (auto& thread : threads) {
            thread.join();
        }
        threads.clear();
        this->SetBnd(b, x, N);
    }
}


// Diffusion
// Notes:
// b: Boundary condition parameter.
// x: The array representing the quantity to be diffused. It gets updated by the diffusion process.
// x0: The input array representing the initial state of the quantity.
// diff: Diffusion rate or viscosity. It determines how fast the quantity diffuses in the fluid.
// dt: Time step. It represents the discrete time interval between simulation steps.
void Physics::Diffuse(int b, float x[], float x0[], float diff, float dt, int iter, int N) {
    // Calculate the diffusion coefficient 'a' based on time step, diffusion rate, and grid size
    float a = dt * diff * (N - 2) * (N - 2);
    // Use the LinSolve function to iteratively solve the diffusion equation
    this->LinSolve(b, x, x0, a, 1 + 6 * a, iter, N);
    // Use Thread to help optimise LinSolve
    //this->LinSolveThreads(b, x, x0, a, 1 + 6 * a, iter, N);
}

// Projection
// NOTES:
// vx and vy: Arrays representing the x and y components of the velocity field.
// p: Array representing the pressure field.
// div: Array representing the divergence of the velocity field.
void Physics::Project(float vx[], float vy[], float p[], float div[], int iter, int N) {
    // Calculate the divergence of the velocity field
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j, N)] = -0.5f * ( vx[IX(i + 1, j, N)] - vx[IX(i - 1, j, N)] + vy[IX(i, j + 1, N)] - vy[IX(i, j - 1, N)]) / N;
            p[IX(i, j, N)] = 0;
        }
    }

    // Enforce boundary conditions on divergence and pressure
    this->SetBnd(0, div, N);
    this->SetBnd(0, p, N);

    // Solve the Poisson equation to obtain a pressure field
    this->LinSolve(0, p, div, 1, 6, iter, N);
    // Use Thread to help optimise LinSolve
    //this->LinSolveThreads(0, p, div, 1, 6, iter, N);


    // Update the velocity field based on the pressure gradient
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            vx[IX(i, j, N)] -= 0.5f * (p[IX(i + 1, j, N)] - p[IX(i - 1, j, N)]) * N;
            vy[IX(i, j, N)] -= 0.5f * (p[IX(i, j + 1, N)] - p[IX(i, j - 1, N)]) * N;
        }
    }

    // Enforce boundary conditions on the updated velocity field
    this->SetBnd(1, vx, N);
    this->SetBnd(2, vy, N);
}

// Advecting (moving the fluid density field)
// Notes:
// b: Boundary condition flag.
// d: Array representing the fluid density field after advection.
// d0: Array representing the original fluid density field.
// vx and vy: Arrays representing the x and y components of the velocity field.
// dt: deltaTime
void Physics::Advect(int b, float d[], float d0[], float vx[], float vy[], float dt, int N) {
    float i0, i1, j0, j1;

    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);

    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;

    float Nfloat = N;
    float ifloat, jfloat;

    int i, j;

    // Iterate over the fluid grid
    for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
        for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            // Calculate the back-traced position using the velocity field
            tmp1 = dtx * vx[IX(i, j, N)];
            tmp2 = dty * vy[IX(i, j, N)];
            x = ifloat - tmp1;
            y = jfloat - tmp2;

            // Ensure the back-traced position is within bounds
            if (x < 0.5f) 
            {
                x = 0.5f;
            }

            if (x > Nfloat + 0.5f)
            {
                x = Nfloat + 0.5f;
            }

            i0 = ::floorf(x);
            i1 = i0 + 1.0f;

            if (y < 0.5f)
            {
                y = 0.5f;
            }

            if (y > Nfloat + 0.5f)
            {
                y = Nfloat + 0.5f;
            }

            j0 = ::floorf(y);
            j1 = j0 + 1.0f;

            // Interpolate the fluid density from the original positions
            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;

            int i0i = i0;
            int i1i = i1;
            int j0i = j0;
            int j1i = j1;

            d[IX(i, j, N)] =
                s0 * (t0 * d0[IX(i0i, j0i, N)] + t1 * d0[IX(i0i, j1i, N)]) +
                s1 * (t0 * d0[IX(i1i, j0i, N)] + t1 * d0[IX(i1i, j1i, N)]);
        }
    }
    // Enforce boundary conditions on the advected fluid density field
    this->SetBnd(b, d, N);
}