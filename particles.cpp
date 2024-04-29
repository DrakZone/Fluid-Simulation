#include "particles.h"
#include<iostream>

// The IX function
int IX(int x, int y, int N);

// Constructor & Destructor
Particles::Particles() : physics(Physics()) {}
Particles::~Particles() {}

// Initialize
Particles::Particles(float dt, float diff, float visc) {
	this->size = SIZE;
	this->dt = dt;
	this->diff = diff;
	this->visc = visc;

	this->InitArr(this->px, SIZE * SIZE);
	this->InitArr(this->py, SIZE * SIZE);
	this->InitArr(this->x, SIZE * SIZE);
	this->InitArr(this->y, SIZE * SIZE);
	this->InitArr(this->previousDensity, SIZE * SIZE);
	this->InitArr(this->density, SIZE * SIZE);
}

// Initializing Arrays
void Particles::InitArr(float arr[], int size) {
	for (int i = 0; i < size; i++) {
		arr[i] = 0;
	}
}

// Add Density
void Particles::AddDensity(float x, float y, float amount) {
	this->density[IX(x, y, this->size)] += amount;
}

// Add Velocity
void Particles::AddVelocity(float x, float y, float px, float py) {
	int index = IX(x, y, this->size);

	this->x[index] += px;
	this->y[index] += py;
}

// Step / Update (For Physics)
void Particles::Step(float dt) {
	// Diffusion (X)
	this->physics.Diffuse(1, this->px, this->x, GetViscosity(), dt, 16, this->size);
	// Diffusion (Y)
	this->physics.Diffuse(2, this->py, this->y, GetViscosity(), dt, 16, this->size);

	// Projection Step (x & y)
	this->physics.Project(this->px, this->py, this->x, this->y, 16, this->size);

	// Advection (X)
	this->physics.Advect(1, this->x, this->px, this->px, this->py, dt, this->size);
	// Advection (Y)
	this->physics.Advect(2, this->y, this->py, this->px, this->py, dt, this->size);

	// Projection Step Again(x & y)
	this->physics.Project(this->x, this->y, this->px, this->py, 16, this->size);

	// Diffusion and Advection of Density
	this->physics.Diffuse(0, this->previousDensity, this->density, GetViscosity(), dt, 16, this->size);
	this->physics.Advect(0, this->density, this->previousDensity, this->x, this->y, dt, this->size);
}

// Mapping a value from one range to another
float Particles::MapToRange(float val, float minIn, float maxIn, float minOut, float maxOut) {
	// Normalize the input value to the range [0, 1]
	float x = (val - minIn) / (maxIn - minIn);
	// Map the normalized value to the output range [minOut, maxOut]
	float result = minOut + (maxOut - minOut) * x;
	// Ensure that the result is within the specified output range
	return (result < minOut) ? minOut : (result > maxOut) ? maxOut : result;
}

// Render
void Particles::Render(sf::RenderWindow& win) {
	win.clear();
	for (int i = 0; i < this->size; i++) {
		for (int j = 0; j < this->size; j++) {
			sf::RectangleShape rect;
			rect.setSize(sf::Vector2f(SCALE, SCALE));
			rect.setPosition(j * SCALE, i * SCALE);
			rect.setFillColor(sf::Color(Color.x, Color.y, Color.z, (this->density[IX(i, j, this->size)] > 255) ? 255 : this->density[IX(i, j, this->size)]));

			win.draw(rect);
		}
	}
}

// Thread Render
void renderRectangles(int startRow, int endRow, int size, const std::vector<float>& density, sf::RenderWindow& win, const float3& color) {
	for (int i = startRow; i < endRow; ++i) {
		for (int j = 0; j < size; ++j) {
			sf::RectangleShape rect;
			rect.setSize(sf::Vector2f(SCALE, SCALE));
			rect.setPosition(j * SCALE, i * SCALE);
			int alpha = density[IX(i, j, size)];
			rect.setFillColor(sf::Color(static_cast<sf::Uint8>(color.x * 255),
				static_cast<sf::Uint8>(color.y * 255),
				static_cast<sf::Uint8>(color.z * 255),
				(alpha > 255) ? 255 : alpha));
			win.draw(rect);
		}
	}
}

void Particles::parallelRender(sf::RenderWindow& win) {
    const int numThreads = std::thread::hardware_concurrency();
    const int chunkSize = this->size / numThreads;
    std::vector<std::thread> threads;

	// Convert density to std::vector<float>
	std::vector<float> densityVec(this->density, this->density + this->size * this->size);

    // Create threads and assign work
    for (int i = 0; i < numThreads; ++i) {
        int startRow = i * chunkSize;
        int endRow = (i == numThreads - 1) ? this->size : (i + 1) * chunkSize;
        threads.emplace_back([=, &win] {
            renderRectangles(startRow, endRow, this->size, densityVec, win, this->Color);
        });
    }

    // Wait for all threads to finish
    for (std::thread& thread : threads) {
        thread.join();
    }
}

// Fade Density (remove added density)
void Particles::FadeDensity(int size) {
	for (int i = 0; i < size; i++) {
		// Density Update
		float d = this->density[i];
		density[i] = (d - 0.05f < 0) ? 0 : d - 0.05f;
	}
}