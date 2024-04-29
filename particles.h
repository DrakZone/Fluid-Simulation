#include <vector>
#include <thread>
#include<SFML/Graphics.hpp>

#include "physics.h"
#include "constants.h"

class Particles {
private:
	Physics physics;

	int size;

	float dt;
	float diff;
	float visc;

	float px[SIZE * SIZE];
	float py[SIZE * SIZE];

	float x[SIZE * SIZE];
	float y[SIZE * SIZE];

	float previousDensity[SIZE * SIZE];
	float density[SIZE * SIZE];

	float3 Color = { 255.0f, 255.0f, 255.0f };

	void InitArr(float arr[], int size);
	float MapToRange(float value, float minIn, float maxIn, float minOut, float maxOut);
public:
	Particles();
	Particles(float dt, float diff, float visc);
	~Particles();

	void AddDensity(float x, float y, float amount);
	void AddVelocity(float x, float y, float px, float py);
	void Step(float dt);
	void Render(sf::RenderWindow& win);
	void parallelRender(sf::RenderWindow& win);
	void FadeDensity(int size);

	// Getter and setter for viscosity
	float GetViscosity() const { return visc; }
	void SetViscosity(float newViscosity) { visc = newViscosity; }

	float3 GetColor() const { return Color; }
	void SetColor(float3 newcolor) { Color = newcolor; }
};