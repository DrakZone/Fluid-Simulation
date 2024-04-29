#include "simulation.h"
#include <iostream>
#include <time.h> // delta time

const int Sim::numParticles = 10;
float visc = 0.000001f;
// Initialize
// Constructor & Destructor
Sim::Sim() : particles(Particles(0.2f, 0, visc)) {
	// Create Window
	this->win.create(sf::VideoMode(SIZE * SCALE, SIZE * SCALE), "2D Fluid Simulation");
	// Seed for DeltaTime
	srand(time_t(static_cast<unsigned>(0)));
}
Sim::~Sim() {}

void Sim::Setup() {}

float3 Sim::ColorConverter(float3 color, bool isConverted)
{
	if (!isConverted)
	{
		color.x /= 255.0f;
		color.y /= 255.0f;
		color.z /= 255.0f;
	}
	else
	{
		color.x *= 255.0f;
		color.y *= 255.0f;
		color.z *= 255.0f;
	}

	return color;
}

// Run / Game loop
void Sim::Run() {
	// Initialize
	this->Setup();
	sf::Vector2i previousMouse = sf::Mouse::getPosition(this->win);
	sf::Vector2i currentMouse = sf::Mouse::getPosition(this->win);
	// Initialize ImGui
	ImGui::SFML::Init(this->win);
	// ImGui Delta Time
	sf::Clock deltaClock;

	// Loop
	while (this->win.isOpen()) {

		// Delta Time
		float deltaTime = clock.getElapsedTime().asSeconds();
		clock.restart().asSeconds();

		// Event handling
		sf::Event event;
		while (this->win.pollEvent(event)) {
			//ImGui Close window
			ImGui::SFML::ProcessEvent(event);
			switch (event.type) {
				//Close window
			case sf::Event::Closed:
				this->win.close();
				break;
			case sf::Event::KeyReleased:
				// up the viscosity of the fluid
				if (event.key.code == sf::Keyboard::Key::Up) {
					visc = this->particles.GetViscosity();
					visc += 0.5f;
					this->particles.SetViscosity(visc);
				}
				// down the viscosity of the fluid
				if (event.key.code == sf::Keyboard::Key::Down) {
					visc = this->particles.GetViscosity();
					visc -= 0.5f;
					this->particles.SetViscosity(visc);
				}
				std::cout << "\rCurrent Viscosity: " << std::fixed << visc << std::flush;
				break;
			default:
				break;
			}
		}

		// ImGui Update
		ImGui::SFML::Update(this->win, deltaClock.restart());

		ImGui::Begin("Settings");
		ImGui::DragFloat("Viscosity", (float*) & visc, 0.1f, -10.0f, 10.0f);
		float3 Color = this->ColorConverter(this->particles.GetColor(), false);
		ImGui::ColorEdit3("Color", (float*)&Color);
		this->particles.SetColor(this->ColorConverter(Color, true));
		ImGui::End();

		// Mouse Input
		if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
			this->particles.AddDensity(currentMouse.y / SCALE, currentMouse.x / SCALE, 200);

		currentMouse = sf::Mouse::getPosition(this->win);

		float amountX = currentMouse.x - previousMouse.x;
		float amountY = currentMouse.y - previousMouse.y;

		this->particles.AddVelocity(currentMouse.y / SCALE, currentMouse.x / SCALE, amountY / 10, amountX / 10);

		previousMouse = currentMouse;

		this->particles.Step(deltaTime);
		this->particles.Render(this->win);
		this->particles.FadeDensity(SIZE * SIZE);

		// ImGui Render
		ImGui::SFML::Render(this->win);

		this->win.display();

		std::cout << "\rTime Elapsed: " << std::fixed << deltaTime << std::flush;
	}

	ImGui::SFML::Shutdown();
}