#include "imgui.h"
#include "imgui-SFML.h"

#include <vector>

#include "particles.h"

class Sim {
private:
	sf::RenderWindow win;
	sf::Clock clock;
	Particles particles;
	void Setup();
	float3 ColorConverter(float3 color, bool isConverted);
public:
	static const int numParticles;
	Sim();
	~Sim();
	void Run();
};