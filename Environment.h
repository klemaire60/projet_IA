#pragma once
#include <SFML/Graphics.hpp>
#include <deque>
#include <vector>
#include <random>

struct Genome;

struct Gate
{
	double x;
	double yMin;
	double yMax;

	Gate(double x, double yMin, double yMax);
};

struct Robot
{
	double x;
	double y;

	double yForce;
	bool isDead;
	double fitness;
	bool hasJump;

	Genome* brain;

	Robot();
	void applyAI(sf::RectangleShape &obstacle);
	void move();
	float getDistanceFrontObstacle();
	
	void setDead(double x);
};

class Environment
{
private:
    sf::RectangleShape rectangleShape;

	sf::RectangleShape obstacle;

    sf::Vector2f position;

	std::random_device device;
	std::deque<Gate> gates;
	double xIncrementPerIteration; // Iteration = 17 ms
	double gravity;

	std::vector<Robot*> actors;
	void generateObstacle();
	double nextGateX = 50;

public:
	Environment()// int width, int height, float x, float y) // 700, 500, 50, 50
	{
		setPosition(sf::Vector2f(50, 50));
		rectangleShape.setSize(sf::Vector2f(700, 500));
		rectangleShape.setPosition(position);
		rectangleShape.setFillColor(sf::Color::Transparent);
		rectangleShape.setOutlineColor(sf::Color::White);
		rectangleShape.setOutlineThickness(5.f);

		xIncrementPerIteration = 0.5;
		gravity = -90.89;

		generateObstacle();
	};
    ~Environment();
    void draw(sf::RenderWindow& window) const;

    sf::Vector2f getPosition() const;
    void setPosition(sf::Vector2f newPosition);

	double getXOffset()
	{
		if (actors.size() > 0)
			return actors[0]->x;

		return 0;
	}

	std::deque<Gate>& getGates() { return gates; }
	std::vector<Robot*>& getActors() { return actors; }
	void addActor(Robot* actor) { actors.push_back(actor); }

	bool haveActorAlive() {
		bool alive = false;
		for (int i = 0; i < actors.size(); i++)
		{
			if (!actors[i]->isDead)
			{
				alive = true;
				break;
			}
		}
		return alive;
	}

	void iterate();
};