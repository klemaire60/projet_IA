#pragma once
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

struct Bird
{
	double x;
	double y;

	double yForce;
	bool isDead;
	double fitness;
	bool hasJump;

	Genome * brain;

	Bird();
	void applyGravity(double gravity);
	void applyAI(Gate nextGate, Gate secondNextGate);
	void jump();
	void applyPhysics(double xIncrementPerIteration);
	void setDead(double x);
};

class Environment
{
	std::random_device device;
	std::deque<Gate> gates;
	double xIncrementPerIteration; // Iteration = 17 ms
	double gravity;

	std::vector<Bird*> actors;
	void generateGates();
	double nextGateX = 50;

public:
	Environment() {
		xIncrementPerIteration = 0.5;
		gravity = -90.89;

		generateGates();
	}

	double getXOffset()
	{
		if (actors.size() > 0)
			return actors[0]->x;

		return 0;
	}

	std::deque<Gate> & getGates() { return gates; }
	std::vector<Bird*> & getActors() { return actors; }
	void addActor(Bird* actor) { actors.push_back(actor); }

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


