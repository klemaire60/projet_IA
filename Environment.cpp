#include "Environment.h"
#include "NEATAlgorithm.h"

Gate::Gate(double x, double yMin, double yMax)
{
	this->x = x;
	this->yMin = yMin;
	this->yMax = yMax;
}

Bird::Bird()
{
	x = 0;
	y = 50;

	yForce = 0;
	isDead = false;
	fitness = 0;
	brain = nullptr;
	hasJump = false;
}

void Bird::applyGravity(double gravity)
{
	yForce -= gravity * 0.017;
}

void Bird::applyAI(Gate nextGate, Gate secondNextGate)
{
	if (brain != nullptr && !isDead)
	{
		double xVectorToGate = nextGate.x - x;
		double yVectorToGateBottom = nextGate.yMin - y;
		double yVectorToGateTop = nextGate.yMax - y;

		double distanceToGateBottom = sqrt((xVectorToGate * xVectorToGate) + (yVectorToGateBottom * yVectorToGateBottom));
		double distanceToGateTop = sqrt((xVectorToGate * xVectorToGate) + (yVectorToGateTop * yVectorToGateTop));

		std::vector<double> inputs;
		inputs.push_back(y);
		inputs.push_back(nextGate.x - x);
		inputs.push_back(secondNextGate.x - x);
		inputs.push_back(nextGate.yMin);
		inputs.push_back(nextGate.yMax);
		inputs.push_back(secondNextGate.yMin);
		inputs.push_back(secondNextGate.yMax);

		std::vector<double> out = NEATAlgorithm::evaluateNetwork(*brain->network, inputs);

		if (out[0] > 0.8)
		{
			jump();
		}
	}
}


void Bird::jump()
{
	yForce = -30;
	hasJump = true;
}

void Bird::applyPhysics(double xIncrementPerIteration)
{
	x += xIncrementPerIteration;
	y += yForce * 0.017;
}

void Bird::setDead(double x)
{
	if (!isDead)
	{
		isDead = true;
		fitness = x + (hasJump ? 10 : 0);
	}
}

void Environment::generateGates()
{
	std::uniform_int_distribution<int> gateIntervalDistribution(30, 50);
	std::uniform_int_distribution<int> gateAltitudeDistribution(35, 65);
	std::uniform_int_distribution<int> gateSizeDistribution(5, 15);

	if (actors.size() > 0)
	{
		double currentX = actors[0]->x;

		while (gates.size() > 0 && gates.front().x + 150 < currentX)
		{
			gates.pop_front();
		}

		while (currentX + 900 > nextGateX)
		{
			double altitude = gateAltitudeDistribution(device);
			Gate gate(nextGateX, altitude, altitude + gateSizeDistribution(device));
			gates.push_back(gate);
			nextGateX += gateIntervalDistribution(device);
		}
	}
}

void Environment::iterate()
{
	generateGates();
	double currentX = actors[0]->x;
	Gate * nextGate = nullptr;
	Gate * secondNextGate = nullptr;
	int nb = 0;
	for (int i = 0; i < gates.size(); i++)
	{
		if (gates[i].x > currentX)
		{
			if (nb == 0)
			{
				nextGate = &gates[i];
			}
			else if (nb == 1)
			{
				secondNextGate = &gates[i];
			}
			nb++;
			if(nb >= 2)
				break;
		}
	}

	for (int i = 0; i < actors.size(); i++)
	{
		Bird * bird = actors[i];
		double xBeforeUpdate = bird->x;
		double yBeforeUpdate = bird->y;

		bird->applyAI(*nextGate, *secondNextGate);
		bird->applyGravity(gravity);
		bird->applyPhysics(xIncrementPerIteration);

		double xAfterUpdate = bird->x;
		double yAfterUpdate = bird->y;

		if (xBeforeUpdate <= nextGate->x && xAfterUpdate >= nextGate->x)
		{
			double meanY = (yBeforeUpdate + yAfterUpdate) / 2.0;
			if (meanY < nextGate->yMin || meanY > nextGate->yMax)
			{
				bird->setDead(nextGate->x);
			}
		}
	}
}
