#include "SimulationJob.h"

SimulationJob::SimulationJob(Simulation * simulation)
{
	this->simulation = simulation;
	done = false;
}

SimulationJob::~SimulationJob()
{
}

void SimulationJob::doJob()
{
	simulation->playAGame();
	done = true;
}
