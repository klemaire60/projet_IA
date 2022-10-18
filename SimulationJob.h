#pragma once

#include "Simulation.h"

class SimulationJob
{
	bool done;
	Simulation * simulation;

public:
	SimulationJob(Simulation * simulation);
	~SimulationJob();
	void doJob();
	inline Simulation * getSimulation() { return simulation; }
	bool isDone() { return done; }
};

