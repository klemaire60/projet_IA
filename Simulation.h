#pragma once

class Environment;

class Simulation
{
	Environment * env;
	int result;	//
	//GameState * state;

	//PlayerInterface * adversary;

public: 
	Simulation(Environment * env);
	~Simulation();
	void playAGame();

	inline int getResult() { return result; }

	void printGameState();
};

