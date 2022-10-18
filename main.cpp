#include <SFML/Graphics.hpp>
#include "Environment.h"
#include "NEATAlgorithm.h"
#include <iostream>


int main(int argc, char **argv)
{
	sf::RenderWindow window(sf::VideoMode(800, 600), "BirdAI");
	window.setFramerateLimit(60);

	NEATAlgorithm algorithm(5, 1);
	algorithm.evaluatePopulation();

	int quit = 0;

	while (quit != 1)
	{
		int choice = -1;

		while (choice < 1 || choice > 3)
		{
			std::cout << "1 - Learn" << std::endl;
			std::cout << "2 - Play best NN" << std::endl;
			std::cout << "3 - Quit" << std::endl;
			std::cin >> choice;
		}

		if (choice == 1)
		{
			int nbGen;
			std::cout << "Nb gen : ";
			std::cin >> nbGen;

			for (int i = 0; i < nbGen; i++)
			{
				Environment env;
				algorithm.newGeneration();
				printf("------------------------------------------\n");
				printf("|              Generation %4d           |\n", algorithm.getGenerationNumber());
				printf("------------------------------------------\n");
				algorithm.evaluatePopulation();
				printf("------------------------------------------\n\n\n");
			}
		}
		else if (choice == 2)
		{
			Environment env;
			Bird * bird = new Bird();
			bird->brain = algorithm.getNeuralNetwork(0);
			env.addActor(bird);
			// on fait tourner le programme jusqu'à ce que la fenêtre soit fermée
			while (env.haveActorAlive())
			{
				// on inspecte tous les évènements de la fenêtre qui ont été émis depuis la précédente itération
				sf::Event event;
				while (window.pollEvent(event))
				{
					// évènement "fermeture demandée" : on ferme la fenêtre
					if (event.type == sf::Event::Closed)
						window.close();
					if (event.type == sf::Event::KeyPressed)
					{
						//bird->jump();
					}
				}

				env.iterate();

				window.clear();

				std::deque<Gate> & gates = env.getGates();
				double xOffset = env.getXOffset() - 100;

				for (int i = 0; i < gates.size(); i++)
				{
					Gate * gate = &gates[i];
					sf::RectangleShape rect;
					rect.setPosition(gate->x - xOffset, gate->yMin);
					rect.setSize(sf::Vector2f(5, gate->yMax - gate->yMin));
					window.draw(rect);
				}


				std::vector<Bird*> & actors = env.getActors();
				for (int i = 0; i < actors.size(); i++)
				{
					Bird * actor = actors[i];
					sf::RectangleShape rect;
					rect.setPosition(actor->x - xOffset, actor->y);
					rect.setSize(sf::Vector2f(5, 5));
					rect.setFillColor(actor->isDead ? sf::Color::Red : sf::Color::Green);
					window.draw(rect);
				}

				window.display();
			}
		}
		else if (choice == 3)
		{
			quit = 1;
		}
	}

	return 0;
}