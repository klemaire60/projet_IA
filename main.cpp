#include <SFML/Graphics.hpp>
#include "Environment.h"
#include "NEATAlgorithm.h"
#include <iostream>


int main(int argc, char **argv)
{
	sf::RenderWindow window(sf::VideoMode(800, 600), "BirdAI");
	window.setFramerateLimit(60);

	NEATAlgorithm algorithm(7, 1);
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
				if (algorithm.getNeuralNetwork(0)->fitness >= 10000)
					break;
				printf("------------------------------------------\n\n\n");
			}
		}
		else if (choice == 2) 
		{
			Environment env;
			Robot* robot = new Robot();
			env.addActor(robot);
			//while (env.haveActorAlive())
			{
				sf::Event event;
				while (window.pollEvent(event))
				{
					// évènement "fermeture demandée" : on ferme la fenêtre
					if (event.type == sf::Event::Closed)
						window.close();
				}
				env.iterate();

				window.clear();

				//std::vector<Robot*>& actors = env.getActors();
				/*
				for (int i = 0; i < actors.size(); i++)
				{
					Robot* actor = actors[i];
					sf::RectangleShape rect;
					rect.setPosition(actor->x, actor->y);
					rect.setSize(sf::Vector2f(5, 5));
					rect.setFillColor(actor->isDead ? sf::Color::Red : sf::Color::Green);
					window.draw(rect);
				}
				*/
				sf::RectangleShape rectangleShape;
				rectangleShape.setSize(sf::Vector2f(500, 500));
				rectangleShape.setPosition(50, 50);
				rectangleShape.setFillColor(sf::Color::Transparent);
				rectangleShape.setOutlineColor(sf::Color::White);
				rectangleShape.setOutlineThickness(5.f);

				sf::RectangleShape rect;
				rect.setPosition(300, 300);
				rect.setSize(sf::Vector2f(15, 15));
				rect.setFillColor(sf::Color::Green);


				window.draw(rect);
				window.draw(rectangleShape);
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