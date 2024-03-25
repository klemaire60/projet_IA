#include "SFML/Graphics.hpp"
#include "Rectangle.h"
#include "Pixel.h"

int main()
{
	//Création de la fenêtre 
	sf::RenderWindow window(sf::VideoMode(800, 600), "Simulation IA");

	//Création du rectangle
	Rectangle rectangle(700, 500, 50, 50);

	//Création du pixel
	Pixel pixel(20, 20, 400, 300, 100);

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}
		window.clear(/*sf::Color::Black*/);

		rectangle.draw(window);
		pixel.draw(window);

		window.display();
	}
}