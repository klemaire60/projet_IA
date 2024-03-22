#include "SFML/Graphics.hpp"
#include "Rectangle.h"

int main()
{
	//Création de la fenêtre 
	sf::RenderWindow window(sf::VideoMode(800, 600), "Simulation IA");

	//Création du rectangle
	Rectangle rectangle(700, 500);
	rectangle.setPosition(50.f, 50.f);
	rectangle.setBorderColor(sf::Color::White);
	rectangle.setBorderThickness(5.f);

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}
		window.clear(/*sf::Color::Black*/);

		window.draw(rectangle);

		window.display();
	}
}