#include <iostream>
#include "SFML/Graphics.hpp"

int main()
{
	sf::RenderWindow window(sf::VideoMode(800, 600), "Simulation IA");

	sf::RectangleShape rectangle(sf::Vector2f(700.f, 500.f));

	rectangle.setFillColor(sf::Color::Transparent);
	rectangle.setOutlineThickness(10.f);
	rectangle.setOutlineColor(sf::Color::White);


	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		window.clear(sf::Color::Black);

		rectangle.setPosition(50.f, 50.f);

		window.draw(rectangle);

		window.display();
	}
}