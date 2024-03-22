#include "Rectangle.h"

Rectangle::Rectangle(float width, float height) : sf::RectangleShape(sf::Vector2f(width,height))
{
	setFillColor(sf::Color::Transparent);
	setOutlineThickness(10.f);
	setOutlineColor(sf::Color::White);
	borderThickness = 2.f;
}

void Rectangle::setBorderColor(sf::Color color)
{
	setOutlineColor(color);
}

void Rectangle::setBorderThickness(float thickness)
{
	borderThickness = thickness;
	setOutlineThickness(thickness);
}