#include "Rectangle.h"

Rectangle::Rectangle(int width, int height, float x, float y)
{
    setPosition(sf::Vector2f(x, y));
    rectangleShape.setSize(sf::Vector2f(width, height));
    rectangleShape.setPosition(position);
    rectangleShape.setFillColor(sf::Color::Transparent);
    rectangleShape.setOutlineColor(sf::Color::White);
    rectangleShape.setOutlineThickness(5.f);
}

Rectangle::~Rectangle()
{
}

void Rectangle::draw(sf::RenderWindow& window) const
{
    window.draw(rectangleShape);
}

sf::Vector2f Rectangle::getPosition() const
{
    return position;
}

void Rectangle::setPosition(sf::Vector2f newPosition)
{
    position = newPosition;
    rectangleShape.setPosition(position);
}
