#include "Pixel.h"

Pixel::Pixel(int width, int height, float x, float y, float speed)
{
    setPosition(sf::Vector2f(x, y));
    pixelShape.setSize(sf::Vector2f(width, height));
    pixelShape.setPosition(position);
    pixelShape.setFillColor(sf::Color::White);

    setSpeed(speed);
    setDirection(sf::Vector2f(1, 0));
}

Pixel::~Pixel()
{
}

void Pixel::move(float ellapsedTime)
{
    oldPosition = position;
    position.x += direction.x * speed * ellapsedTime;
    position.y += direction.y * speed * ellapsedTime;
}

void Pixel::draw(sf::RenderWindow& window) const
{
    window.draw(pixelShape);
}

float Pixel::getSpeed()
{
    return speed;
}

void Pixel::setSpeed(float newSpeed)
{
    speed = newSpeed;
}

void Pixel::setPosition(sf::Vector2f newPosition)
{
    position = newPosition;
}

void Pixel::setDirection(sf::Vector2f newDirection)
{
    direction = newDirection;
}
