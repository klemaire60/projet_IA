#include "Pixel.h"

Pixel::Pixel(int width, int height, float x, float y, float speed, Rectangle& targetRect) : targetRectangle(targetRect)
{
    pixelShape.setSize(sf::Vector2f(width, height));
    pixelShape.setPosition(x, y);
    pixelShape.setFillColor(sf::Color::White);

    this->speed = speed;
    position = sf::Vector2f(x, y);
    setDirection(sf::Vector2f(1, 0));
}

Pixel::~Pixel()
{
}

void Pixel::update(float ellapsedTime)
{
    sf::Vector2f targetPosition = targetRectangle.getPosition();
    sf::Vector2f directionToTarget = targetPosition - position;
    float distanceToTarget = std::sqrt(directionToTarget.x * directionToTarget.x + directionToTarget.y * directionToTarget.y);

    if (distanceToTarget > speed * ellapsedTime)
    {
        direction = directionToTarget / distanceToTarget;
        position.x += direction.x * speed * ellapsedTime;
        position.y += direction.y * speed * ellapsedTime;
        pixelShape.setPosition(position);
    }
    else
    {
        position = targetPosition;
        pixelShape.setPosition(targetPosition);
    }
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

sf::Vector2f Pixel::getPosition() const
{
    return position;
}

void Pixel::setDirection(sf::Vector2f newDirection)
{
    direction = newDirection;
}
