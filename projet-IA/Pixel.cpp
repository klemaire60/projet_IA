#include "Pixel.h"
#include <iostream>

Pixel::Pixel(int width, int height, float x, float y, float speed, Rectangle& targetRect) : targetRectangle(targetRect)
{
    pixelShape.setSize(sf::Vector2f(width, height));
    pixelShape.setPosition(x, y);
    pixelShape.setFillColor(sf::Color::White);

    if (speed <= 0) speed = 1;
    else this->speed = speed;
    position = sf::Vector2f(x, y);
    setDirection(sf::Vector2f(1, 0));
}

Pixel::~Pixel()
{
}

void Pixel::update(float ellapsedTime, Perceptron& perceptron)
{
    sf::Vector2f targetPosition = targetRectangle.getPosition();
    sf::Vector2f directionToTarget = targetPosition - position;

    double input[1] = { std::abs(directionToTarget.x) };
    double predictedSpeedFactor = perceptron.predict(input);

    // Appliquer un facteur de vitesse prédit par le perceptron à la vitesse de base
    double currentSpeed = baseSpeed + predictedSpeedFactor;

    std::cout << baseSpeed << std::endl;
    std::cout << predictedSpeedFactor << std::endl;
    std::cout << currentSpeed << std::endl;

    position.x += currentSpeed * ellapsedTime;
    pixelShape.setPosition(position);
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
