#pragma once
#include <SFML/Graphics.hpp>

#include "Rectangle.h"
#include "Perceptron.h"

class Pixel
{
private:
    sf::RectangleShape pixelShape;

    sf::Vector2f position;
    float speed;
    float baseSpeed = 1;
    sf::Vector2f direction;

    Rectangle& targetRectangle;

public:
    Pixel(int width, int height, float x, float y, float speed, Rectangle& targetRect);
    ~Pixel();
    void update(float ellapsedTime, Perceptron& perceptron);
    void draw(sf::RenderWindow& window) const;
    float getSpeed();
    void setSpeed(float newSpeed);

    sf::Vector2f getPosition() const;
    void setDirection(sf::Vector2f newDirection);
};
