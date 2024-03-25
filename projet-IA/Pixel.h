#pragma once
#include <SFML/Graphics.hpp>

class Pixel
{
private:
    sf::RectangleShape pixelShape;

    sf::Vector2f oldPosition;
    sf::Vector2f position;
    float speed;
    sf::Vector2f direction;

public:
    Pixel(int width, int height, float x, float y, float speed);
    ~Pixel();
    void move(float ellapsedTime);
    void draw(sf::RenderWindow& window) const;
    float getSpeed();
    void setSpeed(float newSpeed);

    void setPosition(sf::Vector2f newPosition);
    void setDirection(sf::Vector2f newDirection);
};