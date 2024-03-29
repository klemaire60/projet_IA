#pragma once
#include <SFML/Graphics.hpp>

class Rectangle
{
private:
    sf::RectangleShape rectangleShape;

    sf::Vector2f position;

public:
    Rectangle(int width, int height, float x, float y);
    ~Rectangle();
    void draw(sf::RenderWindow& window) const;

    sf::Vector2f getPosition() const;
    void setPosition(sf::Vector2f newPosition);
};
