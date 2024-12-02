#include "Environment.h"

//Environment::Environment(int width, int height, float x, float y);

void Environment::generateObstacle()
{
    obstacle.setPosition(sf::Vector2f(300, 300));
    obstacle.setSize(sf::Vector2f(10, 200));
    obstacle.setFillColor(sf::Color::Red);
}

Environment::~Environment()
{};

void Environment::draw(sf::RenderWindow& window) const
{
    window.draw(rectangleShape);
    window.draw(obstacle);
}

sf::Vector2f Environment::getPosition() const
{
    return position;
}

void Environment::setPosition(sf::Vector2f newPosition)
{
    position = newPosition;
    rectangleShape.setPosition(position);
}

void Environment::iterate()
{

}

Robot::Robot()
{
}

void Robot::applyAI(sf::RectangleShape& obstacle)
{
}

void Robot::move()
{
}

float Robot::getDistanceFrontObstacle()
{
    return 0.0f;
}

void Robot::setDead(double x)
{
}
