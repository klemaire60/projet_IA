#pragma once
#include <SFML/Graphics.hpp>

//Position = 

class Rectangle : public sf::RectangleShape
{
private: 
	float borderThickness;

public : 
	Rectangle(float width, float height);

	void setBorderColor(sf::Color color);
	void setBorderThickness(float thickness);
};