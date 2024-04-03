#include "SFML/Graphics.hpp"

#include "Rectangle.h"
#include "Pixel.h"
#include "Perceptron.h"

int main()
{
    // Création de la fenêtre
    sf::RenderWindow window(sf::VideoMode(800, 600), "Simulation IA");

    // Création du rectangle
    Rectangle rectangle(700, 500, 50, 50);

    // Création du pixel
    Pixel pixel(20, 20, 400, 300, 100, rectangle);

    // Création du perceptron
    Perceptron perceptron(0.1, 2);
    
    sf::Clock clock;

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        window.clear(/*sf::Color::Black*/);

        float ellapsedTime = clock.restart().asSeconds();

        // Données d'entraînement
        double trainingData[6][2] = {
            {10, 0.1},
            {20, 0.3},
            {30, 0.6},
            {40, 0.8},
            {50, 1.0},
            {60, 0.5}
        };

        // Entraînement du perceptron
        for (int i = 0; i < 10000; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                perceptron.train(trainingData[j], trainingData[j][1], 0.1);
            }
        }

        pixel.update(ellapsedTime, perceptron);

        rectangle.draw(window);
        pixel.draw(window);

        window.display();
    }
}
