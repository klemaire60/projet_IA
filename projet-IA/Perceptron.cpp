#include <cmath>
#include <cstdlib>
#include <ctime>

#include "Perceptron.h"

Perceptron::Perceptron(double learningRate, int numInputs)
{
    this->learningRate = learningRate;

    // Initialize weights randomly
    srand(time(0));
    for (int i = 0; i < numInputs + 1; i++)
    {
        weights[i] = (double)rand() / RAND_MAX * 2 - 1;
    }
}

void Perceptron::train(double inputs[], double target, double learningRate)
{
    double prediction = predict(inputs);

    double error = target - prediction;

    for (int i = 0; i < 2; i++)
    {
        weights[i] += learningRate * error * inputs[i];
    }

    weights[2] += learningRate * error; // Update bias
}

double Perceptron::predict(double inputs[])
{
    double sum = 0;

    for (int i = 0; i < 2; i++)
    {
        sum += weights[i] * inputs[i];
    }

    sum += weights[2]; // Add bias

    return 1.0 / (1.0 + exp(-sum));
}
