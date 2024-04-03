#ifndef PERCEPTRON_H
#define PERCEPTRON_H

class Perceptron
{
public:
    Perceptron(double learningRate, int numInputs);
    void train(double inputs[], double target, double learningRate);
    double predict(double inputs[]);

private:
    double weights[3]; // 2 inputs + 1 bias
    double learningRate;
};

#endif