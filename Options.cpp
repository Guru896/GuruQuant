#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

double generateRandomNumber() {
    return static_cast<double>(rand()) / RAND_MAX;
}

double simulateStockPrice(double S0, double mu, double sigma, double T, int steps) {
    double dt = T / steps;
    double S = S0;
    for (int i = 0; i < steps; ++i) {
        double dW = sqrt(dt) * generateRandomNumber();
        S *= exp((mu - 0.5 * sigma * sigma) * dt + sigma * dW);
    }
    return S;
}

double europeanCallPayoff(double S, double K) {
    return max(0.0, S - K);
}

double monteCarloOptionPricing(double S0, double K, double T, double mu, double sigma, int steps, int simulations) {
    double sumPayoff = 0.0;
    for (int i = 0; i < simulations; ++i) {
        double S = simulateStockPrice(S0, mu, sigma, T, steps);
        sumPayoff += europeanCallPayoff(S, K);
    }
    return exp(-mu * T) * (sumPayoff / simulations);
}

double normalCDF(double x) {
    return 0.5 * erfc(-x * sqrt(0.5));
}

double blackScholesCallPrice(double S0, double K, double T, double r, double sigma) {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    return S0 * normalCDF(d1) - K * exp(-r * T) * normalCDF(d2);
}

int main() {
    srand(static_cast<unsigned int>(time(0))); // Explicit cast to unsigned int

    double S0 = 100.0; // Initial stock price
    double K = 100.0;  // Strike price
    double T = 1.0;    // Time to maturity (1 year)
    double r = 0.05;   // Risk-free rate
    double sigma = 0.2; // Volatility
    int steps = 100;   // Number of time steps
    int simulations = 10000; // Number of simulations

    double monteCarloPrice = monteCarloOptionPricing(S0, K, T, r, sigma, steps, simulations);
    double blackScholesPrice = blackScholesCallPrice(S0, K, T, r, sigma);

    cout << "European Call Option Price using Monte Carlo: " << monteCarloPrice << endl;
    cout << "European Call Option Price using Black-Scholes: " << blackScholesPrice << endl;

    cin.get(); // Keep the console window open until you press Enter
    return 0;
}
