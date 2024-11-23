#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <sstream>

using namespace std;

vector<double> readCSV(const string& filename) {
    vector<double> data;
    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string value;
        while (getline(ss, value, ',')) {
            data.push_back(stod(value));
        }
    }
    return data;
}

vector<double> calculateReturns(const vector<double>& prices) {
    vector<double> returns;
    for (size_t i = 1; i < prices.size(); ++i) {
        double ret = (prices[i] - prices[i - 1]) / prices[i - 1];
        returns.push_back(ret);
    }
    return returns;
}

double calculateMean(const vector<double>& returns) {
    double sum = accumulate(returns.begin(), returns.end(), 0.0);
    return sum / returns.size();
}

double calculateCovariance(const vector<double>& returns1, const vector<double>& returns2) {
    double mean1 = calculateMean(returns1);
    double mean2 = calculateMean(returns2);
    double covariance = 0.0;
    for (size_t i = 0; i < returns1.size(); ++i) {
        covariance += (returns1[i] - mean1) * (returns2[i] - mean2);
    }
    return covariance / (returns1.size() - 1);
}

vector<vector<double>> calculateCovarianceMatrix(const vector<vector<double>>& returns) {
    size_t n = returns.size();
    vector<vector<double>> covarianceMatrix(n, vector<double>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            covarianceMatrix[i][j] = calculateCovariance(returns[i], returns[j]);
        }
    }
    return covarianceMatrix;
}

double calculatePortfolioVariance(const vector<double>& weights, const vector<vector<double>>& covarianceMatrix) {
    double variance = 0.0;
    for (size_t i = 0; i < weights.size(); ++i) {
        for (size_t j = 0; j < weights.size(); ++j) {
            variance += weights[i] * weights[j] * covarianceMatrix[i][j];
        }
    }
    return variance;
}

double calculatePortfolioReturn(const vector<double>& weights, const vector<double>& means) {
    double portfolioReturn = 0.0;
    for (size_t i = 0; i < weights.size(); ++i) {
        portfolioReturn += weights[i] * means[i];
    }
    return portfolioReturn;
}

void optimizePortfolio(const vector<vector<double>>& covarianceMatrix, const vector<double>& means, vector<double>& weights) {
    size_t numAssets = means.size();
    weights = vector<double>(numAssets, 1.0 / numAssets);
}

int main() {
    vector<double> asset1 = readCSV("AAPL.csv");
    vector<double> asset2 = readCSV("MSFT.csv");
    vector<double> asset3 = readCSV("AMZN.csv");

    vector<double> returns1 = calculateReturns(asset1);
    vector<double> returns2 = calculateReturns(asset2);
    vector<double> returns3 = calculateReturns(asset3);

    vector<vector<double>> returns = { returns1, returns2, returns3 };

    vector<double> means = {
        calculateMean(returns1),
        calculateMean(returns2),
        calculateMean(returns3)
    };

    vector<vector<double>> covarianceMatrix = calculateCovarianceMatrix(returns);

    cout << "Expected Returns:" << endl;
    for (const auto& mean : means) {
        cout << mean << " ";
    }
    cout << endl << "\nCovariance Matrix:" << endl;
    for (const auto& row : covarianceMatrix) {
        for (const auto& value : row) {
            cout << value << " ";
        }
        cout << endl;
    }

    vector<double> weights;
    optimizePortfolio(covarianceMatrix, means, weights);

    double optimizedReturn = calculatePortfolioReturn(weights, means);
    double optimizedVariance = calculatePortfolioVariance(weights, covarianceMatrix);

    cout << "\nOptimized Weights:" << endl;
    for (const auto& weight : weights) {
        cout << weight << " ";
    }
    cout << endl << "Expected Return: " << optimizedReturn << endl;
    cout << "Portfolio Variance: " << optimizedVariance << endl;

    cin.get();
    return 0;
}
