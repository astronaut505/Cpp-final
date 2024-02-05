#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <random>
#include <numeric>
#include <cstring>
#include <filesystem>

// Define a structure to represent an option
struct Option {
    std::string type;          // "call" or "put"
    std::string position;      // "long" or "short"
    double strikePrice;
    double underlyingPrice;    // Price of the underlying asset
    double daysToExpiry;       // Time to expiry in days
    double volatility;         // Annualized volatility (as a percentage)
    double tradeVolume;
    double price;              // Calculated price of the option

    Option(const std::string& t, const std::string& pos, double strike, double days, double volume, double underlying = 5000, double vol = 20)
        : type(t), position(pos), strikePrice(strike), underlyingPrice(underlying), daysToExpiry(days), volatility(vol), tradeVolume(volume), price(0) {}
};

class OptionsPortfolio {
private:
    std::vector<double> payoffs; // Store each option's payoff
    double meanPayoff;           // Mean payoff of the portfolio
    double portfolioVariance;    // Variance of the portfolio based on payoffs
    double standardDeviation;    // Standard deviation of the portfolio based on payoffs
    std::vector<Option> optionsList;   // Container for options
    double meanPrice;                  // Mean price of the portfolio
    std::default_random_engine generator;
    std::normal_distribution<double> distribution;

public:
    OptionsPortfolio() : meanPayoff(0), portfolioVariance(0), standardDeviation(0) {}

    double calculatePayoff(const Option& option) {
    double payoff = 0.0;
    if (option.type == "Call") {
        payoff = std::max(option.underlyingPrice - option.strikePrice, 0.0);
    } else if (option.type == "Put") {
        payoff = std::max(option.strikePrice - option.underlyingPrice, 0.0);
    }
    // Adjust for position; if short, payoff is negated
    if (option.position == "short") {
        payoff = -payoff;
    }
    return payoff;
}
    void addOption(Option option) {
    double price = calculateOptionPrice(option);
    payoffs.push_back(price);
    updateStatistics();
}

    // Modified calculateOptionPrice method to include Monte Carlo simulation
    // Black-Scholes Price calculation
    double blackScholesPrice(const Option& option, double riskFreeRate) {
        double S = option.underlyingPrice;
        double K = option.strikePrice;
        double T = option.daysToExpiry / 365.0; // Convert days to years
        double sigma = option.volatility / 100.0;
        double r = riskFreeRate;

        double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
        double d2 = d1 - sigma * sqrt(T);

        if (option.type == "call") {
            return S * N(d1) - K * exp(-r * T) * N(d2);
        } else { // "put"
            return K * exp(-r * T) * N(-d2) - S * N(-d1);
                }
        }

    double calculateOptionPrice(const Option& option) {
        const double riskFreeRate = 0.05; // Example risk-free rate

        // Calculate option price using Black-Scholes formula
        double optionPrice = blackScholesPrice(option, riskFreeRate);

        // For demonstration purposes, if you still want to simulate and calculate a mean payoff for some reason:
        // This part doesn't make practical sense with BSM as it calculates the price directly, but included for alignment with your structure
        const int num_simulations = 10000; // Arbitrary, as we're not actually simulating with BSM
        std::vector<double> simulatedPrices(num_simulations, optionPrice); // Fill the vector with the calculated BSM price

        double meanPrice = std::accumulate(simulatedPrices.begin(), simulatedPrices.end(), 0.0) / simulatedPrices.size();

        // Assuming you wanted to see how the direct BSM price compares when treated as a mean of simulated constant values
        return meanPrice; // This will just return the BSM price as we're not simulating different outcomes
        }

    double simulateAssetPrice(double initialPrice, double volatility, double riskFreeRate, double time) {
        std::normal_distribution<double> distribution(0.0, 1.0);
        double stochasticComponent = distribution(generator);
        return initialPrice * exp((riskFreeRate - 0.5 * volatility * volatility) * time + volatility * sqrt(time) * stochasticComponent);
    }

    void simulateAndPriceOption(Option& option) {
        const double r = 0.05; // Risk-free rate

        // Simulate underlying price
        double simulatedPrice = simulateAssetPrice(option.underlyingPrice, option.volatility / 100.0, r, option.daysToExpiry / 365.0);
        option.underlyingPrice = simulatedPrice; // Update the option with the simulated price

        // Calculate option price using Black-Scholes
        double optionPrice = calculateOptionPrice(option);
        option.price = optionPrice; // Update the option with the calculated price
    }


    static double N(double x) {
        return 0.5 * (1 + std::erf(x / std::sqrt(2.0)));
    }


    void updateStatistics() {
        if (payoffs.empty()) return;

        // Calculate mean payoff
        meanPayoff = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / payoffs.size();

        // Calculate variance
        double varianceSum = std::accumulate(payoffs.begin(), payoffs.end(), 0.0, [this](double acc, double payoff) {
            return acc + std::pow(payoff - meanPayoff, 2);
        });
        portfolioVariance = varianceSum / payoffs.size();

        // Calculate standard deviation
        standardDeviation = std::sqrt(portfolioVariance);
    }

    void printStatistics() const {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Portfolio Payoff Statistics:\n";
        std::cout << "Mean Payoff: " << meanPayoff << "\n";
        std::cout << "Variance: " << portfolioVariance << "\n";
        std::cout << "Standard Deviation: " << standardDeviation << "\n\n";
    }


   void loadOptionsFromFile(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Unable to open file." << std::endl;
        return;
    }

    std::string line;
    int tradeCount = 0;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string type, position, temp;
        double strikePrice, daysToExpiry, tradeVolume;

        std::getline(ss, type, ',');
        std::getline(ss, position, ',');
        std::getline(ss, temp, ','); strikePrice = std::stod(temp);
        std::getline(ss, temp, ','); daysToExpiry = std::stod(temp);
        std::getline(ss, temp); tradeVolume = std::stod(temp);

        // Adjust the type and position if you're using shorthand notation in the file
        type = (type == "C" ? "call" : "put");
        position = (position == "L" ? "long" : "short");

        Option option(type, position, strikePrice, daysToExpiry, tradeVolume);
        addOption(option);

        // Print statistics after adding each option
        std::cout << "After adding option " << ++tradeCount << ":\n";
        printStatistics();
    }
    file.close();
    }

};


int main() {
    OptionsPortfolio portfolio;
    std::filesystem::current_path(std::filesystem::path(__FILE__).parent_path());
    std::filesystem::path currentPath = std::filesystem::current_path();
    std::cout << "Current working directory: " << currentPath << std::endl;

    // Load options from a CSV file (assuming "sample_large.CSV" is case-sensitive)
    portfolio.loadOptionsFromFile("sample_large.CSV");

    // Print statistics after loading the options
    std::cout << "After loading options from file:" << std::endl;
    portfolio.printStatistics();

    return 0;
}



