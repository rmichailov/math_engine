#include "probability.hpp"
#include <stdexcept>
#include <array>
#include <cmath>
#include <random>

namespace math_engine::probability {    
    int bernoulli_trial(double p) {
        if (p < 0 || p > 1) {
            throw std::invalid_argument("Probability p must be between 0 and 1");
        }

        static std::random_device rand;
        static std::mt19937 generate(rand());
        std::bernoulli_distribution dist(p);

        return dist(generate);
    }
    
    double binomial_pmf(int64_t k, int64_t n, double p) {
        if (p < 0.0 || p > 1.0 || k < 0 || n < 0 || k > n) {
            return 0.0;
        }

        double coeff = 1.0;
        for (int64_t i = 1; i <= k; i++) {
            coeff *= (n - i + 1) / static_cast<double>(i);
        }

        return coeff * std::pow(p, k) * std::pow(1.0 - p, n - k);

    }

    double geometric_pmf(int64_t k, double p) {
        if (p < 0.0 || p > 1.0 || k < 1) {
            return 0.0;
        }

        return std::pow(1.0 - p, k-1) * p;
    }

    double expected_value(const std::vector<double>& outcomes, const std::vector<double>& probabilities) {
        if (outcomes.size() != probabilities.size()) {
            throw std::invalid_argument("Outcomes and probabilities must have same size.");
        }

        double ev = 0.0;
        for (size_t i = 0; i < outcomes.size(); i++) {
            ev += (outcomes[i] * probabilities[i]);
        }

        return ev;
    }

    // Random Simulations
    std::vector<int> simulate_binomial(int64_t n, double p, int trials) {
        if (n < 0 || trials <= 0 || p < 0.0 || p > 1.0) {
            throw std::invalid_argument("Invalid parameters.");
        }

        std::vector<int> res;
        res.reserve(trials);

        for (int64_t i = 0; i < trials; i++) {
            int succ = 0;
            for (int64_t j = 0; j < n; j++) {
                succ += bernoulli_trial(p);
            }
            res.push_back(succ);
        }
        return res;

    }

    std::vector<double> monte_carlo_integration(double (*f)(double), double a, double b, int64_t samples) {
        if (samples <= 0 || a >= b) {
            throw std::invalid_argument("Invalid parameters.");
        }

        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(a, b);

        std::vector<double> est;
        est.reserve(samples);

        double sum = 0.0;
        for (int64_t i = 1; i <= samples; i++) {
            double x = dist(gen);
            sum += f(x);

            double estimate = (b - a) * (sum / i);
            est.push_back(estimate);
        }
        return est;
    }



}