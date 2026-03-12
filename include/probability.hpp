#pragma once
#include <cstdint>
#include <vector>
#include <cmath>

namespace math_engine::probability {
    // Performs one bernoulli trial
    // Returns 1 with probability p of success and
    // probability 1 - p of failure
    int bernoulli_trial(double p);

    // Probability mass function
    // Probability of getting k successes in n independent 
    // Bernoulli trials with success probability p
    double binomial_pmf(int64_t k, int64_t n, double p);

    // Computes probability mass function of a geometric variable
    // k is number of trials before first success
    double geometric_pmf(int64_t k, double p);

    // Finds expected value of a random variable given
    // its outcomes and probabilities
    double expected_value(const std::vector<double>& outcomes, const std::vector<double>& probabilities);

    // Random Simulations
    // Simulates a Binomial(n, p) random variable multiple times 
    // with Bernoulli trials.
    // Each simulation represents the number of successes in n independent trials
    std::vector<int> simulate_binomial(int64_t n, double p, int trials);

    // Estimates definite integral of a function over [a, b]
    // using Monte Carlo sampling. Returns running estimate of the integral
    // as the number of samples increases
    std::vector<double> monte_carlo_integration(double (*f)(double), double a, double b, int64_t samples);

}