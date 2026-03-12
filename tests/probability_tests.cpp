#include <cassert>
#include <cmath>
#include <vector>
#include "probability.hpp"

using namespace math_engine::probability;

int main() {

    // ---- Bernoulli ----
    int result = bernoulli_trial(0.5);
    assert(result == 0 || result == 1);

    // ---- Binomial PMF ----
    double prob = binomial_pmf(2, 4, 0.5);
    // Expected: C(4,2) * (0.5)^4 = 6 * 0.0625 = 0.375
    assert(std::abs(prob - 0.375) < 1e-9);

    // ---- Geometric PMF ----
    double g = geometric_pmf(3, 0.5);
    // (1 - 0.5)^(2) * 0.5 = 0.125
    assert(std::abs(g - 0.125) < 1e-9);

    // ---- Expected Value ----
    std::vector<double> outcomes = {0, 1};
    std::vector<double> probabilities = {0.5, 0.5};
    double ev = expected_value(outcomes, probabilities);
    assert(std::abs(ev - 0.5) < 1e-9);

    // ---- Binomial Simulation ----
    auto sims = simulate_binomial(10, 0.5, 1000);
    assert(sims.size() == 1000);

    // ---- Monte Carlo ----
    auto f = [](double x) { return x * x; };
    auto estimates = monte_carlo_integration(f, 0.0, 1.0, 1000);
    assert(estimates.size() == 1000);

    return 0;
}