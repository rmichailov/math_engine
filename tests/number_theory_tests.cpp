#include <cassert>
#include <cmath>
#include <vector>
#include <tuple>
#include <array>
#include <iostream>
#include "number_theory.hpp"

using namespace math_engine::number_theory;

int main() {

    // -----------------------
    // GCD & LCM
    // -----------------------

    assert(gcd(48, 18) == 6);
    assert(gcd(0, 5) == 5);
    assert(gcd(5, 0) == 5);

    assert(lcm(12, 18) == 36);
    assert(lcm(7, 5) == 35);

    // gcd * lcm identity
    int64_t a = 21, b = 6;
    assert(gcd(a, b) * lcm(a, b) == std::abs(a * b));

    // -----------------------
    // Extended Euclidean
    // -----------------------

    auto [g, x, y] = extended_euclidean(30, 12);
    assert(g == 6);
    assert(30 * x + 12 * y == g);

    // -----------------------
    // Prime Testing
    // -----------------------

    assert(is_prime(2));
    assert(is_prime(13));
    assert(!is_prime(1));
    assert(!is_prime(100));

    // -----------------------
    // Prime Factors
    // -----------------------

    auto pf = prime_factors(60); // 2,2,3,5
    std::vector<int64_t> expected_pf = {2,2,3,5};
    assert(pf == expected_pf);

    auto pff = prime_factorization(60);
    // 60 = 2^2 * 3^1 * 5^1
    assert(pff[0][0] == 2 && pff[0][1] == 2);
    assert(pff[1][0] == 3 && pff[1][1] == 1);
    assert(pff[2][0] == 5 && pff[2][1] == 1);

    // -----------------------
    // Divisors
    // -----------------------

    auto divs = divisors(12);
    std::vector<int64_t> expected_divs = {1,2,3,4,6,12};
    assert(divs == expected_divs);

    assert(count_divisors(12) == 6);
    assert(sum_divisors(12) == 28); // 1+2+3+4+6+12

    // -----------------------
    // Totient
    // -----------------------

    assert(totient(9) == 6);   // 1,2,4,5,7,8
    assert(totient(10) == 4);  // 1,3,7,9

    // Euler's theorem check:
    // if gcd(a,n)=1, then a^phi(n) ≡ 1 (mod n)
    int64_t n = 10;
    int64_t phi = totient(n);
    assert(mod_pow(3, phi, n) == 1);

    // -----------------------
    // Coprime
    // -----------------------

    assert(is_coprime(14, 15));
    assert(!is_coprime(14, 21));

    // -----------------------
    // Modular Arithmetic
    // -----------------------

    assert(mod_add(7, 8, 5) == 0); // 15 mod 5
    assert(mod_sub(3, 5, 7) == 5); // -2 mod 7
    assert(mod_mul(4, 6, 5) == 4); // 24 mod 5

    assert(mod_pow(2, 10, 1000) == 24); // 1024 mod 1000

    // modular inverse
    assert(mod_inverse(3, 11) == 4); // 3*4 = 12 ≡ 1 mod 11

    // -----------------------
    // Combinations & Permutations
    // -----------------------

    assert(nCr(5, 2) == 10);
    assert(nPr(5, 2) == 20);

    assert(nCr(5, 2, 7) == 3); // 10 mod 7
    assert(nPr(5, 2, 7) == 6); // 20 mod 7

    // -----------------------
    // Chinese Remainder Theorem
    // Solve:
    // x ≡ 2 (mod 3)
    // x ≡ 3 (mod 5)
    // solution = 8 mod 15
    // -----------------------

    std::vector<int64_t> remainders = {2, 3};
    std::vector<int64_t> moduli = {3, 5};
    assert(chinese_remainder_theorem(remainders, moduli) == 8);

    // -----------------------
    // Square Divisibility
    // -----------------------

    assert(is_divisible_by_square(12));  // divisible by 4
    assert(!is_divisible_by_square(30)); // squarefree

    // -----------------------
    // Mobius Function
    // -----------------------

    assert(mobius(1) == 1);
    assert(mobius(6) == 1);   // 2*3
    assert(mobius(30) == -1); // 3 primes
    assert(mobius(12) == 0);  // divisible by 4

    // -----------------------
    // Jacobi Symbol
    // -----------------------

    assert(jacobi_symbol(1, 7) == 1);
    assert(jacobi_symbol(2, 7) == 1);

    // -----------------------
    // Perfect Square
    // -----------------------

    assert(is_perfect_square(49));
    assert(!is_perfect_square(50));

    return 0;
}