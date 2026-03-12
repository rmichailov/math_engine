#pragma once
#include <cstdint>
#include <vector>
#include <tuple>

namespace math_engine::number_theory {
    // Compute the greatest common divisor of two integers
    int64_t gcd(int64_t x, int64_t y);

    // Compute the LCM of two integers
    int64_t lcm(int64_t x, int64_t y);

    // Returns the gcd of two numbers as well as the coefficients for Bizout's Identity:
    // ax + by = gcd(x, y)
    std::tuple<int64_t, int64_t, int64_t> extended_euclidean(int64_t x, int64_t y); 

    // Check if an integer is prime
    bool is_prime(int64_t num);

    // Returns the prime factors of a number
    std::vector<int64_t> prime_factors(int64_t num);

    //Returns the prime factorization of a number
    std::vector<std::array<int64_t, 2>> prime_factorization(int64_t num);

    // Returns the divisors of a number
    std::vector<int64_t> divisors(int64_t num);

    // Find the number of divisors an integer has
    int64_t count_divisors(int64_t num);

    // Add the divisors of an integer
    int64_t sum_divisors(int64_t num);    

    // Returns the number of positive integers less 
    // than or equal to num that are relatively prime
    // to num.
    int64_t totient(int64_t num);

    // Checks if two integers are coprime (relatively prime)
    bool is_coprime(int64_t x, int64_t y);

    // Modular Arithmetic
    int64_t mod_add(int64_t x, int64_t y, int64_t m);
    int64_t mod_sub(int64_t x, int64_t y, int64_t m);
    int64_t mod_mul(int64_t x, int64_t y, int64_t m);
    int64_t mod_pow(int64_t x, int64_t y, int64_t m);
    int64_t mod_inverse(int64_t x, int64_t m);

    // Combination Formula [modulo m]
    int64_t nCr(int64_t n, int64_t r);
    int64_t nCr(int64_t n, int64_t r, int64_t m);

    // Permutation Formula [modulo m]
    int64_t nPr(int64_t n, int64_t r);
    int64_t nPr(int64_t n, int64_t r, int64_t m);


    // Chinese remainder theorem
    int64_t chinese_remainder_theorem(const std::vector<int64_t>& x, const std::vector<int64_t>& y);
    

    // Returns whether the number is divisible
    // by a square greater than 1
    bool is_divisible_by_square(int64_t num);


    // Function that assigns values to positive integers 
    // based on prime factorization
    // Always (-1, 0, 1)
    // 1 -> num = 1
    // (-1)^k -> if n is the product of k distinct primes
    // 0 if n is divisible by a square > 1
    // (Used to study primes and their distribution)
    int mobius(int64_t num);

    // 
    // Always (-1, 0, 1)    
    int jacobi_symbol(int64_t x, int64_t num);

    // Checks if a number is a perfect square
    bool is_perfect_square(int64_t num);    


}