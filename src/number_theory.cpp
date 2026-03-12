#include "number_theory.hpp"
#include <stdexcept>
#include <array>
#include <cmath>
#include <iostream>

namespace math_engine::number_theory {
    

    int64_t gcd(int64_t x, int64_t y) {
        while (y != 0) {
            int64_t temp = y;
            y = x % y;
            x = temp;
        }

        return std::abs(x);
    }

    int64_t lcm(int64_t x, int64_t y) {
        if (x == 0 || y == 0) {
            return 0;
        }

        return std::abs(x / gcd(x, y) * y);
    }

    std::tuple<int64_t, int64_t, int64_t> extended_euclidean(int64_t x, int64_t y) {
        if (y == 0) {
            return {x, 1, 0};
        }

        auto [gcd, x1, y1] = extended_euclidean(y, x % y);

        int64_t a = y1;
        int64_t b = x1 - (x / y) * y1;

        return {gcd, a, b};
    }


    bool is_prime(int64_t num) {
        if (num <= 1) {
            return false;
        }

        if (num == 2 || num == 3) {
            return true;
        }

        if ((num & 1) == 0 || num % 3 == 0) {
            return false;
        }

        for (int64_t i = 5; i * i <= num; i += 2) {
            if (num % i == 0) {
                return false;
            }
        }

        return true;
    }

    std::vector<int64_t> prime_factors(int64_t num) {
        std::vector<int64_t> pFactors;

        for (int64_t i = 2; i * i <= num; i++) {
            while (num % i == 0) {
                pFactors.push_back(i);
                num /= i;
            }
        }

        if (pFactors.size() == 0 || num > 1) {
            pFactors.push_back(num);
        }      
        return pFactors;
    }

    std::vector<std::array<int64_t, 2>> prime_factorization(int64_t num) {
        std::vector<std::array<int64_t, 2>> pFactorization;
        for (int64_t i = 2; i * i <= num; i++) {
            int64_t count = 0;
            while (num % i == 0) {
                count++;
                num /= i;
            }
            if (count != 0) {
                pFactorization.push_back({i, count});
            }
            
        }

        if (num > 1) {
            pFactorization.push_back({num, 1});
        }

        return pFactorization;
    }

    std::vector<int64_t> divisors(int64_t num) {
        std::vector<int64_t> d;

        for (int64_t i = 1; i <= num; i++) {
            if (num % i == 0) {
                d.push_back(i);
            }
        }
        return d;
    }



    int64_t count_divisors(int64_t num) {
        std::vector<std::array<int64_t, 2>> pFactorization = prime_factorization(num);
        int64_t count = 1;

        for (const auto& [p, k] : pFactorization) {
            count *= (k + 1);
        }        

        return count;
    }

    
    int64_t sum_divisors(int64_t num) {
        std::vector<int64_t> d = divisors(num);

        int64_t count = 0;

        for (const auto& n : d) {
            count += n;
        }  
        return count;

    }      

    int64_t totient(int64_t num) {
        if (num <= 0) {
            throw std::invalid_argument("Totient function is undefined for num <= 0");
        }
        
        std::vector<std::array<int64_t, 2>> pFactorization = prime_factorization(num);

        int64_t res = 1;

        for (const auto& [p, k] : pFactorization) {
            res *= (std::pow(p, k-1) * (p - 1));
        }
        return res;
    }

    bool is_coprime(int64_t x, int64_t y) {
        return gcd(x, y) == 1;
    }

    int64_t mod_add(int64_t x, int64_t y, int64_t m) {
        return (x % m + y % m + m) % m;
    }

    int64_t mod_sub(int64_t x, int64_t y, int64_t m) {
        return (x % m - y % m + m) % m;
    }

    int64_t mod_mul(int64_t x, int64_t y, int64_t m) {
        return ((x % m) * (y % m) + m) % m;
    }

    int64_t mod_pow(int64_t x, int64_t y, int64_t m) {
        int64_t res = 1;
        x = x % m;
        while (y > 0) {
            if (y & 1) {
                res = ((res * x) % m);
            }
            y = y >> 1;
            x = (x * x) % m;
        }
        return res;
    }


    int64_t mod_inverse(int64_t x, int64_t m) {
        auto [g, a, b] = extended_euclidean(x, m);
        if (g != 1) {
            throw std::invalid_argument("Modular Inverse doesn't exist.");
        }
        return (a % m + m) % m;
    }

    int64_t nCr(int64_t n, int64_t r) {
        if (r < 0 || r > n) {
            return 0;
        }
        int64_t ans = 1;
        int64_t min;
        if (r < n - r) {
            min = r;
        }
        else {
            min = n - r;
        }

        for (int64_t i = 1; i <= min; i++) {
            ans *= (n - i + 1)/i;
        }
        return ans;
    }

    int64_t nCr(int64_t n, int64_t r, int64_t m) {
        if (r < 0 || r > n) {
            return 0;
        }
        int64_t ans = 1;
        int64_t min;
        if (r < n - r) {
            min = r;
        }
        else {
            min = n - r;
        }

        for (int64_t i = 1; i <= min; i++) {
            ans = (ans * (n - i + 1)) % m;
            ans = (ans * mod_inverse(i, m)) % m;
        }
        return ans;
    }

    int64_t nPr(int64_t n, int64_t r) {
        if (r < 0 || r > n) {
            return 0;
        }        
        int64_t ans = 1;
        for (int64_t i = n; i > n - r; i--) {
            ans *= i;
        }

        return ans;
    }

    int64_t nPr(int64_t n, int64_t r, int64_t m) {
        if (r < 0 || r > n) {
            return 0;
        }        
        int64_t ans = 1;
        for (int64_t i = n; i > n - r; i--) {
            ans *= (i) % m;
            ans %= m;
        }

        return ans;
    }

    int64_t chinese_remainder_theorem(const std::vector<int64_t>& x, const std::vector<int64_t>& y) {
        int64_t prod = 1;
        for (auto v : y) prod *= v;

        int64_t res = 0;
        for (size_t i = 0; i < x.size(); i++) {
            int64_t pp = prod / y[i];
            res += x[i] * mod_inverse(pp, y[i]) * pp;
        }
        return (res % prod + prod) % prod;
    }
       

    bool is_divisible_by_square(int64_t num) {
        // int64_t diff = 3;
        // for (int64_t i = 4; i * i <= num; i += diff) {
        //     if (num % i == 0) {
        //         return true;
        //     }
        //     diff += 2;
        // }
        // return false;

        // std::vector<std::array<int64_t, 2>> pFactorization = prime_factorization(num);

        // for (const auto& [p, k] : pFactorization) {
        //     if (k % 2 == 0) {
        //         return true;
        //     }
        // }
        // return false;

        if (num < 2) {
            return false;
        }

        for (int64_t i = 2; i * i <= num; i++) {
            if (num % (i * i) == 0) {
                return true;
            } 
        }
        return false;

    }

    int mobius(int64_t num) {
        if (num == 1) {
            return 1;
        }

        if (is_divisible_by_square(num)) {
            return 0;
        }

        std::vector<int64_t> pFactors = prime_factors(num);

        if (pFactors.size() % 2 == 1) {
            return -1;
        }
        else {
            return 1;
        }
    }

    int jacobi_symbol(int64_t x, int64_t num) {
        if (num <= 0 || num % 2 == 0) {
            throw std::invalid_argument("Num must be odd and positive");
        } 
        x = x % num;
        int result = 1;

        while (x != 0) {
            while (x % 2 == 0) {
                x /= 2;
                int r = num % 8;
                if (r == 3 || r == 5) {
                    result = -result;
                }
            }
            std::swap(x, num);
            if (x % 4 == 3 && num % 4 == 3) {
                result = -result;
            }
            x = x % num;
        }
        return (num == 1) ? result : 0;
    }  

    // Checks if a number is a perfect square
    bool is_perfect_square(int64_t num) {
        if (num < 0) {
            return false;
        }

        for (int64_t i = 1; i * i <= num; i++) {
            if (num == (i * i)) {
                return true;
            } 
        }
        return false;        
    } 
}