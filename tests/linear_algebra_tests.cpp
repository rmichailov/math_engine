#include <cassert>
#include <cmath>
#include <vector>

#include "linear_algebra.hpp"

using namespace math_engine::linear_algebra;

constexpr double EPS = 1e-9;

bool approx(double a, double b) {
    return std::abs(a - b) < EPS;
}

bool matrix_equal(const Matrix& A, const Matrix& B) {
    if (A.size() != B.size()) return false;
    if (A[0].size() != B[0].size()) return false;

    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < A[0].size(); ++j)
            if (!approx(A[i][j], B[i][j]))
                return false;

    return true;
}

int main() {

    // -----------------------
    // Vector Operations
    // -----------------------

    Vector v = {1, 2, 3};
    Vector w = {4, 5, 6};

    auto sum = vector_add(v, w);
    assert(sum == Vector({5, 7, 9}));

    auto diff = vector_sub(w, v);
    assert(diff == Vector({3, 3, 3}));

    auto scaled = vector_scale(v, 2.0);
    assert(scaled == Vector({2, 4, 6}));

    assert(approx(vector_dot(v, w), 32.0));

    assert(approx(vector_norm({3, 4}), 5.0));

    auto normed = normalized_vector({3, 4});
    assert(approx(normed[0], 0.6));
    assert(approx(normed[1], 0.8));

    // -----------------------
    // Matrix Addition/Subtraction
    // -----------------------
    Matrix A = {
        {1, 2},
        {3, 4}
    };

    Matrix B = {
        {5, 6},
        {7, 8}
    };

    Matrix C = matrix_add(A, B);
    Matrix expected_add = {
        {6, 8},
        {10, 12}
    };
    assert(matrix_equal(C, expected_add));

    Matrix D = matrix_sub(B, A);
    Matrix expected_sub = {
        {4, 4},
        {4, 4}
    };
    assert(matrix_equal(D, expected_sub));
    
    // -----------------------
    // Matrix Scaling
    // -----------------------

    Matrix scaled_A = matrix_scale(A, 2.0);
    Matrix expected_scaled = {
        {2, 4},
        {6, 8}
    };
    assert(matrix_equal(scaled_A, expected_scaled));
    
    // -----------------------
    // Matrix Multiplication
    // -----------------------

    Matrix M1 = {
        {1, 2},
        {3, 4}
    };

    Matrix M2 = {
        {2, 0},
        {1, 2}
    };

    Matrix product = matrix_mul(M1, M2);
    Matrix expected_product = {
        {4, 4},
        {10, 8}
    };

    assert(matrix_equal(product, expected_product));
    
    // -----------------------
    // Matrix-Vector Multiplication
    // -----------------------

    Vector vec = {1, 1};
    Vector mv = matrix_vector_mult(M1, vec);
    assert(mv == Vector({3, 7}));

    // -----------------------
    // Determinant
    // -----------------------

    assert(approx(matrix_determinant_recursion(A), -2.0));

    Matrix A3 = {
        {1, 2, 3},
        {0, 4, 5},
        {1, 0, 6}
    };

    assert(approx(matrix_determinant_recursion(A3), 22.0));

    // -----------------------
    // Transpose
    // -----------------------

    Matrix At = matrix_transpose(A);
    Matrix expected_transpose = {
        {1, 3},
        {2, 4}
    };
    assert(matrix_equal(At, expected_transpose));
    
    // -----------------------
    // Inverse
    // -----------------------

    Matrix invA = matrix_inverse(A);

    Matrix identity = matrix_mul(A, invA);

    Matrix expected_identity = {
        {1, 0},
        {0, 1}
    };

    assert(matrix_equal(identity, expected_identity));
    
    // -----------------------
    // RREF
    // -----------------------

    Matrix R = {
        {1, 2, 1},
        {2, 4, 2},
        {3, 6, 3}
    };

    Matrix rref = matrix_rref(R);

    // First row should be pivot row [1, 2, 1]
    assert(approx(rref[0][0], 1.0));
    
    return 0;
}
