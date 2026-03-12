#pragma once
#include <cstdint>
#include <vector>
#include <cmath>

namespace math_engine::linear_algebra {
    // Vector Operations
    using Vector = std::vector<double>;
    Vector vector_add(const Vector& v, const Vector& w);
    Vector vector_sub(const Vector& v, const Vector& w);
    Vector vector_scale(const Vector& v, double scale);
    double vector_dot(const Vector& v, const Vector& w);
    double vector_norm(const Vector& v);
    Vector normalized_vector(const Vector& v);


    // Matrix operations
    using Matrix = std::vector<std::vector<double>>;
    Matrix matrix_add(const Matrix& A, const Matrix& B);
    Matrix matrix_sub(const Matrix& A, const Matrix& B);
    Matrix matrix_scale(const Matrix& A, double scale);
    Matrix matrix_mul(const Matrix& A, const Matrix& B);
    Vector matrix_vector_mult(const Matrix& A, const Vector& v);

    double matrix_determinant_recursion(const Matrix& A);
    Matrix matrix_transpose(const Matrix& A);
    Matrix matrix_inverse(const Matrix& A);
    Matrix matrix_rref(const Matrix& A);

    int64_t matrix_rank(const Matrix& A);
}