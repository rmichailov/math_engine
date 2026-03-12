#include "linear_algebra.hpp"
#include <stdexcept>
#include <array>
#include <cmath>

namespace math_engine::linear_algebra {    
    // Vectors
    Vector vector_add(const Vector& v, const Vector& w) {
        if (v.size() != w.size()) {
            throw std::invalid_argument("Vectors must be same size");
        }
        Vector res;
        for (size_t i = 0; i < v.size(); i++) {
            res.push_back(v[i] + w[i]);
        }
        return res;
    }

    Vector vector_sub(const Vector& v, const Vector& w) {
        if (v.size() != w.size()) {
            throw std::invalid_argument("Vectors must be same size");
        }
        Vector res;
        for (size_t i = 0; i < v.size(); i++) {
            res.push_back(v[i] - w[i]);
        }
        return res;
    }

    Vector vector_scale(const Vector& v, double scale) {
        Vector res;
        for (size_t i = 0; i < v.size(); i++) {
            res.push_back(v[i] * scale);
        }
        return res;
    } 

    double vector_dot(const Vector& v, const Vector& w) {
        if (v.size() != w.size()) {
            throw std::invalid_argument("Vectors must be same size");
        }
        double res = 0;
        for (size_t i = 0; i < v.size(); i++) {
            res += (v[i] * w[i]);
        }
        return res;
    }

    double vector_norm(const Vector& v) {
        double sum = 0;
        for (size_t i = 0; i < v.size(); i++) {
            sum += (v[i] * v[i]);
        }
        return std::sqrt(sum);
    }

    Vector normalized_vector(const Vector& v) {
        double scale = vector_norm(v);
        return vector_scale(v, 1.0 / scale);
    }   


    // Matrices
    Matrix matrix_add(const Matrix& A, const Matrix& B) {
        if (A.size() != B.size() || A[0].size() != B[0].size()) {
            throw std::invalid_argument("Matrices must be same dimensions");
        }
        size_t rows = A.size();
        size_t cols = A[0].size();
        Matrix result(rows, std::vector<double>(cols));

        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                result[i][j] = (A[i][j] + B[i][j]);
            }
        }
        return result;   
    }

    Matrix matrix_sub(const Matrix& A, const Matrix& B) {
        if (A.size() != B.size() || A[0].size() != B[0].size()) {
            throw std::invalid_argument("Matrices must be same dimensions");
        }
        size_t rows = A.size();
        size_t cols = A[0].size();
        Matrix result(rows, std::vector<double>(cols));

        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                result[i][j] = (A[i][j] - B[i][j]);
            }
        }         
        return result;
    }

    Matrix matrix_scale(const Matrix& A, double scale) {
        size_t rows = A.size();
        size_t cols = A[0].size();
        Matrix result(rows, std::vector<double>(cols));
        
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                result[i][j] = (A[i][j] * scale);
            }
        }
        return result;            
    }  

    Matrix matrix_mul(const Matrix& A, const Matrix& B) {
        if (A[0].size() != B.size()) {
            throw std::invalid_argument("Matrices must be same dimensions");
        }
        size_t rows = A.size();
        size_t cols = B[0].size();
        Matrix res(rows, std::vector<double>(cols));
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                double sum = 0;
                for (size_t k = 0; k < B.size(); k++) {
                    sum +=(A[i][k] * B[k][j]);
                }
                res[i][j] = sum;
            }
        }
        return res;
    }
    Vector matrix_vector_mult(const Matrix& A, const Vector& v) {
        if (A[0].size() != v.size()) {
            throw std::invalid_argument("Matrix and vector dimensions must be same");
        }
        Vector res;
        for (size_t i = 0; i < A.size(); i++) {
            double sum = 0;
            for (size_t j = 0; j < v.size(); j++) {
                sum += (A[i][j] * v[j]);
            }
            res.push_back(sum);
        }
        return res;         
    }

    double matrix_determinant_recursion(const Matrix& A) {
        size_t n = A.size();
        if (n != A[0].size()) {
            throw std::invalid_argument("Matrix must be square");
        }

        if (n == 1) {
            return A[0][0];
        }
        if (n == 2) {
            return (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
        }

        double d = 0;
        for (size_t col = 0; col < n; col++) {
            Matrix sub_matrix;
            for (size_t i = 1; i < n; i++) {
                Vector row;
                for (size_t j = 0; j < n; j++) {
                    if (j != col) {
                        row.push_back(A[i][j]);
                    }
                }
                sub_matrix.push_back(row);
            }
            double sign = 1;
            if (col % 2 == 1) {
                sign = -1;
            }
            d += (sign * A[0][col] * matrix_determinant_recursion(sub_matrix));
        }

        return d;

    }

    
    Matrix matrix_transpose(const Matrix& A) {
        size_t rows = A.size();
        size_t cols = A[0].size();
        Matrix res(rows, std::vector<double>(cols));

        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                res[i][j] = A[j][i];
            }
        }
        return res;
    }

    
    Matrix matrix_inverse(const Matrix& A) {
        if (A.size() != A[0].size())
            throw std::invalid_argument("Matrix must be square");
        int64_t n = A.size();
        Matrix mat = A;
        Matrix inv(n, Vector(n, 0.0));
        for (int64_t i = 0; i < n; i++) {
            inv[i][i] = 1.0;
        }


        for (int64_t i = 0; i < n; i++) {
            int64_t pivot = i;
            for (int64_t j = i+1; j < n; j++) {
                if (std::abs(mat[j][i]) > std::abs(mat[pivot][i])) {
                    pivot = j;
                }
            }
            if (std::abs(mat[pivot][i]) < 1e-12) {
                throw std::runtime_error("Matrix is singular and cannot be inverted");
            }
            if (i != pivot) {
                std::swap(mat[i], mat[pivot]);
                std::swap(inv[i], inv[pivot]);
            }

            double factor = mat[i][i];
            for (int64_t j = 0; j < n; j++) {
                mat[i][j] /= factor;
                inv[i][j] /= factor;
            }
            for (int64_t j = 0; j < n; j++) {
                if (j == i) continue;
                double f = mat[j][i];
                for (int64_t k = 0; k < n; k++) {
                    mat[j][k] -= f * mat[i][k];
                    inv[j][k] -= f * inv[i][k];
                }
            }
        }
        return inv;        
    }

    Matrix matrix_rref(const Matrix& A) {
        const double epsilon = 1e-12;
        int64_t row = A.size();
        int64_t col = A[0].size();

        Matrix res = A;

        int64_t lead = 0;
        for (int64_t r = 0; r < row && lead < col; r++) {
            int64_t i = r;

            while (i < row && std::abs(A[i][lead]) < epsilon) {
                i++;
            }

            if (i == row) {
                lead++;
                r--;
                continue;
            }

            Vector temp = A[i];
            res[i] = A[r];
            res[r] = temp;

            double pivot = A[r][lead];
            for (int64_t j = 0; j < col; j++) {
                res[r][j] /= pivot;
            }

            for (int64_t j = 0; j < row; j++) {
                if (j == r) {
                    continue;
                }
                double factor = A[j][lead];
                for (int64_t k = 0; k < col; k++) {
                    res[j][k] -= (factor * A[r][k]);
                }
            }
            lead ++;
        }
        for (int64_t i = 0; i < row; i++) {
            for (int64_t j = 0; j < col; j++) {
                if (std::abs(res[i][j]) < epsilon) {
                    res[i][j] = 0;
                }
            }
        }

        return res;
    }

}