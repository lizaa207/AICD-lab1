#include <complex>
#include <cstddef>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>

double AbsValue(int value) { return std::abs(static_cast<double>(value)); }
double AbsValue(float value) { return static_cast<double>(std::abs(value)); }
double AbsValue(double value) { return std::abs(value); }

template <typename T>
double AbsValue(const std::complex<T>& value) {
    return static_cast<double>(std::abs(value));
}

int GetRandom(std::mt19937& generator, int lower, int upper) {
    if (lower > upper) {
        throw std::invalid_argument("Lower bound is greater than upper bound");
    }
    std::uniform_int_distribution<int> distribution(lower, upper);
    return distribution(generator);
}

float GetRandom(std::mt19937& generator, float lower, float upper) {
    if (lower > upper) {
        throw std::invalid_argument("Lower bound is greater than upper bound");
    }
    std::uniform_real_distribution<float> distribution(lower, upper);
    return distribution(generator);
}

double GetRandom(std::mt19937& generator, double lower, double upper) {
    if (lower > upper) {
        throw std::invalid_argument("Lower bound is greater than upper bound");
    }
    std::uniform_real_distribution<double> distribution(lower, upper);
    return distribution(generator);
}

template <typename T>
std::complex<T> GetRandom(std::mt19937& generator, const std::complex<T>& lower,
    const std::complex<T>& upper) {
    T real = GetRandom(generator, lower.real(), upper.real());
    T imag = GetRandom(generator, lower.imag(), upper.imag());
    return std::complex<T>(real, imag);
}

template <typename T>
class Matrix {
public:
    static constexpr double kEpsilon = 1e-6;

    Matrix(std::size_t rows, std::size_t cols, const T& value)
        : rows_(rows), cols_(cols), data_(nullptr) {
        CheckSize();
        data_ = new T[rows_ * cols_];
        for (std::size_t i = 0; i < rows_ * cols_; ++i) {
            data_[i] = value;
        }
    }

    Matrix(std::size_t rows, std::size_t cols, const T& lower, const T& upper)
        : rows_(rows), cols_(cols), data_(nullptr) {
        CheckSize();
        data_ = new T[rows_ * cols_];
        std::random_device device;
        std::mt19937 generator(device());
        try {
            for (std::size_t i = 0; i < rows_ * cols_; ++i) {
                data_[i] = GetRandom(generator, lower, upper);
            }
        }
        catch (...) {
            delete[] data_;
            throw;
        }
    }

    Matrix(const Matrix& other)
        : rows_(other.rows_),
        cols_(other.cols_),
        data_(new T[other.rows_ * other.cols_]) {
        for (std::size_t i = 0; i < rows_ * cols_; ++i) {
            data_[i] = other.data_[i];
        }
    }

    Matrix& operator=(const Matrix& other) {
        if (this == &other) {
            return *this;
        }
        T* new_data = new T[other.rows_ * other.cols_];
        for (std::size_t i = 0; i < other.rows_ * other.cols_; ++i) {
            new_data[i] = other.data_[i];
        }
        delete[] data_;
        data_ = new_data;
        rows_ = other.rows_;
        cols_ = other.cols_;
        return *this;
    }

    ~Matrix() { delete[] data_; }

    std::size_t GetRows() const { return rows_; }
    std::size_t GetCols() const { return cols_; }

    T& operator()(std::size_t row, std::size_t col) {
        CheckIndex(row, col);
        return data_[row * cols_ + col];
    }

    const T& operator()(std::size_t row, std::size_t col) const {
        CheckIndex(row, col);
        return data_[row * cols_ + col];
    }
        Matrix & operator+=(const Matrix & other) {
        CheckSameSize(other);
        for (std::size_t i = 0; i < rows_ * cols_; ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }

    Matrix& operator-=(const Matrix& other) {
        CheckSameSize(other);
        for (std::size_t i = 0; i < rows_ * cols_; ++i) {
            data_[i] -= other.data_[i];
        }
        return *this;
    }

    Matrix& operator*=(const T& scalar) {
        for (std::size_t i = 0; i < rows_ * cols_; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }

    Matrix& operator/=(const T& scalar) {
        if (scalar == T()) {
            throw std::invalid_argument("Division by zero");
        }
        for (std::size_t i = 0; i < rows_ * cols_; ++i) {
            data_[i] /= scalar;
        }
        return *this;
    }

    Matrix operator+(const Matrix& other) const {
        Matrix result(*this);
        result += other;
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        Matrix result(*this);
        result -= other;
        return result;
    }

    Matrix operator*(const T& scalar) const {
        Matrix result(*this);
        result *= scalar;
        return result;
    }

    friend Matrix operator*(const T& scalar, const Matrix& matrix) {
        return matrix * scalar;
    }

    Matrix operator/(const T& scalar) const {
        Matrix result(*this);
        result /= scalar;
        return result;
    }
     private:
         void CheckSize() const {
             if (rows_ == 0 || cols_ == 0) {
                 throw std::invalid_argument("Matrix size must be positive");
             }
         }

         void CheckIndex(std::size_t row, std::size_t col) const {
             if (row >= rows_ || col >= cols_) {
                 throw std::out_of_range("Matrix index out of range");
             }
         }
         std::size_t rows_;
         std::size_t cols_;
         T* data_;
};
void CheckSameSize(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("Matrices have different sizes");
    }
}
Matrix operator*(const Matrix& other) const {
    if (cols_ != other.rows_) {
        throw std::invalid_argument("Wrong sizes for matrix multiplication");
    }
    Matrix result(rows_, other.cols_, T());
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < other.cols_; ++j) {
            for (std::size_t k = 0; k < cols_; ++k) {
                result(i, j) += (*this)(i, k) * other(k, j);
            }
        }
    }
    return result;
}

T Trace() const {
    if (rows_ != cols_) {
        throw std::invalid_argument("Matrix is not square");
    }
    T sum = T();
    for (std::size_t i = 0; i < rows_; ++i) {
        sum += (*this)(i, i);
    }
    return sum;
}

bool operator==(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        return false;
    }
    for (std::size_t i = 0; i < rows_ * cols_; ++i) {
        if (AbsValue(data_[i] - other.data_[i]) > kEpsilon) {
            return false;
        }
    }
    return true;
}

bool operator!=(const Matrix& other) const { return !(*this == other); }
template <typename T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& matrix) {
    for (std::size_t i = 0; i < matrix.GetRows(); ++i) {
        for (std::size_t j = 0; j < matrix.GetCols(); ++j) {
            out << matrix(i, j) << '\t';
        }
        out << '\n';
    }
    return out;
}