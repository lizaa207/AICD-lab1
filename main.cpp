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

    Matrix operator/(const T& scalar) const {
        Matrix result(*this);
        result /= scalar;
        return result;
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

    void CheckSameSize(const Matrix& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Matrices have different sizes");
        }
    }

    std::size_t rows_;
    std::size_t cols_;
    T* data_;
};

    
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
template <typename T>
void DecomposeLu(const Matrix<T>& matrix, Matrix<T>& lower, Matrix<T>& upper) {
    static_assert(!std::is_integral<T>::value,
        "LU decomposition requires non-integer type");
    if (matrix.GetRows() != matrix.GetCols()) {
        throw std::invalid_argument("Matrix is not square");
    }
    const std::size_t size = matrix.GetRows();
    lower = Matrix<T>(size, size, T());
    upper = matrix;
    for (std::size_t i = 0; i < size; ++i) {
        lower(i, i) = T(1);
    }

    for (std::size_t k = 0; k < size; ++k) {
        const T pivot = upper(k, k);
        if (AbsValue(pivot) < Matrix<T>::kEpsilon) {
            throw std::domain_error("Zero pivot, decomposition is impossible");
        }
        for (std::size_t i = k + 1; i < size; ++i) {
            const T factor = upper(i, k) / pivot;
            lower(i, k) = factor;
            for (std::size_t j = k; j < size; ++j) {
                upper(i, j) -= factor * upper(k, j);
            }
            upper(i, k) = T();
        }
    }
}
template <typename T>
void DemonstrateMatrix(const std::string& name, const T& lower,
    const T& upper) {
    std::cout << "===== " << name << " =====\n";

    Matrix<T> filled(2, 3, T(1));
    Matrix<T> a(3, 3, lower, upper);
    Matrix<T> b(3, 3, lower, upper);
    std::cout << "Filled:\n" << filled << "A:\n" << a << "B:\n" << b;

    Matrix<T> copy(a);
    copy(0, 0) = upper;
    std::cout << "Copy of A with changed (0, 0):\n" << copy;
    copy = b;
    std::cout << "Copy after assignment of B:\n" << copy;

    std::cout << "A + B:\n" << a + b << "A - B:\n" << a - b;
    std::cout << "A * B:\n" << a * b;

    const T scalar = T(2);
    std::cout << "A * 2:\n" << a * scalar << "2 * A:\n" << scalar * a;
    std::cout << "A / 2:\n" << a / scalar;
    std::cout << "Trace of A: " << a.Trace() << '\n';
    std::cout << "A == A: " << (a == a) << ", A != B: " << (a != b) << '\n';

    try {
        Matrix<T> sum = filled + a;
    }
    catch (const std::invalid_argument& error) {
        std::cout << "Error: " << error.what() << '\n';
    }
    try {
        std::cout << filled(5, 5) << '\n';
    }
    catch (const std::out_of_range& error) {
        std::cout << "Error: " << error.what() << '\n';
    }
    try {
        std::cout << filled.Trace() << '\n';
    }
    catch (const std::invalid_argument& error) {
        std::cout << "Error: " << error.what() << '\n';
    }
    try {
        Matrix<T> quotient = a / T();
    }
    catch (const std::invalid_argument& error) {
        std::cout << "Error: " << error.what() << '\n';
    }
    std::cout << '\n';
}
template <typename T>
void DemonstrateLu(const std::string& name, const T& lower_bound,
    const T& upper_bound) {
    std::cout << "===== LU, " << name << " =====\n";
    Matrix<T> a(4, 4, lower_bound, upper_bound);
    Matrix<T> lower(1, 1, T());
    Matrix<T> upper(1, 1, T());
    std::cout << "A:\n" << a;
    try {
        DecomposeLu(a, lower, upper);
        std::cout << "L:\n" << lower << "U:\n" << upper;
        std::cout << "L * U:\n" << lower * upper;
        std::cout << "L * U == A: " << (lower * upper == a) << '\n';
    }
    catch (const std::domain_error& error) {
        std::cout << "Error: " << error.what() << '\n';
    }
    std::cout << '\n';
}

int main() {
    DemonstrateMatrix<int>("int", 1, 9);
    DemonstrateMatrix<float>("float", -5.0F, 5.0F);
    DemonstrateMatrix<double>("double", -5.0, 5.0);
    DemonstrateMatrix<std::complex<float>>("complex<float>",
        std::complex<float>(-5.0F, -5.0F),
        std::complex<float>(5.0F, 5.0F));
    DemonstrateMatrix<std::complex<double>>("complex<double>",
        std::complex<double>(-5.0, -5.0),
        std::complex<double>(5.0, 5.0));

    DemonstrateLu<double>("double", -5.0, 5.0);
    DemonstrateLu<std::complex<double>>("complex<double>",
        std::complex<double>(-5.0, -5.0),
        std::complex<double>(5.0, 5.0));
    return 0;
}