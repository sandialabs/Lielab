#ifndef LIELAB_DOMAIN_LIEALGEBRAS_cn_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_cn_HPP

#include "glc.hpp"
#include "LieAlgebra.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <complex>
#include <limits>

namespace Lielab::domain
{

class cn : public glc
{
    public:

    // Storage and typing
    using data_t = Eigen::VectorXcd;
    data_t data;

    // Lie Algebra class information
    bool is_abelian() const override;
    std::string to_string() const override;

    // Constructors and destructors
    cn();
    cn(const size_t n);
    template<typename OtherDerived> cn(const Eigen::MatrixBase<OtherDerived> & other);
    static cn basis(const ptrdiff_t i, const size_t n);
    static cn from_shape(const size_t shape);

    // Lie algebra object information
    size_t get_dimension() const;

    // IO and data manipulation
    cn::matrix_t get_matrix() const;
    Eigen::VectorXd get_vector() const;
    void set_vector(const Eigen::VectorXd& vector);
    void set_vector(std::initializer_list<double> vector);
    
    double operator()(const ptrdiff_t index) const;
    std::complex<double> operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;
    std::complex<double> operator[](const ptrdiff_t index) const;

    // Lie algebra math operations
    cn operator+(const cn & other) const;
    cn & operator+=(const cn & other);
    cn operator-(const cn & other) const;
    cn & operator-=(const cn & other);
    cn operator-() const;
    cn operator*(const double other) const;
    friend cn operator*(const double other, const cn & rhs);
    cn operator*(const std::complex<int> other) const;
    cn operator*(const std::complex<double> other) const;
    friend cn operator*(const std::complex<int> other, const cn & rhs);
    friend cn operator*(const std::complex<double> other, const cn & rhs);
    cn & operator*=(const double other);
    cn & operator*=(const std::complex<int> other);
    cn & operator*=(const std::complex<double> other);
    cn operator/(const double other) const;
    cn operator/(const std::complex<int> other) const;
    cn operator/(const std::complex<double> other) const;
    cn & operator/=(const double other);
    cn & operator/=(const std::complex<int> other);
    cn & operator/=(const std::complex<double> other);

    // Extra methods
    static cn from_vector(const Eigen::VectorXd & other);
    static cn from_vector(std::initializer_list<double> other);
    static cn from_complex_vector(const Eigen::VectorXcd & other);
    static cn from_complex_vector(std::initializer_list<std::complex<double>> other);
    Eigen::VectorXcd to_complex_vector() const;

    static Eigen::MatrixXcd project(const Eigen::MatrixXcd & other);

    friend std::ostream& operator<<(std::ostream & os, const cn & other);
};

}

#include "cn.tpp"

#endif
