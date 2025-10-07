#ifndef LIELAB_DOMAIN_LIEALGEBRAS_glc_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_glc_HPP

#include "LieAlgebra.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <complex>

namespace Lielab::domain
{

class glc : public LieAlgebra<std::complex<double>>
{
    public:

    static constexpr bool abelian = false;

    size_t _shape = 0;

    std::string to_string() const override;
    glc();
    glc(const size_t n);
    template<typename OtherDerived> glc(const Eigen::MatrixBase<OtherDerived> & other);
    static glc basis(const ptrdiff_t i, const size_t n);
    static glc from_shape(const size_t shape);

    size_t get_dimension() const;
    size_t get_shape() const;
    Eigen::VectorXd get_vector() const;
    void set_vector(const Eigen::VectorXd& vector);
    void set_vector(std::initializer_list<double> vector);
    glc::matrix_t get_matrix() const;

    double operator()(const ptrdiff_t index) const;
    std::complex<double> operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    glc operator+(const glc & other) const;
    glc & operator+=(const glc & other);
    glc operator-(const glc & other) const;
    glc & operator-=(const glc & other);
    glc operator-() const;
    glc operator*(const double other) const;
    friend glc operator*(const double other, const glc & rhs);
    glc operator*(const std::complex<double> other) const;
    friend glc operator*(const std::complex<double> other, const glc & rhs);
    glc & operator*=(const double other);
    glc & operator*=(const std::complex<double> other);
    glc operator/(const double other) const;
    glc operator/(const std::complex<double> other) const;
    glc & operator/=(const double other);
    glc & operator/=(const std::complex<double> other);

    static glc from_vector(const Eigen::VectorXd & other);
    static glc from_vector(std::initializer_list<double> other);
    static glc from_complex_vector(const Eigen::VectorXcd & other);
    static glc from_complex_vector(std::initializer_list<std::complex<double>> other);

    static Eigen::MatrixXcd project(const Eigen::MatrixXcd& other);

    friend std::ostream& operator<<(std::ostream & os, const glc & other);
};

}

#include "glc.tpp"

#endif
