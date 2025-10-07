#ifndef LIELAB_DOMAIN_LIEALGEBRAS_su_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_su_HPP

#include "LieAlgebra.hpp"
#include "glc.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{
/*!
    * The su class.
    *
    * Defined by skew-Hermitian matrices.
    */
class su : public glc
{
    public:

    static constexpr bool abelian = false;

    std::string to_string() const override;
    su();
    su(const size_t n);
    template<typename OtherDerived> su(const Eigen::MatrixBase<OtherDerived> & other);
    static su basis(const ptrdiff_t i, const size_t n);
    static su from_shape(const size_t shape);

    size_t get_dimension() const;
    Eigen::VectorXd get_vector() const;
    void set_vector(const Eigen::VectorXd& vector);
    void set_vector(std::initializer_list<double> vector);
    su::matrix_t get_matrix() const;

    double operator()(const ptrdiff_t index) const;
    std::complex<double> operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    su operator+(const su & other) const;
    su & operator+=(const su & other);
    su operator-(const su & other) const;
    su & operator-=(const su & other);
    su operator-() const;
    su operator*(const double other) const;
    su operator*(const std::complex<double> other) const;
    su & operator*=(const double other);
    friend su operator*(const double other, const su & rhs);
    su & operator*=(const std::complex<double> other);
    friend su operator*(const std::complex<double> other, const su & rhs);
    su operator/(const double other) const;
    su operator/(const std::complex<double> other) const;
    su & operator/=(const double other);
    su & operator/=(const std::complex<double> other);

    static su from_vector(const Eigen::VectorXd& vector);
    static su from_vector(std::initializer_list<double> vector);

    // TODO: Project function

    friend std::ostream & operator<<(std::ostream& os, const su & other);
};

}

#include "su.tpp"

#endif
