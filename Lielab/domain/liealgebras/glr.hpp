#ifndef LIELAB_DOMAIN_LIEALGEBRAS_glr_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_glr_HPP

#include "LieAlgebra.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class glr : public LieAlgebra<double>
{
    public:

    static constexpr bool abelian = false;

    size_t _shape = 0;

    std::string to_string() const override;
    glr();
    glr(const size_t n);
    template<typename OtherDerived> glr(const Eigen::MatrixBase<OtherDerived> & other);
    static glr basis(const ptrdiff_t i, const size_t n);
    static glr from_shape(const size_t shape);

    static Eigen::MatrixXd project(const Eigen::MatrixXd & other);

    size_t get_dimension() const;
    size_t get_shape() const;
    Eigen::VectorXd get_vector() const;
    void set_vector(const Eigen::VectorXd & vector);
    void set_vector(std::initializer_list<double> vector);
    glr::matrix_t get_matrix() const;

    double operator()(const ptrdiff_t index) const;
    double operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    glr operator+(const glr & other) const;
    glr & operator+=(const glr & other);
    glr operator-(const glr & other) const;
    glr & operator-=(const glr & other);
    glr operator-() const;
    glr operator*(const double other) const;
    friend glr operator*(const double other, const glr & rhs);
    glr & operator*=(const double other);
    glr operator/(const double other) const;
    glr & operator/=(const double other);

    static glr from_vector(const Eigen::VectorXd & other);
    static glr from_vector(std::initializer_list<double> other);

    friend std::ostream& operator<<(std::ostream & os, const glr & other);
};

}

#include "glr.tpp"

#endif
