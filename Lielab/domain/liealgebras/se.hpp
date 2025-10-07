#ifndef LIELAB_DOMAIN_LIEALGEBRAS_se_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_se_HPP

#include "LieAlgebra.hpp"
#include "glr.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class se : public glr
{
    public:

    static constexpr bool abelian = false;

    std::string to_string() const override;
    se();
    se(const size_t n);
    template<typename OtherDerived> se(const Eigen::MatrixBase<OtherDerived> & other);
    static se basis(const ptrdiff_t i, const size_t n);
    static se from_shape(const size_t shape);

    size_t get_dimension() const;
    Eigen::VectorXd get_vector() const;
    void set_vector(const Eigen::VectorXd& vector);
    void set_vector(std::initializer_list<double> vector);
    se::matrix_t get_matrix() const;

    double operator()(const ptrdiff_t index) const;
    double operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    se operator+(const se& other) const;
    se& operator+=(const se& other);
    se operator-(const se& other) const;
    se& operator-=(const se& other);
    se operator-() const;
    se operator*(const double other) const;
    friend se operator*(const double other, const se& rhs);
    se& operator*=(const double other);
    se operator/(const double other) const;
    se& operator/=(const double other);

    static se from_vector(const Eigen::VectorXd& other);
    static se from_vector(std::initializer_list<double> other);

    static Eigen::MatrixXd project(const Eigen::MatrixXd& other);

    friend std::ostream& operator<<(std::ostream& os, const se& other);
};

}

#include "se.tpp"

#endif
