#ifndef LIELAB_DOMAIN_LIEALGEBRAS_rn_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_rn_HPP

#include "glr.hpp"
#include "LieAlgebra.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class rn : public glr
{
    public:

    static constexpr bool abelian = true;

    using data_t = Eigen::VectorXd;

    data_t data;

    std::string to_string() const override;
    rn();
    rn(const size_t n);
    template<typename OtherDerived> rn(const Eigen::MatrixBase<OtherDerived> & other);
    static rn basis(const ptrdiff_t i, const size_t n);
    static rn from_shape(const size_t shape);

    size_t get_dimension() const;
    Eigen::VectorXd get_vector() const;
    void set_vector(const Eigen::VectorXd& vector);
    void set_vector(std::initializer_list<double> vector);
    rn::matrix_t get_matrix() const;

    double operator()(const ptrdiff_t index) const;
    double operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    rn operator+(const rn & other) const;
    rn & operator+=(const rn & other);
    rn operator-(const rn & other) const;
    rn & operator-=(const rn & other);
    rn operator-() const;
    rn operator*(const double other) const;
    friend rn operator*(const double other, const rn & rhs);
    rn & operator*=(const double other);
    rn operator/(const double other) const;
    rn & operator/=(const double other);

    static rn from_vector(const Eigen::VectorXd& other);
    static rn from_vector(std::initializer_list<double> other);

    static Eigen::MatrixXd project(const Eigen::MatrixXd & other);

    friend std::ostream& operator<<(std::ostream & os, const rn& other);
};

}

#include "rn.tpp"

#endif
