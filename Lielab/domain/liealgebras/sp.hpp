#ifndef LIELAB_DOMAIN_LIEALGEBRAS_sp_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_sp_HPP

#include "LieAlgebra.hpp"
#include "glr.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class sp : public glr
{
    public:

    static constexpr bool abelian = false;

    std::string to_string() const override;
    sp();
    sp(const size_t n);
    template<typename OtherDerived> sp(const Eigen::MatrixBase<OtherDerived> & other);
    static sp basis(const ptrdiff_t i, const size_t n);
    static sp from_shape(const size_t shape);

    size_t get_dimension() const;
    Eigen::VectorXd get_vector() const;
    void set_vector(const Eigen::VectorXd& vector);
    void set_vector(std::initializer_list<double> vector);
    sp::matrix_t get_matrix() const;

    double operator()(const ptrdiff_t index) const;
    double operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    sp operator+(const sp & other) const;
    sp & operator+=(const sp & other);
    sp operator-(const sp & other) const;
    sp & operator-=(const sp & other);
    sp operator-() const;
    sp operator*(const double other) const;
    friend sp operator*(const double other, const sp & rhs);
    sp & operator*=(const double other);
    sp operator/(const double other) const;
    sp & operator/=(const double other);

    static sp from_vector(const Eigen::VectorXd & other);
    static sp from_vector(std::initializer_list<double> other);

    static Eigen::MatrixXd project(const Eigen::MatrixXd & other);

    friend std::ostream& operator<<(std::ostream & os, const sp & other);
};

}

#include "sp.tpp"

#endif
