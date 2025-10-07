#ifndef LIELAB_DOMAIN_LIEALGEBRAS_so_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_so_HPP

#include "LieAlgebra.hpp"
#include "glr.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{
class so : public glr
{
    public:

    static constexpr bool abelian = false;

    std::string to_string() const override;
    so();
    so(const size_t n);
    template<typename OtherDerived> so(const Eigen::MatrixBase<OtherDerived> & other);
    static so basis(const ptrdiff_t i, const size_t n);
    static so from_shape(const size_t shape);

    size_t get_dimension() const;
    Eigen::VectorXd get_vector() const;
    void set_vector(const Eigen::VectorXd& vector);
    void set_vector(std::initializer_list<double> vector);
    so::matrix_t get_matrix() const;

    double operator()(const ptrdiff_t index) const;
    double operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    so operator+(const so& other) const;
    so& operator+=(const so& other);
    so operator-(const so& other) const;
    so& operator-=(const so& other);
    so operator-() const;
    so operator*(const double other) const;
    friend so operator*(const double other, const so& rhs);
    so& operator*=(const double other);
    so operator/(const double other) const;
    so& operator/=(const double other);

    static so from_vector(const Eigen::VectorXd& other);
    static so from_vector(std::initializer_list<double> other);

    static Eigen::MatrixXd project(const Eigen::MatrixXd& other);

    friend std::ostream& operator<<(std::ostream& os, const so& other);
};

}

#include "so.tpp"

#endif
