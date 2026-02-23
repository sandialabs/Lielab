#ifndef LIELAB_DOMAIN_LIEALGEBRAS_rn_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_rn_HPP

#include "../VirtualManifolds.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class rn : public VirtualLieAlgebra<double>
{
    public:
    // Manifold typing
    using point_t = Eigen::VectorXd;

    // Manifold storage
    point_t point;

    // rn storage
    int _shape = 0;

    // Manifold constructors
    rn();
    // ~rn()

    // Lie algebra constructors
    rn(const matrix_t& matrix);
    static rn basis(const int index, const int shape);
    static rn zero(const int shape);
    static rn from_vector(const Eigen::VectorXd& other);
    static rn from_vector(std::initializer_list<double> other);
    static rn project(const matrix_t& other);

    // rn constructors
    rn(const int n);

    // Manifold information
    std::string to_string() const override;
    int get_dimension() const override;
    int get_size() const override;

    // Lie algebra information
    bool is_abelian() const override;
    int get_shape() const override;

    // Manifold IO
    point_t get_point() const;
    Eigen::VectorXd serialize() const override;
    void unserialize(const Eigen::VectorXd& serialized) override;
    void unserialize(std::initializer_list<double> serialized) override;

    // Lie algebra IO
    matrix_t get_matrix() const;
    Eigen::VectorXd get_vector() const override;
    void set_vector(const Eigen::VectorXd& vector) override;
    void set_vector(std::initializer_list<double> vector) override;
    double operator()(const int index) const;
    field_t operator()(const int index1, const int index2) const;

    // rn IO
    const field_t& operator[](const int index) const;
    field_t& operator[](const int index);

    // Lie algebra math
    rn operator+(const rn& other) const;
    rn& operator+=(const rn& other);
    rn operator-(const rn& other) const;
    rn& operator-=(const rn& other);
    rn operator-() const;
    rn operator*(const double other) const;
    friend rn operator*(const double other, const rn& rhs);
    rn& operator*=(const double other);
    rn operator/(const double other) const;
    rn& operator/=(const double other);
};

}

#endif
