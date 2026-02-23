#ifndef LIELAB_DOMAIN_LIEALGEBRAS_su_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_su_HPP

#include "../VirtualManifolds.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{
/*!
    * The su class.
    *
    * Defined by skew-Hermitian matrices.
    */
class su : public VirtualLieAlgebra<std::complex<double>>
{
    public:
    // Manifold typing
    using point_t = Eigen::MatrixXcd;

    // Manifold storage
    point_t point;

    // Manifold constructors
    su();
    // ~su();

    // Lie algebra constructors
    su(const matrix_t& matrix);
    static su basis(const int index, const int shape);
    static su zero(const int shape);
    static su from_vector(const Eigen::VectorXd& vector);
    static su from_vector(std::initializer_list<double> vector);
    static su project(const matrix_t& other);

    // su constructors
    su(const int n);

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
    void unserialize(const Eigen::VectorXd& serialized_rep) override;
    void unserialize(std::initializer_list<double> serialized) override;

    // Lie algebra IO
    matrix_t get_matrix() const;
    Eigen::VectorXd get_vector() const override;
    void set_vector(const Eigen::VectorXd& vector) override;
    void set_vector(std::initializer_list<double> vector) override;
    double operator()(const int index) const;
    field_t operator()(const int index1, const int index2) const;

    // Lie algebra math
    su operator+(const su& other) const;
    su& operator+=(const su& other);
    su operator-(const su& other) const;
    su& operator-=(const su& other);
    su operator-() const;
    su operator*(const double other) const;
    su& operator*=(const double other);
    friend su operator*(const double other, const su& rhs);
    su operator/(const double other) const;
    su& operator/=(const double other);
};

}

#endif
