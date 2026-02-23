#ifndef LIELAB_DOMAIN_LIEALGEBRAS_sp_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_sp_HPP

#include "../VirtualManifolds.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class sp : public VirtualLieAlgebra<double>
{
    public:
    // Manifold typing
    using point_t = Eigen::MatrixXd;

    // Manifold storage
    point_t point;

    // Manifold constructors
    sp();
    // ~sp();

    // Lie algebra constructors
    sp(const matrix_t& matrix);
    static sp basis(const int index, const int shape);
    static sp zero(const int shape);
    static sp from_vector(const Eigen::VectorXd& other);
    static sp from_vector(std::initializer_list<double> other);
    static sp project(const matrix_t& matrix);

    // cn constructors
    sp(const int n);

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

    // Lie algebra math
    sp operator+(const sp& other) const;
    sp& operator+=(const sp& other);
    sp operator-(const sp& other) const;
    sp& operator-=(const sp& other);
    sp operator-() const;
    sp operator*(const double other) const;
    friend sp operator*(const double other, const sp& rhs);
    sp& operator*=(const double other);
    sp operator/(const double other) const;
    sp& operator/=(const double other);
};

}

#endif
