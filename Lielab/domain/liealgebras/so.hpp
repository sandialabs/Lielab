#ifndef LIELAB_DOMAIN_LIEALGEBRAS_so_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_so_HPP

#include "../VirtualManifolds.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{
class so : public VirtualLieAlgebra<double>
{
    public:
    // Manifold typing
    using point_t = Eigen::MatrixXd;

    // Manifold storage
    point_t point;

    // Manifold constructors
    so();
    // ~so();

    // Lie algebra constructors
    so(const matrix_t& other);
    static so basis(const int index, const int shape);
    static so zero(const int shape);
    static so from_vector(const Eigen::VectorXd& other);
    static so from_vector(std::initializer_list<double> other);
    static so project(const matrix_t& other);

    // cn constructors
    so(const int n);

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
};

}

#endif
