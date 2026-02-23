#ifndef LIELAB_DOMAIN_LIEALGEBRAS_se_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_se_HPP

#include "../VirtualManifolds.hpp"

#include "rn.hpp"
#include "so.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <tuple>

namespace Lielab::domain
{

class se : public VirtualLieAlgebra<double>
{
    public:
    // Manifold typing
    using point_t = std::tuple<rn, so>;

    // Manifold storage
    point_t point;

    // se storage
    int _shape = 0;

    // Manifold constructors
    se();
    // ~se();

    // Lie algebra constructors
    se(const matrix_t& matrix);
    static se basis(const int index, const int shape);
    static se zero(const int shape);
    static se from_vector(const Eigen::VectorXd& other);
    static se from_vector(std::initializer_list<double> other);
    static se project(const matrix_t& matrix);

    // se constructors
    se(const int n);
    se(const rn& rn_c, const so& so_c);

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
};

}

#endif
