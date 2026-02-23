#ifndef LIELAB_DOMAIN_LIEALGEBRAS_cn_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_cn_HPP

#include "../VirtualManifolds.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <complex>
#include <limits>

namespace Lielab::domain
{

class cn : public VirtualLieAlgebra<std::complex<double>>
{
    public:
    // Manifold typing
    using point_t = Eigen::VectorXcd;
    
    // Manifold storage
    point_t point;

    // cn storage
    int _shape;

    // Manifold constructors
    cn();
    // ~cn();

    // Lie algebra constructors
    cn(const matrix_t& matrix);
    static cn basis(const int index, const int shape);
    static cn zero(const int shape);
    static cn from_vector(const Eigen::VectorXd& other);
    static cn from_vector(std::initializer_list<double> other);
    static cn project(const matrix_t& matrix);

    // cn constructors
    cn(const int n);
    static cn from_complex_vector(const Eigen::VectorXcd& other);
    static cn from_complex_vector(std::initializer_list<std::complex<double>> other);

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

    // cn IO
    Eigen::VectorXcd to_complex_vector() const; // TODO: Remove in favor of get_point()
    const field_t& operator[](const int index) const;
    field_t& operator[](const int index);

    // Lie algebra math
    cn operator+(const cn& other) const;
    cn& operator+=(const cn& other);
    cn operator-(const cn& other) const;
    cn& operator-=(const cn& other);
    cn operator-() const;
    cn operator*(const double other) const;
    friend cn operator*(const double other, const cn& rhs);
    cn& operator*=(const double other);
    cn operator/(const double other) const;
    cn& operator/=(const double other);

    // cn math
    cn operator*(const std::complex<int> other) const;
    cn operator*(const std::complex<double> other) const;
    friend cn operator*(const std::complex<int> other, const cn& rhs);
    friend cn operator*(const std::complex<double> other, const cn& rhs);
    cn& operator*=(const std::complex<int> other);
    cn& operator*=(const std::complex<double> other);
    cn operator/(const std::complex<int> other) const;
    cn operator/(const std::complex<double> other) const;
    cn& operator/=(const std::complex<int> other);
    cn& operator/=(const std::complex<double> other);
};

}

#endif
