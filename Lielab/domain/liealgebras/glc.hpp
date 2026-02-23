#ifndef LIELAB_DOMAIN_LIEALGEBRAS_glc_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_glc_HPP

#include "../VirtualManifolds.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <complex>

namespace Lielab::domain
{

class glc : public VirtualLieAlgebra<std::complex<double>>
{
    public:
    // Manifold typing
    using point_t = Eigen::MatrixXcd;

    // Manifold storage
    point_t point;

    // Manifold constructors
    glc();
    // ~glc();

    // Lie algebra constructors
    glc(const matrix_t& matrix);
    static glc basis(const int index, const int shape);
    static glc zero(const int shape);
    static glc from_vector(const Eigen::VectorXd& other);
    static glc from_vector(std::initializer_list<double> other);
    static glc project(const matrix_t& matrix);

    // glc constructors
    glc(const int n);
    static glc from_complex_vector(const Eigen::VectorXcd& other);
    static glc from_complex_vector(std::initializer_list<std::complex<double>> other);

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
    glc operator+(const glc& other) const;
    glc& operator+=(const glc& other);
    glc operator-(const glc& other) const;
    glc& operator-=(const glc& other);
    glc operator-() const;
    glc operator*(const double other) const;
    friend glc operator*(const double other, const glc& rhs);
    glc& operator*=(const double other);
    glc operator/(const double other) const;
    glc& operator/=(const double other);

    // glc math
    glc operator*(const std::complex<double> other) const;
    friend glc operator*(const std::complex<double> other, const glc& rhs);
    glc& operator*=(const std::complex<double> other);
    glc operator/(const std::complex<double> other) const;
    glc& operator/=(const std::complex<double> other);
};

}

#endif
