#ifndef LIELAB_DOMAIN_CN_HPP
#define LIELAB_DOMAIN_CN_HPP

#include "../VirtualManifolds.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class CN : public VirtualLieGroup<std::complex<double>>
{
    public:
    // Manifold typing
    using point_t = Eigen::VectorXcd;

    // Manifold storage
    point_t point;

    // CN storage
    int _shape = 0;

    // Manifold constructors
    CN();
    // ~CN();

    // Lie group constructors
    CN(const matrix_t& matrix);
    static CN identity(const int shape);
    static CN project(const matrix_t& matrix);

    // CN constructors
    CN(const int n);
    static CN from_vector(const Eigen::VectorXd& other); // TODO: remove?
    static CN from_vector(std::initializer_list<double> other); // TODO: remove?
    static CN from_complex_vector(const Eigen::VectorXcd& other);
    static CN from_complex_vector(const std::initializer_list<std::complex<double>> other);

    // Manifold information
    std::string to_string() const override;
    int get_dimension() const override;
    int get_size() const override;

    // Lie group information
    bool is_abelian() const override;
    int get_shape() const override;

    // Manifold IO
    point_t get_point() const;
    Eigen::VectorXd serialize() const override;
    void unserialize(const Eigen::VectorXd& serialized) override;
    void unserialize(std::initializer_list<double> serialized) override;
    
    // Lie group IO
    matrix_t get_matrix() const;
    field_t operator()(const int index1, const int index2) const;

    // CN IO
    Eigen::VectorXcd to_complex_vector() const;
    const field_t& operator[](const int index) const;
    field_t& operator[](const int index);

    // Lie group math
    CN operator*(const CN& other) const;
    CN& operator*=(const CN& other);
    CN inverse() const;

    // CN math
    // CN operator+(const CN& other) const;
    // CN& operator+=(const CN& other);
    // CN operator-(const CN& other) const;
    // CN& operator-=(const CN& other);
    // and doubles, std::complex, etc...
};

}

#endif
