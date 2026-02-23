#ifndef LIELAB_DOMAIN_RN_HPP
#define LIELAB_DOMAIN_RN_HPP

#include "../VirtualManifolds.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{
class RN : public VirtualLieGroup<double>
{
    public:
    // Manifold typing
    using point_t = Eigen::VectorXd;

    // Manifold storage
    point_t point;

    // RN storage
    int _shape = 0;

    // Manifold constructors
    RN();
    // ~RN();

    // Lie group constructors
    RN(const matrix_t& matrix);
    static RN identity(const int shape);
    static RN project(const matrix_t& matrix);

    // RN constructors
    RN(const int n);
    static RN from_vector(const Eigen::VectorXd& other);
    static RN from_vector(const std::initializer_list<double> other);

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

    // RN IO
    const field_t& operator[](const int index) const;
    field_t& operator[](const int index);

    // Lie group math
    RN operator*(const RN& other) const;
    RN& operator*=(const RN& other);
    RN inverse() const;
};

}

#endif
