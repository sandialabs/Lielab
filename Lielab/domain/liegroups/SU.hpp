#ifndef LIELAB_DOMAIN_SU_HPP
#define LIELAB_DOMAIN_SU_HPP

#include "../VirtualManifolds.hpp"
#include "SO.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <complex>

namespace Lielab::domain
{

class SO;

class SU : public VirtualLieGroup<std::complex<double>>
{
    public:
    // Manifold typing
    using point_t = Eigen::MatrixXcd;
    
    // Manifold storage
    point_t point;

    // Manifold constructors
    SU();
    // ~SU();

    // Lie group constructors
    SU(const matrix_t& other);
    static SU identity(const int shape);
    // TODO: Project

    // SU constructors
    SU(const int shape);
    static SU from_quaternion(const double e0, const double e1, const double e2, const double e3);
    static SU from_SO3(const SO& dcm);
    
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

    // SU IO
    std::array<double, 4> to_quaternion() const;

    // Lie group math
    SU operator*(const SU& other) const;
    SU& operator*=(const SU& other);
    SU inverse() const;
};

}

#endif
