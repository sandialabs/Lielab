#ifndef LIELAB_DOMAIN_GLC_HPP
#define LIELAB_DOMAIN_GLC_HPP

#include "../VirtualManifolds.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class GLC : public VirtualLieGroup<std::complex<double>>
{
    public:
    // Manifold typing
    using point_t = Eigen::MatrixXcd;

    // Manifold storage
    point_t point;

    // Manifold constructors
    GLC();
    // ~GLC();

    // Lie group constructors
    GLC(const matrix_t& matrix);
    static GLC identity(const int shape);
    static GLC project(const matrix_t& matrix);

    // GLC constructors
    GLC(const int shape);

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

    // Lie group math
    GLC operator*(const GLC& other) const;
    GLC& operator*=(const GLC& other);
    GLC inverse() const;
};

}

#endif
