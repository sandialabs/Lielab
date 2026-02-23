#ifndef LIELAB_DOMAIN_GL_HPP
#define LIELAB_DOMAIN_GL_HPP

#include "../VirtualManifolds.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class GLR : public VirtualLieGroup<double>
{
    public:
    // Manifold typing
    using point_t = Eigen::MatrixXd;

    // Manifold storage
    point_t point;

    // Manifold constructors
    GLR();
    // ~GLR();

    // Lie group constructors
    GLR(const matrix_t& matrix);
    static GLR identity(const int shape);
    static GLR project(const matrix_t& matrix);

    // GLR constructors
    GLR(const int shape);

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
    GLR operator*(const GLR& other) const;
    GLR& operator*=(const GLR& other);
    GLR inverse() const;
};

}

#endif
