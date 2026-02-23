#ifndef LIELAB_DOMAIN_SP_HPP
#define LIELAB_DOMAIN_SP_HPP

#include "../VirtualManifolds.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{
class SP : public VirtualLieGroup<double>
{
    public:
    // Manifold typing
    using point_t = Eigen::MatrixXd;
    
    // Manifold storage
    point_t point;

    // Manifold constructors
    SP();
    // ~SP();

    // Lie group constructors
    SP(const matrix_t& matrix);
    static SP identity(const int shape);
    // TODO: project

    // SP constructors
    SP(const int shape);

    // Manifold information
    std::string to_string() const override;
    int get_dimension() const override;
    int get_size() const override;

    // Lie group informatioin
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
    SP operator*(const SP& other) const;
    SP& operator*=(const SP& other);
    SP inverse() const;
};

}

#endif
