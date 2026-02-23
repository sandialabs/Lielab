#ifndef LIELAB_DOMAIN_SE_HPP
#define LIELAB_DOMAIN_SE_HPP

#include "../VirtualManifolds.hpp"
#include "RN.hpp"
#include "SO.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include "tuple"

namespace Lielab::domain
{
class SE : public VirtualLieGroup<double>
{
    public:
    // Manifold typing
    using point_t = std::tuple<RN, SO>;

    // Manifold storage
    point_t point;

    // SE storage
    int _shape = 0;

    // Manifold constructors
    SE();
    // ~SE();

    // Lie group constructors
    SE(const matrix_t& matrix);
    static SE identity(const int shape);
    static SE project(const matrix_t& matrix);

    // SE constructors
    SE(const int shape);
    SE(const RN& RN_c, const SO& SO_c);

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
    SE operator*(const SE& other) const;
    SE& operator*=(const SE& other);
    SE inverse() const;
};

}

#endif
