#ifndef LIELAB_DOMAIN_COMPOSITEGROUP_HPP
#define LIELAB_DOMAIN_COMPOSITEGROUP_HPP

#include "../VirtualManifolds.hpp"

#include "CN.hpp"
#include "GLC.hpp"
#include "GLR.hpp"
#include "RN.hpp"
#include "SE.hpp"
#include "SO.hpp"
#include "SP.hpp"
#include "SU.hpp"

#include "Lielab/utils.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <complex>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <variant>

namespace Lielab::domain
{

typedef std::variant<CN, GLC, GLR, RN, SE, SO, SP, SU> CompositeGroupTYPES;

class CompositeGroup : public VirtualLieGroup<std::complex<double>>, virtual public VirtualComposite<CompositeGroupTYPES>
{
    public:
    // Manifold typing
    using point_t = std::vector<CompositeGroupTYPES>;

    // CompositeGroup typing
    using TYPES = CompositeGroupTYPES;
    static constexpr size_t INDEX_CN  = 0;
    static constexpr size_t INDEX_GLC = 1;
    static constexpr size_t INDEX_GLR = 2;
    static constexpr size_t INDEX_RN  = 3;
    static constexpr size_t INDEX_SE  = 4;
    static constexpr size_t INDEX_SO  = 5;
    static constexpr size_t INDEX_SP  = 6;
    static constexpr size_t INDEX_SU  = 7;

    struct dataproxy
    {
        TYPES& var;

        template <class T>
        operator T() const
        {
            return std::get<T>(var);
        }

        template <class T>
        dataproxy& operator=(const T& other)
        {
            var = other;
            return *this;
        }
    };    
    
    // Manifold storage
    point_t point;

    // Manifold constructors
    CompositeGroup();

    // Lie group constructors
    CompositeGroup(const Eigen::MatrixXcd& matrix);
    static CompositeGroup identity(const int shape);
    static CompositeGroup project(const matrix_t& matrix);

    // CompositeGroup constructors
    CompositeGroup(const int n);
    CompositeGroup(std::initializer_list<TYPES> others);
    CompositeGroup(const std::vector<TYPES>& others);

    // Manifold information
    std::string to_string() const override;
    int get_dimension() const override;
    int get_size() const override;

    // Lie group information
    bool is_abelian() const override;
    int get_shape() const override;

    // Composite information
    std::vector<TYPES>::iterator begin() override;
    std::vector<TYPES>::iterator end() override;
    std::vector<TYPES>::const_iterator begin() const override;
    std::vector<TYPES>::const_iterator end() const override;

    // CompositeGroup information
    std::vector<int> get_dimensions() const;
    std::vector<int> get_sizes() const;
    std::vector<int> get_shapes() const;

    // Manifold IO
    point_t get_point() const;
    Eigen::VectorXd serialize() const override;
    void unserialize(const Eigen::VectorXd& serialized) override;
    void unserialize(std::initializer_list<double> serialized) override;

    // Lie group IO
    CompositeGroup::matrix_t get_matrix() const;
    field_t operator()(const int index1, const int index2) const;

    // CompositeGroup IO
    // TODO: get_matrices() (plural)
    const dataproxy operator[](const int index) const;
    dataproxy operator[](const int index);

    // Lie Group math ops
    CompositeGroup operator*(const CompositeGroup& other) const;
    CompositeGroup& operator*=(const CompositeGroup& other);
    CompositeGroup inverse() const;
};

}

#endif
