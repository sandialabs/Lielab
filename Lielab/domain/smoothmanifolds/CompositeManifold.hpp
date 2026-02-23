#ifndef LIELAB_DOMAIN_COMPOSITEMANIFOLD_HPP
#define LIELAB_DOMAIN_COMPOSITEMANIFOLD_HPP

#include "../VirtualManifolds.hpp"

#include "../liealgebras.hpp"
#include "../liegroups.hpp"
#include "Grassmannian.hpp"

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

typedef std::variant<CN, GLC, GLR, RN, SE, SO, SP, SU,
                     cn, glc, glr, rn, se, so, sp, su,
                     Grassmannian> CompositeManifoldTYPES;

class CompositeManifold : virtual public VirtualComposite<CompositeManifoldTYPES>
{
    public:
    // Manifold typing
    using field_t = std::complex<double>;
    using point_t = std::vector<CompositeManifoldTYPES>;

    // CompositeManifold typing
    using TYPES = CompositeManifoldTYPES;
    static constexpr size_t INDEX_CN  = 0;
    static constexpr size_t INDEX_GLC = 1;
    static constexpr size_t INDEX_GLR = 2;
    static constexpr size_t INDEX_RN  = 3;
    static constexpr size_t INDEX_SE  = 4;
    static constexpr size_t INDEX_SO  = 5;
    static constexpr size_t INDEX_SP  = 6;
    static constexpr size_t INDEX_SU  = 7;
    static constexpr size_t INDEX_cn  = 8;
    static constexpr size_t INDEX_glc = 9;
    static constexpr size_t INDEX_glr = 10;
    static constexpr size_t INDEX_rn  = 11;
    static constexpr size_t INDEX_se  = 12;
    static constexpr size_t INDEX_so  = 13;
    static constexpr size_t INDEX_sp  = 14;
    static constexpr size_t INDEX_su  = 15;
    static constexpr size_t INDEX_Grassmannian = 16;

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
    CompositeManifold();

    // CompositeManifold constructors
    CompositeManifold(const int n);
    CompositeManifold(std::initializer_list<TYPES> others);
    CompositeManifold(const std::vector<TYPES>& others);

    // Manifold information
    std::string to_string() const;
    int get_dimension() const;
    int get_size() const;

    // Composite information
    std::vector<TYPES>::iterator begin() override;
    std::vector<TYPES>::iterator end() override;
    std::vector<TYPES>::const_iterator begin() const override;
    std::vector<TYPES>::const_iterator end() const override;

    // CompositeManifold information
    std::vector<int> get_dimensions() const;
    std::vector<int> get_sizes() const;

    // Manifold IO
    point_t get_point() const;
    Eigen::VectorXd serialize() const;
    void unserialize(const Eigen::VectorXd& vec);
    void unserialize(std::initializer_list<double> vec);

    // CompositeManifold IO
    const dataproxy operator[](const int index) const;
    dataproxy operator[](const int index);
};

}

#endif
