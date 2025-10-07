#ifndef LIELAB_DOMAIN_COMPOSITEGROUP_HPP
#define LIELAB_DOMAIN_COMPOSITEGROUP_HPP

#include "CN.hpp"
#include "GLR.hpp"
#include "GLC.hpp"
#include "RN.hpp"
#include "SE.hpp"
#include "SO.hpp"
#include "SP.hpp"
#include "SU.hpp"

#include "Lielab/domain/liealgebras.hpp"
#include "Lielab/utils.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <variant>

namespace Lielab::domain
{

class CompositeGroup : public LieGroup<std::complex<double>>
{
    public:

    // Storage and typing
    static constexpr size_t INDEX_CN  = 0;
    static constexpr size_t INDEX_GLR = 1;
    static constexpr size_t INDEX_GLC = 2;
    static constexpr size_t INDEX_RN  = 3;
    static constexpr size_t INDEX_SE  = 4;
    static constexpr size_t INDEX_SO  = 5;
    static constexpr size_t INDEX_SP  = 6;
    static constexpr size_t INDEX_SU  = 7;

    typedef std::variant<Lielab::domain::CN,
                         Lielab::domain::GLR,
                         Lielab::domain::GLC,
                         Lielab::domain::RN,
                         Lielab::domain::SE,
                         Lielab::domain::SO,
                         Lielab::domain::SP,
                         Lielab::domain::SU> TYPES;

    std::vector<TYPES> space;

    // Lie Group class information
    std::string to_string() const;

    // Constructors and destructors
    CompositeGroup();
    CompositeGroup(const size_t n);
    static CompositeGroup from_shape(const size_t shape);
    CompositeGroup(std::initializer_list<TYPES> others);
    CompositeGroup(const std::vector<TYPES>& others);

    // Object information
    size_t get_dimension() const;
    std::vector<size_t> get_dimensions() const;
    size_t get_shape() const;
    std::vector<size_t> get_shapes() const;
    size_t get_size() const;
    std::vector<size_t> get_sizes() const;

    // Object IO and data manipulation
    Eigen::MatrixXcd get_matrix() const;
    // TODO: get_matrices() (plural)
    Eigen::VectorXd serialize() const;
    void unserialize(const Eigen::VectorXd& vec);
    void unserialize(std::initializer_list<double> vec);

    // double operator()(const ptrdiff_t index) const;
    std::complex<double> operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;
    TYPES operator[](const ptrdiff_t index) const;

    // Lie Group math ops
    CompositeGroup operator*(const CompositeGroup & other) const;
    CompositeGroup & operator*=(const CompositeGroup & other);
    CompositeGroup inverse() const;

    friend std::ostream & operator<<(std::ostream & os, const CompositeGroup & other);
};

}

#endif
