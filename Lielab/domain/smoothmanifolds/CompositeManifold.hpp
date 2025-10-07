#ifndef LIELAB_DOMAIN_COMPOSITEMANIFOLD_HPP
#define LIELAB_DOMAIN_COMPOSITEMANIFOLD_HPP

#include "../liealgebras.hpp"
#include "../liegroups.hpp"

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

class CompositeManifold
{
    public:

    static constexpr size_t INDEX_CN  = 0;
    static constexpr size_t INDEX_GLR = 1;
    static constexpr size_t INDEX_GLC = 2;
    static constexpr size_t INDEX_RN  = 3;
    static constexpr size_t INDEX_SE  = 4;
    static constexpr size_t INDEX_SO  = 5;
    static constexpr size_t INDEX_SP  = 6;
    static constexpr size_t INDEX_SU  = 7;
    static constexpr size_t INDEX_cn  = 8;
    static constexpr size_t INDEX_glr = 9;
    static constexpr size_t INDEX_glc = 10;
    static constexpr size_t INDEX_rn  = 11;
    static constexpr size_t INDEX_se  = 12;
    static constexpr size_t INDEX_so  = 13;
    static constexpr size_t INDEX_sp  = 14;
    static constexpr size_t INDEX_su  = 15;

    typedef std::variant<Lielab::domain::CN,
                         Lielab::domain::GLR,
                         Lielab::domain::GLC,
                         Lielab::domain::RN,
                         Lielab::domain::SE,
                         Lielab::domain::SO,
                         Lielab::domain::SP,
                         Lielab::domain::SU,
                         Lielab::domain::cn,
                         Lielab::domain::glr,
                         Lielab::domain::glc,
                         Lielab::domain::rn,
                         Lielab::domain::se,
                         Lielab::domain::so,
                         Lielab::domain::sp,
                         Lielab::domain::su> TYPES;

    std::vector<TYPES> space;

    std::string to_string() const;

    CompositeManifold();
    CompositeManifold(const size_t n);
    static CompositeManifold from_shape(const size_t shape);
    CompositeManifold(std::initializer_list<TYPES> others);
    CompositeManifold(const std::vector<TYPES>& others);

    size_t get_dimension() const;
    std::vector<size_t> get_dimensions() const;
    size_t get_size() const;
    std::vector<size_t> get_sizes() const;
    size_t get_shape() const;
    std::vector<size_t> get_shapes() const;

    Eigen::VectorXd serialize() const;
    void unserialize(const Eigen::VectorXd& vec);
    void unserialize(std::initializer_list<double> vec);
    Eigen::MatrixXcd get_matrix() const;

    std::complex<double> operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;
    TYPES operator[](const ptrdiff_t index) const;
    friend std::ostream & operator<<(std::ostream & os, const CompositeManifold & other);
};

}

#endif
