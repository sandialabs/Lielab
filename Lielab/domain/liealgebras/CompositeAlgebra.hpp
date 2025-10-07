#ifndef LIELAB_DOMAIN_LIEALGEBRAS_COMPOSITEALGEBRA_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_COMPOSITEALGEBRA_HPP

#include "cn.hpp"
#include "glr.hpp"
#include "glc.hpp"
#include "rn.hpp"
#include "se.hpp"
#include "so.hpp"
#include "sp.hpp"
#include "su.hpp"

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

class CompositeAlgebra : public LieAlgebra<std::complex<double>>
{
    public:

    // Storage and typing
    static constexpr size_t INDEX_cn  = 0;
    static constexpr size_t INDEX_glr = 1;
    static constexpr size_t INDEX_glc = 2;
    static constexpr size_t INDEX_rn  = 3;
    static constexpr size_t INDEX_se  = 4;
    static constexpr size_t INDEX_so  = 5;
    static constexpr size_t INDEX_sp  = 6;
    static constexpr size_t INDEX_su  = 7;

    typedef std::variant<Lielab::domain::cn,
                         Lielab::domain::glr,
                         Lielab::domain::glc,
                         Lielab::domain::rn,
                         Lielab::domain::se,
                         Lielab::domain::so,
                         Lielab::domain::sp,
                         Lielab::domain::su> TYPES;

    std::vector<TYPES> space;

    // Lie Algebra class information
    std::string to_string() const override;

    // Constructors and destructors
    CompositeAlgebra();
    CompositeAlgebra(const size_t n);
    static CompositeAlgebra basis(const ptrdiff_t i, const size_t n);
    static CompositeAlgebra from_shape(const size_t shape);
    CompositeAlgebra(std::initializer_list<TYPES> others);
    CompositeAlgebra(const std::vector<TYPES>& others);

    // Object information
    size_t get_dimension() const;
    std::vector<size_t> get_dimensions() const;
    size_t get_shape() const;
    std::vector<size_t> get_shapes() const;

    // Object IO and data manipulation
    CompositeAlgebra::matrix_t get_matrix() const;
    // TODO: std::vector<Eigen::MatrixBase> get_matrices() (plural)
    Eigen::VectorXd get_vector() const;
    std::vector<Eigen::VectorXd> get_vectors() const;
    void set_vector(const Eigen::VectorXd& vec);
    void set_vector(std::initializer_list<double> vec);
    // TODO: set_vectors() (plural)

    double operator()(const ptrdiff_t index) const;
    std::complex<double> operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;
    TYPES operator[](const ptrdiff_t index) const;

    // Lie Algebra math ops
    CompositeAlgebra operator+(const CompositeAlgebra & other) const;
    CompositeAlgebra & operator+=(const CompositeAlgebra & other);
    CompositeAlgebra operator-(const CompositeAlgebra & other) const;
    CompositeAlgebra & operator-=(const CompositeAlgebra & other);
    CompositeAlgebra operator-() const;
    CompositeAlgebra operator*(const double other) const;
    friend CompositeAlgebra operator*(const double other, const CompositeAlgebra & rhs);
    CompositeAlgebra & operator*=(const double other);
    CompositeAlgebra operator/(const double other) const;
    CompositeAlgebra & operator/=(const double other);

    friend std::ostream & operator<<(std::ostream & os, const CompositeAlgebra & other);
};

// inline std::ostream & operator<<(std::ostream & os, const CompositeAlgebra & other)
// {
//     /*!
//     * Overloads the "<<" stream insertion operator.
//     */
//     // TODO: Get rid of inline
//     os << other.to_string();
//     return os;
// }

// inline CompositeAlgebra operator*(const double other, const CompositeAlgebra & rhs)
// {
//     // TODO: Get rid of inline
//     return rhs*other;
// }

}

#endif
