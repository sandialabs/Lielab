#ifndef LIELAB_DOMAIN_LIEALGEBRAS_COMPOSITEALGEBRA_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_COMPOSITEALGEBRA_HPP

#include "../VirtualManifolds.hpp"

#include "cn.hpp"
#include "glc.hpp"
#include "glr.hpp"
#include "rn.hpp"
#include "se.hpp"
#include "so.hpp"
#include "sp.hpp"
#include "su.hpp"

#include "Lielab/utils.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <complex>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <memory>
#include <variant>

namespace Lielab::domain
{

typedef std::variant<cn, glc, glr, rn, se, so, sp, su> CompositeAlgebraTYPES;

class CompositeAlgebra : public VirtualLieAlgebra<std::complex<double>>, virtual public VirtualComposite<CompositeAlgebraTYPES>
{
    public:
    // Manifold typing
    using point_t = std::vector<CompositeAlgebraTYPES>;

    // CompositeAlgebra typing
    using TYPES = CompositeAlgebraTYPES;
    static constexpr size_t INDEX_cn  = 0;
    static constexpr size_t INDEX_glc = 1;
    static constexpr size_t INDEX_glr = 2;
    static constexpr size_t INDEX_rn  = 3;
    static constexpr size_t INDEX_se  = 4;
    static constexpr size_t INDEX_so  = 5;
    static constexpr size_t INDEX_sp  = 6;
    static constexpr size_t INDEX_su  = 7;

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
    CompositeAlgebra();

    // Lie algebra constructors
    CompositeAlgebra(const matrix_t& other);
    static CompositeAlgebra basis(const int index, const int shape);
    static CompositeAlgebra zero(const int shape);
    // from_vector() TODO:
    // from_vector() TODO:
    static CompositeAlgebra project(const matrix_t& matrix);

    // CompositeAlgebra constructors
    CompositeAlgebra(const int n);
    CompositeAlgebra(std::initializer_list<TYPES> others);
    CompositeAlgebra(const std::vector<TYPES>& others);

    // Manifold information
    std::string to_string() const override;
    int get_dimension() const override;
    int get_size() const override;

    // Lie algebra information
    bool is_abelian() const override;
    int get_shape() const override;

    // Composite information
    std::vector<TYPES>::iterator begin() override;
    std::vector<TYPES>::iterator end() override;
    std::vector<TYPES>::const_iterator begin() const override;
    std::vector<TYPES>::const_iterator end() const override;

    // CompositeAlgebra information
    std::vector<int> get_dimensions() const;
    std::vector<int> get_sizes() const;
    std::vector<int> get_shapes() const;
    
    // Manifold IO
    point_t get_point() const;
    Eigen::VectorXd serialize() const override;
    void unserialize(const Eigen::VectorXd& serialized) override;
    void unserialize(std::initializer_list<double> serialized) override;

    // Lie algebra IO
    matrix_t get_matrix() const;
    Eigen::VectorXd get_vector() const override;
    void set_vector(const Eigen::VectorXd& vec) override;
    void set_vector(std::initializer_list<double> vec) override;
    double operator()(const int index) const;
    field_t operator()(const int index1, const int index2) const;

    // CompositeAlgebra IO
    // TODO: std::vector<Eigen::MatrixBase> get_matrices() (plural)
    std::vector<Eigen::VectorXd> get_vectors() const;
    // TODO: set_vectors() (plural)
    const dataproxy operator[](const int index) const;
    dataproxy operator[](const int index);

    // Lie Algebra math ops
    CompositeAlgebra operator+(const CompositeAlgebra& other) const;
    CompositeAlgebra& operator+=(const CompositeAlgebra& other);
    CompositeAlgebra operator-(const CompositeAlgebra& other) const;
    CompositeAlgebra& operator-=(const CompositeAlgebra& other);
    CompositeAlgebra operator-() const;
    CompositeAlgebra operator*(const double other) const;
    friend CompositeAlgebra operator*(const double other, const CompositeAlgebra& rhs);
    CompositeAlgebra& operator*=(const double other);
    CompositeAlgebra operator/(const double other) const;
    CompositeAlgebra& operator/=(const double other);
};

}

#endif
