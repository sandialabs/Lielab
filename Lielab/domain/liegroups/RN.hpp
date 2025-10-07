#ifndef LIELAB_DOMAIN_RN_HPP
#define LIELAB_DOMAIN_RN_HPP

#include "LieGroup.hpp"
#include "GLR.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{
class RN : public GLR
{
    /*!
    * The RN class.
    */
    public:
    static constexpr bool abelian = true;

    using data_t = Eigen::VectorXd;
    data_t data;

    std::string to_string() const;

    // Initialization methods

    RN();
    RN(const size_t n);
    static RN from_shape(const size_t shape);

    template<typename OtherDerived>
    RN(const Eigen::MatrixBase<OtherDerived> & other);

    template<typename OtherDerived>
    RN & operator=(const Eigen::MatrixBase<OtherDerived> & other);

    static Eigen::MatrixXd project(const Eigen::MatrixXd & other);

    size_t get_dimension() const;
    size_t get_shape() const;

    size_t get_size() const;

    Eigen::MatrixXd get_matrix() const;

    RN inverse() const;

    // Data representation

    Eigen::VectorXd serialize() const override;

    void unserialize(const Eigen::VectorXd &vec) override;
    void unserialize(std::initializer_list<double> vec);

    double operator()(const ptrdiff_t index) const;
    double operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    RN operator*(const RN & other) const;

    RN & operator*=(const RN & other);

    static RN from_vector(const Eigen::VectorXd & other);
    static RN from_vector(const std::initializer_list<double> other);

    friend std::ostream & operator<<(std::ostream& os, const RN & other);
};

}

#include "RN.tpp"

#endif
