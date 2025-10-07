#ifndef LIELAB_DOMAIN_SP_HPP
#define LIELAB_DOMAIN_SP_HPP

#include "LieGroup.hpp"
#include "GLR.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{
class SP : public GLR
{
    /*!
    * The SP class.
    */
    public:
    static constexpr bool abelian = false;

    std::string to_string() const;

    // Initialization methods

    SP();
    SP(const size_t shape);
    static SP from_shape(const size_t shape);

    template<typename OtherDerived>
    SP(const Eigen::MatrixBase<OtherDerived> & other);

    template<typename OtherDerived>
    SP & operator=(const Eigen::MatrixBase<OtherDerived> & other);

    // TODO: project

    size_t get_dimension() const;
    size_t get_shape() const;

    size_t get_size() const;

    Eigen::MatrixXd get_matrix() const;

    SP inverse() const;

    Eigen::VectorXd serialize() const override;

    void unserialize(const Eigen::VectorXd &vec) override;
    void unserialize(std::initializer_list<double> vec);

    double operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    SP operator*(const SP & other) const;

    SP & operator*=(const SP & other);

    friend std::ostream & operator<<(std::ostream& os, const SP & other);
};

}

#include "SP.tpp"

#endif
