#ifndef LIELAB_DOMAIN_SE_HPP
#define LIELAB_DOMAIN_SE_HPP

#include "LieGroup.hpp"
#include "GLR.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{
class SE : public GLR
{
    /*!
    * The SE class.
    */
    public:
    static constexpr bool abelian = false;

    std::string to_string() const;

    SE();
    SE(const size_t shape);
    static SE from_shape(const size_t shape);

    template<typename OtherDerived>
    SE(const Eigen::MatrixBase<OtherDerived> & other);

    template<typename OtherDerived>
    SE & operator=(const Eigen::MatrixBase<OtherDerived> & other);

    size_t get_dimension() const;
    size_t get_shape() const;

    size_t get_size() const;

    Eigen::MatrixXd get_matrix() const;

    SE inverse() const;

    Eigen::VectorXd serialize() const override;

    void unserialize(const Eigen::VectorXd &vec) override;
    void unserialize(std::initializer_list<double> vec);

    double operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    SE operator*(const SE & other) const;

    SE & operator*=(const SE & other);

    friend std::ostream & operator<<(std::ostream & os, const SE & other);
};

}

#include "SE.tpp"

#endif
