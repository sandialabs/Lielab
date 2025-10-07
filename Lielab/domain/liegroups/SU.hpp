#ifndef LIELAB_DOMAIN_SU_HPP
#define LIELAB_DOMAIN_SU_HPP

#include "LieGroup.hpp"
#include "GLC.hpp"
#include "SO.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class SO;

class SU : public GLC
{
    /*!
    * The SU class.
    */
    public:
    static constexpr bool abelian = false;
    
    std::string to_string() const;

    // Initialization methods

    SU();
    SU(const size_t shape);
    static SU from_shape(const size_t shape);

    template<typename OtherDerived>
    SU(const Eigen::MatrixBase<OtherDerived> & other);

    template<typename OtherDerived>
    SU & operator=(const Eigen::MatrixBase<OtherDerived> & other);

    // TODO: Project

    size_t get_dimension() const;
    size_t get_shape() const;

    size_t get_size() const;

    Eigen::MatrixXcd get_matrix() const;

    SU inverse() const;

    Eigen::VectorXd serialize() const override;

    void unserialize(const Eigen::VectorXd &vec) override;
    void unserialize(std::initializer_list<double> vec);

    std::complex<double> operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    SU operator*(const SU & other) const;

    SU & operator*=(const SU & other);

    friend std::ostream & operator<<(std::ostream& os, const SU & other);

    /*
     * Additional static initializers. Not a part of the core Lie group, but are convenient.
     */
    template <typename T>
    static SU from_quaternion(const T e0, const T e1, const T e2, const T e3);

    static SU from_SO3(const SO & dcm);

    std::array<double, 4> to_quaternion() const;
};

}

#include "SU.tpp"

#endif
