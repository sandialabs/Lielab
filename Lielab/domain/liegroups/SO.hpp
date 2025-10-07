#ifndef LIELAB_DOMAIN_SO_HPP
#define LIELAB_DOMAIN_SO_HPP

#include "GLR.hpp"
#include "SU.hpp"
#include "LieGroup.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class SU;

class SO : public GLR
{
    /*!
    * The SO class.
    */
    public:
    static constexpr bool abelian = false;
    
    std::string to_string() const;
    SO();
    SO(const size_t shape);
    static SO from_shape(const size_t shape);

    template<typename OtherDerived>
    SO(const Eigen::MatrixBase<OtherDerived> & other);

    template<typename OtherDerived>
    SO & operator=(const Eigen::MatrixBase<OtherDerived> & other);

    static Eigen::MatrixXd project(const Eigen::MatrixXd & other);

    size_t get_dimension() const;
    size_t get_shape() const;

    size_t get_size() const;

    Eigen::MatrixXd get_matrix() const;

    SO inverse() const;

    Eigen::VectorXd serialize() const override;

    void unserialize(const Eigen::VectorXd &vec) override;
    void unserialize(std::initializer_list<double> vec);

    double operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    SO operator*(const SO & other) const;

    SO & operator*=(const SO & other);

    friend std::ostream & operator<<(std::ostream & os, const SO & other);

    /*
     * Additional static initializers. Not a part of the core Lie group, but are convenient.
     */
    template <typename T>
    static SO from_eulerangles_body123(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_body231(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_body312(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_body132(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_body213(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_body321(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_body121(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_body131(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_body212(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_body232(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_body313(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_body323(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_space123(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_space231(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_space312(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_space132(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_space213(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_space321(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_space121(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_space131(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_space212(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_space232(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_space313(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_eulerangles_space323(const T theta1, const T theta2, const T theta3);

    template <typename T>
    static SO from_quaternion(const T e0, const T e1, const T e2, const T e3);

    template <typename T>
    static SO from_rodriguesvector(const T g1, const T g2, const T g3);

    static SO from_SU2(const SU & quaternion);

    /*
    * Additional outputs. Not a part of the core Lie group, but are convenient.
    */
    
    template <typename T>
    std::array<T, 4> to_quaternion() const;

    template <typename T>
    std::array<T, 3> to_gibbs() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_body123() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_body231() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_body312() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_body132() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_body213() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_body321() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_body121() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_body131() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_body212() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_body232() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_body313() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_body323() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_space123() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_space231() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_space312() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_space132() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_space213() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_space321() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_space121() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_space131() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_space212() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_space232() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_space313() const;

    template <typename T>
    std::array<T, 3> to_eulerangles_space323() const;
};

}

#include "SO.tpp"

#endif
