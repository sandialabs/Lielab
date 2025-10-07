#ifndef LIELAB_DOMAIN_GL_HPP
#define LIELAB_DOMAIN_GL_HPP

#include "LieGroup.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class GLR : public LieGroup<double>
{
    /*!
    * The GLR class.
    */
    public:
    static constexpr bool abelian = false;
    size_t _shape = 0;

    std::string to_string() const;

    // Initialization methods

    GLR();

    GLR(const size_t shape);
    static GLR from_shape(const size_t shape);

    template<typename OtherDerived>
    GLR(const Eigen::MatrixBase<OtherDerived> & other);

    template<typename OtherDerived>
    GLR & operator=(const Eigen::MatrixBase<OtherDerived> & other);

    static Eigen::MatrixXd project(const Eigen::MatrixXd & other);

    size_t get_dimension() const;
    size_t get_shape() const;

    size_t get_size() const;

    Eigen::MatrixXd get_matrix() const;

    GLR inverse() const;

    // Data representation

    Eigen::VectorXd serialize() const override;

    void unserialize(const Eigen::VectorXd &vec) override;
    void unserialize(std::initializer_list<double> vec);

    double operator()(const ptrdiff_t index1, const ptrdiff_t index2) const;

    GLR operator*(const GLR & other) const;

    GLR & operator*=(const GLR & other);

    friend std::ostream & operator<<(std::ostream& os, const GLR & other);
};

}

#include "GLR.tpp"

#endif
