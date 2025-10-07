#ifndef LIELAB_DOMAIN_LIEALGEBRAS_LIEALGEBRA_HPP
#define LIELAB_DOMAIN_LIEALGEBRAS_LIEALGEBRA_HPP

#include <Eigen/Core>

#include <string>

namespace Lielab::domain
{

/*!
 * Lie algebra base class.
 */
template <typename Field>
class LieAlgebra
{
    public:

    // Storage and typing
    using field_t = Field;
    using matrix_t = Eigen::Matrix<Field, Eigen::Dynamic, Eigen::Dynamic>;
    using data_t = Eigen::Matrix<Field, Eigen::Dynamic, Eigen::Dynamic>;
    data_t data = data_t::Zero(0, 0);

    // Lie Algebra class information
    virtual bool is_abelian() const;
    virtual bool is_complex() const;
    virtual std::string to_string() const;

    // Constructors and destructors
    LieAlgebra();
    LieAlgebra(const size_t n);
    LieAlgebra(const LieAlgebra::matrix_t& other);
    ~LieAlgebra();

    // Object information
    virtual size_t get_shape() const;
    virtual size_t get_dimension() const = 0;

    // Object IO and data manipulation
    virtual data_t get_data() const;
    virtual matrix_t get_matrix() const = 0;
    virtual Eigen::VectorXd get_vector() const = 0;
    virtual void set_vector(const Eigen::VectorXd &vec) = 0;

    // operator() and []'s here

    // Lie Algebra math ops
};

}

#include "LieAlgebra.tpp"

#endif
