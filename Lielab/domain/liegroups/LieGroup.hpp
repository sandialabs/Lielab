#ifndef LIELAB_DOMAIN_LIEGROUP_HPP
#define LIELAB_DOMAIN_LIEGROUP_HPP

#include <Eigen/Core>

#include <cassert>
#include <string>

namespace Lielab::domain
{

/*!
 * Lie group base class.
 */
template <typename Field>
class LieGroup
{
    public:

    // Storage and typing
    using field_t = Field;
    using matrix_t = Eigen::Matrix<Field, Eigen::Dynamic, Eigen::Dynamic>;
    using data_t = Eigen::Matrix<Field, Eigen::Dynamic, Eigen::Dynamic>;
    data_t data = data_t::Zero(0, 0);

    // Lie Group class information
    virtual bool is_abelian() const;
    virtual bool is_complex() const;
    virtual std::string to_string() const;

    // Constructors and destructors
    LieGroup();
    LieGroup(const size_t n);
    LieGroup(const LieGroup::matrix_t& other);
    ~LieGroup();
    
    // Object information
    virtual size_t get_dimension() const = 0;
    virtual size_t get_shape() const = 0;
    virtual size_t get_size() const = 0;

    // Object IO and data manipulation
    // virtual matrix_t get_matrix() = 0;
    virtual Eigen::VectorXd serialize() const = 0;
    virtual void unserialize(const Eigen::VectorXd& vec) = 0;
    
    // operator() and []'s here

    // Lie Group math ops
};

}

#include "LieGroup.tpp"

#endif
