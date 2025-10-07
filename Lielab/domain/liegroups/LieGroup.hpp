#ifndef LIELAB_DOMAIN_LIEGROUP_HPP
#define LIELAB_DOMAIN_LIEGROUP_HPP

#include <Eigen/Core>

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

    using field_t = Field;
    using matrix_t = Eigen::Matrix<Field, Eigen::Dynamic, Eigen::Dynamic>;
    using data_t = Eigen::Matrix<Field, Eigen::Dynamic, Eigen::Dynamic>;

    // Lie group base info
    // bool abelian = false;
    // virtual bool is_abelian() {return this->abelian;}
    // bool complex = true;
    // virtual bool is_complex() {return this->complex;}

    data_t data = data_t::Zero(0, 0);

    // std::string name = "LieAlgebra";
    // ptrdiff_t shape = 0;
    // virtual ptrdiff_t get_shape() {return this->shape;}

    // data_t data = data_t::Zero(0, 0);

    // virtual matrix_t get_matrix() {return this->data;}
    // virtual data_t get_data() {return this->data;}

    LieGroup() { }
    virtual Eigen::VectorXd serialize() const = 0;
    virtual void unserialize(const Eigen::VectorXd&) = 0;
    virtual size_t get_shape() const = 0;
    virtual size_t get_size() const = 0;
    // LieGroup(const matrix_t &other)
    // {
        // assert(other.rows() == other.cols());

        // this->data = other;
        // this->shape = other.rows();
    // }

    // virtual ~LieGroup() {}

    // virtual size_t get_dimension() const = 0;
};

}

#include "LieGroup.tpp"

#endif
