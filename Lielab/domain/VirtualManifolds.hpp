#ifndef LIELAB_DOMAIN_VIRTUALMANIFOLDS_HPP
#define LIELAB_DOMAIN_VIRTUALMANIFOLDS_HPP

#include <Eigen/Core>

#include <concepts>
#include <string>
#include <vector>

namespace Lielab::domain
{

class VirtualManifold
{
    public:
    // Manifold typing
    // using point_t = ?

    // Manifold storage
    // point_t point;

    // Manifold constructors
    // virtual VirtualManifold() = default;
    // virtual ~VirtualManifold() = default;

    // Manifold information
    virtual std::string to_string() const = 0;
    virtual int get_dimension() const = 0;
    virtual int get_size() const = 0;

    // Manifold IO
    // virtual point_t get_point const() = 0;
    virtual Eigen::VectorXd serialize() const = 0;
    virtual void unserialize(const Eigen::VectorXd& serialized) = 0;
    virtual void unserialize(std::initializer_list<double> serialized) = 0;
};

template <typename Field>
class VirtualLieAlgebra : virtual public VirtualManifold
{
    public:
    // Lie algebra typing
    using field_t = Field;
    using matrix_t = Eigen::Matrix<Field, Eigen::Dynamic, Eigen::Dynamic>;

    // Lie algebra information
    virtual bool is_abelian() const = 0;

    // Lie algebra constructors
    // virtual VirtualLieAlgebra(const matrix_t& other) = 0;
    // virtual static VirtualLieAlgebra basis(const int index, const int shape) = 0;
    // virtual static VirtualLieAlgebra zero(const int shape) = 0;
    // virtual static VirtualLieAlgebra from_vector(const Eigen::VectorXd& other);
    // virtual static VirtualLieAlgebra from_vector(std::initializer_list<double> other);
    // virtual static VirtualLieAlgebra project(const matrix_t& other);

    // Lie algebra information
    virtual int get_shape() const = 0;

    // Lie algebra IO
    // virtual matrix_t get_matrix() const = 0;
    virtual Eigen::VectorXd get_vector() const = 0;
    virtual void set_vector(const Eigen::VectorXd& vec) = 0;
    virtual void set_vector(std::initializer_list<double> vec) = 0;
    // double operator()(const int index) const = 0;
    // field_t operator()(const int index1, const int index2) const = 0;

    // Lie algebra math
    // virtual VirtualLieAlgebra operator+(const VirtualLieAlgebra& other) const = 0;
    // virtual VirtualLieAlgebra& operator+=(const VirtualLieAlgebra& other) = 0;
    // virtual VirtualLieAlgebra operator-(const VirtualLieAlgebra& other) const = 0;
    // virtual VirtualLieAlgebra& operator-=(const VirtualLieAlgebra& other) = 0;
    // virtual VirtualLieAlgebra operator-() const = 0;
    // virtual VirtualLieAlgebra operator*(const double other) const = 0;
    // virtual friend VirtualLieAlgebra operator*(const double other, const VirtualLieAlgebra& rhs) = 0;
    // virtual VirtualLieAlgebra operator*(const double other) const = 0;
    // virtual friend VirtualLieAlgebra operator*(const double other, const VirtualLieAlgebra& rhs) = 0;
    // virtual VirtualLieAlgebra& operator*=(const double other) = 0;
    // virtual VirtualLieAlgebra operator/(const double other) const = 0;
    // virtual VirtualLieAlgebra& operator/=(const double other) = 0;
};

template <typename Field>
class VirtualLieGroup : virtual public VirtualManifold
{
    public:
    // Lie group typing
    using field_t = Field;
    using matrix_t = Eigen::Matrix<Field, Eigen::Dynamic, Eigen::Dynamic>;

    // Lie group information
    virtual bool is_abelian() const = 0;

    // Lie group constructors
    // virtual VirtualLieGroup(const matrix_t& other) = 0;
    // virtual static VirtualLieGroup identity(const int shape) = 0;
    // virtual static VirtualLieGroup project(const matrix_t& other);

    // Lie group information
    virtual int get_shape() const = 0;

    // Lie group IO
    // virtual matrix_t get_matrix() const = 0;
    // field_t operator()(const int index1, const int index2) const = 0;

    // Lie group math
    // VirtualLieGroup operator*(const VirtualLieGroup& other) const;
    // VirtualLieGroup& operator*=(const VirtualLieGroup& other);
    // VirtualLieGroup inverse() const;

};

template <typename Types>
class VirtualComposite
{
    public:

    virtual std::vector<Types>::iterator begin() = 0;
    virtual std::vector<Types>::iterator end() = 0;
    virtual std::vector<Types>::const_iterator begin() const = 0;
    virtual std::vector<Types>::const_iterator end() const = 0;
};

}

#endif
