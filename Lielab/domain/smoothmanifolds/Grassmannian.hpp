#ifndef LIELAB_DOMAIN_Grassmannian_HPP
#define LIELAB_DOMAIN_Grassmannian_HPP

#include "Lielab/domain/liealgebras/rn.hpp"

#include <Eigen/Core>

namespace Lielab::domain
{

class Grassmannian
{
    public:
    // Manifold typing
    using field_t = double;
    using point_t = Eigen::VectorXd;

    // Manifold storage
    point_t point;
    
    // Grassmannian storage
    Eigen::MatrixXd axes;
    Eigen::MatrixXd point_projector;

    // Manifold constructors
    Grassmannian();
    // ~Grassmannian();

    // Grassmannian constructors
    Grassmannian(const int k, const int n);
    Grassmannian(const Eigen::VectorXd& point, const Eigen::MatrixXd& axes);
    static Grassmannian project(const Eigen::VectorXd& point, const Eigen::MatrixXd& axes);

    // Manifold information
    std::string to_string() const;
    int get_dimension() const;
    int get_size() const;

    // Manifold IO
    point_t get_point() const;
    Eigen::VectorXd serialize() const;
    void unserialize(const Eigen::VectorXd& vec);
    void unserialize(std::initializer_list<double> vec);

    // Grassmannian IO
    const double& operator[](const int index) const;
    double& operator[](const int index);

    // Other misc methods
    Eigen::VectorXd project_point(const Eigen::VectorXd& other_point) const;
    Eigen::VectorXd project_vector_onto_tangent_space(const Eigen::VectorXd& vector) const;
    Eigen::VectorXd project_vector_onto_normal_space(const Eigen::VectorXd& vector) const;
    Eigen::VectorXd axes_intersection(const Eigen::VectorXd& other_point, const Eigen::MatrixXd& other_axes) const;
};

}

#endif
