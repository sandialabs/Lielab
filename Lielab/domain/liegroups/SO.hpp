#ifndef LIELAB_DOMAIN_SO_HPP
#define LIELAB_DOMAIN_SO_HPP

#include "../VirtualManifolds.hpp"
#include "SU.hpp"

#include "Lielab/domain/liealgebras/so.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

class SU;

class SO : public VirtualLieGroup<double>
{
    public:
    // Manifold typing
    using point_t = Eigen::MatrixXd;
    
    // Manifold storage
    point_t point;

    // Manifold constructors
    SO();
    // ~SO();

    // Lie group constructors
    SO(const matrix_t& other);
    static SO identity(const int shape);
    static SO project(const matrix_t& other);

    // SO constructors
    SO(const int shape);
    static SO from_eulerangles_body123(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_body231(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_body312(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_body132(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_body213(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_body321(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_body121(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_body131(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_body212(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_body232(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_body313(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_body323(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_space123(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_space231(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_space312(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_space132(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_space213(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_space321(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_space121(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_space131(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_space212(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_space232(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_space313(const double theta1, const double theta2, const double theta3);
    static SO from_eulerangles_space323(const double theta1, const double theta2, const double theta3);
    static SO from_quaternion(const double e0, const double e1, const double e2, const double e3);
    static SO from_rodriguesvector(const double g1, const double g2, const double g3);
    static SO from_SU2(const SU& quaternion);

    // Manifold information
    std::string to_string() const override;
    int get_dimension() const override;
    int get_size() const override;

    // Lie group information
    bool is_abelian() const override;
    int get_shape() const override;

    // Manifold IO
    point_t get_point() const;
    Eigen::VectorXd serialize() const override;
    void unserialize(const Eigen::VectorXd& serialized) override;
    void unserialize(std::initializer_list<double> serialized) override;

    // Lie group IO
    matrix_t get_matrix() const;
    field_t operator()(const int index1, const int index2) const;

    // SO IO
    std::array<double, 3> to_eulerangles_body123() const;
    std::array<double, 3> to_eulerangles_body231() const;
    std::array<double, 3> to_eulerangles_body312() const;
    std::array<double, 3> to_eulerangles_body132() const;
    std::array<double, 3> to_eulerangles_body213() const;
    std::array<double, 3> to_eulerangles_body321() const;
    std::array<double, 3> to_eulerangles_body121() const;
    std::array<double, 3> to_eulerangles_body131() const;
    std::array<double, 3> to_eulerangles_body212() const;
    std::array<double, 3> to_eulerangles_body232() const;
    std::array<double, 3> to_eulerangles_body313() const;
    std::array<double, 3> to_eulerangles_body323() const;
    std::array<double, 3> to_eulerangles_space123() const;
    std::array<double, 3> to_eulerangles_space231() const;
    std::array<double, 3> to_eulerangles_space312() const;
    std::array<double, 3> to_eulerangles_space132() const;
    std::array<double, 3> to_eulerangles_space213() const;
    std::array<double, 3> to_eulerangles_space321() const;
    std::array<double, 3> to_eulerangles_space121() const;
    std::array<double, 3> to_eulerangles_space131() const;
    std::array<double, 3> to_eulerangles_space212() const;
    std::array<double, 3> to_eulerangles_space232() const;
    std::array<double, 3> to_eulerangles_space313() const;
    std::array<double, 3> to_eulerangles_space323() const;
    std::array<double, 4> to_quaternion() const;
    std::array<double, 3> to_gibbs() const;

    // Lie group math
    SO operator*(const SO& other) const;
    SO& operator*=(const SO& other);
    SO inverse() const;

    
    // TODO: Organize these functions
    so project_onto_tangent_space(const Eigen::MatrixXd& vector);
};

}

#endif
