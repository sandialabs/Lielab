#ifndef LIELAB_INTEGRATE_IVPSETTINGS_HPP
#define LIELAB_INTEGRATE_IVPSETTINGS_HPP

#include "Lielab/domain.hpp"
#include "Lielab/functions.hpp"

#include <Eigen/Core>

#include <functional>

namespace Lielab::integrate
{

enum class IVPMethod
{
    Undefined = 0,
    RungeKutta = 1,
    MuntheKaas = 2,
};

struct IVPOptions
{
    public:
    
    // IVP Meta-options
    IVPMethod method = IVPMethod::Undefined;

    // Multi-method options
    double dt_initial = std::numeric_limits<double>::quiet_NaN();
    double dt_min = 1e-4;
    double dt_max = 10.0;
    bool variable_time_step = true;

    double reltol = 1e-4;
    double abstol = 1e-8;

    double small = 0.2;
    double large = 10.0;
    double pessimist = 0.9;

    // Munthe-Kaas specific options
    // bool rebase_every_step = true; # TODO:
};

using EuclideanIVP_event_t = std::function<double(const double, const Eigen::VectorXd&)>;
using EuclideanIVP_vectorfield_t = std::function<Eigen::VectorXd(const double, const Eigen::VectorXd&)>;

struct EuclideanIVPSystem
{
    public:

    EuclideanIVP_event_t event = [](const double t, const Eigen::VectorXd& y){return std::numeric_limits<double>::signaling_NaN();};
    EuclideanIVP_vectorfield_t vectorfield = [](const double t, const Eigen::VectorXd& y){return Eigen::VectorXd::Ones(1)*std::numeric_limits<double>::signaling_NaN();};
};

using HomogeneousIVP_action_t = std::function<Lielab::domain::CompositeManifold(const Lielab::domain::CompositeGroup&, const Lielab::domain::CompositeManifold&)>;
using HomogeneousIVP_connection_t = std::function<Lielab::domain::CompositeAlgebra(const Lielab::domain::CompositeAlgebra&, const Lielab::domain::CompositeAlgebra&)>;
using HomogeneousIVP_coordinates_t = std::function<Lielab::domain::CompositeGroup(const Lielab::domain::CompositeAlgebra&)>;
using HomogeneousIVP_event_t = std::function<double(const double, const Lielab::domain::CompositeManifold&)>;
using HomogeneousIVP_vectorfield_t = std::function<Lielab::domain::CompositeAlgebra(const double, const Lielab::domain::CompositeManifold&)>;

struct HomogeneousIVPSystem
{
    public:

    HomogeneousIVP_action_t action = Lielab::functions::left_Lie_group_action;
    HomogeneousIVP_connection_t connection = [](const Lielab::domain::CompositeAlgebra& a, const Lielab::domain::CompositeAlgebra& b){return Lielab::functions::dexpinv(a, b, 5);};
    HomogeneousIVP_coordinates_t coordinates = Lielab::functions::exp<Lielab::domain::CompositeAlgebra>;
    HomogeneousIVP_event_t event = [](const double t, const Lielab::domain::CompositeManifold& y){return std::numeric_limits<double>::signaling_NaN();};
    HomogeneousIVP_vectorfield_t vectorfield = [](const double t, const Lielab::domain::CompositeManifold& y){return Lielab::domain::CompositeAlgebra({Lielab::domain::rn::from_vector({std::numeric_limits<double>::signaling_NaN()})});};
};

}

#endif
