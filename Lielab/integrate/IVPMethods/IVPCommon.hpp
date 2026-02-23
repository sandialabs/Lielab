#ifndef LIELAB_INTEGRATE_IVPCOMMON_HPP
#define LIELAB_INTEGRATE_IVPCOMMON_HPP

#include "Coefficients.hpp"

#include "Lielab/domain.hpp"
#include "Lielab/functions.hpp"
#include "Lielab/utils.hpp"

#include <Eigen/Core>

#include <limits>
#include <functional>

namespace Lielab::integrate
{

enum class IVPMethod
{
    Undefined = 0,
    RungeKutta = 1,
    CrouchGrossman = 2,
    MuntheKaas = 3,
};

struct IVPOptions
{
    public:
    
    // IVP Meta-options
    IVPMethod method = IVPMethod::Undefined;

    // Multi-method options
    double dt = std::numeric_limits<double>::quiet_NaN();
    double dt_min = 1e-4;
    double dt_max = 10.0;
    bool variable_time_step = true;

    double reltol = 1e-4;
    double abstol = 1e-8;

    double small = 0.2;
    double large = 10.0;
    double pessimist = 0.9;

    int max_iterations = -1;

    // Runge-Kutta specific options
    RungeKuttaCoefficients coefficients = RungeKuttaCoefficients::RKV87r;

    // Crouch-Grossman specific options
    CrouchGrossmanCoefficients crouch_grossman_coefficients = CrouchGrossmanCoefficients::CG23;

    // Munthe-Kaas specific options
    bool rebase_every_step = true;
};

enum class IVPStatus
{
    ERROR = -99,
    ERROR_MAX_ITERATIONS = -6,
    ERROR_INFS_IN_EVENT = -5,
    ERROR_NANS_IN_EVENT = -4,
    ERROR_INFS_IN_VF = -3,
    ERROR_NANS_IN_VF = -2,
    RUNNING = -1,
    SUCCESS = 0,
    SUCCESS_EVENT = 1,
    SUCCESS_BUT_TOL = 2,
    SUCCESS_EVENT_BUT_TOL = 3,
};

class IVPSolution
{
    public:

    // Metadata
    bool success = false;
    IVPStatus status = IVPStatus::RUNNING;
    std::string message = "";
    double time_to_solution = std::numeric_limits<double>::quiet_NaN();

    // Chunking algorithm params
    ptrdiff_t chunk_size = 2048;
    ptrdiff_t current_index = 0;
    
    // Solution main data
    Eigen::VectorXd t = Eigen::VectorXd::Zero(0);
    std::vector<Lielab::domain::CompositeManifold> y;
    Eigen::MatrixXd ybar = Eigen::MatrixXd::Zero(0, 0);
    std::vector<Lielab::domain::CompositeAlgebra> theta;
    Eigen::MatrixXd thetabar = Eigen::MatrixXd::Zero(0, 0);

    // Debug variables
    Eigen::MatrixXd debug;

    std::string to_string() const;

    IVPSolution();
    IVPSolution(const IVPSolution& other);
    IVPSolution(const size_t num_eoms);
    IVPSolution& operator=(const IVPSolution& other);
    void trim_chunk();
    void trim_chunk(const ptrdiff_t last_index);
    void add_chunk();
    void add_data(const double t_add, const Eigen::VectorXd& ybar_add);
    void add_data(const double t_add, const Eigen::VectorXd& ybar_add, const Eigen::VectorXd& thetabar_add);
};

using EuclideanIVP_event_t = std::function<double(const double, const Eigen::VectorXd&)>;
using EuclideanIVP_vectorfield_t = std::function<Eigen::VectorXd(const double, const Eigen::VectorXd&)>;

struct EuclideanIVPSystem
{
    public:

    EuclideanIVP_event_t event = [](const double t, const Eigen::VectorXd& y){return std::numeric_limits<double>::signaling_NaN();};
    EuclideanIVP_vectorfield_t vectorfield = [](const double t, const Eigen::VectorXd& y){return Eigen::VectorXd::Ones(1)*std::numeric_limits<double>::signaling_NaN();};

    EuclideanIVPSystem(EuclideanIVP_vectorfield_t vf);
};

using HomogeneousIVP_action_t = std::function<Lielab::domain::CompositeManifold(const Lielab::domain::CompositeGroup&, const Lielab::domain::CompositeManifold&)>;
using HomogeneousIVP_connection_t = std::function<Lielab::domain::CompositeAlgebra(const Lielab::domain::CompositeAlgebra&, const Lielab::domain::CompositeAlgebra&)>;
using HomogeneousIVP_coordinates_t = std::function<Lielab::domain::CompositeGroup(const Lielab::domain::CompositeAlgebra&)>;
using HomogeneousIVP_event_t = std::function<double(const double, const Lielab::domain::CompositeManifold&)>;
using HomogeneousIVP_generator_t = std::function<Lielab::domain::CompositeAlgebra(const double, const Lielab::domain::CompositeManifold&)>;

struct HomogeneousIVPSystem
{
    public:

    HomogeneousIVP_action_t action = Lielab::functions::left_Lie_group_action;
    HomogeneousIVP_connection_t connection = [](const Lielab::domain::CompositeAlgebra& a, const Lielab::domain::CompositeAlgebra& b){return Lielab::functions::dexpinv(a, b, 5);};
    HomogeneousIVP_coordinates_t coordinates = Lielab::functions::exp<Lielab::domain::CompositeAlgebra>;
    HomogeneousIVP_event_t event = [](const double t, const Lielab::domain::CompositeManifold& y){return std::numeric_limits<double>::signaling_NaN();};
    HomogeneousIVP_generator_t generator = [](const double t, const Lielab::domain::CompositeManifold& y){return Lielab::domain::CompositeAlgebra({Lielab::domain::rn::from_vector({std::numeric_limits<double>::signaling_NaN()})});};

    HomogeneousIVPSystem(HomogeneousIVP_generator_t vf);
};

using EuclideanClassicHamiltonianIVP_dTdp_t = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;
using EuclideanClassicHamiltonianIVP_dVdq_t = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;
using EuclideanClassicHamiltonianIVP_event_t = std::function<double(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>;

struct EuclideanClassicHamiltonianIVPSystem
{
    public:

    EuclideanClassicHamiltonianIVP_dVdq_t dVdq = [](const Eigen::VectorXd& q){return Lielab::utils::to_VectorXd({std::numeric_limits<double>::signaling_NaN()});};
    EuclideanClassicHamiltonianIVP_dTdp_t dTdp = [](const Eigen::VectorXd& p){return Lielab::utils::to_VectorXd({std::numeric_limits<double>::signaling_NaN()});};
    EuclideanClassicHamiltonianIVP_event_t event = [](const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& p){return std::numeric_limits<double>::signaling_NaN();};

    EuclideanClassicHamiltonianIVPSystem(EuclideanClassicHamiltonianIVP_dVdq_t dVdq_, EuclideanClassicHamiltonianIVP_dTdp_t dTdp_);
};

}

#endif
