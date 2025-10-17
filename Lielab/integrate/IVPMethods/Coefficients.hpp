#ifndef LIELAB_INTEGRATE_COEFFICIENTS_HPP
#define LIELAB_INTEGRATE_COEFFICIENTS_HPP

#include <Eigen/Core>

#include <array>
#include <tuple>

namespace Lielab::integrate
{

enum class Coefficients {FE1, RK3, RK4a, RK4b, RK5a, RK5b, // Runge-Kutta classics
    RKF12a, RKF12b, RKF23a, RKF23b, RKF34a, RKF34b, RKF45a, RKF45b, RKF56, RKF67, RKF78, RKF8, // Fehlberg's methods
    RKDP54_7M, // Dormand-Prince
    RKV65e, RKV65r, RKV76e, RKV76r, RKV87e, RKV87r, RKV98e, RKV98r, // Verner's RK methods
    CG4a, CG5a, // Crouch-Grossman methods
    BE1}; // Implicit RK classics

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, int, int, bool, bool> get_butcher_tableau(const Coefficients method);


}

#endif
