#ifndef LIELAB_INTEGRATE_COEFFICIENTS_HPP
#define LIELAB_INTEGRATE_COEFFICIENTS_HPP

#include <Eigen/Core>

#include <array>
#include <tuple>

namespace Lielab::integrate
{

constexpr int COEFFICIENTS_STORAGE = 20;

enum class RungeKuttaCoefficients {FE1, RK3, RK4a, RK4b, RK5a, RK5b, // Runge-Kutta classics
    RKF12a, RKF12b, RKF23a, RKF23b, RKF34a, RKF34b, RKF45a, RKF45b, RKF56, RKF67, RKF78, RKF8, // Fehlberg's methods
    RKDP54_7M, // Dormand-Prince
    RKV65e, RKV65r, RKV76e, RKV76r, RKV87e, RKV87r, RKV98e, RKV98r, // Verner's RK methods
    BE1, // Implicit RK classics
    LG2, LG4, LG4s, LG6, LG6s, // Legendre-Gauss methods
    Lobatto3A2, Lobatto3A4, Lobatto3A6}; // Gauss-Lobatto methods

std::tuple<std::array<std::array<double, COEFFICIENTS_STORAGE>, COEFFICIENTS_STORAGE>, std::array<double, COEFFICIENTS_STORAGE>, std::array<double, COEFFICIENTS_STORAGE>, std::array<double, COEFFICIENTS_STORAGE>, std::array<double, COEFFICIENTS_STORAGE>, int, int, bool, bool> get_butcher_tableau(const RungeKuttaCoefficients method);

enum class CrouchGrossmanCoefficients {CG23,
    CG4a, CG5a};

std::tuple<std::array<std::array<double, COEFFICIENTS_STORAGE>, COEFFICIENTS_STORAGE>, std::array<double, COEFFICIENTS_STORAGE>, std::array<double, COEFFICIENTS_STORAGE>, std::array<double, COEFFICIENTS_STORAGE>, std::array<double, COEFFICIENTS_STORAGE>, int, int, bool, bool> get_crouch_grossman_coefficients(const CrouchGrossmanCoefficients method);

}

#endif
