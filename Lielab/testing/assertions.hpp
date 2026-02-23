#ifndef LIELAB_TESTING_ASSERTIONS_HPP
#define LIELAB_TESTING_ASSERTIONS_HPP

#include "Lielab/domain.hpp"

namespace Lielab::testing
{

bool check_topology(const Lielab::domain::CompositeManifold& a, const Lielab::domain::CompositeManifold& b);
bool check_topology(const Lielab::domain::CompositeGroup& a, const Lielab::domain::CompositeGroup& b);
bool check_topology(const Lielab::domain::CompositeAlgebra& a, const Lielab::domain::CompositeAlgebra& b);

bool check_almost_equal_tol(const double a, const double b, const double abstol = 1.0e-14, const double reltol = 1.0e-14);
bool check_almost_equal_tol(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const double abstol = 1.0e-14, const double reltol = 1.0e-14);
bool check_almost_equal_tol(const Eigen::MatrixXcd& a, const Eigen::MatrixXcd& b, const double abstol = 1.0e-14, const double reltol = 1.0e-14);
bool check_almost_equal_tol(const Lielab::domain::CompositeManifold& a, const Lielab::domain::CompositeManifold& b, const double abstol = 1.0e-14, const double reltol = 1.0e-14);
bool check_almost_equal_tol(const Lielab::domain::CompositeGroup& a, const Lielab::domain::CompositeGroup& b, const double abstol = 1.0e-14, const double reltol = 1.0e-14);
bool check_almost_equal_tol(const Lielab::domain::CompositeAlgebra& a, const Lielab::domain::CompositeAlgebra& b, const double abstol = 1.0e-14, const double reltol = 1.0e-14);

bool check_almost_equal_nulp(const double a, const double b, const int nulp = 1, const bool gate = false);
bool check_almost_equal_nulp(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const int nulp = 1, const bool gate = false);
bool check_almost_equal_nulp(const Eigen::MatrixXcd& a, const Eigen::MatrixXcd& b, const int nulp = 1, const bool gate = false);
bool check_almost_equal_nulp(const Lielab::domain::CompositeManifold& a, const Lielab::domain::CompositeManifold& b, const int nulp = 1, const bool gate = false);
bool check_almost_equal_nulp(const Lielab::domain::CompositeGroup& a, const Lielab::domain::CompositeGroup& b, const int nulp = 1, const bool gate = false);
bool check_almost_equal_nulp(const Lielab::domain::CompositeAlgebra& a, const Lielab::domain::CompositeAlgebra& b, const int nulp = 1, const bool gate = false);

void assert_handler(const char* expr, const std::string why, const char* file, const char* func, const int line);

}

#ifdef LIELAB_INCLUDE_ASSERTS
#define lielab_assert(expr, why) ((expr) ? (void)0 : Lielab::testing::assert_handler(#expr, why, __FILE__, __func__, __LINE__))
#else
#define lielab_assert(expr, why) ((void)0)
#endif

#endif
