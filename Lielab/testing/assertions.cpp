#include "assertions.hpp"

#include "Lielab/domain.hpp"

#include "fmt/core.h"

#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Lielab::testing
{

bool check_topology(const Lielab::domain::CompositeManifold& a, const Lielab::domain::CompositeManifold& b)
{
    if (a.point.size() != b.point.size())
    {
        return false;
    }

    const auto adims = a.get_dimensions();
    const auto bdims = b.get_dimensions();

    for (size_t ii = 0; ii < a.point.size(); ii++)
    {
        const size_t inda = a.point[ii].index();
        const size_t indb = b.point[ii].index();
        if (inda != indb)
        {
            return false;
        }

        if (adims[ii] != bdims[ii])
        {
            return false;
        }
    }

    return true;
}

bool check_topology(const Lielab::domain::CompositeGroup& a, const Lielab::domain::CompositeGroup& b)
{
    if (a.point.size() != b.point.size())
    {
        return false;
    }

    const auto adims = a.get_dimensions();
    const auto bdims = b.get_dimensions();

    for (size_t ii = 0; ii < a.point.size(); ii++)
    {
        const size_t inda = a.point[ii].index();
        const size_t indb = b.point[ii].index();
        if (inda != indb)
        {
            return false;
        }

        if (adims[ii] != bdims[ii])
        {
            return false;
        }
    }

    return true;
}

bool check_topology(const Lielab::domain::CompositeAlgebra& a, const Lielab::domain::CompositeAlgebra& b)
{
    if (a.point.size() != b.point.size())
    {
        return false;
    }

    const auto adims = a.get_dimensions();
    const auto bdims = b.get_dimensions();

    for (size_t ii = 0; ii < a.point.size(); ii++)
    {
        const size_t inda = a.point[ii].index();
        const size_t indb = b.point[ii].index();
        if (inda != indb)
        {
            return false;
        }

        if (adims[ii] != bdims[ii])
        {
            return false;
        }
    }

    return true;
}

bool check_almost_equal_tol(const double a, const double b, const double abstol, const double reltol)
{
    constexpr double inf = std::numeric_limits<double>::infinity();

    if (a == b) return true;

    // Any nans always return false
    if (std::isnan(a) || std::isnan(b)) return false;

    // Infs are true only if the signage is the same
    if (a == inf && b == inf) return true;
    if (a == inf && b == -inf) return false;
    if (a == -inf && b == inf) return false;
    if (a == -inf && b == -inf) return true;
    if (std::isinf(a) || std::isinf(b)) return false;

    return !(std::abs(a - b) > (abstol + reltol*std::abs(b)));
}

bool check_almost_equal_tol(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const double abstol, const double reltol)
{
    if (a.rows() != b.rows()) return false;
    if (a.cols() != b.cols()) return false;

    for (int ii = 0; ii < a.rows(); ii++)
    {
        for (int jj = 0; jj < a.cols(); jj++)
        {
            const bool res = check_almost_equal_tol(a(ii,jj), b(ii,jj), abstol, reltol);
            if (!res) return false;
        }
    }

    return true;
}

bool check_almost_equal_tol(const Eigen::MatrixXcd& a, const Eigen::MatrixXcd& b, const double abstol, const double reltol)
{
    if (a.rows() != b.rows()) return false;
    if (a.cols() != b.cols()) return false;

    for (int ii = 0; ii < a.rows(); ii++)
    {
        for (int jj = 0; jj < a.cols(); jj++)
        {
            const bool res = check_almost_equal_tol(a(ii,jj).real(), b(ii,jj).real(), abstol, reltol);
            if (!res) return false;
            const bool resj = check_almost_equal_tol(a(ii,jj).imag(), b(ii,jj).imag(), abstol, reltol);
            if (!resj) return false;
        }
    }

    return true;
}

bool check_almost_equal_tol(const Lielab::domain::CompositeManifold& a, const Lielab::domain::CompositeManifold& b, const double abstol, const double reltol)
{
    if (!check_topology(a, b)) return false;

    const Eigen::VectorXd abar = a.serialize();
    const Eigen::VectorXd bbar = b.serialize();

    if (abar.size() != bbar.size())
    {
        // This should never be called. It would have already been captured
        // by the check_topology() above.
        return false;
    }

    for (ptrdiff_t ii = 0; ii < abar.size(); ii++)
    {
        const double ai = abar(ii);
        const double bi = bbar(ii);

        if (!check_almost_equal_tol(ai, bi, abstol, reltol)) return false;
    }

    return true;
}

bool check_almost_equal_tol(const Lielab::domain::CompositeGroup& a, const Lielab::domain::CompositeGroup& b, const double abstol, const double reltol)
{
    if (!check_topology(a, b)) return false;

    const Eigen::VectorXd abar = a.serialize();
    const Eigen::VectorXd bbar = b.serialize();

    if (abar.size() != bbar.size())
    {
        // This should never be called. It would have already been captured
        // by the check_topology() above.
        return false;
    }

    for (ptrdiff_t ii = 0; ii < abar.size(); ii++)
    {
        const double ai = abar(ii);
        const double bi = bbar(ii);

        if (!check_almost_equal_tol(ai, bi, abstol, reltol)) return false;
    }

    return true;
}

bool check_almost_equal_tol(const Lielab::domain::CompositeAlgebra& a, const Lielab::domain::CompositeAlgebra& b, const double abstol, const double reltol)
{
    if (!check_topology(a, b)) return false;

    const Eigen::VectorXd abar = a.get_vector();
    const Eigen::VectorXd bbar = b.get_vector();

    if (abar.size() != bbar.size())
    {
        // This should never be called. It would have already been captured
        // by the check_topology() above.
        return false;
    }

    for (ptrdiff_t ii = 0; ii < abar.size(); ii++)
    {
        const double ai = abar(ii);
        const double bi = bbar(ii);

        if (!check_almost_equal_tol(ai, bi, abstol, reltol)) return false;
    }

    return true;
}

bool check_almost_equal_nulp(const double a, const double b, const int nulp, const bool gate)
{
    constexpr double inf = std::numeric_limits<double>::infinity();

    if (a == b) return true;

    // Any nans always return false
    if (std::isnan(a) || std::isnan(b)) return false;

    // Infs are true only if the signage is the same
    if (a == inf && b == inf) return true;
    if (a == inf && b == -inf) return false;
    if (a == -inf && b == inf) return false;
    if (a == -inf && b == -inf) return true;
    if (std::isinf(a) || std::isinf(b)) return false;

    const double aabs = std::abs(a);
    const double babs = std::abs(b);

    double max_ab = std::max(aabs, babs);
    if (gate) max_ab = std::max(max_ab, 1.0);

    const double nextu = std::nextafter(max_ab, std::numeric_limits<double>::infinity());
    const double nextd = std::nextafter(max_ab, -std::numeric_limits<double>::infinity());
    const double ulp = std::min(std::abs(nextu - max_ab), std::abs(max_ab - nextd));

    return !(std::abs(a - b) > nulp*ulp);
}

bool check_almost_equal_nulp(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const int nulp, const bool gate)
{
    if (a.rows() != b.rows()) return false;
    if (a.cols() != b.cols()) return false;

    for (int ii = 0; ii < a.rows(); ii++)
    {
        for (int jj = 0; jj < a.cols(); jj++)
        {
            const bool res = check_almost_equal_nulp(a(ii,jj), b(ii,jj), nulp, gate);
            if (!res) return false;
        }
    }

    return true;
}

bool check_almost_equal_nulp(const Eigen::MatrixXcd& a, const Eigen::MatrixXcd& b, const int nulp, const bool gate)
{
    if (a.rows() != b.rows()) return false;
    if (a.cols() != b.cols()) return false;

    for (int ii = 0; ii < a.rows(); ii++)
    {
        for (int jj = 0; jj < a.cols(); jj++)
        {
            const bool res = check_almost_equal_nulp(a(ii,jj).real(), b(ii,jj).real(), nulp, gate);
            if (!res) return false;
            const bool resj = check_almost_equal_nulp(a(ii,jj).imag(), b(ii,jj).imag(), nulp, gate);
            if (!resj) return false;
        }
    }

    return true;
}

bool check_almost_equal_nulp(const Lielab::domain::CompositeManifold& a, const Lielab::domain::CompositeManifold& b, const int nulp, const bool gate)
{
    if (!check_topology(a, b)) return false;

    const Eigen::VectorXd abar = a.serialize();
    const Eigen::VectorXd bbar = b.serialize();

    if (abar.size() != bbar.size())
    {
        // This should never be called. It would have already been captured
        // by the check_topology() above.
        return false;
    }

    for (ptrdiff_t ii = 0; ii < abar.size(); ii++)
    {
        const double ai = abar(ii);
        const double bi = bbar(ii);

        if (!check_almost_equal_nulp(ai, bi, nulp, gate)) return false;
    }

    return true;
}

bool check_almost_equal_nulp(const Lielab::domain::CompositeGroup& a, const Lielab::domain::CompositeGroup& b, const int nulp, const bool gate)
{
    if (!check_topology(a, b)) return false;

    const Eigen::VectorXd abar = a.serialize();
    const Eigen::VectorXd bbar = b.serialize();

    if (abar.size() != bbar.size())
    {
        // This should never be called. It would have already been captured
        // by the check_topology() above.
        return false;
    }

    for (ptrdiff_t ii = 0; ii < abar.size(); ii++)
    {
        const double ai = abar(ii);
        const double bi = bbar(ii);

        if (!check_almost_equal_nulp(ai, bi, nulp, gate)) return false;
    }

    return true;
}

bool check_almost_equal_nulp(const Lielab::domain::CompositeAlgebra& a, const Lielab::domain::CompositeAlgebra& b, const int nulp, const bool gate)
{
    if (!check_topology(a, b)) return false;

    const Eigen::VectorXd abar = a.get_vector();
    const Eigen::VectorXd bbar = b.get_vector();

    if (abar.size() != bbar.size())
    {
        // This should never be called. It would have already been captured
        // by the check_topology() above.
        return false;
    }

    for (ptrdiff_t ii = 0; ii < abar.size(); ii++)
    {
        const double ai = abar(ii);
        const double bi = bbar(ii);

        if (!check_almost_equal_nulp(ai, bi, nulp, gate)) return false;
    }

    return true;
}

void assert_handler(const char* expr, const std::string why, const char* file, const char* func, const int line)
{
    const std::string exprstr = expr;
    std::filesystem::path filepath(file);
    const std::string filename = filepath.filename().string();

    std::string msg = "\x1b[31mLielab runtime error\n";
    msg += fmt::format("Pass condition: {}\n", expr);
    msg += fmt::format("Reason: {}\n", why);
    msg += fmt::format("Where: {}() located in {} (Line {})", func, filename, line);
    msg += "\x1b[0m";
    throw std::runtime_error(msg);
}

}
