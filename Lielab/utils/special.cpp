#include "special.hpp"

#include <cmath>
#include <limits>
#include <vector>

namespace Lielab::utils
{

double bernoulli(const int num, const int sign)
{
    /*!
    * Returns the nth Bernoulli number.
    * 
    * @param[in] num The Bernoulli number to calculate
    * @param[in] sign The sign convention to use for num=1. Defaults to -1.
    */

    if (num == 1)
    {
        return sign*1.0/2.0;
    }

    std::vector<double> A(num+1);
    for (int ii = 0; ii <= num; ii++)
    {
        A[ii] = 1.0/(ii + 1.0);
        for (int jj = ii; jj > 0; jj--)
        {
            A[jj-1] = jj*(A[jj-1] - A[jj]);
        }
    }
    return A[0];
}

int sign(const double x)
{
    if (x == 0.0) return 0;
    const double inf = std::numeric_limits<double>::infinity();
    if (x == -inf) return -1;
    if (x == inf) return 1;
    if (x < 0.0) return -1;
    if (x > 0.0) return 1;

    // This should never get called.
    return 0;
}

double sinc(const double x)
{
    const double absx = std::abs(x);
    if (absx < 1e-10) return 1.0;
    return std::sin(x)/x;
}

}

