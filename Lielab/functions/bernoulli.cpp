#include "bernoulli.hpp"

#include <vector>

namespace Lielab::functions
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

}

