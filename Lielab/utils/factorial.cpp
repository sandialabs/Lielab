#include "factorial.hpp"

namespace Lielab::utils
{

unsigned int factorial(const unsigned int n)
{
    /*!
    * An inefficient factorial function.
    */
    if (n == 0)
    {
        return 1;
    }

    return (factorial(n - 1)*n);
}

// template <std::size_t Length>
// constexpr auto factorial_precomputed = precompute_<Length>(factorial);

}
