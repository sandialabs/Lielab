#ifndef LIELAB_CONSTANTS_CONSTANTS_ALL_HPP_
#define LIELAB_CONSTANTS_CONSTANTS_ALL_HPP_

#include <array>

namespace Lielab
{
namespace constants
{
/*!
 * pi out to 64 decimal places.
 */

template <typename T>
constexpr T PI = 3.1415926535897932384626433832795028841971693993751058209749445923;

/*!
 * sqrt of 3 out to 64 decimal places.
 */

template <typename T>
constexpr T SQRT3 = 1.7320508075688772935274463415058723669428052538103806280558069795;
}
}

#endif
