#ifndef _LIELAB_DOMAIN_LIEIII_HPP
#define _LIELAB_DOMAIN_LIEIII_HPP

#include <type_traits>

namespace Lielab
{
namespace domain
{
/*! \f{equation}{ \mathfrak{g} \leftrightarrow G \f}
* Implements Lie's third theorem as a template macro.
*/
using std::conditional;
using std::is_same;
template<class T>
using lieiii = typename conditional<is_same<T, gl>::value, GL,
               typename conditional<is_same<T, GL>::value, gl,
               typename conditional<is_same<T, rn>::value, RN,
               typename conditional<is_same<T, RN>::value, rn,
               typename conditional<is_same<T, se>::value, SE,
               typename conditional<is_same<T, SE>::value, se,
               typename conditional<is_same<T, so>::value, SO,
               typename conditional<is_same<T, SO>::value, so,
               typename conditional<is_same<T, sp>::value, SP,
               typename conditional<std::is_same<T, SP>::value, sp,
               typename conditional<std::is_same<T, su>::value, SU,
               typename conditional<std::is_same<T, SU>::value, su,
               Eigen::MatrixXcd>::type>::type>::type>::type>::type>::type>::type>::type>::type>::type>::type>::type; // Default false is complex matrix, throw compile time error instead??
}
}

#endif
