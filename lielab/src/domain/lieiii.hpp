#ifndef _LIELAB_DOMAIN_LIEIII_HPP
#define _LIELAB_DOMAIN_LIEIII_HPP

#include <type_traits>

namespace lielab
{
    namespace domain
    {
        /*! \f{equation}{ \mathfrak{g} \leftrightarrow G \f}
        * Implements Lie's third theorem as a template macro.
        */
        template<class T>
        using lieiii = typename std::conditional<std::is_same<T, gl>::value, GL,
                       typename std::conditional<std::is_same<T, GL>::value, gl,
                       typename std::conditional<std::is_same<T, rn>::value, RN,
                       typename std::conditional<std::is_same<T, RN>::value, rn,
                       typename std::conditional<std::is_same<T, so>::value, SO,
                       typename std::conditional<std::is_same<T, SO>::value, so,
                       typename std::conditional<std::is_same<T, sp>::value, SP,
                       typename std::conditional<std::is_same<T, SP>::value, sp,
                       typename std::conditional<std::is_same<T, su>::value, SU,
                       typename std::conditional<std::is_same<T, SU>::value, su,
                       Eigen::MatrixXcd>::type>::type>::type>::type>::type>::type>::type>::type>::type>::type; // Default false is complex matrix, throw compile time error instead??
    }
}

#endif
