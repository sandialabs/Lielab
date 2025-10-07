#ifndef LIELAB_DOMAIN_LIEIII_HPP
#define LIELAB_DOMAIN_LIEIII_HPP

#include "liealgebras.hpp"
#include "liegroups.hpp"

#include <Eigen/Core>

#include <type_traits>

namespace Lielab::domain
{

template <typename T>
struct _LieIII;

template <>
struct _LieIII<cn>{using type = CN;};
template <>
struct _LieIII<CN>{using type = cn;};
template <>
struct _LieIII<glr>{using type = GLR;};
template <>
struct _LieIII<GLR>{using type = glr;};
template <>
struct _LieIII<glc>{using type = GLC;};
template <>
struct _LieIII<GLC>{using type = glc;};
template <>
struct _LieIII<rn>{using type = RN;};
template <>
struct _LieIII<RN>{using type = rn;};
template <>
struct _LieIII<se>{using type = SE;};
template <>
struct _LieIII<SE>{using type = se;};
template <>
struct _LieIII<so>{using type = SO;};
template <>
struct _LieIII<SO>{using type = so;};
template <>
struct _LieIII<sp>{using type = SP;};
template <>
struct _LieIII<SP>{using type = sp;};
template <>
struct _LieIII<su>{using type = SU;};
template <>
struct _LieIII<SU>{using type = su;};

template <>
struct _LieIII<CompositeAlgebra>{using type = CompositeGroup;};
template <>
struct _LieIII<CompositeGroup>{using type = CompositeAlgebra;};

template <typename T>
struct _LieIII{using type = Eigen::MatrixXcd;};

/*! \f{equation}{ \mathfrak{g} \leftrightarrow G \f}
* Implements Lie's third theorem as a template macro.
*/
template <typename T>
using LieIII = typename _LieIII<T>::type;

}

#endif
