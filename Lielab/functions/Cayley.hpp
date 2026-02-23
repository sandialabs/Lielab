#ifndef LIELAB_FUNCTIONS_CAYLEY_HPP
#define LIELAB_FUNCTIONS_CAYLEY_HPP

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

// cay
template <typename LA> Lielab::domain::LieIII<LA> cay(const LA& x);

// cayinv
template <typename LG> Lielab::domain::LieIII<LG> cayinv(const LG& g);

// dcay
// template <typename LA> Lielab::domain::glr dcay(const LA& x);
template <typename LA> LA dcay(const LA& x, const LA& y);

// dcayinv
// template <typename LA> Lielab::domain::glr dcayinv(const LA& x);
template <typename LA> LA dcayinv(const LA& x, const LA& y);

// cay2
template <typename LA> Lielab::domain::LieIII<LA> cay2(const LA& x);

// Composite overloads
template <> Lielab::domain::CompositeGroup cay(const Lielab::domain::CompositeAlgebra& x);
template <> Lielab::domain::CompositeAlgebra cayinv(const Lielab::domain::CompositeGroup& g);
template <> Lielab::domain::CompositeAlgebra dcay(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y);
template <> Lielab::domain::CompositeAlgebra dcayinv(const Lielab::domain::CompositeAlgebra& x, const Lielab::domain::CompositeAlgebra& y);

}

#include "Cayley.tpp"

#endif
