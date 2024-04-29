#ifndef LIELAB_ABSTRACT_ABSTRACT_ALL_HPP_
#define LIELAB_ABSTRACT_ABSTRACT_ALL_HPP_

#include <complex>
#include <concepts>

namespace Lielab
{
namespace abstract
{
template <typename T>
concept Real = std::integral<T> || std::floating_point<T>;

template <typename T>
struct is_complex : std::false_type {};

template <std::floating_point T>
struct is_complex<std::complex<T>> : std::true_type {};

template <typename T>
concept Imaginary = Real<T> || is_complex<T>::value;

template <typename T>
concept Additive = requires (T x)
{
    {x + x};
    {x - x};
};

template <typename T>
concept Multiplicative = requires (T x)
{
    {x * x};
};

template <typename T>
concept ScalarMultiplicative = requires (T x)
{
    {1 * x};
    {x * 1};
    {x / 1};
};

template <typename T>
concept Invertible = requires (T x)
{
    {x.inverse()};
};

template <typename T>
concept Algebra = Additive<T> && Multiplicative<T> && ScalarMultiplicative<T> && requires (T x)
{
    {x.abelian};
    {x.basis};
    {x.get_dimension()};
    {x.get_vector()};
};

template <typename T>
concept LieAlgebra = Algebra<T> && requires (T x)
{
    {x.shape};
    {x.get_ados_representation()};
};

template <typename T>
concept Group = Invertible<T> && Multiplicative<T> && requires (T x)
{
    {x.abelian};
};

template <typename T>
concept LieGroup = Group<T> && requires (T x)
{
    {x.shape};
    {x.get_ados_representation()};
};
}
}

#endif
