#ifndef _LIELAB_CONSTANTS_HPP
#define _LIELAB_CONSTANTS_HPP

#include <array>

namespace lielab
{
    namespace constants
    {
        /*!
        * pi out to 64 decimal places.
        */

        template <typename T>
        constexpr T PI = 3.1415926535897932384626433832795028841971693993751058209749445923;

        /*!
        * Pauli matrices.
        *
        * Indexing:
        *     (0) -> (0, 0)
        *     (1) -> (0, 1)
        *     (2) -> (1, 0)
        *     (3) -> (1, 1)
        */
        
        template <typename T>
        constexpr std::array<std::complex<T>, 4> pauli1 = {std::complex<T>(0.0, 0.0),  std::complex<T>(1.0, 0.0),
                                                           std::complex<T>(1.0, 0.0),  std::complex<T>(0.0, 0.0)};
        
        template <typename T>
        constexpr std::array<std::complex<T>, 4> pauli2 = {std::complex<T>(0.0, 0.0), -std::complex<T>(0.0, 1.0),
                                                           std::complex<T>(0.0, 1.0),  std::complex<T>(0.0, 0.0)};
        
        template <typename T>
        constexpr std::array<std::complex<T>, 4> pauli3 = {std::complex<T>(1.0, 0.0),  std::complex<T>(0.0, 0.0),
                                                           std::complex<T>(0.0, 0.0), -std::complex<T>(1.0, 0.0)};
        
        /*!
        * Gell-Mann matrices.
        *
        * Indexing:
        *     (0) -> (0, 0)
        *     (1) -> (0, 1)
        *     (2) -> (0, 2)
        *     (3) -> (1, 0)
        *     (4) -> (1, 1)
        *     (5) -> (1, 2)
        *     (6) -> (2, 0)
        *     (7) -> (2, 1)
        *     (8) -> (2, 2)
        */

        template <typename T>
        constexpr std::array<std::complex<T>, 9> gellmann1 = {std::complex<T>( 0.0, 0.0), std::complex<T>( 1.0, 0.0), std::complex<T>( 0.0, 0.0),
                                                              std::complex<T>( 1.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0),
                                                              std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0)};
        
        template <typename T>
        constexpr std::array<std::complex<T>, 9> gellmann2 = {std::complex<T>( 0.0, 0.0),-std::complex<T>( 0.0, 1.0), std::complex<T>( 0.0, 0.0),
                                                              std::complex<T>( 0.0, 1.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0),
                                                              std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0)};

        template <typename T>
        constexpr std::array<std::complex<T>, 9> gellmann3 = {std::complex<T>( 1.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0),
                                                              std::complex<T>( 0.0, 0.0),-std::complex<T>( 1.0, 0.0), std::complex<T>( 0.0, 0.0),
                                                              std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0)};
        
        template <typename T>
        constexpr std::array<std::complex<T>, 9> gellmann4 = {std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 1.0, 0.0),
                                                              std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0),
                                                              std::complex<T>( 1.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0)};
        
        template <typename T>
        constexpr std::array<std::complex<T>, 9> gellmann5 = {std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0,-1.0),
                                                              std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0),
                                                              std::complex<T>( 0.0, 1.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0)};
        
        template <typename T>
        constexpr std::array<std::complex<T>, 9> gellmann6 = {std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0),
                                                              std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 1.0, 0.0),
                                                              std::complex<T>( 0.0, 0.0), std::complex<T>( 1.0, 0.0), std::complex<T>( 0.0, 0.0)};

        template <typename T>
        constexpr std::array<std::complex<T>, 9> gellmann7 = {std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0),
                                                              std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0,-1.0),
                                                              std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 1.0), std::complex<T>( 0.0, 0.0)};
        
        /*!
        * sqrt of 3 out to 64 decimal places.
        */
        
        template <typename T>
        constexpr T SQRT3 = 1.7320508075688772935274463415058723669428052538103806280558069795;

        template <typename T>
        constexpr std::array<std::complex<T>, 9> gellmann8 = {std::complex<T>( 1.0, 0.0) / SQRT3<T>, std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0),
                                                              std::complex<T>( 0.0, 0.0), std::complex<T>( 1.0, 0.0) / SQRT3<T>, std::complex<T>( 0.0, 0.0),
                                                              std::complex<T>( 0.0, 0.0), std::complex<T>( 0.0, 0.0), std::complex<T>(-2.0, 0.0) / SQRT3<T>};
    }
}

#endif
