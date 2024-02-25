#ifndef _LIELAB_FUNCTIONS_LOG_HPP
#define _LIELAB_FUNCTIONS_LOG_HPP

namespace lielab
{
    namespace functions
    {
        template<typename LG>
        lielab::domain::lieiii<LG> _log(const LG & G)
        {
            /*!
            * This is the main computationally expensive log
            * function.
            */
            
            return (G.get_ados_representation()).log();
        }

        template<typename LG>
        lielab::domain::lieiii<LG> log(const LG & G, const bool optimize = false)
        {
            /*!
            * This is the logarithm function.
            */
            
            return _log(G);
        }

        template <>
        lielab::domain::rn log(const lielab::domain::RN & G, const bool optimize)
        {
            /*!
            * RN log overload
            */

            lielab::domain::rn out(G.shape);
            out._data = G._data;
            return out;
        }

        template <>
        lielab::domain::so log(const lielab::domain::SO & G, const bool optimize)
        {
            /*!
            * SO logarithm overload
            *
            * Uses Engø's formula for SO(3), otherwise
            * calls the numerical method for SO(n), n != 3.
            * 
            * Source: Engø, Kenth (2001), "On the BCH-formula in so(3)", BIT Numerical Mathematics, 41 (3): 629–632
            * 
            * Engø's formula is a local solution. Do not use when theta >= pi/2
            *
            * TODO: This directly uses _data. Replace with get_ados_representation().
            */

            if (optimize && G.shape == 3)
            {
                lielab::domain::so out(3);
                Eigen::MatrixXd A = 1.0/2.0*(G._data - G._data.transpose());
                double normA = std::sqrt(-1.0/2.0 * (A*A).trace());
                out._data = std::asin(normA)/normA*A;
                return out;
            }

            return _log(G);
        }
    }
}

#endif
