#ifndef _LIELAB_FUNCTIONS_COAD_HPP
#define _LIELAB_FUNCTIONS_COAD_HPP

// namespace Lielab
// {
//     namespace functions
//     {
//         template<Lielab::abstract::LieGroup LG>
//         Lielab::domain::cola _coAd(const LG & G, const Lielab::domain::cola & mu)
//         {
//             // TODO: This can be sped up with .get_vector()
//             // TODO: coAd needs a rework overall. No more cola's.
//             std::vector<domain::lieiii<LG>> basis;
//             const size_t d1 = G.get_dimension();
//             const size_t d2 = mu.get_dimension();
//             double norm = 0.0;

//             if (d1 != d2)
//             {
//                 throw SizeError("Input data malformed: dimension " + std::to_string(d1) + " != " + std::to_string(d2) + ".");
//             }
            
//             Eigen::VectorXd out(d1);

//             const size_t s = G.shape;

//             for (int ii = 0; ii < d1; ii++)
//             {
//                 basis.push_back(domain::lieiii<LG>::basis(ii, s));
//             }

//             for (int ii = 0; ii < d1; ii++)
//             {
//                 out(ii) = copair(mu, Ad(G.inverse(), basis[ii]));
//             }

//             return out;
//         }

//         template<Lielab::abstract::LieGroup LG>
//         Lielab::domain::cola coAd(const LG & G, const Lielab::domain::cola & mu)
//         {
//             /*! \f{equation*}{ (G, \mathfrak{g}^*) \rightarrow \mathfrak{g}^* \f}
//             *
//             * CoAdjoint function.
//             */

//             return _coAd(G, mu);
//         }

//         template<>
//         Lielab::domain::cola coAd(const Lielab::domain::GL & G, const Lielab::domain::cola & mu)
//         {
//             /*! \f{equation*}{ (GL, \mathfrak{gl}^*) \rightarrow \mathfrak{gl}^*}
//             *
//             * coAd overload for known shortcut formulas
//             */
            
//             if (G.shape == mu.get_dimension())
//             {
//                 return G._data*mu._data;
//             }
            
//             // TODO: Find out why the following line doesn't compile and fix
//             // return _coAd(G, mu);
//         }

//         template<>
//         Lielab::domain::cola coAd(const Lielab::domain::SO & G, const Lielab::domain::cola & mu)
//         {
//             /*! \f{equation*}{ (SO, \mathfrak{so}^*) \rightarrow \mathfrak{so}^*}
//             *
//             * coAd overload for known shortcut formulas
//             */
            
//             if (G.shape == mu.get_dimension())
//             {
//                 return G._data*mu._data;
//             }
            
//             return _coAd(G, mu);
//         }
//     }
// }

#endif
