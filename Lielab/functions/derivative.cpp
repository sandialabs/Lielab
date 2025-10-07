#include "derivative.hpp"

#include "log.hpp"

#include "Lielab/domain.hpp"

namespace Lielab::functions
{

Lielab::domain::CompositeAlgebra forward_difference(std::function<Lielab::domain::CompositeGroup(const double)> & fun, const double t, const double dt)
{
    /*!

    Foward difference based on Solà's concept of a "right" Jacobian. [1]

    Arguments
    ---------
    @param[in] fun Function mapping R -> G.
    @param[in] t Location to take derivative at.
    @param[in] dt Finite difference distance.
    @param[out] out Derivative of input function at T_t G = g.

    References
    ----------
    [1] - Solà, Joan, Jeremie Deray, and Dinesh Atchuthan. "A micro lie theory for state estimation in robotics." arXiv preprint arXiv:1812.01537 (2018).

    */
    
    using namespace Lielab::domain;
    using namespace Lielab::functions;

    const CompositeGroup f0 = fun(t);
    const CompositeGroup f1 = fun(t+dt);
    const CompositeAlgebra out = log(f0.inverse()*f1)/dt;
    return out;
}

Lielab::domain::CompositeAlgebra backward_difference(std::function<Lielab::domain::CompositeGroup(const double)> & fun, const double t, const double dt)
{
    /*!

    Backward difference based on Solà's concept of a "right" Jacobian. [1]

    Arguments
    ---------
    @param[in] fun Function mapping R -> G.
    @param[in] t Location to take derivative at.
    @param[in] dt Finite difference distance.
    @param[out] out Derivative of input function at T_t G = g.

    References
    ----------
    [1] - Solà, Joan, Jeremie Deray, and Dinesh Atchuthan. "A micro lie theory for state estimation in robotics." arXiv preprint arXiv:1812.01537 (2018).

    */
    
    using namespace Lielab::domain;
    using namespace Lielab::functions;

    const CompositeGroup f0 = fun(t-dt);
    const CompositeGroup f1 = fun(t);
    const CompositeAlgebra out = log(f0.inverse()*f1)/dt;
    return out;
}

Lielab::domain::CompositeAlgebra central_difference(std::function<Lielab::domain::CompositeGroup(const double)> & fun, const double t, const double dt)
{
    /*!

    Central difference based on Solà's concept of a "right" Jacobian. [1]

    Arguments
    ---------
    @param[in] fun Function mapping R -> G.
    @param[in] t Location to take derivative at.
    @param[in] dt Finite difference distance.
    @param[out] out Derivative of input function at T_t G = g.

    References
    ----------
    [1] - Solà, Joan, Jeremie Deray, and Dinesh Atchuthan. "A micro lie theory for state estimation in robotics." arXiv preprint arXiv:1812.01537 (2018).

    */
    
    using namespace Lielab::domain;
    using namespace Lielab::functions;

    const CompositeGroup f0 = fun(t-dt);
    const CompositeGroup f1 = fun(t+dt);
    const CompositeAlgebra out = log(f0.inverse()*f1)/(2.0*dt);
    return out;
}

}
