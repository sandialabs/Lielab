#include <Lielab.hpp>
#include <string>
#include <Eigen/Core>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "pybind_domain.hpp"
#include "domain/liealgebras/pybind_cn.hpp"
#include "domain/liealgebras/pybind_glc.hpp"
#include "domain/liealgebras/pybind_glr.hpp"
#include "domain/liealgebras/pybind_rn.hpp"
#include "domain/liealgebras/pybind_se.hpp"
#include "domain/liealgebras/pybind_so.hpp"
#include "domain/liealgebras/pybind_sp.hpp"
#include "domain/liealgebras/pybind_su.hpp"
#include "domain/liealgebras/pybind_CompositeAlgebra.hpp"

#include "domain/liegroups/pybind_CN.hpp"
#include "domain/liegroups/pybind_GLC.hpp"
#include "domain/liegroups/pybind_GLR.hpp"
#include "domain/liegroups/pybind_RN.hpp"
#include "domain/liegroups/pybind_SE.hpp"
#include "domain/liegroups/pybind_SO.hpp"
#include "domain/liegroups/pybind_SP.hpp"
#include "domain/liegroups/pybind_SU.hpp"
#include "domain/liegroups/pybind_CompositeGroup.hpp"

#include "domain/smoothmanifolds/pybind_Grassmannian.hpp"
#include "domain/smoothmanifolds/pybind_CompositeManifold.hpp"

namespace py = pybind11;

void bind_domain(py::module& m_domain)
{
    bind_cn(m_domain);
    bind_glc(m_domain);
    bind_glr(m_domain);
    bind_rn(m_domain);
    bind_se(m_domain);
    bind_so(m_domain);
    bind_sp(m_domain);
    bind_su(m_domain);
    bind_CompositeAlgebra(m_domain);

    bind_CN(m_domain);
    bind_GLC(m_domain);
    bind_GLR(m_domain);
    bind_RN(m_domain);
    bind_SE(m_domain);
    bind_SO(m_domain);
    bind_SP(m_domain);
    bind_SU(m_domain);
    bind_CompositeGroup(m_domain);

    bind_Grassmannian(m_domain);
    bind_CompositeManifold(m_domain);
}
