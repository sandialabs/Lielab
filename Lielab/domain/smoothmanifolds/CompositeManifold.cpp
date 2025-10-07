#include "CompositeManifold.hpp"

#include "../liealgebras.hpp"
#include "../liegroups.hpp"

#include "Lielab/utils.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <exception>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <variant>

namespace Lielab::domain
{

std::string CompositeManifold::to_string() const
{
    std::string out = "";
    const size_t sz = this->space.size();
    std::vector<size_t> shapes = this->get_shapes();

    for (size_t ii = 0; ii < sz; ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeManifold::INDEX_CN)
        {
            out += std::get<CN>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_GLR)
        {
            out += std::get<GLR>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_GLC)
        {
            out += std::get<GLC>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_RN)
        {
            out += std::get<RN>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_SE)
        {
            out += std::get<SE>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_SO)
        {
            out += std::get<SO>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_SP)
        {
            out += std::get<SP>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_SU)
        {
            out += std::get<SU>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_cn)
        {
            out += std::get<cn>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_glr)
        {
            out += std::get<glr>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_glc)
        {
            out += std::get<glc>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_rn)
        {
            out += std::get<rn>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_se)
        {
            out += std::get<se>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_so)
        {
            out += std::get<so>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_sp)
        {
            out += std::get<sp>(this->space[ii]).to_string();
        }
        else if (ind == CompositeManifold::INDEX_su)
        {
            out += std::get<su>(this->space[ii]).to_string();
        }

        if (ii < sz-1)
        {
            out += " ";
            out += "x";
            out += " ";
        }
    }

    return out;
}

CompositeManifold::CompositeManifold()
{
    
}

CompositeManifold::CompositeManifold(const size_t n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{CompositeManifold} \f}
    *
    * Constructor instantiating a \f$\mathfrak{CompositeManifold}\f$ as a single glc object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::CompositeManifold x(3), y(4), z(5);
    * 
    * @param[in] n The shape of the data matrix.
    */

    this->space.push_back(glc(n));
}

CompositeManifold CompositeManifold::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{CompositeManifold} \f}
    *
    * Returns an identity group with given shape.
    * 
    * @param[in] shape The shape of the group.
    * @param[out] out The CompositeManifold element. 
    */

    CompositeManifold out;
    out.space.push_back(glc::from_shape(shape));
    return out;
}

CompositeManifold::CompositeManifold(std::initializer_list<TYPES> others)
{
    /*!
    * Initializer list constructor for CompositeManifold
    *
    * Enables construction like:
    *     Lielab::domain::RN R(3);
    *     Lielab::domain::SO O(3);
    *     Lielab::domain::CompositeManifold M{R, O};
    */

    this->space = std::vector<TYPES>{std::move(others)};
}

CompositeManifold::CompositeManifold(const std::vector<TYPES>& others)
{
    /*!
    * Vector list constructor for CompositeManifold
    *
    * Not needed for C++, but enables construction in Python like:
    *     R = lielab.domain.RN(3)
    *     O = lielab.domain.SO(3)
    *     M = lielab.domain.CompositeManifold([R, O])
    */

    this->space = others;
}

size_t CompositeManifold::get_dimension() const
{
    const std::vector<size_t> dims = this->get_dimensions();
    return std::accumulate(dims.begin(), dims.end(), size_t(0));
}

std::vector<size_t> CompositeManifold::get_dimensions() const
{
    std::vector<size_t> dims(this->space.size());

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeManifold::INDEX_CN)
        {
            dims[ii] = std::get<CN>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_GLR)
        {
            dims[ii] = std::get<GLR>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_GLC)
        {
            dims[ii] = std::get<GLC>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_RN)
        {
            dims[ii] = std::get<RN>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_SE)
        {
            dims[ii] = std::get<SE>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_SO)
        {
            dims[ii] = std::get<SO>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_SP)
        {
            dims[ii] = std::get<SP>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_SU)
        {
            dims[ii] = std::get<SU>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_cn)
        {
            dims[ii] = std::get<cn>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_glr)
        {
            dims[ii] = std::get<glr>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_glc)
        {
            dims[ii] = std::get<glc>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_rn)
        {
            dims[ii] = std::get<rn>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_se)
        {
            dims[ii] = std::get<se>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_so)
        {
            dims[ii] = std::get<so>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_sp)
        {
            dims[ii] = std::get<sp>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_su)
        {
            dims[ii] = std::get<su>(this->space[ii]).get_dimension();
        }
    }

    return dims;
}

size_t CompositeManifold::get_size() const
{
    const std::vector<size_t> sizes = this->get_sizes();
    return std::accumulate(sizes.begin(), sizes.end(), size_t(0));
}

std::vector<size_t> CompositeManifold::get_sizes() const
{
    std::vector<size_t> sizes(this->space.size());

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeManifold::INDEX_CN)
        {
            sizes[ii] = std::get<CN>(this->space[ii]).get_size();
        }
        else if (ind == CompositeManifold::INDEX_GLR)
        {
            sizes[ii] = std::get<GLR>(this->space[ii]).get_size();
        }
        else if (ind == CompositeManifold::INDEX_GLC)
        {
            sizes[ii] = std::get<GLC>(this->space[ii]).get_size();
        }
        else if (ind == CompositeManifold::INDEX_RN)
        {
            sizes[ii] = std::get<RN>(this->space[ii]).get_size();
        }
        else if (ind == CompositeManifold::INDEX_SE)
        {
            sizes[ii] = std::get<SE>(this->space[ii]).get_size();
        }
        else if (ind == CompositeManifold::INDEX_SO)
        {
            sizes[ii] = std::get<SO>(this->space[ii]).get_size();
        }
        else if (ind == CompositeManifold::INDEX_SP)
        {
            sizes[ii] = std::get<SP>(this->space[ii]).get_size();
        }
        else if (ind == CompositeManifold::INDEX_SU)
        {
            sizes[ii] = std::get<SU>(this->space[ii]).get_size();
        }
        else if (ind == CompositeManifold::INDEX_cn)
        {
            sizes[ii] = std::get<cn>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_glr)
        {
            sizes[ii] = std::get<glr>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_glc)
        {
            sizes[ii] = std::get<glc>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_rn)
        {
            sizes[ii] = std::get<rn>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_se)
        {
            sizes[ii] = std::get<se>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_so)
        {
            sizes[ii] = std::get<so>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_sp)
        {
            sizes[ii] = std::get<sp>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeManifold::INDEX_su)
        {
            sizes[ii] = std::get<su>(this->space[ii]).get_dimension();
        }
    }

    return sizes;
}

size_t CompositeManifold::get_shape() const
{
    const std::vector<size_t> shapes = this->get_shapes();
    return std::accumulate(shapes.begin(), shapes.end(), size_t(0));
}

std::vector<size_t> CompositeManifold::get_shapes() const
{
    std::vector<size_t> out(space.size());

    for (size_t ii = 0; ii < space.size(); ii++)
    {
        const size_t ind = space[ii].index();
        if (ind == INDEX_CN)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::CN>(space[ii]).get_shape());
        }
        else if (ind == INDEX_GLR)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::GLR>(space[ii]).get_shape());
        }
        else if (ind == INDEX_GLC)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::GLC>(space[ii]).get_shape());
        }
        else if (ind == INDEX_RN)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::RN>(space[ii]).get_shape());
        }
        else if (ind == INDEX_SE)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::SE>(space[ii]).get_shape());
        }
        else if (ind == INDEX_SO)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::SO>(space[ii]).get_shape());
        }
        else if (ind == INDEX_SP)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::SP>(space[ii]).get_shape());
        }
        else if (ind == INDEX_SU)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::SU>(space[ii]).get_shape());
        }
        else if (ind == INDEX_cn)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::cn>(space[ii]).get_shape());
        }
        else if (ind == INDEX_glr)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::glr>(space[ii]).get_shape());
        }
        else if (ind == INDEX_glc)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::glc>(space[ii]).get_shape());
        }
        else if (ind == INDEX_rn)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::rn>(space[ii]).get_shape());
        }
        else if (ind == INDEX_se)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::se>(space[ii]).get_shape());
        }
        else if (ind == INDEX_so)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::so>(space[ii]).get_shape());
        }
        else if (ind == INDEX_sp)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::sp>(space[ii]).get_shape());
        }
        else if (ind == INDEX_su)
        {
            out[ii] = static_cast<int>(std::get<Lielab::domain::su>(space[ii]).get_shape());
        }
    }

    return out;
}

Eigen::VectorXd CompositeManifold::serialize() const
{

    if (space.size() == 0)
    {
        return Eigen::VectorXd(0);
    }
    
    std::vector<Eigen::VectorXd> serials;
    const size_t sz = this->space.size();

    for (size_t ii = 0; ii < sz; ii++)
    {
        const size_t ind = space[ii].index();

        if (ind == INDEX_CN)
        {
            serials.push_back(std::get<Lielab::domain::CN>(space[ii]).serialize());
        }
        else if (ind == INDEX_GLR)
        {
            serials.push_back(std::get<Lielab::domain::GLR>(space[ii]).serialize());
        }
        else if (ind == INDEX_GLC)
        {
            serials.push_back(std::get<Lielab::domain::GLC>(space[ii]).serialize());
        }
        else if (ind == INDEX_RN)
        {
            serials.push_back(std::get<Lielab::domain::RN>(space[ii]).serialize());
        }
        else if (ind == INDEX_SE)
        {
            serials.push_back(std::get<Lielab::domain::SE>(space[ii]).serialize());
        }
        else if (ind == INDEX_SO)
        {
            serials.push_back(std::get<Lielab::domain::SO>(space[ii]).serialize());
        }
        else if (ind == INDEX_SP)
        {
            serials.push_back(std::get<Lielab::domain::SP>(space[ii]).serialize());
        }
        else if (ind == INDEX_SU)
        {
            serials.push_back(std::get<Lielab::domain::SU>(space[ii]).serialize());
        }
        else if (ind == INDEX_cn)
        {
            serials.push_back(std::get<Lielab::domain::cn>(space[ii]).get_vector());
        }
        else if (ind == INDEX_glr)
        {
            serials.push_back(std::get<Lielab::domain::glr>(space[ii]).get_vector());
        }
        else if (ind == INDEX_glc)
        {
            serials.push_back(std::get<Lielab::domain::glc>(space[ii]).get_vector());
        }
        else if (ind == INDEX_rn)
        {
            serials.push_back(std::get<Lielab::domain::rn>(space[ii]).get_vector());
        }
        else if (ind == INDEX_se)
        {
            serials.push_back(std::get<Lielab::domain::se>(space[ii]).get_vector());
        }
        else if (ind == INDEX_so)
        {
            serials.push_back(std::get<Lielab::domain::so>(space[ii]).get_vector());
        }
        else if (ind == INDEX_sp)
        {
            serials.push_back(std::get<Lielab::domain::sp>(space[ii]).get_vector());
        }
        else if (ind == INDEX_su)
        {
            serials.push_back(std::get<Lielab::domain::su>(space[ii]).get_vector());
        }
    }
    
    Eigen::VectorXd out = Lielab::utils::concatenate(serials);
    return out;
}

void CompositeManifold::unserialize(const Eigen::VectorXd &vec)
{
    size_t jj = 0;
    size_t sz = 0;

    for (auto & M : this->space)
    {
        const size_t ind = M.index();
        jj += sz;

        if (ind == INDEX_CN)
        {
            sz = std::get<Lielab::domain::CN>(M).get_size();
            std::get<Lielab::domain::CN>(M).unserialize(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_GLR)
        {
            sz = std::get<Lielab::domain::GLR>(M).get_size();
            std::get<Lielab::domain::GLR>(M).unserialize(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_GLC)
        {
            sz = std::get<Lielab::domain::GLC>(M).get_size();
            std::get<Lielab::domain::GLC>(M).unserialize(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_RN)
        {
            sz = std::get<Lielab::domain::RN>(M).get_size();
            std::get<Lielab::domain::RN>(M).unserialize(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_SE)
        {
            sz = std::get<Lielab::domain::SE>(M).get_size();
            std::get<Lielab::domain::SE>(M).unserialize(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_SO)
        {
            sz = std::get<Lielab::domain::SO>(M).get_size();
            std::get<Lielab::domain::SO>(M).unserialize(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_SP)
        {
            sz = std::get<Lielab::domain::SP>(M).get_size();
            std::get<Lielab::domain::SP>(M).unserialize(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_SU)
        {
            sz = std::get<Lielab::domain::SU>(M).get_size();
            std::get<Lielab::domain::SU>(M).unserialize(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_cn)
        {
            sz = std::get<Lielab::domain::cn>(M).get_dimension();
            std::get<Lielab::domain::cn>(M).set_vector(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_glr)
        {
            sz = std::get<Lielab::domain::glr>(M).get_dimension();
            std::get<Lielab::domain::glr>(M).set_vector(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_glc)
        {
            sz = std::get<Lielab::domain::glc>(M).get_dimension();
            std::get<Lielab::domain::glc>(M).set_vector(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_rn)
        {
            sz = std::get<Lielab::domain::rn>(M).get_dimension();
            std::get<Lielab::domain::rn>(M).set_vector(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_se)
        {
            sz = std::get<Lielab::domain::se>(M).get_dimension();
            std::get<Lielab::domain::se>(M).set_vector(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_so)
        {
            sz = std::get<Lielab::domain::so>(M).get_dimension();
            std::get<Lielab::domain::so>(M).set_vector(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_sp)
        {
            sz = std::get<Lielab::domain::sp>(M).get_dimension();
            std::get<Lielab::domain::sp>(M).set_vector(vec(Eigen::seqN(jj, sz)));
        }
        else if (ind == INDEX_su)
        {
            sz = std::get<Lielab::domain::su>(M).get_dimension();
            std::get<Lielab::domain::su>(M).set_vector(vec(Eigen::seqN(jj, sz)));
        }
    }
}

void CompositeManifold::unserialize(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->unserialize(Eigen::VectorXd{std::move(vector)});
}

Eigen::MatrixXcd CompositeManifold::get_matrix() const
{
    const std::vector<size_t> shapes = this->get_shapes();
    const size_t shape = this->get_shape();

    Eigen::MatrixXcd out = Eigen::MatrixXcd::Zero(shape, shape);
    ptrdiff_t kk = 0;

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == INDEX_CN)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::CN>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_GLR)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::GLR>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_GLC)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::GLC>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_RN)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::RN>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_SE)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::SE>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_SO)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::SO>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_SP)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::SP>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_SU)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::SU>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_cn)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::cn>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_glr)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::glr>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_glc)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::glc>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_rn)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::rn>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_se)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::se>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_so)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::so>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_sp)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::sp>(this->space[ii]).get_matrix();
        }
        else if (ind == INDEX_su)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<Lielab::domain::su>(this->space[ii]).get_matrix();
        }
        kk += shapes[ii];
    }

    return out;
}

std::complex<double> CompositeManifold::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t shape = this->get_shape();
    if (shape == 0) return std::complex<double>(nan, nan);

    if (index1 >= static_cast<ptrdiff_t>(shape)) return std::complex<double>(nan, nan);
    if (index2 >= static_cast<ptrdiff_t>(shape)) return std::complex<double>(nan, nan);
    if (std::abs(index1) > static_cast<ptrdiff_t>(shape)) return std::complex<double>(nan, nan);
    if (std::abs(index2) > static_cast<ptrdiff_t>(shape)) return std::complex<double>(nan, nan);

    size_t _index1;
    if (index1 < 0)
    {
        _index1 = static_cast<size_t>(static_cast<ptrdiff_t>(shape) + index1);
    }
    else
    {
        _index1 = static_cast<size_t>(index1);
    }

    size_t _index2;
    if (index2 < 0)
    {
        _index2 = static_cast<size_t>(static_cast<ptrdiff_t>(shape) + index2);
    }
    else
    {
        _index2 = static_cast<size_t>(index2);
    }
    
    const std::vector<size_t> shapes = this->get_shapes();

    size_t relind = 0;
    for (size_t ii = 0; ii < shapes.size(); ii++)
    {
        if ((_index1 - relind) < shapes[ii] && (_index2 - relind) < shapes[ii])
        {
            const size_t ind = this->space[ii].index();
            if (ind == CompositeManifold::INDEX_CN)
            {
                return std::get<CN>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_GLR)
            {
                return std::get<GLR>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_GLC)
            {
                return std::get<GLC>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_RN)
            {
                return std::get<RN>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_SE)
            {
                return std::get<SE>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_SO)
            {
                return std::get<SO>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_SP)
            {
                return std::get<SP>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_SU)
            {
                return std::get<SU>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_cn)
            {
                return std::get<cn>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_glr)
            {
                return std::get<glr>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_glc)
            {
                return std::get<glc>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_rn)
            {
                return std::get<rn>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_se)
            {
                return std::get<se>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_so)
            {
                return std::get<so>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_sp)
            {
                return std::get<sp>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeManifold::INDEX_su)
            {
                return std::get<su>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
        }
        relind += shapes[ii];
    }

    return std::complex<double>(0.0, 0.0);
}

CompositeManifold::TYPES CompositeManifold::operator[](const ptrdiff_t index) const
{
    const size_t sz = this->space.size();
    if (index >= static_cast<ptrdiff_t>(sz)) return Lielab::domain::glc();

    if (index >= static_cast<ptrdiff_t>(sz)) return Lielab::domain::glc();
    if (std::abs(index) > static_cast<ptrdiff_t>(sz)) return Lielab::domain::glc();

    size_t _index;
    if (index < 0)
    {
        _index = static_cast<size_t>(static_cast<ptrdiff_t>(sz) + index);
    }
    else
    {
        _index = static_cast<size_t>(index);
    }
    
    const size_t ind = space[_index].index();
    if (ind == INDEX_CN)
    {
        return std::get<Lielab::domain::CN>(space[_index]);
    }
    else if (ind == INDEX_GLR)
    {
        return std::get<Lielab::domain::GLR>(space[_index]);
    }
    else if (ind == INDEX_GLC)
    {
        return std::get<Lielab::domain::GLC>(space[_index]);
    }
    else if (ind == INDEX_RN)
    {
        return std::get<Lielab::domain::RN>(space[_index]);
    }
    else if (ind == INDEX_SE)
    {
        return std::get<Lielab::domain::SE>(space[_index]);
    }
    else if (ind == INDEX_SO)
    {
        return std::get<Lielab::domain::SO>(space[_index]);
    }
    else if (ind == INDEX_SP)
    {
        return std::get<Lielab::domain::SP>(space[_index]);
    }
    else if (ind == INDEX_SU)
    {
        return std::get<Lielab::domain::SU>(space[_index]);
    }
    else if (ind == INDEX_cn)
    {
        return std::get<Lielab::domain::cn>(space[_index]);
    }
    else if (ind == INDEX_glr)
    {
        return std::get<Lielab::domain::glr>(space[_index]);
    }
    else if (ind == INDEX_glc)
    {
        return std::get<Lielab::domain::glc>(space[_index]);
    }
    else if (ind == INDEX_rn)
    {
        return std::get<Lielab::domain::rn>(space[_index]);
    }
    else if (ind == INDEX_se)
    {
        return std::get<Lielab::domain::se>(space[_index]);
    }
    else if (ind == INDEX_so)
    {
        return std::get<Lielab::domain::so>(space[_index]);
    }
    else if (ind == INDEX_sp)
    {
        return std::get<Lielab::domain::sp>(space[_index]);
    }
    else if (ind == INDEX_su)
    {
        return std::get<Lielab::domain::su>(space[_index]);
    }

    // This should never be called
    return Lielab::domain::glc();
}

std::ostream & operator<<(std::ostream & os, const CompositeManifold & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */
    
    os << other.to_string();
    return os;
}

}
