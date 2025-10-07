#include "CompositeGroup.hpp"

#include "CN.hpp"
#include "GLR.hpp"
#include "GLC.hpp"
#include "RN.hpp"
#include "SE.hpp"
#include "SO.hpp"
#include "SP.hpp"
#include "SU.hpp"

#include "Lielab/domain/liealgebras.hpp"
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

std::string CompositeGroup::to_string() const
{
    std::string out = "";
    const size_t sz = this->space.size();
    std::vector<size_t> shapes = this->get_shapes();

    for (size_t ii = 0; ii < sz; ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeGroup::INDEX_CN)
        {
            out += std::get<CN>(this->space[ii]).to_string();
        }
        else if (ind == CompositeGroup::INDEX_GLR)
        {
            out += std::get<GLR>(this->space[ii]).to_string();
        }
        else if (ind == CompositeGroup::INDEX_GLC)
        {
            out += std::get<GLC>(this->space[ii]).to_string();
        }
        else if (ind == CompositeGroup::INDEX_RN)
        {
            out += std::get<RN>(this->space[ii]).to_string();
        }
        else if (ind == CompositeGroup::INDEX_SE)
        {
            out += std::get<SE>(this->space[ii]).to_string();
        }
        else if (ind == CompositeGroup::INDEX_SO)
        {
            out += std::get<SO>(this->space[ii]).to_string();
        }
        else if (ind == CompositeGroup::INDEX_SP)
        {
            out += std::get<SP>(this->space[ii]).to_string();
        }
        else if (ind == CompositeGroup::INDEX_SU)
        {
            out += std::get<SU>(this->space[ii]).to_string();
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

CompositeGroup::CompositeGroup()
{
    
}

CompositeGroup::CompositeGroup(const size_t n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{CompositeGroup} \f}
    *
    * Constructor instantiating a \f$\mathfrak{CompositeGroup}\f$ as a single GLC object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::CompositeGroup x(3), y(4), z(5);
    * 
    * @param[in] n The shape of the data matrix.
    */

    this->space.push_back(GLC(n));
}

CompositeGroup CompositeGroup::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{CompositeGroup} \f}
    *
    * Returns an identity group with given shape.
    * 
    * @param[in] shape The shape of the group.
    * @param[out] out The CompositeGroup element. 
    */

    CompositeGroup out;
    out.space.push_back(GLC::from_shape(shape));
    return out;
}

CompositeGroup::CompositeGroup(std::initializer_list<TYPES> others)
{
    /*!
    * Initializer list constructor for CompositeGroup
    *
    * Enables construction like:
    *     Lielab::domain::RN R(3);
    *     Lielab::domain::SO O(3);
    *     Lielab::domain::CompositeGroup M{R, O};
    */

    this->space = std::vector<TYPES>{std::move(others)};
}

CompositeGroup::CompositeGroup(const std::vector<TYPES>& others)
{
    /*!
    * Vector list constructor for CompositeGroup
    *
    * Not needed for C++, but enables construction in Python like:
    *     R = lielab.domain.RN(3)
    *     O = lielab.domain.SO(3)
    *     M = lielab.domain.CompositeGroup([R, O])
    */

    this->space = others;
}

size_t CompositeGroup::get_dimension() const
{
    const std::vector<size_t> dims = this->get_dimensions();
    return std::accumulate(dims.begin(), dims.end(), size_t(0));
}

std::vector<size_t> CompositeGroup::get_dimensions() const
{
    std::vector<size_t> dims(this->space.size());

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeGroup::INDEX_CN)
        {
            dims[ii] = std::get<CN>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeGroup::INDEX_GLR)
        {
            dims[ii] = std::get<GLR>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeGroup::INDEX_GLC)
        {
            dims[ii] = std::get<GLC>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeGroup::INDEX_RN)
        {
            dims[ii] = std::get<RN>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeGroup::INDEX_SE)
        {
            dims[ii] = std::get<SE>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeGroup::INDEX_SO)
        {
            dims[ii] = std::get<SO>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeGroup::INDEX_SP)
        {
            dims[ii] = std::get<SP>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeGroup::INDEX_SU)
        {
            dims[ii] = std::get<SU>(this->space[ii]).get_dimension();
        }
    }

    return dims;
}

size_t CompositeGroup::get_size() const
{
    const std::vector<size_t> sizes = this->get_sizes();
    return std::accumulate(sizes.begin(), sizes.end(), size_t(0));
}

std::vector<size_t> CompositeGroup::get_sizes() const
{
    std::vector<size_t> sizes(this->space.size());

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeGroup::INDEX_CN)
        {
            sizes[ii] = std::get<CN>(this->space[ii]).get_size();
        }
        else if (ind == CompositeGroup::INDEX_GLR)
        {
            sizes[ii] = std::get<GLR>(this->space[ii]).get_size();
        }
        else if (ind == CompositeGroup::INDEX_GLC)
        {
            sizes[ii] = std::get<GLC>(this->space[ii]).get_size();
        }
        else if (ind == CompositeGroup::INDEX_RN)
        {
            sizes[ii] = std::get<RN>(this->space[ii]).get_size();
        }
        else if (ind == CompositeGroup::INDEX_SE)
        {
            sizes[ii] = std::get<SE>(this->space[ii]).get_size();
        }
        else if (ind == CompositeGroup::INDEX_SO)
        {
            sizes[ii] = std::get<SO>(this->space[ii]).get_size();
        }
        else if (ind == CompositeGroup::INDEX_SP)
        {
            sizes[ii] = std::get<SP>(this->space[ii]).get_size();
        }
        else if (ind == CompositeGroup::INDEX_SU)
        {
            sizes[ii] = std::get<SU>(this->space[ii]).get_size();
        }
    }

    return sizes;
}

size_t CompositeGroup::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    const std::vector<size_t> shapes = this->get_shapes();
    return std::accumulate(shapes.begin(), shapes.end(), size_t(0));
}

std::vector<size_t> CompositeGroup::get_shapes() const
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
    }

    return out;
}

Eigen::VectorXd CompositeGroup::serialize() const
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
    }
    
    Eigen::VectorXd out = Lielab::utils::concatenate(serials);
    return out;
}

void CompositeGroup::unserialize(const Eigen::VectorXd &vec)
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
    }
}

void CompositeGroup::unserialize(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->unserialize(Eigen::VectorXd{std::move(vector)});
}

Eigen::MatrixXcd CompositeGroup::get_matrix() const
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
        kk += shapes[ii];
    }

    return out;
}

// double CompositeGroup::operator()(const ptrdiff_t index) const
// {
//     /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
//     *
//     * Gets a value in the vector representation.
//     */
    
//     const double nan = std::numeric_limits<double>::quiet_NaN();
//     const size_t dim = this->get_dimension();

//     if (index >= static_cast<ptrdiff_t>(dim)) return nan;
//     if (std::abs(index) > static_cast<ptrdiff_t>(dim)) return nan;

//     size_t _index;
//     if (index < 0)
//     {
//         _index = static_cast<size_t>(static_cast<ptrdiff_t>(dim) + index);
//     }
//     else
//     {
//         _index = static_cast<size_t>(index);
//     }
    
//     const std::vector<size_t> sizes = this->get_sizes();

//     size_t base_index = 0;
//     for (size_t ii = 0; ii < sizes.size(); ii++)
//     {
//         if ((_index - base_index) < sizes[ii])
//         {
//             const size_t ind = this->space[ii].index();
//             const size_t relind = _index - base_index;
//             if (ind == CompositeGroup::INDEX_CN)
//             {
//                 return std::get<CN>(this->space[ii]).operator()(relind);
//             }
//             else if (ind == CompositeGroup::INDEX_GLR)
//             {
//                 return std::get<GLR>(this->space[ii]).operator()(relind);
//             }
//             else if (ind == CompositeGroup::INDEX_GLC)
//             {
//                 return std::get<GLC>(this->space[ii]).operator()(relind);
//             }
//             else if (ind == CompositeGroup::INDEX_RN)
//             {
//                 return std::get<RN>(this->space[ii]).operator()(relind);
//             }
//             else if (ind == CompositeGroup::INDEX_SE)
//             {
//                 return std::get<SE>(this->space[ii]).operator()(relind);
//             }
//             else if (ind == CompositeGroup::INDEX_SO)
//             {
//                 return std::get<SO>(this->space[ii]).operator()(relind);
//             }
//             else if (ind == CompositeGroup::INDEX_SP)
//             {
//                 return std::get<SP>(this->space[ii]).operator()(relind);
//             }
//             else if (ind == CompositeGroup::INDEX_SU)
//             {
//                 return std::get<SU>(this->space[ii]).operator()(relind);
//             }
//         }
//         base_index += sizes[ii];
//     }

//     // This should never be returned.
//     return nan;
// }

std::complex<double> CompositeGroup::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
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
            if (ind == CompositeGroup::INDEX_CN)
            {
                return std::get<CN>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeGroup::INDEX_GLR)
            {
                return std::get<GLR>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeGroup::INDEX_GLC)
            {
                return std::get<GLC>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeGroup::INDEX_RN)
            {
                return std::get<RN>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeGroup::INDEX_SE)
            {
                return std::get<SE>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeGroup::INDEX_SO)
            {
                return std::get<SO>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeGroup::INDEX_SP)
            {
                return std::get<SP>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeGroup::INDEX_SU)
            {
                return std::get<SU>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
        }
        relind += shapes[ii];
    }

    return std::complex<double>(0.0, 0.0);
}

CompositeGroup CompositeGroup::operator*(const CompositeGroup & other) const
{
    CompositeGroup out;

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == INDEX_CN)
        {
            out.space.push_back(std::get<Lielab::domain::CN>(this->space[ii]) * std::get<Lielab::domain::CN>(other.space[ii]));
        }
        else if (ind == INDEX_GLR)
        {
            out.space.push_back(std::get<Lielab::domain::GLR>(this->space[ii]) * std::get<Lielab::domain::GLR>(other.space[ii]));
        }
        else if (ind == INDEX_GLC)
        {
            out.space.push_back(std::get<Lielab::domain::GLC>(this->space[ii]) * std::get<Lielab::domain::GLC>(other.space[ii]));
        }
        else if (ind == INDEX_RN)
        {
            out.space.push_back(std::get<Lielab::domain::RN>(this->space[ii]) * std::get<Lielab::domain::RN>(other.space[ii]));
        }
        else if (ind == INDEX_SE)
        {
            out.space.push_back(std::get<Lielab::domain::SE>(this->space[ii]) * std::get<Lielab::domain::SE>(other.space[ii]));
        }
        else if (ind == INDEX_SO)
        {
            out.space.push_back(std::get<Lielab::domain::SO>(this->space[ii]) * std::get<Lielab::domain::SO>(other.space[ii]));
        }
        else if (ind == INDEX_SP)
        {
            out.space.push_back(std::get<Lielab::domain::SP>(this->space[ii]) * std::get<Lielab::domain::SP>(other.space[ii]));
        }
        else if (ind == INDEX_SU)
        {
            out.space.push_back(std::get<Lielab::domain::SU>(this->space[ii]) * std::get<Lielab::domain::SU>(other.space[ii]));
        }
    }

    return out;
}

CompositeGroup & CompositeGroup::operator*=(const CompositeGroup & other)
{
    assert(this->space.size() == other.space.size());

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == INDEX_CN)
        {
            std::get<Lielab::domain::CN>(this->space[ii]) *= std::get<Lielab::domain::CN>(other.space[ii]);
        }
        else if (ind == INDEX_GLR)
        {
            std::get<Lielab::domain::GLR>(this->space[ii]) *= std::get<Lielab::domain::GLR>(other.space[ii]);
        }
        else if (ind == INDEX_GLC)
        {
            std::get<Lielab::domain::GLC>(this->space[ii]) *= std::get<Lielab::domain::GLC>(other.space[ii]);
        }
        else if (ind == INDEX_RN)
        {
            std::get<Lielab::domain::RN>(this->space[ii]) *= std::get<Lielab::domain::RN>(other.space[ii]);
        }
        else if (ind == INDEX_SE)
        {
            std::get<Lielab::domain::SE>(this->space[ii]) *= std::get<Lielab::domain::SE>(other.space[ii]);
        }
        else if (ind == INDEX_SO)
        {
            std::get<Lielab::domain::SO>(this->space[ii]) *= std::get<Lielab::domain::SO>(other.space[ii]);
        }
        else if (ind == INDEX_SP)
        {
            std::get<Lielab::domain::SP>(this->space[ii]) *= std::get<Lielab::domain::SP>(other.space[ii]);
        }
        else if (ind == INDEX_SU)
        {
            std::get<Lielab::domain::SU>(this->space[ii]) *= std::get<Lielab::domain::SU>(other.space[ii]);
        }
    }

    return *this;
}

CompositeGroup CompositeGroup::inverse() const
{
    CompositeGroup out;

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == INDEX_CN)
        {
            out.space.push_back(std::get<Lielab::domain::CN>(this->space[ii]).inverse());
        }
        if (ind == INDEX_GLR)
        {
            out.space.push_back(std::get<Lielab::domain::GLR>(this->space[ii]).inverse());
        }
        if (ind == INDEX_GLC)
        {
            out.space.push_back(std::get<Lielab::domain::GLC>(this->space[ii]).inverse());
        }
        else if (ind == INDEX_RN)
        {
            out.space.push_back(std::get<Lielab::domain::RN>(this->space[ii]).inverse());
        }
        else if (ind == INDEX_SE)
        {
            out.space.push_back(std::get<Lielab::domain::SE>(this->space[ii]).inverse());
        }
        else if (ind == INDEX_SO)
        {
            out.space.push_back(std::get<Lielab::domain::SO>(this->space[ii]).inverse());
        }
        else if (ind == INDEX_SP)
        {
            out.space.push_back(std::get<Lielab::domain::SP>(this->space[ii]).inverse());
        }
        else if (ind == INDEX_SU)
        {
            out.space.push_back(std::get<Lielab::domain::SU>(this->space[ii]).inverse());
        }
    }

    return out;
}

CompositeGroup::TYPES CompositeGroup::operator[](const ptrdiff_t index) const
{
    const size_t sz = this->space.size();
    if (index >= static_cast<ptrdiff_t>(sz)) return Lielab::domain::GLC();

    if (std::abs(index) > static_cast<ptrdiff_t>(sz)) return Lielab::domain::GLC();

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

    // This should never be called.
    return Lielab::domain::GLC();
}

std::ostream & operator<<(std::ostream & os, const CompositeGroup & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */
    
    os << other.to_string();
    return os;
}

}
