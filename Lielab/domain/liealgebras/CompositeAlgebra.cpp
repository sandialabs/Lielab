#include "CompositeAlgebra.hpp"

#include "cn.hpp"
#include "glr.hpp"
#include "glc.hpp"
#include "rn.hpp"
#include "se.hpp"
#include "so.hpp"
#include "sp.hpp"
#include "su.hpp"

#include <Eigen/Core>

#include <initializer_list>
#include <numeric>

namespace Lielab::domain
{

std::string CompositeAlgebra::to_string() const
{
    std::string out = "";
    const size_t sz = this->space.size();
    std::vector<size_t> shapes = this->get_shapes();

    for (size_t ii = 0; ii < sz; ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            out += std::get<cn>(this->space[ii]).to_string();
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            out += std::get<glr>(this->space[ii]).to_string();
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            out += std::get<glc>(this->space[ii]).to_string();
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            out += std::get<rn>(this->space[ii]).to_string();
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            out += std::get<se>(this->space[ii]).to_string();
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            out += std::get<so>(this->space[ii]).to_string();
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            out += std::get<sp>(this->space[ii]).to_string();
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            out += std::get<su>(this->space[ii]).to_string();
        }

        if (ii < sz-1)
        {
            out += " ";
            out += "⊕";
            out += " ";
        }
    }

    return out;
}

CompositeAlgebra::CompositeAlgebra()
{
    /*! \f{equation*}{() \rightarrow \mathfrak{CompositeAlgebra} \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::CompositeAlgebra x, y, z;
    * 
    */
}

CompositeAlgebra::CompositeAlgebra(const size_t n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{CompositeAlgebra} \f}
    *
    * Constructor instantiating a \f$\mathfrak{CompositeAlgebra}\f$ as a single glc object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::CompositeAlgebra x(3), y(4), z(5);
    * 
    * @param[in] n The shape of the data matrix.
    */

    this->space.push_back(glc(n));
}

CompositeAlgebra CompositeAlgebra::basis(const ptrdiff_t i, const size_t n)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{CompositeAlgebra} \f}
    *
    * Returns the i'th basis element of the CompositeAlgebra algebra.
    * 
    * @param[in] i The basis vector.
    * @param[in] n The shape of the algebra.
    * @param[out] out The CompositeAlgebra element.
    */

    CompositeAlgebra out;
    out.space.push_back(glc::basis(i, n));
    return out;
}

CompositeAlgebra CompositeAlgebra::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{CompositeAlgebra} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The CompositeAlgebra element. 
    */

    CompositeAlgebra out;
    out.space.push_back(glc::from_shape(shape));
    return out;
}

CompositeAlgebra::CompositeAlgebra(std::initializer_list<CompositeAlgebra::TYPES> others)
{
    /*!
    * Initializer list constructor for CompositeAlgebra
    *
    * Enables construction like:
    *     Lielab::domain::rn R(3);
    *     Lielab::domain::so O(3);
    *     Lielab::domain::CompositeAlgebra M{R, O};
    */

    this->space = std::vector<CompositeAlgebra::TYPES>{std::move(others)};
}

CompositeAlgebra::CompositeAlgebra(const std::vector<CompositeAlgebra::TYPES>& others)
{
    /*!
    * Vector list constructor for CompositeAlgebra
    *
    * Not needed for C++, but enables construction in Python like:
    *     R = lielab.domain.rn(3)
    *     O = lielab.domain.so(3)
    *     M = lielab.domain.CompositeAlgebra([R, O])
    */

    this->space = others;
}

size_t CompositeAlgebra::get_dimension() const
{
    const std::vector<size_t> dims = this->get_dimensions();
    return std::accumulate(dims.begin(), dims.end(), size_t(0));
}

Eigen::VectorXd CompositeAlgebra::get_vector() const
{
    return Lielab::utils::concatenate(this->get_vectors());
}

void CompositeAlgebra::set_vector(const Eigen::VectorXd& vec)
{
    size_t jj = 0;

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            const size_t dim = std::get<cn>(this->space[ii]).get_dimension();
            std::get<cn>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
            jj += dim;
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            const size_t dim = std::get<glr>(this->space[ii]).get_dimension();
            std::get<glr>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
            jj += dim;
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            const size_t dim = std::get<glc>(this->space[ii]).get_dimension();
            std::get<glc>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
            jj += dim;
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            const size_t dim = std::get<rn>(this->space[ii]).get_dimension();
            std::get<rn>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
            jj += dim;
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            const size_t dim = std::get<se>(this->space[ii]).get_dimension();
            std::get<se>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
            jj += dim;
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            const size_t dim = std::get<so>(this->space[ii]).get_dimension();
            std::get<so>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
            jj += dim;
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            const size_t dim = std::get<sp>(this->space[ii]).get_dimension();
            std::get<sp>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
            jj += dim;
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            const size_t dim = std::get<su>(this->space[ii]).get_dimension();
            std::get<su>(this->space[ii]).set_vector(vec(Eigen::seqN(jj, dim)));
            jj += dim;
        }
    }
}

void CompositeAlgebra::set_vector(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->set_vector(Eigen::VectorXd{std::move(vector)});
}

CompositeAlgebra::matrix_t CompositeAlgebra::get_matrix() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times n} \f}
    * 
    * Returns a matrix representation.
    * 
    * Formerly called "get_ados_representation()".
    * 
    * Ado, Igor D. "Note on the representation of finite continuous groups by
    *               means of linear substitutions, Izv. Fiz." Mat. Obsch.(Kazan)
    *               7.1 (1935): 935.
    * 
    * Ado, Igor D. "The representation of Lie algebras by matrices." Uspekhi
    *               Matematicheskikh Nauk 2.6 (1947): 159-173.
    */
    
    const size_t shape = this->get_shape();
    const std::vector<size_t> shapes = this->get_shapes();

    Eigen::MatrixXcd out = Eigen::MatrixXcd::Zero(shape, shape);
    ptrdiff_t kk = 0;

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<cn>(this->space[ii]).get_matrix();
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<glr>(this->space[ii]).get_matrix();
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<glc>(this->space[ii]).get_matrix();
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<rn>(this->space[ii]).get_matrix();
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<se>(this->space[ii]).get_matrix();
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<so>(this->space[ii]).get_matrix();
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<sp>(this->space[ii]).get_matrix();
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            out(Eigen::seqN(kk, shapes[ii]), Eigen::seqN(kk, shapes[ii])) = std::get<su>(this->space[ii]).get_matrix();
        }
        kk += shapes[ii];
    }

    return out;
}

size_t CompositeAlgebra::get_shape() const
{
    const std::vector<size_t> shapes = this->get_shapes();
    return std::accumulate(shapes.begin(), shapes.end(), size_t(0));
}

std::vector<size_t> CompositeAlgebra::get_dimensions() const
{
    std::vector<size_t> dims(this->space.size());

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            dims[ii] = std::get<cn>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            dims[ii] = std::get<glr>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            dims[ii] = std::get<glc>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            dims[ii] = std::get<rn>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            dims[ii] = std::get<se>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            dims[ii] = std::get<so>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            dims[ii] = std::get<sp>(this->space[ii]).get_dimension();
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            dims[ii] = std::get<su>(this->space[ii]).get_dimension();
        }
    }

    return dims;
}

std::vector<Eigen::VectorXd> CompositeAlgebra::get_vectors() const
{
    std::vector<Eigen::VectorXd> vectors = std::vector<Eigen::VectorXd>(this->space.size());

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            vectors[ii] = std::get<cn>(this->space[ii]).get_vector();
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            vectors[ii] = std::get<glr>(this->space[ii]).get_vector();
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            vectors[ii] = std::get<glc>(this->space[ii]).get_vector();
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            vectors[ii] = std::get<rn>(this->space[ii]).get_vector();
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            vectors[ii] = std::get<se>(this->space[ii]).get_vector();
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            vectors[ii] = std::get<so>(this->space[ii]).get_vector();
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            vectors[ii] = std::get<sp>(this->space[ii]).get_vector();
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            vectors[ii] = std::get<su>(this->space[ii]).get_vector();
        }
    }

    return vectors;
}

std::vector<size_t> CompositeAlgebra::get_shapes() const
{
    std::vector<size_t> out(this->space.size());

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == INDEX_cn)
        {
            out[ii] = static_cast<int>(std::get<cn>(this->space[ii]).get_shape());
        }
        else if (ind == INDEX_glr)
        {
            out[ii] = static_cast<int>(std::get<glr>(this->space[ii]).get_shape());
        }
        else if (ind == INDEX_glc)
        {
            out[ii] = static_cast<int>(std::get<glc>(this->space[ii]).get_shape());
        }
        else if (ind == INDEX_rn)
        {
            out[ii] = static_cast<int>(std::get<rn>(this->space[ii]).get_shape());
        }
        else if (ind == INDEX_se)
        {
            out[ii] = static_cast<int>(std::get<se>(this->space[ii]).get_shape());
        }
        else if (ind == INDEX_so)
        {
            out[ii] = static_cast<int>(std::get<so>(this->space[ii]).get_shape());
        }
        else if (ind == INDEX_sp)
        {
            out[ii] = static_cast<int>(std::get<sp>(this->space[ii]).get_shape());
        }
        else if (ind == INDEX_su)
        {
            out[ii] = static_cast<int>(std::get<su>(this->space[ii]).get_shape());
        }
    }

    return out;
}

double CompositeAlgebra::operator()(const ptrdiff_t index) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the vector representation.
    */
    
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t dim = this->get_dimension();

    if (index >= static_cast<ptrdiff_t>(dim)) return nan;
    if (std::abs(index) > static_cast<ptrdiff_t>(dim)) return nan;

    size_t _index;
    if (index < 0)
    {
        _index = static_cast<size_t>(static_cast<ptrdiff_t>(dim) + index);
    }
    else
    {
        _index = static_cast<size_t>(index);
    }
    
    const std::vector<size_t> dims = this->get_dimensions();

    size_t base_index = 0;
    for (size_t ii = 0; ii < dims.size(); ii++)
    {
        if ((_index - base_index) < dims[ii])
        {
            const size_t ind = this->space[ii].index();
            const size_t relind = _index - base_index;
            if (ind == CompositeAlgebra::INDEX_cn)
            {
                return std::get<cn>(this->space[ii]).operator()(relind);
            }
            else if (ind == CompositeAlgebra::INDEX_glr)
            {
                return std::get<glr>(this->space[ii]).operator()(relind);
            }
            else if (ind == CompositeAlgebra::INDEX_glc)
            {
                return std::get<glc>(this->space[ii]).operator()(relind);
            }
            else if (ind == CompositeAlgebra::INDEX_rn)
            {
                return std::get<rn>(this->space[ii]).operator()(relind);
            }
            else if (ind == CompositeAlgebra::INDEX_se)
            {
                return std::get<se>(this->space[ii]).operator()(relind);
            }
            else if (ind == CompositeAlgebra::INDEX_so)
            {
                return std::get<so>(this->space[ii]).operator()(relind);
            }
            else if (ind == CompositeAlgebra::INDEX_sp)
            {
                return std::get<sp>(this->space[ii]).operator()(relind);
            }
            else if (ind == CompositeAlgebra::INDEX_su)
            {
                return std::get<su>(this->space[ii]).operator()(relind);
            }
        }
        base_index += dims[ii];
    }

    // This should never be returned.
    return nan;
}

std::complex<double> CompositeAlgebra::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
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
            if (ind == CompositeAlgebra::INDEX_cn)
            {
                return std::get<cn>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeAlgebra::INDEX_glr)
            {
                return std::get<glr>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeAlgebra::INDEX_glc)
            {
                return std::get<glc>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeAlgebra::INDEX_rn)
            {
                return std::get<rn>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeAlgebra::INDEX_se)
            {
                return std::get<se>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeAlgebra::INDEX_so)
            {
                return std::get<so>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeAlgebra::INDEX_sp)
            {
                return std::get<sp>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
            else if (ind == CompositeAlgebra::INDEX_su)
            {
                return std::get<su>(this->space[ii]).operator()(_index1 - relind, _index2 - relind);
            }
        }
        relind += shapes[ii];
    }

    return std::complex<double>(0.0, 0.0);
}

CompositeAlgebra CompositeAlgebra::operator+(const CompositeAlgebra & other) const
{
    assert(this->space.size() == other.space.size());

    CompositeAlgebra out;

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(std::get<cn>(this->space[ii]) + std::get<cn>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(std::get<glr>(this->space[ii]) + std::get<glr>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(std::get<glc>(this->space[ii]) + std::get<glc>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(std::get<rn>(this->space[ii]) + std::get<rn>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(std::get<se>(this->space[ii]) + std::get<se>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(std::get<so>(this->space[ii]) + std::get<so>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(std::get<sp>(this->space[ii]) + std::get<sp>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(std::get<su>(this->space[ii]) + std::get<su>(other.space[ii]));
        }
    }

    return out;
}

CompositeAlgebra & CompositeAlgebra::operator+=(const CompositeAlgebra & other)
{
    assert(this->space.size() == other.space.size());

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            std::get<cn>(this->space[ii]) += std::get<cn>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            std::get<glr>(this->space[ii]) += std::get<glr>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            std::get<glc>(this->space[ii]) += std::get<glc>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            std::get<rn>(this->space[ii]) += std::get<rn>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            std::get<se>(this->space[ii]) += std::get<se>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            std::get<so>(this->space[ii]) += std::get<so>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            std::get<sp>(this->space[ii]) += std::get<sp>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            std::get<su>(this->space[ii]) += std::get<su>(other.space[ii]);
        }
    }

    return *this;
}

CompositeAlgebra CompositeAlgebra::operator-(const CompositeAlgebra & other) const
{
    CompositeAlgebra out;

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(std::get<cn>(this->space[ii]) - std::get<cn>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(std::get<glr>(this->space[ii]) - std::get<glr>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(std::get<glc>(this->space[ii]) - std::get<glc>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(std::get<rn>(this->space[ii]) - std::get<rn>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(std::get<se>(this->space[ii]) - std::get<se>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(std::get<so>(this->space[ii]) - std::get<so>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(std::get<sp>(this->space[ii]) - std::get<sp>(other.space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(std::get<su>(this->space[ii]) - std::get<su>(other.space[ii]));
        }
    }

    return out;
}

CompositeAlgebra & CompositeAlgebra::operator-=(const CompositeAlgebra & other)
{
    assert(this->space.size() == other.space.size());

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            std::get<cn>(this->space[ii]) -= std::get<cn>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            std::get<glr>(this->space[ii]) -= std::get<glr>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            std::get<glc>(this->space[ii]) -= std::get<glc>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            std::get<rn>(this->space[ii]) -= std::get<rn>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            std::get<se>(this->space[ii]) -= std::get<se>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            std::get<so>(this->space[ii]) -= std::get<so>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            std::get<sp>(this->space[ii]) -= std::get<sp>(other.space[ii]);
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            std::get<su>(this->space[ii]) -= std::get<su>(other.space[ii]);
        }
    }

    return *this;
}

CompositeAlgebra CompositeAlgebra::operator-() const
{
    CompositeAlgebra out;

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(-std::get<cn>(this->space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(-std::get<glr>(this->space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(-std::get<glc>(this->space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(-std::get<rn>(this->space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(-std::get<se>(this->space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(-std::get<so>(this->space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(-std::get<sp>(this->space[ii]));
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(-std::get<su>(this->space[ii]));
        }
    }

    return out;
}

CompositeAlgebra CompositeAlgebra::operator*(const double other) const
{
    CompositeAlgebra out;

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            out.space.push_back(std::get<cn>(this->space[ii]) * other);
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            out.space.push_back(std::get<glr>(this->space[ii]) * other);
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            out.space.push_back(std::get<glc>(this->space[ii]) * other);
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            out.space.push_back(std::get<rn>(this->space[ii]) * other);
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            out.space.push_back(std::get<se>(this->space[ii]) * other);
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            out.space.push_back(std::get<so>(this->space[ii]) * other);
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            out.space.push_back(std::get<sp>(this->space[ii]) * other);
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            out.space.push_back(std::get<su>(this->space[ii]) * other);
        }
    }

    return out;
}

CompositeAlgebra operator*(const double other, const CompositeAlgebra & rhs)
{
    return rhs*other;
}

CompositeAlgebra & CompositeAlgebra::operator*=(const double other)
{
    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == CompositeAlgebra::INDEX_cn)
        {
            std::get<cn>(this->space[ii]) *= other;
        }
        else if (ind == CompositeAlgebra::INDEX_glr)
        {
            std::get<glr>(this->space[ii]) *= other;
        }
        else if (ind == CompositeAlgebra::INDEX_glc)
        {
            std::get<glc>(this->space[ii]) *= other;
        }
        else if (ind == CompositeAlgebra::INDEX_rn)
        {
            std::get<rn>(this->space[ii]) *= other;
        }
        else if (ind == CompositeAlgebra::INDEX_se)
        {
            std::get<se>(this->space[ii]) *= other;
        }
        else if (ind == CompositeAlgebra::INDEX_so)
        {
            std::get<so>(this->space[ii]) *= other;
        }
        else if (ind == CompositeAlgebra::INDEX_sp)
        {
            std::get<sp>(this->space[ii]) *= other;
        }
        else if (ind == CompositeAlgebra::INDEX_su)
        {
            std::get<su>(this->space[ii]) *= other;
        }
    }

    return *this;
}

CompositeAlgebra CompositeAlgebra::operator/(const double other) const
{
    CompositeAlgebra out;

    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == INDEX_cn)
        {
            out.space.push_back(std::get<cn>(this->space[ii]) / other);
        }
        else if (ind == INDEX_glr)
        {
            out.space.push_back(std::get<glr>(this->space[ii]) / other);
        }
        else if (ind == INDEX_glc)
        {
            out.space.push_back(std::get<glc>(this->space[ii]) / other);
        }
        else if (ind == INDEX_rn)
        {
            out.space.push_back(std::get<rn>(this->space[ii]) / other);
        }
        else if (ind == INDEX_se)
        {
            out.space.push_back(std::get<se>(this->space[ii]) / other);
        }
        else if (ind == INDEX_so)
        {
            out.space.push_back(std::get<so>(this->space[ii]) / other);
        }
        else if (ind == INDEX_sp)
        {
            out.space.push_back(std::get<sp>(this->space[ii]) / other);
        }
        else if (ind == INDEX_su)
        {
            out.space.push_back(std::get<su>(this->space[ii]) / other);
        }
    }

    return out;
}

CompositeAlgebra & CompositeAlgebra::operator/=(const double other)
{
    for (size_t ii = 0; ii < this->space.size(); ii++)
    {
        const size_t ind = this->space[ii].index();
        if (ind == INDEX_cn)
        {
            std::get<cn>(this->space[ii]) /= other;
        }
        else if (ind == INDEX_glr)
        {
            std::get<glr>(this->space[ii]) /= other;
        }
        else if (ind == INDEX_glc)
        {
            std::get<glc>(this->space[ii]) /= other;
        }
        else if (ind == INDEX_rn)
        {
            std::get<rn>(this->space[ii]) /= other;
        }
        else if (ind == INDEX_se)
        {
            std::get<se>(this->space[ii]) /= other;
        }
        else if (ind == INDEX_so)
        {
            std::get<so>(this->space[ii]) /= other;
        }
        else if (ind == INDEX_sp)
        {
            std::get<sp>(this->space[ii]) /= other;
        }
        else if (ind == INDEX_su)
        {
            std::get<su>(this->space[ii]) /= other;
        }
    }

    return *this;
}

CompositeAlgebra::TYPES CompositeAlgebra::operator[](const ptrdiff_t index) const
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

    const size_t ind = this->space[_index].index();
    if (ind == CompositeAlgebra::INDEX_cn)
    {
        return std::get<cn>(this->space[_index]);
    }
    else if (ind == CompositeAlgebra::INDEX_glr)
    {
        return std::get<glr>(this->space[_index]);
    }
    else if (ind == CompositeAlgebra::INDEX_glc)
    {
        return std::get<glc>(this->space[_index]);
    }
    else if (ind == CompositeAlgebra::INDEX_rn)
    {
        return std::get<rn>(this->space[_index]);
    }
    else if (ind == CompositeAlgebra::INDEX_se)
    {
        return std::get<se>(this->space[_index]);
    }
    else if (ind == CompositeAlgebra::INDEX_so)
    {
        return std::get<so>(this->space[_index]);
    }
    else if (ind == CompositeAlgebra::INDEX_sp)
    {
        return std::get<sp>(this->space[_index]);
    }
    else if (ind == CompositeAlgebra::INDEX_su)
    {
        return std::get<su>(this->space[_index]);
    }

    // This should never be called
    return Lielab::domain::glc();
}

std::ostream & operator<<(std::ostream & os, const CompositeAlgebra & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */

    os << other.to_string();
    return os;
}

}
