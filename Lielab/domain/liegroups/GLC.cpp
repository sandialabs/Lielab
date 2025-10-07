#include "GLC.hpp"

#include "Lielab/utils/Error.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string GLC::to_string() const
{
    return "GL(" + std::to_string(this->_shape) + ", C)";
}

GLC::GLC() : GLC(0)
{
    /*! \f{equation*}{ () \rightarrow GLC \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::GLC x, y, z;
    * 
    */

}

GLC::GLC(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow GLC \f}
    *
    * Constructor instantiating an \f$GLC\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::GLC x(2), y(3), z(4);
    * 
    * @param[in] shape The shape of the data matrix.
    */
    
    this->data = data_t::Identity(shape, shape);
    this->_shape = shape;
}

GLC GLC::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{GLC} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The GLC element.
    */

    return GLC(shape);
}

Eigen::MatrixXcd GLC::project(const Eigen::MatrixXcd & other)
{
    /*! \f{equation*}{ (\mathbb{C}^{n \times n}) \rightarrow \mathbb{C}^{n \times n} \in GLC \f}
    *
    * Projects a matrix suitable for data.
    */

    const size_t shape = std::min(other.rows(), other.cols());
    return other(Eigen::seqN(0, shape), Eigen::seqN(0, shape));
}

size_t GLC::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    return static_cast<size_t>(2*std::pow(this->_shape, 2));
}

size_t GLC::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return this->_shape;
}

size_t GLC::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
        *
        * Gets the size of the data representation.
        */

    return static_cast<size_t>(2*std::pow(this->_shape, 2));
}

Eigen::MatrixXcd GLC::get_matrix() const
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

    return this->data;
}

GLC GLC::inverse() const
{
    /*! \f{equation*}{ (GLC) \rightarrow GLC \f}
    * 
    * Returns the inverse.
    */

    return this->data.inverse();
}

// Data representation

Eigen::VectorXd GLC::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */
    
    const Eigen::MatrixXcd A = this->get_matrix();

    Eigen::VectorXd out = Eigen::VectorXd::Zero(this->get_dimension());
    size_t kk = 0;
    for (size_t ii = 0; ii < this->_shape; ii++)
    {
        for (size_t jj = 0; jj < this->_shape; jj++)
        {
            out(kk) = std::real(A(ii,jj));
            out(kk+1) = std::imag(A(ii,jj));
            kk = kk + 2;
        }
    }

    return out;
}

void GLC::unserialize(const Eigen::VectorXd& vector)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the GL object from a serialized vector.
    */
    
    const size_t vdim = vector.size();
    const size_t max_ind = std::min(this->get_dimension(), vdim);
    
    for (size_t vind = 0; vind < max_ind; vind++)
    {
        const size_t rem = vind % 2;
        const size_t row = static_cast<size_t>(std::floor((vind/2) / this->_shape));
        const size_t col = (vind/2) % this->_shape;
        
        if (rem == 0)
        {
            this->data(row, col).real(vector(vind));
        }
        else
        {
            this->data(row, col).imag(vector(vind));
        }
    }
}

void GLC::unserialize(std::initializer_list<double> vector)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the GL object from a serialized vector.
    */
    
    this->unserialize(Eigen::VectorXd{std::move(vector)});
}

std::complex<double> GLC::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
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
    
    return this->data(_index1, _index2);
}

GLC GLC::operator*(const GLC & other) const
{
    /*! \f{equation*}{ (GLC, GLC) \rightarrow GLC \f}
    *
    * Group product.
    */

    assert(this->_shape == other.get_shape());

    return GLC(this->data * other.data);
}

GLC & GLC::operator*=(const GLC & other)
{
    /*! \f{equation*}{ (GLC, GLC) \rightarrow GLC \f}
    *
    * In place group product.
    */

    assert(this->_shape == other.get_shape());
    this->data *= other.data;
    return *this;
}

std::ostream & operator<<(std::ostream& os, const GLC & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */

    os << static_cast<const Eigen::MatrixXcd>(other.data);
    return os;
}

}
