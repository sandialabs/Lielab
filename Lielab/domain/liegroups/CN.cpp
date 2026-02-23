#include "CN.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

CN::CN() : CN(0)
{
    /*! \f{equation*}{ () \rightarrow CN \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::CN x, y, z;
    * 
    */

}

// CN:~CN();

CN::CN(const Eigen::MatrixXcd& other)
{
    /*! \f{eqnarray*}{(\mathbb{C}^{n \times n}) &\rightarrow& CN \\ (\mathbb{C}^{n \times 1}) &\rightarrow& CN \f}
    *
    * Constructor instantiating an \f$CN\f$ object from either an
    * \f$n \times n\f$ imaginary matrix or \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    lielab_assert(other.rows() == other.cols(), "Input matrix must be square.");

    this->_shape = static_cast<int>(other.rows());

    if (this->_shape == 0)
    {
        this->point.noalias() = Eigen::VectorXcd::Zero(0);
        return;
    }

    this->point = Eigen::VectorXcd::Zero(this->_shape - 1);
    for (int ii = 0; ii < this->_shape - 1; ii++)
    {
        this->point(ii) = other(ii, this->_shape - 1);
    }
}

CN CN::identity(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{CN} \f}
    *
    * Returns an identity group with given shape.
    * 
    * @param[in] shape The shape of the group.
    * @param[out] out The CN element. 
    */

    if (shape == 0)
    {
        CN out;
        out._shape = 0;
        return out;
    }

    return CN(shape-1);
}

CN CN::project(const CN::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{C}^{n \times n}) \rightarrow \mathbb{C}^{n \times n} \in CN \f}
    *
    */

    const size_t shape = std::min(matrix.rows(), matrix.cols());
    return CN(matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape)));
}

CN::CN(const int n)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow CN \f}
    *
    * Constructor instantiating an \f$CN\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::CN x(2), y(3), z(4);
    * 
    * @param[in] shape The shape of the data matrix.
    */
    
    this->_shape = n + 1;
    this->point.noalias() = Eigen::VectorXcd::Zero(n);
}

CN CN::from_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{CN} \f}
    *
    * Constructor instantiating an \f$\mathfrak{CN}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] vector The object to instantiate from as an imaginary vector.
    */

    const int shape = static_cast<int>(std::ceil(vector.size()/2.0)) + 1;
    CN out = CN::identity(shape);
    out.unserialize(vector);

    return out;
}

CN CN::from_vector(std::initializer_list<double> vector)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{CN} \f}
    *
    * Constructor instantiating an \f$\mathfrak{CN}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary vector.
    */

    return CN::from_vector(Eigen::VectorXd{std::move(vector)});
}

CN CN::from_complex_vector(const Eigen::VectorXcd& other)
{
    /*! \f{equation}{(\mathbb{C}^{n \times 1}) \rightarrow \mathfrak{CN} \f}
    *
    * Constructor instantiating an \f$\mathfrak{CN}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary vector.
    */

    const int n = static_cast<int>(other.size());
    CN out(n);
    out.point = other;
    return out;
}

CN CN::from_complex_vector(const std::initializer_list<std::complex<double>> other)
{
    /*! \f{equation}{(\mathbb{C}^{n \times 1}) \rightarrow \mathfrak{CN} \f}
    *
    * Constructor instantiating an \f$\mathfrak{CN}\f$ object from either a
    * \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as an imaginary vector.
    */

    return CN::from_complex_vector(Eigen::VectorXcd{std::move(other)});
}

std::string CN::to_string() const
{
    const size_t shape = this->get_shape();
    if (shape == 0) return "C^nan";
    return "C^" + std::to_string(shape-1);
}

int CN::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    if (this->_shape == 0) return 0; // TODO: Return nan?
    return 2*(this->_shape - 1);
}

int CN::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
        *
        * Gets the size of the data representation.
        */

    if (this->_shape == 0) return 0;
    return 2*(this->_shape - 1);
}

bool CN::is_abelian() const
{
    return true;
}

int CN::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return this->_shape;
}

CN::point_t CN::get_point() const
{
    return this->point;
}

Eigen::VectorXd CN::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */

    const size_t dim = this->get_dimension();
    Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);

    size_t kk = 0;
    for (size_t ii = 0; ii < dim/2; ii++)
    {
        out(kk) = std::real(this->point(ii));
        kk += 1;
        out(kk) = std::imag(this->point(ii));
        kk += 1;
    }

    return out;
}

void CN::unserialize(const Eigen::VectorXd& serialized)
{
    /*! \f{equation*}{ (\mathbb{C}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the CN object from a serialized vector.
    */
    
    const int vdim = static_cast<int>(serialized.size());
    const int max_ind = std::min(this->get_dimension(), vdim);

    int dind = 0;
    for (int vind = 0; vind < max_ind; vind++)
    {
        const int rem = vind % 2;

        if (rem == 0)
        {
            this->point(dind).real(serialized(vind));
        }
        else
        {
            this->point(dind).imag(serialized(vind));
            dind += 1;
        }
    }
}

void CN::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
}

CN::matrix_t CN::get_matrix() const
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

    CN::matrix_t out = CN::matrix_t::Identity(this->_shape, this->_shape);

    if (this->_shape == 0) return out;

    for (int ii = 0; ii < this->_shape - 1; ii++)
    {
        out(ii, this->_shape - 1) = this->point(ii);
    }

    return out;
}

CN::field_t CN::operator()(const int index1, const int index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int shape = this->get_shape();
    if (shape == 0) return std::complex<double>(nan, nan);

    // If input index is negative, index from the back of the array
    const int _index1 = (index1 < 0) ? shape + index1 : index1;
    const int _index2 = (index2 < 0) ? shape + index2 : index2;

    // Error check for out of bounds
    if (_index1 < 0) return std::complex<double>(nan, nan);
    if (_index1 >= shape) return std::complex<double>(nan, nan);
    if (_index2 < 0) return std::complex<double>(nan, nan);
    if (_index2 >= shape) return std::complex<double>(nan, nan);


    if (_index1 == _index2) return std::complex<double>(1.0, 0.0);
    if (_index1 == shape - 1) return std::complex<double>(0.0, 0.0);
    if (_index2 != shape - 1) return std::complex<double>(0.0, 0.0);

    return this->point(_index1);
}

Eigen::VectorXcd CN::to_complex_vector() const
{
    /*! \f{equation}{ () \rightarrow (\mathbb{C}^{n \times 1}) \f}
    *
    * @param[out] vec The object as a complex vector.
    */

    return Eigen::VectorXcd(this->point);
}

// std::complex<double>& CN::operator()(const size_t index)
// {
//     /*! \f{equation*}{ CN(\mathbb{Z}) := \mathbb{C} \f}
//     *
//     * Assignment of a value in the column vector representation.
//     */

//     return this->point(index);
// }

const CN::field_t& CN::operator[](const int index) const
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{C} \f}
    *
    */

    const int dim = this->get_dimension();

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? dim/2 + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < dim/2), "Index " + std::to_string(index) + " is out of bounds for CN of length " + std::to_string(dim/2));

    return this->point(_index);
}

CN::field_t& CN::operator[](const int index)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{C} \f}
    *
    */

    const int dim = this->get_dimension();

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? dim/2 + index : index;

    // Error check for out of bounds
    lielab_assert((_index >= 0) && (_index < dim/2), "Index " + std::to_string(index) + " is out of bounds for CN of length " + std::to_string(dim/2));

    return this->point(_index);
}

CN CN::operator*(const CN& other) const
{
    /*! \f{equation*}{ (CN, CN) \rightarrow CN \f}
    *
    * Group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return CN::from_complex_vector(this->point + other.point);
}

CN& CN::operator*=(const CN& other)
{
    /*! \f{equation*}{ (CN, CN) \rightarrow CN \f}
    *
    * In place group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point += other.point;
    return *this;
}

CN CN::inverse() const
{
    /*! \f{equation*}{ (CN) \rightarrow CN \f}
    * 
    * Returns the inverse.
    */

    return CN::from_complex_vector(-this->point);
}

}
