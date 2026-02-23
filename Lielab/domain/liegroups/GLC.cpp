#include "GLC.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>


namespace Lielab::domain
{

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

// GLC::~GLC();

GLC::GLC(const GLC::matrix_t& matrix)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& GLC \\ (\mathbb{R}^{n \times 1}) &\rightarrow& GLC \f}
    *
    * Constructor instantiating an \f$GLC\f$ object from either an
    * \f$n \times n\f$ imaginary matrix or \f$n \times 1\f$ imaginary vector.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->point.noalias() = matrix;
}

GLC GLC::identity(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{GLC} \f}
    *
    * Returns an identity group with given shape.
    * 
    * @param[in] shape The shape of the group.
    * @param[out] out The GLC element.
    */

    return GLC(shape);
}

GLC GLC::project(const GLC::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{C}^{n \times n}) \rightarrow \mathbb{C}^{n \times n} \in GLC \f}
    *
    */

    const size_t shape = std::min(matrix.rows(), matrix.cols());
    return GLC(matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape)));
}

GLC::GLC(const int shape)
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
    
    this->point.noalias() = point_t::Identity(shape, shape);
}

std::string GLC::to_string() const
{
    return "GL(" + std::to_string(this->get_shape()) + ", C)";
}

int GLC::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    return 2*static_cast<int>(std::pow(this->get_shape(), 2));
}

int GLC::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
        *
        * Gets the size of the data representation.
        */

    return 2*static_cast<int>(std::pow(this->get_shape(), 2));
}

bool GLC::is_abelian() const
{
    return false;
}

int GLC::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return static_cast<int>(this->point.rows());
}

GLC::point_t GLC::get_point() const
{
    return this->point;
}

Eigen::VectorXd GLC::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */
    
    const Eigen::MatrixXcd A = this->get_matrix();

    Eigen::VectorXd out = Eigen::VectorXd::Zero(this->get_dimension());
    int kk = 0;
    for (int ii = 0; ii < this->get_shape(); ii++)
    {
        for (int jj = 0; jj < this->get_shape(); jj++)
        {
            out(kk) = std::real(A(ii,jj));
            out(kk+1) = std::imag(A(ii,jj));
            kk = kk + 2;
        }
    }

    return out;
}

void GLC::unserialize(const Eigen::VectorXd& serialized)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the GL object from a serialized vector.
    */
    
    const int vdim = static_cast<int>(serialized.size());
    const int max_ind = std::min(this->get_dimension(), vdim);
    
    for (int vind = 0; vind < max_ind; vind++)
    {
        const int rem = vind % 2;
        const int row = (vind / 2) / this->get_shape();
        const int col = (vind / 2) % this->get_shape();
        
        if (rem == 0)
        {
            this->point(row, col).real(serialized(vind));
        }
        else
        {
            this->point(row, col).imag(serialized(vind));
        }
    }
}

void GLC::unserialize(std::initializer_list<double> serialized)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the GL object from a serialized vector.
    */
    
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
}

GLC::matrix_t GLC::get_matrix() const
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

    return this->point;
}

GLC::field_t GLC::operator()(const int index1, const int index2) const
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
    
    return this->point(_index1, _index2);
}

GLC GLC::operator*(const GLC& other) const
{
    /*! \f{equation*}{ (GLC, GLC) \rightarrow GLC \f}
    *
    * Group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return GLC(this->point * other.point);
}

GLC& GLC::operator*=(const GLC& other)
{
    /*! \f{equation*}{ (GLC, GLC) \rightarrow GLC \f}
    *
    * In place group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point *= other.point;
    return *this;
}

GLC GLC::inverse() const
{
    /*! \f{equation*}{ (GLC) \rightarrow GLC \f}
    * 
    * Returns the inverse.
    */

    return GLC(this->point.inverse());
}

}
