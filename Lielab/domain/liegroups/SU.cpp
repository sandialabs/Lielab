#include "SU.hpp"

#include "SO.hpp"

#include "Lielab/testing.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

SU::SU() : SU(0)
{
    /*! \f{equation*}{ () \rightarrow SU \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::SU x, y, z;
    * 
    */

}

// SU::~SU();

SU::SU(const SU::matrix_t& matrix)
{
    /*! \f{equation}{(\mathbb{C}^{n \times n}) \rightarrow SU \f}
    *
    * Constructor instantiating an \f$SU\f$ object from an
    * \f$n \times n\f$ imaginary matrix.
    *
    */

    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->point.noalias() = matrix;
}

SU SU::identity(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SU \f}
    *
    * Returns an identity group with given shape.
    * 
    * @param[in] shape The shape of the group.
    * @param[out] out The SU element. 
    */

    return SU(shape);
}

// TODO: project

SU::SU(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SU \f}
    *
    * Constructor instantiating an \f$SU\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::SU x(2), y(4), z(6);
    * 
    * @param[in] shape The shape of the data matrix.
    */

    this->point.noalias() = Eigen::MatrixXcd::Identity(shape, shape);
}

SU SU::from_quaternion(const double e0, const double e1, const double e2, const double e3)
{
    /*! \f{equation*}{ (\mathbb{R}^4) \rightarrow SU \f}
     *
     * Constructor instantiating a Quaternion as an \f$SU\f$ object.
     * 
     * Enables instatiation like:
     * 
     *     Lielab::domain::SU Quaternion0 = Lielab::domain::SU::from_quaternion(1.0, 0.0, 0.0, 0.0);
     * 
     * @param[out] quaternion An SU object representing the Quaternion.
     */

    constexpr std::complex<double> j(0.0, 1.0);
    SU qout = SU(2);

    qout.point(0,0) = e0 + e1*j;
    qout.point(1,1) = e0 - e1*j;
    qout.point(0,1) = -e2 + e3*j;
    qout.point(1,0) = e2 + e3*j;

    return qout;
}

SU SU::from_SO3(const SO& dcm)
{
    /*!
    * Transforms an SO(3) object into an SU(2) object.
    *
    * @param[in] dcm A direction cosine as an SO object of shape 3.
    * @param[out] q A quaternion as an SU object of shape 2.
    */

    lielab_assert(dcm.get_shape(), "Expected input shape 3. Got " + std::to_string(dcm.get_shape()) + ".");

    const auto [q0, q1, q2, q3] = dcm.to_quaternion();

    return SU::from_quaternion(q0, q1, q2, q3);

}

std::string SU::to_string() const
{
    return "SU(" + std::to_string(this->get_shape()) + ")";
}

int SU::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    if (this->get_shape() == 0) return 0;
    return this->get_shape() * this->get_shape() - 1;
}

int SU::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
        *
        * Gets the size of the data representation.
        */

    return 2*static_cast<int>(std::pow(this->get_shape(), 2));
}

bool SU::is_abelian() const
{
    return false;
}

int SU::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return static_cast<int>(this->point.rows());
}

SU::point_t SU::get_point() const
{
    return this->point;
}

Eigen::VectorXd SU::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */
    
    if (this->get_shape() == 0) return Eigen::VectorXd::Zero(0);

    const Eigen::MatrixXcd A = this->get_matrix();

    Eigen::VectorXd out = Eigen::VectorXd::Zero(this->get_size());
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

void SU::unserialize(const Eigen::VectorXd& serialized)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow SU \f}
    * 
    * Sets the SP object from a serialized vector.
    */

    const int vdim = static_cast<int>(serialized.size());
    const int max_ind = std::min(this->get_size(), vdim);
    
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

void SU::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
}

SU::matrix_t SU::get_matrix() const
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

SU::field_t SU::operator()(const int index1, const int index2) const
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

std::array<double, 4> SU::to_quaternion() const
{
    /*! \f{equation*}{ (SU) \rightarrow \mathbb{R}^4 \f}
     *
     * Method returning a quaternion as an \f$\mathbb{R}^4\f$ object.
     * 
     * Enables outputs like:
     * 
     *     std::array<double, 4> Quaternion0 = Qobj.to_quaternion();
     * 
     * @param[out] quaternion An array representing the Quaternion.
     */

    if (this->get_shape() != 2)
    {
        throw std::domain_error("SU::to_quaternion: Expected input shape 2. Got " + std::to_string(this->get_shape()) + ".");
    }

    const double q0 = (this->point(0,0).real() + this->point(1,1).real())/2.0;
    const double q1 = (this->point(0,0).imag() - this->point(1,1).imag())/2.0;
    const double q2 = (-this->point(0,1).real() + this->point(1,0).real())/2.0;
    const double q3 = (this->point(0,1).imag() + this->point(1,0).imag())/2.0;

    return {q0, q1, q2, q3};
}

SU SU::operator*(const SU& other) const
{
    /*! \f{equation*}{ (SU, SU) \rightarrow SU \f}
    *
    * Group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    return SU(this->point * other.point);
}

SU& SU::operator*=(const SU& other)
{
    /*! \f{equation*}{ (SU, SU) \rightarrow SU \f}
    *
    * In place group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    this->point *= other.point;
    return *this;
}

SU SU::inverse() const
{
    /*! \f{equation*}{ () \rightarrow SU \f}
    * 
    * Returns the inverse.
    */

    return SU(this->point.inverse());
}

}
