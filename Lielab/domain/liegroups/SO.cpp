#include "SO.hpp"

#include "SU.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab::domain
{

std::string SO::to_string() const
{
    return "SO(" + std::to_string(this->_shape) + ")";
}

SO::SO() : SO(0)
{
    /*! \f{equation*}{ () \rightarrow SO \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::SO x, y, z;
    * 
    */

}

SO::SO(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SO \f}
    *
    * Constructor instantiating an \f$SO\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::SO x(2), y(3), z(4);
    * 
    * @param[in] shape The shape of the data matrix.
    */

    this->data = data_t::Identity(shape, shape);
    this->_shape = shape;
}

SO SO::from_shape(const size_t shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{SO} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The SO element. 
    */

    return SO(shape);
}

Eigen::MatrixXd SO::project(const Eigen::MatrixXd & other)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{SO} \f}
    *
    * Projects a matrix suitable for data.
    */

    const size_t shape = std::min(other.rows(), other.cols());
    const Eigen::MatrixXd other_square = other(Eigen::seqN(0, shape), Eigen::seqN(0, shape));

    Eigen::HouseholderQR<Eigen::MatrixXd> qr(other_square);
    Eigen::MatrixXd Q = qr.householderQ();
    Eigen::MatrixXd R = qr.matrixQR();
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(shape, shape);

    for (size_t ii = 0; ii < shape; ii++)
    {
        if (R(ii, ii) < 0)
        {
            P(ii, ii) = -1.0;
        }
        else
        {
            P(ii, ii) = 1.0;
        }
    }

    return Q*P;
}

size_t SO::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    return this->_shape * (this->_shape - 1) / 2;
}

size_t SO::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return this->_shape;
}

size_t SO::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
    *
    * Gets the size of the data representation.
    */

    return static_cast<size_t>(std::pow(this->_shape, 2));
}

Eigen::MatrixXd SO::get_matrix() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times n} \f}
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

SO SO::inverse() const
{
    /*! \f{equation*}{ (SO) \rightarrow SO \f}
    * 
    * Returns the inverse.
    */

    return SO(this->data.transpose());
}

Eigen::VectorXd SO::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */

    return this->data.reshaped<Eigen::RowMajor>();
}

void SO::unserialize(const Eigen::VectorXd &vec)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the GL object from a serialized vector.
    */

    const size_t vdim = vec.size();
    const size_t max_ind = std::min(this->get_size(), vdim);

    for (size_t vind = 0; vind < max_ind; vind++)
    {
        const size_t row = static_cast<size_t>(std::floor(vind / this->_shape));
        const size_t col = vind % this->_shape;
        this->data(row, col) = vec(vind);
    }
}

void SO::unserialize(std::initializer_list<double> vec)
{
    this->unserialize(Eigen::VectorXd{std::move(vec)});
}

double SO::operator()(const ptrdiff_t index1, const ptrdiff_t index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t shape = this->get_shape();
    if (shape == 0) return nan;

    if (index1 >= static_cast<ptrdiff_t>(shape)) return nan;
    if (index2 >= static_cast<ptrdiff_t>(shape)) return nan;
    if (std::abs(index1) > static_cast<ptrdiff_t>(shape)) return nan;
    if (std::abs(index2) > static_cast<ptrdiff_t>(shape)) return nan;

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

SO SO::operator*(const SO & other) const
{
    /*! \f{equation*}{ (SO, SO) \rightarrow SO \f}
    *
    * Group product.
    */

    assert(this->_shape == other.get_shape());
    return this->data * other.data;
}

SO & SO::operator*=(const SO & other)
{
    /*! \f{equation*}{ (SO, SO) \rightarrow SO \f}
    *
    * In place group product.
    */

    assert(this->_shape == other.get_shape());
    this->data *= other.data;
    return *this;
}

std::ostream & operator<<(std::ostream& os, const SO & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */

    os << static_cast<const Eigen::MatrixXd>(other.data);
    return os;
}

template SO SO::from_eulerangles_body123<double>(const double, const double, const double);
template SO SO::from_eulerangles_body231<double>(const double, const double, const double);
template SO SO::from_eulerangles_body312<double>(const double, const double, const double);
template SO SO::from_eulerangles_body132<double>(const double, const double, const double);
template SO SO::from_eulerangles_body213<double>(const double, const double, const double);
template SO SO::from_eulerangles_body321<double>(const double, const double, const double);
template SO SO::from_eulerangles_body121<double>(const double, const double, const double);
template SO SO::from_eulerangles_body131<double>(const double, const double, const double);
template SO SO::from_eulerangles_body212<double>(const double, const double, const double);
template SO SO::from_eulerangles_body232<double>(const double, const double, const double);
template SO SO::from_eulerangles_body313<double>(const double, const double, const double);
template SO SO::from_eulerangles_body323<double>(const double, const double, const double);
template SO SO::from_eulerangles_space123<double>(const double, const double, const double);
template SO SO::from_eulerangles_space231<double>(const double, const double, const double);
template SO SO::from_eulerangles_space312<double>(const double, const double, const double);
template SO SO::from_eulerangles_space132<double>(const double, const double, const double);
template SO SO::from_eulerangles_space213<double>(const double, const double, const double);
template SO SO::from_eulerangles_space321<double>(const double, const double, const double);
template SO SO::from_eulerangles_space121<double>(const double, const double, const double);
template SO SO::from_eulerangles_space131<double>(const double, const double, const double);
template SO SO::from_eulerangles_space212<double>(const double, const double, const double);
template SO SO::from_eulerangles_space232<double>(const double, const double, const double);
template SO SO::from_eulerangles_space313<double>(const double, const double, const double);
template SO SO::from_eulerangles_space323<double>(const double, const double, const double);
template SO SO::from_quaternion<double>(const double, const double, const double, const double);
template SO SO::from_rodriguesvector<double>(const double, const double, const double);

SO SO::from_SU2(const SU & other)
{
    /*!
        * Transforms an SU(2) object to a Direction Cosine object.
        *
        * @param[in] q A quaternion as an SU object of shape 2.
        * @param[out] dcm A direction cosine as an SO object of shape 3.
        */

    if (other.get_shape() != 2)
    {
        throw std::domain_error("SO::from_SU2: Expected input shape 2. Got " + std::to_string(other.get_shape()) + ".");
    }

    const auto [q0, q1, q2, q3] = other.to_quaternion();

    SO dcm(3);
    dcm.data(0,0) = 1.0 - 2.0*q2*q2 - 2.0*q3*q3;
    dcm.data(0,1) = 2.0*q1*q2 - 2.0*q3*q0;
    dcm.data(0,2) = 2.0*q1*q3 + 2.0*q2*q0;
    dcm.data(1,0) = 2.0*q1*q2 + 2.0*q3*q0;
    dcm.data(1,1) = 1.0 - 2.0*q1*q1 - 2.0*q3*q3;
    dcm.data(1,2) = 2.0*q2*q3 - 2.0*q1*q0;
    dcm.data(2,0) = 2.0*q1*q3 - 2.0*q2*q0;
    dcm.data(2,1) = 2.0*q2*q3 + 2.0*q1*q0;
    dcm.data(2,2) = 1.0 - 2.0*q1*q1 - 2.0*q2*q2;

    return dcm;
}

template std::array<double, 4> SO::to_quaternion<double>() const;
template std::array<double, 3> SO::to_gibbs<double>() const;

template std::array<double, 3> SO::to_eulerangles_body123<double>() const;
template std::array<double, 3> SO::to_eulerangles_body231<double>() const;
template std::array<double, 3> SO::to_eulerangles_body312<double>() const;
template std::array<double, 3> SO::to_eulerangles_body132<double>() const;
template std::array<double, 3> SO::to_eulerangles_body213<double>() const;
template std::array<double, 3> SO::to_eulerangles_body321<double>() const;
template std::array<double, 3> SO::to_eulerangles_body121<double>() const;
template std::array<double, 3> SO::to_eulerangles_body131<double>() const;
template std::array<double, 3> SO::to_eulerangles_body212<double>() const;
template std::array<double, 3> SO::to_eulerangles_body232<double>() const;
template std::array<double, 3> SO::to_eulerangles_body313<double>() const;
template std::array<double, 3> SO::to_eulerangles_body323<double>() const;
template std::array<double, 3> SO::to_eulerangles_space123<double>() const;
template std::array<double, 3> SO::to_eulerangles_space231<double>() const;
template std::array<double, 3> SO::to_eulerangles_space312<double>() const;
template std::array<double, 3> SO::to_eulerangles_space132<double>() const;
template std::array<double, 3> SO::to_eulerangles_space213<double>() const;
template std::array<double, 3> SO::to_eulerangles_space321<double>() const;
template std::array<double, 3> SO::to_eulerangles_space121<double>() const;
template std::array<double, 3> SO::to_eulerangles_space131<double>() const;
template std::array<double, 3> SO::to_eulerangles_space212<double>() const;
template std::array<double, 3> SO::to_eulerangles_space232<double>() const;
template std::array<double, 3> SO::to_eulerangles_space313<double>() const;
template std::array<double, 3> SO::to_eulerangles_space323<double>() const;

}
