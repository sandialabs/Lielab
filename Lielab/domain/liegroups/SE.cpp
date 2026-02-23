#include "SE.hpp"

#include "RN.hpp"
#include "SO.hpp"

#include "Lielab/testing.hpp"
#include "Lielab/utils/eigentools.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <tuple>

namespace Lielab::domain
{

SE::SE() : SE(0)
{
    /*! \f{equation*}{ () \rightarrow SE \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::SE x, y, z;
    * 
    */

}

// SE::~SE();

SE::SE(const SE::matrix_t& matrix)
{
    /*! \f{equation}{ (\mathbb{R}^{n \times n}) \rightarrow SE \f}
    *
    * Constructor instantiating an \f$SE\f$ object from an
    * \f$n \times n\f$ real matrix.
    *
    */

    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->_shape = static_cast<int>(matrix.rows());

    if (this->_shape == 0)
    {
        this->point = std::make_tuple(RN::identity(0), SO::identity(0));
    }
    else
    {
        const RN d1 = RN::from_vector(matrix(Eigen::seqN(0, this->_shape - 1), this->_shape - 1));
        const SO d2 = SO(matrix(Eigen::seqN(0, this->_shape - 1), Eigen::seqN(0, this->_shape - 1)));

        this->point = std::make_tuple(d1, d2);
    }
}

SE SE::identity(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{SE} \f}
    *
    * Returns an identity group with given shape.
    * 
    * @param[in] shape The shape of the group.
    * @param[out] out The SE element. 
    */

    if (shape == 0)
    {
        SE out;
        out._shape = 0;
        out.point = std::make_tuple(RN::identity(0), SO::identity(0));
        return out;
    }

    return SE(shape - 1);
}

SE SE::project(const SE::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{SE} \f}
    *
    */

    const size_t shape = std::min(matrix.rows(), matrix.cols());

    if (shape == 0) return SE::identity(0);

    const SE::matrix_t matrix_sq = matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape));
    const Eigen::MatrixXd SO_component = matrix_sq(Eigen::seqN(0, shape - 1), Eigen::seqN(0, shape - 1));
    const Eigen::MatrixXd RN_component = matrix_sq(Eigen::seqN(0, shape - 1), shape - 1);

    return SE(RN::from_vector(RN_component), SO::project(SO_component));
}

SE::SE(const int n)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow SE \f}
    *
    * Constructor instantiating an \f$SE\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::SE x(2), y(3), z(4);
    * 
    * @param[in] shape The shape of the data matrix.
    */

    this->_shape = n + 1;
    this->point = std::make_tuple(RN(n), SO(n));
}

SE::SE(const RN& RN_c, const SO& SO_c)
{
    const int shape = std::min(RN_c.get_shape(), SO_c.get_shape() + 1);
    this->_shape = shape;

    if (this->_shape == 0)
    {
        this->point = std::make_tuple(RN::identity(0), SO::identity(0));
        return;
    }

    const Eigen::VectorXd RN_cbar = RN_c.serialize()(Eigen::seqN(0, shape - 1));
    const Eigen::MatrixXd SO_chat = SO_c.get_matrix()(Eigen::seqN(0, shape - 1), Eigen::seqN(0, shape - 1));

    this->point = std::make_tuple(RN::from_vector(RN_cbar), SO(SO_chat));
}

std::string SE::to_string() const
{
    if (this->_shape == 0) return "SE(0)";
    return "SE(" + std::to_string(this->_shape-1) + ")";
}

int SE::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the group.
    */

    return (this->_shape - 1) * (this->_shape - 2) / 2 + this->_shape - 1;
}

int SE::get_size() const
{
    /*! \f{quation*}{ () \rightarrow \mathbb{Z} \f}
        *
        * Gets the size of the data representation.
        */

    return std::get<RN>(this->point).get_size() + std::get<SO>(this->point).get_size();
}

bool SE::is_abelian() const
{
    return false;
}

int SE::get_shape() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the shape of the group.
    */

    return this->_shape;
}

SE::point_t SE::get_point() const
{
    return this->point;
}

Eigen::VectorXd SE::serialize() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns a serialized representation.
    */

    return Lielab::utils::concatenate({std::get<RN>(this->point).serialize(), std::get<SO>(this->point).serialize()});
}

void SE::unserialize(const Eigen::VectorXd& serialized)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times 1}) \rightarrow () \f}
    * 
    * Sets the SE object from a serialized vector.
    */

    const size_t vdim = serialized.size();
    const size_t rnsize = std::get<RN>(this->point).get_size();
    
    std::get<RN>(this->point).unserialize(serialized);

    if (vdim <= rnsize) return;

    std::get<SO>(this->point).unserialize(serialized(Eigen::seqN(rnsize, vdim - rnsize)));
}

void SE::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
}

SE::matrix_t SE::get_matrix() const
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

    const size_t shape = this->get_shape();
    SE::matrix_t out = SE::matrix_t::Identity(shape, shape);
    if (shape == 0 || shape == 1) return out;
    
    out(Eigen::seqN(0, shape - 1), shape - 1) = std::get<RN>(this->point).serialize();
    out(Eigen::seqN(0, shape - 1), Eigen::seqN(0, shape - 1)) = std::get<SO>(this->point).get_matrix();

    return out;
}

SE::field_t SE::operator()(const int index1, const int index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int shape = this->get_shape();
    if (shape == 0) return nan;

    // If input index is negative, index from the back of the array
    const int _index1 = (index1 < 0) ? shape + index1 : index1;
    const int _index2 = (index2 < 0) ? shape + index2 : index2;

    // Error check for out of bounds
    if (_index1 < 0) return nan;
    if (_index1 >= shape) return nan;
    if (_index2 < 0) return nan;
    if (_index2 >= shape) return nan;

    if (_index1 == shape - 1 && _index2 == shape - 1)
    {
        return 1.0;
    }

    if (_index1 == shape - 1)
    {
        return 0.0;
    }

    if (_index2 < shape - 1)
    {
        return std::get<SO>(this->point)(_index1, _index2);
    }

    return std::get<RN>(this->point)(_index1, _index2);
}

SE SE::operator*(const SE& other) const
{
    /*! \f{equation*}{ (SE, SE) \rightarrow SE \f}
    *
    * Group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    SE out = SE::identity(this->get_shape());
    std::get<RN>(out.point).point.noalias() = std::get<SO>(this->point).get_matrix()*std::get<RN>(other.point).serialize() + std::get<RN>(this->point).serialize();
    std::get<SO>(out.point) = std::get<SO>(this->point)*std::get<SO>(other.point);

    return out;
}

SE& SE::operator*=(const SE& other)
{
    /*! \f{equation*}{ (SE, SE) \rightarrow SE \f}
    *
    * In place group product.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    std::get<RN>(this->point).point.noalias() = std::get<SO>(this->point).get_matrix()*std::get<RN>(other.point).serialize() + std::get<RN>(this->point).serialize();
    std::get<SO>(this->point) *= std::get<SO>(other.point);

    return *this;
}

SE SE::inverse() const
{
    /*! \f{equation*}{ (SE) \rightarrow SE \f}
    * 
    * Returns the inverse.
    */

    SE out = SE::identity(this->get_shape());
    std::get<SO>(out.point) = std::get<SO>(this->point).inverse();
    std::get<RN>(out.point) = RN::from_vector(-std::get<SO>(out.point).get_matrix()*std::get<RN>(this->point).serialize());

    return out;
}

}
