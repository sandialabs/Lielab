#include "se.hpp"

#include "rn.hpp"
#include "so.hpp"

#include "Lielab/testing.hpp"
#include "Lielab/utils/eigentools.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <tuple>

namespace Lielab::domain
{

se::se() : se(0)
{
    /*! \f{equation*}{() \rightarrow \mathfrak{se} \f}
    * 
    * Empty initialization function. Enables instantiation like:
    * 
    *     Lielab::domain::se x, y, z;
    * 
    */

}

// se::~se();

se::se(const se::matrix_t& matrix)
{
    /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& \mathfrak{se} \f}
    *
    * Constructor instantiating an \f$\mathfrak{se}\f$ object from either an
    * \f$n \times n\f$ real matrix or \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real matrix.
    */

    lielab_assert(matrix.rows() == matrix.cols(), "Input matrix must be square.");

    this->_shape = static_cast<int>(matrix.rows());

    if (this->_shape == 0)
    {
        this->point = std::make_tuple(rn::zero(0), so::zero(0));
    }
    else
    {
        const rn d1 = rn::from_vector(matrix(Eigen::seqN(0, this->_shape - 1), this->_shape - 1));
        const so d2 = so(matrix(Eigen::seqN(0, this->_shape - 1), Eigen::seqN(0, this->_shape - 1)));

        this->point = std::make_tuple(d1, d2);
    }
}

se se::basis(const int index, const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{se} \f}
    *
    * Returns the i'th basis element of the se algebra.
    * 
    * @param[in] index The basis vector.
    * @param[in] shape The shape of the algebra.
    * @param[out] out The se element.
    */

    se out(shape);
    if (index < 0) return out;

    const int dim = out.get_dimension();
    if (index >= dim) return out;

    Eigen::VectorXd v = Eigen::VectorXd::Zero(dim);
    v(index) = 1.0;
    out.set_vector(v);
    return out;
}

se se::zero(const int shape)
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{se} \f}
    *
    * Returns a zero algebra with given shape.
    * 
    * @param[in] shape The shape of the algebra.
    * @param[out] out The se element. 
    */

    if (shape == 0)
    {
        se out;
        out.point = std::make_tuple(rn::zero(0), so::zero(0));
        out._shape = shape;
        return out;
    }

    return se(shape - 1);
}

se se::from_vector(const Eigen::VectorXd& other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{se} \f}
    *
    * Constructor instantiating an \f$\mathfrak{se}\f$ object from either a
    * \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real vector.
    */

    const int len = static_cast<int>(other.size());
    const int shape = static_cast<int>(std::ceil(std::sqrt(8.0*len + 1.0)/2.0 + 0.5));

    se out = se::zero(shape);
    out.set_vector(other);
    return out;
}

se se::from_vector(std::initializer_list<double> other)
{
    /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{se} \f}
    *
    * Constructor instantiating an \f$\mathfrak{se}\f$ object from either a
    * \f$n \times 1\f$ real vector.
    *
    * @param[in] other The object to instantiate from as a real vector.
    */

    return se::from_vector(Eigen::VectorXd{std::move(other)});
}

se se::project(const se::matrix_t& matrix)
{
    /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{se} \f}
    *
    * Projects a matrix suitable for data.
    */

    const size_t shape = std::min(matrix.rows(), matrix.cols());

    if (shape == 0) return se::zero(0);

    const se::matrix_t matrix_sq = matrix(Eigen::seqN(0, shape), Eigen::seqN(0, shape));
    const Eigen::MatrixXd so_component = matrix_sq(Eigen::seqN(0, shape - 1), Eigen::seqN(0, shape - 1));
    const Eigen::MatrixXd rn_component = matrix_sq(Eigen::seqN(0, shape - 1), shape - 1);

    return se(rn::from_vector(rn_component), so::project(so_component));
}

se::se(const int n)
{
    /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{se} \f}
    *
    * Constructor instantiating an \f$\mathfrak{se}\f$ object.
    * 
    * Enables instantiation like:
    * 
    *     Lielab::domain::se x(3), y(4), z(5);
    * 
    * @param[in] shape The shape of the data matrix.
    */
    
    this->point = std::make_tuple(rn(n), so(n));
    this->_shape = n + 1;
}

se::se(const rn& rn_c, const so& so_c)
{
    this->_shape = std::min(rn_c.get_shape(), so_c.get_shape() + 1);

    if (this->_shape == 0)
    {
        this->point = std::make_tuple(rn::zero(0), so::zero(0));
        return;
    }

    const Eigen::VectorXd rn_cbar = rn_c.get_vector();
    const Eigen::MatrixXd so_chat = so_c.get_matrix();

    this->point = std::make_tuple(rn::from_vector(rn_cbar), so(so_chat));
}

std::string se::to_string() const
{
    if (this->get_shape() == 0) return "se(nan)";
    return "se(" + std::to_string(this->get_shape()-1) + ")";
}

int se::get_dimension() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
    * 
    * Gets the dimension of the algebra.
    */

    if (this->get_shape() == 0) return 0; // Return nan?
    if (this->get_shape() == 1) return 0; // Return nan?

    return (this->get_shape() - 1) * (this->get_shape() - 2) / 2 + this->get_shape() - 1;
}

int se::get_size() const
{
    return this->get_dimension();
}

bool se::is_abelian() const
{
    return false;
}

int se::get_shape() const
{
    return this->_shape;
}

se::point_t se::get_point() const
{
    /*!
    */

    return this->point;
}

Eigen::VectorXd se::serialize() const
{
    return this->get_vector();
}

void se::unserialize(const Eigen::VectorXd& serialized)
{
    this->set_vector(serialized);
}

void se::unserialize(std::initializer_list<double> serialized)
{
    this->unserialize(Eigen::VectorXd{std::move(serialized)});
}

se::matrix_t se::get_matrix() const
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
    
    const int shape = this->get_shape();
    se::matrix_t out = se::matrix_t::Zero(shape, shape);
    if (shape == 0 || shape == 1) return out;

    out(Eigen::seqN(0, shape - 1), shape - 1) = std::get<rn>(this->point).get_vector();
    out(Eigen::seqN(0, shape - 1), Eigen::seqN(0, shape - 1)) = std::get<so>(this->point).get_matrix();

    return out;
}

Eigen::VectorXd se::get_vector() const
{
    /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
    * 
    * Returns the vector representation.
    */

    const size_t dim = this->get_dimension();
    Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);
    if (this->get_shape() == 0) return out;

    return Lielab::utils::concatenate({std::get<rn>(this->point).get_vector(), std::get<so>(this->point).get_vector()});
}

void se::set_vector(const Eigen::VectorXd& vector)
{
    /*! \f{equation*}{ \mathfrak{se} := \mathbb{R}^{n \times 1} \f}
    * 
    * @param[in] vector An Eigen::VectorXd to assign.
    */

    const size_t vdim = vector.size();
    const size_t rndim = std::get<rn>(this->point).get_dimension();

    std::get<rn>(this->point).set_vector(vector);

    if (vdim <= rndim) return;

    std::get<so>(this->point).set_vector(vector(Eigen::seqN(rndim, vector.size() - rndim)));
}

void se::set_vector(std::initializer_list<double> vector)
{
    /*!
    *
    * @param[in] vector
    */
   
    this->set_vector(Eigen::VectorXd{std::move(vector)});
}

double se::operator()(const int index) const
{
    /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the column vector representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int dim = this->get_dimension();

    // If input index is negative, index from the back of the array
    const int _index = (index < 0) ? dim + index : index;

    // Error check for out of bounds
    if (_index < 0) return nan;
    if (_index >= dim) return nan;

    const Eigen::VectorXd vector = this->get_vector();
    return vector(_index);
}

se::field_t se::operator()(const int index1, const int index2) const
{
    /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
    *
    * Gets a value in the square matrix representation.
    */

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int shape = this->get_shape();
    if (shape <= 1) return nan;

    // If input index is negative, index from the back of the array
    const int _index1 = (index1 < 0) ? shape + index1 : index1;
    const int _index2 = (index2 < 0) ? shape + index2 : index2;

    // Error check for out of bounds
    if (_index1 < 0) return nan;
    if (_index1 >= shape) return nan;
    if (_index2 < 0) return nan;
    if (_index2 >= shape) return nan;

    if (_index1 == shape - 1)
    {
        return 0.0;
    }

    if (_index2 < shape - 1)
    {
        return std::get<so>(this->point)(_index1, _index2);
    }

    return std::get<rn>(this->point)(_index1);
}

se se::operator+(const se& other) const
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * Addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    se out = se::zero(this->get_shape());
    std::get<rn>(out.point).point.noalias() = std::get<rn>(this->point).point + std::get<rn>(other.point).point;
    std::get<so>(out.point).point.noalias() = std::get<so>(this->point).point + std::get<so>(other.point).point;
    return out;
}

se& se::operator+=(const se& other)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * In place addition of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    std::get<rn>(this->point) += std::get<rn>(other.point);
    std::get<so>(this->point) += std::get<so>(other.point);
    return *this;
}

se se::operator-(const se& other) const
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * Subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    se out = se::zero(this->get_shape());
    std::get<rn>(out.point).point.noalias() = std::get<rn>(this->point).point - std::get<rn>(other.point).point;
    std::get<so>(out.point).point.noalias() = std::get<so>(this->point).point - std::get<so>(other.point).point;
    return out;
}

se& se::operator-=(const se& other)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * In place subtraction of two vectors in the algebra.
    */

    lielab_assert(this->get_shape() == other.get_shape(), "Shapes must be equal.");
    std::get<rn>(this->point) -= std::get<rn>(other.point);
    std::get<so>(this->point) -= std::get<so>(other.point);
    return *this;
}

se se::operator-() const
{
    /*! \f{equation*}{ (\mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * Unary negative of the vector.
    */

    se out = se::zero(this->get_shape());
    std::get<rn>(out.point).point.noalias() = -std::get<rn>(this->point).point;
    std::get<so>(out.point).point.noalias() = -std::get<so>(this->point).point;
    return out;
}

se se::operator*(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
    *
    * Scalar product.
    */

    se out = se::zero(this->get_shape());
    std::get<rn>(out.point) = std::get<rn>(this->point)*other;
    std::get<so>(out.point) = std::get<so>(this->point)*other;
    return out;
}

se operator*(const double other, const se& rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

se& se::operator*=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
    *
    * In place scalar multiplication.
    */

    std::get<rn>(this->point) *= other;
    std::get<so>(this->point) *= other;
    return *this;
}

se se::operator/(const double other) const
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
    *
    * Scalar division.
    */

    se out = se::zero(this->get_shape());
    std::get<rn>(out.point) = std::get<rn>(this->point) / other;
    std::get<so>(out.point) = std::get<so>(this->point) / other;
    return out;
}

se& se::operator/=(const double other)
{
    /*! \f{equation*}{ (\mathfrak{se}, \mathbb{RR}) \rightarrow \mathfrak{se} \f}
    *
    * In place scalar division.
    */

    std::get<rn>(this->point) = std::get<rn>(this->point) / other;
    std::get<so>(this->point) = std::get<so>(this->point) / other;
    return *this;
}

}

