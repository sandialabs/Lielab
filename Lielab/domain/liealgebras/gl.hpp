#ifndef _LIELAB_DOMAIN_LIEALGEBRAS_gl_HPP
#define _LIELAB_DOMAIN_LIEALGEBRAS_gl_HPP

#include "../../abstract/abstract_all.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
namespace domain
{

using Lielab::abstract::Real;
using Lielab::abstract::Imaginary;

class gl
{
    public:

    /*! Flag specifying whether or not the algebra is abelian. */
    static constexpr bool abelian = false;
    size_t shape = 0;
    

    /*! Internal data for the algebra. */
    Eigen::MatrixXd _data;

    gl() : gl(0)
    {
        /*! \f{equation*}{() \rightarrow \mathfrak{gl} \f}
        * 
        * Empty initialization function. Enables instantiation like:
        * 
        *     Lielab::domain::gl x, y, z;
        * 
        */

    }

    gl(const size_t shape)
    {
        /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{gl} \f}
        *
        * Constructor instantiating an \f$\mathfrak{gl}\f$ object.
        * 
        * Enables instantiation like:
        * 
        *     Lielab::domain::gl x(3), y(4), z(5);
        * 
        * @param[in] shape The shape of the data matrix.
        */

        this->_data = Eigen::MatrixXd::Zero(shape,shape);
        this->shape = shape;
    }

    gl(std::initializer_list<double> other) 
    {
        /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{gl} \f}
        *
        * Constructor instantiating an \f$\mathfrak{gl}\f$ object from either a
        * \f$n \times 1\f$ real vector.
        *
        * @param[in] other The object to instantiate from as a real vector.
        */

        this->set_vector(Eigen::VectorXd{std::move(other)});
    }

    template<typename OtherDerived>
    gl(const Eigen::MatrixBase<OtherDerived> & other)
    {
        /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& \mathfrak{gl} \\ (\mathbb{R}^{n \times 1}) &\rightarrow& \mathfrak{gl} \f}
        *
        * Constructor instantiating an \f$\mathfrak{gl}\f$ object from either an
        * \f$n \times n\f$ real matrix.
        *
        * @param[in] other The object to instantiate from as a real matrix.
        */

        if (other.rows() == other.cols())
        {
            shape = other.rows();
            _data = other; 
        }
        else
        {
            throw Errorx("Input data to gl malformed.");
        }
    }

    template<typename OtherDerived>
    gl & operator=(const Eigen::MatrixBase<OtherDerived> & other)
    {
        /*! \f{eqnarray*}{ \mathfrak{gl} &:= \mathbb{R}^{n \times n}\f}
        * 
        * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `gl`.
        */

        if (other.rows() == other.cols())
        {
            shape = other.rows();
            _data = other;
        }
        else
        {
            throw Errorx("Input data to gl malformed.");
        }
        return *this;
    }

    static gl basis(const size_t i, const size_t shape) 
    {
        /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{gl} \f}
        *
        * Returns the i'th basis element of the gl algebra.
        * 
        * @param[in] i The basis vector.
        * @param[in] shape The shape of the algebra.
        * @param[out] out The gl element.
        */

        gl out(shape);
        const size_t _dim = out.get_dimension();
        Eigen::VectorXd _u(_dim);
        _u = Eigen::VectorXd::Zero(_dim);
        _u(i) = 1.0;
        out.set_vector(_u);
        return out;
    }

    static Eigen::MatrixXd project(const Eigen::MatrixXd & other)
    {
        /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{gl} \f}
        *
        * Projects a matrix suitable for data.
        */

        if (other.rows() != other.cols())
        {
            throw Errorx("Size of the matrix must be square.");
        }
        return other;
    }

    size_t get_dimension() const
    {
        /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
        * 
        * Gets the dimension of the algebra.
        */

        return static_cast<size_t>(std::pow(this->shape, 2));
    }

    Eigen::VectorXd get_vector() const
    {
        /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
        * 
        * Returns the vector representation.
        */
        
        return Eigen::Map<const Eigen::VectorXd>(_data.data(), _data.size());
    }

    Eigen::MatrixXd get_matrix() const
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

        return _data;
    }

    void set_vector(const Eigen::VectorXd & vector) 
    {
        /*! \f{equation*}{ \mathfrak{gl} := \mathbb{R}^{n \times 1} \f}
        * 
        * @param[in] vector An Eigen::VectorXd to assign.
        */

        this->shape = static_cast<size_t>(std::sqrt(vector.size()));
        this->_data = Eigen::Map<const Eigen::MatrixXd>(vector.data(), shape, shape);
    }

    double operator()(const size_t index)
    {
        /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
        *
        * Gets a value in the column vector representation.
        */

        Eigen::VectorXd vector = get_vector();
        return vector(index);
    }

    // double & operator()(const size_t index) // TODO [not important]
    // {
    //     /*! \f{equation*}{ \mathfrak{gl}(\mathbb{Z}) := \mathbb{R} \f}
    //     *
    //     * Assignment of a value in the column vector representation.
    //     */
        
    //     return _data(index);
    // }

    double operator()(const size_t index1, const size_t index2) const
    {
        /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
        *
        * Gets a value in the square matrix representation.
        */
        
        return _data(index1, index2);
    }

    double & operator()(const size_t index1, const size_t index2)
    {
        /*! \f{equation*}{ \mathfrak{gl}(\mathbb{Z}, \mathbb{Z}) := \mathbb{R} \f}
        *
        * Assignment of a value in the square matrix representation.
        */
        
        return _data(index1, index2);
    }

    gl & operator+=(const gl & other)
    {
        /*! \f{equation*}{ (\mathfrak{gl}, \mathfrak{gl}) \rightarrow \mathfrak{gl} \f}
        *
        * In place addition of two vectors in the algebra.
        */

        assert(this->shape == other.shape);
        this->_data += other._data;
        return *this;
    }

    gl & operator-=(const gl & other)
    {
        /*! \f{equation*}{ (\mathfrak{gl}, \mathfrak{gl}) \rightarrow \mathfrak{gl} \f}
        *
        * In place subtraction of two vectors in the algebra.
        */

        assert(this->shape == other.shape);
        this->_data -= other._data;
        return *this;
    }

    gl operator-() const
    {
        /*! \f{equation*}{ (\mathfrak{gl}) \rightarrow \mathfrak{gl} \f}
        *
        * Unary negative of the vector.
        */

        return -this->_data;
    }

    gl operator*(const Real auto other) const
    {
        /*! \f{equation*}{ (\mathfrak{gl}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
        *
        * Scalar product.
        */

        return this->_data * other;
    }

    gl & operator*=(const Real auto other)
    {
        /*! \f{equation*}{ (\mathfrak{gl}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
        *
        * In place scalar product.
        */

        this->_data *= other;
        return *this;
    }

    gl & operator*=(const gl & other)
    {
        /*! \f{equation*}{ (\mathfrak{gl}, \mathfrak{gl}) \rightarrow \mathbb{R}^{n \times n} \notin \mathfrak{gl} \f}
        *
        * In place product of two vectors in the algebra.
        */

        assert(this->shape == other.shape);
        this->_data *= other._data;
        return *this;
    }

    gl operator/(const Real auto other) const
    {
        /*! \f{equation*}{ (\mathfrak{gl}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
        *
        * Scalar division.
        */

        return this->_data / other;
    }

    gl & operator/=(const Real auto other)
    {
        /*! \f{equation*}{ (\mathfrak{gl}, \mathbb{R}) \rightarrow \mathfrak{gl} \f}
        *
        * In place scalar division.
        */

        this->_data /= other;
        return *this;
    }

    friend std::ostream& operator<<(std::ostream & os, const gl & other);
};

gl operator*(const Real auto other, const gl & rhs)
{
    /*! \f{equation*}{ (\mathbb{R}, \mathfrak{gl}) \rightarrow \mathfrak{gl} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

std::ostream& operator<<(std::ostream & os, const gl & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */
    
    os << static_cast<const Eigen::MatrixXd>(other._data);
    return os;
}
}
}

#endif
