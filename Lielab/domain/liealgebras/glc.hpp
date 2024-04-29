#ifndef _LIELAB_DOMAIN_LIEALGEBRAS_glc_HPP
#define _LIELAB_DOMAIN_LIEALGEBRAS_glc_HPP

#include "../../abstract/abstract_all.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
namespace domain
{

using Lielab::abstract::Imaginary;

class glc
{
    public:

    /*! Flag specifying whether or not the algebra is abelian. */
    static constexpr bool abelian = false;
    size_t shape = 0;
    

    /*! Internal data for the algebra. */
    Eigen::MatrixXcd _data;

    glc() : glc(0)
    {
        /*! \f{equation*}{() \rightarrow \mathfrak{glc} \f}
        * 
        * Empty initialization function. Enables instantiation like:
        * 
        *     Lielab::domain::glc x, y, z;
        * 
        */

    }

    glc(const size_t shape)
    {
        /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{glc} \f}
        *
        * Constructor instantiating an \f$\mathfrak{glc}\f$ object.
        * 
        * Enables instantiation like:
        * 
        *     Lielab::domain::glc x(3), y(4), z(5);
        * 
        * @param[in] shape The shape of the data matrix.
        */

        this->_data = Eigen::MatrixXcd::Zero(shape,shape);
        this->shape = shape;
    }

    glc(std::initializer_list<std::complex<double>> other) 
    {
        /*! \f{equation}{(\mathbb{C}^{n \times 1}) \rightarrow \mathfrak{glc} \f}
        *
        * Constructor instantiating a \f$\mathfrak{glc}\f$ object from either an
        * \f$n \times 1\f$ imaginary vector.
        *
        * @param[in] other The object to instantiate from as an imaginary vector.
        */

        const Eigen::VectorXcd vector_imag = Eigen::VectorXcd{std::move(other)};
        const ptrdiff_t sz = vector_imag.size();
        Eigen::VectorXd vector_real = Eigen::VectorXd::Zero(2*sz);
        ptrdiff_t kk = 0;        
        for (ptrdiff_t ii = 0; ii < sz; ii++)
        {
            vector_real(kk) = std::real(vector_imag(ii));
            vector_real(kk+1) = std::imag(vector_imag(ii));
            kk += 2;
        }

        this->set_vector(vector_real);
    }

    template<typename OtherDerived>
    glc(const Eigen::MatrixBase<OtherDerived> & other)
    {
        /*! \f{eqnarray*}{(\mathbb{C}^{n \times n}) &\rightarrow& \mathfrak{glc} \\ (\mathbb{C}^{n \times 1}) &\rightarrow& \mathfrak{glc} \f}
        *
        * Constructor instantiating an \f$\mathfrak{glc}\f$ object from an
        * \f$n \times n\f$ imaginary matrix.
        *
        * @param[in] other The object to instantiate from as an imaginary matrix.
        */

        if (other.rows() == other.cols())
        {
            shape = other.rows();
            _data = other; 
        }
        else
        {
            throw Errorx("Input data to glc malformed.");
        }
    }

    template<typename OtherDerived>
    glc & operator=(const Eigen::MatrixBase<OtherDerived> & other)
    {
        /*! \f{eqnarray*}{ \mathfrak{glc} &:= \mathbb{C}^{n \times n}\f}
        * 
        * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `glc`.
        */

        if (other.rows() == other.cols())
        {
            shape = other.rows();
            _data = other;
        }
        else
        {
            throw Errorx("Input data to u malformed.");
        }
        return *this;
    }

    static glc basis(const size_t i, const size_t shape) 
    {
        /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{glc} \f}
        *
        * Returns the i'th basis element of the glc algebra.
        * 
        * @param[in] i The basis vector.
        * @param[in] shape The shape of the algebra.
        * @param[out] out The glc element.
        */

        glc out(shape);
        const size_t _dim = out.get_dimension();
        Eigen::VectorXd _u(_dim);
        _u = Eigen::VectorXd::Zero(_dim);
        _u(i) = 1.0;
        out.set_vector(_u);
        return out;
    }

    static Eigen::MatrixXcd project(const Eigen::MatrixXcd & other)
    {
        /*! \f{equation*}{ (\mathbb{C}^{n \times n}) \rightarrow \mathbb{C}^{n \times n} \in \mathfrak{glc} \f}
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

        return static_cast<size_t>(2*std::pow(this->shape, 2));
    }

    Eigen::VectorXd get_vector() const
    {
        /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
        * 
        * Returns the vector representation.
        */
        
        const ptrdiff_t dim = this->get_dimension();
        const Eigen::MatrixXcd A = this->get_ados_representation();

        Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);
        ptrdiff_t kk = 0;
        for (ptrdiff_t ii = 0; ii < this->shape; ii++)
        {
            for (ptrdiff_t jj = 0; jj < this->shape; jj++)
            {
                out(kk) = std::real(A(ii,jj));
                out(kk+1) = std::imag(A(ii,jj));
                kk = kk + 2;
            }
        }

        return out;
    }

    Eigen::MatrixXcd get_ados_representation() const
    {
        /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times n} \f}
        * 
        * Returns a matrix representation.
        */

        return _data;
    }

    void set_vector(const Eigen::VectorXd & vector) 
    {
        /*! \f{equation*}{ \mathfrak{glc} := \mathbb{R}^{n \times 1} \f}
        * 
        * @param[in] vector An Eigen::VectorXd to assign.
        */

        const ptrdiff_t dim = vector.size();
        this->shape = static_cast<size_t>(std::sqrt(dim/2));
        this->_data = Eigen::MatrixXcd::Zero(this->shape, this->shape);
        
        ptrdiff_t kk = 0;
        for (ptrdiff_t ii = 0; ii < this->shape; ii++)
        {
            for (ptrdiff_t jj = 0; jj < this->shape; jj++)
            {
                this->_data(ii, jj) = std::complex<double>(vector(kk), vector(kk+1));
                kk = kk + 2;
            }
        }
    }

    double operator()(const size_t index)
    {
        /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
        *
        * Gets a value in the column vector representation.
        */

        const Eigen::VectorXd vector = get_vector();
        return vector(index);
    }

    // double & operator()(const size_t index) // TODO [not important]
    // {
    //     /*! \f{equation*}{ \mathfrak{glc}(\mathbb{Z}) := \mathbb{R} \f}
    //     *
    //     * Assignment of a value in the column vector representation.
    //     */
        
    //     return _data(index);
    // }

    std::complex<double> operator()(const size_t index1, const size_t index2) const
    {
        /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
        *
        * Gets a value in the square matrix representation.
        */
        
        return _data(index1, index2);
    }

    std::complex<double> & operator()(const size_t index1, const size_t index2)
    {
        /*! \f{equation*}{ \mathfrak{glc}(\mathbb{Z}, \mathbb{Z}) := \mathbb{C} \f}
        *
        * Assignment of a value in the square matrix representation.
        */
        
        return _data(index1, index2);
    }

    glc operator+(const glc & other) const
    {
        /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
        *
        * Addition of two vectors in the algebra.
        */

        assert(this->shape == other.shape);
        return this->_data + other._data;
    }

    glc & operator+=(const glc & other)
    {
        /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
        *
        * In place addition of two vectors in the algebra.
        */

        assert(this->shape == other.shape);
        this->_data += other._data;
        return *this;
    }

    glc operator-(const glc & other) const
    {
        /*! \f{equation*}{ (\mathfrak{glc}}, \mathfrak{glc}) \rightarrow \mathfrak{glc}} \f}
        *
        * Subtraction of two vectors in the algebra.
        */
        
        assert(this->shape == other.shape);
        return this->_data - other._data;
    }

    glc & operator-=(const glc & other)
    {
        /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
        *
        * In place subtraction of two vectors in the algebra.
        */

        assert(this->shape == other.shape);
        this->_data -= other._data;
        return *this;
    }

    glc operator-() const
    {
        /*! \f{equation*}{ (\mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
        *
        * Unary negative of the vector.
        */

        return -this->_data;
    }

    glc operator*(const Imaginary auto other) const
    {
        /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{C}) \rightarrow \mathfrak{glc} \f}
        *
        * Scalar product.
        */

        return this->_data * other;
    }

    glc & operator*=(const Imaginary auto other)
    {
        /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{C}) \rightarrow \mathfrak{glc} \f}
        *
        * In place scalar product.
        */

        this->_data *= other;
        return *this;
    }

    glc operator*(const glc & other) const
    {
        /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
        *
        * Vector product.
        */

        assert(this->shape == other.shape);
        return this->_data * other._data;
    }

    glc & operator*=(const glc & other)
    {
        /*! \f{equation*}{ (\mathfrak{glc}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
        *
        * In place product of two vectors in the algebra.
        */

        assert(this->shape == other.shape);
        this->_data *= other._data;
        return *this;
    }

    glc operator/(const Imaginary auto other) const
    {
        /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{C}) \rightarrow \mathfrak{glc} \f}
        *
        * Scalar division.
        */

        return this->_data / other;
    }

    glc & operator/=(const Imaginary auto other)
    {
        /*! \f{equation*}{ (\mathfrak{glc}, \mathbb{C}) \rightarrow \mathfrak{glc} \f}
        *
        * In place scalar division.
        */

        this->_data /= other;
        return *this;
    }

    friend std::ostream& operator<<(std::ostream & os, const glc & other);
};

glc operator*(const Imaginary auto other, const glc & rhs)
{
    /*! \f{equation*}{ (\mathbb{C}, \mathfrak{glc}) \rightarrow \mathfrak{glc} \f}
    *
    * Scalar product.
    */

    return rhs*other;
}

std::ostream& operator<<(std::ostream & os, const glc & other)
{
    /*!
    * Overloads the "<<" stream insertion operator.
    */
    
    os << static_cast<const Eigen::MatrixXcd>(other._data);
    return os;
}
}
}

#endif
