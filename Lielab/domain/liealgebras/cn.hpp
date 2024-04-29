#ifndef _LIELAB_DOMAIN_LIEALGEBRAS_cn_HPP
#define _LIELAB_DOMAIN_LIEALGEBRAS_cn_HPP

#include "../../abstract/abstract_all.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
    namespace domain
    {
        using Lielab::abstract::Imaginary;

        class cn
        {
            public:

            /*! Flag specifying whether or not the algebra is abelian. */
            static constexpr bool abelian = true;
            size_t shape = 0;
            
            /*! Internal data for the algebra. */
            Eigen::VectorXcd _data;

            cn() : cn(0)
            {
                /*! \f{equation*}{() \rightarrow \mathfrak{cn} \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     Lielab::domain::cn x, y, z;
                * 
                */

            }

            cn(const size_t shape)
            {
                /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{cn} \f}
                *
                * Constructor instantiating an \f$\mathfrak{cn}\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     Lielab::domain::cn x(3), y(4), z(5);
                * 
                * @param[in] shape The shape of the data matrix.
                */

                this->shape = shape;
                this->_data = Eigen::VectorXcd::Zero(std::max((int)shape-1, 0));
            }

            cn(std::initializer_list<std::complex<double>> other)
            {
                /*! \f{equation}{(\mathbb{C}^{n \times 1}) \rightarrow \mathfrak{cn} \f}
                *
                * Constructor instantiating an \f$\mathfrak{cn}\f$ object from either a
                * \f$n \times 1\f$ imaginary vector.
                *
                * @param[in] other The object to instantiate from as an imaginary vector.
                */

                this->_data = Eigen::VectorXcd{std::move(other)};
                this->shape = this->_data.size() + 1;
            }

            template<typename OtherDerived>
            cn(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{(\mathbb{C}^{n \times n}) &\rightarrow& \mathfrak{cn} \\ (\mathbb{C}^{n \times 1}) &\rightarrow& \mathfrak{cn} \f}
                *
                * Constructor instantiating an \f$\mathfrak{cn}\f$ object from either an
                * \f$n \times n\f$ imaginary matrix or \f$n \times 1\f$ imaginary vector.
                *
                * @param[in] other The object to instantiate from as an imaginary matrix.
                */

                if (other.rows() == other.cols() && other.rows() != 1)
                {
                    this->shape = other.rows();
                    this->_data = Eigen::VectorXcd::Zero(shape-1);
                    for (size_t ii = 0; ii < shape-1; ii++)
                    {
                        _data(ii) = other(ii, shape-1);
                    }
                }
                else if (other.cols() == 1)
                {
                    this->shape = other.rows()+1;
                    this->_data = other;
                }
                else
                {
                    throw Errorx("Input data to cn malformed.");
                }
            }

            template<typename OtherDerived>
            cn & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{ \mathfrak{cn} &:= \mathbb{C}^{n \times n} \\ &:= \mathbb{C}^{n \times 1} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `cn`.
                */

                if (other.rows() == other.cols() && other.rows() != 1)
                {
                    this->shape = other.rows();
                    this->_data = Eigen::VectorXcd::Zero(shape-1);
                    for (int ii = 0; ii < shape-1; ii++)
                    {
                        _data(ii) = other(ii, shape-1);
                    }
                }
                else if (other.cols() == 1)
                {
                    this->shape = other.rows()+1;
                    this->_data = other;
                }
                else
                {
                    throw Errorx("Input data to cn malformed.");
                }
                return *this;
            }

            static cn basis(const size_t i, const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{cn} \f}
                *
                * Returns the i'th basis element of the cn algebra.
                * 
                * @param[in] i The basis vector.
                * @param[in] shape The shape of the algebra.
                * @param[out] out The cn element. 
                */

                cn out(shape);
                const size_t _dim = out.get_dimension();
                Eigen::VectorXd _u(_dim);
                _u = Eigen::VectorXd::Zero(_dim);
                _u(i) = 1.0;
                out.set_vector(_u);
                return out;
            }

            static Eigen::MatrixXcd project(const Eigen::MatrixXcd & other)
            {
                /*! \f{equation*}{ (\mathbb{C}^{n \times n}) \rightarrow \mathbb{C}^{n \times n} \in \mathfrak{cn} \f}
                *
                * Projects a matrix suitable for data.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }

                Eigen::MatrixXcd out = other;

                for (int ii = 0; ii < out.rows(); ii++)
                {
                    for (int jj = 0; jj < out.cols(); jj++)
                    {
                        if (jj != other.cols()-1)
                        {
                            out(ii,jj) = 0.0;
                        }
                    }
                }
                out(out.rows()-1, out.cols()-1) = 0.0;
                return out;
            }

            size_t get_dimension() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
                * 
                * Gets the dimension of the algebra.
                */

                return 2*(this->shape - 1);
            }

            Eigen::VectorXd get_vector() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times 1} \f}
                * 
                * Returns the vector representation.
                */

                const size_t dim = this->get_dimension();
                Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);

                ptrdiff_t kk = 0;
                for (ptrdiff_t ii = 0; ii < dim/2; ii++)
                {
                    out(kk) = std::real(this->_data(ii));
                    kk += 1;
                    out(kk) = std::imag(this->_data(ii));
                    kk += 1;
                }

                return out;
            }

            Eigen::MatrixXcd get_ados_representation() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{C}^{n \times n} \f}
                * 
                * Returns a matrix representation.
                */

                Eigen::MatrixXcd out = Eigen::MatrixXcd::Zero(shape, shape);

                for (int ii = 0; ii < shape-1; ii++)
                {
                    out(ii, shape-1) = _data(ii);
                }

                return out;
            }

            void set_vector(const Eigen::VectorXd & vector)
            {
                /*! \f{equation*}{ \mathfrak{cn} := \mathbb{C}^{n \times 1} \f}
                * 
                * @param[in] vector An Eigen::VectorXd to assign.
                */

                const size_t dim = vector.size();
                Eigen::VectorXcd temp = Eigen::VectorXcd::Zero(static_cast<ptrdiff_t>(dim/2.0));

                ptrdiff_t kk = 0;
                for (ptrdiff_t ii = 0; ii < static_cast<ptrdiff_t>(dim/2); ii++)
                {
                    temp(ii) = std::complex<double>(vector(kk), vector(kk+1));
                    kk += 2;
                }

                this->_data = temp;
                this->shape = this->_data.size() + 1;
            }

            double operator()(const size_t index) const
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{C} \f}
                *
                * Gets a value in the column vector representation.
                */

                return get_vector()(index);
            }

            std::complex<double> & operator()(const size_t index)
            {
                /*! \f{equation*}{ \mathfrak{cn}(\mathbb{Z}) := \mathbb{C} \f}
                *
                * Assignment of a value in the column vector representation.
                */
                
                return _data(index);
            }

            std::complex<double> operator()(const size_t index1, const size_t index2) const
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{C} \f}
                *
                * Gets a value in the square matrix representation.
                */
                
                Eigen::MatrixXcd mat = get_ados_representation();
                return mat(index1, index2);
            }

            cn operator+(const cn & other) const
            {
                /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
                *
                * Addition of two vectors in the algebra.
                */
                
                assert(this->shape == other.shape);
                return this->_data + other._data;
            }

            cn & operator+=(const cn & other)
            {
                /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
                *
                * In place addition of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                this->_data += other._data;
                return *this;
            }

            cn operator-(const cn & other) const
            {
                /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
                *
                * Subtraction of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                return this->_data - other._data;
            }

            cn & operator-=(const cn & other)
            {
                /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
                *
                * In place subtraction of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                this->_data -= other._data;
                return *this;
            }

            cn operator-() const
            {
                /*! \f{equation*}{ (\mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
                *
                * Unary negative of the vector.
                */

                return -this->_data;
            }

            cn operator*(const Imaginary auto other) const
            {
                /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
                *
                * Scalar product.
                */

                return this->_data * other;
            }

            cn & operator*=(const Imaginary auto other)
            {
                /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
                *
                * In place scalar product.
                */

                this->_data *= other;
                return *this;
            }

            cn operator*(const cn & other) const
            {
                /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathbb{C}^{n \times n} \notin \mathfrak{cn} \f}
                *
                * Vector product.
                */

                assert(this->shape == other.shape);
                return Eigen::VectorXcd::Zero(shape-1);
            }

            cn & operator*=(const cn & other)
            {
                /*! \f{equation*}{ (\mathfrak{cn}, \mathfrak{cn}) \rightarrow \mathbb{C}^{n \times n} \notin \mathfrak{cn} \f}
                *
                * In place vector product.
                */

                assert(this->shape == other.shape);
                this->_data *= 0.0;
                return *this;
            }

            cn operator/(const Imaginary auto other) const
            {
                /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
                *
                * Scalar division.
                */

                return this->_data / other;
            }

            cn & operator/=(const Imaginary auto other)
            {
                /*! \f{equation*}{ (\mathfrak{cn}, \mathbb{C}) \rightarrow \mathfrak{cn} \f}
                *
                * In place scalar division.
                */

                this->_data /= other;
                return *this;
            }

            friend std::ostream& operator<<(std::ostream & os, const cn & other);
        };

        cn operator*(const Imaginary auto other, const cn & rhs)
        {
            /*! \f{equation*}{ (\mathbb{C}, \mathfrak{cn}) \rightarrow \mathfrak{cn} \f}
            *
            * Scalar product.
            */

            return rhs*other;
        }

        std::ostream& operator<<(std::ostream & os, const cn & other)
        {
            /*!
            * Overloads the "<<" stream insertion operator.
            */
            
            os << static_cast<const Eigen::VectorXcd>(other._data);
            return os;
        }
    }
}

#endif
