#ifndef _LIELAB_DOMAIN_LIEALGEBRAS_so_HPP
#define _LIELAB_DOMAIN_LIEALGEBRAS_so_HPP

#include "../../abstract/abstract_all.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
    namespace domain
    {
        class so
        {
            /*!
            * The special orthonormal (so) algebra class.
            */

            public:
            static constexpr bool abelian = false;
            size_t shape = 0;
            
            Eigen::MatrixXd _data;

            so() : so(0)
            {
                /*! \f{equation*}{() \rightarrow \mathfrak{so} \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     Lielab::domain::so x, y, z;
                * 
                */

            }

            so(const size_t shape)
            {
                /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{so} \f}
                *
                * Constructor instantiating an \f$\mathfrak{so}\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     Lielab::domain::so x(3), y(4), z(5);
                * 
                * @param[in] shape The shape of the data matrix.
                */
                
                this->_data = Eigen::MatrixXd::Zero(shape, shape);
                this->shape = shape;
            }

            so(std::initializer_list<double> other)
            {
                /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{so} \f}
                *
                * Constructor instantiating an \f$\mathfrak{so}\f$ object from either a
                * \f$n \times 1\f$ real vector.
                *
                * @param[in] other The object to instantiate from as a real vector.
                */

                this->set_vector(Eigen::VectorXd{std::move(other)});
            }

            template<typename OtherDerived>
            so(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& \mathfrak{so} \\ (\mathbb{R}^{n \times 1}) &\rightarrow& \mathfrak{so} \f}
                *
                * Constructor instantiating an \f$\mathfrak{so}\f$ object from either an
                * \f$n \times n\f$ real matrix or \f$n \times 1\f$ real vector.
                *
                * @param[in] other The object to instantiate from as a real matrix.
                */

                if (other.rows() == other.cols())
                {
                    this->_data = Eigen::MatrixXd(other);
                    this->shape = _data.rows();
                }
                else if (other.cols() == 1)
                {
                    this->set_vector(other);
                }
                else
                {
                    throw Errorx("Input data to so malformed.");
                }

            }

            template<typename OtherDerived>
            so & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{ \mathfrak{so} &:= \mathbb{R}^{n \times n} \\ &:= \mathbb{R}^{n \times 1} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `so`.
                */

                if (other.rows() == other.cols())
                {
                    this->_data = Eigen::MatrixXd(other);
                    this->shape = _data.rows();
                }
                else if (other.cols() == 1)
                {
                    this->set_vector(other);
                }
                else
                {
                    throw Errorx("Input data to so malformed.");
                }

                return *this;
            }

            static so basis(const size_t i, const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{so} \f}
                *
                * Returns the i'th basis element of the so algebra.
                * 
                * @param[in] i The basis vector.
                * @param[in] shape The shape of the algebra.
                * @param[out] out The so element.
                */

                so out(shape);
                const size_t _dim = out.get_dimension();
                Eigen::VectorXd _u(_dim);
                _u = Eigen::VectorXd::Zero(_dim);
                _u(i) = 1.0;
                out.set_vector(_u);
                return out;
            }

            static Eigen::MatrixXd project(const Eigen::MatrixXd & other)
            {
                /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{so} \f}
                *
                * Projects a matrix suitable for data.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }

                return (other - other.transpose())/2;
            }

            size_t get_dimension() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
                * 
                * Gets the dimension of the algebra.
                */

                return this->shape * (this->shape - 1) / 2;
            }

            Eigen::VectorXd get_vector() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
                * 
                * Returns the vector representation.
                */

                const size_t dim = this->get_dimension();
                Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);
                int k = 0;

                for (size_t ii = shape-1; ii > 0; ii--)
                {
                    for (size_t jj = shape; jj > ii; jj--)
                    {
                        out(k) = _data(ii-1, jj-1)/std::pow(-1.0, ii+jj);
                        k++;
                    }
                }

                return out;
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
                /*! \f{equation*}{ \mathfrak{so} := \mathbb{R}^{n \times 1} \f}
                * 
                * @param[in] vector An Eigen::VectorXd to assign.
                */

                this->shape = static_cast<size_t>(std::sqrt(2.0*vector.size() + 0.25) + 0.5);
                Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(shape, shape);
                int k = 0;

                for (size_t ii = shape-1; ii > 0; ii--)
                {
                    for (size_t jj = shape; jj > ii; jj--)
                    {
                        temp(ii-1, jj-1) = std::pow(-1, (ii+jj))*vector(k);
                        k++;
                    }
                }

                this->_data = (temp - temp.transpose());
            }

            double operator()(const size_t index) const
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
                *
                * Gets a value in the column vector representation.
                */

                Eigen::VectorXd vector = get_vector();
                return vector(index);
            }

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
                /*! \f{equation*}{ \mathfrak{so}(\mathbb{Z}, \mathbb{Z}) := \mathbb{R} \f}
                *
                * Assignment of a value in the square matrix representation.
                */

                return _data(index1, index2);
            }

            so & operator+=(const so & other)
            {
                /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}) \rightarrow \mathfrak{so} \f}
                *
                * In place addition of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                this->_data += other._data;
                return *this;
            }

            so & operator-=(const so & other)
            {
                /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}) \rightarrow \mathfrak{so} \f}
                *
                * In place subtraction of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                this->_data -= other._data;
                return *this;
            }

            so operator-() const
            {
                /*! \f{equation*}{ (\mathfrak{so}) \rightarrow \mathfrak{so} \f}
                *
                * Unary negative of the vector.
                */
                
                return -this->_data;
            }

            so operator*(const Real auto other) const
            {
                /*! \f{equation*}{ (\mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{so} \f}
                *
                * Scalar product.
                */

                so out = this->_data * other;
                return out;
            }

            so & operator*=(const Real auto other)
            {
                /*! \f{equation*}{ (\mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{so} \f}
                *
                * In place scalar multiplication.
                */

                this->_data *= other;
                return *this;
            }

            so & operator*=(const so & other)
            {
                /*! \f{equation*}{ (\mathfrak{so}, \mathfrak{so}) \rightarrow \mathbb{R}^{n \times n} \notin \mathfrak{so} \f}
                *
                * In place vector product.
                */

                assert(this->shape == other.shape);
                this->_data *= other._data;
                return *this;
            }

            so operator/(const Real auto other) const
            {
                /*! \f{equation*}{ (\mathfrak{so}, \mathbb{R}) \rightarrow \mathfrak{so} \f}
                *
                * Scalar division.
                */

                Eigen::MatrixXd out = _data / other;
                return out;
            }

            so & operator/=(const Real auto other)
            {
                /*! \f{equation*}{ (\mathfrak{so}, \mathbb{RR}) \rightarrow \mathfrak{so} \f}
                *
                * In place scalar division.
                */

                this->_data /= other;
                return *this;
            }

            friend std::ostream & operator<<(std::ostream& os, const so & other);
        };

        so operator*(const Real auto other, const so & rhs)
        {
            /*! \f{equation*}{ (\mathbb{R}, \mathfrak{so}) \rightarrow \mathfrak{so} \f}
            *
            * Scalar product.
            */

            return rhs*other;
        }

        std::ostream & operator<<(std::ostream & os, const so & other)
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
