#ifndef _LIELAB_DOMAIN_se_HPP
#define _LIELAB_DOMAIN_se_HPP

#include "./rn.hpp"
#include "./so.hpp"

namespace lielab
{
    namespace domain
    {
        class se
        {
            /*!
            * The special Euclidean (se) algebra class.
            *
            * Defined as the interacting product rn x so
            */

            public:
            static constexpr bool abelian = false;
            size_t shape = 0;
            
            Eigen::MatrixXd _data;

            se() : se(0)
            {
                /*! \f{equation*}{() \rightarrow \mathfrak{se} \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     lielab::domain::se x, y, z;
                * 
                */

            }

            se(const size_t shape)
            {
                /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{se} \f}
                *
                * Constructor instantiating an \f$\mathfrak{se}\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     lielab::domain::se x(3), y(4), z(5);
                * 
                * @param[in] shape The shape of the data matrix.
                */
                
                this->_data = Eigen::MatrixXd::Zero(shape, shape);
                this->shape = shape;
            }

            se(std::initializer_list<double> other)
            {
                /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{se} \f}
                *
                * Constructor instantiating an \f$\mathfrak{se}\f$ object from either a
                * \f$n \times 1\f$ real vector.
                *
                * @param[in] other The object to instantiate from as a real vector.
                */

                this->set_vector(Eigen::VectorXd{std::move(other)});
            }

            template<typename OtherDerived>
            se(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& \mathfrak{se} \\ (\mathbb{R}^{n \times 1}) &\rightarrow& \mathfrak{se} \f}
                *
                * Constructor instantiating an \f$\mathfrak{se}\f$ object from either an
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
                    throw Errorx("Input data to se malformed.");
                }

            }

            template<typename OtherDerived>
            se & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{ \mathfrak{se} &:= \mathbb{R}^{n \times n} \\ &:= \mathbb{R}^{n \times 1} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `se`.
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
                    throw Errorx("Input data to se malformed.");
                }

                return *this;
            }

            static se basis(const size_t i, const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{se} \f}
                *
                * Returns the i'th basis element of the se algebra.
                * 
                * @param[in] i The basis vector.
                * @param[in] shape The shape of the algebra.
                * @param[out] out The se element.
                */

                se out(shape);
                const size_t _dim = out.get_dimension();
                Eigen::VectorXd _u(_dim);
                _u = Eigen::VectorXd::Zero(_dim);
                _u(i) = 1.0;
                out.set_vector(_u);
                return out;
            }

            // static Eigen::MatrixXd project(const Eigen::MatrixXd & other)
            // {
            //     /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{se} \f}
            //     *
            //     * Projects a matrix suitable for data.
            //     */

            //     if (other.rows() != other.cols())
            //     {
            //         throw Errorx("Size of the matrix must be square.");
            //     }

            //     return (other - other.transpose())/2;
            // }

            size_t get_dimension() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
                * 
                * Gets the dimension of the algebra.
                */

                return (this->shape - 1) * (this->shape - 2) / 2 + this->shape - 1;
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

                for (size_t ii = 0; ii < this->shape - 1; ii++)
                {
                    out(k) = this->_data(ii, this->shape - 1);
                    k++;
                }

                for (size_t ii = shape-2; ii > 0; ii--)
                {
                    for (size_t jj = shape-1; jj > ii; jj--)
                    {
                        out(k) = _data(ii-1, jj-1)/std::pow(-1.0, ii+jj);
                        k++;
                    }
                }

                return out;
            }

            Eigen::MatrixXd get_ados_representation() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times n} \f}
                * 
                * Returns a matrix representation.
                */
               
                return _data;
            }

            void set_vector(const Eigen::VectorXd & vector)
            {
                /*! \f{equation*}{ \mathfrak{se} := \mathbb{R}^{n \times 1} \f}
                * 
                * @param[in] vector An Eigen::VectorXd to assign.
                */

                this->shape = static_cast<size_t>(std::sqrt(8.0*vector.size() + 1.0)/2.0 + 0.5);
                Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(shape, shape);
                int k = 0;

                for (size_t ii = 0; ii < this->shape - 1; ii++)
                {
                    temp(ii, this->shape-1) = vector(k);
                    k++;
                }

                for (size_t ii = shape-2; ii > 0; ii--)
                {
                    for (size_t jj = shape-1; jj > ii; jj--)
                    {
                        temp(ii-1, jj-1) =  std::pow(-1, (ii+jj))*vector(k);
                        temp(jj-1, ii-1) = -std::pow(-1, (ii+jj))*vector(k);
                        k++;
                    }
                }

                this->_data = temp;
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
                /*! \f{equation*}{ \mathfrak{se}(\mathbb{Z}, \mathbb{Z}) := \mathbb{R} \f}
                *
                * Assignment of a value in the square matrix representation.
                */

                return _data(index1, index2);
            }

            se operator+(const se & other) const
            {
                /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
                *
                * Addition of two vectors in the algebra.
                */
                
                assert(this->shape == other.shape);
                se out = this->_data + other._data;
                return out;
            }

            se & operator+=(const se & other)
            {
                /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
                *
                * In place addition of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                this->_data += other._data;
                return *this;
            }

            se operator-(const se & other) const
            {
                /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
                *
                * Subtraction of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                Eigen::MatrixXd out = this->_data - other._data;
                return out;
            }

            se & operator-=(const se & other)
            {
                /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
                *
                * In place subtraction of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                this->_data -= other._data;
                return *this;
            }

            se operator-() const
            {
                /*! \f{equation*}{ (\mathfrak{se}) \rightarrow \mathfrak{se} \f}
                *
                * Unary negative of the vector.
                */
                
                return -this->_data;
            }

            se operator*(const Real auto other) const
            {
                /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
                *
                * Scalar product.
                */

                se out = this->_data * other;
                return out;
            }

            se & operator*=(const Real auto other)
            {
                /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
                *
                * In place scalar multiplication.
                */

                this->_data *= other;
                return *this;
            }

            se operator*(const se & other) const
            {
                /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathbb{R}^{n \times n} \notin \mathfrak{se} \f}
                *
                * Vector product.
                */

                assert(this->shape == other.shape);
                se out = this->_data * other._data;
                return out;
            }

            se & operator*=(const se & other)
            {
                /*! \f{equation*}{ (\mathfrak{se}, \mathfrak{se}) \rightarrow \mathbb{R}^{n \times n} \notin \mathfrak{se} \f}
                *
                * In place vector product.
                */

                assert(this->shape == other.shape);
                this->_data *= other._data;
                return *this;
            }

            se operator/(const Real auto other) const
            {
                /*! \f{equation*}{ (\mathfrak{se}, \mathbb{R}) \rightarrow \mathfrak{se} \f}
                *
                * Scalar division.
                */

                Eigen::MatrixXd out = _data / other;
                return out;
            }

            se & operator/=(const Real auto other)
            {
                /*! \f{equation*}{ (\mathfrak{se}, \mathbb{RR}) \rightarrow \mathfrak{se} \f}
                *
                * In place scalar division.
                */

                this->_data /= other;
                return *this;
            }

            friend std::ostream & operator<<(std::ostream& os, const se & other);
        };

        se operator*(const Real auto other, const se & rhs)
        {
            /*! \f{equation*}{ (\mathbb{R}, \mathfrak{se}) \rightarrow \mathfrak{se} \f}
            *
            * Scalar product.
            */

            return rhs*other;
        }

        std::ostream & operator<<(std::ostream & os, const se & other)
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