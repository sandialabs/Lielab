#ifndef _LIELAB_DOMAIN_rn_HPP
#define _LIELAB_DOMAIN_rn_HPP

using lielab::abstract::Real;
using lielab::abstract::Imaginary;

namespace lielab
{
    namespace domain
    {
        class rn
        {
            public:

            /*! Flag specifying whether or not the algebra is abelian. */
            static constexpr bool abelian = true;
            size_t shape = 0;
            
            /*! Internal data for the algebra. */
            Eigen::VectorXd _data;

            rn() : rn(0)
            {
                /*! \f{equation*}{() \rightarrow \mathfrak{rn} \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     lielab::domain::rn x, y, z;
                * 
                */

            }

            rn(const size_t shape)
            {
                /*! \f{equation*}{(\mathbb{Z}) \rightarrow \mathfrak{rn} \f}
                *
                * Constructor instantiating an \f$\mathfrak{rn}\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     lielab::domain::rn x(3), y(4), z(5);
                * 
                * @param[in] shape The shape of the data matrix.
                */

                this->shape = shape;
                this->_data = Eigen::VectorXd::Zero(std::max((int)shape-1, 0));
            }

            rn(std::initializer_list<double> other)
            {
                /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{rn} \f}
                *
                * Constructor instantiating an \f$\mathfrak{rn}\f$ object from either a
                * \f$n \times 1\f$ real vector.
                *
                * @param[in] other The object to instantiate from as a real vector.
                */

                this->set_vector(Eigen::VectorXd{std::move(other)});
            }

            template<typename OtherDerived>
            rn(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{(\mathbb{R}^{n \times n}) &\rightarrow& \mathfrak{rn} \\ (\mathbb{R}^{n \times 1}) &\rightarrow& \mathfrak{rn} \f}
                *
                * Constructor instantiating an \f$\mathfrak{rn}\f$ object from either an
                * \f$n \times n\f$ real matrix or \f$n \times 1\f$ real vector.
                *
                * @param[in] other The object to instantiate from as a real matrix.
                */

                if (other.rows() == other.cols() && other.rows() != 1)
                {
                    this->shape = other.rows();
                    this->_data = Eigen::VectorXd::Zero(shape-1);
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
                    throw Errorx("Input data to rn malformed.");
                }
            }

            template<typename OtherDerived>
            rn & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{eqnarray*}{ \mathfrak{rn} &:= \mathbb{R}^{n \times n} \\ &:= \mathbb{R}^{n \times 1} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `rn`.
                */

                if (other.rows() == other.cols() && other.rows() != 1)
                {
                    this->shape = other.rows();
                    this->_data = Eigen::VectorXd::Zero(shape-1);
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
                    throw Errorx("Input data to rn malformed.");
                }
                return *this;
            }

            static rn basis(const size_t i, const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{rn} \f}
                *
                * Returns the i'th basis element of the rn algebra.
                * 
                * @param[in] i The basis vector.
                * @param[in] shape The shape of the algebra.
                * @param[out] out The rn element. 
                */

                rn out(shape);
                const size_t _dim = out.get_dimension();
                Eigen::VectorXd _u(_dim);
                _u = Eigen::VectorXd::Zero(_dim);
                _u(i) = 1.0;
                out.set_vector(_u);
                return out;
            }

            static Eigen::MatrixXd project(const Eigen::MatrixXd & other)
            {
                /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{rn} \f}
                *
                * Projects a matrix suitable for data.
                */

                if (other.rows() != other.cols())
                {
                    throw Errorx("Size of the matrix must be square.");
                }

                Eigen::MatrixXd out = other;

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

                return this->shape - 1;
            }

            Eigen::VectorXd get_vector() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
                * 
                * Returns the vector representation.
                */

                return _data;
            }

            Eigen::MatrixXd get_ados_representation() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times n} \f}
                * 
                * Returns a matrix representation.
                */

                Eigen::MatrixXd out = Eigen::MatrixXd::Zero(shape, shape);

                for (int ii = 0; ii < shape-1; ii++)
                {
                    out(ii, shape-1) = _data(ii);
                }

                return out;
            }

            void set_vector(const Eigen::VectorXd & vector)
            {
                /*! \f{equation*}{ \mathfrak{rn} := \mathbb{R}^{n \times 1} \f}
                * 
                * @param[in] vector An Eigen::VectorXd to assign.
                */

                this->_data = vector;
                this->shape = this->_data.size() + 1;
            }

            double operator()(const size_t index) const
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
                *
                * Gets a value in the column vector representation.
                */

                Eigen::VectorXd vector = get_vector();
                return _data(index);
            }

            double & operator()(const size_t index)
            {
                /*! \f{equation*}{ \mathfrak{rn}(\mathbb{Z}) := \mathbb{R} \f}
                *
                * Assignment of a value in the column vector representation.
                */
                
                return _data(index);
            }

            double operator()(const size_t index1, const size_t index2) const
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathbb{R} \f}
                *
                * Gets a value in the square matrix representation.
                */
                
                Eigen::MatrixXd mat = get_ados_representation();
                return mat(index1, index2);
            }

            rn operator+(const rn & other) const
            {
                /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
                *
                * Addition of two vectors in the algebra.
                */
                
                assert(this->shape == other.shape);
                return this->_data + other._data;
            }

            rn & operator+=(const rn & other)
            {
                /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
                *
                * In place addition of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                this->_data += other._data;
                return *this;
            }

            rn operator-(const rn & other) const
            {
                /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
                *
                * Subtraction of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                return this->_data - other._data;
            }

            rn & operator-=(const rn & other)
            {
                /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
                *
                * In place subtraction of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                this->_data -= other._data;
                return *this;
            }

            rn operator-() const
            {
                /*! \f{equation*}{ (\mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
                *
                * Unary negative of the vector.
                */

                return -this->_data;
            }

            rn operator*(const Real auto other) const
            {
                /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
                *
                * Scalar product.
                */

                return this->_data * other;
            }

            rn & operator*=(const Real auto other)
            {
                /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
                *
                * In place scalar product.
                */

                this->_data *= other;
                return *this;
            }

            rn operator*(const rn & other) const
            {
                /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathbb{R}^{n \times n} \notin \mathfrak{rn} \f}
                *
                * Vector product.
                */

                assert(this->shape == other.shape);
                return Eigen::VectorXd::Zero(shape-1);
            }

            rn & operator*=(const rn & other)
            {
                /*! \f{equation*}{ (\mathfrak{rn}, \mathfrak{rn}) \rightarrow \mathbb{R}^{n \times n} \notin \mathfrak{rn} \f}
                *
                * In place vector product.
                */

                assert(this->shape == other.shape);
                this->_data *= 0.0;
                return *this;
            }

            rn operator/(const Real auto other) const
            {
                /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
                *
                * Scalar division.
                */

                return this->_data / other;
            }

            rn & operator/=(const Real auto other)
            {
                /*! \f{equation*}{ (\mathfrak{rn}, \mathbb{R}) \rightarrow \mathfrak{rn} \f}
                *
                * In place scalar division.
                */

                this->_data /= other;
                return *this;
            }

            friend std::ostream& operator<<(std::ostream & os, const rn & other);
        };

        rn operator*(const Real auto other, const rn & rhs)
        {
            /*! \f{equation*}{ (\mathbb{R}, \mathfrak{rn}) \rightarrow \mathfrak{rn} \f}
            *
            * Scalar product.
            */

            return rhs*other;
        }

        std::ostream& operator<<(std::ostream & os, const rn & other)
        {
            /*!
            * Overloads the "<<" stream insertion operator.
            */
            
            os << static_cast<const Eigen::VectorXd>(other._data);
            return os;
        }
    }
}

#endif
