#ifndef _LIELAB_DOMAIN_LIEALGEBRAS_sp_HPP
#define _LIELAB_DOMAIN_LIEALGEBRAS_sp_HPP

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
    namespace domain
    {
        class sp
        {
            public:

            /*! Flag specifying whether or not the algebra is abelian. */
            static constexpr bool abelian = false;
            size_t shape = 0;

            Eigen::MatrixXd _data;

            sp() : sp(0)
            {
                /*! \f{equation*}{ () \rightarrow \mathfrak{sp} \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     Lielab::domain::sp x, y, z;
                * 
                */

            }

            sp(const size_t n)
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{sp} \f}
                *
                * Constructor instantiating an \f$\mathfrak{sp}\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     Lielab::domain::sp x(2), y(4), z(6);
                * 
                * @param[in] shape The shape of the data matrix.
                */

                if (n % 2 != 0)
                {
                    throw Errorx("Shape of sp must be even dimensional.");
                }

                this->_data = Eigen::MatrixXd::Zero(n,n);
                this->shape = n;
            }

            sp(std::initializer_list<double> other)
            {
                /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{sp} \f}
                *
                * Constructor instantiating an \f$\mathfrak{sp}\f$ object from either a
                * \f$n \times 1\f$ real vector.
                *
                * @param[in] other The object to instantiate from as a real vector.
                */

                this->set_vector(Eigen::VectorXd{std::move(other)});
            }

            template<typename OtherDerived>
            sp(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{equation}{ (\mathbb{R}^{n \times n}) \rightarrow \mathfrak{sp} \f}
                *
                * Constructor instantiating an \f$\mathfrak{sp}\f$ object from an
                * \f$n \times n\f$ real matrix.
                *
                * @param[in] other The object to instantiate from.
                */
                
                if (other.rows() % 2 != 0)
                {
                    throw Errorx("Shape of sp must be even dimensional.");
                }

                if (other.rows() != other.cols())
                {
                    throw Errorx("Shape of sp must be square.");
                }

                this->_data = Eigen::MatrixXd(other);
                this->shape = _data.rows();
            }

            template<typename OtherDerived>
            sp & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{equation*}{ \mathfrak{sp} := \mathbb{R}^{n \times n} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `sp`.
                */

                if (other.rows() % 2 != 0)
                {
                    throw Errorx("Shape of sp must be even dimensional.");
                }

                if (other.rows() != other.cols())
                {
                    throw Errorx("Shape of sp must be square.");
                }
                
                this->_data = Eigen::MatrixXd(other);
                this->shape = _data.rows();
                return *this;
            }

            static sp basis(const size_t i, const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{sp} \f}
                *
                * Returns the i'th basis element of the sp algebra.
                * 
                * @param[in] i The basis vector.
                * @param[in] shape The shape of the algebra.
                * @param[out] out The rn element.
                */

                sp out(shape);
                const size_t _dim = out.get_dimension();
                Eigen::VectorXd _u(_dim);
                _u = Eigen::VectorXd::Zero(_dim);
                _u(i) = 1.0;
                out.set_vector(_u);
                return out;
            }

            static Eigen::MatrixXd project(const Eigen::MatrixXd & other)
            {
                /*! \f{equation*}{ (\mathbb{R}^{n \times n}) \rightarrow \mathbb{R}^{n \times n} \in \mathfrak{sp} \f}
                *
                * Projects a matrix suitable for data.
                */

                if (other.rows() % 2 != 0)
                {
                    throw Errorx("Shape of sp must be even dimensional.");
                }

                if (other.rows() != other.cols())
                {
                    throw Errorx("Shape of sp must be square.");
                }

                const size_t half_shape = other.rows() / 2;
                Eigen::MatrixXd J = Eigen::MatrixXd::Zero(other.rows(), other.cols());
                J.block(0, half_shape, half_shape, half_shape) = -Eigen::MatrixXd::Identity(half_shape, half_shape);
                J.block(half_shape, 0, half_shape, half_shape) = Eigen::MatrixXd::Identity(half_shape, half_shape);
                Eigen::MatrixXd temp = -J*other;
                return J*(temp + temp.transpose())/2.0;
            }

            size_t get_dimension() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
                * 
                * Gets the dimension of the algebra.
                */

                return this->shape * (this->shape + 1) / 2;
            }

            Eigen::VectorXd get_vector() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
                * 
                * Returns the vector representation.
                */

                Eigen::VectorXd out(get_dimension());
                int k = 0;

                for (int ii = 0; ii < shape/2; ii++)
                {
                    for (int jj = 0; jj < shape/2; jj++)
                    {
                        out(k) = _data(ii,jj);
                        k++;
                    }
                }

                for (int ii = 0; ii < shape/2; ii++)
                {
                    for (int jj = 0; jj < shape/2; jj++)
                    {
                        if (jj >= ii)
                        {
                            out(k) = _data(ii, jj+shape/2);
                            k++;
                        }
                    }
                }

                for (int ii = 0; ii < shape/2; ii++)
                {
                    for (int jj = 0; jj < shape/2; jj++)
                    {
                        if (jj >= ii)
                        {
                            out(k) = _data(ii+shape/2, jj);
                            k++;
                        }
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
                /*! \f{equation*}{ \mathfrak{sp} := \mathbb{R}^{n \times 1} \f}
                * 
                * @param[in] vector An Eigen::VectorXd to assign.
                * TODO: Set shape.
                */

                Eigen::MatrixXd A(shape/2, shape/2), B(shape/2, shape/2), C(shape/2, shape/2);
                int k = 0;

                for (int ii = 0; ii < shape/2; ii++)
                {
                    for (int jj = 0; jj < shape/2; jj++)
                    {
                        A(ii, jj) = vector(k);
                        k++;
                    }
                }

                for (int ii = 0; ii < shape/2; ii++)
                {
                    for (int jj = 0; jj < shape/2; jj++)
                    {
                        if (jj >= ii)
                        {
                            B(ii, jj) = vector(k);
                            k++;
                        }

                        if (ii != jj)
                        {
                            B(jj, ii) = B(ii, jj);
                        }
                    }
                }

                for (int ii = 0; ii < shape/2; ii++)
                {
                    for (int jj = 0; jj < shape/2; jj++)
                    {
                        if (jj >= ii)
                        {
                            C(ii, jj) = vector(k);
                            k++;
                        }

                        if (ii != jj)
                        {
                            C(jj, ii) = C(ii, jj);
                        }
                    }
                }

                _data.block(0,0,shape/2,shape/2) = A;
                _data.block(0, shape/2, shape/2, shape/2) = B;
                _data.block(shape/2, 0, shape/2, shape/2) = C;
                _data.block(shape/2, shape/2, shape/2, shape/2) = -A.transpose();
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
                /*! \f{equation*}{ \mathfrak{sp}(\mathbb{Z}, \mathbb{Z}) := \mathbb{R} \f}
                *
                * Assignment of a value in the square matrix representation.
                */

                return _data(index1, index2);
            }

            sp & operator+=(const sp & other)
            {
                /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
                *
                * In place addition of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                this->_data += other._data;
                return *this;
            }

            sp & operator-=(const sp & other)
            {
                /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
                *
                * In place subtraction of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                this->_data -= other._data;
                return *this;
            }

            sp operator-() const
            {
                /*! \f{equation*}{ (\mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
                *
                * Unary negative of the vector.
                */

                Eigen::MatrixXd out = - _data;
                return out;
            }

            sp operator*(const Real auto other) const
            {
                /*! \f{equation*}{ (\mathfrak{sp}, \mathbb{R}) \rightarrow \mathfrak{sp} \f}
                *
                * Scalar product.
                */

                Eigen::MatrixXd out = _data * other;
                return out;
            }

            sp & operator*=(const Real auto other)
            {
                /*! \f{equation*}{ (\mathfrak{sp}, \mathbb{R}) \rightarrow \mathfrak{sp} \f}
                *
                * In place scalar product.
                */

                this->_data *= other;
                return *this;
            }

            sp & operator*=(const sp & other)
            {
                /*! \f{equation*}{ (\mathfrak{sp}, \mathfrak{sp}) \rightarrow \mathbb{R}^{n \times n} \notin \mathfrak{sp} \f}
                *
                * In place vector product.
                */

                assert(this->shape == other.shape);
                this->_data *= other._data;
                return *this;
            }

            sp operator/(const Real auto other) const
            {
                /*! \f{equation*}{ (\mathfrak{sp}, \mathbb{R}) \rightarrow \mathfrak{sp} \f}
                *
                * Scalar division.
                */

                Eigen::MatrixXd out = _data / other;
                return out;
            }

            sp & operator/=(const Real auto other)
            {
                /*! \f{equation*}{ (\mathfrak{sp}, \mathbb{R}) \rightarrow \mathfrak{sp} \f}
                *
                * In place scalar division.
                */

                this->_data /= other;
                return *this;
            }

            friend std::ostream& operator<<(std::ostream & os, const sp & other);
        };

        sp operator*(const Real auto other, const sp & rhs)
        {
            /*! \f{equation*}{ (\mathbb{R}, \mathfrak{sp}) \rightarrow \mathfrak{sp} \f}
            *
            * Scalar product.
            */

            return rhs*other;
        }

        std::ostream & operator<<(std::ostream & os, const sp & other)
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
