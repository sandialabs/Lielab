#ifndef _LIELAB_DOMAIN_LIEALGEBRAS_su_HPP
#define _LIELAB_DOMAIN_LIEALGEBRAS_su_HPP

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace Lielab
{
    namespace domain
    {
        class su
        {
            public:
            static constexpr bool abelian = false;
            size_t shape = 0;

            Eigen::MatrixXcd _data;

            su() : su(0)
            {
                /*! \f{equation*}{ () \rightarrow \mathfrak{su} \f}
                * 
                * Empty initialization function. Enables instantiation like:
                * 
                *     Lielab::domain::su x, y, z;
                * 
                */

            }

            su(const size_t n)
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathfrak{su} \f}
                *
                * Constructor instantiating an \f$\mathfrak{su}\f$ object.
                * 
                * Enables instantiation like:
                * 
                *     Lielab::domain::su x(2), y(3), z(4);
                * 
                * @param[in] shape The shape of the data matrix.
                * 
                * TODO: Only shape 2 is implemented.
                */
                
                this->_data = Eigen::MatrixXcd::Zero(n,n);
                this->shape = n;
            }

            su(std::initializer_list<double> other)
            {
                /*! \f{equation}{(\mathbb{R}^{n \times 1}) \rightarrow \mathfrak{su} \f}
                *
                * Constructor instantiating an \f$\mathfrak{su}\f$ object from either a
                * \f$n \times 1\f$ real vector.
                *
                * @param[in] other The object to instantiate from as a real vector.
                */

                this->set_vector(Eigen::VectorXd{std::move(other)});
            }

            template<typename OtherDerived>
            su(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{equation}{ (\mathbb{C}^{n \times n}) \rightarrow \mathfrak{su} \f}
                *
                * Constructor instantiating an \f$\mathfrak{su}\f$ object from an
                * \f$n \times n\f$ real matrix.
                *
                * @param[in] other The object to instantiate from.
                */
                
                this->_data = Eigen::MatrixXcd(other);
                this->shape = _data.rows();
            }

            template<typename OtherDerived>
            su & operator=(const Eigen::MatrixBase<OtherDerived> & other)
            {
                /*! \f{equation*}{ \mathfrak{su} := \mathbb{C}^{n \times n} \f}
                * 
                * Overload of the assignment operator. Allows Eigen Matrix data to be directly assigned to `su`.
                */
                
                this->_data = Eigen::MatrixXcd(other);
                this->shape = _data.rows();
                return *this;
            }

            static su basis(const size_t i, const size_t shape)
            {
                /*! \f{equation*}{ (\mathbb{Z}, \mathbb{Z}) \rightarrow \mathfrak{su} \f}
                *
                * Returns the i'th basis element of the su algebra.
                * 
                * @param[in] i The basis vector.
                * @param[in] shape The shape of the algebra.
                * @param[out] out The su element.
                */

                su out(shape);
                const size_t _dim = out.get_dimension();
                Eigen::VectorXd _u(_dim);
                _u = Eigen::VectorXd::Zero(_dim);
                _u(i) = 1.0;
                out.set_vector(_u);
                return out;
            }

            // TODO: Project function

            size_t get_dimension() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{Z} \f}
                * 
                * Gets the dimension of the algebra.
                */

                return this->shape * this->shape - 1;
            }

            Eigen::VectorXd get_vector() const
            {
                /*! \f{equation*}{ () \rightarrow \mathbb{R}^{n \times 1} \f}
                * 
                * Returns the vector representation.
                */

                if (this->shape <= 1)
                {
                    return Eigen::VectorXd::Zero(0);
                }

                const size_t dim = this->get_dimension();
                Eigen::VectorXd out = Eigen::VectorXd::Zero(dim);

                // General case of su(2+). Use's the Generalized Gell-Mann matrices
                size_t k = 0;

                // Symmetric
                for (size_t jj = 1; jj < this->shape; jj++)
                {
                    for (size_t ii = 0; ii < jj; ii++)
                    {
                        out(k) = std::imag(this->_data(ii, jj) + this->_data(jj, ii))/2.0;
                        k++;
                    }
                }

                // Anti-symmetric
                for (size_t jj = 1; jj < this->shape; jj++)
                {
                    for (size_t ii = 0; ii < jj; ii++)
                    {
                        out(k) = std::real(this->_data(jj, ii) - this->_data(ii, jj))/2.0;
                        k++;
                    }
                }

                // Diagonal
                size_t zz = this->shape;
                size_t yy = this->shape - 1;
                k = out.size() - 1;
                Eigen::MatrixXcd temp = this->get_matrix();
                for (size_t yy = this->shape - 1; yy >= 1; yy--)
                {
                    const double multiplier = std::sqrt(2.0/((zz-1)*(zz)));
                    out(k) = -std::imag(temp(yy, yy))/(multiplier*(zz - 1));

                    for (size_t ii = 0; ii < yy; ii++)
                    {
                        temp(ii, ii) -= std::imag(multiplier*out(k));
                    }

                    zz--;
                    k--;
                }

                return out;
            }

            Eigen::MatrixXcd get_matrix() const
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
                /*! \f{equation*}{ \mathfrak{su} := \mathbb{R}^{n \times 1} \f}
                * 
                * @param[in] vector An Eigen::VectorXd to assign.
                * TODO: Set shape.
                */
                
                // Error check for vector size
                assert((this->get_dimension() == vector.size()) && "su.set_vector(): Input dimension mismatch.");

                this->_data = Eigen::MatrixXcd::Zero(this->shape, this->shape);

                // su(0) and su(1) are 0-dimensional. Return 0's.
                if (this->shape <= 1)
                {
                    return;
                }

                // General case of su(2+). Use's the Generalized Gell-Mann matrices
                size_t k = 0;

                // Symmetric
                for (size_t jj = 1; jj < this->shape; jj++)
                {
                    for (size_t ii = 0; ii < jj; ii++)
                    {
                        this->_data(ii, jj) = std::complex<double>(0.0, vector(k));
                        this->_data(jj, ii) = std::complex<double>(0.0, vector(k));
                        k++;
                    }
                }

                // Anti-symmetric
                for (size_t jj = 1; jj < this->shape; jj++)
                {
                    for (size_t ii = 0; ii < jj; ii++)
                    {
                        this->_data(ii, jj) += std::complex<double>(-vector(k), 0.0);
                        this->_data(jj, ii) += std::complex<double>(vector(k), 0.0);
                        k++;
                    }
                }

                // Diagonal
                size_t zz = 2;
                while (k < vector.size())
                {
                    const double multiplier = std::sqrt(2.0/((zz-1)*(zz)));

                    for (size_t ii = 0; ii < zz; ii++)
                    {
                        if (ii == (zz - 1))
                        {
                            this->_data(ii, ii) -= std::complex<double>(0.0, multiplier*(zz - 1)*vector(k));
                        }
                        else
                        {
                            this->_data(ii, ii) += std::complex<double>(0.0, multiplier*vector(k));
                        }
                    }

                    zz++;
                    k++;
                }
            }

            double operator()(const size_t index) const
            {
                /*! \f{equation*}{ (\mathbb{Z}) \rightarrow \mathbb{R} \f}
                *
                * Gets a value in the column vector representation.
                */

                return get_vector()(index);
            }

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
                /*! \f{equation*}{ \mathfrak{su}(\mathbb{Z}, \mathbb{Z}) := \mathbb{C} \f}
                *
                * Assignment of a value in the square matrix representation.
                */

                return _data(index1, index2);
            }

            su & operator+=(const su & other)
            {
                /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
                *
                * In place addition of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                this->_data += other._data;
                return *this;
            }

            su & operator-=(const su & other)
            {
                /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
                *
                * In place subtraction of two vectors in the algebra.
                */

                assert(this->shape == other.shape);
                this->_data -= other._data;
                return *this;
            }

            su operator-() const
            {
                /*! \f{equation*}{ (\mathfrak{su}) \rightarrow \mathfrak{su} \f}
                *
                * Unary negative of the vector.
                */

                Eigen::MatrixXcd out = - _data;
                return out;
            }

            su operator*(const Imaginary auto other) const
            {
                /*! \f{equation*}{ (\mathfrak{su}, \mathbb{C}) \rightarrow \mathfrak{su} \f}
                *
                * Scalar product.
                */

                Eigen::MatrixXcd out = _data * other;
                return out;
            }

            su & operator*=(const Imaginary auto other)
            {
                /*! \f{equation*}{ (\mathfrak{su}, \mathbb{C}) \rightarrow \mathfrak{su} \f}
                *
                * In place scalar product.
                */

                this->_data *= other;
                return *this;
            }

            su & operator*=(const su & other)
            {
                /*! \f{equation*}{ (\mathfrak{su}, \mathfrak{su}) \rightarrow \mathbb{C}^{n \times n} \notin \mathfrak{su} \f}
                *
                * In place vector product.
                */

                assert(this->shape == other.shape);
                this->_data *= other._data;
                return *this;
            }

            su operator/(const Imaginary auto other) const
            {
                /*! \f{equation*}{ (\mathfrak{su}, \mathbb{C}) \rightarrow \mathfrak{su} \f}
                *
                * Scalar division.
                */

                Eigen::MatrixXcd out = _data / other;
                return out;
            }

            su & operator/=(const Imaginary auto other)
            {
                /*! \f{equation*}{ (\mathfrak{su}, \mathbb{C}) \rightarrow \mathfrak{su} \f}
                *
                * In place scalar division.
                */

                this->_data /= other;
                return *this;
            }

            friend std::ostream & operator<<(std::ostream& os, const su & other);
        };

        su operator*(const Imaginary auto other, const su & rhs)
        {
            /*! \f{equation*}{ (\mathbb{C}, \mathfrak{su}) \rightarrow \mathfrak{su} \f}
            *
            * Scalar product.
            */

            return rhs*other;
        }

        std::ostream & operator<<(std::ostream & os, const su & other)
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
